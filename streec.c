/**********  simulate_coalescent.c **********************************
*
*	This subroutine uses a Monte Carlo algorithm described in
*	Hudson,R. 1983. Theor.Pop.Biol. 23: 183-201, to produce
*	a history of a random sample of gametes under a neutral
*	Wright-Fisher model with recombination and geographic structure.
*	Input parameters
*	are the sample size (nsam), the number of sites between
*	which recombination can occur (nsites), and the recombination
*	rate between the ends of the gametes (r). The function returns
*	nsegs, the number of segments the gametes were broken into
*	in tracing back the history of the gametes.  The histories of
*	these segments are passed back to the calling function in the
*	array of structures seglst[]. An element of this array,  seglst[i],
* 	consists of three parts: (1) beg, the starting point of
*	of segment i, (2) ptree, which points to the first node of the
*	tree representing the history of the segment, (3) next, which
*	is the index number of the next segment.
*	     A tree is a contiguous set of 2*nsam nodes. The first nsam
*	nodes are the tips of the tree, the sampled gametes.  The other
*	nodes are the nodes ancestral to the sampled gametes. Each node
*	consists of an "abv" which is the number of the node ancestral to
*	that node, an "ndes", which is not used or assigned to in this routine,
*	and a "time", which is the time (in units of 4N generations) of the
*	node. For the tips, time equals zero.
*	Returns a pointer to an array of segments, seglst.

**************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ms.h"

#define NL putchar('\n')
#define size_t unsigned

#define MIN(x, y) ((x) < (y) ? (x) : (y))

#define ERROR(message) fprintf(stderr, message), NL, exit(1)

#define SEGINC 80
#define CHROMOSOME_INC 20

int nchrom, begs, nsegs;
long nchrom_links;
static int *nnodes = NULL;
double t, cleft, pc, lnpc;

static unsigned seglimit = SEGINC;
static unsigned maxchr;

static struct seg *buf = NULL ;
static int buf_cap = 0 ;

struct seg {
    int beg;
    int end;
    int desc;
};

struct chromo {
    int nseg;
    int pop;
    struct seg *pseg;
};

static struct chromo *chrom = NULL;

struct node *ptree1, *ptree2, *ptree0;

static struct segl *seglst = NULL;

/* ---- funciones externas (rand1.c, ms.c) ---- */
extern double ran1(void);

struct segl *simulate_coalescent(const struct c_params *cp, int *pnsegs);
static void allocate_buffers(int nsam);
static void init_chromo(int ind, int desc, int pop, int nsites);
static void check_isolated_populations(const int *config, int npop, double mig, const struct devent *nextevent);
static void migration_event(int *config, double *const *migration_matrix, int npop, double mig);
static int recombination_event(int nsam);
static int conv_left(int nsam);
static int conv_inside(int nsam, int nsites);  // Conversion left y conversion inside
static int find_crossover_chrom(long spot, int *pis, int *pend);
static int crossover(int nsam, int ic, int is);
static void insert_new_tree(int i, int nsam);
static int merge_chrom(int nsam, int nsites, int c1, int c2);
static void remove_chrom_links(int c);
static void add_chrom_links(int c);
static void move_last_chrom_to(int c);
static int chromo_has_segment(int start, int c, int *psg);
static void pick_2_chrom(int pop, const int config[], int *pc1, int *pc2);
static int next_chrom_in_pop(int pop, int target, int i, int count);
static int pick_2_random_numbers(int n, int *i, int *j);
static int chrom_links(int c);

/*
 *	@brief Es la función generadora de arboles coalescentes.
 *
 *	En el bucle determina de manera probabilística cual será el proximo evento a ejecutar (a excepción de los
 *  eventos demográficos que son determinístas). Entre los eventos se encuentran los demográficos, recombinación, 
 *  conversión génica, coalescencia y eventos migratorios. El bucle sigue mientras quede más de un linaje, contando
 *  las muestras de ADN antiguo que todavía no entraron en la simulación. Finaliza cuando queda un solo linaje.
 * 
 *  Al terminar la ejecución se obtiene como resultado la lista de arboles genealógicos correspondientes a cada 
 *  segmento del locus.
 * 
 *	@param[in] cp : Son los parámetros de coalescencia que incluyen también la lista de eventos demográficos. 
 *	@param[out]	pnsegs : Es la cantidad de segmentos en los que se dividió el locus.
 *	@return	seglst : Es la lista de segmentos del locus, cada uno con su árbol genealógico asociado de 2*nsam nodos.
 *
 */
struct segl *simulate_coalescent(const struct c_params *cp, int *pnsegs)
{
    int i, j, k, seg, dec, pop, pop2, c1, c2, ind, rchrom, nsamin, numanc, num;
    int *config;
    double ttemp, rft, clefta, tmin, p;
    double prec, cin, prect, mig, ran, coal_prob, rdum, arg;
    char event;
    int eflag, cpop, ic;
    int nsam, npop, nsites ;
    const int *inconfig ;
    double r, f, rf, track_len, **migration_matrix;
    double *size, *alphag, *tlast;
    const struct devent *nextevent;

    nsam = cp->nsam;                   // Número total de muestras
    nsamin = cp->nsamin; /* adna */    // Número de muestras en el presente cuando t = 0
    numanc = nsam - nsamin; /* adna */ // Número de muestras de ADN antiguo cuando t = x
                                       // especificado por el usuario con -eA
    npop = cp->npop;                   // Número de subpoblaciones
    nsites = cp->nsites;               // Número de sitios donde pueden ocurrir recombinación
    inconfig = cp->config;             // Configuración de muestras por subpoblación
    r = cp->r;                         // Tasa de recombinación
    f = cp->f;                         // Otra tasa
    track_len = cp->track_len;         // Longitud promedio del track

    migration_matrix = (double **)malloc((unsigned)npop * sizeof(double *));

    for (i = 0; i < npop; i++) {
        migration_matrix[i] = (double *)malloc((unsigned)npop * sizeof(double));
        for (j = 0; j < npop; j++)
            migration_matrix[i][j] = (cp->mig_mat)[i][j];
    }
    nextevent = cp->deventlist;

    allocate_buffers(nsam);

    config = (int *)malloc((unsigned)((npop + 1) * sizeof(int)));
    if (config == NULL)
        perror("malloc error. segtre.");
    size = (double *)malloc((unsigned)((npop) * sizeof(double)));
    if (size == NULL)
        perror("malloc error. segtre.");
    alphag = (double *)malloc((unsigned)((npop) * sizeof(double)));
    if (alphag == NULL)
        perror("malloc error. segtre.");
    tlast = (double *)malloc((unsigned)((npop) * sizeof(double)));
    if (alphag == NULL)
        perror("malloc error. segtre.");

    for (pop = 0; pop < npop; pop++) {
        config[pop] = inconfig[pop];
        size[pop] = (cp->size)[pop];
        alphag[pop] = (cp->alphag)[pop];
        tlast[pop] = 0.0;
    }

    for (pop = ind = 0; pop < npop; pop++)          // Itera por población y en cada población, por individuo
        for (j = 0; j < inconfig[pop]; j++, ind++)
            init_chromo(ind, ind, pop, nsites);

    seglst[0].beg = 0;
    if (!(seglst[0].ptree = (struct node *)calloc((unsigned)(2 * nsam), sizeof(struct node))))
        perror("calloc error. se2");
    nnodes[0] = nsam - 1;

    nchrom = nsamin;                          /* adna */
    nchrom_links = ((long)(nsamin)) * (nsites - 1); /* adna */
    nsegs = 1;
    t = 0.;
    r /= (nsites - 1);
    if (f > 0.0)
        pc = (track_len - 1.0) / track_len;
    else
        pc = 1.0;
    lnpc = log(pc);
    cleft = nsamin * (1.0 - pow(pc, (double)(nsites - 1))); /* adna */
    if (r > 0.0)
        rf = r * f;
    else
        rf = f / (nsites - 1);
    rft = rf * track_len;

    /* Main loop */

    while ((nchrom + numanc) > 1) { /* adna */
        prec = nchrom_links * r;
        cin = nchrom_links * rf;
        clefta = cleft * rft;
        prect = prec + cin + clefta;
        mig = 0.0; // Migración

        for (i = 0; i < npop; i++)
            mig += config[i] * migration_matrix[i][i];
        check_isolated_populations(config, npop, mig, nextevent);
        eflag = 0;

        if (prect > 0.0) { /* cross-over or gene conversion */
            while ((rdum = ran1()) == 0.0) ;
            ttemp = -log(rdum) / prect;
            if ((eflag == 0) || (ttemp < tmin)) {
                tmin = ttemp;
                event = 'r';
                eflag = 1;
            }
        }
        if (mig > 0.0) { /* migration   */
            while ((rdum = ran1()) == 0.0) ;
            ttemp = -log(rdum) / mig;
            if ((eflag == 0) || (ttemp < tmin)) {
                tmin = ttemp;
                event = 'm';
                eflag = 1;
            }
        }

        for (pop = 0; pop < npop; pop++) { /* coalescent */
            coal_prob = ((double)config[pop]) * (config[pop] - 1.);
            if (coal_prob > 0.0) {
                while ((rdum = ran1()) == .0) ;
                if (alphag[pop] == 0) {
                    ttemp = -log(rdum) * size[pop] / coal_prob;
                    if ((eflag == 0) || (ttemp < tmin)) {
                        tmin = ttemp;
                        event = 'c';
                        eflag = 1;
                        cpop = pop;
                    }
                } else {
                    arg = 1. - alphag[pop] * size[pop] * exp(-alphag[pop] * (t - tlast[pop])) * log(rdum) / coal_prob;
                    if (arg > 0.0) { /*if arg <= 0,  no coalescent within interval */
                        ttemp = log(arg) / alphag[pop];
                        if ((eflag == 0) || (ttemp < tmin)) {
                            tmin = ttemp;
                            event = 'c';
                            eflag = 1;
                            cpop = pop;
                        }
                    }
                }
            }
        }

        if ((eflag == 0) && (nextevent == NULL)) {
            fprintf(stderr,
                    " infinite time to next event. Negative growth rate in last time interval or "
                    "non-communicating subpops.\n");
            exit(0);
        }
        if (((eflag == 0) && (nextevent != NULL)) ||
            ((nextevent != NULL) && ((t + tmin) >= nextevent->time))) {
            t = nextevent->time;
            switch (nextevent->detype) {
            case 'N':
                for (pop = 0; pop < npop; pop++) {
                    size[pop] = nextevent->paramv;
                    alphag[pop] = 0.0;
                }
                nextevent = nextevent->nextde;
                break;
            case 'n':
                size[nextevent->popi] = nextevent->paramv;
                alphag[nextevent->popi] = 0.0;
                nextevent = nextevent->nextde;
                break;
            case 'G':
                for (pop = 0; pop < npop; pop++) {
                    size[pop] = size[pop] * exp(-alphag[pop] * (t - tlast[pop]));
                    alphag[pop] = nextevent->paramv;
                    tlast[pop] = t;
                }
                nextevent = nextevent->nextde;
                break;
            case 'g':
                pop = nextevent->popi;
                size[pop] = size[pop] * exp(-alphag[pop] * (t - tlast[pop]));
                alphag[pop] = nextevent->paramv;
                tlast[pop] = t;
                nextevent = nextevent->nextde;
                break;
            case 'M':
                for (pop = 0; pop < npop; pop++)
                    for (pop2 = 0; pop2 < npop; pop2++)
                        migration_matrix[pop][pop2] = (nextevent->paramv) / (npop - 1.0);
                for (pop = 0; pop < npop; pop++)
                    migration_matrix[pop][pop] = nextevent->paramv;
                nextevent = nextevent->nextde;
                break;
            case 'a':
                for (pop = 0; pop < npop; pop++)
                    for (pop2 = 0; pop2 < npop; pop2++)
                        migration_matrix[pop][pop2] = (nextevent->mat)[pop][pop2];
                nextevent = nextevent->nextde;
                break;
            case 'm':
                i = nextevent->popi;
                j = nextevent->popj;
                migration_matrix[i][i] += nextevent->paramv - migration_matrix[i][j];
                migration_matrix[i][j] = nextevent->paramv;
                nextevent = nextevent->nextde;
                break;
            case 'j': /* merge pop i into pop j  (join) */
                i = nextevent->popi;
                j = nextevent->popj;
                config[j] += config[i];
                config[i] = 0;
                for (ic = 0; ic < nchrom; ic++)
                    if (chrom[ic].pop == i)
                        chrom[ic].pop = j;
                /*  the following was added 19 May 2007 */
                for (k = 0; k < npop; k++) {
                    if (k != i) {
                        migration_matrix[k][k] -= migration_matrix[k][i];
                        migration_matrix[k][i] = 0.;
                    }
                }
                /* end addition */
                nextevent = nextevent->nextde;
                break;
            case 'A': /* ancient dna sample */
                num = nextevent->num;
                nchrom_links += ((long)(num)) * (nsites - 1);
                cleft += num * (1.0 - pow(pc, (double)(nsites - 1)));

                if ( (unsigned int)(nchrom + num) >= maxchr) {
                    maxchr += CHROMOSOME_INC;
                    chrom =
                        (struct chromo *)realloc(chrom, (unsigned)(maxchr * sizeof(struct chromo)));
                    if (chrom == NULL)
                        perror("malloc error. segtre2");
                }

                for (ind = nchrom, i = nsam - numanc; ind < nchrom + num; ind++, i++)
                    init_chromo(ind, i, nextevent->popi, nsites);
                ptree1 = seglst[0].ptree;
                for (i = nsam - numanc; i < nsam - numanc + num; i++) {
                    (ptree1 + i)->time = t;
                }
                nchrom += num;
                config[nextevent->popi] += num;
                numanc -= num;
                nextevent = nextevent->nextde;
                break;
            case 's': /*split  pop i into two;p is the proportion from pop i, and 1-p from pop
                         n+1  */
                i = nextevent->popi;
                p = nextevent->paramv;
                npop++;
                config = (int *)realloc(config, (unsigned)(npop * sizeof(int)));
                size = (double *)realloc(size, (unsigned)(npop * sizeof(double)));
                alphag = (double *)realloc(alphag, (unsigned)(npop * sizeof(double)));
                tlast = (double *)realloc(tlast, (unsigned)(npop * sizeof(double)));
                tlast[npop - 1] = t;
                size[npop - 1] = 1.0;
                alphag[npop - 1] = 0.0;
                migration_matrix = (double **)realloc(migration_matrix, (unsigned)(npop * sizeof(double *)));
                for (j = 0; j < npop - 1; j++)
                    migration_matrix[j] = (double *)realloc(migration_matrix[j], (unsigned)(npop * sizeof(double)));
                migration_matrix[npop - 1] = (double *)malloc((unsigned)(npop * sizeof(double)));
                for (j = 0; j < npop; j++)
                    migration_matrix[npop - 1][j] = migration_matrix[j][npop - 1] = 0.0;
                config[npop - 1] = 0;
                config[i] = 0;
                for (ic = 0; ic < nchrom; ic++) {
                    if (chrom[ic].pop == i) {
                        if (ran1() < p)
                            config[i]++;
                        else {
                            chrom[ic].pop = npop - 1;
                            config[npop - 1]++;
                        }
                    }
                }
                nextevent = nextevent->nextde;
                break;
            }
        } else {
            t += tmin;
            if (event == 'r') {
                if ((ran = ran1()) < (prec / prect)) { /*recombination */
                    rchrom = recombination_event(nsam);
                    config[chrom[rchrom].pop] += 1;
                } else if (ran < (prec + clefta) / (prect)) { /*  cleft event */
                    rchrom = conv_left(nsam);
                    config[chrom[rchrom].pop] += 1;
                } else { /* cin event */
                    rchrom = conv_inside(nsam, nsites);
                    if (rchrom >= 0)
                        config[chrom[rchrom].pop] += 1;
                }
            } else if (event == 'm') { /* migration event */
                migration_event(config, migration_matrix, npop, mig);
            } else { /* coalescent event */
                /* pick the two, c1, c2  */
                pick_2_chrom(cpop, config, &c1, &c2); /* c1 and c2 are chrom's to coalesce */
                dec = merge_chrom(nsam, nsites, c1, c2);
                config[cpop] -= dec;
            }
        }
    }
    /* adna */
    if (nsegs > 1) {
        ptree0 = seglst[0].ptree;
        for (seg = seglst[0].next, k = 1; k < nsegs; seg = seglst[seg].next, k++) {
            ptree1 = seglst[seg].ptree;
            for (i = nsamin; i < nsam; i++)
                (ptree1 + i)->time = (ptree0 + i)->time;
        }
    }

    *pnsegs = nsegs;
    free(config);
    free(size);
    free(alphag);
    free(tlast);
    for (i = 0; i < npop; i++)
        free(migration_matrix[i]);
    free(migration_matrix);
    return (seglst);
}

/*
 *	@brief Reserva la memoria inicial de las estructuras chrom, nnodes y seglst. 
 *
 *  Esta función se ejecuta solamente si el puntero correspondiente es NULL.
 * 
 *  - chrom es el arreglo de linajes o cromosomas, con nsam + CHROMOSOME_INC sitios. 
 * 
 *  - nnodes es el contador de nodos utilizado en el árbol de cada semento con seglimit
 *  lugares.
 * 
 *  - seglst es la lista de segmentos del locus con seglimit lugares.
 * 
 *	@param[in] nsam : Es el tamaño de la muestra. Se utiliza para alocar la memoria 
 *  inicial de chrom.
 *	@return	void : El resultado son las estructuras con un espacio inicial suficiente
 *  para almacenar los linajes en el caso de chrom, el contador de nodos del arbol 
 *  general en el caso de nnodes y la lista de segmentos (seglst).
 */
static void allocate_buffers(int nsam)
{
    /* Initialization */
    if (chrom == NULL) {
        maxchr = nsam + CHROMOSOME_INC;
        chrom = (struct chromo *)malloc((unsigned)(maxchr * sizeof(struct chromo)));
        if (chrom == NULL)
            perror("malloc error. segtre");
    }
    if (nnodes == NULL) {
        nnodes = (int *)malloc((unsigned)(seglimit * sizeof(int)));
        if (nnodes == NULL)
            perror("malloc error. simulate_coalescent");
    }
    if (seglst == NULL) {
        seglst = (struct segl *)malloc((unsigned)(seglimit * sizeof(struct segl)));
        if (seglst == NULL)
            perror("malloc error. simulate_coalescent.c 2");
    }
}

/*
 *	@brief Inicializa el cromosoma (linaje) de la posición ind en el arreglo chrom. 
 *
 *  Se reserva memoria para un solo segmento (nseg = 1) que cubre todo el locus que va desde 
 *  beg = 0 hasta end = nsites - 1. 
 * 
 *	@param[in] ind : Es el indice del linaje a inicializar en la estructura chrom
 *	@param[in] desc : Representa a qué nodo pertenece el segmento. Al principio 
 *  coincide con el índice del individuo 
 *	@param[in] pop : Indica a qué subpoblación pertenece el linaje
 *	@param[in] nsites : Es la cantidad de sitios que tiene el segmento
 *	@return	void : La posición "ind" de chrom inicializado
 */
static void init_chromo(int ind, int desc, int pop, int nsites)
{
    chrom[ind].nseg = 1;
    if (!(chrom[ind].pseg = (struct seg *)malloc((unsigned)sizeof(struct seg))))
        ERROR("calloc error. se1");
    (chrom[ind].pseg)->beg = 0;
    (chrom[ind].pseg)->end = nsites - 1;
    (chrom[ind].pseg)->desc = desc;
    chrom[ind].pop = pop;
}

// Termina el programa si hay linajes en más de una subpoblación sin migración entre ellas ni eventos demográficos pendientes (nunca coalescerían) 
static void check_isolated_populations(const int *config, int npop, double mig, const struct devent *nextevent)
{
    int i, j;

    if ((npop > 1) && (mig == 0.0) && (nextevent == NULL)) {
        i = 0;
        for (j = 0; j < npop; j++)
            if (config[j] > 0)
                i++;
        if (i > 1) {
            fprintf(stderr, " Infinite coalescent time. No migration.\n");
            exit(1);
        }
    }
}

// Evento de migración: elige un linaje al azar (proporcional a la tasa de migración de su subpoblación) y lo mueve a otra subpoblación
static void migration_event(int *config, double *const *migration_matrix, int npop, double mig)
{
    int i, migrant, source_pop;
    double x, sum;

    x = mig * ran1();
    sum = 0.0;
    for (i = 0; i < nchrom; i++) {
        sum += migration_matrix[chrom[i].pop][chrom[i].pop];
        if (x < sum)
            break;
    }
    migrant = i;
    x = ran1() * migration_matrix[chrom[i].pop][chrom[i].pop];
    sum = 0.0;
    for (i = 0; i < npop; i++) {
        if (i != chrom[migrant].pop) {
            sum += migration_matrix[chrom[migrant].pop][i];
            if (x < sum)
                break;
        }
    }
    source_pop = i;
    config[chrom[migrant].pop] -= 1;
    config[source_pop] += 1;
    chrom[migrant].pop = source_pop;
}

/*
 *	@brief Es la subrutina que ejecuta la recombinación : Elige un cromosoma y lo divide en dos partes. 
 *  Si el punto de corte es un lugar nuevo, un nuevo segmento se agrega a seglst y se asigna un nuevo arbol 
 *  para ese segmento.
 *	
 *	@param[in]	nsam : Es el tamaño de la muestra
 *	@return	ic : Es el índice del cromosoma modificado
 *
 */
static int recombination_event(int nsam)
{
    int ic, is, endic;

    /* First generate a random x-over spot, then locate it as to chrom and seg. */

    ic = find_crossover_chrom((long)(nchrom_links * ran1() + 1.), &is, &endic);
    crossover(nsam, ic, is);
    return (ic);
}

static int conv_left(int nsam)
{
    struct seg *pseg;
    int ic, is;
    double x, sum, len;

    while ((x = cleft * ran1()) == 0.0) ;
    sum = 0.0;
    ic = -1;
    while (sum < x) {
        sum += 1.0 - pow(pc, chrom_links(++ic));
    }
    pseg = chrom[ic].pseg;
    len = chrom_links(ic);
    is = pseg->beg + floor(1.0 + log(1.0 - (1.0 - pow(pc, len)) * ran1()) / lnpc) - 1;
    crossover(nsam, ic, is);
    return (ic);
}

/*
 *	@brief 
 *	
 *  @param[in] spot : Es el número de enlace aleatorio que representa la posición global a lo largo de 
 *                    todos los cromosomas donde se hace el corte.
 *  @param[in, out] pis : Puntero a la posición física donde se debe hacer el corte.
 *	@param[in, out] pend : Puntero al final del segmento del cromosoma seleccionado.
 *	@return	ic : Es el índice del cromosoma a ser modificado.
 *
 */
static int conv_inside(int nsam, int nsites)
{
    int len, ic, is, endic;

    /* First generate a random x-over spot, then locate it as to chrom and seg. */

    ic = find_crossover_chrom((long)(nchrom_links * ran1() + 1.), &is, &endic);
    crossover(nsam, ic, is);

    len = floor(1.0 + log(ran1()) / lnpc);
    if (is + len >= endic)
        return (ic);
    if (is + len < (chrom[nchrom - 1].pseg)->beg) {
        merge_chrom(nsam, nsites, ic, nchrom - 1);
        return (-1);
    }
    crossover(nsam, nchrom - 1, is + len);
    merge_chrom(nsam, nsites, ic, nchrom - 1);
    return (ic);
}

/*
 *	@brief Busca la posición del segmento a cortar (spot) y guarda la posición física del punto y el 
 *         final del segmento del cromosoma que se seleccionó.
 *	
 *  @param[in] spot : Es el número de enlace aleatorio que representa la posición global a lo largo de 
 *                    todos los cromosomas donde se hace el corte.
 *  @param[in, out] pis : Puntero a la posición física donde se debe hacer el corte.
 *	@param[in, out] pend : Puntero al final del segmento del cromosoma seleccionado.
 *	@return	ic : Es el índice del cromosoma a ser modificado.
 *
 */
static int find_crossover_chrom(long spot, int *pis, int *pend)
{
    struct seg *pseg;
    int el, lsgm1, ic;

    /* get chromosome # (ic)  */

    for (ic = 0; ic < nchrom; ic++) {
        lsgm1 = chrom[ic].nseg - 1;
        pseg = chrom[ic].pseg;
        el = ((pseg + lsgm1)->end) - (pseg->beg);
        if (spot <= el)
            break;
        spot -= el;
    }
    *pis = pseg->beg + spot - 1;
    *pend = (pseg + lsgm1)->end;
    return (ic);
}

/*
 *	@brief Es la función que efectúa el corte del cromosomas, se encarga de modificar las estructuras
 *         involucradas.
 * 
 *  
 *	
 *  @param[in] nsam : Es el tamaño de la muestra.
 *  @param[in] ic : Es el índice del cromosoma que se debe modificar.
 *	@param[in] is : Es la posición física donde se debe hacer el corte.
 *	@return	ic : Es el índice del cromosoma modificado.
 *
 */
static int crossover(int nsam, int ic, int is)
{
    struct seg *pseg, *pseg2;
    int i, lsg, lsgm1, newsg, jseg, k, in;
    double len;

    pseg = chrom[ic].pseg;
    lsg = chrom[ic].nseg;
    len = (pseg + lsg - 1)->end - pseg->beg;
    cleft -= 1 - pow(pc, len);
    /* get seg # (jseg)  */

    for (jseg = 0; is >= (pseg + jseg)->end; jseg++);

    if (is >= (pseg + jseg)->beg)
        in = 1;
    else
        in = 0;
    newsg = lsg - jseg;

    /* copy last part of chrom to nchrom  */

    nchrom++;
    if ( (unsigned int)nchrom >= maxchr) {
        maxchr += CHROMOSOME_INC;
        chrom = (struct chromo *)realloc(chrom, (unsigned)(maxchr * sizeof(struct chromo)));
        if (chrom == NULL)
            perror("malloc error. segtre2");
    }
    if (!(pseg2 = chrom[nchrom - 1].pseg =
              (struct seg *)calloc((unsigned)newsg, sizeof(struct seg))))
        ERROR(" alloc error. re1");
    chrom[nchrom - 1].nseg = newsg;
    chrom[nchrom - 1].pop = chrom[ic].pop;
    pseg2->end = (pseg + jseg)->end;
    if (in) {
        pseg2->beg = is + 1;
        (pseg + jseg)->end = is;
    } else
        pseg2->beg = (pseg + jseg)->beg;
    pseg2->desc = (pseg + jseg)->desc;
    for (k = 1; k < newsg; k++) {
        (pseg2 + k)->beg = (pseg + jseg + k)->beg;
        (pseg2 + k)->end = (pseg + jseg + k)->end;
        (pseg2 + k)->desc = (pseg + jseg + k)->desc;
    }

    lsg = chrom[ic].nseg = lsg - newsg + in;
    lsgm1 = lsg - 1;
    nchrom_links -= pseg2->beg - (pseg + lsgm1)->end;
    len = (pseg + lsgm1)->end - (pseg->beg);
    cleft += 1.0 - pow(pc, len);
    len = (pseg2 + newsg - 1)->end - pseg2->beg;
    cleft += 1.0 - pow(pc, len);
    if (!(chrom[ic].pseg =
              (struct seg *)realloc(chrom[ic].pseg, (unsigned)(lsg * sizeof(struct seg)))))
        perror(" realloc error. re2");
    if (in) {
        begs = pseg2->beg;
        for (i = 0, k = 0; (k < nsegs - 1) && (begs > seglst[seglst[i].next].beg - 1);
             i = seglst[i].next, k++)
            ;
        if (begs != seglst[i].beg) /* new tree  */
            insert_new_tree(i, nsam);
    }
    return (ic);
}

/*
 *	@brief Agrega a seglst un segmento nuevo que empieza en begs ubicado después del segmento i 
 *  con una copia del árbol del segmento i.
 * 
 *  
 *  @param[in] i : .  
 *  @param[in] nsam : Es el tamaño de la muestra.
 *	@return	void
 *
 */
static void insert_new_tree(int i, int nsam)
{
    int k;

    if ( (unsigned int)nsegs >= seglimit) {
        seglimit += SEGINC;
        nnodes = (int *)realloc(nnodes, (unsigned)(sizeof(int) * seglimit));
        if (nnodes == NULL)
            perror("realloc error. 1. simulate_coalescent.c");
        seglst = (struct segl *)realloc(seglst, (unsigned)(sizeof(struct segl) * seglimit));
        if (seglst == NULL)
            perror("realloc error. 2. simulate_coalescent.c");
        /*  printf("seglimit: %d\n",seglimit);  */
    }
    seglst[nsegs].next = seglst[i].next;
    seglst[i].next = nsegs;
    seglst[nsegs].beg = begs;
    if (!(seglst[nsegs].ptree =
              (struct node *)calloc((unsigned)(2 * nsam), sizeof(struct node))))
        perror("calloc error. re3.");
    nnodes[nsegs] = nnodes[i];
    ptree1 = seglst[i].ptree;
    ptree2 = seglst[nsegs].ptree;
    nsegs++;
    for (k = 0; k <= nnodes[i]; k++) {
        (ptree2 + k)->abv = (ptree1 + k)->abv;
        (ptree2 + k)->time = (ptree1 + k)->time;
    }
}

//   common ancestor subroutine
//   Pick two chromosomes and merge them. Update trees if necessary.

/*
 *	@brief Subrutina que ejecuta el evento de coalescencia o ancestro común. Actualiza arboles si es 
 *  necesario.
 *  
 *  Recorre todos los segmentos
 * 
 *	@param[in] nsam : Tamaño de la muestra
 *	@param[in] nsites : Cantidad de sitios
 *	@param[in] c1 : Indice del cromosoma c1
 *	@param[in] c2 : Indice del cromosoma c2
 *	@return	int : 2 linajes menos 
 *                1 linaje menos
 *
 */
static int merge_chrom(int nsam, int nsites, int c1, int c2)
{
    int yes1, yes2, seg1, seg2, seg;
    int tseg, start, end, desc, k;
    struct seg *pseg;
    struct node *ptree;

    seg1 = 0;   // Punteros que guardan el segmento de cada cromosoma, se van   
    seg2 = 0;   // incrementando de a uno en la función chromo_has_segment()

    if (nsegs > buf_cap){
        free(buf) ;
        buf_cap = nsegs + SEGINC ; 
        if (!(buf = (struct seg *)malloc((unsigned)buf_cap * sizeof(struct seg))))
            perror("Error malloc") ;
    }
    pseg = buf ;

    tseg = -1;

    for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++) {  // Recorre todos los segmentos de seglst usando start que
        start = seglst[seg].beg;                                    // apunta al inicio de cada segmento de seglst
        yes1 = chromo_has_segment(start, c1, &seg1);
        yes2 = chromo_has_segment(start, c2, &seg2);
        if (yes1 || yes2) {
            tseg++;
            (pseg + tseg)->beg = seglst[seg].beg;
            end = (k < nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : nsites - 1);
            (pseg + tseg)->end = end;

            if (yes1 && yes2) {
                nnodes[seg]++;
                if (nnodes[seg] >= (2 * nsam - 2))  // Llegamos a la raíz?
                    tseg--;
                else
                    (pseg + tseg)->desc = nnodes[seg];
                ptree = seglst[seg].ptree;
                desc = (chrom[c1].pseg + seg1)->desc;
                (ptree + desc)->abv = nnodes[seg];
                desc = (chrom[c2].pseg + seg2)->desc;
                (ptree + desc)->abv = nnodes[seg];
                (ptree + nnodes[seg])->time = t;

            } else {
                (pseg + tseg)->desc =
                    (yes1 ? (chrom[c1].pseg + seg1)->desc : (chrom[c2].pseg + seg2)->desc);
            }
        }
    }
    remove_chrom_links(c1);
    free(chrom[c1].pseg);
    if (tseg < 0) {
        move_last_chrom_to(c1);
        if (c2 == nchrom - 1)
            c2 = c1;
        nchrom--;
    } else {
        struct seg *nuevo = (struct seg *)malloc((unsigned)(tseg + 1) * sizeof(struct seg)) ;
        if (nuevo == NULL)
            perror("Error malloc") ;
        memcpy(nuevo, pseg, (unsigned)(tseg + 1) * sizeof(struct seg));
        chrom[c1].pseg = nuevo ;
        chrom[c1].nseg = tseg + 1 ;
        add_chrom_links(c1);
    }
    remove_chrom_links(c2);
    free(chrom[c2].pseg);
    move_last_chrom_to(c2);
    nchrom--;
    if (tseg < 0)
        return (2); /* decrease of nchrom is two */
    else
        return (1);
}

// Descuenta del total los enlaces del cromosoma c y su aporte a cleft 
static void remove_chrom_links(int c)
{
    nchrom_links -= chrom_links(c);
    cleft -= 1.0 - pow(pc, (double)chrom_links(c));
}

// Suma al total los enlaces del cromosoma c y su aporte a cleft 
static void add_chrom_links(int c)
{
    nchrom_links += chrom_links(c);
    cleft += 1.0 - pow(pc, (double)chrom_links(c));
}

/*
 *	@brief Desplaza el último cromosoma del arreglo hacia la posición "c"
 *  
 *	@param[in] c : Posición a la que se quiere desplazar el último cromosoma
 *	@return	void : Como resultado, el último cromosoma se desplaza a la posición "c" 
 */
static void move_last_chrom_to(int c)
{
    chrom[c].pseg = chrom[nchrom - 1].pseg;
    chrom[c].nseg = chrom[nchrom - 1].nseg;
    chrom[c].pop = chrom[nchrom - 1].pop;
}

// chromo_has_segment: Does chromosome c contain the segment on seglst which starts at
//            start? *psg is the segment of chrom[c] at which one is to begin looking.  

/*
 *	@brief Responde a la pregunta ¿el cromosoma c contiene el segmento de seglst que empieza en start?
 *  *psg es el segmento de chrom[c] que está siendo buscado.
 *  
 *  El bucle principal verifica si el sitio start está acotado por end. Si esto es así, entonces el 
 *  cromosoma tiene el sitio start y retorna verdadero, caso contrario, no lo tiene y retorna falso.
 *  
 *	@param[in] start : Posición que queremos ver si el cromosoma solapa
 *  @param[in] c : Es el índice del cromosoma que se quiere verificar 
 *  @param[in] psg : Puntero auxiliar para recorrer los segmentos de manera eficiente
 *	@return	int : 1 para indicar encontrado
 *                0 para indicar no encontrado
 */
static int chromo_has_segment(int start, int c, int *psg)
{
    int ns;
    struct seg *pseg;

    ns = chrom[c].nseg;
    pseg = chrom[c].pseg;

    /*  changed order of test conditions in following line on 6 Dec 2004 */
    for (; ((*psg) < ns) && ((pseg + (*psg))->beg <= start); ++(*psg))
        if ((pseg + (*psg))->end >= start)
            return (1);
    return (0);
}

static void pick_2_chrom(int pop, const int config[], int *pc1, int *pc2)
{
    int c1, c2, cs, cb;

    pick_2_random_numbers(config[pop], &c1, &c2);
    cs = (c1 > c2) ? c2 : c1;
    cb = (c1 > c2) ? c1 : c2;
    *pc1 = next_chrom_in_pop(pop, cs, 0, 0);
    *pc2 = next_chrom_in_pop(pop, cb, *pc1 + 1, cs + 1);
}

// Avanza desde el cromosoma i (que es el número count de la subpoblación pop) hasta el cromosoma target de esa subpoblación, y devuelve su índice
static int next_chrom_in_pop(int pop, int target, int i, int count)
{
    for (;;) {
        while (chrom[i].pop != pop)
            i++;
        if (count == target)
            return (i);
        count++;
        i++;
    }
}

/* pick_2_random_numbers()  */

static int pick_2_random_numbers(int n, int *i, int *j)
{
    *i = n * ran1();
    while ((*j = n * ran1()) == *i);

    return (0);
}

/****  chrom_links(c): returns the number of chrom_links between beginning and end of chrom **/
// Retorna el número de enlaces del cromosoma c que van desde el inicio y el final del cromosoma
// (chrom[c].pseg + ns) apunta al último seg y (chrom[c].pseg) al primero
static int chrom_links(int c)
{
    int ns;

    ns = chrom[c].nseg - 1;

    return ((chrom[c].pseg + ns)->end - (chrom[c].pseg)->beg);
}