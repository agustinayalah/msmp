/***** ms.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.
        usage: ms nsam howmany  -t  theta  [options]
                or
               ms nsam howmany -s segsites  [options]

           nsam is the number of gametes per sample.
           howmany is the number of samples to produce.
           With -t the numbers of segregating sites will randomly vary
                from one sample to the next.
           with -s segsites,  the number of segregating sites will be
                segsites in each sample.

           Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


*	  Arguments of the options are explained here:

           npop:  Number of subpopulations which make up the total population
           ni:  the sample size from the i th subpopulation (all must be
                specified.) The output will have the gametes in order such that
                the first n1 gametes are from the first island, the next n2 are
                from the second island, etc.
           nsites: number of sites between which recombination can occur.
           theta: 4No times the neutral mutation rate
           rho: recombination rate between ends of segment times 4No
           f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
           track_len:  mean length of conversion track in units of sites.  The
                       total number of sites is nsites, specified with the -r option.
           mig_rate: migration rate: the fraction of each subpop made up of
                 migrants times 4No.
           howmany: howmany samples to generate.

        Note:  In the above definition, No is the total diploid population if
                npop is one, otherwise, No is the diploid population size of each
                subpopulation.
        A seed file called "seedms" will be created  if it doesn't exist. The
                seed(s) in this file will be modified by the program.
                So subsequent runs
                will produce new output.  The initial contents of seedms will be
                printed on the second line of the output.
        Output consists of one line with the command line arguments and one
                line with the seed(s).
                The samples appear sequentially following that line.
                Each sample begins with "//", then the number of segregating sites, the positions
                of the segregating sites (on a scale of 0.0 - 1.0). On the following
                lines are the sampled gametes, with mutants alleles represented as
                ones and ancestral alleles as zeros.
        To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
                or:  cc -o ms ms.c streec.c rand2.c -lm
         (Of course, gcc would be used instead of cc on some machines.  And -O3 or
                some other optimization switches might be usefully employed with some
                compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
        Modifications to add gene conversion 5 Mar 2001
        Added demographic options -d  13 Mar 2001
        Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
        Changed cleftr() to check for zero rand() .13 June 2001
        Move seed stuff to subroutine seedit()  11 July 2001
        Modified streec.c to handle zero length demographic intervals 9 Aug 2001
        Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
        Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
        Changed many command line options, allowing arbitrary migration matrix, and subpopulation
           sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
           split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
        !!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used.
9 Dec 2003 Fixed bug which resulted in incorrect results for the case where rho = 0.0 and gene
conversion rate > 0.0. This case was not handled correctly in early versions of the program. 
5 Apr 2004.  (Thanks to Vincent Plagnol for pointing out this problem.) Fixed bug in print_tree().  Earlier
versions may have produced garbage when the -T option was used. 
1 Jul 2004. Fixed bug in -e. options that caused problems with -f option  
13 Aug 2004. Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 
6 Nov. 2004 Added -F option: -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  
11 Nov. 2004. Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks to Melissa Jane
Todd-Hubisz for finding and reporting this bug.) Added -seeds option 
4 Nov 2006 Added "tbs" arguments feature 
4 Nov 2006 Added -L option.  
10 May 2007 Changed -ej option to set Mki = 0 pastward of the ej event.  See msdoc.pdf.  
May 19 2007. fixed bug with memory allocation when using -I option. This caused problems expecially on Windows machines.  
Thanks to several people, including Vitor Sousa and Stephane De Mita for help on this one. 
Oct. 17, 2007. Modified pick_branch() and pick_branch_freq() to eliminate rare occurrence of fixed mutations 
Thanks to J. Liechty and K. Thornton.
4 Feb 2010. Added -p n  switch to allow position output to be higher precision.  
10 Nov 2012. Changed spot from int to long in the function re().  
29 July 2013. ( Thanks to Yuri D'Elia for this suggestion.) Changed function definitions, and misc other things to comply with c99
     requirements.  4 Mar 2014.
***************************************************************************/

#include "ms.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SITESINC 10       
#define TBS_PARAMS_LENGTH 30
#define PARAM_FILE_ARG_LENGTH 20
#define MIN_ARGS 4                  // Representa la cantidad mínima de argumentos
#define FIRST_OPTION_ARG 3          // Representa la posición de la primera opción en los parámetros
#define OUTPUT_PRECISION 4          // Representa la cantidad de decimales para las posiciones
#define MIN_NSITES 2
#define MIN_MFREQ 2

unsigned maxsites = SITESINC;

double *posit;
double *agevec;
double alleleage;
double segfac;
int count, ntbs, nseeds;
struct params pars;

// Prototypes

static int generate_sample(char **list, double *pprobss, double *ptmrca, double *pttot);
static void count_nodes(struct node *ptree, int nsam);
static void bigger_list(int nsam, char **list);
static char **cmatrix(int nsam, int len);
static void locate(int n, double beg, double len, double *ptr);

static void getpars(int argc, char *argv[], int *phowmany);
static void argcheck(int arg, int argc, char *const argv[]);
static void usage(void);
static void insert_event_sorted(struct devent *pt);
static char **allocate_buffers(void);
static void print_sample(FILE *pf, char *const *list, int segsites, double probss);
static int segment_end(const struct segl *seglst, int seg, int k, int nsegs, int nsites);
static double segment_tree_time(const struct node *ptree, int nsam, int mfreq);
static void print_segment_trees(const struct segl *seglst, int nsegs);
static void compute_tree_times(const struct segl *seglst, int nsegs, double *ptmrca, double *pttot);
static int mutate_with_theta(const struct segl *seglst, int nsegs, char **list);
static int mutate_fixed_segsites(const struct segl *seglst, int nsegs, char **list, double *pprobss);
static void ensure_sites_capacity(int needed, int nsam, char **list);
static void set_default_params(void);
static int parse_island_model(int arg, int argc, char *const argv[], int *pnpop);
static int parse_demographic_event(int arg, int argc, char *const argv[]);
static void calc_matrix_diagonal(double **mat, int npop);
static double branch_length(const struct node *ptree, int i); 
static double random_age_in_branch(const struct node *ptree, int i);  
static void free_event_list(struct devent *pt, int npop);
static void make_gametes(int nsam, int mfreq, const struct node *ptree, double tt, int newsites, int ns, char **list);
static double tree_time(const struct node *ptree, int nsam);
static double tree_time_mfreq(const struct node *ptree, int nsam, int mfreq);
static void print_tree(const struct node *ptree, int nsam);
static void print_parenthesis(const struct node *ptree, const int *left_child, const int *right_child, int noden);
static int pick_branch(int nsam, const struct node *ptree, double tt);
static int pick_branch_freq(int nsam, int mfreq, const struct node *ptree, double tt);
static int is_descendant(const struct node *ptree, int tip, int node);
static void ordran(int n, double pbuf[]);
static void multinomial(int n, int nclass, const double p[], int rv[]);
static void order(int n, double pbuf[]);
static void ranvec(int n, double pbuf[]);
static int poisso(double u);
static double gasdev(double m, double v);

int main(int argc, char *argv[])
{
    int i, k, howmany, segsites;
    char **list, **tbsparamstrs;
    FILE *pf;
    double probss, tmrca, ttot;

    /* these next few lines are for reading in parameters from a file (for each sample) */
    ntbs = 0;
    tbsparamstrs = (char **)malloc(argc * sizeof(char *));

    for (i = 0; i < argc; i++)
        printf("%s ", argv[i]);
    for (i = 0; i < argc; i++)
        tbsparamstrs[i] = (char *)malloc(TBS_PARAMS_LENGTH * sizeof(char));
    for (i = 1; i < argc; i++)
        if (strcmp(argv[i], "tbs") == 0)
            argv[i] = tbsparamstrs[ntbs++];

    count = 0;

    if (ntbs > 0)
        for (k = 0; k < ntbs; k++)
            scanf(" %s", tbsparamstrs[k]);

    getpars(argc, argv, &howmany);

    /* results are stored in global variable, pars */

    if (!pars.commandlineseedflag)
        seedit("s");               // Si no hay parametros, seteamos el seed en "s"

    pf = stdout;

    list = allocate_buffers();

    while (howmany - count++) {
        if ((ntbs > 0) && (count > 1)) {
            for (k = 0; k < ntbs; k++) {
                if (scanf(" %s", tbsparamstrs[k]) == EOF) {
                    if (!pars.commandlineseedflag)
                        seedit("end");
                    exit(0);
                }
            }
            getpars(argc, argv, &howmany);
        }

        fprintf(pf, "\n//");
        if (ntbs > 0) {
            for (k = 0; k < ntbs; k++)
                printf("\t%s", tbsparamstrs[k]);
        }
        printf("\n");

        segsites = generate_sample(list, &probss, &tmrca, &ttot);

        if (pars.mp.timeflag)                                      
            fprintf(pf, "time:\t%lf\t%lf\n", tmrca, ttot);

        print_sample(pf, list, segsites, probss);
    }

    if (!pars.commandlineseedflag)
        seedit("end");
}

/*
 *	@brief Genera los arboles de coalescencia y aplica las mutaciones.
 *
 *	Invoca a simulate_coalescent() que guarda en seglst el arbol de coalescencia. Luego en función de las
 *opciones ingresadas, muestra el arbol en formato Newick o muestra el tiempo que le tomó a cada
 *rama hasta coalescer. Por último aplica las mutaciones sobre las ramas de manera aleatoria con
 *make_gametes().
 *
 *	@param[out]	list : Lista de segmentos de ADN de los individuos a completar
 *	@param[out]	pprobss : La probailidad de que bajo una tasa de mutación theeta se den n
 *segsites fijos
 *	@param[out]	ptmrca : Tiempo en llegar al most recent common ancester
 *	@param[out]	pttot : Tiempo total para la coalescencia del arbol
 *	@return	ns : Cantidad de sitios segregantes
 *
 */
static int generate_sample(char **list, double *pprobss, double *ptmrca, double *pttot)
{
    int nsegs, i, ns = 0;
    struct segl *seglst;

    seglst = simulate_coalescent(&(pars.cp), &nsegs);

    if (pars.mp.treeflag)
        print_segment_trees(seglst, nsegs);

    if (pars.mp.timeflag)
        compute_tree_times(seglst, nsegs, ptmrca, pttot);

    if ((pars.mp.segsitesin == 0) && (pars.mp.theta > 0.0)) // Si se pasa theta sin segsites fijos
        ns = mutate_with_theta(seglst, nsegs, list);
    else if (pars.mp.segsitesin > 0) // Si se pasa con segsites fijos
        ns = mutate_fixed_segsites(seglst, nsegs, list, pprobss);

    for (i = 0; i < pars.cp.nsam; i++)
        list[i][ns] = '\0';
    return (ns);
}

// Imprime el árbol, se activa con la opción -T 
static void print_segment_trees(const struct segl *seglst, int nsegs)
{
    int seg, k, start, end, len;
    int nsites = pars.cp.nsites, nsam = pars.cp.nsam;

    for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++) {
        if ((pars.cp.r > 0.0) || (pars.cp.f > 0.0)) {
            end = segment_end(seglst, seg, k, nsegs, nsites);
            start = seglst[seg].beg;
            len = end - start + 1;
            fprintf(stdout, "[%d]", len);
        }
        print_tree(seglst[seg].ptree, nsam);
        if ((pars.mp.segsitesin == 0) && (pars.mp.theta == 0.0) && (pars.mp.timeflag == 0))
            free(seglst[seg].ptree);
    }
}

// Usando la opción -L calcula el tiempo al MRCA y la longitud total del árbol 
static void compute_tree_times(const struct segl *seglst, int nsegs, double *ptmrca, double *pttot)
{
    int seg, k, start, end, len;
    int nsites = pars.cp.nsites, nsam = pars.cp.nsam, mfreq = pars.mp.mfreq;
    double tseg, tt = 0.0;

    for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++) {
        if (mfreq > 1)
            count_nodes(seglst[seg].ptree, nsam);
        end = segment_end(seglst, seg, k, nsegs, nsites);
        start = seglst[seg].beg;
        if ((nsegs == 1) || ((start <= nsites / 2) && (end >= nsites / 2)))
            *ptmrca = (seglst[seg].ptree + 2 * nsam - 2)->time;
        len = end - start + 1;
        tseg = len / (double)nsites;
        tt += segment_tree_time(seglst[seg].ptree, nsam, mfreq) * tseg;
        if ((pars.mp.segsitesin == 0) && (pars.mp.theta == 0.0))
            free(seglst[seg].ptree);
    }
    *pttot = tt;
}

// Aplica mutaciones con tasa theta (opción -t) y devuelve la cantidad de sitios segregantes 
static int mutate_with_theta(const struct segl *seglst, int nsegs, char **list)
{
    int seg, k, start, end, len, segsit, ns = 0;
    int nsites = pars.cp.nsites, nsam = pars.cp.nsam, mfreq = pars.mp.mfreq;
    double nsinv = 1. / nsites, theta = pars.mp.theta, tseg, tt;

    for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++) {
        if (mfreq > 1)
            count_nodes(seglst[seg].ptree, nsam);
        end = segment_end(seglst, seg, k, nsegs, nsites);
        start = seglst[seg].beg;
        len = end - start + 1;
        tseg = len * (theta / nsites);
        tt = segment_tree_time(seglst[seg].ptree, nsam, mfreq);
        segsit = poisso(tseg * tt);
        ensure_sites_capacity(segsit + ns, nsam, list);
        make_gametes(nsam, mfreq, seglst[seg].ptree, tt, segsit, ns, list);
        free(seglst[seg].ptree);
        locate(segsit, start * nsinv, len * nsinv, posit + ns);
        ns += segsit;
    }
    return (ns);
}

// Aplica una cantidad fija de sitios segregantes entre los segmentos (opción -s) y retorna los segsites
static int mutate_fixed_segsites(const struct segl *seglst, int nsegs, char **list, double *pprobss)
{
    int seg, k, start, end, len, ns;
    int nsites = pars.cp.nsites, nsam = pars.cp.nsam, mfreq = pars.mp.mfreq;
    int segsitesin = pars.mp.segsitesin;
    double nsinv = 1. / nsites, theta = pars.mp.theta, tseg, tt, es;
    double *pk;
    int *ss;

    pk = (double *)malloc((unsigned)(nsegs * sizeof(double)));
    ss = (int *)malloc((unsigned)(nsegs * sizeof(int)));
    if ((pk == NULL) || (ss == NULL))
        perror("malloc error. generate_sample.2");

    tt = 0.0;
    for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++) {
        if (mfreq > 1)
            count_nodes(seglst[seg].ptree, nsam);
        end = segment_end(seglst, seg, k, nsegs, nsites);
        start = seglst[seg].beg;
        len = end - start + 1;
        tseg = len / (double)nsites;
        pk[k] = segment_tree_time(seglst[seg].ptree, nsam, mfreq) * tseg;
        tt += pk[k];
    }
    if (theta > 0.0) {
        es = theta * tt;
        *pprobss = exp(-es) * pow(es, (double)segsitesin) / segfac;
    }
    if (tt > 0.0) {
        for (k = 0; k < nsegs; k++)
            pk[k] /= tt;
        multinomial(segsitesin, nsegs, pk, ss);
    } else
        for (k = 0; k < nsegs; k++)
            ss[k] = 0;
    ns = 0;
    for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++) {
        end = segment_end(seglst, seg, k, nsegs, nsites);
        start = seglst[seg].beg;
        len = end - start + 1;
        tseg = len / (double)nsites;
        make_gametes(nsam, mfreq, seglst[seg].ptree, tt * pk[k] / tseg, ss[k], ns, list);

        free(seglst[seg].ptree);
        locate(ss[k], start * nsinv, len * nsinv, posit + ns);
        ns += ss[k];
    }
    free(pk);
    free(ss);
    return (ns);
}

static int segment_end(const struct segl *seglst, int seg, int k, int nsegs, int nsites)
{
    return (k < nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : nsites - 1);
}

static double segment_tree_time(const struct node *ptree, int nsam, int mfreq)
{
    if (mfreq == 1)
        return tree_time(ptree, nsam);
    else
        return tree_time_mfreq(ptree, nsam, mfreq);
}

static void ensure_sites_capacity(int needed, int nsam, char **list)
{
    if ( (unsigned int)needed >= maxsites) {
        maxsites = needed + SITESINC;
        posit = (double *)realloc(posit, maxsites * sizeof(double));
        agevec = (double *)realloc(agevec, maxsites * sizeof(double));
        bigger_list(nsam, list);
    }
}

//  
/*
 *	@brief Reserva espacio para la matriz de gametos y los vectores de posiciones y edades
 *
 *	Si no se fija la cantidad de sitios segregantes con -s, entonces reserva con un máximo predefinido 
 *  maxsegsites. Caso contrario, entonces se reserva un espacio fijo de tamaño segsites, adicionalmente
 *  se calcula el factor de segregación "segfac" que se utiliza para calcular la probabilidad de obtener
 *  segsites cantidad de mutaciones para el tamaño de muestra elegido.
 *
 *	@return	list : Matríz de gametos inicializada con tamaño maxsites o segsites
 *
 */
static char **allocate_buffers(void)
{
    int i;
    char **list;

    if (pars.mp.segsitesin == 0) {                  // Tasa de mutación con -t theta
        list = cmatrix(pars.cp.nsam, maxsites + 1);
        posit = (double *)malloc((unsigned)(maxsites * sizeof(double)));    // posit y agevec son los metadatos de los segmentos como la
        agevec = (double *)malloc((unsigned)(maxsites * sizeof(double)));   // posición donde ocurrió la mutación (rango de 0 a 1) o la edad
                                                                            // en la que ocurrió la mutación (hacia atrás) 
    } else {                                        // Sitios segregantes fijos -s segsites
        list = cmatrix(pars.cp.nsam, pars.mp.segsitesin + 1);
        posit = (double *)malloc((unsigned)(pars.mp.segsitesin * sizeof(double)));
        agevec = (double *)malloc((unsigned)(pars.mp.segsitesin * sizeof(double)));
        if (pars.mp.theta > 0.0) {
            segfac = 1.0;                           // segfac = segsitesin!
            for (i = pars.mp.segsitesin; i > 1; i--)
                segfac *= i;
        }
    }
    return (list);
}

// Imprime los resultados de una réplica como los segsites, posiciones, edades y gametos 
static void print_sample(FILE *pf, char *const *list, int segsites, double probss)
{
    int i, k, afreq;

    if ((segsites > 0) || (pars.mp.theta > 0.0)) {
        if ((pars.mp.segsitesin > 0) && (pars.mp.theta > 0.0))
            fprintf(pf, "prob: %g\n", probss);
        fprintf(pf, "segsites: %d\n", segsites);
        if (segsites > 0)
            fprintf(pf, "positions: ");
        for (i = 0; i < segsites; i++)
            fprintf(pf, "%6.*lf ", pars.output_precision, posit[i]);
        fprintf(pf, "\n");
        if ((segsites > 0) && pars.mp.ageflag) {
            fprintf(pf, "allele ages: ");
            for (i = 0; i < segsites; i++)
                fprintf(pf, "%6.*lf ", pars.output_precision, agevec[i]);
            fprintf(pf, "\n");
            fprintf(pf, "allele freqs: ");
            for (i = 0; i < segsites; i++) {
                for (k = afreq = 0; k < pars.cp.nsam; k++)
                    afreq += ((list[k][i] == '1') ? 1 : 0);
                fprintf(pf, "%d ", afreq);
            }
            fprintf(pf, "\n");
        }
        if (segsites > 0)
            for (i = 0; i < pars.cp.nsam; i++) {
                fprintf(pf, "%s\n", list[i]);
            }
    }
}

static void count_nodes(struct node *ptree, int nsam)
{
    int i;
    for (i = 0; i < nsam; i++)
        (ptree + i)->ndes = 1;
    for (i = nsam; i < 2 * nsam - 1; i++)
        (ptree + i)->ndes = 0;
    for (i = 0; i < 2 * nsam - 2; i++)
        (ptree + ((ptree + i)->abv))->ndes += (ptree + i)->ndes;
}

static void bigger_list(int nsam, char **list)
{
    int i;
    /*  fprintf(stderr,"maxsites: %d\n",maxsites);  */
    for (i = 0; i < nsam; i++) {
        list[i] = (char *)realloc(list[i], maxsites * sizeof(char));
        if (list[i] == NULL)
            perror("realloc error. bigger");
    }
}

/*
 *	@brief Reserva memoria de manera dinámica para la matríz que guarda las secuencias de ADN
 *  
 *  Asigna espacio para un arreglo de punteros donde la cantidad de filas son la cantidad de 
 *  individuos y la cantidad de columnas la cantidad de sitios segregantes o mutaciones que
 *  surgen para cada secuencia.
 * 
 *	@param[in]	nsam : Tamaño de la muestra
 *	@param[in]	len : Longitud de cada cadena. Puede ser prefijado con nsegsites (opción -s)
 *                    o puede usar una cantidad maxima fijada de antemano (maxsites)
 *  @return matrix : Puntero a la matríz de caracteres
 */
static char **cmatrix(int nsam, int len)
{
    int i;
    char **matrix;

    if (!(matrix = (char **)malloc((unsigned)nsam * sizeof(char *))))
        perror("alloc error in cmatrix");
    for (i = 0; i < nsam; i++) {
        if (!(matrix[i] = (char *)malloc((unsigned)len * sizeof(char))))
            perror("alloc error in cmatric. 2");
    }
    return (matrix);
}

static void locate(int n, double beg, double len, double *ptr)
{
    int i;
    ordran(n, ptr);
    for (i = 0; i < n; i++)
        ptr[i] = beg + ptr[i] * len;
}

/*
 *	@brief Obtiene los parámetros de la simulación, los valida y procesa.
 *
 *  Analiza los parámetros validando la cantidad de argumentos que requiere el programa. Los guarda 
 *  en la estructura global de parámetros "pars" y valida que sean correctos. Soporta el pasaje de
 *  parámetros en diferido o mediante un archivo con la opción -f.
 * 
 *	@param[in]	argc : Cantidad de parámetros pasados al ejecutar el programa
 *	@param[in]	argv : Lista de parámetros pasados al ejecutar el programa
 *	@param[in]	phowmany : Cantidad de repeticiones sobre un mismo conjunto de muestras a simular
 *
 */
static void getpars(int argc, char *argv[], int *phowmany)
{
    int arg, i, j, sum, pop, argstart, npop, pop2;
    double mij, psize, palpha;
    FILE *pf;

    if (count == 0) {
        if (argc < MIN_ARGS) {
            fprintf(stderr, "Too few command line arguments\n");
            usage();
        }
        /* adna */

        pars.cp.nsamin = atoi(argv[1]);
        pars.cp.nsam = pars.cp.nsamin;

        if (pars.cp.nsam <= 0) {
            fprintf(stderr, "First argument error. nsam <= 0. \n");
            usage();
        }

        *phowmany = atoi(argv[2]);

        if (*phowmany <= 0) {
            fprintf(stderr, "Second argument error. howmany <= 0. \n");
            usage();
        }

        set_default_params();
        npop = pars.cp.npop;
    } else {
        npop = pars.cp.npop;
        free_event_list(pars.cp.deventlist, npop);
    }

    pars.cp.deventlist = NULL;

    arg = FIRST_OPTION_ARG ;

    while (arg < argc) {
        if (argv[arg][0] != '-') {
            fprintf(stderr, " argument should be -%s ?\n", argv[arg]);
            usage();
        }
        switch (argv[arg][1]) {
        case 'f':
            if (ntbs > 0) {
                fprintf(stderr, " can't use tbs args and -f option.\n");
                exit(1);
            }
            arg++;
            argcheck(arg, argc, argv);
            pf = fopen(argv[arg], "r");
            if (pf == NULL) {
                fprintf(stderr, " no parameter file %s\n", argv[arg]);
                exit(0);
            }
            arg++;
            argc++;
            argv = (char **)malloc((unsigned)(argc + 1) * sizeof(char *));
            argv[arg] = (char *)malloc((unsigned)(PARAM_FILE_ARG_LENGTH * sizeof(char)));
            argstart = arg;
            while (fscanf(pf, " %s", argv[arg]) != EOF) {
                arg++;
                argc++;
                argv = (char **)realloc(argv, (unsigned)argc * sizeof(char *));
                argv[arg] = (char *)malloc((unsigned)(PARAM_FILE_ARG_LENGTH * sizeof(char)));
            }
            fclose(pf);
            argc--;
            arg = argstart;
            break;
        case 'r':
            arg++;
            argcheck(arg, argc, argv);
            pars.cp.r = atof(argv[arg++]);
            argcheck(arg, argc, argv);
            pars.cp.nsites = atoi(argv[arg++]);
            if (pars.cp.nsites < MIN_NSITES) {
                fprintf(stderr, "with -r option must specify both rec_rate and nsites>1\n");
                usage();
            }
            break;
        case 'p':
            arg++;
            argcheck(arg, argc, argv);
            pars.output_precision = atoi(argv[arg++]);
            break;
        case 'c':
            arg++;
            argcheck(arg, argc, argv);
            pars.cp.f = atof(argv[arg++]);
            argcheck(arg, argc, argv);
            pars.cp.track_len = atof(argv[arg++]);
            if (pars.cp.track_len < 1.) {
                fprintf(stderr, "with -c option must specify both f and track_len>0\n");
                usage();
            }
            break;
        case 't':
            arg++;
            argcheck(arg, argc, argv);
            pars.mp.theta = atof(argv[arg++]);
            break;
        case 's':
            arg++;
            argcheck(arg, argc, argv);
            if (argv[arg - 1][2] == 'e') { /* command line seeds */
                pars.commandlineseedflag = 1;
                if (count == 0)
                    nseeds = commandlineseed(argv + arg);
                arg += nseeds;
            } else {
                pars.mp.segsitesin = atoi(argv[arg++]);
            }
            break;
        case 'F':
            arg++;
            argcheck(arg, argc, argv);
            pars.mp.mfreq = atoi(argv[arg++]);
            if ((pars.mp.mfreq < MIN_MFREQ) || (pars.mp.mfreq > pars.cp.nsam / 2)) {
                fprintf(stderr, " mfreq must be >= 2 and <= nsam/2.\n");
                usage();
            }
            break;
        case 'T':
            pars.mp.treeflag = 1;
            arg++;
            break;
        case 'L':
            pars.mp.timeflag = 1;
            arg++;
            break;
        case 'a':
            pars.mp.ageflag = 1;
            arg++;
            break;
        case 'I':
            arg = parse_island_model(arg, argc, argv, &npop);
            break;
        case 'm':
            if (npop < 2) {
                fprintf(stderr, "Must use -I option first.\n");
                usage();
            }
            if (argv[arg][2] == 'a') {
                arg++;
                for (pop = 0; pop < npop; pop++)
                    for (pop2 = 0; pop2 < npop; pop2++) {
                        argcheck(arg, argc, argv);
                        pars.cp.mig_mat[pop][pop2] = atof(argv[arg++]);
                    }
                calc_matrix_diagonal(pars.cp.mig_mat, npop);
            } else {
                arg++;
                argcheck(arg, argc, argv);
                i = atoi(argv[arg++]) - 1;
                argcheck(arg, argc, argv);
                j = atoi(argv[arg++]) - 1;
                argcheck(arg, argc, argv);
                mij = atof(argv[arg++]);
                pars.cp.mig_mat[i][i] += mij - pars.cp.mig_mat[i][j];
                pars.cp.mig_mat[i][j] = mij;
            }
            break;
        case 'n':
            if (npop < 2) {
                fprintf(stderr, "Must use -I option first.\n");
                usage();
            }
            arg++;
            argcheck(arg, argc, argv);
            pop = atoi(argv[arg++]) - 1;
            argcheck(arg, argc, argv);
            psize = atof(argv[arg++]);
            pars.cp.size[pop] = psize;
            break;
        case 'g':
            if (npop < 2) {
                fprintf(stderr, "Must use -I option first.\n");
                usage();
            }
            arg++;
            argcheck(arg, argc, argv);
            pop = atoi(argv[arg++]) - 1;
            if (arg >= argc) {
                fprintf(stderr, "Not enough arg's after -G.\n");
                usage();
            }
            palpha = atof(argv[arg++]);
            pars.cp.alphag[pop] = palpha;
            break;
        case 'G':
            arg++;
            if (arg >= argc) {
                fprintf(stderr, "Not enough arg's after -G.\n");
                usage();
            }
            palpha = atof(argv[arg++]);
            for (i = 0; i < pars.cp.npop; i++)
                pars.cp.alphag[i] = palpha;
            break;
        case 'e':
            arg = parse_demographic_event(arg, argc, argv);
            break;
        default:
            fprintf(stderr, " option default\n");
            usage();
        }
    }

    if ((pars.mp.theta == 0.0) && (pars.mp.segsitesin == 0) && (pars.mp.treeflag == 0) &&
        (pars.mp.timeflag == 0)) {
        fprintf(stderr, " either -s or -t or -T option must be used. \n");
        usage();
        exit(1);
    }

    sum = 0;
    for (i = 0; i < pars.cp.npop; i++)
        sum += (pars.cp.config)[i];

    if (sum != pars.cp.nsamin) { /* adna */
        fprintf(stderr, " sum sample sizes != nsam\n");
        usage();
        exit(1);
    }
}

static void set_default_params(void)
{
    pars.commandlineseedflag = 0;
    pars.output_precision = OUTPUT_PRECISION;
    pars.cp.r = pars.mp.theta = pars.cp.f = 0.0;
    pars.cp.track_len = 0.;
    pars.cp.npop = 1;
    pars.cp.mig_mat = (double **)malloc((unsigned)sizeof(double *));
    pars.cp.mig_mat[0] = (double *)malloc((unsigned)sizeof(double));
    pars.cp.mig_mat[0][0] = 0.0;
    pars.mp.segsitesin = 0;
    pars.mp.treeflag = 0;
    pars.mp.timeflag = 0;
    pars.mp.ageflag = 0;
    pars.mp.mfreq = 1;
    pars.cp.config = (int *)malloc((unsigned)((pars.cp.npop + 1) * sizeof(int)));
    /* adna */
    (pars.cp.config)[0] = pars.cp.nsamin;
    pars.cp.size = (double *)malloc((unsigned)(pars.cp.npop * sizeof(double)));
    (pars.cp.size)[0] = 1.0;
    pars.cp.alphag = (double *)malloc((unsigned)((pars.cp.npop) * sizeof(double)));
    (pars.cp.alphag)[0] = 0.0;
    pars.cp.nsites = MIN_NSITES;
}

static void calc_matrix_diagonal(double **mat, int npop)
{
    int pop, pop2;

    for (pop = 0; pop < npop; pop++) {
        mat[pop][pop] = 0.0;
        for (pop2 = 0; pop2 < npop; pop2++) {
            if (pop2 != pop)
                mat[pop][pop] += mat[pop][pop2];
        }
    }
}

// Procesa la opción -I npop n1 n2 ... [mig_rate] del modelo de islas, recibe en arg la posición de "-I" y devuelve la posición del siguiente argumento
static int parse_island_model(int arg, int argc, char *const argv[], int *pnpop)
{
    int i, j;
    double migr;

    arg++;
    if (count == 0) {
        argcheck(arg, argc, argv);
        pars.cp.npop = atoi(argv[arg]);
        pars.cp.config =
            (int *)realloc(pars.cp.config, (unsigned)(pars.cp.npop * sizeof(int)));
        *pnpop = pars.cp.npop;
    }
    arg++;
    for (i = 0; i < pars.cp.npop; i++) {
        argcheck(arg, argc, argv);
        pars.cp.config[i] = atoi(argv[arg++]);
    }
    if (count == 0) {
        pars.cp.mig_mat = (double **)realloc(pars.cp.mig_mat,
                                             (unsigned)(pars.cp.npop * sizeof(double *)));

        pars.cp.mig_mat[0] = (double *)realloc(pars.cp.mig_mat[0],
                                               (unsigned)(pars.cp.npop * sizeof(double)));

        for (i = 1; i < pars.cp.npop; i++)
            pars.cp.mig_mat[i] =
                (double *)malloc((unsigned)(pars.cp.npop * sizeof(double)));

        pars.cp.size =
            (double *)realloc(pars.cp.size, (unsigned)(pars.cp.npop * sizeof(double)));

        pars.cp.alphag =
            (double *)realloc(pars.cp.alphag, (unsigned)(pars.cp.npop * sizeof(double)));

        for (i = 1; i < pars.cp.npop; i++) {
            (pars.cp.size)[i] = (pars.cp.size)[0];
            (pars.cp.alphag)[i] = (pars.cp.alphag)[0];
        }
    }
    if ((arg < argc) && (argv[arg][0] != '-')) {
        argcheck(arg, argc, argv);
        migr = atof(argv[arg++]);
    } else
        migr = 0.0;

    for (i = 0; i < pars.cp.npop; i++)
        for (j = 0; j < pars.cp.npop; j++)
            pars.cp.mig_mat[i][j] = migr / (pars.cp.npop - 1);

    for (i = 0; i < pars.cp.npop; i++)
        pars.cp.mig_mat[i][i] = migr;
    return (arg);
}

static int parse_demographic_event(int arg, int argc, char *const argv[])
{
    int i, pop, npop2;
    struct devent *pt;
    char ch3;

    pt = (struct devent *)malloc(sizeof(struct devent));
    pt->detype = argv[arg][2];
    ch3 = argv[arg][3];
    arg++;
    argcheck(arg, argc, argv);
    pt->time = atof(argv[arg++]);
    pt->nextde = NULL;

    insert_event_sorted(pt);

    switch (pt->detype) {
    case 'N':
        argcheck(arg, argc, argv);
        pt->paramv = atof(argv[arg++]);
        break;
    case 'G':
        if (arg >= argc) {
            fprintf(stderr, "Not enough arg's after -eG.\n");
            usage();
        }
        pt->paramv = atof(argv[arg++]);
        break;
    case 'M':
        argcheck(arg, argc, argv);
        pt->paramv = atof(argv[arg++]);
        break;
    case 'n':
        argcheck(arg, argc, argv);
        pt->popi = atoi(argv[arg++]) - 1;
        argcheck(arg, argc, argv);
        pt->paramv = atof(argv[arg++]);
        break;
    case 'g':
        argcheck(arg, argc, argv);
        pt->popi = atoi(argv[arg++]) - 1;
        if (arg >= argc) {
            fprintf(stderr, "Not enough arg's after -eg.\n");
            usage();
        }
        pt->paramv = atof(argv[arg++]);
        break;
    case 's':
        argcheck(arg, argc, argv);
        pt->popi = atoi(argv[arg++]) - 1;
        argcheck(arg, argc, argv);
        pt->paramv = atof(argv[arg++]);
        break;
    case 'm':
        if (ch3 == 'a') {
            pt->detype = 'a';
            argcheck(arg, argc, argv);
            npop2 = atoi(argv[arg++]);
            pt->mat = (double **)malloc((unsigned)npop2 * sizeof(double *));
            for (pop = 0; pop < npop2; pop++) {
                (pt->mat)[pop] = (double *)malloc((unsigned)npop2 * sizeof(double));
                for (i = 0; i < npop2; i++) {
                    if (i == pop) {
                        arg++;
                    } else {
                        argcheck(arg, argc, argv);
                        (pt->mat)[pop][i] = atof(argv[arg++]);
                    }
                }
            }
            calc_matrix_diagonal(pt->mat, npop2);
        } else {
            argcheck(arg, argc, argv);
            pt->popi = atoi(argv[arg++]) - 1;
            argcheck(arg, argc, argv);
            pt->popj = atoi(argv[arg++]) - 1;
            argcheck(arg, argc, argv);
            pt->paramv = atof(argv[arg++]);
        }
        break;
    case 'j':
        argcheck(arg, argc, argv);
        pt->popi = atoi(argv[arg++]) - 1;
        argcheck(arg, argc, argv);
        pt->popj = atoi(argv[arg++]) - 1;
        break;
        /* adna */
    case 'A':
        argcheck(arg, argc, argv);
        pt->popi = atoi(argv[arg++]) - 1;
        argcheck(arg, argc, argv);
        pt->num = atoi(argv[arg++]);
        pars.cp.nsam += pt->num;
        break;
    default:
        fprintf(stderr, "e event\n");
        usage();
    }
    return (arg);
}

static void argcheck(int arg, int argc, char *const argv[])
{
    if ((arg >= argc) || (argv[arg][0] == '-')) {
        fprintf(stderr, "not enough arguments after %s\n", argv[arg - 1]);
        fprintf(stderr, "For usage type: ms<return>\n");
        exit(0);
    }
}

static void usage()
{
    fprintf(stderr, "usage: ms nsam howmany \n");
    fprintf(stderr, "  Options: \n");
    fprintf(stderr, "\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
    fprintf(stderr, "\t -s segsites   ( fixed number of segregating sites)\n");
    fprintf(stderr, "\t -T          (Output gene tree.)\n");
    fprintf(stderr, "\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
    fprintf(stderr, "\t -r rho nsites     (rho here is 4Nc)\n");
    fprintf(stderr, "\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is "
                    "mean length.) \n");
    fprintf(stderr, "\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n");
    fprintf(stderr, "\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");
    fprintf(
        stderr,
        "\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");
    fprintf(stderr, "\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n");
    fprintf(stderr, "\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of "
                    "migration matrix.)\n");
    fprintf(stderr, "\t\t -n i size_i   (popi has size set to size_i*N0 \n");
    fprintf(stderr, "\t\t -g i alpha_i  (If used must appear after -M option.)\n");
    fprintf(stderr, "\t   The following options modify parameters at the time 't' specified as the "
                    "first argument:\n");
    fprintf(stderr, "\t -eG t alpha  (Modify growth rate of all pop's.)\n");
    fprintf(stderr, "\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");
    fprintf(stderr,
            "\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n");
    fprintf(stderr, "\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n");
    fprintf(stderr, "\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements "
                    "of migration matrix.)\n");
    fprintf(stderr, "\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");
    fprintf(stderr,
            "\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
    fprintf(stderr,
            "\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");
    fprintf(stderr, "\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are "
                    "admixt. proport.\n");
    fprintf(stderr, "\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i "
                    "are unchanged.\n");
    fprintf(stderr, "\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
    fprintf(stderr, "\t\t  size, alpha and M are unchanged.\n");
    fprintf(stderr, "\t -eA t popi num ( num Ancient samples from popi at time t\n");
    fprintf(stderr, "\t  -f filename     ( Read command line arguments from file filename.)\n");
    fprintf(stderr, "\t  -p n ( Specifies the precision of the position output.  n is the number "
                    "of digits after the decimal.)\n");
    fprintf(stderr, "\t  -a (Print out allele ages.)\n");
    fprintf(stderr, " See msdoc.pdf for explanation of these parameters.\n");

    exit(1);
}

static void insert_event_sorted(struct devent *pt)
{
    struct devent *plast, *pevent;

    if ((pars.cp.deventlist == NULL) || (pt->time < pars.cp.deventlist->time)) {
        pt->nextde = pars.cp.deventlist;
        pars.cp.deventlist = pt;
        return;
    }
    pevent = pars.cp.deventlist;
    while ((pevent != NULL) && (pevent->time <= pt->time)) {
        plast = pevent;
        pevent = pevent->nextde;
    }
    plast->nextde = pt;
    pt->nextde = pevent;
}

static void free_event_list(struct devent *pt, int npop)
{
    struct devent *next;
    int pop;

    while (pt != NULL) {
        next = pt->nextde;
        if (pt->detype == 'a') {
            for (pop = 0; pop < npop; pop++)
                free((pt->mat)[pop]);
            free(pt->mat);
        }
        free(pt);
        pt = next;
    }
}

/************ make_gametes.c  *******************************************
 *
 *
 *****************************************************************************/

#define STATE1 '1' // Mutación
#define STATE2 '0' // Conservación del gen

static void make_gametes(int nsam, int mfreq, const struct node *ptree, double tt, int newsites, int ns, char **list)
{
    int tip, j, node;
    for (j = ns; j < ns + newsites; j++) {
        if (mfreq == 1)
            node = pick_branch(nsam, ptree, tt);
        else
            node = pick_branch_freq(nsam, mfreq, ptree, tt);
        agevec[j] = alleleage;
        for (tip = 0; tip < nsam; tip++) {
            if (is_descendant(ptree, tip, node))
                list[tip][j] = STATE1;
            else
                list[tip][j] = STATE2;
        }
    }
}

/***  tree_time.c : Returns the total time in the tree, *ptree, with nsam tips. **/

static double tree_time(const struct node *ptree, int nsam)
{
    double t;
    int i;
    t = (ptree + 2 * nsam - 2)->time;
    for (i = nsam; i < 2 * nsam - 1; i++)
        t += (ptree + i)->time;
    /* adna */
    for (i = 0; i < nsam; i++)
        t -= (ptree + i)->time;

    return (t);
}

static double tree_time_mfreq(const struct node *ptree, int nsam, int mfreq)
{
    double t;
    int i;

    t = 0.;
    for (i = 0; i < 2 * nsam - 2; i++) {
        if (((ptree + i)->ndes >= mfreq) && ((ptree + i)->ndes <= nsam - mfreq)) {
            t += branch_length(ptree, i);
        }
    }

    return (t);
}

static void print_tree(const struct node *ptree, int nsam)
{
    int i, *left_child, *right_child;

    left_child = (int *)malloc((unsigned)(2 * nsam - 1) * sizeof(int));
    right_child = (int *)malloc((unsigned)(2 * nsam - 1) * sizeof(int));
    for (i = 0; i < 2 * nsam - 1; i++) {
        left_child[i] = right_child[i] = -1;
    }

    for (i = 0; i < 2 * nsam - 2; i++) {
        if (left_child[(ptree + i)->abv] == -1) {
            left_child[(ptree + i)->abv] = i;
        } else {
            right_child[(ptree + i)->abv] = i;
        }
    }

    print_parenthesis(ptree, left_child, right_child, 2 * nsam - 2);
    free(left_child);
    free(right_child);
}

// Función para imprimir los arboles en formato de Newick
static void print_parenthesis(const struct node *ptree, const int *left_child, const int *right_child, int noden)
{
    if (left_child[noden] == -1) {
        printf("%d:%5.3lf", noden + 1, branch_length(ptree, noden)); /* adna */
    } else {
        printf("(");
        print_parenthesis(ptree, left_child, right_child, left_child[noden]);
        printf(",");
        print_parenthesis(ptree, left_child, right_child, right_child[noden]);
        if ((ptree + noden)->abv == 0) {
            printf(");\n");
        } else {
            printf("):%5.3lf", branch_length(ptree, noden));
        }
    }
}

/* Duración de la rama que une al nodo i con su nodo padre */
static double branch_length(const struct node *ptree, int i)
{
    return (ptree + (ptree + i)->abv)->time - (ptree + i)->time;
}

/* Instante elegido al azar (uniforme) dentro de la rama del nodo i */
static double random_age_in_branch(const struct node *ptree, int i)
{
    return (ptree + i)->time + ran1() * branch_length(ptree, i);
}

/***  pick_branch : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
              time in tree.   ****/

static int pick_branch(int nsam, const struct node *ptree, double tt)
{
    double x, y;
    int i;

    x = ran1() * tt;

    for (i = 0, y = 0; i < 2 * nsam - 2; i++) {
        y += branch_length(ptree, i);
        if (y >= x) {
            alleleage = random_age_in_branch(ptree, i);
            return (i);
        }
    }

    alleleage = random_age_in_branch(ptree, 2 * nsam - 3);

    return (2 * nsam - 3); /* changed 4 Feb 2010 */
}

static int pick_branch_freq(int nsam, int mfreq, const struct node *ptree, double tt)
{
    double x, y;
    int i, lastbranch = 0;

    x = ran1() * tt;

    for (i = 0, y = 0; i < 2 * nsam - 2; i++) {
        if (((ptree + i)->ndes >= mfreq) && ((ptree + i)->ndes <= nsam - mfreq)) {
            y += branch_length(ptree, i);
            lastbranch = i; /* changed 4 Feb 2010 */
        }
        if (y >= x) {
            alleleage = random_age_in_branch(ptree, i);
            return (i);
        }
    }

    alleleage = random_age_in_branch(ptree, lastbranch);

    return (lastbranch); /*  changed 4 Feb 2010 */
}

/****  is_descendant : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/

static int is_descendant(const struct node *ptree, int tip, int node)
{
    int k;

    for (k = tip; k < node; k = (ptree + k)->abv)
        ;

    if (k == node) {
        return (1);
    } else {
        return (0);
    }
}

/**** ordran.c  ***/

static void ordran(int n, double pbuf[])
{
    ranvec(n, pbuf);
    order(n, pbuf);
    return;
}

static void multinomial(int n, int nclass, const double p[], int rv[])
{
    double x, s;
    int i, j;

    for (i = 0; i < nclass; i++) {
        rv[i] = 0;
    }
    for (i = 0; i < n; i++) {
        x = ran1();
        j = 0;
        s = p[0];
        while ((x > s) && (j < (nclass - 1)))
            s += p[++j];
        rv[j]++;
    }
    return;
}

static void order(int n, double pbuf[])
{
    int gap, i, j;
    double temp;

    for (gap = n / 2; gap > 0; gap /= 2)
        for (i = gap; i < n; i++)
            for (j = i - gap; j >= 0 && pbuf[j] > pbuf[j + gap]; j -= gap) {
                temp = pbuf[j];
                pbuf[j] = pbuf[j + gap];
                pbuf[j + gap] = temp;
            }
    return;
}

static void ranvec(int n, double pbuf[])
{
    int i;

    for (i = 0; i < n; i++) {
        pbuf[i] = ran1();
    }

    return;
}

static int poisso(double u)
{
    double cump, ru, p;
    int i = 1;

    if (u > 30.) {
        i = (int)(0.5 + gasdev(u, u));
        if (i < 0)
            return (0);
        else
            return (i);
    }

    ru = ran1();
    p = exp(-u);
    if (ru < p)
        return (0);
    cump = p;

    while (ru > (cump += (p *= u / i))) {
        i++;
    }
    return (i);
}

/* a slight modification of crecipes version */

static double gasdev(double m, double v)
{
    static int iset = 0;
    static float gset;
    float fac, r, v1, v2;

    if (iset == 0) {
        do {
            v1 = 2.0 * ran1() - 1.0;
            v2 = 2.0 * ran1() - 1.0;
            r = v1 * v1 + v2 * v2;
        } while (r >= 1.0);

        fac = sqrt(-2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;

        return (m + sqrt(v) * v2 * fac);
    } else {
        iset = 0;
        return (m + sqrt(v) * gset);
    }
}
