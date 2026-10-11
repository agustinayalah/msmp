struct devent {
    double time;
    int popi;
    int popj;
    double paramv;
    int num; /* adna */
    double **mat;
    char detype;
    struct devent *nextde;
};

struct c_params { // Parámetros de coalescencia
    int npop;
    int nsam;   /* total sample size including ancient samples */
    int nsamin; /*  sample size at present time,  does not include ancient samples */
    int *config;
    double **mig_mat;
    double r;
    int nsites;
    double f;
    double track_len;
    double *size;
    double *alphag;
    struct devent *deventlist;
};

struct m_params { // Parámetros de mutación
    double theta;
    int segsitesin; // Para fijar cantidad de segsites
    int treeflag;
    int timeflag;
    int mfreq;
    int ageflag;
};

struct params {
    struct c_params cp;
    struct m_params mp;
    int commandlineseedflag;
    int output_precision;
};

struct node {
    int abv;
    int ndes;
    float time;
};

struct segl {
    int beg;
    struct node *ptree;
    int next;
};

/*KRT -- prototypes added*/
void ordran(int n, double pbuf[]);
void ranvec(int n, double pbuf[]);
void order(int n, double pbuf[]);

void bigger_list(int nsam, char **list);
int poisso(double u);
void locate(int n, double beg, double len, double *ptr);
void multinomial(int n, int nclass, double p[], int rv[]);
void usage(void);
int is_descendant(struct node *ptree, int tip, int node);
int pick_2_random_numbers(int n, int *i, int *j);
int crossover(int nsam, int ic, int is);
int chrom_links(int c);
void seedit(const char *flag);
struct segl *simulate_coalescent(struct c_params *cp, int *pnsegs);
int commandlineseed(char **seeds);
double ran1(void);