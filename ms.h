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

#include <stdlib.h>

#define LIKELY(x) __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)

extern unsigned short ms_xsubi[3];
#pragma omp threadprivate(ms_xsubi)

static inline __attribute__((always_inline)) double ran1(void)
{
    return erand48(ms_xsubi);
}

void ms_seed_thread(void);
void seedit(const char *flag);
int commandlineseed(char *const *seeds);
struct segl *simulate_coalescent(const struct c_params *cp, int *pnsegs);
