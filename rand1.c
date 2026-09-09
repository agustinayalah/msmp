/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

double ran1(void);
void seedit(char *flag);
int commandlineseed(char **seeds);
void ms_seed_thread(void);

static unsigned short master_seed[3] = { 3579, 27011, 59243 } ;
static unsigned short ms_xsubi[3] = { 3579, 27011, 59243 } ;

#pragma omp threadprivate(ms_xsubi)

double ran1(void)
{
    return (erand48(ms_xsubi));
}

void ms_seed_thread(void)
{
#ifdef _OPENMP
    unsigned tid = (unsigned) omp_get_thread_num();
#else 
    unsigned tid = 0 ;
#endif
	ms_xsubi[0] = (unsigned short)( master_seed[0] ^ (unsigned short)(tid * 0x9E37u + 0x1u) );
	ms_xsubi[1] = (unsigned short)( master_seed[1] ^ (unsigned short)(tid * 0x85EBu + 0x3u) );
	ms_xsubi[2] = (unsigned short)( master_seed[2] ^ (unsigned short)(tid * 0xC2B2u + 0x5u) );
}

void seedit(char *flag)
{
	FILE *fopen(), *pfseed;
	unsigned short seedv2[3] ;
	int i;
 
	if( flag[0] == 's' ) {
		pfseed = fopen("seedms","r");
		if( pfseed == NULL ) {
			master_seed[0] = 3579 ; master_seed[1] = 27011; master_seed[2] = 59243;
		}
		else {
			seedv2[0] = 3579; seedv2[1] = 27011; seedv2[2] = 59243;
			for(i=0;i<3;i++){
				if(  fscanf(pfseed," %hd",master_seed+i) < 1 )
					master_seed[i] = seedv2[i] ;
			}
			fclose( pfseed);
		}
		/* 
            Seed the calling thread too (used e.g. for the -f/tbs
		   sequential pre-pass, which runs before the parallel region). 
        */
		ms_seed_thread();
 
		printf("\n%d %d %d\n", master_seed[0], master_seed[1], master_seed[2] );
	}
	else {
		pfseed = fopen("seedms","w");
		if( pfseed != NULL ){
			fprintf(pfseed,"%d %d %d\n", master_seed[0], master_seed[1], master_seed[2] );
			fclose(pfseed);
		}
	}
}
 
int commandlineseed(char **seeds)
{
 
	master_seed[0] = (unsigned short) atoi( seeds[0] );
	master_seed[1] = (unsigned short) atoi( seeds[1] );
	master_seed[2] = (unsigned short) atoi( seeds[2] );
	printf("\n%d %d %d\n", master_seed[0], master_seed[1], master_seed[2] );
 
	ms_seed_thread();
	return(3);
}
 