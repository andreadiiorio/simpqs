#include "SIMPQS.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <CONFIGURATION.h>
#include "utils/utils.h"
#include "utils/gmp_patch.h"
#include "sievingSIMPQS.h"

mpz_t   N;
SIEVE_ARRAY_BLOCK SieveArrayBlock;              ///array block in memory
struct Configuration Config;

///polinomial family coeff
struct polynomial_actual* ActualPolynomial;
void initializationVars(char** argv){
    /*
     generated Pol Coeff a: 16079919820099 b: 19000050618295 c: -6218915954954063976
 N: 100000030925519250968982645360649 B: 17330 M: 7330

     */
    if(!(ActualPolynomial=malloc(sizeof(struct polynomial_actual)))){
        fprintf(stderr,"malloc failed for actual polynomio instantiation");
        exit(EXIT_FAILURE);
    }
    //TODO MOCKED ARGV
    //char* _N=argv[1];
    //char* _M=argv[2];
    //char* _B=argv[3];
    //char* _a=argv[4];
    Config.M=7096;
    Config.B=10000;
    Config.SIEVING_THREAD_NUM=_SIEVING_THREAD_NUM;
    Config.ARRAY_IN_MEMORY_MAX_SIZE=_ARRAY_IN_MEMORY_MAX_SIZE;
    char* _N="100000030925519250968982645360649";
    char* _a="16079919820099";
    char* _b="19000050618295";
    char* _c="-6218915954954063976";
    mpz_init(N);
    strToMpz(N,_N)
    mpz_init((*ActualPolynomial).a);
    strToMpz((*ActualPolynomial).a,_a)
    mpz_init((*ActualPolynomial).b);
    strToMpz((*ActualPolynomial).b,_b)
    mpz_init((*ActualPolynomial).c);
    strToMpz((*ActualPolynomial).c,_c)
}
int main(int argc, char** argv){
    //FORKED WORKER PROCESS TO SIEVE BSMOOTH RELATION AND PARTIAL RELATION WITH POLINOMIO FAMILIES FROM a
    if(argc < 5){
        printf("USAGE: N,M,B,a, ....\n"); //todo append config
        //exit(EXIT_FAILURE);
    }
    initializationVars(argv);
    struct Precomputes* precomputations=preComputations(&((*ActualPolynomial).a));
    gmp_printf("\n DONE precomputation for polynomial: a: %Zd b: %Zd factorizing N: %Zd \n",ActualPolynomial->a,ActualPolynomial->b,N);
#ifdef VERBOSE
    printPrecomputations(Precomputations,5);
#endif
    /// sieve array allocate
    if(!(SieveArrayBlock=malloc(Config.ARRAY_IN_MEMORY_MAX_SIZE* sizeof(*SieveArrayBlock)))){   //TODO DEBUG ALLOCATED SIZE
        fprintf(stderr,"sieve array block in mem malloc failed\n");
        exit(EXIT_FAILURE);
    }
    if(Sieve(&Config,precomputations,SieveArrayBlock)==EXIT_FAILURE){
        fprintf(stderr,"sieve error \n");
        exit(EXIT_FAILURE);
    } //TODO ITERATE UNTIL MATRIX IS FILLED --> EVENTUALLY CHANGE POLYNOMIAL

    //TODO MATRIX STAGE CALLBACK TO MASTER
}
