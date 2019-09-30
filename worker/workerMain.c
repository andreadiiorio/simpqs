#include "SIMPQS.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <CONFIGURATION.h>
#include <matrix/matrix.h>
#include <factorization/factorizerQuick.h>
#include "utils/utils.h"
#include "utils/gmp_patch.h"
#include "sievingSIMPQS.h"

mpz_t   N;
SIEVE_ARRAY_BLOCK SieveArrayBlock;              ///array block in memory
struct Configuration Config;
const char* n_str="100000030925519250968982645360649";

int main(int argc, char** argv){
    //FORKED WORKER PROCESS TO SIEVE BSMOOTH RELATION AND PARTIAL RELATION WITH POLINOMIO FAMILIES FROM a
    if(argc < 5){
        printf("USAGE: N,M,B,a, ....\n"); //todo append config
        //exit(EXIT_FAILURE);
    }
    CONFIGURATION *configuration = initConfiguration(n_str, 0, 0, 0, 0);
    //// get first polynomial:
    struct polynomial pol;
    PRECOMPUTES *precomputes = preComputations(configuration, &pol);
    if (!precomputes )
        exit(EXIT_FAILURE);
    gmp_printf("\n DONE precomputation for polynomial: a: %Zd b: %Zd factorizing N: %Zd \n", pol.a, pol.b, N);

    // init local worker matrix, used to aggregate relations founded in the various sieving iteration
    //TODO AT POLYNOMIAL FAMILY END MATRIX WILL BE SERIALIZED AND SENT TO MASTER
    MATRIX matrixWorker;
    init_matrix(&matrixWorker,1,precomputes->primes.vectorSize);

#ifdef VERBOSE
    printPrecomputations(Precomputations,5);
#endif
    /// Sieve
    REPORTS* polynomialFamilyReports=calloc(1, sizeof(*polynomialFamilyReports));
    if (!polynomialFamilyReports){
        fprintf(stderr,"Out of MEM in pol family reports");
        exit(EXIT_FAILURE);
    }
    //TODO ITERATE CONCURRENTLY OVER BLOCK OF POLYNOMIALS in family passed by MASTER in configuration IN TASK
    if(!(SieveArrayBlock=malloc(Config.ARRAY_IN_MEMORY_MAX_SIZE* sizeof(*SieveArrayBlock)))){
        fprintf(stderr,"sieve array block in mem malloc failed\n");
        exit(EXIT_FAILURE);
    }
    REPORTS* polynomialReports=Sieve(&Config, precomputes, SieveArrayBlock, &pol);
    if (!polynomialReports){
        fprintf(stderr,"sieve error occurred\n");
        exit(EXIT_FAILURE);
    }
    if(mergeReports(polynomialFamilyReports,polynomialReports)==EXIT_FAILURE){
        fprintf(stderr,"REPORTS MERGE ERROR\n");
    }
}
