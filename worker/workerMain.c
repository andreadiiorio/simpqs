#include "SIMPQS.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <CONFIGURATION.h>
#include <matrix/matrix.h>
#include <factorization/factorizerQuick.h>
#include <unistd.h>
#include <wait.h>
#include "utils/utils.h"
#include "utils/gmp_patch.h"
#include "sievingSIMPQS.h"
#include "../matrix/matrix.h"
SIEVE_ARRAY_BLOCK SieveArrayBlock;              ///array block in memory
const char* n_str="100000030925519250968982645360649";
CONFIGURATION* configuration;
DYNAMIC_VECTOR* primes_B;
#define TEST_1_POL_FAMILY
#ifdef TEST_1_POL_FAMILY
int main(){
#else
int finalStepWrap(){
#endif

    configuration = initConfiguration(n_str, 0, 0, 0, 0);
    struct polynomial pol;
    PRECOMPUTES *precomputes = preComputations(configuration, &pol,NULL);
    char** reportsPaths=findReportsLocally(1,REPORTS_POLYNOMIAL_FAMILY_FILENAME_SUFFIX);
    REPORTS* reports=loadReports(*reportsPaths);
    if(reports->relationsNum<(precomputes->primes.vectorSize+1)) {
        fprintf(stderr, "not founded enough reports for linear algebra phase\n");
        exit(EXIT_FAILURE);
    }
    MATRIX matrix;
    initMatrixFromRows(&matrix,reports->relationsNum,reports->bsmoothEntries,precomputes->primes.vectorSize+1);
    gauss_elimination_over_GF2(&matrix);
//    print_matrix_matrix(&matrix);
//    print_matrix_identity(&matrix);
    return quadraticRelationTry(reports,&matrix);
}
#ifdef TEST_1_POL_FAMILY
int MAIN(){
#else
int main(int argc, char** argv){
    //FORKED WORKER PROCESS TO SIEVE BSMOOTH RELATION AND PARTIAL RELATION WITH POLINOMIO FAMILIES FROM a
    if(argc < 5){           //TODO MOCKED ARGv
        printf("USAGE: N,M,B,a, ....\n");
        //exit(EXIT_FAILURE);
    }
#endif
    /// init configuration, with argv, on un setted configuration defalut setting will be used
    configuration = initConfiguration(n_str, 0, 0, 0, 0);
    struct polynomial pol;
    /// init precomputation for polynomial family   todo next version computation of a coeff ... <- MASTER COMUNICATION
    PRECOMPUTES *precomputes = preComputations(configuration, &pol,NULL);
    if (!precomputes ){
        free(configuration);
        exit(EXIT_FAILURE);
    }
    gmp_printf("\n DONE precomputation for polynomial: a: %Zd b: %Zd factorizing N: %Zd \n", pol.a, pol.b, *(pol.N));

    // init local worker matrix, used to aggregate relations founded in the various sieving iteration
//    MATRIX matrixWorker;
//    init_matrix(&matrixWorker,1,precomputes->primes.vectorSize);

#ifdef VERBOSE
    printPrecomputations(Precomputations,5);
#endif
    /*
     * SIEVE CONCURRENTLY ON DIFFERENT PROCESSES
     * Generated Mongomery Polynomial Coefficients belonging to polynomial family of Master passed coeff a
     * Sieving will be Done with each polynomial on a different process
     * Each Process will concurrently sieve on a block of the sieving array  with pthreads
     * and factorization of each LikellyBSmooth entry (identified by log sieving) will be done with 2 level of parallelism by FactorizerThreadGroupsPool
     * before each Process termination, the founded (partial) reports will be serialized to file with filename REPORTS_FILENAME_PREFIX N a b #rel #partialRel REPORTS_FILENAME_SUFFIX
     */
    const unsigned int polynoamilFamilySize=(1 << (configuration->a_coefficient.a_factors_num - 1));
//    int* polynomial_processes_pid=malloc(polynoamilFamilySize* sizeof(*polynomial_processes_pid));
//    if(!polynomial_processes_pid){
//        free(configuration);free(precomputes);
//        exit(EXIT_FAILURE);
//    }
    //// polynomials family generation:
#ifdef DEBUG
    REPORTS *reportsFounded = Sieve(configuration, precomputes, &pol);
    saveReports(reportsFounded, precomputes->factorbaseSize, true, &pol);
        exit(!reportsFounded);
#endif
    //  TODO DEBUG FORKED CHILD: set follow-fork-mode child \n set detach-on-fork off
    for (u_int j = 1; j  < polynoamilFamilySize; ++j) {
        gmp_printf("\n\n\n\npolynomial:%d\ta=%Zd;\tb=%Zd;\n", j, pol.a, pol.b);
        if (!fork()) {
            REPORTS *reportsFounded = Sieve(configuration, precomputes, &pol);
            int result = EXIT_SUCCESS;
            if (!reportsFounded) {
                fprintf(stderr, "SIEVING ERROR\n");
                result = EXIT_FAILURE;
            }
            if ((saveReports(reportsFounded, precomputes->factorbaseSize, false,false, &pol)) == EXIT_FAILURE) {
                fprintf(stderr, "REPORTS SERIALIZING ERROR\n");
                result = EXIT_FAILURE;
            }
            free(reportsFounded);
            exit(result);
        }
        nextPolynomial_b_i(&(pol.b), j, precomputes);                           //change polynomial
    }

    for (u_int j = 1; j  < polynoamilFamilySize; ++j) {
        int workerPolynomialRes=0;
        wait(&workerPolynomialRes);
        printf("polynomial process returned: %s\n",workerPolynomialRes==EXIT_SUCCESS?"Exit success":"Exit failure");
    }
    fflush(0);
    REPORTS* polynomialsReportsAggregated=aggregateSieversWorkers(polynoamilFamilySize);
    if(!polynomialsReportsAggregated){
        fprintf(stderr,"ERR DURING POLYNOMIAL REPORTS AGGREGATION");
        exit(EXIT_FAILURE);
    }
    printf("\n\n\n\n\n\nSERIALIZING AGGREGATED REPORTS\n");fflush(0);
    if (saveReports(polynomialsReportsAggregated,precomputes->primes.vectorSize,false,true,&pol)==EXIT_FAILURE){
        fprintf(stderr, "REPORTS  AGGREGATED SERIALIZING ERROR\n");
        exit(EXIT_FAILURE);
    }


    //TODO MATRIX STAGE HERE --> NEXT MOVE TO MASTER
#ifndef TEST_1_POL_FAMILY
    return finalStepWrap();
#else
    exit(EXIT_SUCCESS);
#endif
}

REPORTS* aggregateSieversWorkers(const unsigned int polynomialN) {
    int result=EXIT_SUCCESS;
    char** reportsLocalFilenames=findReportsLocally(polynomialN,REPORTS_FILENAME_SUFFIX);
    REPORTS** reportsAll=malloc(sizeof(*reportsAll)*polynomialN);       //will hold reports pointers
    REPORTS* reportsAllMerged=calloc(1,sizeof(*reportsAllMerged));       //will hold reports pointers
    if(!reportsLocalFilenames || !reportsAll || !reportsAllMerged){
        fprintf(stderr,"Out of Mem or reports paths retrive err\n");
        result =EXIT_FAILURE;goto exit;
    }
    char* reportFilePath=*(reportsLocalFilenames);
    unsigned int report_i;
    u_int64_t cumulativeReportsN=0,cumulativePartialReportsN=0;
    for (report_i = 0; report_i < polynomialN && reportFilePath; reportFilePath=reportsLocalFilenames[++report_i]) {
        reportsAll[report_i]= loadReports(reportFilePath);              //on error reports_i will set to NULL
        if(!reportsAll[report_i]){
            fprintf(stderr,"invalid read\n");
            result =EXIT_FAILURE;goto exit;

        }
        cumulativeReportsN+=reportsAll[report_i]->relationsNum;
        cumulativePartialReportsN+=reportsAll[report_i]->partialRelationsNum;
        if(mergeReports(reportsAllMerged,reportsAll[report_i])==EXIT_FAILURE){
            fprintf(stderr,"merge reports error\n");
            result =EXIT_FAILURE;goto exit;
        }
    }
    //// allocate enough space for merged reports
    /*reportsAllMerged->bsmoothEntries=malloc(sizeof(*(reportsAllMerged->bsmoothEntries))*cumulativeReportsN);
    reportsAllMerged->largePrimesEntries=malloc(sizeof(*(reportsAllMerged->largePrimesEntries))*cumulativePartialReportsN);
    if(!reportsAllMerged->bsmoothEntries ||  !reportsAllMerged->largePrimesEntries){
        fprintf(stderr,"Out of Mem at reports aggregation allocation\n");
        result =EXIT_FAILURE;goto exit;
    }
    //// aggregate all founded reports in 1! reports struct
    for (u_int i = 0; i < report_i; ++i) {
        REPORTS* reports=reportsAll[i];
        if(reports){
            if(mergeReports(reportsAllMerged, reports) == EXIT_FAILURE){
                fprintf(stderr,"error at aggregate founded reports\n");
                result =EXIT_FAILURE;goto exit;
            }
        }
    }*/
    //// try to pair partial reports in extra reports
    if(pairPartialReports(reportsAllMerged)==EXIT_FAILURE){
        fprintf(stderr,"error during partial reports pairing\n");
        result =EXIT_FAILURE;goto exit;                                 ///TODO SKIPPABLE EXIT
    }
    exit:
        free(reportsAll);free(reportsLocalFilenames);//TODO FREE REPORTS FILE PATHS NESTED POINTERS
        if(result==EXIT_FAILURE){
            if(reportsAllMerged){
                free(reportsAllMerged->largePrimesEntries);free(reportsAllMerged->bsmoothEntries);
                return NULL;
            }
            free(reportsAllMerged);
        }
    return reportsAllMerged;
}
