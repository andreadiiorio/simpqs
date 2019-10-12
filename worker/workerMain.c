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

#define POL_FAMILY_TRIES 150
///TODO CLEAN REPORTS FILES:  find -iname "reports_*" | xargs -d "\n" rm
SIEVE_ARRAY_BLOCK SieveArrayBlock;              ///array block in memory
//const char* n_str=  "100000030925519250968982645360649"; //OLD
//const char* n_str="100000030925519650969044496394369";
//const char* n_str="10000000000251715795601347229089999344259";
const char* n_str="100000000000000028598093420000002011524934548107677";

//TOO LARGE TRYS
//const char* n_str="100000000000000000000000000000024260646520000000000000000000001360395958358684667";
//const char* n_str="1000000000000000000002859809340000000000002011524921863497539";
CONFIGURATION* Configuration;
//#define TEST_1_POL_FAMILY

#ifdef TEST_1_POL_FAMILY
int main(){
#else
int finalStepWrap(){
#endif

    Configuration =initConfiguration(n_str, 0, 0, 0, 0);
    struct polynomial pol;
    PRECOMPUTES *precomputes = preComputations(Configuration, &pol, true);
    char** reportsPaths=findReportsLocally(1,REPORTS_POLYNOMIAL_FAMILY_FILENAME_SUFFIX);
    if(!(*reportsPaths)){
        fprintf(stderr,"error in reports deserialization\n");
        exit(EXIT_FAILURE);
    }
    REPORTS* reports=loadReports(*reportsPaths);
    if(reports->relationsNum<(precomputes->factorbaseSize)+1) {
        fprintf(stderr, "not founded enough reports for linear algebra phase\n");
        exit(EXIT_FAILURE);
    }
    MATRIX matrix;
    initMatrixFromRows(&matrix,reports->relationsNum,reports->bsmoothEntries,precomputes->factorbaseSize+1);
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
    _deleteLocalReports(false);          //TODO DEBUG RESET OLD SERIALIZED REPORTS
    if(argc < 5){           //TODO MOCKED ARGv
        printf("USAGE: N,M,B,a, ....\n");
        //exit(EXIT_FAILURE);
    }
#endif
    /// init configuration, with argv, on un setted configuration defalut setting will be used
    CONFIGURATION* configuration= initConfiguration(n_str, 0, 0, 0, 0);
    Configuration=configuration;
    struct polynomial pol;
    /// init precomputation for polynomial family   todo next version computation of a coeff ... <- MASTER COMUNICATION
    PRECOMPUTES *precomputes = preComputations(configuration, &pol, false);
    if (!precomputes ){
        free(Configuration);
        exit(EXIT_FAILURE);
    }

    gmp_printf("\n DONE precomputation for polynomial: a: %Zd b: %Zd factorizing N: %Zd \n", pol.a.a, pol.b, *(pol.N));
    // init local worker matrix, used to aggregate relations founded in the various sieving iteration
//    MATRIX matrixWorker;
//    init_matrix(&matrixWorker,1,precomputes->primes.vectorSize);

    REPORTS *polynomialsReportsAggregated;
    A_COEFF *polynomialFamilyCoefficients;
    polynomialFamilyCoefficients = genPolynomialFamilies_a(POL_FAMILY_TRIES, Configuration, precomputes, configuration->a_factors_indexes_FB_families);
    if(!polynomialFamilyCoefficients)
        exit(EXIT_FAILURE);


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

    REPORTS* reportsAllWorker=calloc(1, sizeof(*reportsAllWorker));
    if(!reportsAllWorker) {
        fprintf(stderr,"Out of mem at reports all calloc\n");
        exit(EXIT_FAILURE);
    }
    unsigned int polynoamilFamilySize;
    int polFamily_i = 1;

    polynomial_family_sieve:
    polynoamilFamilySize = (1 << (pol.a.a_factors_num - 1));
    char polynomial_str[440];
//    fflush(0); printf("starting siever process in polynomial family of size:_%d\n",polynoamilFamilySize);fflush(0);
    //  TODO DEBUG FORKED CHILD: set follow-fork-mode child \n set detach-on-fork off
    for (u_int j = 1; j  < polynoamilFamilySize; ++j) {
        gmp_snprintf(polynomial_str,440," polynomial:%d\ta=%Zd;\tb=%Zd;\n", j, pol.a.a, pol.b);
        if (!fork()) {
            REPORTS *reportsFounded = Sieve(Configuration, precomputes, &pol);
            int result = EXIT_SUCCESS;
            if (!reportsFounded) {
                fprintf(stderr, "SIEVING ERROR\n");
                result = EXIT_FAILURE;
            }
            if(reportsFounded->relationsNum || reportsFounded->partialRelationsNum ) {  //serialize reports only on need
                if ((saveReports(reportsFounded, precomputes->factorbaseSize, false, false, &pol)) == EXIT_FAILURE) {
                    fprintf(stderr, "REPORTS SERIALIZING ERROR\n");
                    result = EXIT_FAILURE;
                }
            }
            fprintf(stderr,"EXITING WITH CODE %d SIEVER PROCESS founded full reports:%lu\tpartial reports:%lu\tOF POLYNOMIAL:\t%s\n",result,reportsFounded->relationsNum,reportsFounded->partialRelationsNum, polynomial_str);
            free(reportsFounded);
//            freePrecomputations(precomputes);
            exit(result);
        }
        nextPolynomial_b_i(&(pol.b), j, precomputes);                           //change polynomial
    }

    for (u_int j = 1; j  < polynoamilFamilySize; ++j) {
        int workerPolynomialRes=0;
        wait(&workerPolynomialRes);
        if(workerPolynomialRes!=EXIT_SUCCESS)
            exit(workerPolynomialRes);
    }
    polynomialsReportsAggregated = aggregateSieversWorkers(polynoamilFamilySize, false,reportsAllWorker);
    if(!polynomialsReportsAggregated){
        fprintf(stderr,"ERR DURING POLYNOMIAL REPORTS AGGREGATION");
        exit(EXIT_FAILURE);
    }
    printf("\n\naggreagated reports:%lu	partialReports:%lu vs FB SIZE: %lu \n", polynomialsReportsAggregated->relationsNum,polynomialsReportsAggregated->partialRelationsNum,precomputes->factorbaseSize);
    if(polynomialsReportsAggregated->relationsNum < precomputes->factorbaseSize - 88 && (polFamily_i) < POL_FAMILY_TRIES){
        changePolynomialFamily(precomputes, polynomialFamilyCoefficients + polFamily_i++, &pol);
        gmp_printf("Changed to polynomial family :%d\t new a=%Zd;\tb=%Zd;\n",polFamily_i,pol.a.a,pol.b);
        goto polynomial_family_sieve;
    }
    if(pairPartialReports(reportsAllWorker)==EXIT_FAILURE)
        exit(EXIT_FAILURE);

    printf("\n\n\n\n\n\nSERIALIZING POLYNOMIAL FAMILIES AGGREGATED REPORTS\n");fflush(0);
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


REPORTS *aggregateSieversWorkers(const unsigned int polynomialN, bool pairLargePrimeEntries, REPORTS *reportsAllMerged) {
    //aggregate polynomial siver outputs, if reportsAllMerged is null will allocated space for it otherwise reports will be concat there
    //if pairLargePrimeEntries is true, large prime entries will be aggregated in new relations
    //returned reportsAllMerged (allocated if needed)


    int result=EXIT_SUCCESS;
    char** reportsLocalFilenames=findReportsLocally(polynomialN,REPORTS_FILENAME_SUFFIX);
    REPORTS** reportsAll=malloc(sizeof(*reportsAll)*polynomialN);
    if(!reportsAllMerged) reportsAllMerged=calloc(1, sizeof(*reportsAllMerged)); //allocate all reports if has been null passed
    if(!reportsLocalFilenames || !reportsAll || !reportsAllMerged){
        fprintf(stderr,"Out of Mem or reports paths retrive err\n");
        result =EXIT_FAILURE;goto exit;
    }
    char* reportFilePath=*(reportsLocalFilenames);          //filename of report file (init with first fo find cmd output
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
        if(mergeReports(reportsAllMerged, reportsAll[report_i], true) == EXIT_FAILURE){
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
    if(pairLargePrimeEntries && pairPartialReports(reportsAllMerged)==EXIT_FAILURE){
        fprintf(stderr,"error during partial reports pairing\n");
        result =EXIT_FAILURE;goto exit;                                 ///TODO SKIPPABLE EXIT
    }
    exit:
    free(reportsAll);free(reportsLocalFilenames);//TODO FREE REPORTS FILE PATHS NESTED POINTERS
    if(result==EXIT_FAILURE){
        if(reportsAllMerged){
//            free(reportsAllMerged->largePrimesEntries);free(reportsAllMerged->bsmoothEntries);
            freeReports(reportsAllMerged);
            return NULL;
        }
        free(reportsAllMerged);
    }
    _deleteLocalReports(true);
    return reportsAllMerged;
}
