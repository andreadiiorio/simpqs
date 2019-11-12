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

// developped by andysnake96
// Highly Concurrent SIMPQS with Large Prime variation
// based on Contini phD for SIMPQS Algo
// see readme for algo information, see CONFIGURATION.h for algorithm parameter setting (widely used c
#define POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS 5
#define POLYNOMIAL_FAMILIES_PER_SIEVER 7
#define LARGE_PRIMES_SURPLUS_EXPECTED 0
int POL_FAMILY_TRIES=950;


char* n;
//const char* n_str=  "100000030925519250968982645360649"; //OLD
//const char* n_str="100000030925519650969044496394369";
//const char* n_str="10000000000251715795601347229089999344259";   //~40
char* n_str= "100000000000000028598093420000002011524934548107677";//~50
//char* n_str= "1000000000000000000002426064040000000000001360395223980847119";//~60
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

    Configuration =initConfiguration(n, 0, 0, 0, 0);
    struct polynomial pol;
    PRECOMPUTES *precomputes = preComputations(Configuration, &pol, true);
    char** reportsPaths=findReportsLocally(1,REPORTS_POLYNOMIAL_FAMILY_FILENAME_SUFFIX);
    if(!(*reportsPaths)){
        fprintf(stderr,"error in reports deserialization\n");
        exit(EXIT_FAILURE);
    }
    REPORTS* reports=loadReports(*reportsPaths);
    if(reports->relationsNum<(precomputes->factorbaseSize)+1) {
        fprintf(stderr, "not founded enough reports for linear algebra phase \n founded:%lu vs %lu\n",reports->relationsNum,precomputes->factorbaseSize);
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

extern bool auditExtra;
int main(int argc, char** argv){
    //FORKED WORKER PROCESS TO SIEVE BSMOOTH RELATION AND PARTIAL RELATION WITH POLINOMIO FAMILIES FROM a
    _deleteLocalReports(false);          //TODO DEBUG RESET OLD SERIALIZED REPORTS
    n=n_str;
    if(argc > 1){
        n=argv[1];
    }
#endif
    /// init configuration, with argv, on un setted configuration defalut setting will be used
    CONFIGURATION* configuration= initConfiguration(n, 0, 0, 0, 0);
    Configuration=configuration;
    struct polynomial pol;
    /// init precomputation for polynomial family   todo next version computation of a coeff ... <- MASTER COMUNICATION
    PRECOMPUTES *precomputes = preComputations(configuration, &pol, false);
    if (!precomputes ){
        free(Configuration);
        exit(EXIT_FAILURE);
    }
//    printPrecomputations(precomputes,10);
    gmp_printf("\n DONE precomputation for polynomial: a: %Zd b: %Zd factorizing N: %Zd \n", pol.a.a, pol.b, *(pol.N));
    REPORTS *polynomialsReportsAggregated;
    A_COEFF *polynomialFamilyCoefficients,*polynomialFamily_a;
    polynomialFamilyCoefficients = genPolynomialFamilies_a(&POL_FAMILY_TRIES, Configuration, precomputes, configuration->a_factors_indexes_FB_families);
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
    int polFamilies_series_i = 0,sievedPolynomials=0;
    bool firstPolynomialFamily=true,polynomialFamiliesAvaible=false;

    int polFamilyOffset;
    auditExtra=false;
    polynomial_families_sieve:
    printf("concurrent multiple polynomial families sieving :%d\t on processes start\t%d %d\n", polFamilies_series_i, sievedPolynomials, POL_FAMILY_TRIES);
    sievedPolynomials+= POLYNOMIAL_FAMILIES_PER_SIEVER * POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS;
    polynomialFamiliesAvaible= sievedPolynomials < POL_FAMILY_TRIES;
    if(polFamilies_series_i>20) {auditExtra=true;printf("enabled auditing\n");}  //TODO DEBUG

    polFamilyOffset = polFamilies_series_i * POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS * POLYNOMIAL_FAMILIES_PER_SIEVER;
    for (int sieverID = 0; sieverID < POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS && polynomialFamiliesAvaible; ++sieverID) { //c hold polynomial sieved globally for a cleaner exit
        if(!fork()) {           //SIEVER TASK WILL SEARCH ITERATIVELLY ON MULTIPLE POLYNOMIAL FAMILIES
            REPORTS* sieverTaskReports;
            REPORTS sieverTaskReportsAggregated; memset(&sieverTaskReportsAggregated,0, sizeof(sieverTaskReportsAggregated));
            for (int polFamilyIndx = sieverID + polFamilyOffset,j=0; j < POLYNOMIAL_FAMILIES_PER_SIEVER && polFamilyIndx < POL_FAMILY_TRIES ; polFamilyIndx+=POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS,++j) {
                if(firstPolynomialFamily){       //one time only do first polynomial family sieve
                    firstPolynomialFamily=false;
                    if(!(sieverTaskReports=Sieve(configuration,precomputes,&pol))){
                        fprintf(stderr,"SIEVE ERROR\n");
                        exit(EXIT_FAILURE);
                    }
                    if(mergeReports(&sieverTaskReportsAggregated,sieverTaskReports,true)==EXIT_FAILURE)
                        exit(EXIT_FAILURE);
                }
                polynomialFamily_a= polynomialFamilyCoefficients + polFamilyIndx;       //get process unique polynomial family a coeff.
                if (changePolynomialFamily(precomputes, polynomialFamily_a, &pol) == EXIT_FAILURE){
                    fprintf(stderr,"CHANGE POLYNOMIAL ERROR\n");
                    exit(EXIT_FAILURE);
                }
                polynoamilFamilySize = (1 << (pol.a.a_factors_num - 1));
                for (unsigned int p = 1; p <= polynoamilFamilySize; ++p) {
                    if(!(sieverTaskReports=Sieve(configuration,precomputes,&pol))){
                        fprintf(stderr,"SIEVE ERROR\n");
                        exit(EXIT_FAILURE);
                    }
                    if(mergeReports(&sieverTaskReportsAggregated,sieverTaskReports,true)==EXIT_FAILURE) {
                        fprintf(stderr,"Merge reports error\n");
                        exit(EXIT_FAILURE);
                    }
                    nextPolynomial_b_i(&(pol.b), p, precomputes);   //change polynomial inside the family
                }
            }
            //// after all polynomial families search serialize aggregated reports -> serialization cost amortized per polynomial is low  now
            if(sieverTaskReportsAggregated.partialRelationsNum || sieverTaskReportsAggregated.relationsNum) {
                if (saveReports(&sieverTaskReportsAggregated, precomputes->factorbaseSize, false, false, &pol) == EXIT_FAILURE)
                    exit(EXIT_FAILURE);
            }
            fflush(0);fprintf(stderr,"SIEVER PROCESS founded full reports:%lu\tpartial reports:%lu\tOF %d  POLYNOMIAL FAMILIES\n",sieverTaskReportsAggregated.relationsNum,sieverTaskReportsAggregated.partialRelationsNum,POLYNOMIAL_FAMILIES_PER_SIEVER );
            exit(EXIT_SUCCESS);
        }
        if(firstPolynomialFamily)
            firstPolynomialFamily=false;
    }
    /// wait all polynomial families processes
    for (u_int j = 0; j  < POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS && polynomialFamiliesAvaible; ++j) {
        int workerPolynomialRes=0;
        wait(&workerPolynomialRes);
        if(workerPolynomialRes!=EXIT_SUCCESS)
            exit(workerPolynomialRes);
    }
    /// aggregate all reports founded until here (returned same pointer in last param) re iterate if needed
    polynomialsReportsAggregated = aggregateSieversWorkers(POLYNOMIAL_FAMILIES_CONCURRENT_SIEVERS, false,reportsAllWorker);
    if(!polynomialsReportsAggregated){
        fprintf(stderr,"ERR DURING POLYNOMIAL REPORTS AGGREGATION");
        exit(EXIT_FAILURE);
    }
    fflush(0);
    printf("\n\naggreagated reports:%lu	partialReports:%lu vs FB SIZE: %lu \n", polynomialsReportsAggregated->relationsNum,polynomialsReportsAggregated->partialRelationsNum,precomputes->factorbaseSize);
    polFamilies_series_i++;
    if((polynomialsReportsAggregated->relationsNum < precomputes->factorbaseSize - LARGE_PRIMES_SURPLUS_EXPECTED) && polynomialFamiliesAvaible){
        goto polynomial_families_sieve;
    }
    if(polynomialsReportsAggregated->relationsNum<precomputes->factorbaseSize) {
        if (pairPartialReports(reportsAllWorker) == EXIT_FAILURE)
            exit(EXIT_FAILURE);
    }
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
