#include <stdlib.h>
#include <stdio.h>
#include <utils/utils.h>
#include "sievingSIMPQS.h"
#include "factorization/factorizerQuick.h"
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <stdbool.h>
#include <mpfr.h>
#include <unistd.h>

#define OVERFLOW_CHECK_POLYNOMIAL   //perform overflow check
//SIEVE_ARRAY_BLOCK  SieveArrayBlock;

bool auditExtra;
//#define SIEVE_JUMP_CHECK
struct polynomial* polRef;
void polynomialValueQuick(int64_t j, mpz_t outputVal, struct polynomial polynomial){
    //return a*j^2+2*b*j+c in outputVal
    //outputVal assumed  to be already initiated at call time
//#ifdef OVERFLOW_CHECK_POLYNOMIAL TODO OLD CHECK FOR DIFFERENT VERSION
//    if(LONG_MUL_OVERFLOW_P(j,j)){
//        fprintf(stderr,"overflow on j^2 computation");
//        exit(EXIT_FAILURE);
//    }
//#endif
    //a*j
    mpz_mul_si(outputVal, polynomial.a.a, j );
    //a*j+b
    mpz_add(outputVal,outputVal,polynomial.b);
    //(a*j+b)^2
    mpz_pow_ui(outputVal, outputVal, 2);

    mpz_sub(outputVal,outputVal,*(polynomial.N));

}


#define DEBUG_CHECK
void sieveSubArrayForPrime(SIEVE_ARRAY_BLOCK subArray, u_int64_t subArrayLen, u_int64_t prime, int64_t sol_p) {
    //sieve sub array for multiple of prime at primeMultipleLocation by adding log p to array entries divisible by p
    //later entries witch will be more than LOG_THREASHOLD will be tried factorized

    long int log_p=lrint(ceil(log(prime)));    //TODO PRECOMPUTE ????
//    double log_p=log(prime);    //TODO PRECOMPUTE ????
    int64_t firstSieveDeltha;
    int64_t j_start=subArray->j;
    firstSieveDeltha = ((sol_p - j_start) % (int64_t )prime); //deltha firstIndx<->startSubArray simplyfied expression
    firstSieveDeltha<0?firstSieveDeltha=prime+firstSieveDeltha:firstSieveDeltha;         //eventually handle negative sphasement with prime complement
    for(u_int64_t i=firstSieveDeltha;i<subArrayLen;i+=prime){
        subArray[i].logSieveCumulative+=log_p;              //add log p to log cumulative
#ifdef SIEVE_JUMP_CHECK
        mpz_t tmpDebug,tmp;mpz_inits(tmpDebug,tmp,NULL);
        POLYNOMIAL_VAL_COMPUTE(i + j_start, tmpDebug,polRef);
        if(!mpz_divisible_ui_p(tmpDebug,prime)) {
            gmp_fprintf(stderr, "INVALID LOG SIEVING PRIME:%lu ARRAY ELEMENT:%Zd j:%ld", prime, tmpDebug,i+j_start);fflush(0);
            pause();
        }
#endif
    }
}

///siver pthread args wrap :)
typedef struct sieverThreadArg{
    u_int64_t arrayShareSize;
//    u_int64_t arrayStartIndx;     //TODO USELESS IF ARRAYPNTR POINT TO FIRST SUB ARRAY ELEMENT FOR THREAD
    int64_t j_start;                //polynomialf(j) j value related to first array share
    SIEVE_ARRAY_BLOCK arrayBlockPntr;
    struct Precomputes* precomputes;
    struct Configuration* configuration;
    FACTORIZE_JOB_QUEUE* factorizeJobQueue;  //R-W safe queue to read and write factorization job
    //TODO OTHER USEFUL POINTER DECOUPLINGIZERS
    struct polynomial actualPol;
} SIEVER_THREAD_ARG;

void* siever_thread_logic(void* arg){
    /*
     * sieve a fair array portion for Bsmooth values passed by caller
     *  ->add log(p) with p in FB to sub array elements divisible for p
     *  ->factorize array elements that are above LOG_THREASHOLD --> writing  factorize job on queue passed
     *  ->ANY factorized entries building (partial)reports
     *
     *  ->join to caller main thread returning founded reports
     */
    SIEVER_THREAD_ARG sieverArg= *((SIEVER_THREAD_ARG*) arg);
    SIEVE_ARRAY_BLOCK subArray=sieverArg.arrayBlockPntr;
    //probablyBSmooth ref inside array in a dynamicVector
    const int BLOCK_SIZE=PROBABLY_BSMOOTH_ARRAY_BLOCK;
    struct ArrayEntry** ProbablyBsmoothArrayEntries = malloc(BLOCK_SIZE * sizeof(struct ArrayEntry*));
    u_int64_t probablyBsmoothsNum=0,dynamicSizeLikellyBsmooth = BLOCK_SIZE;
    //BSmooth*LargePrime ref inside array
    struct ArrayEntry* largePrimeEntries=malloc(BLOCK_SIZE * sizeof(*largePrimeEntries));
    u_int dynamicSizeLargePrimes=BLOCK_SIZE,foundedLargePrimes=0;
    //BSmooth entries founded in sieving
    struct ArrayEntry* bsmoothEntries=malloc(BLOCK_SIZE * sizeof(*bsmoothEntries));
    u_int dynamicSizeBsmoothEntries=BLOCK_SIZE,foundedBsmoothEntries=0;
    REPORTS* reportsOut=malloc(sizeof(*reportsOut));          //output reports produced by siever to return to worker process
    if(!largePrimeEntries || ! ProbablyBsmoothArrayEntries || !reportsOut){
        free(largePrimeEntries);free(ProbablyBsmoothArrayEntries);free(reportsOut);
        fprintf(stderr,"out of mem at siever thread initialziation\n");
        return (void*)EXIT_FAILURE;
    }
    /// initialize array block share with indexes
    struct ArrayEntry*  bsmoothToFactPntr;
    u_int64_t i;    //array index
    int64_t j;      //array polynomial j for f(j)'s

    for(i = 0,j=sieverArg.j_start; i < sieverArg.arrayShareSize; i++,j++){
        subArray[i].j=j;
        subArray[i].logSieveCumulative=0;
        /*mpz_init(subArray[i].element);
        //init exp vector of element with space to hold all primes vector exp in factorbase
        mpz_init2(subArray[i].exp_vector,sieverArg.precomputes->primes.vectorSize);
        polynomialValueQuick(j,(subArray[i].element),sieverArg.actualPol);    //write polynomial value f(j) in correct array location*/ //TODO LAZY ARRAY VALUES COMPUTATION ONLY FOR PROBABLY BSMOTH ENTRIES IDENTIFIED BY SIEVING SEE BELLOW  ***
    }

    /////// LOG SIEVE FOR PRIMES IN FACTOR BASE EXPLOITING SIEVING JUMPS IN PRECOMPUTATION
    // added log(p) to each array elements divisible by p
    polRef=&(sieverArg.actualPol);
    struct Precomputes* precomputations=sieverArg.precomputes;
    u_int64_t p=0;
    int64_t sol1_p,sol2_p;
    for(i=0;i<sieverArg.precomputes->factorbaseSize;i++){
        p=precomputations->factorbase[i];
        //sieve jumps for prime p -->
        sol1_p=(precomputations->polPrecomputedData.sol1p[i]);
        sol2_p=(precomputations->polPrecomputedData.sol2p[i]);
        ///array location divisible for p will have logSieveCumulative incremented by log p
        sieveSubArrayForPrime(subArray, sieverArg.arrayShareSize, p, sol1_p);
        sieveSubArrayForPrime(subArray, sieverArg.arrayShareSize, p, sol2_p);
    }

    // find probably Bsmooth exploiting LOG SIEVING array locations in probablyBsmoothArrayLocations factorizing polinomials values
    //probably BSmooth entries setted in an array of pointer ( next to minimun size ) dynamically reallocated  on need
    mpfr_t LOG_THREASHOLD,Ncpy; mpfr_inits(LOG_THREASHOLD,Ncpy,NULL); mpfr_set_z(Ncpy,*(sieverArg.actualPol.N),MPFR_RNDN); //TODO OPTIMIZE LOG THRESHOLD
    mpfr_sqrt(Ncpy,Ncpy,MPFR_RNDN);
    mpfr_log(LOG_THREASHOLD,Ncpy,MPFR_RNDN);
    long double LogThreashold=mpfr_get_ld(LOG_THREASHOLD,MPFR_RNDN);
    LogThreashold+=log(sieverArg.configuration->M);                         //log(sqrt(N)*M)
    LogThreashold-=DELTA_LOG_SIEVE_TOLLERATION;
    long int logThresh=lrintl(LogThreashold);
    for(i = 0; i < sieverArg.arrayShareSize; i++){
        if(subArray[i].logSieveCumulative>logThresh){
//            printf("probablyBsmoothEntry at j:%ld log cumul %Lf\n",subArray[i].j,subArray[i].logSieveCumulative);
            ProbablyBsmoothArrayEntries[probablyBsmoothsNum++]=&(subArray[i]);         //NOTE LIKELLY B SMOOTH ARRAY ELEMENT INDEX
            REALLOC_WRAP(probablyBsmoothsNum, dynamicSizeLikellyBsmooth, ProbablyBsmoothArrayEntries, PROBABLY_BSMOOTH_ARRAY_BLOCK)
                    return (void*) EXIT_FAILURE;
            }}
        }
    }
//    fprintf(stderr,"founded %lu likelly to be BSmooth array entries \tvs\t array share of size %lu \t starting from:%ld\n",probablyBsmoothsNum,sieverArg.arrayShareSize,sieverArg.j_start);fflush(0);
    ProbablyBsmoothArrayEntries[probablyBsmoothsNum]=NULL;              //set end dynamically allocated array end

    AUDIT_EXTRA
    if (!probablyBsmoothsNum) gmp_printf("no likelly BSmooth for LogSieve->useless array chunk from :%lld of %lu \t pol: %Zd - %Zd\n",sieverArg.j_start,sieverArg.arrayShareSize,sieverArg.actualPol.a.a,sieverArg.actualPol.b);
    mpfr_clears(LOG_THREASHOLD,Ncpy,NULL);
/////// FACTORIZE LIKELLY BSMOOTH ENTRIES WITH MULTIPLE CONCURRENT PARALLEL FACTORIZATIONS
    /*
     * each likelly to be Bsmooth entry will be appended in thread safe queue in a block of entries
     * each factorize thread manager will dequeue these entries 1 by 1 and try factorizing it with other factorize thread
     * Factorization, exponent vector and eventual Large prime founded will be added in place in the entry by pointer container structure ArrayEntryList
     */
    u_int64_t blockIndex=0;
    struct ArrayEntryList* factorizeJobsBlock[FACTORIZE_JOB_BLOCK_APPEND+1];    //+1 for fast next block assignement
    memset(factorizeJobsBlock,0, sizeof(*factorizeJobsBlock)*(FACTORIZE_JOB_BLOCK_APPEND+1));
    int errCode=0;
    for(u_int64_t w=0; w<probablyBsmoothsNum; w++) {
        bsmoothToFactPntr = ProbablyBsmoothArrayEntries[w];
        /*because only entry that are likelly to be BSmooth (see log sieving) will be factorized and eventually added to founded reports
         * only that one will hold actual values allocated and initiated x(j) and a*f(j) and exp vector
         */
        mpz_inits(bsmoothToFactPntr->element,bsmoothToFactPntr->x,NULL);
        //init exp vector of element with space to hold all primes vector exp in factorbase
//        mpz_init2(bsmoothToFactPntr->exp_vector,sieverArg.precomputes->primes.vectorSize+1);
        POLYNOMIAL_VAL_COMPUTE(bsmoothToFactPntr->j,(bsmoothToFactPntr->element),polRef)    //write polynomial value f(j) in correct array location
        POLYNOMIAL_VAL_COMPUTE_X(bsmoothToFactPntr->j,(bsmoothToFactPntr->x),polRef)    //write polynomial value f(j) in correct array location
        //// sync append (block of )job to job queue
#ifdef FACTORIZE_JOB_BLOCK_APPEND
        blockIndex=w%FACTORIZE_JOB_BLOCK_APPEND;
        if(blockIndex==0) {
            if (w > 0) {           //if we are not in start block iteration flush buffered block to queue
                LOCK_MUTEX( &(sieverArg.factorizeJobQueue->mutex),errCode)
                appendBlockJobs(sieverArg.factorizeJobQueue, factorizeJobsBlock[0], factorizeJobsBlock[FACTORIZE_JOB_BLOCK_APPEND-1]);
                UNLOCK_MUTEX( &(sieverArg.factorizeJobQueue->mutex),errCode)
            }
            ///malloc new blocks -> freed by Factorize Managers on dequeue
            for (int z = 0; z < FACTORIZE_JOB_BLOCK_APPEND; z++) {
                factorizeJobsBlock[z] = malloc(sizeof(*factorizeJobsBlock[z]));
            }
        }
        factorizeJobsBlock[blockIndex]->arrayEntry=bsmoothToFactPntr;
        factorizeJobsBlock[blockIndex]->nextArrayEntry=factorizeJobsBlock[blockIndex+1];        //link to next element to fill
#else
        if(!(newFactorizeJob = malloc(sizeof(*newFactorizeJob)))) {
            fprintf(stderr, "newJob malloc errd\n");
            return (void *) EXIT_FAILURE;                    // TODO CHECK FREE NEEDED WITHOUT GOTO FOR BETTER GCC OPTIMIZATION
        }
        newFactorizeJob->arrayEntry=bsmoothToFactPntr;          //set job data
        if((errCode= pthread_mutex_lock(&(sieverArg.factorizeJobQueue->mutex)))) {
            fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
            return (void *) EXIT_FAILURE;                    // TODO CHECK FREE NEEDED WITHOUT GOTO FOR BETTER GCC OPTIMIZATION
        }
        appendJob(sieverArg.factorizeJobQueue,newFactorizeJob); //// safely pop a job from job queue
        if((errCode= pthread_mutex_unlock(&(sieverArg.factorizeJobQueue->mutex)))) {
            fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
            return (void *) EXIT_FAILURE;                    // TODO CHECK FREE NEEDED WITHOUT GOTO FOR BETTER GCC OPTIMIZATION
        }
#endif
    }
#ifdef FACTORIZE_JOB_BLOCK_APPEND
    if(probablyBsmoothsNum && (blockIndex=(probablyBsmoothsNum-1)%FACTORIZE_JOB_BLOCK_APPEND!=0)) {
        LOCK_MUTEX( &(sieverArg.factorizeJobQueue->mutex),errCode)
        appendBlockJobs(sieverArg.factorizeJobQueue, factorizeJobsBlock[0], factorizeJobsBlock[blockIndex]);
        UNLOCK_MUTEX( &(sieverArg.factorizeJobQueue->mutex),errCode)
    }
#endif
    //HERE ALL JOB HAS BEEN CORRECTLY ENQUEUED
    //mark as terminated producer's queue end and wait until consumers thread will finish to factorize jobs (posix condvar wait)
    LOCK_MUTEX(&(sieverArg.factorizeJobQueue->mutex),errCode)
    sieverArg.factorizeJobQueue->producersEnded++;       //end of the job production for this thread
//    fprintf(stderr ,"waiting queue job ended from siever %lu\t condvar at %p mutex at %p\n",pthread_self(),&(sieverArg.factorizeJobQueue->emptyAndClosedQueue),&(sieverArg.factorizeJobQueue->mutex));fflush(0);
    if(pthread_cond_wait(&(sieverArg.factorizeJobQueue->emptyAndClosedQueue),&(sieverArg.factorizeJobQueue->mutex))!=0) {
        fprintf(stderr, "invalid cond wait on conditional var\n");
        return (void *) EXIT_FAILURE;
    }
    UNLOCK_MUTEX(&(sieverArg.factorizeJobQueue->mutex),errCode)

    //////   CHECK JOB RESULTS AND BUILD REPORTS
    int uselessEntry=0;
    struct ArrayEntry* entry;
    for(u_int64_t w=0; w<probablyBsmoothsNum; w++) {
        entry=ProbablyBsmoothArrayEntries[w];
        if(entry->factors) {
            if(entry->largePrime->_mp_d) {                //true if large prime is set TODO CALLOC AND CLASSIC NON ZER CHECK
//                printf("partial relation\n");
                REALLOC_WRAP(foundedLargePrimes, dynamicSizeLargePrimes, largePrimeEntries, BLOCK_SIZE)
                        fprintf(stderr, "realloc failed on large primes ref resize\n");
                }}
                arrayEntryCopy(&(largePrimeEntries[foundedLargePrimes++]),entry);
            } else {
//                printf("relation\n");
                REALLOC_WRAP(foundedBsmoothEntries, dynamicSizeBsmoothEntries, bsmoothEntries, BLOCK_SIZE)
                        fprintf(stderr, "realloc failed on bsmoothentries ref resize\n");
                }}
                arrayEntryCopy(&(bsmoothEntries[foundedBsmoothEntries++]), entry);
            }
        } else
            uselessEntry++;

    }

    //TODO DEBUG
    AUDIT_EXTRA        if(auditExtra)     printf("founded useless rel:%d\trelations:%d\t,partialRelation:%d vs likellYBSmooth :%d\n", uselessEntry, foundedBsmoothEntries, foundedLargePrimes,probablyBsmoothsNum);
    // 1)useless entry => skip
    // 2) BSmooth entry => matrix row compute &store -> localMatrixRow list append -> later aggregated
    // 3) "partial" BSmooth relation=> localLargePrimes list append -> later aggregated
#ifdef FINAL_RESIZE
    if(  foundedBsmoothEntries && !(bsmoothEntries=realloc(bsmoothEntries, foundedBsmoothEntries * sizeof(*bsmoothEntries)))) {
        fprintf(stderr, "last resize on large primes failed\n");
        return (void *) EXIT_FAILURE;
    }
    if(foundedLargePrimes && !(largePrimeEntries=realloc(largePrimeEntries, foundedLargePrimes * sizeof(*largePrimeEntries)))) {
        fprintf(stderr, "last resize on large primes failed\n");
        return (void *) EXIT_FAILURE;
    }
#endif

    ///return to main founded (partial) reports
    reportsOut->bsmoothEntries=bsmoothEntries;
    reportsOut->relationsNum=foundedBsmoothEntries;
    reportsOut->largePrimesEntries=largePrimeEntries;
    reportsOut->partialRelationsNum=foundedLargePrimes;

    //deallocate no more usefull suff
    free(ProbablyBsmoothArrayEntries);
    return (void*) reportsOut;
}

REPORTS* Sieve(struct Configuration *config, struct Precomputes *precomputes, struct polynomial* actualPol) {
    /*sieve concurrently an array of 1 Mongomery polynomial values divided it in fixed sized blocks sieveArrayBlock
      precomputations will be used to find quickly array location divisible per primes in factorbase
      log(p) will be added to shadow location of array A[j] for each p that divide A[j]
      array entries that will be more than a specific threashold will be factorized (hopefully quickly) concurrently with multiple thread group reading jobs from a locked Queue

      Large primes identified will produce "partial" relation TODO C-HASHMAP VS SUPER ENTRY LIST, SORTED FOR LARGE PRIME MATCHING
      Bsmooth values will produce relation
      relation will produce matrix rows for matrix step
    */
    if(auditExtra)     gmp_printf("start Sieving process on polynomial %Zd - %Zd \n",actualPol->a.a,actualPol->b);fflush(0);
    int result=EXIT_SUCCESS;
    pthread_t *factorizersThreadManagers=NULL;
    pthread_t* sieversThreads= calloc(config->SIEVING_THREAD_NUM, sizeof(*sieversThreads));
    SIEVER_THREAD_ARG* sieverThreadArgs = malloc(config->SIEVING_THREAD_NUM* sizeof(*sieverThreadArgs));
    //large primes pointers to pointers list for partial relation exploit

    REPORTS* polynomialReports=calloc(1,sizeof(*polynomialReports));
    u_int64_t arrayInMemSize=MIN(2 * config->M, config->ARRAY_IN_MEMORY_MAX_SIZE);
    SIEVE_ARRAY_BLOCK sieveArrayBlock=calloc(1, arrayInMemSize* sizeof(*sieveArrayBlock));
    if(!sieversThreads || !sieverThreadArgs || !polynomialReports || !sieveArrayBlock){
        fprintf(stderr,"out of memory in siever thread initialization\n");
        result=EXIT_FAILURE;
        goto exit;
    }
    FACTORIZE_JOB_QUEUE* factorizeJobQueue= initFactorizeJobQueue(config->B, precomputes, config->SIEVING_THREAD_NUM, NUM_FACTORIZER_GROUPS);
    if(!factorizeJobQueue) {
        result=EXIT_FAILURE;
        goto exit;
    }
    /////iterate trough j in [-M,M] with a block per time so a fixed and configurable ammount of MEM will be used for the array

    for(int64_t j= config->M * (-1),upLimit=config->M;  j <= upLimit; j+=config->ARRAY_IN_MEMORY_MAX_SIZE) { //move sieve array index j of block jumps to keep in mem only an array block of fixed size
#ifdef VERBOSE
        printf("---\tnew array block loading of %lu elements starting from index j:%ld\n",config->ARRAY_IN_MEMORY_MAX_SIZE, j);
#endif
        u_int64_t blockToAssignShare = (arrayInMemSize) / config->SIEVING_THREAD_NUM;//fair share of array block to assign to each siever thread
        u_int64_t blockToAssignShareReminder = (arrayInMemSize) % config->SIEVING_THREAD_NUM;//last thread takes reminder too
        //// division of array block in subset to assign to siever threads
        for (int t = 0; t < config->SIEVING_THREAD_NUM; t++) {
            sieverThreadArgs[t].arrayBlockPntr =sieveArrayBlock + (blockToAssignShare * t);  //point to siever's array block subset
            sieverThreadArgs[t].j_start = j++ + (blockToAssignShare * t);
            sieverThreadArgs[t].arrayShareSize = blockToAssignShare;
            sieverThreadArgs[t].factorizeJobQueue = factorizeJobQueue;
            sieverThreadArgs[t].configuration = config;
            sieverThreadArgs[t].precomputes = precomputes;
            sieverThreadArgs[t].actualPol = *actualPol;

            if (t == config->SIEVING_THREAD_NUM - 1)    //last thread takes the reminder too
                sieverThreadArgs[t].arrayShareSize += blockToAssignShareReminder;
            //// start siever threads
            if (pthread_create(sieversThreads + t, NULL, siever_thread_logic, (void *) (sieverThreadArgs + t)) != 0) {
                fprintf(stderr, " siever creataing error");
                result = EXIT_FAILURE;
                goto exit;
            }
#ifdef VERBOSE
            printf(">>>\tcreated thread with array start address:%p\t j in [%ld,%ld]\n", sieverThreadArgs[t].arrayBlockPntr,sieverThreadArgs[t].j_start,sieverThreadArgs[t].j_start+sieverThreadArgs[t].arrayShareSize);fflush(0);
#endif
        }
        ///start factorizers thread manager --> subsequentially workers :=)
        factorizersThreadManagers = StartFactorizerThreadGroups(factorizeJobQueue, NUM_FACTORIZER_GROUPS);
        if (!factorizersThreadManagers) {
            fprintf(stderr, "factorize thread group init faild\n");
            result = EXIT_FAILURE;
            goto exit;
        }
        /// join factorizer threads
        if((result=JoinFactorizerThreadGroups(factorizersThreadManagers,NUM_FACTORIZER_GROUPS))!=EXIT_SUCCESS){
            goto exit;
        }   //now siever will see updated array with exp vector of BSmooth entry and LargePrime one

        //// join sievers threads
        REPORTS *reportsFounded=NULL;
        for (int t = 0; t < config->SIEVING_THREAD_NUM; t++) {
            int threadRetErrCode;
#ifdef VERBOSE
            printf("joining siever thread\t:%d \t %lu \n",t,sieversThreads[t]);fflush(0);
#endif
            if ((threadRetErrCode = pthread_join(sieversThreads[t], (void **) &reportsFounded))) {
                fprintf(stderr, "siever JOIN ERR\n %s__\n", strerror(threadRetErrCode));
                result = EXIT_FAILURE;
                goto exit;
            }
//            print_reports(reportsFounded, precomputes->primes.vectorSize);
            if (mergeReports(polynomialReports, (REPORTS *) reportsFounded, true) == EXIT_FAILURE) {
                fprintf(stderr, "MERGE POLYNOMIAL REPORTS ERR\n");
                result = EXIT_FAILURE;
                goto exit;
            }
        }
        resetFactorizeJobQueue(factorizeJobQueue, config->SIEVING_THREAD_NUM, NUM_FACTORIZER_GROUPS);
    }

    exit:
//    printf("Done with polynomial\n");
    free(sieversThreads);free(sieverThreadArgs);free(sieveArrayBlock);
    pthread_cond_destroy(&factorizeJobQueue->emptyAndClosedQueue);
    pthread_mutex_destroy(&factorizeJobQueue->mutex);
    free(factorizeJobQueue);
    if(result==EXIT_FAILURE){
        fprintf(stderr,"error on polynomial\n");
        return NULL;
    }
    return polynomialReports;
}


