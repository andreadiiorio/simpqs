#include <stdlib.h>
#include <stdio.h>
#include <utils/utils.h>
#include "sievingSIMPQS.h"
#include "factorization/factorizerQuick.h"
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <factorization/factorizerQuick.h>
#include <stdbool.h>
#include <mpfr.h>

#define OVERFLOW_CHECK_POLYNOMIAL   //perform overflow check
#define VERBOSE
#define DELTA_LOG_SIEVE_TOLLERATION 7       //tolleration constant for log threashold comparison for factorize job mark
mpz_t N,tmp;
//SIEVE_ARRAY_BLOCK  SieveArrayBlock;

struct polynomial_actual* ActualPolynomial;

//#undef DEBUG
void sieveSubArrayForPrime(SIEVE_ARRAY_BLOCK subArray, u_int64_t subArrayLen, u_int64_t prime, int64_t sol_p) {
    //sieve sub array for multiple of prime at primeMultipleLocation by adding log p to array entries divisible by p
    //later entries witch will be more than LOG_THREASHOLD will be tried factorized

    double log_p=log(prime);    //TODO PRECOMPUTE ????
    int64_t firstSieveDeltha;
    int64_t j_start=subArray->j;
//    firstSieveDeltha = j_start < 0 ? firstSieveDeltha=((sol_p-j_start) % prime): firstSieveDeltha=((sol_p - j_start) % prime); //deltha firstIndx<->startSubArray
    firstSieveDeltha = ((sol_p - j_start) % (int64_t )prime); //deltha firstIndx<->startSubArray simplyfied expression
//    sol_p < j_start?firstSieveDeltha=prime+firstSieveDeltha:firstSieveDeltha;          //reverse first index if from sol_1 -> positive jumps to j_start (
    firstSieveDeltha<0?firstSieveDeltha=prime-firstSieveDeltha:firstSieveDeltha;         //eventually handle negative sphasement with prime complement
    for(u_int64_t i=firstSieveDeltha;i<subArrayLen;i+=prime){
        subArray[i].logSieveCumulative+=log_p;              //add log p to log cumulative
#ifdef DEBUG_CHECK
        if(!mpz_divisible_ui_p(subArray[i].element,prime)) {
            gmp_fprintf(stderr, "INVALID LOG SIEVING PRIME:%lu ARRAY ELEMENT:%Zd j:%ld\n N:%Zd", prime, subArray[i].element,subArray[i].j,N);
            exit(EXIT_FAILURE);
        }
#endif
    }
}
int mergeLargePrimeRefs(struct ArrayEntry **foundedPartialRelation, u_int64_t foundedPartialRelNum,struct ArrayEntry **allPartialRelations, int *dinamicAllPartialRelSize,int cumulPartialRelNum) ;

void polynomialValueQuick(int64_t j, mpz_t outputVal, struct polynomial_actual polynomial){
    //return a*j^2+2*b*j+c in outputVal
    //outputVal assumed  to be already initiated at call time
//#ifdef OVERFLOW_CHECK_POLYNOMIAL TODO OLD CHECK FOR DIFFERENT VERSION
//    if(LONG_MUL_OVERFLOW_P(j,j)){
//        fprintf(stderr,"overflow on j^2 computation");
//        exit(EXIT_FAILURE);
//    }
//#endif
    //a*j
    mpz_mul_si(outputVal, polynomial.a, j );
    //a*j+b
    mpz_add(outputVal,outputVal,polynomial.b);
    //(a*j+b)^2
    mpz_pow_ui(outputVal, outputVal, 2);

    mpz_sub(outputVal,outputVal,N);

}
void polynomialValue(u_int64_t j,mpz_t outputVal,struct polynomial_actual polynomial){
    //return a*j^2+2*b*j+c in outputVal
    //outputVal assumed  to be already initiated at call time
    //2*b*j
    mpz_mul_ui(tmp,polynomial.b,j);
    mpz_mul_ui(tmp,polynomial.b,2);
    //a*j^2
    mpz_mul_ui(outputVal,polynomial.a,j);
    mpz_mul_ui(outputVal,outputVal,j);
    //a*j^2+2*b*j
    mpz_add(outputVal,outputVal,tmp);
    //a*j^2+2*b*j+c
    mpz_add(outputVal,outputVal,polynomial.c);
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

} SIEVER_THREAD_ARG;
typedef struct sieverThreadReturn{
    //TODO MATRIX ROWS FROM "direct relation
    DYNAMIC_VECTOR largePrimesEntries;              //partial relations
    u_int64_t largePrimesNum;
} SIEVER_THREAD_RETURN;

void* siever_thread_logic(void* arg){
    /*
     * sieve a fair array portion for Bsmooth values
     *  ->add log(p) with p in FB to sub array elements divisible for p
     *  ->directly factorize array elements that are above LOG_THREASHOLD --> write factorize job on queue passed
     *  ->join to caller main thread
     */
    SIEVER_THREAD_ARG sieverArg= *((SIEVER_THREAD_ARG*) arg);
    SIEVE_ARRAY_BLOCK subArray=sieverArg.arrayBlockPntr;
    ///probablyBSmooth ref inside array in a dynamicVector
    const int SIEVE_OUTS_BLOCK_SIZE=PROBABLY_BSMOOTH_ARRAY_BLOCK;
    struct ArrayEntry** ProbablyBsmoothArrayEntries = malloc(PROBABLY_BSMOOTH_ARRAY_BLOCK * sizeof(struct ArrayEntry*));
    u_int64_t probablyBsmoothsNum=0; int dynamicSizeLikellyBsmooth = SIEVE_OUTS_BLOCK_SIZE;
    ///BSmooth*LargePrime ref inside array --> merged at a next step
    struct ArrayEntry** LargePrimeArrayEntries=malloc(PROBABLY_BSMOOTH_ARRAY_BLOCK * sizeof(struct ArrayEntry*));
    int dynamicSizeLargePrimes=SIEVE_OUTS_BLOCK_SIZE,foundedLargePrimes=0;
    if(!LargePrimeArrayEntries || ! ProbablyBsmoothArrayEntries){
        free(LargePrimeArrayEntries);free(ProbablyBsmoothArrayEntries);
        fprintf(stderr,"out of mem at siever thread initialziation\n");
        return (void*)EXIT_FAILURE;
    }
    struct ArrayEntry*  bsmoothToFactPntr;

    /// initialize array block share with polynomial elements
    u_int64_t i;    //array index
    int64_t j;      //array polynomial j for f(j)'s

    //TODO MEM PASSEDDEBUG BLOCK
//    char myArrayStartPntr[44] ;
//    snprintf(myArrayStartPntr,44,"%p",subArray);
//    printf("MY_START PNTR:%s\n",myArrayStartPntr);

    for(i = 0,j=sieverArg.j_start; i < sieverArg.arrayShareSize; i++,j++){
        //TODO LAZY ARRAY VALUES COMPUTATION BELLOW --> SEE ***
        subArray[i].j=j;
        mpz_init(subArray[i].element);
        polynomialValueQuick(j,(subArray[i].element),*ActualPolynomial);    //write polynomial value f(j) in correct array location
        subArray[i].logSieveCumulative=0;
    }
#ifdef VERBOSE
    printf("\n --- \tINITIATED SUB ARRAY at %lu\t --- \n",pthread_self());
#endif

    /// sieve for primes in factor base exploiting sieving jumps in precomputation
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
    /// find probably Bsmooth array locations in probablyBsmoothArrayLocations
    //TODO OPTIMIZE LOG THRESHOLD..---...---..---..---..---...---...---
    mpfr_t LOG_THREASHOLD,Ncpy; mpfr_inits(LOG_THREASHOLD,Ncpy,NULL); mpfr_set_z(Ncpy,N,MPFR_RNDN);
    mpfr_sqrt(Ncpy,Ncpy,MPFR_RNDN);
    mpfr_log(LOG_THREASHOLD,Ncpy,MPFR_RNDN);
    long double LogThreashold=mpfr_get_ld(LOG_THREASHOLD,MPFR_RNDN);
    LogThreashold+=log(sieverArg.configuration->M);
    LogThreashold-=DELTA_LOG_SIEVE_TOLLERATION;
    fflush(0);

#ifdef VERBOSE
    printf("\n ---\tLOG SIEVING START at %lu\t---\n "
           "\t\tlogSieve threshold :%Lf reduced of %d\n",pthread_self() ,LogThreashold,DELTA_LOG_SIEVE_TOLLERATION);
#endif
    for(i = 0; i < sieverArg.arrayShareSize; i++){
        if(subArray[i].logSieveCumulative>LogThreashold){
//            printf("probablyBsmoothEntry at j:%ld log cumul %Lf\n",subArray[i].j,subArray[i].logSieveCumulative);
            ProbablyBsmoothArrayEntries[probablyBsmoothsNum++]=&(subArray[i]);         //NOTE LIKELLY B SMOOTH ARRAY ELEMENT INDEX
            REALLOC_WRAP(probablyBsmoothsNum, dynamicSizeLikellyBsmooth, ProbablyBsmoothArrayEntries, PROBABLY_BSMOOTH_ARRAY_BLOCK)
                    return (void*) EXIT_FAILURE;
                }}
        }
//        else
//            fprintf(stderr,"bellow log threshold at j:%ld of %Lf\n",subArray[i].j,LogThreashold -subArray[i].logSieveCumulative);
        fflush(0);
    }
#ifdef VERBOSE
    fprintf(stderr,"founded %lu likelly to be BSmooth array entries\n",probablyBsmoothsNum);fflush(0);
#endif
    ProbablyBsmoothArrayEntries[probablyBsmoothsNum++]=NULL;              //set end dynamically allocated array end

    /// try factorizing (probably) BSmooth array  elements
    u_int64_t blockIndex=0;
    struct ArrayEntryList* factorizeJobsBlock[FACTORIZE_JOB_BLOCK_APPEND+1];    //+1 for fast next block assignement
    int errCode=0;
    for(u_int64_t w=0; ProbablyBsmoothArrayEntries[w] != NULL; w++) {
        //TODO *** ARRAY ELEMENT COMPUTATION ULTRA LAZY HERE ONLY ON LIKELLY BSMOOTH ELEMENTS!!!
        bsmoothToFactPntr = ProbablyBsmoothArrayEntries[w];
        //// sync append job to job queue
#ifdef FACTORIZE_JOB_BLOCK_APPEND
        blockIndex=w%FACTORIZE_JOB_BLOCK_APPEND;
        if(blockIndex==0) {
            if (w > 0) {           //if we are not in start block iteration flush buffered block to queue
                if ((errCode = pthread_mutex_lock(&(sieverArg.factorizeJobQueue->mutex)))) {
                    fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
                    return (void *) EXIT_FAILURE;                    // TODO CHECK FREE NEEDED WITHOUT GOTO FOR BETTER GCC OPTIMIZATION
                }
                appendBlockJobs(sieverArg.factorizeJobQueue, factorizeJobsBlock[0], factorizeJobsBlock[FACTORIZE_JOB_BLOCK_APPEND-1]);
                if ((errCode = pthread_mutex_unlock(&(sieverArg.factorizeJobQueue->mutex)))) {
                    fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
                    return (void *) EXIT_FAILURE;                    // TODO CHECK FREE NEEDED WITHOUT GOTO FOR BETTER GCC OPTIMIZATION
                }
            }
            ///malloc new block
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
    }   //here all job has been correctly enqueued
#ifdef FACTORIZE_JOB_BLOCK_APPEND
        //TODO FLUSH LAST BLOCK OR IGNORE??
#endif
    ///set queue end for this siever thread
    if((errCode= pthread_mutex_lock(&(sieverArg.factorizeJobQueue->mutex)))) {
        fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
        free(LargePrimeArrayEntries);free(ProbablyBsmoothArrayEntries);
        return (void *) EXIT_FAILURE;                    // TODO CHECK FREE NEEDED WITHOUT GOTO FOR BETTER GCC OPTIMIZATION
    }
    sieverArg.factorizeJobQueue->producersEnded++;       //end of the job production for this thread
    if((errCode= pthread_mutex_unlock(&(sieverArg.factorizeJobQueue->mutex)))) {
        fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
        free(LargePrimeArrayEntries);free(ProbablyBsmoothArrayEntries);
        return (void *) EXIT_FAILURE;                    // TODO CHECK FREE NEEDED WITHOUT GOTO FOR BETTER GCC OPTIMIZATION
    }

    //// wait until factorizer thread finish to factorize jobs
    LOCK_MUTEX(&(sieverArg.factorizeJobQueue->mutex),errCode)
    fprintf(stderr ,"waiting queue job ended from siever %lu\t condvar at %p mutex at %p\n",
            pthread_self(),&(sieverArg.factorizeJobQueue->emptyAndClosedQueue),&(sieverArg.factorizeJobQueue->mutex));fflush(0);
    if(pthread_cond_wait(&(sieverArg.factorizeJobQueue->emptyAndClosedQueue),&(sieverArg.factorizeJobQueue->mutex))!=0) {
        fprintf(stderr, "invalid cond wait on conditional var\n");
        return (void *) EXIT_FAILURE;
    }
    UNLOCK_MUTEX(&(sieverArg.factorizeJobQueue->mutex),errCode)
    //////         check job results
    struct ArrayEntry* entry;
    for(u_int64_t w=0; ProbablyBsmoothArrayEntries[w] != NULL; w++) {
        entry=ProbablyBsmoothArrayEntries[w];
        if(entry->factors) {
            if(entry->largePrime->_mp_d) {                //true if large prime is set
                printf("L\n");
                LargePrimeArrayEntries[foundedLargePrimes++] = entry;
                REALLOC_WRAP(foundedLargePrimes, dynamicSizeLargePrimes, LargePrimeArrayEntries,
                             PROBABLY_BSMOOTH_ARRAY_BLOCK)
                        fprintf(stderr, "realloc failed on large primes ref resize\n");
                    }}
            }
            else {
                printf("B\n");  //TODO SET IN MATRIX ROW
            }
        }
        // 1)useless entry => skip
        // 2) BSmooth entry => matrix row compute &store -> localMatrixRow list append -> later aggregated
        // 3) "partial" BSmooth relation=> localLargePrimes list append -> later aggregated
    }
#ifdef FINAL_RESIZE
    if(!(LargePrimeArrayEntries=realloc(LargePrimeArrayEntries,foundedLargePrimes* sizeof(LargePrimeArrayEntries)))&&foundedLargePrimes) {
        fprintf(stderr, "last resize on large primes failed\n");
        return (void *) EXIT_FAILURE;
    }
#endif
    ///return to main thread BSmooth largePrime locations easily discriminable
    SIEVER_THREAD_RETURN* sieveOutput=(SIEVER_THREAD_RETURN*) malloc(sizeof(SIEVER_THREAD_RETURN));
    sieveOutput->largePrimesEntries=(DYNAMIC_VECTOR){.pntr =LargePrimeArrayEntries,.vectorSize=dynamicSizeLargePrimes};
    //TODO sieveOutput->matrixRows=rows
    free(ProbablyBsmoothArrayEntries);
    return (void*) sieveOutput;
}



int Sieve(struct Configuration *config, struct Precomputes *precomputes, SIEVE_ARRAY_BLOCK sieveArrayBlock) {
    /*sieve concurrently an array of Mongomery polynomial values divided it in fixed sized blocks sieveArrayBlock
      sieveArrayBlock expected to be already allocated
      precomputations will be used to find quickly array location divisible per primes in factorbase
      log(p) will be added to shadow location of array A[j] for each p that divide A[j]
      array entries that will be more than a specific threashold will be factorized (hopefully quickly) TODO TIMEOUTED SKIP FACT TRY

      Large primes identified will produce "partial" relation TODO C-HASHMAP VS SUPER ENTRY LIST, SORTED FOR LARGE PRIME MATCHING
      Bsmooth values will produce relation
      relation will produce matrix rows for matrix step
    */
    int result=EXIT_SUCCESS;
    pthread_t *factorizersThreadManagers;
    pthread_t* sieversThreads= calloc(config->SIEVING_THREAD_NUM, sizeof(*sieversThreads));
    SIEVER_THREAD_ARG* sieverThreadArgs = malloc(config->SIEVING_THREAD_NUM* sizeof(SIEVER_THREAD_ARG));
    //large primes pointers to pointers list for partial relation exploit
    int largePrimesAllDyncamicSize=PARTIAL_RELATIONS_BLOCK;
    struct ArrayEntry** largePrimesAll = malloc(sizeof(*largePrimesAll)*largePrimesAllDyncamicSize);        //TODO CONFRONT WITH HASHMAP?

    if(!sieversThreads || !sieverThreadArgs || !largePrimesAll){
        fprintf(stderr,"out of memory in siever thread initialization\n");
        result=EXIT_FAILURE;
        goto exit;
    }
    FACTORIZE_JOB_QUEUE* factorizeJobQueue=initFactorizeJobQueue(config->B,precomputes,config->SIEVING_THREAD_NUM);
    if(!factorizeJobQueue) {
        result=EXIT_FAILURE;
        goto exit;
    }
    ///// sieving array block in MEM for --> new array block loded in mem selecting a range of j
    for(int64_t j= config->M * (-1);  j < config->M; j+=config->ARRAY_IN_MEMORY_MAX_SIZE){           //move sieve array index j of block jumps to keep in mem only an array block of fixed size
        printf("---\tnew array block loading of %lu elements starting from index j:%ld\n",config->ARRAY_IN_MEMORY_MAX_SIZE,j);
        u_int64_t blockToAssignShare= (config->ARRAY_IN_MEMORY_MAX_SIZE) / config ->SIEVING_THREAD_NUM;                 //fair share of array block to assign to each siever thread
        u_int64_t blockToAssignShareReminder= (config->ARRAY_IN_MEMORY_MAX_SIZE) % config->SIEVING_THREAD_NUM;         //last thread takes reminder too
        factorizeJobQueue->producersEnded=0;
        factorizeJobQueue->endedQueue=false;
        //// division of array block in subset to assign to siever threads
        for(int t=0;t<config->SIEVING_THREAD_NUM;t++){
            sieverThreadArgs[t].arrayBlockPntr=sieveArrayBlock+(blockToAssignShare*t);  //will point to siever's array block subset
            sieverThreadArgs[t].j_start= j+(blockToAssignShare*t);
            sieverThreadArgs[t].arrayShareSize=blockToAssignShare;
            sieverThreadArgs[t].factorizeJobQueue=factorizeJobQueue;
            sieverThreadArgs[t].configuration=config;
            sieverThreadArgs[t].precomputes=precomputes;
            if(t==config->SIEVING_THREAD_NUM-1)                                                     //last thread takes the reminder too
                sieverThreadArgs[t].arrayShareSize+=blockToAssignShareReminder;
            printf(">>>\tcreated thread with array start address:%p\n",sieverThreadArgs[t].arrayBlockPntr);fflush(0);
            //// start siever threads
            if( pthread_create(sieversThreads+t,NULL,siever_thread_logic, (void *) (sieverThreadArgs+t))!=0){
                fprintf(stderr," siever creataing error");
                result= EXIT_FAILURE;
                goto exit;
            }
        }
        ///start factorizers thread manager --> subsequentially workers :=)
        factorizersThreadManagers = StartFactorizerThreadGroups(factorizeJobQueue, NUM_FACTORIZER_GROUPS);
        if(!factorizersThreadManagers){
            fprintf(stderr,"factorize thread group init faild\n");
            result=EXIT_FAILURE;
            goto exit;
        }

        //// join sievers threads
        for(int t=0;t<config->SIEVING_THREAD_NUM;t++){
            SIEVER_THREAD_RETURN* retvalGenThread;
            int threadRetErrCode;
            //joining...
            if((threadRetErrCode= pthread_join(sieversThreads[t],(void**) &retvalGenThread))) {
                fprintf(stderr, "siever JOIN ERR\n %s__\n", strerror(threadRetErrCode));
                result= EXIT_FAILURE;
                goto exit;
            }
            printf("<<<\tjoined siever thread with retval at %p\n",retvalGenThread);fflush(0);
//            if(mergeLargePrimeRefs((*retvalGenThread).largePrimesEntries, retvalGenThread->largePrimesNum, largePrimesAll,
//                                &largePrimesAllCumulNum)==EXIT_FAILURE){
//                result= EXIT_FAILURE;
//                goto exit;
//            }

        }
    }


    exit:
    free(sieversThreads);
    free(sieverThreadArgs);
    free(factorizersThreadManagers);
//    free(largePrimesAll);         //TODO LATER RE CUMULATION ???
    return result;
}

    int mergeLargePrimeRefs(struct ArrayEntry **foundedPartialRelation, u_int64_t foundedPartialRelNum,
                         struct ArrayEntry **allPartialRelations, int *dinamicAllPartialRelSize,int cumulPartialRelNum) {
    //TODO CHECK PERF WITH HASHMAP
    REALLOC_WRAP(cumulPartialRelNum+foundedPartialRelNum,(*dinamicAllPartialRelSize),allPartialRelations,PARTIAL_RELATIONS_BLOCK)
            return EXIT_FAILURE;
        }
    }
    memcpy(allPartialRelations+cumulPartialRelNum,foundedPartialRelation,foundedPartialRelNum);//TODO foundedPartialRelNum*sizeof(...) ????
    free(foundedPartialRelation);
    return EXIT_SUCCESS;
}

