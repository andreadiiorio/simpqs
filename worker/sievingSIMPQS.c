#include <stdlib.h>
#include <stdio.h>
#include <utils/utils.h>
#include "sievingSIMPQS.h"
#include "factorizerQuick.h"
#include <math.h>
#include <pthread.h>
#include <string.h>
struct polynomial_actual* ActualPolynomial;

void sieveSubArrayForPrime(SIEVE_ARRAY_BLOCK subArray, u_int64_t subArrayLen, u_int64_t prime,
                           u_int64_t primeMultipleLocations, int64_t i);



mpz_t tmp;
#define OVERFLOW_CHECK_POLYNOMIAL   //perform overflow check
#define PROBABLY_BSMOOTH_ARRAY_BLOCK 256

#define PARTIAL_RELATIONS_BLOCK 1024

void polynomialValueQuick(int64_t j, mpz_t outputVal, struct polynomial_actual polynomial){
    //return a*j^2+2*b*j+c in outputVal
    //outputVal assumed  to be already initiated at call time
    //TODO CONTRACTED 2 MUL OPERATION OVERFLOW MAY OCCUR AT j>2^32-1
#ifdef OVERFLOW_CHECK_POLYNOMIAL
    if(LONG_MUL_OVERFLOW_P(j,j)){
        fprintf(stderr,"overflow on j^2 computation");
        exit(EXIT_FAILURE);
    }
#endif
    //2*b*j
    mpz_mul_ui(tmp, polynomial.b, j * 2);
    //a*j^2
    mpz_mul_ui(outputVal, polynomial.a, j * j);
    //a*j^2+2*b*j
    mpz_add(outputVal,outputVal,tmp);
    //a*j^2+2*b*j+c
    mpz_add(outputVal, outputVal, polynomial.c);
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

//// array block
SIEVE_ARRAY_BLOCK SieveArrayBlock;

int exit_failure=EXIT_FAILURE;

int sieveArrayBlockAllocate(struct Configuration* config, SIEVE_ARRAY_BLOCK* sieveArrayBlockPntr){
    //allocate array block of blockSized returning error to caller
    *sieveArrayBlockPntr=calloc(config->ARRAY_IN_MEMORY_MAX_SIZE, sizeof(struct ArrayEntry));
    if(!(*sieveArrayBlockPntr)){
        fprintf(stderr,"array block calloc error \n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
typedef struct sieverThreadArg{
    u_int64_t arrayShareSize;
//    u_int64_t arrayStartIndx;     //TODO USELESS IF ARRAYPNTR POINT TO FIRST SUB ARRAY ELEMENT FOR THREAD
    int64_t j_start;                //polynomialf(j) j value related to first array share
    SIEVE_ARRAY_BLOCK arrayBlockPntr;
    struct Precomputes* precomputes;
    struct Configuration* configuration;
    //TODO OTHER USEFUL POINTER DECOUPLINGIZERS

} SIEVER_THREAD_ARG;

typedef struct sieverThreadReturn{
    //TODO MATRIX ROWS FROM "direct relation
    struct ArrayEntry** largePrimesEntries;              //partial relations
    u_int64_t largePrimesNum;
} SIEVER_THREAD_RETURN;

void* siever_thread_logic(void* arg){
    /*
     * sieve a fair array portion for Bsmooth values
     *  ->add log(p) with p in FB to sub array elements divisible for p
     *  ->directly factorize array elements that are above LOG_THREASHOLD TODO**
     *  ->join to caller main thread
     *  TODO**
     *      1) subThread groups of Factorizers with dispatching with semaphore posix
     *      2) locked comunication (pipe) to SINGLE FACTORIZER THREADGROUP( threading inside it) NB NEEDED (DE)SERIALIZE mpz_t OVERHEAD
     *      3) locked pipe to Factorizer Process ( threading inside it) NB NEEDED (DE)SERIALIZE mpz_t OVERHEAD
     */
    SIEVER_THREAD_ARG sieverArg= (*(SIEVER_THREAD_ARG*) arg);
    SIEVE_ARRAY_BLOCK subArray=sieverArg.arrayBlockPntr;
    //after log sieve annotate probably Bsmooth array elements in dinamically resized array
    struct ArrayEntry** ProbablyBsmoothArrayEntries;        //likelly Bsmooth elements array pointers
    struct ArrayEntry** LargePrimeArrayEntries;             //large primes array elements pointers
    struct ArrayEntry*  bsmoothToFactPntr;
    u_int64_t k=0;                                          //counter in that array filling
    if(!(ProbablyBsmoothArrayEntries=malloc(PROBABLY_BSMOOTH_ARRAY_BLOCK * sizeof(struct ArrayEntry*)))
            || !(LargePrimeArrayEntries= malloc(PROBABLY_BSMOOTH_ARRAY_BLOCK * sizeof(struct ArrayEntry*)))){
        fprintf(stderr,"out of mem in probblyBsmoothlocations array of indexes \n");
        return (void*) &exit_failure;
    }


    /// initialize array with polynomial elements
    printf("\n --- INIT SUB ARRAY ---\n");
    u_int64_t i;    //array index
    int64_t j;      //array polynomial j for f(j)'s
    for(i = 0,j=sieverArg.j_start; i < sieverArg.arrayShareSize; i++,j++){
        //TODO LAZY ARRAY VALUES COMPUTATION --> SEE ***
        mpz_init(subArray[i].element);
        polynomialValueQuick(j,(subArray[i].element),*ActualPolynomial);    //write polynomial value f(j) in correct array location
        subArray[i].logSieveCumulative=0;
    }
    
    /// sieve for primes in factor base exploiting sieving jumps in precomputation
    printf("\n ---LOG SIEVING---\n");
    struct Precomputes* precomputations=sieverArg.precomputes;
    u_int64_t p=0;
//    u_int64_t sol1_p,sol2_p;
    u_int64_t sol1_p,sol2_p;
    for(i=0;i<sieverArg.precomputes->factorbaseSize;i++){
        p=precomputations->factorbase[i];
        //sieve jumps for prime p -->
        sol1_p=(precomputations->polPrecomputedData.sol1p[i]);
        sol2_p=(precomputations->polPrecomputedData.sol2p[i]);
        sieveSubArrayForPrime(subArray, sieverArg.arrayShareSize, p, sol1_p, sieverArg.j_start);
        sieveSubArrayForPrime(subArray, sieverArg.arrayShareSize, p, sol2_p, sieverArg.j_start);
        ///array location divisible for p will be in array entries
    }
    /// find probably Bsmooth array locations in probablyBsmoothArrayLocations
    printf("\n---identify probably Bsmooth array location---\n");
    int dynamicSizeLikellyBsmooth = PROBABLY_BSMOOTH_ARRAY_BLOCK;
    for(i = 0; i < sieverArg.arrayShareSize; i++){
        if(subArray[i].logSieveCumulative>sieverArg.configuration->LOG_SIEVING_THREASHOLD){
            if(k>dynamicSizeLikellyBsmooth){        //likelly Bsmooth indexes needed a new block reallocation
                dynamicSizeLikellyBsmooth+=PROBABLY_BSMOOTH_ARRAY_BLOCK;
                if(!realloc(ProbablyBsmoothArrayEntries, dynamicSizeLikellyBsmooth * (sizeof(struct ArrayEntry**)))){
                    fprintf(stderr,"out of mem in probblyBsmoothlocations array of indexes realloc \n");
                    free(ProbablyBsmoothArrayEntries);
                    return (void*) &exit_failure;
                }
            }
            ProbablyBsmoothArrayEntries[k++]=&(subArray[i]);         //NOTE LIKELLY B SMOOTH ARRAY ELEMENT INDEX
        }
    }
    ProbablyBsmoothArrayEntries[k++]=NULL;              //set end dynamically allocated array end

    /// try factorizing (probably) BSmooth array  elements
    printf("\n--- trial division for Bsmoothness check ---\n");
    //TODO *** ARRAY ELEMENT COMPUTATION ULTRA LAZY HERE ONLY ON LIKELLY BSMOOTH ELEMENTS!!!
    for(u_int64_t w=0; ProbablyBsmoothArrayEntries[w] != NULL; w++) {
        bsmoothToFactPntr = ProbablyBsmoothArrayEntries[w];
        factorizeTrialDivide(*bsmoothToFactPntr->element);
    }   //here all job has been correctly enqueued
    //todo lock until job queue empty => all factorized
    for(u_int64_t w=0; ProbablyBsmoothArrayEntries[w] != NULL; w++) {
        bsmoothToFactPntr = ProbablyBsmoothArrayEntries[w];
        //TODO CHECK IF FACTORIZE RESULT WRITTEN IN PLACE IS EITHER:
        // 1)useless entry => skip
        // 2) BSmooth entry => matrix row compute &store -> localMatrixRow list append -> later aggregated
        // 3) "partial" BSmooth relation=> localLargePrimes list append -> later aggregated
    }

    }
    ///return to main thread BSmooth largePrime locations easily discriminable
    SIEVER_THREAD_RETURN* sieveOutput=malloc(sizeof(SIEVER_THREAD_RETURN));
    sieveOutput->largePrimesEntries=LargePrimeArrayEntries;
    //TODO sieveOutput->matrixRows=rows
    return (void*) sieveOutput;
}

void sieveSubArrayForPrime(SIEVE_ARRAY_BLOCK subArray, u_int64_t subArrayLen, u_int64_t prime,
                           u_int64_t sol_p, int64_t j_start) {
    //sieve sub array for multiple of prime at primeMultipleLocation
    double log_p=log(prime);    //TODO CAST prime to double
    int64_t firstSieveDeltha;
    firstSieveDeltha = j_start < 0 ? (j_start + sol_p) % prime : (sol_p - j_start) % prime; //deltha firstIndx<->startSubArray
    firstSieveDeltha= firstSieveDeltha < 0 ? firstSieveDeltha * (-1) : firstSieveDeltha;    //abs deltha
    for(u_int64_t i=firstSieveDeltha;i<subArrayLen;i+=prime){
        subArray[i].logSieveCumulative+=log_p;              //add log p to log cumulative
        if(!mpz_divisible_ui_p(subArray[i].element,prime))  //TODO DEBUG CHECK
            gmp_fprintf(stderr,"INVALID LOG SIEVING PRIME:%lu ARRAY ELEMENT:%lu \n",prime,subArray[i].element);
    }
}


int Sieve(struct Configuration *config, struct Precomputes *precomputes, SIEVE_ARRAY_BLOCK sieveArrayBlock) {
    /*sieve concurrently an array of Mongomery polynomial values divided it in fixed sized blocks sieveArrayBlock
      todo sieveArrayBlock expected to be already allocated
      precomputations will be used to find quickly array location divisible per primes in factorbase
      log(p) will be added to shadow location of array A[j] for each p that divide A[j]
      array entries that will be more than a specific threashold will be factorized (hopefully quickly) TODO TIMEOUTED SKIP FACT TRY

      Large primes identified will produce "partial" relation TODO C-HASHMAP VS SUPER ENTRY LIST, SORTED FOR LARGE PRIME MATCHING
      Bsmooth values will produce relation
      relation will produce matrix rows for matrix step
    */
    int result=EXIT_SUCCESS;
    pthread_t* sieversThreads= calloc(config->SIEVING_THREAD_NUM, sizeof(*sieversThreads));
    SIEVER_THREAD_ARG* sieverThreadArgs = malloc(config->SIEVING_THREAD_NUM* sizeof(SIEVER_THREAD_ARG));
    //large primes pointers to pointers list for partial relation exploit
    int largePrimesAllCumulNum=0;
    int largePrimesAllDyncamicSize=PARTIAL_RELATIONS_BLOCK;
    struct ArrayEntry** largePrimesAll = malloc(sizeof(*largePrimesAll)*largePrimesAllDyncamicSize);        //TODO CONFRONT WITH HASHMAP?

    if(!sieversThreads || !sieverThreadArgs || !largePrimesAll){
        fprintf(stderr,"out of memory in siever thread initialization\n");
        if(sieverThreadArgs)
            free(largePrimesAll);                   //free eventually last (most significative) allocated struct
        return EXIT_FAILURE;
    }

    for(u_int64_t j= config->M * (-1); j < config->M; j+=config->ARRAY_IN_MEMORY_MAX_SIZE){           //move sieve array index j of block jumps to keep in mem only an array block of fixed size
        u_int64_t blockToAssignShare= (config->M - j) / config->SIEVING_THREAD_NUM;                 //fair share of array block to assign to each siever thread
        u_int64_t blockToAssignShareReminder= (config->M - j) % config->SIEVING_THREAD_NUM;         //last thread takes reminder too
        for(int t=0;t<config->SIEVING_THREAD_NUM;t++){
            sieverThreadArgs[t].configuration=config;
            sieverThreadArgs[t].precomputes=precomputes;
            sieverThreadArgs[t].arrayBlockPntr=sieveArrayBlock;
            sieverThreadArgs[t].j_start= j+(blockToAssignShare*t);
            sieverThreadArgs[t].arrayShareSize=blockToAssignShare;
            if(t==config->SIEVING_THREAD_NUM-1)                                                     //last thread takes the reminder too
                sieverThreadArgs[t].arrayShareSize+=blockToAssignShareReminder;

            //// start siever threads
            if( pthread_create(sieversThreads+t,NULL,siever_thread_logic, (void *) sieverThreadArgs+t)!=0){
                fprintf(stderr," siever creataing error");
                result= EXIT_FAILURE;
                goto exit;
            }
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
            if(mergeLargePrimeRefs((*retvalGenThread).largePrimesEntries, retvalGenThread->largePrimesNum, largePrimesAll,
                                &largePrimesAllCumulNum)==EXIT_FAILURE){
                result= EXIT_FAILURE;
                goto exit;
            }

                                ;     //unificate founded partial relation for later aggreagation
        }
    }


    exit:
    free(sieversThreads);
    free(sieverThreadArgs);
//    free(largePrimesAll);         //TODO LATER RE CUMULATION ???
    return result;
}

int mergeLargePrimeRefs(struct ArrayEntry **foundedPartialRelation, u_int64_t foundedPartialRelNum,
                        struct ArrayEntry **allPartialRelations, int *dinamicAllPartialRelSize,int cumulPartialRelNum) ;
    int mergeLargePrimeRefs(struct ArrayEntry **foundedPartialRelation, u_int64_t foundedPartialRelNum,
                         struct ArrayEntry **allPartialRelations, int *dinamicAllPartialRelSize,int cumulPartialRelNum) {
    //TODO CHECK PERF WITH HASHMAP
    if(cumulPartialRelNum+foundedPartialRelNum>*dinamicAllPartialRelSize){
        ///large primes all resize up needed
        *dinamicAllPartialRelSize+=PARTIAL_RELATIONS_BLOCK;
        if(!realloc(allPartialRelations,*dinamicAllPartialRelSize * sizeof(*allPartialRelations))){
            fprintf(stderr,"reallocation of partial relation cumulative fail \n");
            return EXIT_FAILURE;
        }
    }
    memcpy(allPartialRelations+cumulPartialRelNum,foundedPartialRelation,foundedPartialRelNum);//TODO foundedPartialRelNum*sizeof(...) ????
    free(foundedPartialRelation);
    return EXIT_SUCCESS;
}
