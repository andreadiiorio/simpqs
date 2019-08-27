#include <unistd.h>
#include "factorizerQuick.h"

//TODO DEBUG MACRO
#define JOBS_ON_STACK                       //quick factorization module debug putting jobs on stack --> no free
#define TEST_CUNCURRENT_FACTORIZATION
#define MAX_PRIMES_TRY_PER_ITERATION 100000
#define TIMEBOXED_FACTORIZATION_ITERATION_POLLING
mpz_t LARGE_PRIME_THREASHOLD;

void appendJob(FACTORIZE_JOB_QUEUE *jobQueue, struct ArrayEntryList* newJob) {
    newJob->nextArrayEntry=NULL;        //set tail
    if(!jobQueue->queueHead) {          ////empty list
        jobQueue->queueHead=newJob;
        jobQueue->queueTail=newJob;
        return;
    }
    jobQueue->queueTail->nextArrayEntry=newJob;           //append new job
    jobQueue->queueTail=newJob;
}
void appendBlockJobs(FACTORIZE_JOB_QUEUE *jobQueue, struct ArrayEntryList* firstJobsinBlock, int block_size) {
    struct ArrayEntryList* lastBlockEntry=firstJobsinBlock+block_size-1;
    lastBlockEntry->nextArrayEntry=NULL;
    if(!jobQueue->queueHead) {          ////empty list
        jobQueue->queueHead=firstJobsinBlock;
        jobQueue->queueTail=lastBlockEntry;
        return;
    }
    jobQueue->queueTail->nextArrayEntry=firstJobsinBlock;           //append new job
    jobQueue->queueTail=lastBlockEntry;
}
struct ArrayEntry* popFirstJob(FACTORIZE_JOB_QUEUE *jobQueue) {
    ///pop first job
    struct ArrayEntryList* jobEntry=jobQueue->queueHead;
    if(!jobEntry)                                   //empty list
        return NULL;
    if(jobEntry->nextArrayEntry!=NULL){
        jobQueue->queueHead = jobEntry->nextArrayEntry ;        //eventually move the head pntr if list is not empty after pop
    }
    else {                                                      //single element (to pop) in queue ==> update head and tail to NULL
        jobQueue->queueTail=NULL;
        jobQueue->queueHead=NULL;
    }
    struct ArrayEntry* job=jobEntry->arrayEntry;
#ifndef JOBS_ON_STACK
    free(jobEntry);
#endif
    return job;
}
FACTORIZE_JOB_QUEUE* initFactorizeJobQueue(u_int64_t B,DYNAMIC_VECTOR primes){
    FACTORIZE_JOB_QUEUE* factorizeJobQueue=(FACTORIZE_JOB_QUEUE*)calloc(1,sizeof(FACTORIZE_JOB_QUEUE));
    if(!factorizeJobQueue){
        fprintf(stderr,"job queue init errd\n");
        return NULL;
    }
    if(pthread_mutex_init(&(factorizeJobQueue->mutex),NULL)==-1){
        fprintf(stderr,"mutex init errd\n");
        free(factorizeJobQueue);
        return NULL;
    }
//    if(!(factorizeJobQueue->queueHead=calloc(1,sizeof(struct ArrayEntryList)))){
//        fprintf(stderr,"job list head calloc errd\n");
//        free(factorizeJobQueue);
//        return NULL;
//    }
//    factorizeJobQueue->queueTail=factorizeJobQueue->queueHead;  //empty List have tail and head coincident
    //// other usefull links
    factorizeJobQueue->closedQueue=false;
    factorizeJobQueue->B=B;
    factorizeJobQueue->primesLink=primes.pntr;
    factorizeJobQueue->primesN=primes.vectorSize;
    return factorizeJobQueue;
}

void* FactorizeTrialDivide(void* args) {
    //trial divide factorize logic
    //lockstep factorization in a group of factorizer thread
    //each one will try a subset of primes up to Bound B
    //after a cumulation of founded factors or timeout the group will synchronize on barrier
    // the manager will aggregate all founded factors and reduce N (eliminating factors from it)
    // job done on either: 1! largePrime founded, compleate factorization over BSmooth numbers, max iteration reached (fail)

    // todo LAST FACTOR (NON B SMOOTH) identified with primality test vs dict search over primeList??


    struct factorizerArgs *factorizerArgs = (struct factorizerArgs *) args;
    int factorizerID = factorizerArgs->factorizerID;
    int exp;                        //tmp exp for factorization
    int barrierWaitReturn;
    ///timer vars for timed factorizations ierations
    STOPWATCH_INIT;
    struct timeval FACTORIZE_BARRIER_SYNC_TIMEOUT = (struct timeval) {.tv_sec=FACTORIZE_BARRIER_SYNC_TIMEOUT_SEC,
            .tv_usec=FACTORIZE_BARRIER_SYNC_TIMEOUT_USEC};
    ///vars accessed only by thread manager
    int totalFoundedFactors = 0;                //total founded factors in all iterations
    int factorizeIteration = 0;                 //iteration elapsed since factorization start
    ///output factors dynamic array malloc
    FACTOR *factors = NULL;
    int factorsDynamicSize = FACTORS_IN_SIEVE_BLOCK_REALLOC;
    if (factorizerID == FACTORIZER_MANAGER_ID) {        ///only manager allocate space for final factors
        if (!(factors = malloc(sizeof(*factors) * factorsDynamicSize))) {
            fprintf(stderr, "factors malloc fail\n");
            return (void *) EXIT_FAILURE;
        }
    }
    for (;;) {                  ///end less loop for worker thread --> they will break exitFlag setted by manager thread
// ----------------------------------  ~~~~~~  threads sync 1 ~~~~~~    ------------------------------------------------
        //// synchronize thread on N (factorize job) value <=== setted by thread manager ( last to join at the barrier)
        barrierWaitReturn = pthread_barrier_wait(factorizerArgs->barrier);
        if (barrierWaitReturn != PTHREAD_BARRIER_SERIAL_THREAD && barrierWaitReturn) {
            fprintf(stderr, "barrier error code:%d\n", barrierWaitReturn);
            return (void *) EXIT_FAILURE;
        }//NOW ALL THREAD CONCORDE WITH N VALUE (SETTED BY THREAD GROUP MANAGER)
#ifdef  DEBUG   //posix pthread barrier memory barrier check
        if(!factorizerID)
            gmp_printf("factorizing %Zd, from thread:%d \n", *(factorizerArgs->N),factorizerID);fflush(0);
#endif
        /////////// factorize group end condition <-- setted by manager on queue end and closed
        if (*(factorizerArgs->exitCondition)) {   //exit cond from manager
#ifdef DEBUG
            printf("exiting from concurrent factorization \t thread:%lu\n", pthread_self());
#endif
            break;
        }

        ////primes for trial divide worker share will fetched from precomputed list distribuiting fairly primes share by primeIndx
        u_int64_t primeIndx, prime,primeTmp;                    //index of prime for trial divide and actual prime
        u_int64_t j;                                       //counter of primes already tried to divide N
        int factorsFoundedIteration = 0;                       //num of founded factors during this factorization iteration

        /*
         * factorize tries may continue for a fiexed time or until an ammount of prime has been tried to divide N
         */
#ifdef TIMEBOXED_FACTORIZATION_ITERATION_POLLING
        STARTSTOPWACHT(STOPWATCH_START_VAR)                    //start timer for timed factorization
#endif
        ////search for N factors in thread primes subset until K factors founded or timeout
        for (primeIndx = factorizerID,j=0;  //for initializations
                    ////for loop conditions selected with a macro
                    (primeIndx < factorizerArgs->primesNum && factorsFoundedIteration < FACTORS_ASYNC_PER_ITERATION)
#ifndef TIMEBOXED_FACTORIZATION_ITERATION_POLLING
                                            && (++j)<MAX_PRIMES_TRY_PER_ITERATION;
#else
                                                                                 ;
#endif
                                                                            primeIndx += FACTORIZER_THREAD_GROUP_SIZE) {
            prime = factorizerArgs->primes[primeIndx];
#ifdef TIMEBOXED_FACTORIZATION_ITERATION_POLLING
            if (++j % TIMEOUT_CHECK_INTERVAL_FACTORIZE_WORKER == 0) {      //timed factorization iterations
                STOPWACHT_STOP(STOPWATCH_START_VAR, STOPWATCH_STOP_VAR, STOPWATCH_DELTA_VAR)
                if (timercmp(&STOPWATCH_DELTA_VAR, &FACTORIZE_BARRIER_SYNC_TIMEOUT, >)) {
                    break;                                       //factorize job has timeout
                }
            }
#endif
            ///actual trial divide for fetched prime in worker primes subset
            exp = 0; primeTmp = prime; int u = 0;
            while ((u = mpz_divisible_ui_p(*(factorizerArgs->N),primeTmp))) {
                primeTmp *= prime;                   //TODO MORE EFFICENT WAY ~~ N COPY & LOCAL DIVISION
                factorizerArgs->factorsTempStorage[factorizerID * FACTORS_ASYNC_PER_ITERATION + factorsFoundedIteration]
                                        = (FACTOR) {.factor=prime, .exp=++exp};  //set founded factors to new location
                //TODO 1! FACTOR SET WITH IF ... WHILE .exp++
            }
            if (exp)         //ANALOG TO SET FLAG IN DIVISION WHILE AND INCREMENT FACTOR INDEX HERE ! ONLY ONCE !
                factorsFoundedIteration++;

        } ///for end K factors founded or timeout expired


// ----------------------------------  ~~~~~~  threads sync 2 ~~~~~~    ------------------------------------------------
        //// synchronize factorizer thread to reduce N
        barrierWaitReturn = pthread_barrier_wait(factorizerArgs->barrier);
        if (barrierWaitReturn != PTHREAD_BARRIER_SERIAL_THREAD && barrierWaitReturn) {
            fprintf(stderr, "barrier error code:%d\n", barrierWaitReturn);
            return (void *) EXIT_FAILURE;
        }                               //NOW ALL THREAD CONCORDE ON FOUNDED FACTORS (in particular the manager)
        FACTOR factorFounded;
        /////// MANAGER AGGREAGTION
        if (factorizerID == FACTORIZER_MANAGER_ID) {

            for (int i = 0; i < FACTORIZER_THREAD_GROUP_SIZE * FACTORS_ASYNC_PER_ITERATION; i++) {
                factorFounded = factorizerArgs->factorsTempStorage[i];
                if (factorFounded.factor != 0) {    //aggregate factors locations only where workers has setted
                    factors[totalFoundedFactors++] = factorFounded;
#ifdef DEBUG
                    printf("founded factor %lu ^ %d indx:%d iteration%d\n",factorFounded.factor,factorFounded.exp,i,factorizeIteration);fflush(0);
#endif
                    //TODO CHECK EFFICIENCY DIFFERENCES WITH EXP RYSE AND DIVIDE IN 2! OP.
                    for (int k = 0; k < factorFounded.exp; ++k) {
                        ///reduce N eliminating founded factors
                        mpz_tdiv_q_ui(*(factorizerArgs->N), *(factorizerArgs->N),factorFounded.factor);
                    }
                    factorizerArgs->factorsTempStorage[i].factor=0;     //reset factor entry to void
                    REALLOC_WRAP(totalFoundedFactors, factorsDynamicSize, factors, FACTORS_IN_SIEVE_BLOCK_REALLOC)
                            return (void *) EXIT_FAILURE;
                        }}
                }
            }

            ///check if N factorization is done
            //either N fully factorized (now is 1) or N residue is a large prime (less then ~20*B and probably a prime)
            int N_has_been_fully_factorized=mpz_cmp_ui(*(factorizerArgs->N),1)==0;
//            int N_quotient_is_large_prime= mpz_cmp_ui(*(factorizerArgs->N),factorizerArgs->B*LARGE_PRIME_FACTOR_THREASHOLD ) <0 &&  mpz_probab_prime_p(*(factorizerArgs->N), PROB_PRIME_REPS);
            int N_quotient_is_likelly_prime= mpz_probab_prime_p(*(factorizerArgs->N), PROB_PRIME_REPS);
            if (N_has_been_fully_factorized || N_quotient_is_likelly_prime){
#ifdef FINAL_RESIZE
                //realloc factors to actual needed size
                if (!realloc(factors, sizeof(*factors) * totalFoundedFactors)) {
                    fprintf(stderr, "final realloc of founded factors has failed\n");
                    return (void *) EXIT_FAILURE;
                }
#endif
                ///set in place of array entry factors founded
                factorizerArgs->arrayEntryJob->factorsNum = totalFoundedFactors;
                factorizerArgs->arrayEntryJob->factors = factors;
                if (N_quotient_is_likelly_prime && mpz_cmp(*(factorizerArgs->N),LARGE_PRIME_THREASHOLD)){
                    mpz_init_set(factorizerArgs->arrayEntryJob->largePrime,*(factorizerArgs->N)); //TODO CHECK DOUBLE INIT AT LARGE PRIMEs
                }
#ifdef DEBUG
                gmp_printf("N quotient residue %Zd\tfounded %s \n",*(factorizerArgs->N),N_has_been_fully_factorized?"BSmooth number":"Large Prime mul BSmooth number");fflush(0);
#endif
                break;
            }

            ////manager exit condition
            if (++factorizeIteration == FACTORIZATION_ITERATION_MAX) {
                fprintf(stderr, "MAX FACT ITERATIONS REACHED -->FACTORIZATION FAILED\n");
                return (void *) EXIT_FAILURE;
            }
        }
#ifdef DEBUG
        if(!factorizerID)printf("new factorize iteration worker:%d\n",factorizerID);fflush(0);
#endif
    }
    return (void *) EXIT_SUCCESS;
}



void* FactorizeTrialDivideThreadGroupManager(void*factorizeJobQueueArg) {
    /*
     * factorizer thread group manager code, will get jobs from the queue until NULL jobs readed (end of queue )
     */
    FACTORIZE_JOB_QUEUE* factorizeJobQueue=(FACTORIZE_JOB_QUEUE*)factorizeJobQueueArg;
    ////init factorization stuff
    bool exitCond=false;                                        //thread worker shared exit condition
    pthread_barrier_t* barrier=malloc(sizeof(*barrier));
    if(!barrier){
        fprintf(stderr,"out of mem barrier\n");
        return (void*) EXIT_FAILURE;
    }
    if(pthread_barrier_init(barrier,NULL,FACTORIZER_THREAD_GROUP_SIZE)<0){
        fprintf(stderr,"barrier init err \n");
        return (void*) EXIT_FAILURE;
    }
    pthread_t* factorizeWorkers=calloc(FACTORIZER_THREAD_GROUP_SIZE, sizeof(*factorizeWorkers));
    if(!factorizeWorkers){
        fprintf(stderr,"workers malloc errd\n");
        return (void*)EXIT_FAILURE;
    }
    FACTOR* factorizersTempStorage=malloc(sizeof(FACTOR)*FACTORIZER_THREAD_GROUP_SIZE*FACTORS_ASYNC_PER_ITERATION);
    if(!factorizersTempStorage){
        fprintf(stderr,"temp storage malloc fail\n");
        return (void*) EXIT_FAILURE;
    }
    mpz_t N_cp;  mpz_init(N_cp);        //job N copy for in place division(only by manager)
    //set large prime threashold
    mpz_init_set_ui(LARGE_PRIME_THREASHOLD,factorizeJobQueue->B);
    mpz_mul_ui(LARGE_PRIME_THREASHOLD,LARGE_PRIME_THREASHOLD,LARGE_PRIME_THREASHOLD_FACTOR);

    //factorizer arg template --> only few change needed in copyes 5 for workers
    struct factorizerArgs factorizerArgTemplate=(struct factorizerArgs){
            .N=&N_cp,
            .B=factorizeJobQueue->B,
            .primes=factorizeJobQueue->primesLink,
            .primesNum=factorizeJobQueue->primesN,
            .barrier=barrier,
            .factorsTempStorage=factorizersTempStorage,
            .exitCondition=&exitCond
    };
    //// start factorizers worker as detached thread  --> they will be stopped manually on jobs end
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
    struct factorizerArgs factorizerArgs[FACTORIZER_THREAD_GROUP_SIZE];
    for (int i = 0; i < FACTORIZER_THREAD_GROUP_SIZE-1; i++) {  //-1 because the manager is the caller so already created
//        memcpy(factorizerArgs,&factorizerArgTemplate, sizeof(factorizerArgTemplate));   //set template arg on factorizer i arg
        factorizerArgs[i]=factorizerArgTemplate;    //TODO CHECK value set
        factorizerArgs[i].factorizerID=i;
        if (pthread_create(factorizeWorkers + i, &attr, FactorizeTrialDivide, (void *) (factorizerArgs+i)) != 0){
            fprintf(stderr, " factorizer worker creataing error at %d\n",i);
            exit(EXIT_FAILURE);
        }
//        printf("created factorizer worker %d\n",i);
    }
    //////// job fetch and workers comunication MUTEX FETCHed
    int errCode=0;
    for(;!exitCond;) {
        if ((errCode = pthread_mutex_lock(&(factorizeJobQueue->mutex)))) {
            fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
            return (void*)EXIT_FAILURE;
        }
        struct ArrayEntry *job = popFirstJob(factorizeJobQueue); //// safely pop a job from job queue
        if ((errCode = pthread_mutex_unlock(&(factorizeJobQueue->mutex)))) {
            fprintf(stderr, "mutex lock err\n err :%s\n", strerror(errCode));
            return (void*)EXIT_FAILURE;
        }
        if(!job){
            //it may means too fast job consuming || closed queue and queue ended
           if(factorizeJobQueue->closedQueue){
               ///job queue is closed setting exit condition
               printf("job queue empty\n");
               exitCond=true;
               break;
           }
           /// if queue is not closed it means that thread managers consume jobs to fastly or bad luck in thread scheduling
           /// so wait a bit and retry
           usleep(QUEUE_FILLING_POLLING_USLEEP);
           continue;
        }

        mpz_set(*(factorizerArgTemplate.N),job->element);   //set a copy of N for local computation on factorizer
        factorizerArgTemplate.arrayEntryJob=job;            //set array entry link for in place output write (only by manager)

        /////// sync comunication of newly fetched job to worker --> by mem barrier of barrier wait -> wakeup because manager entered
        //2 settings above will be visible to worker on barrier wait wake up because of memory barrier
        factorizerArgTemplate.factorizerID=FACTORIZER_MANAGER_ID;
        int factorizeJobRes=(int) FactorizeTrialDivide(&factorizerArgTemplate);               //manager will partecipate to concurrent factorization with his special ID
#ifdef DEBUG_MANAGER
        if(factorizeJobRes==EXIT_SUCCESS){
            fflush(0);gmp_printf("%Zd founded factors\t", (job->element));
            for (u_int i = 0; i < job->factorsNum; ++i) {
                printf("%lu^%d\t",job->factors[i].factor,job->factors[i].exp);
            }
            printf("\n\n-&:$-\n\n");fflush(0);
        }
#endif
    }

    //// comunicate setted exit condition to workers
    factorizerArgTemplate.factorizerID=FACTORIZER_MANAGER_ID;
    FactorizeTrialDivide(&factorizerArgTemplate);               //manager will partecipate to concurrent factorization with his special ID

    ///mem clean
    pthread_barrier_destroy(barrier); free(barrier); free(factorizeWorkers); free(factorizersTempStorage);

    return (void*)EXIT_SUCCESS;
}


/*
 * TODO test job queue append / pop sequentially
 *      test factorization of 10 random jobs +2
 *      tune experiment show:
 *              -with  little timeout of factorizer worker --> very fast founding little factors
 *              -large num of max iteration useful only with medium/large timeout
 *                  -NB LARGE TIMEOUT => SLOW DOWN LITTLE JOBS!!!
 *               -maybe good max iteration num ~ 5 (enhanced timeout polling interval)
 *                  -supposition: too mutch threading=> ???
 */


#ifdef TEST_CUNCURRENT_FACTORIZATION
int main() {
#else
int _main() {
#endif
    u_int64_t B=400000;
    DYNAMIC_VECTOR primes=ReadPrimes(PRIMES_32B_PATH,B);
    printf("----\t\treaded %lu primes\t\t----\n",primes.vectorSize);
    if(!primes.pntr){
        fprintf(stderr, "primes init fail\n");
        exit(EXIT_FAILURE);
    }
    FACTORIZE_JOB_QUEUE *factorizeJobQueue = initFactorizeJobQueue(B,primes);
    if (!factorizeJobQueue) {
        fprintf(stderr, "job queue init fail\n");
        exit(EXIT_FAILURE);
    }
    /////   random gen of job entries
    gmp_randstate_t RANDSTATE;
    gmp_randinit_default(RANDSTATE);
    mp_bitcnt_t RANDOM_SIZE = 69;
    /////   first dflt entries from mongomery polynomial
    struct ArrayEntry arrayEntry1,arrayEntry2;
    mpz_init_set_str(arrayEntry1.element,"1206727674889036177328",10);
    mpz_init_set_str(arrayEntry2.element,"504260813376616547",10);
    struct ArrayEntryList entryJob1=(struct ArrayEntryList){.arrayEntry=&arrayEntry1};
    struct ArrayEntryList entryJob2=(struct ArrayEntryList){.arrayEntry=&arrayEntry2};
    /// define block of entries
    const int JOBS_NUM_BLOCK=10;
    struct ArrayEntryList jobs[JOBS_NUM_BLOCK];
    struct ArrayEntry jobsEntries[JOBS_NUM_BLOCK];
    //gmp random initialization
    for (int i = 0; i < JOBS_NUM_BLOCK; ++i) {
        //set link in list
        jobs[i].nextArrayEntry=jobs+i+1;
        jobs[i].arrayEntry=jobsEntries+i;
        mpz_init(jobs[i].arrayEntry->element);
        mpz_urandomb(jobs[i].arrayEntry->element,RANDSTATE,RANDOM_SIZE);        //fill element with a random number
//        gmp_printf("appended job %Zd %lu \n",jobs[i].arrayEntry->element,RANDOM_SIZE);fflush(0);
    }
    jobs[JOBS_NUM_BLOCK-1].nextArrayEntry=NULL;                   //set tail of the linked list
    /// mixed append the entries block
    appendJob(factorizeJobQueue,&entryJob1);
    appendBlockJobs(factorizeJobQueue,jobs,JOBS_NUM_BLOCK);
    appendJob(factorizeJobQueue,&entryJob2);

    factorizeJobQueue->closedQueue=true;
    pthread_t manager;

    STOPWATCH_INIT;
//    memset(&STOPWATCH_DELTA_VAR,0,sizeof(STOPWATCH_DELTA_VAR));
    STARTSTOPWACHT(STOPWATCH_START_VAR)

    if (pthread_create(&manager, NULL, FactorizeTrialDivideThreadGroupManager,(void *)factorizeJobQueue)!=0){
        fprintf(stderr, " siever creataing error");
        exit(EXIT_FAILURE);
    }
    int threadRetErrCode,retvalGenThread;
    if ((threadRetErrCode = pthread_join(manager, (void **) &retvalGenThread))) {
        fprintf(stderr, "siever JOIN ERR\n %s__\n", strerror(threadRetErrCode));
        exit(EXIT_FAILURE);
    }
    printf("ret manger %d \n",retvalGenThread);
    STOPWACHT_STOP_DFLT
    return EXIT_SUCCESS;
}
