#include <unistd.h>
#include "factorizerQuick.h"

//#define TEST_CUNCURRENT_FACTORIZATION //TODO DEBUG MACRO
#ifdef TEST_CUNCURRENT_FACTORIZATION
#define JOBS_ON_STACK                       //quick factorization module debug putting jobs on stack --> no free
#endif
#define FACTORIZE_END_POLLING_USEC 333
//#define DEBUG_MANAGER


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
void appendBlockJobs(FACTORIZE_JOB_QUEUE *jobQueue, struct ArrayEntryList *firstJobsinBlock,struct ArrayEntryList *lastBlockEntry) {
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
    else {                                                      //present only 1 element (to pop) in queue ==> update head and tail to NULL
        jobQueue->queueTail=NULL;
        jobQueue->queueHead=NULL;
    }
    struct ArrayEntry* job=jobEntry->arrayEntry;
#ifndef JOBS_ON_STACK
    free(jobEntry);
#endif
    return job;
}
void saveExponentVector(mpz_t *expVectorDst, FACTOR *factors,int numFactors,u_int64_t numPrimes) {
    //save into an mpz_t the exponent vector of factors with coeff. in GF(2)

    mpz_init2(*expVectorDst,numPrimes);
    for (int i = 0; i < numFactors; i++) {
        if (factors[i].exp%2==1)                              //mark factor in exp vector with 1 if exponent is odd
            mpz_setbit(*expVectorDst,factors[i].factor_indx);
    }
}

FACTORIZE_JOB_QUEUE * initFactorizeJobQueue(u_int64_t B, struct Precomputes *precomputes, int producersNum, int consumersNum) {
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
    if(pthread_cond_init(&(factorizeJobQueue->emptyAndClosedQueue),NULL)!=0) {
        fprintf(stderr, "cond var initialization error\n");
    }
    fprintf(stderr,"initialized condvar at %p mutex at %p\n",&(factorizeJobQueue->emptyAndClosedQueue),&(factorizeJobQueue->mutex));fflush(0);
//    if(!(factorizeJobQueue->queueHead=calloc(1,sizeof(struct ArrayEntryList)))){
//        fprintf(stderr,"job list head calloc errd\n");
//        free(factorizeJobQueue);
//        return NULL;
//    }
//    factorizeJobQueue->queueTail=factorizeJobQueue->queueHead;  //empty List have tail and head coincident
    ///termination vars
    factorizeJobQueue->producersEnded=0;
    factorizeJobQueue->producersNum=producersNum;
    factorizeJobQueue->consumersEnded=0;
    factorizeJobQueue->consumersNum=consumersNum;
    factorizeJobQueue->endedQueue = false;
    factorizeJobQueue->B=B;
    //// other usefull links
    factorizeJobQueue->precomputes=precomputes;
    return factorizeJobQueue;
}

void resetFactorizeJobQueue(FACTORIZE_JOB_QUEUE* factorizeJobQueue,int producersNum, int consumersNum)  {
    //reset default parameters for job queue
    factorizeJobQueue->producersEnded=0;
    factorizeJobQueue->producersNum=producersNum;
    factorizeJobQueue->consumersEnded=0;
    factorizeJobQueue->consumersNum=consumersNum;
    factorizeJobQueue->endedQueue = false;
}

void* FactorizeTrialDivide(void* args) {
    //trial divide factorize logic
    //lockstep factorization in a group of factorizer thread
    //each one will try a subset of primes up to Bound B
    //after a cumulation of founded factors or timeout (ON MACROs CONFIGURATION) the threads group will synchronize on barrier
    // the manager will aggregate all founded factors and reduce N (eliminating factors from it)
    // job done on either: 1! largePrime founded, compleate factorization over BSmooth numbers, max iteration reached (fail)
    // on each arrayEntry passed in the thread group pool exponent vector  (and eventually explicit factors) will be setted in place
    struct factorizerArgs *factorizerArgs = (struct factorizerArgs *) args;
    int factorizerID = factorizerArgs->factorizerID;
    int exp;    int barrierWaitReturn;
    ///timer vars for timed factorizations ierations
    STOPWATCH_INIT;  struct timeval FACTORIZE_BARRIER_SYNC_TIMEOUT =
            (struct timeval) {.tv_sec=FACTORIZE_BARRIER_SYNC_TIMEOUT_SEC,.tv_usec=FACTORIZE_BARRIER_SYNC_TIMEOUT_USEC};
    ///vars accessed only by thread manager
    int totalFoundedFactors = 0;                //total founded factors in all iterations
    int factorizeIteration = 0;                 //iteration elapsed since factorization start
    FACTOR *factors = NULL;  int factorsDynamicSize = FACTORS_IN_SIEVE_BLOCK_REALLOC;
    if (factorizerID == FACTORIZER_MANAGER_ID) {
        if (!(factors = malloc(sizeof(*factors) * factorsDynamicSize))) {
            fprintf(stderr, "factors malloc fail\n");
            return (void *) EXIT_FAILURE;
        }
    }
    ////primes for trial divide worker share will fetched from precomputed list distribuiting fairly primes share by primeIndx
    u_int64_t* primes=(u_int64_t*) factorizerArgs->precomputes->primes.pntr;
    u_int64_t primes_num=factorizerArgs->precomputes->primes.vectorSize;
    u_int64_t primeIndx, prime,primeTmp;                    //index of prime and related computation vars
    u_int64_t j;                                            //counter of primes tired during factorizing iteration
    int factorsFoundedIteration=0; //num of founded factors during this factorization iteration
    for (;;) {                  ///end less loop for worker thread --> they will break exitFlag setted by manager thread
//---------------------------------  ~~~~~~ threads sync 1 ~~~~~~    --------------------------------------------
        //// synchronize thread on N (factorize job) value <=== setted by thread manager ( last to join at the barrier)
        barrierWaitReturn = pthread_barrier_wait(factorizerArgs->barrier);
        if (barrierWaitReturn != PTHREAD_BARRIER_SERIAL_THREAD && barrierWaitReturn) {
            fprintf(stderr, "barrier error code:%d\n", barrierWaitReturn);
            return (void *) EXIT_FAILURE;
        }//NOW ALL THREAD CONCORDE WITH N VALUE (SETTED BY THREAD GROUP MANAGER) and eventual setted exit condition flag
#ifdef  DEBUG
        if(factorizerID==FACTORIZER_MANAGER_ID) {
            gmp_printf("factorizing %Zd, from thread:%lu iteration %d \n", *(factorizerArgs->N),pthread_self(),factorizeIteration );
            fflush(0);
        }
#endif
        /////////// factorize group end condition <-- setted by manager on queue end and closed
        if (*(factorizerArgs->exitCondition)) {   //exit cond setted by manager
            break;
        }
        /*
         * factorize tries may continue for a fiexed time or until an ammount of prime has been tried to divide N
         * these 2 different versions are selectable (not) defining macro TIMEBOXED_FACTORIZATION_ITERATION_POLLING
         */
        factorsFoundedIteration = 0;
#ifdef TIMEBOXED_FACTORIZATION_ITERATION_POLLING
        STARTSTOPWACHT(STOPWATCH_START_VAR)                    //start timer for timed factorization
#endif
        ////search for N factors in thread primes subset until K factors founded or timeout
        for (primeIndx = factorizerID,j=0;(primeIndx < primes_num && factorsFoundedIteration < FACTORS_ASYNC_PER_ITERATION)
#ifndef TIMEBOXED_FACTORIZATION_ITERATION_POLLING // --> fixed amount of tries--> exit condition with fixed set of max num of primes try
                                            && (++j)<MAX_PRIMES_TRY_PER_ITERATION;
#else
                                                                    ;
#endif
                                                                    primeIndx += FACTORIZER_THREAD_GROUP_SIZE) {
            prime = primes[primeIndx];
//            printf("new tring prime %lu from thread %lu\n",prime,pthread_self());fflush(0); //TODO EXTRA DEBUG

#ifdef TIMEBOXED_FACTORIZATION_ITERATION_POLLING       //--> exit condition of fixed max time per iteration
            if (++j % TIMEOUT_CHECK_INTERVAL_FACTORIZE_WORKER == 0) {      //timed factorization iterations
                STOPWACHT_STOP(STOPWATCH_START_VAR, STOPWATCH_STOP_VAR, STOPWATCH_DELTA_VAR)
                if (timercmp(&STOPWATCH_DELTA_VAR, &FACTORIZE_BARRIER_SYNC_TIMEOUT, >)) {
                    break;                                       //factorize job has timeout
                }
            }
#endif
            ///actual trial divide for fetched prime in worker primes subset
            exp = 0; primeTmp = prime;
            char buf[44];gmp_snprintf(buf,44,"%Zd",*(factorizerArgs->N));
            while ( mpz_divisible_ui_p(*(factorizerArgs->N),primeTmp)) {
                primeTmp *= prime;
                factorizerArgs->factorsTempStorage[factorizerID * FACTORS_ASYNC_PER_ITERATION + factorsFoundedIteration]
                                        = (FACTOR) {.factor=prime, .factor_indx=primeIndx, .exp=++exp};  //set founded factors to new location
            }
            if (exp)         //ANALOG TO SET FLAG IN DIVISION WHILE AND INCREMENT FACTOR INDEX HERE ! ONLY ONCE !
                factorsFoundedIteration++;

        } ///for end K factors founded or timeout expired


// ----------------------------------  ~~~~~~  threads sync 2 ~~~~~~    ------------------------------------------------
        //// synchronized factorizer thread to reduce N
        barrierWaitReturn = pthread_barrier_wait(factorizerArgs->barrier);
        if (barrierWaitReturn != PTHREAD_BARRIER_SERIAL_THREAD && barrierWaitReturn) {
            fprintf(stderr, "barrier error code:%d\n", barrierWaitReturn);
            return (void *) EXIT_FAILURE;
        } //NOW ALL THREAD CONCORDE ON FOUNDED FACTORS (in particular the manager)
        FACTOR factorFounded;
        /////// MANAGER AGGREAGTION
        if (factorizerID == FACTORIZER_MANAGER_ID) {
            for (int i = 0; i < FACTORIZER_THREAD_GROUP_SIZE * FACTORS_ASYNC_PER_ITERATION; i++) {
                factorFounded = factorizerArgs->factorsTempStorage[i];
                if (factorFounded.factor != 0) {    //aggregate factors locations only where thread has setted them
                    factors[totalFoundedFactors++] = factorFounded;
#ifdef DEBUG
                    printf("founded factor %lu ^ %d indx:%d iteration%d\n",factorFounded.factor,factorFounded.exp,i,factorizeIteration);fflush(0);
#endif
                    //TODO CHECK EFFICIENCY DIFFERENCES WITH EXP RYSE AND DIVIDE IN 2! OP.
                    for (int k = 0; k < factorFounded.exp; ++k) {
                        ///reduce N eliminating founded factors
                        mpz_tdiv_q_ui(*(factorizerArgs->N), *(factorizerArgs->N),factorFounded.factor);
                    }
                    factorizerArgs->factorsTempStorage[i].factor=0;     //reset factor entry to void for next iteration
                    REALLOC_WRAP(totalFoundedFactors, factorsDynamicSize, factors, FACTORS_IN_SIEVE_BLOCK_REALLOC)
                            return (void *) EXIT_FAILURE;
                        }}
                }
            }

            ///check if N factorization is done
            //either N fully factorized (now is 1) or N residue is a large prime (less then ~20*B and probably a prime)
            int N_has_been_fully_factorized=mpz_cmpabs_ui(*(factorizerArgs->N),1)==0;
            int N_is_large_prime= mpz_probab_prime_p(*(factorizerArgs->N), PROB_PRIME_REPS)  && mpz_cmp(*(factorizerArgs->N), *(factorizerArgs->largePrimeThreshold));
            if (N_has_been_fully_factorized || N_is_large_prime){
#ifdef VERBOSE

                printf(".\n");fflush(0);
#endif
#ifdef FINAL_RESIZE
                //realloc factors to actual needed size
                if (!(factors=realloc(factors, sizeof(*factors) * totalFoundedFactors))) {
                    fprintf(stderr, "final realloc of founded factors has failed\n");
                    return (void *) EXIT_FAILURE;
                }
#endif
                ///set in place into array entry the founded factors
#ifdef SAVE_FACTORS
                factorizerArgs->arrayEntryJob->factorsNum = totalFoundedFactors;
                factorizerArgs->arrayEntryJob->factors = factors;
#endif
                //save in place exp vector for matrix stage if this array entry yield to (partial) relation
                saveExponentVector(&(factorizerArgs->arrayEntryJob->exp_vector),factors,totalFoundedFactors,factorizerArgs->precomputes->primes.vectorSize);
                /// set large prime for partial relation
                if (N_is_large_prime){
                    mpz_init_set(factorizerArgs->arrayEntryJob->largePrime,*(factorizerArgs->N));
                }
#ifdef DEBUG
                gmp_printf("N quotient residue %Zd\tfounded %s \n",*(factorizerArgs->N),N_has_been_fully_factorized?"RELATION founded":"PARTIAL RELATION founded");fflush(0);
#endif
                break;
            }

            ////manager exit condition
            if (++factorizeIteration == FACTORIZATION_ITERATION_MAX) {
                free(factors);
#ifdef VERBOSE
              printf("MAX FACT ITERATIONS REACHED -->FACTORIZATION FAILED\n");
#endif
                return (void *) EXIT_FAILURE;
            }
        }
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
        free(barrier);
        return (void*) EXIT_FAILURE;
    }
    pthread_t* factorizeWorkers=calloc(FACTORIZER_THREAD_GROUP_SIZE, sizeof(*factorizeWorkers));
    FACTOR* factorizersTempStorage=calloc(FACTORIZER_THREAD_GROUP_SIZE*FACTORS_ASYNC_PER_ITERATION, sizeof(*factorizersTempStorage));
    if(!factorizeWorkers || !factorizersTempStorage){
        fprintf(stderr,"workers malloc errd\n");
        free(factorizeWorkers);free(factorizersTempStorage);
        return (void*)EXIT_FAILURE;
    }
    mpz_t N_cp;  mpz_init(N_cp);        //job N copy for in place division(only by manager)
    //set large prime threashold
    mpz_t LARGE_PRIME_THREASHOLD;
    mpz_init_set_ui(LARGE_PRIME_THREASHOLD,factorizeJobQueue->B);
    mpz_mul_ui(LARGE_PRIME_THREASHOLD,LARGE_PRIME_THREASHOLD,LARGE_PRIME_THREASHOLD_FACTOR);

    //factorizer arg template --> only few change needed in copyes 5 for workers
    struct factorizerArgs factorizerArgTemplate=(struct factorizerArgs){
            .N=&N_cp,
            .B=factorizeJobQueue->B,
            .precomputes=factorizeJobQueue->precomputes,
            .barrier=barrier,
            .factorsTempStorage=factorizersTempStorage,
            .largePrimeThreshold=&LARGE_PRIME_THREASHOLD,
            .exitCondition=&exitCond
    };

    //// start factorizers worker as detached thread  --> they will be stopped manually on jobs end
//    pthread_attr_t attr; pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
    struct factorizerArgs factorizerArgs[FACTORIZER_THREAD_GROUP_SIZE];
    for (int i = 0; i < FACTORIZER_THREAD_GROUP_SIZE-1; i++) {  //-1 because the manager is the caller so already created
//        memcpy(factorizerArgs,&factorizerArgTemplate, sizeof(factorizerArgTemplate));   //set template arg on factorizer i arg
        factorizerArgs[i]=factorizerArgTemplate;
        factorizerArgs[i].factorizerID=i;
        if (pthread_create(factorizeWorkers + i, NULL, FactorizeTrialDivide, (void *) (factorizerArgs+i)) != 0){
            fprintf(stderr, " factorizer worker creataing error at %d\n",i);
            exit(EXIT_FAILURE);
        }
//        printf("created factorizer worker %d\n",i);
    }
    //////// job fetch and workers comunication MUTEX FETCHed
    struct ArrayEntry *job; int errCode=0,closedQueue;//non zero if queue has been closed by producers
    bool endedQueue=false;
    for(;!exitCond;) {
        LOCK_MUTEX(&(factorizeJobQueue->mutex),errCode)
        job = popFirstJob(factorizeJobQueue); //// safely pop a job from job queue
        closedQueue = factorizeJobQueue->producersEnded == factorizeJobQueue->producersNum;
        if(!job && closedQueue){    //ended queue
            if(!factorizeJobQueue->endedQueue) {
                //--> notify all producer in wait for results only if not already did
                factorizeJobQueue->endedQueue=true;     //set that the queue has been ended=> no more broadcast needed
                endedQueue=true;                        //only this thread will have this flag set-->used to quit other threads
            }
            factorizeJobQueue->consumersEnded++;        //mark ended consumer
            exitCond=true;
        }
        UNLOCK_MUTEX(&(factorizeJobQueue->mutex),errCode)
        if(!job && !closedQueue){    // here it means too fast job consuming
           /// if queue is not closed it means that thread managers consume jobs to fastly or bad luck in thread scheduling
           /// so wait a bit and retry
           usleep(QUEUE_FILLING_POLLING_USLEEP);
           continue;
        }
        if(job) {
            mpz_set(*(factorizerArgTemplate.N), job->element);   //set a copy of N for local computation on factorizer
            factorizerArgTemplate.arrayEntryJob = job;            //set array entry link for in place output write (only by manager)
            if (!mpz_cmp_ui(job->element,1)|| !mpz_cmp_ui(job->element,0))                //TODO ROBUSTNESS----???----
                continue;
#ifdef DEBUG_MANAGER
            gmp_printf("picked job:%Zd\n",job->element);
#endif
        }    // otherwhise comunicate setted exit condition to workers with exitCond Setted
        /////// sync comunication of newly fetched job to worker --> by mem barrier of barrier wait -> wakeup because manager entered
        //2 settings above will be visible to worker on barrier wait wake up because of memory barrier
        factorizerArgTemplate.factorizerID=FACTORIZER_MANAGER_ID;
        int factorizeJobRes=(int) FactorizeTrialDivide((void*)&factorizerArgTemplate);               //manager will participate to concurrent factorization with his special ID
/*#ifdef DEBUG_MANAGER
        if(factorizeJobRes==EXIT_SUCCESS && job){
            fflush(0);gmp_printf("%Zd founded factors\t", (job->element));
            for (u_int i = 0; i < job->factorsNum; ++i) {
                printf("%lu^%d\t",job->factors[i].factor,job->factors[i].exp);
            }
            printf("\n\n-&:$-\n\n");fflush(0);
        }
#endif*/
    }
    //join factorizer thread for safe barrier destroy
    for (int j = 0; j < FACTORIZER_THREAD_GROUP_SIZE; ++j) {
        if(pthread_join(factorizeWorkers[j],NULL)<0){
            fprintf(stderr,"INVALID JOIN OF WORKER THREAD");
        }
    }
    //// only marked thread manager will safelly unlock sievers with condvar when each factorize group has ended
    if (endedQueue) {
        //// spinlock on consuemrs ended counter (hopefully short wait)
        while(factorizeJobQueue->consumersEnded<factorizeJobQueue->consumersNum){
            usleep(FACTORIZE_END_POLLING_USEC);
        }
        // now all job has been taked from queue and ended so siever will se updated array
        fprintf(stderr,"END OF QUEUE on condvar %p and mutex %p \n",&(factorizeJobQueue->emptyAndClosedQueue),&(factorizeJobQueue->mutex));fflush(0);
        if (pthread_cond_broadcast(&(factorizeJobQueue->emptyAndClosedQueue))!=0) {
            fprintf(stderr, "broadcast on condvar err\n");
            return (void *) EXIT_FAILURE;
        }
    }
    ///mem clean
    pthread_barrier_destroy(barrier);free(barrier);
    free(factorizeWorkers); free(factorizersTempStorage);
#ifdef VERBOSE
    printf("GROUP_MANAGER_EXIT\t:%lu\n",pthread_self());fflush(0);
#endif
    return (void*)EXIT_SUCCESS;
}

pthread_t* StartFactorizerThreadGroups(FACTORIZE_JOB_QUEUE* factorizeJobQueue,int numGroupsFactorizers){
    // start numGroupFactorizers group of factorizeManger --> workers simply starting only group managers
    pthread_t* pthreads=malloc(sizeof(*pthreads)*numGroupsFactorizers);
    if(!pthreads){
        fprintf(stderr,"out of mem for thread managers\n");
        return NULL;
    }
//    start thread as detached
//    pthread_attr_t attr;     pthread_attr_init(&attr);     pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
    for (int i = 0; i <numGroupsFactorizers ; ++i) {
        if(pthread_create(pthreads+i,NULL,FactorizeTrialDivideThreadGroupManager,factorizeJobQueue)<0){
            fprintf(stderr,"creation of factorizer manager thread errd\n");
            free(pthreads);
            return NULL;
        }
    }
    return pthreads;
}
int JoinFactorizerThreadGroups(pthread_t* threadManager,int numManager){
    //wait for manager to exit --> that will happend when all factorize job has compleated and 1 manager will notify siever of that event (condvar)
    //simply join thread and free mem
    for (int i = 0; i < numManager; ++i) {
        int threadRetErrCode;
        if ((threadRetErrCode = pthread_join(threadManager[i], NULL   ))) {
            fprintf(stderr, "factorizer jERR\n %s__\n", strerror(threadRetErrCode));
            return threadRetErrCode;
        }
    }
    free(threadManager);
    return EXIT_SUCCESS;
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
    const char* n_str="100000030925519250968982645360649";
    CONFIGURATION *configuration = initConfiguration(n_str, 0, 0, 0, 0);
    //// get first polynomial:
    struct polynomial pol;
    PRECOMPUTES *precomputes = preComputations(configuration, &pol,NULL);
    u_int64_t B=400000;
    FACTORIZE_JOB_QUEUE *factorizeJobQueue = initFactorizeJobQueue(B, precomputes, configuration->SIEVING_THREAD_NUM,
                                                                   0);
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
    const int JOBS_NUM_BLOCK=100;
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
    appendBlockJobs(factorizeJobQueue, jobs, jobs+JOBS_NUM_BLOCK-1);
    appendJob(factorizeJobQueue,&entryJob2);

    factorizeJobQueue->producersEnded++;
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
