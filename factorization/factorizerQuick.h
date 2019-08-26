
#ifndef SIMPQS_FACTORIZERQUICK_H
#define SIMPQS_FACTORIZERQUICK_H

#include "SIMPQS.h"
#include <stdlib.h>
#include "utils/utils.h"
#include <CONFIGURATION.h>
#include <stdio.h>
#include <pthread.h>
#include <gmp.h>
#include <worker/sievingSIMPQS.h>
#include <stdbool.h>

////factorize job queue for job dispatching to threadGroup Pool
struct ArrayEntryList{  //double linkedList
    struct ArrayEntry* arrayEntry;
    struct ArrayEntryList* nextArrayEntry;

};
typedef struct FactorizeJobQueue{
    struct ArrayEntryList* queueHead;    //list head    --> if null empty queue
    struct ArrayEntryList* queueTail;    //list head
    pthread_mutex_t mutex;                  //syncronize access to job queue
    u_int64_t B;
    u_int64_t* primesLink;
    u_int64_t  primesN;
    bool closedQueue;                     //true if the queue has been closed
} FACTORIZE_JOB_QUEUE;
void appendJob(FACTORIZE_JOB_QUEUE *jobQueue,struct ArrayEntryList* newJob);
struct ArrayEntry* popFirstJob(FACTORIZE_JOB_QUEUE *jobQueue);
void appendBlockJobs(FACTORIZE_JOB_QUEUE *jobQueue, struct ArrayEntryList* firstJobsinBlock,struct ArrayEntryList* lastJobInBlockPntr) ;
FACTORIZE_JOB_QUEUE* initFactorizeJobQueue(u_int64_t B,DYNAMIC_VECTOR primes);

struct factorizerArgs{
    mpz_t* N;                        //num to factorize pntr
    u_int64_t B;
    u_int64_t  primesNum;
    u_int64_t* primes;              // primes<= B
    int factorizerID;
    ////fized size factors mem write safe for workers thread (each one has his set of factors location separated)
    ///for each worker thread in the group there're  FACTORS_ASYNC_PER_ITERATION FACTORS  write sa
    ///FACTOR* factorizersTempStorage=malloc(sizeof(FACTOR)*FACTORIZER_THREAD_GROUP_SIZE*FACTORS_ASYNC_PER_ITERATION);
    FACTOR* factorsTempStorage;
    pthread_barrier_t* barrier;     //barrier for lockstep factorization
    //job queue
    struct ArrayEntry* arrayEntryJob;  //job to do
    bool* exitCondition;
};

void* FactorizeTrialDivide(void* args);
#endif //SIMPQS_FACTORIZERQUICK_H
