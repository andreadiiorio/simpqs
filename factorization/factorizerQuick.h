
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
    struct Precomputes* precomputes;
    /// queue termianation vars
    int producersNum;                     //number of producer thread on the queue
    int producersEnded;                   //number of producer thread  that has finished the queue
    int consumersNum;
    int consumersEnded;                   //number of factorize thread manager that has ended
    pthread_cond_t emptyAndClosedQueue;   //closed and flushed job queue --> setted from consumers at queue end
    bool endedQueue;                      //true if condvar broadcast already called
} FACTORIZE_JOB_QUEUE;
void appendJob(FACTORIZE_JOB_QUEUE *jobQueue,struct ArrayEntryList* newJob);
struct ArrayEntry* popFirstJob(FACTORIZE_JOB_QUEUE *jobQueue);
void appendBlockJobs(FACTORIZE_JOB_QUEUE *jobQueue, struct ArrayEntryList *firstJobsinBlock,struct ArrayEntryList *lastBlockEntry);
FACTORIZE_JOB_QUEUE *
initFactorizeJobQueue(u_int64_t B, struct Precomputes *precomputes, int producersNum, int consumersNum);
void resetFactorizeJobQueue(FACTORIZE_JOB_QUEUE* factorizeJobQueue,int producersNum, int consumersNum);


struct factorizerArgs{
    mpz_t* N;                        //num to factorize pntr
    u_int64_t B;
    struct Precomputes* precomputes;
    int factorizerID;
    ////fized size factors mem write safe for workers thread (each one has his set of factors location separated)
    ///for each worker thread in the group there're  FACTORS_ASYNC_PER_ITERATION FACTORS  write sa
    ///FACTOR* factorizersTempStorage=malloc(sizeof(FACTOR)*FACTORIZER_THREAD_GROUP_SIZE*FACTORS_ASYNC_PER_ITERATION);
    FACTOR* factorsTempStorage;
    mpz_t* largePrimeThreshold;
    pthread_barrier_t* barrier;     //barrier for lockstep factorization
    //job queue
    struct ArrayEntry* arrayEntryJob;  //job to do
    bool* exitCondition;
};

void* FactorizeTrialDivide(void* args);
pthread_t* StartFactorizerThreadGroups(FACTORIZE_JOB_QUEUE* factorizeJobQueues,int numGroupsFactorizers);
int JoinFactorizerThreadGroups(pthread_t* threadManager,int numManager);

void arrayEntryCopy(struct ArrayEntry *destEntry, struct ArrayEntry *entry);
int mergeReports(REPORTS *dstReports, const REPORTS *new_reports);
int mergeReportsFast(REPORTS *dstReports, const REPORTS *new_reports); //TODO DEBUG MEMMOVE CORRUPT EVERYTHING
#endif //SIMPQS_FACTORIZERQUICK_H
