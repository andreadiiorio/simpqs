#ifndef SIMPQS_CONFIGURATION_H
#define SIMPQS_CONFIGURATION_H

//// configuration of SIMPQS
//dflt config of main algo parameters
#define _ARRAY_IN_MEMORY_MAX_SIZE 5000      //max num of array elements holded in memory
#define _M 7096
#define _B 10000
#define _SIEVING_THREAD_NUM  3              //siever threads
#define FINAL_RESIZE 1                       //EVERY DINAMIC ARRAY WILL BE RESIZED TO ACTUAL SIZE NEEDED

//// dynamically resized array configurations
#define PROBABLY_BSMOOTH_ARRAY_BLOCK 1024
#define PARTIAL_RELATIONS_BLOCK 1024
#define FACTORIZE_JOB_BLOCK_APPEND 32       //block of factorize job to atomically append on mutex protected job queue

//// factorization at end of sieving stage configs
#define FACTORIZER_THREAD_GROUP_SIZE    5   //num of thread used to factorize array elements probably BSmooth (included the manager)
#define NUM_FACTORIZER_GROUPS 11            /// number of factorizer groups
    //// concurrent trial divide factorization MODES
//iteration mode timeboxed trigger (if non defined fixed set of try mode enabled)
#define TIMEBOXED_FACTORIZATION_ITERATION_POLLING   ///ITERATION MODE TRI
#define MAX_PRIMES_TRY_PER_ITERATION 100000

#define TIMEOUT_CHECK_INTERVAL_FACTORIZE_WORKER 5969                 //number of trial division to do before check if sub job has timeout
#define FACTORS_ASYNC_PER_ITERATION 6                              //max num of factors to find per iteration
#define FACTORIZER_MANAGER_ID (FACTORIZER_THREAD_GROUP_SIZE -1)
#define FACTORS_IN_SIEVE_BLOCK_REALLOC 96
#define PROB_PRIME_REPS 44                                          //number of repetition of miller rabin pseudo primality test
#define FACTORIZATION_ITERATION_MAX 3                               //max number of factorization iteration in a factorization thread group
#define LARGE_PRIME_THREASHOLD_FACTOR 20                            //will be accepted primes up to B*THIS_MACRO
///timeouts
//factorization iteration timeout
#define FACTORIZE_BARRIER_SYNC_TIMEOUT_SEC 0
#define FACTORIZE_BARRIER_SYNC_TIMEOUT_USEC 6000
#define QUEUE_FILLING_POLLING_USLEEP 1000                            //u sleep to wait job is filled by a producer of job in the queue

#define PRIMES_32B_PATH "../primes.32b"
/// debug flag
//#define DEBUG                                                       //debug print of
//#define DEBUG_MANAGER                                               //factorization thread manager debug print
#endif //SIMPQS_CONFIGURATION_H
#define VERBOSE_0
