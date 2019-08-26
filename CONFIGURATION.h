#ifndef SIMPQS_CONFIGURATION_H
#define SIMPQS_CONFIGURATION_H

//// configuration of SIMPQS
#define _ARRAY_IN_MEMORY_MAX_SIZE 5000      //max num of array elements holded in memory
#define _SIEVING_THREAD_NUM 7               //siever threads
#define FINAL_RESIZE 1                       //EVERY DINAMIC ARRAY WILL BE RESIZED TO ACTUAL SIZE NEEDED
//// factorization at end of sieving stage configs
#define FACTORIZER_THREAD_GROUP_SIZE    30   //num of thread used to factorize array elements probably BSmooth (included the manager)
#define FACTORIZE_THREAD_GROUPS_NUM     5   //num of factorize thread group (equal to number of factorizer manager)
#define TIMEOUT_CHECK_INTERVAL_FACTORIZE_WORKER 596                 //number of trial division to do before check if sub job has timeout
#define FACTORS_ASYNC_PER_ITERATION 60                              //max num of factors to find per iteration
#define FACTORIZER_MANAGER_ID (FACTORIZER_THREAD_GROUP_SIZE -1)
#define FACTORS_IN_SIEVE_BLOCK_REALLOC 96
#define PROB_PRIME_REPS 44                                          //number of repetition of miller rabin pseudo primality test
#define FACTORIZATION_ITERATION_MAX 16                               //max number of factorization iteration in a factorization thread group
///timeouts
//factorization iteration timeout
#define FACTORIZE_BARRIER_SYNC_TIMEOUT_SEC 0
#define FACTORIZE_BARRIER_SYNC_TIMEOUT_USEC 1000

#define PRIMES_32B_PATH "/home/andysnake/Desktop/tenPrj/primes.32b"
#endif //SIMPQS_CONFIGURATION_H
