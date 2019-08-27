//
// Created by andysnake on 19/08/19.
//

#ifndef QUADRATIC_SIEVE_UTILS_H
#define QUADRATIC_SIEVE_UTILS_H


#include "SIMPQS.h"
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/time.h>
#include <gmp.h>
#include "CONFIGURATION.h"
//initialize multiple precision int from sourceString in dest
#define strToMpz(dest,sourceStr)\
    if (mpz_init_set_str(dest,sourceStr, 10) == -1) {\
    printf("Cannot load N %s\n", argv[1]);\
    exit(2);\
}

#define LONG_MUL_OVERFLOW_P(a, b) \
   __builtin_mul_overflow_p (a, b, (__typeof__ ((a) * (b))) 0)

void printPrecomputations(struct Precomputes* precomputes,int blockPrint);
typedef struct DynamicVector{
    void* pntr;                         //to cast pntr
    u_int64_t vectorSize;               //actual vector size
}DYNAMIC_VECTOR ;
//read primes from precomputed primes list file until reached smoothnessBound
DYNAMIC_VECTOR ReadPrimes(char *primesListPath, u_int64_t smoothnessBound );
//read primes until bound reached and filter them for them with legendre (N/p)==1
DYNAMIC_VECTOR ReadFactorBase(char *primesListPath, u_int64_t smoothnessBound, mpz_t N);

//// wrap resizing evalutation and realloc for a dynamic list
#define REALLOC_WRAP(cumulativeCount,dynamicVectSize,dynamicVectPntr,DYNAMIC_VECT_BLOCK_SIZE)\
    if ( cumulativeCount> dynamicVectSize) {   \
        dynamicVectSize+= DYNAMIC_VECT_BLOCK_SIZE; \
        if (!(dynamicVectPntr=realloc(dynamicVectPntr, dynamicVectSize * sizeof(*dynamicVectPntr)))) {  \
            fprintf(stderr, "realloc error\n"); \
            free(dynamicVectPntr);
            //TODO CONTINUE MACRO ON CALL
            //  needed }}
/*
 * stopwatch MACROs to get elapsed time of an operation with simply call START AND STOP
 * meaning by name of used timevals vars in STOPWATCH_INIT
 * result will be printed on stdout
 */
///TMP VAR NAMES MACROS
#define STOPWATCH_START_VAR startime
#define STOPWATCH_STOP_VAR stoptime
#define STOPWATCH_DELTA_VAR delta
///STOP WATCH BASIC OP MACROS
#define STOPWATCH_INIT struct timeval STOPWATCH_START_VAR,STOPWATCH_STOP_VAR,STOPWATCH_DELTA_VAR
#define STARTSTOPWACHT(startime)\
    gettimeofday(&(startime),NULL);

#define STOPWACHT_STOP(startime,stoptime,delta)\
    gettimeofday(&(stoptime),NULL);\
    timersub(&(stoptime),&(startime),&(delta));

#define STOPWACHT_STOP_DFLT\
    gettimeofday(&(STOPWATCH_STOP_VAR ),NULL);\
    timersub(&(STOPWATCH_STOP_VAR),&(STOPWATCH_START_VAR),&(STOPWATCH_DELTA_VAR));  \
    fprintf(stderr,"\n\n Elapsed secs and micros: %ld , %ld \n\n",(delta).tv_sec,(delta).tv_usec);fflush(0);


#endif //QUADRATIC_SIEVE_UTILS_H
