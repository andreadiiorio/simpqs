#include <stdlib.h>
#include <stdio.h>
#include <utils/utils.h>
#include "sievingSIMPQS.h"

mpz_t tmp;
#define OVERFLOW_CHECK_POLYNOMIAL   //perform overflow check
void polynomialValueQuick(int64_t j,mpz_t outputVal,struct polynomial_actual polynomial){
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
int sieveArrayBlockAllocate(struct Configuration* config,SIEVE_ARRAY_BLOCK* sieveArrayBlockPntr){
    //allocate array block of blockSized returning error to caller
    *sieveArrayBlockPntr=calloc(config->ARRAY_IN_MEMORY_MAX_SIZE, sizeof(mpz_t));
    if(!(*sieveArrayBlockPntr)){
        fprintf(stderr,"array block calloc error \n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

void* sieve_thread_logic(void* arg){
    /*
     * sieve a fair array portion for Bsmooth values
     *  ->add log(p) with p in FB to sub array elements divisible for p
     *  ->directly factorize array elements that are above LOG_THREASHOLD TODO**
     *  ->join to caller main thread
     *  TODO**
     *      1) subThread Factorizers create and exec
     *      2) locked pipe to Factorizer Process ( threading inside it)
     *      3) return array locations above LOG_THREASHOLD
     */

}
int sieve(struct Configuration* config, struct Precomputes* precomputes,SIEVE_ARRAY_BLOCK sieveArrayBlock){
    /*sieve concurrently an array of Mongomery polynomial values divided it in fixed sized blocks sieveArrayBlock
      todo sieveArrayBlock expected to be already allocated
      precomputations will be used to find quickly array location divisible per primes in factorbase
      log(p) will be added to shadow location of array A[j] for each p that divide A[j]
      array entries that will be more than a specific threashold will be factorized (hopefully quickly) TODO TIMEOUTED SKIP FACT TRY

      Large primes identified will produce "partial" relation TODO C-HASHMAP VS SUPER ENTRY LIST, SORTED FOR LARGE PRIME MATCHING
      Bsmooth values will produce relation
      relation will produce matrix rows for matrix step
    */



    return EXIT_SUCCESS;
}
