#ifndef SIMPQS_SIEVINGSIMPQS_H
#define SIMPQS_SIEVINGSIMPQS_H

#include "SIMPQS.h"
#include <mpfr.h>
struct ArrayEntry{
    int64_t j;
    mpz_t element;                  //f(j)
    ///log sieving
    ////////TODO FIND FASTER ALTERNATIVE
    double logSieveCumulative;      //space saver, faster comparison,
//    u_int64_t logSieveCumulative;   //rounding  double version

    u_int64_t* factors;              //non NULL element has logSieveCumulative above threashold
    ///large prime
    mpz_t* largePrime;              //not NULL if element has 1! Large prime in his factorization
    //TODO DISCRIMINATE ARRAYENTRIES:
    //  factors && largePrime == NULL ==> useless entry ==> skippable
    //  largePrime !=NULL (=> factors not NULL) ==> largePrime to accopiate
    //  largePrime==NULL && factors !=NULL  ==> BSmooth entry => direct matrix row
};
typedef struct ArrayEntry* SIEVE_ARRAY_BLOCK;


int Sieve(struct Configuration* config, struct Precomputes* precomputes,SIEVE_ARRAY_BLOCK sieveArrayBlock);
void* siever_thread_logic(void* arg);       //siever pthread func
#endif //SIMPQS_SIEVINGSIMPQS_H
