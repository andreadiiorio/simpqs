#ifndef SIMPQS_SIEVINGSIMPQS_H
#define SIMPQS_SIEVINGSIMPQS_H

#include "SIMPQS.h"
#include <mpfr.h>
#include <string.h>

typedef struct factor{
    u_int64_t factor;
    int         exp;
}FACTOR;
struct ArrayEntry{
    int64_t j;
    mpz_t element;                  //f(j)
    ///log sieving
    ////////TODO FIND FASTER ALTERNATIVE
    long double logSieveCumulative;      //space saver, faster comparison,
//    u_int64_t logSieveCumulative;   //rounding  double version

    FACTOR* factors;              //non NULL element has logSieveCumulative above threashold
    u_int factorsNum;
    ///large prime
    mpz_t largePrime;              //not NULL if element has 1! Large prime in his factorization
    //TODO DISCRIMINATE ARRAYENTRIES:
    //  factors && largePrime == NULL ==> useless entry ==> skippable
    //  largePrime !=NULL (=> factors not NULL) ==> largePrime to accopiate
    //  largePrime==NULL && factors !=NULL  ==> BSmooth entry => direct matrix row
};
typedef struct ArrayEntry* SIEVE_ARRAY_BLOCK;


/////////global
extern mpz_t   N;
extern SIEVE_ARRAY_BLOCK SieveArrayBlock;              ///array block in memory
extern struct Configuration Config;
int Sieve(struct Configuration* config, struct Precomputes* precomputes,SIEVE_ARRAY_BLOCK sieveArrayBlock);
void* siever_thread_logic(void* arg);       //siever pthread func
#endif //SIMPQS_SIEVINGSIMPQS_H
