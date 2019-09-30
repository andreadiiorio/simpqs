#ifndef SIMPQS_SIEVINGSIMPQS_H
#define SIMPQS_SIEVINGSIMPQS_H

#include "SIMPQS.h"
#include <mpfr.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define SAVE_FACTORS  //save factorization of entry other then exponent vector

typedef struct factor{
    u_int64_t factor;
    u_int64_t factor_indx;      //index of factor in primes list
    int         exp;            //factor exponent
}FACTOR;
struct ArrayEntry{
    int64_t j;
    mpz_t element;                  //f(j)
    ///log sieving
    ////////TODO FIND FASTER ALTERNATIVE
    long double logSieveCumulative;      //space saver, faster comparison,
//    u_int64_t logSieveCumulative;   //rounding  double version
#ifdef SAVE_FACTORS
    FACTOR* factors;              //non NULL element has logSieveCumulative above threashold
#endif
    /*
     * will hold exponent vector of element if it produce (partial) relation
     * "bit" of this variable has same indexing of loaded primes of FactorBase
     *  1 at bit i <=> element has odd power of primes[i]
     *  used in matrix stage
     */
    mpz_t exp_vector;
    u_int factorsNum;
    ///large prime
    mpz_t largePrime;              //not NULL if element has 1! Large prime in his factorization
    //TODO DISCRIMINATE ARRAYENTRIES:
    //  factors && largePrime == NULL ==> useless entry ==> skippable
    //  largePrime !=NULL (=> factors not NULL) ==> largePrime to accopiate
    //  largePrime==NULL && factors !=NULL  ==> BSmooth entry => direct matrix row
};
typedef struct ArrayEntry* SIEVE_ARRAY_BLOCK;

typedef struct report{
    ///report infos from sieving
    struct ArrayEntry* bsmoothEntries;
    u_int64_t          relationsNum;             //number of BSmoothEntries founded
    ///partial reports
    struct ArrayEntry* largePrimesEntries;       //array entries holding a BsmoothVal times a largePrime --> exploitable to generate more reports
    u_int64_t          partialRelationsNum;      //number of LargePrimes entries founded
} REPORTS;

void print_reports(REPORTS *reports, u_int64_t colsN, bool printMatrix);

/////////global
extern mpz_t   N;
extern SIEVE_ARRAY_BLOCK SieveArrayBlock;              ///array block in memory
extern struct Configuration Config;
REPORTS* Sieve(struct Configuration *config, struct Precomputes *precomputes, SIEVE_ARRAY_BLOCK sieveArrayBlock,struct polynomial* actualPol);
void* siever_thread_logic(void* arg);       //siever pthread func
#endif //SIMPQS_SIEVINGSIMPQS_H
