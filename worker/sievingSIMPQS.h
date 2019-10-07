#ifndef SIMPQS_SIEVINGSIMPQS_H
#define SIMPQS_SIEVINGSIMPQS_H

#include "SIMPQS.h"
#include <mpfr.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#define DEBUG_CHECK
#define SAVE_FACTORS  //save factorization of entry other then exponent vector

typedef struct factor{
    u_int64_t factor;
    u_int64_t factor_indx;      //index of factor in primes list
    int         exp;            //factor exponent
}FACTOR;
struct ArrayEntry{
    int64_t j;
    mpz_t element;                  //mongomery pol value a*f(j)=(a*j^2+2*b*j+c)==(a*j+b)^2-n
    mpz_t x;                        //(a*j+b) --> it hold Mod(x^2,n)=Mod(a*f(j),n)
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
    mpz_t              n;                       //convenience ref to target N factor
} REPORTS;

REPORTS* Sieve(struct Configuration *config, struct Precomputes *precomputes, struct polynomial* actualPol);
void* siever_thread_logic(void* arg);       //siever pthread func
///// REPORTS
int saveReports(REPORTS *reports, u_int64_t colsN, bool printReports, bool polynomialFamilyReports, struct polynomial* polynomial);
REPORTS *loadReports(char *filePath);
void print_reports(REPORTS *reports, u_int64_t colsN, bool printMatrix);
void sortPartialReports(REPORTS* reports);
int pairPartialReports(REPORTS* reports);
void arrayEntryCopy(struct ArrayEntry *destEntry, struct ArrayEntry *entry);
int mergeReports(REPORTS *dstReports, const REPORTS *new_reports);
int checkReports(REPORTS *reports, bool aggregatedLargePrimeCheck);
void CHECK_X_SQURARE_CONGRUENT_Y_MOD_N(struct ArrayEntry* arrayEntry, mpz_t tmp, mpz_t tmp2, bool largePrimeAggregatedEntryCheck);
REPORTS* aggregateSieversWorkers(const unsigned int polynomialN);

#ifndef REPORTS_PATHS
#define REPORTS_PATHS
#define FIND_REPORTS_POL_FAMILY_BASH_CMD_LINUX "find \"$(pwd -P)\" -iname \"reports_*\""
#define REPORTS_FILENAME_PREFIX  "reports_"
#define REPORTS_FILENAME_SUFFIX  ".reportslist"
#define REPORTS_POLYNOMIAL_FAMILY_FILENAME_SUFFIX  ".reportsfamilylist"
#endif
char** findReportsLocally(unsigned int numReports,const char* reportSuffix);
void _deleteLocalReports();

#endif //SIMPQS_SIEVINGSIMPQS_H
