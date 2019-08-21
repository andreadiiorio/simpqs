#ifndef SIMPQS_SIEVINGSIMPQS_H
#define SIMPQS_SIEVINGSIMPQS_H

#include "SIMPQS.h"
#include <mpfr.h>
/*TODO ALTERNATIVES FOR SIEVE ARRAY
 *  1) ONLY f(x) VALUES AND SHADOW ARRAY        --> ??? QUICKER ARRAY MOVE ??? && SMALLEST FOOTPRINT VS INDEX SYNCRONOSNESS
 *  3) LARGE PRIME RECORD SEPARATED             -->LIGHTER MEM FOOT PRINT IN SIEVE VS MEM COPY OVERHEAD
 *  2) ENTRY ARRAY                              --> "
 */
typedef mpz_t* SIEVE_ARRAY_BLOCK ;
typedef struct Shadow_Array{
    u_int64_t j;
    ///log sieving
    ////////TODO FIND FASTER ALTERNATIVE
    mpfr_t logSieveCumulative;      //most space consumer, NO MPFR CONVERSION NEEDED --> quicker sieve
//    double logSieveCumulative;      //space saver, faster comparison,
//    u_int64_t logSieveCumulative;   //rounding  double version
    u_int64_t* factors;
    ///large prime  TODO vedi  *ALTERNATIVA
    mpz_t* largePrime;              //not NULL if element has 1! Large prime in his factorization
//    char* largePrime;
} SIEVE_ARRAY_BLOCK_SHADOW;
typedef mpz_t* SIEVE_ARRAY_BLOCK ;
struct ArrayEntry{
    u_int64_t j;
    mpz_t element;                  //f(j)
    ///log sieving
    ////////TODO FIND FASTER ALTERNATIVE
    mpfr_t logSieveCumulative;      //most space consumer, NO MPFR CONVERSION NEEDED --> quicker sieve
//    double logSieveCumulative;      //space saver, faster comparison,
//    u_int64_t logSieveCumulative;   //rounding  double version
    ///array entry factorization for BSmooth values
    u_int64_t* factors;              //non NULL element has logSieveCumulative above threashold
    ///large prime  TODO vedi  *ALTERNATIVA
    mpz_t* largePrime;              //not NULL if element has 1! Large prime in his factorization
//    char* largePrime;
};

//struct largePrimeBucket{ //TODO *ALTERNATIVA --> quicker largePrime init/set vs lighter memory in sieve
//    ///BUCKET KEY
////    mpz_t largePrime;
//    char* largePrime;
//    ///BUCKET VALUE
//    u_int64_t* factors;         //factors for matrix row
//    u_int64_t j;                //element identifier
//};



#endif //SIMPQS_SIEVINGSIMPQS_H
