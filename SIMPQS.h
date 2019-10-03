#ifndef SIMPQS_SIMPQS_H
#define SIMPQS_SIMPQS_H

#include <stdlib.h>
#include <sys/types.h>
#include <gmp.h>
typedef struct DynamicVector{
    void* pntr;                         //to cast pntr
    u_int64_t vectorSize;               //actual vector size
}DYNAMIC_VECTOR ;

///precomputated vars
struct precomputatios_polynomial{   //precomputations polynomial dependent -> will change (slightly) at each polynomial change
    //same indexing of FactorBase
    int64_t* sol1p;                //   sieving constant jumps
    int64_t* sol2p ;               //   sieving constant jumps
    mpz_t* B_ainv_2Bj_p;           //   2*Bj*(a^-1) mod p  indexed by (cols) j in 1..s (rows) in p1,..,pn (2B_j1_p1,2Bj2_p1....) in 1 D
    mpz_t* B_l;                    //   generators of b by GreyCode; indexing 1<l<s where a=p1*..*ps
    int s;                          //number of B_l --> num of factors of a
};
/// compute f(x)=a*(a*j^2+2*b*j+c) //TODO j^2 OVERFOW POSSIBLE!!!
#define POLYNOMIAL_VAL_COMPUTE_A_F(j, outputVal,tmp, polynomial)\
    mpz_mul_si(outputVal, (polynomial->a), (j*j) );    \
    mpz_mul_si(tmp, (polynomial->b), (2*j) );    \
    mpz_add(outputVal,outputVal,tmp);\
    mpz_add(outputVal,outputVal,(polynomial->c));\
    mpz_mul(outputVal,outputVal,(polynomial->a));



/// compute (aj+b)
#define POLYNOMIAL_VAL_COMPUTE_X(j, outputVal, polynomial)\
    mpz_mul_si(outputVal, (polynomial->a), j );    \
    mpz_add(outputVal,outputVal,(polynomial->b));

/// compute (aj+b)^2-n
#define POLYNOMIAL_VAL_COMPUTE(j, outputVal, polynomial)\
    mpz_mul_si(outputVal, (polynomial->a), j );    \
    mpz_add(outputVal,outputVal,(polynomial->b)); \
    mpz_pow_ui(outputVal, outputVal, 2); \
    mpz_sub(outputVal,outputVal,(*(polynomial->N)));

typedef struct a_mongomeryPol_coeff{
    mpz_t* a;       //point to stored a coeff
    u_int* a_factors_indexes_FB;    //indexes of factors of a in precomputes->factorbase (also passed from master)
    u_int a_factors_num;                 //number of factors of a
} A_COEFF;
typedef struct Precomputes{ //precomputation for the SELF INIT of SIMPQS
    //factor base, prime up to B for witch N is a quadratic residue
    //because of limited size of B (for several reason) it's enough 64 bit integers
    DYNAMIC_VECTOR primes;
    u_int64_t* factorbase;        //all primes up to B that make N quadratic  residue
    u_int64_t factorbaseSize;
    //sqrt(N) mod p for each p
    //list with same inxeing of factobase
    mpz_t*   sqrtN_mod_p;       //tmem_p in contini phD
    mpz_t*   a_inv_mod_p;       //a^-1 mod p  TODO NB 0 SETTED ON NOT EXISTENT INVERSE (only on p s.t. p  is in a factorization)
    struct precomputatios_polynomial polPrecomputedData;
} PRECOMPUTES;                              //TODO MOVE ALL MALLOC ON START --> GCC OPTIMIZATION & EASIER FREE ON MALLOC'ed

/// worker constants params
typedef struct Configuration {
    mpz_t N;
    u_int64_t B;        //FactorBase threshold
    u_int64_t M;        //sieve array of size 2*M
    A_COEFF a_coefficient;
    /// memory configuration
    u_int64_t ARRAY_IN_MEMORY_MAX_SIZE;
    /// concurrency configuration
    int SIEVING_THREAD_NUM;
} CONFIGURATION;



struct polynomial{   //Mongomery polynomial to sieving
    mpz_t a;
    mpz_t b;
    mpz_t c;
    mpz_t*N;
};
#endif //SIMPQS_SIMPQS_H
