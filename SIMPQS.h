//
// Created by andysnake on 20/08/19.
//

#ifndef SIMPQS_SIMPQS_H
#define SIMPQS_SIMPQS_H

#include <stdlib.h>
#include <sys/types.h>
#include <gmp.h>

///precomputated vars
//TODO ALL WITH SAME INDEXING OF F.B. and mod p may be in uint64 see *
struct precomputatios_polynomial{   //precomputations polynomial dependent -> will change (slightly) at each polynomial change
    mpz_t* B_l;                 //generators of b by GreyCode; indexing 1<l<s where a=p1*..*ps
    //sieving with multiple polynomial vars
    //same indexing of FactorBase
    u_int64_t* sol1p;                 // sieving constant jumps
    u_int64_t* sol2p ;               // sieving constant jumps
    mpz_t* B_ainv_2Bj_p;            //2*Bj*(a^-1) mod p //TODO UPDATE JUMPs ON POLYNOMIAL CHANGE
};
struct Precomputes{ //precomputation for the SELF INIT of SIMPQS
    //factor base, prime up to B for witch N is a quadratic residue
    //because of limited size of B (for several reason) it's enough 64 bit integers
    u_int64_t factorbaseSize;
    u_int64_t* factorbase;        //all primes up to B that make N quadratic  residue
    //sqrt(N) mod p for each p
    //list with same inxeing of factobase
    mpz_t*   sqrtN_mod_p;       //tmem_p in contini phD //TODO * SHANKS TONELLI CODE WANT IN MPZ
    mpz_t*   a_inv_mod_p;       //a^-1 mod p        //TODO *avoid varius mpz extraction<->recomposition if it's in mpz
    struct precomputatios_polynomial polPrecomputedData;
} Precomputations;

/// worker constants params
struct Configuration {
    u_int64_t B;        //FactorBase threshold
    u_int64_t M;        //sieve array of size 2*M
    int s;          //num of coeff. a's factors TODO redundant
    mpz_t* a;       //sieving polynomial families identified by a coeff. passed from master
    double LOG_SIEVING_THREASHOLD;   //threshold that above it an array element A(j)=log(p1)+..+log(pn) is likely to be Bsmooth
    /// memory configuration
    u_int64_t ARRAY_IN_MEMORY_MAX_SIZE;
    /// concurrency configuration
    int SIEVING_THREAD_NUM;
};



struct polynomial_actual{   //actual polynomial in sieving
    mpz_t a;
    mpz_t b;
    mpz_t c;
};
extern struct polynomial_actual* ActualPolynomial;
#endif //SIMPQS_SIMPQS_H
