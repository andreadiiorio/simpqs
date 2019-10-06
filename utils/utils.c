//
// Created by andysnake on 19/08/19.
//

#include "utils.h"


#define DIGITNUM_PRINT 5
#define PRIMES_LIST_WORD_SIZE sizeof(int)
#define DEBUG_CHECK
void printSievingJumps(struct Precomputes* precomputes, int blockPrint) {
    printf("\n sieve jumps <sol1_p,sol2_p> \n");
    __uint64_t sizeFB=precomputes->factorbaseSize;
    for(__uint64_t i=0;i<sizeFB-blockPrint;i+=blockPrint){
        printf("%lu:\t\t",i);
        for(int j=0;j<blockPrint;j++){
            printf(" <%*ld,%*ld> ",DIGITNUM_PRINT,precomputes->polPrecomputedData.sol1p[i+j],DIGITNUM_PRINT,precomputes->polPrecomputedData.sol2p[i+j]);
        }
        printf("\n");
    }
}
int checkSieveJumps(PRECOMPUTES* precomputes,struct polynomial* polynomial) {
    mpz_t polVal,tmp;mpz_inits(polVal,tmp,NULL);
    u_int64_t prime,sieveJump;
    int result=EXIT_SUCCESS;
    for (u_int i = 0; i < precomputes->factorbaseSize; ++i) {
        prime = precomputes->factorbase[i];
        sieveJump=precomputes->polPrecomputedData.sol1p[i];
        POLYNOMIAL_VAL_COMPUTE(sieveJump,polVal,polynomial);
        if(!mpz_divisible_ui_p(polVal,prime)) {
            fprintf(stderr, "INVALID SIEVE JUMP AT INDEX %d\n", i);
            result = EXIT_FAILURE;
        }
        sieveJump=precomputes->polPrecomputedData.sol2p[i];
        POLYNOMIAL_VAL_COMPUTE(sieveJump,polVal,polynomial);
        if(!mpz_divisible_ui_p(polVal,prime)) {
            fprintf(stderr, "INVALID SIEVE JUMP AT INDEX %d\n", i);
            result = EXIT_FAILURE;
        }
    }
    return result;
}
void printPrecomputations(struct Precomputes* precomputes, int blockPrint){
    printf("\n precomputations \n");
    __uint64_t sizeFB=precomputes->factorbaseSize;
    printf("---factor Base p_i of %lu primes for witch legandre (n/p)=1---\n",sizeFB);
    for(__uint64_t i=0;i<sizeFB-blockPrint;i+=blockPrint){
        printf("%lu:\t\t",i);
        for(int j=0;j<blockPrint;j++){

            printf(" %*lu ",DIGITNUM_PRINT,precomputes->factorbase[i+j]);
        }
        printf("\n");
    }
    printf("\n---a^-1 mod p_i---\n");
    for(__uint64_t i=0;i<sizeFB-blockPrint;i+=blockPrint){
        printf("%lu:\t\t",i);
        for(int j=0;j<blockPrint;j++){
            gmp_printf(" %*Zd ",DIGITNUM_PRINT,precomputes->a_inv_mod_p[i+j]);
        }
        gmp_printf("\n");
    }

    printf("\n---sqrt(N) mod p_i---\n");
    for(__uint64_t i=0;i<sizeFB-blockPrint;i+=blockPrint){
        printf("%lu:\t\t",i);
        for(int j=0;j<blockPrint;j++){
            gmp_printf(" %*Zd ",DIGITNUM_PRINT,precomputes->sqrtN_mod_p[i+j]);
        }
        gmp_printf("\n");
    }
    printSievingJumps(precomputes,blockPrint);
    printf("\n---B_l vals ---\n");
    for (int k = 0; k < precomputes->polPrecomputedData.s; ++k) {
//        printf("B_%d\t",k);
        gmp_printf("%Zd\t",precomputes->polPrecomputedData.B_l[k]);
    }
    printf("\n\n---2*B*jainv_mod_p vals---\n");
    for (u_int k = 0; k < (precomputes->polPrecomputedData.s * precomputes->factorbaseSize)-blockPrint;k+=blockPrint) {
        for (int i = 0; i < blockPrint; ++i) {
            gmp_printf("%Zd\t",precomputes->polPrecomputedData.B_ainv_2Bj_p[k+i]);
        }
        gmp_printf("\n");
    }
}
#define FACTOR_BASE_BLOCK_REALLOC_N 1024

DYNAMIC_VECTOR* ReadPrimes(char *primesListPath, u_int64_t smoothnessBound) {

    DYNAMIC_VECTOR* outVect=malloc(sizeof(*outVect));
    if(!outVect) {
        fprintf(stderr, "Out of mem at primes allocation\n");
        return NULL;
    }
    outVect->pntr=NULL;outVect->vectorSize=0;
    //read primes list in file
    u_int64_t *result;
    FILE* primesListPrecomputed=fopen(primesListPath,"r");
    if(!primesListPrecomputed){
        perror(primesListPath);
        return outVect;
    }
    u_int64_t factorBaseDynamicN=FACTOR_BASE_BLOCK_REALLOC_N;
    u_int64_t* primes=malloc(factorBaseDynamicN * sizeof(*primes)); //factor base of primes up to bound that make N quadratic residue
    if(!primes){
        fprintf(stderr,"factor base malloc error\n");
        return outVect;
    }
    u_int64_t prime=2;
    u_int64_t primeIndx=0;
    primes[primeIndx++] = prime;              //TODO PATCH PRIME LIST MISSING 2
    size_t rd;
    while (prime<smoothnessBound) {
        if((rd = fread(&prime, PRIMES_LIST_WORD_SIZE, 1, primesListPrecomputed)) != 1) {
            printf("no more primes to read at %lu...\n", primeIndx);
        }
        /// evaluate reallocation
        REALLOC_WRAP(primeIndx,factorBaseDynamicN,primes,FACTOR_BASE_BLOCK_REALLOC_N)
                result = NULL;
                goto exit;
        }}
        primes[primeIndx++] = prime;              //save prime
        if (ferror(primesListPrecomputed)) {
            fprintf(stderr, "fread error occured in primes reading\t rd:%lu \n", rd);
            free(primes);
            result = NULL;
            goto exit;
        }
    }
#ifdef FINAL_RESIZE
    ////final realloc of primes list for actual space needed
    if(!(primes=realloc(primes,(--primeIndx)* sizeof(*primes)))){
        fprintf(stderr,"final realloc failed on factorbase\n");
        free(primes);
        return outVect;
    }
#endif
    result=primes;
#ifdef DEBUG
    fprintf(stderr,"primes addr %p\n",primes);
#endif
    exit:
    fclose(primesListPrecomputed);
    printf("correctly loaded %lu primes up to %lu \n",primeIndx,smoothnessBound);
    outVect->vectorSize=primeIndx;
    outVect->pntr=result;
    return outVect;
}

DYNAMIC_VECTOR ReadFactorBase(DYNAMIC_VECTOR primes, mpz_t N){
    DYNAMIC_VECTOR outVect=(DYNAMIC_VECTOR){.pntr=NULL,.vectorSize=0};
    size_t dynamicFactorBaseSize = FACTOR_BASE_BLOCK_REALLOC_N;
    u_int64_t* factorBase=malloc(sizeof(*factorBase) * dynamicFactorBaseSize);
    if (!factorBase)
        return outVect;
    u_int64_t  factorBaseIndx=0;
    mpz_t primeMpz;
    u_int64_t prime;
    mpz_init(primeMpz);
    u_int64_t* primesList=(u_int64_t *)primes.pntr;
    for(u_int64_t i=0;i<primes.vectorSize;i++){
        prime=primesList[i];
        mpz_set_ui(primeMpz,prime);
        if(mpz_legendre(N,primeMpz)==1) {            //prime ok for FactorBase
            REALLOC_WRAP(factorBaseIndx,dynamicFactorBaseSize,factorBase,FACTOR_BASE_BLOCK_REALLOC_N)
                    free(primes.pntr);
                    return outVect;
                }}
        ///set prime in factor base
        factorBase[factorBaseIndx++]=prime;
        }
    }
#ifdef FINAL_RESIZE
    ///final resize for actual space needed
    if(!(factorBase=realloc(factorBase,factorBaseIndx* sizeof(*factorBase)))){
        fprintf(stderr,"final realloc failed on factorbase\n");
        free(factorBase);
        free(primes.pntr);
        return outVect;
    }
#endif
    outVect.pntr=(void*)factorBase;
    outVect.vectorSize=factorBaseIndx;
    printf("factor base of %lu primes \n",factorBaseIndx);
    return outVect;
}

