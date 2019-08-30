#include "utils.h"
#include "gmp_patch.h"

mpz_t N;
struct Configuration  Config;

int genSqrtN_modp(struct Precomputes *precomputation) {
    //generate sqrt(N) mod p where p is in FactorBase
    //Shanks-Tonelli algoritm used

    if(!(precomputation->sqrtN_mod_p=calloc(precomputation->factorbaseSize, sizeof(mpz_t))))
        return EXIT_FAILURE;
    mpz_t tmem_p,primeTmp;           //sqrt(N) mod prime
    mpz_init(tmem_p);
    mpz_init(primeTmp);
    int result=EXIT_SUCCESS;
    for(u_int64_t i=0; i<precomputation->factorbaseSize;i++){
        mpz_set_ui(primeTmp, precomputation->factorbase[i]);
#ifdef LEGANDRE_EXTRA_CHECK_SQRTM
        if(mpz_sqrtm(tmem_p,N,primeTmp)==EXIT_FAILURE){
            result = EXIT_FAILURE;
            goto exit;
        }
#else
        mpz_sqrtm(tmem_p,N,primeTmp);
#endif
        mpz_set(precomputation->sqrtN_mod_p[i],tmem_p); //set the sqrt of N mod p in his field
    }
    exit:
    mpz_clears(tmem_p,primeTmp,NULL);
    return result;
}
int genSievingJumps(struct Precomputes *precomputations,const mpz_t* b) {
    //GEN sol1,2_p for each p in FB for the sieving phase for actual polynomial
    //allocate enough space for sol1,2_p
    if(!(precomputations->polPrecomputedData.sol1p=calloc(precomputations->factorbaseSize, sizeof(int64_t)))){
        fprintf(stderr,"calloc fail for sol1,p\n");
        return EXIT_FAILURE;
    }
    if(!(precomputations->polPrecomputedData.sol2p=calloc(precomputations->factorbaseSize, sizeof(int64_t)))){
        fprintf(stderr,"calloc fail for sol2,p\n");
        return EXIT_FAILURE;
    }
    //pointers to p_th sieving jumps in precomputations
    mpz_t sol1_p_tmp, sol2_p_tmp,sub,add;
    mpz_inits(sol1_p_tmp,sol2_p_tmp,sub,add,NULL);
    u_int64_t p;
    for(u_int64_t i=0;i<precomputations->factorbaseSize;i++){
        p = precomputations->factorbase[i];
        mpz_sub(sub,precomputations->sqrtN_mod_p[i],*b);    //tmem_p-b
        mpz_add(add,precomputations->sqrtN_mod_p[i],*b);
        mpz_neg(add,add);                                   //-tmem_p-b
        //mul
        mpz_mul(sol1_p_tmp,precomputations->a_inv_mod_p[i],sub);
        mpz_mul(sol2_p_tmp,precomputations->a_inv_mod_p[i],add);
        //reduce mod p; saving in 64bit integer
        precomputations->polPrecomputedData.sol1p[i]=mpz_mod_ui(sol1_p_tmp,sol1_p_tmp,p);
        precomputations->polPrecomputedData.sol2p[i]=mpz_mod_ui(sol2_p_tmp,sol2_p_tmp,p);
        //save in precomputations
    }
    mpz_clears(sol1_p_tmp,sol2_p_tmp,add,sub,NULL);
    return EXIT_SUCCESS;
}
int gen_a_inv_modp(const mpz_t* a, struct Precomputes *precomputations) {
    //generate a ^-1 mod p for each p in FB
    //for each factor a will not exist inverse mod factor --> will be setted 0
    if(!(precomputations->a_inv_mod_p =calloc(precomputations->factorbaseSize, sizeof(mpz_t)))){
        fprintf(stderr,"Out of mem a^inv precomputes\n");
        return EXIT_FAILURE;
    }
    //tmp i_th prime in FB
    mpz_t primeTmp;
    mpz_init(primeTmp);
    for(u_int64_t i=0;i<precomputations->factorbaseSize;i++){
        mpz_init(precomputations->a_inv_mod_p[i]);
        mpz_set_ui(primeTmp,precomputations->factorbase[i]);
        if(!mpz_invert(precomputations->a_inv_mod_p[i],*a,primeTmp)){
            gmp_fprintf(stderr,"not existent inverse %Zd mod %Zd \n",*a,primeTmp);
            mpz_set_ui(precomputations->a_inv_mod_p[i],0);
//            mpz_clear(primeTmp);
//            return EXIT_FAILURE;
        }
    }

    mpz_clear(primeTmp);
    return EXIT_SUCCESS;
}
int genB(struct Configuration *config, struct Precomputes *precomputations) {
    //TODO GEN
    //  B_l
    //  2*B*ainv_p
    return EXIT_SUCCESS;
}
struct Precomputes* preComputations(const mpz_t* a){
    //smart precomputations stored to reduce computational cost for each sieve iteration on each polynomial
    struct Precomputes* precomputations;
    if(!(precomputations=malloc(sizeof(*precomputations)))){
        fprintf(stderr,"malloc fail on precomputes\n");
        return NULL;
    }
    DYNAMIC_VECTOR primes=ReadPrimes(PRIMES_32B_PATH,Config.B);
    if(!primes.pntr){
        fprintf(stderr,"primes read error\n");
        return NULL;
    }
    DYNAMIC_VECTOR factorBase=ReadFactorBase(primes,N);
    if(!factorBase.pntr){
        fprintf(stderr,"factor base gen error\n");
        free(primes.pntr);
        return NULL;
    }
    precomputations->primes=primes;
    precomputations->factorbase=factorBase.pntr; precomputations->factorbaseSize=factorBase.vectorSize;
    //SQRT(N) MOD P GEN
    if(genSqrtN_modp(precomputations)==EXIT_FAILURE){
        fprintf(stderr,"PRECOMPUTATIONS ERROR ON SQRT(N) FOR EACH p IN FACTOR BASE\n");
        free(primes.pntr);free(factorBase.pntr);
        return NULL;
    }
    printf("sqrt of N mod p in primes generated\n");
    //a^-1
    if(gen_a_inv_modp(a,precomputations)==EXIT_FAILURE){
        fprintf(stderr,"a^-1 PRECOMPUTATION ERROR");
        free(primes.pntr);free(factorBase.pntr);
        return NULL;
    }
    printf("a^-1 mod p in primes generated\n");

    //B_l GENERATION
    if(genB(&Config ,precomputations)==EXIT_FAILURE){
        fprintf(stderr,"error in B generation\n");
        free(primes.pntr);free(factorBase.pntr);
        return NULL;
    }

    //soln1,2_p GENERATION ---> SIEVING JUMPS
    if(genSievingJumps(precomputations,&(ActualPolynomial->b))==EXIT_FAILURE){
        fprintf(stderr,"error in Sieving jumps precomputations \n");
        free(primes.pntr);free(factorBase.pntr);
        return NULL;
    }
    printf("sieving jumps computed\n");
    return precomputations;
}

