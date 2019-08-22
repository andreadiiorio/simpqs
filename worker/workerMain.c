#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
#include <CONFIGURATION.h>
#include "SIMPQS.h"
#include "utils/utils.h"
#include "utils/gmp_patch.h"
#include "sievingSIMPQS.h"

mpz_t N;
SIEVE_ARRAY_BLOCK SieveArrayBlock;              ///array block in memory

///polinomial family coeff
struct polynomial_actual* ActualPolynomial;
struct Configuration Config;
int genFactorBase(mpz_t B, struct Precomputes* precomputes){
    //generate factor base of prime up to B that make n quadratic residue mod p
    //TODO FARE MEGLIO DEL MARIANO

    unsigned long int uBase = mpz_get_ui(B);
    int64_t nb_primes = sieve_primes((int64_t) (uBase));
    precomputes->factorbase=calloc(nb_primes, sizeof(int64_t));
    if(!(precomputes->factorbase))
        return EXIT_FAILURE;
    int nb_primes_quadr_res = fill_primes(Precomputations.factorbase, N);
    precomputes->factorbaseSize=nb_primes_quadr_res;
    return EXIT_SUCCESS;
}
int genSqrtN_modp(struct Precomputes *precomputation) {
    //generate sqrt(N) mod p where p is in FactorBase
    //Shanks-Tonelli algoritm used

    if(!(precomputation->sqrtN_mod_p=calloc(precomputation->factorbaseSize, sizeof(mpz_t))))
        return EXIT_FAILURE;
    mpz_t tmem_p,primeTmp;           //sqrt(N) mod prime
    mpz_init(tmem_p);
    mpz_init(primeTmp);
    int result=EXIT_SUCCESS;
    for(int64_t i=0; i<precomputation->factorbaseSize;i++){
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
        mpz_clear(tmem_p);
        return result;
}


int genSievingJumps(struct Precomputes *precomputations,const mpz_t* b) {
    //GEN sol1,2_p for each p in FB for the sieving phase for actual polynomial
    //allocate enough space for sol1,2_p
    if(!(precomputations->polPrecomputedData.sol1p=calloc(precomputations->factorbaseSize, sizeof(u_int64_t)))){
        fprintf(stderr,"calloc fail for sol1,p\n");
        return EXIT_FAILURE;
    }
    if(!(precomputations->polPrecomputedData.sol2p=calloc(precomputations->factorbaseSize, sizeof(u_int64_t)))){
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
            fprintf(stderr,"not existent inverse "); //TODO FORMATTED OUTPUT ADD
            gmp_fprintf(stderr,"%Zd mod %Zd \n",*a,primeTmp);
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

void mockConfiguration(){
    //TODO MOCKED CONFIGURATION SET  ==> next config file ?
    Config.ARRAY_IN_MEMORY_MAX_SIZE=_ARRAY_IN_MEMORY_MAX_SIZE;
    Config.SIEVING_THREAD_NUM= _SIEVING_THREAD_NUM;
}

void initializationVars(char** argv){
    /*
     generated Pol Coeff a: 16079919820099 b: 19000050618295 c: -6218915954954063976
 N: 100000030925519250968982645360649 B: 17330 M: 7330

     */
    if(!(ActualPolynomial=malloc(sizeof(struct polynomial_actual)))){
        fprintf(stderr,"malloc failed for actual polynomio instantiation");
        exit(EXIT_FAILURE);
    }

    //TODO MOCKED ARGV
    //char* _N=argv[1];
    //char* _M=argv[2];
    //char* _B=argv[3];
    //char* _a=argv[4];
    char* _N="100000030925519250968982645360649";
    char* _M="7330";
    char* _B="17330";
    char* _a="16079919820099";
    char* _b="19000050618295";
    char* _c="-6218915954954063976";
    mpz_init(N);
    strToMpz(N,_N)
    mpz_init(Config.M);
    strToMpz(Config.M,_M)
    mpz_init(Config.B);
    strToMpz(Config.B,_B)
    mpz_init((*ActualPolynomial).a);
    strToMpz((*ActualPolynomial).a,_a)
    mpz_init((*ActualPolynomial).b);
    strToMpz((*ActualPolynomial).b,_b)
    mpz_init((*ActualPolynomial).c);
    strToMpz((*ActualPolynomial).c,_c)
}
int preComputations(const mpz_t* a){
    //smart precomputations stored to reduce computational cost for each sieve iteration on each polynomial

    //FACTOR BASE GENERATION
    if(genFactorBase(Config.B,&Precomputations)==EXIT_FAILURE){
        fprintf(stderr,"PRECOMPUTATIONS ERROR ON FACTORBASE GEN\n");
        return EXIT_FAILURE;
    }
    //SQRT(N) MOD P GEN
    if(genSqrtN_modp(&Precomputations)==EXIT_FAILURE){
        fprintf(stderr,"PRECOMPUTATIONS ERROR ON SQRT(N) FOR EACH p IN FACTOR BASE\n");
        return EXIT_FAILURE;
    }

    //a^-1
    if(gen_a_inv_modp(a,&Precomputations)==EXIT_FAILURE){
        fprintf(stderr,"a^-1 PRECOMPUTATION ERROR");
        return EXIT_FAILURE;
    }

    //B_l GENERATION
    if(genB(&Config ,&Precomputations)==EXIT_FAILURE){
        fprintf(stderr,"error in B generation\n");
        return EXIT_FAILURE;
    }

    //soln1,2_p GENERATION ---> SIEVING JUMPS
    if(genSievingJumps(&Precomputations,&(ActualPolynomial.b))==EXIT_FAILURE){
        fprintf(stderr,"error in Sieving jumps precomputations \n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int main(int argc, char** argv){
    //FORKED WORKER PROCESS TO SIEVE BSMOOTH RELATION AND PARTIAL RELATION WITH POLINOMIO FAMILIES FROM a
    if(argc < 5){
        printf("USAGE: N,M,B,a, ....\n"); //todo append config
        //exit(EXIT_FAILURE);
    }
    initializationVars(argv);
    preComputations(&((*ActualPolynomial).a));
    gmp_printf("\n DONE precomputation for polynomial: a: %Zd b: %Zd factorizing N: %Zd \n",ActualPolynomial->a,ActualPolynomial->b,N);
    printPrecomputations(&Precomputations,5);

    /// sieve array allocate
    if(!(SieveArrayBlock=malloc(Config.ARRAY_IN_MEMORY_MAX_SIZE* sizeof(*SieveArrayBlock)))){   //TODO DEBUG ALLOCATED SIZE
        fprintf(stderr,"sieve array block in mem malloc failed\n");
        exit(EXIT_FAILURE);
    }
    if(Sieve(&Config,&Precomputations,SieveArrayBlock)==EXIT_FAILURE){
        fprintf(stderr,"sieve error \n");
        exit(EXIT_FAILURE);
    } //TODO ITERATE UNTIL MATRIX IS FILLED --> EVENTUALLY CHANGE POLYNOMIAL

    //TODO MATRIX STAGE CALLBACK TO MASTER
}
