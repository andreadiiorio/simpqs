#include "utils.h"
#include "gmp_patch.h"

struct Configuration  Config;
#define  MIN_FACTOR_A_COEFF 2069 //8 primes after 2000
//#define  GEN_A_UNDER_THRESHOLD


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
        if(mpz_sqrtm(tmem_p,Config.N,primeTmp)==EXIT_FAILURE){
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
int gen_a_inv_modp(A_COEFF* a, struct Precomputes *precomputations) {
    //generate a ^-1 mod p and 2Bl*(a^-1)modp for each p in FB
    //NB: for each factor a will not exist inverse mod factor --> will be setted 0
    if(!(precomputations->a_inv_mod_p =calloc(precomputations->factorbaseSize, sizeof(mpz_t)))){
        fprintf(stderr,"Out of mem a^inv precomputes\n");
        return EXIT_FAILURE;
    }
    //tmp i_th prime in FB
    mpz_t primeTmp;
    mpz_t* a_inv_p_i;                   //point to location for (a^-1)mod p
    mpz_t* B_j;                         // point to B_j location in precomputes
    mpz_t* B_j2a_inv;                     //point to location holding 2bjainv
    if(!(precomputations->polPrecomputedData.B_ainv_2Bj_p=malloc(a->a_factors_num*precomputations->factorbaseSize* sizeof(mpz_t)))){
        fprintf(stderr,"out of mem at B_ainv_2Bj mod p matrix\n");
        free(precomputations->a_inv_mod_p);
        return EXIT_FAILURE;
    }
    mpz_init(primeTmp);
    for(u_int64_t p=0; p < precomputations->factorbaseSize; p++){
        ///compute a^-1 mod p
        a_inv_p_i=&(precomputations->a_inv_mod_p[p]);
        mpz_init(*a_inv_p_i);
        mpz_set_ui(primeTmp,precomputations->factorbase[p]);
        if(!mpz_invert(*a_inv_p_i,*(a->a),primeTmp)){
            gmp_fprintf(stderr,"not existent inverse %Zd mod %Zd \n",*a,primeTmp);
            mpz_set_ui(precomputations->a_inv_mod_p[p], 0);
            continue;
        }
        ///compute Bj*2*ainv
        for (u_int j = 0; j < a->a_factors_num; ++j) {
            B_j=&(precomputations->polPrecomputedData.B_l[j]);
            B_j2a_inv=&(precomputations->polPrecomputedData.B_ainv_2Bj_p[j+ (p*(a->a_factors_num))]); //take pointer to dest location of 2Bkainv indexing by a factor and prime
            mpz_init_set(*B_j2a_inv,*B_j);
            mpz_mul(*B_j2a_inv,*B_j2a_inv,*a_inv_p_i);
            mpz_mul_ui(*B_j2a_inv,*B_j2a_inv,2);    //computed 2*B_j*a_inv on correct location
            mpz_mod(*B_j2a_inv,*B_j2a_inv,primeTmp);
        }
    }
    mpz_clear(primeTmp);
    return EXIT_SUCCESS;
}
mpz_t* genB(struct Configuration *config, struct Precomputes *precomputations) {
    /// gen B_l  and first b coefficient for first polynomial
    //multiply computed B_l in first b coefficient for first polynomial

    mpz_t* b_first=malloc(sizeof(*b_first));               //first b coeff returned
    if(!b_first || !(precomputations->polPrecomputedData.B_l=calloc(config->a_coefficient.a_factors_num, sizeof(*(precomputations->polPrecomputedData.B_l))))){
        fprintf(stderr,"B_l CALLOC FAILED \n");
        free(precomputations->polPrecomputedData.B_l );free(b_first);
        return NULL;
    }
    precomputations->polPrecomputedData.s=config->a_coefficient.a_factors_num;
    mpz_init_set_ui(*b_first,0);    //initialize first b coff with 1 (sucessive multiplication will cumulate B_ls)
    u_int64_t  ql;                  //l_th prime factor of a
    u_int64_t  ql_indx;             //index of ql in factor base array
    mpz_t*      t_mem_ql;           //sqrt of N mod ql
    mpz_t*      B_l;                //point to allocated mpz_t for B_l
    mpz_t       gamma,ql_tmp,a_ql; mpz_inits(gamma,ql_tmp,a_ql,NULL);
    void* result=b_first;
    for (int l = 0; l < precomputations->polPrecomputedData.s; l++) {
        //GET l th prime factor of a using a factors indexes in factor base
        ql_indx=config->a_coefficient.a_factors_indexes_FB[l];            //index of prime factor of a ql and sqrt N mod ql
        ql=precomputations->factorbase[ql_indx];            //prime factor of a
        t_mem_ql=&(precomputations->sqrtN_mod_p[ql_indx]); // sqrt of N mod ql has same indexing of factor base
        // compute (a/ql)^-1 mod ql
        mpz_set_ui(ql_tmp,ql);
        mpz_divexact_ui(a_ql,*(config->a_coefficient.a),ql);                  // a/ql
        if(!mpz_invert(gamma,a_ql,ql_tmp)){                //(a/ql)^-1 mod ql
            gmp_fprintf(stderr,"INVALID INVERSION OF GAMMA %Zd mod %Zd\n",gamma,ql_tmp);
            result=NULL;
            goto exit;
        }
        mpz_mul(gamma,gamma,*(t_mem_ql));
        mpz_mod_ui(gamma,gamma,ql);                     //finally gamma= tmem_p*(a/ql)^-1 mod ql
        if(mpz_cmp_ui(gamma,(ql/2))>0) {
            mpz_sub(gamma, ql_tmp, gamma);                    // on gamma >ql/2 replace with ql-gamma
        }
        B_l=&(precomputations->polPrecomputedData.B_l[l]);
        mpz_init(*B_l);
        mpz_mul(*B_l,a_ql,gamma);                            //finally get B_l
        //cumulate computed Bl in first B (gray code 111...11)
        gmp_printf("gamma :%Zd \t a_ql:%Zd\tql:%lu B_l:%Zd\n",gamma,a_ql,ql,*B_l);fflush(0);
        mpz_add(*b_first,*b_first,*B_l);

    }

    exit:
        mpz_clears(ql_tmp,gamma,a_ql,NULL);
        if(result==NULL){
            mpz_clear(*b_first);free(b_first);  // on non memory every fault free malloced b
        }
        return result;
}


CONFIGURATION* initConfiguration(const char* n, int arrayInMemMaxSize, int64_t M, u_int64_t B, int sieverThreadsNum) {
    //init main SIMPQS paramters with passed paramenters if non zero values has been set
    //on parameter passed at 0 default value will be used
    //returned heap allocated config struct, global config also setted
    CONFIGURATION* config=calloc(1, sizeof(*config));
    if(!config){
        fprintf(stderr,"out of mem in config calloc\n");
        return NULL;
    }
    strToMpz(config->N,n);
    config->ARRAY_IN_MEMORY_MAX_SIZE=_ARRAY_IN_MEMORY_MAX_SIZE;
    if(arrayInMemMaxSize!=0)
        config->ARRAY_IN_MEMORY_MAX_SIZE=arrayInMemMaxSize;
    config->M=_M;
    if(M!=0)
        config->M=M;
    config->B=_B;
    if(B!=0)
        config->B=B;
    config->SIEVING_THREAD_NUM=_SIEVING_THREAD_NUM;
    if (sieverThreadsNum!=0)
        config->SIEVING_THREAD_NUM=sieverThreadsNum;
    Config=*config;
    return config;
}
A_COEFF* gen_a_centered(const u_int64_t* factorbase, u_int64_t factorBaseSize, int s,struct Configuration *configuration){
    //generate a coefficient s.t. a=p1* ... * ps=~(Sqrt(2*N)/M)
    //simply aggregate factors until ideall value reached, aggregation factors ceneterd at ideal^(1/s)
    //a coeff setted in configuration
    A_COEFF* aCoeff=calloc(1, sizeof(*aCoeff));
    mpz_t*a=calloc(1, sizeof(*a));
    u_int* aFactorsIndexes=calloc(s, sizeof(*aFactorsIndexes));
    if(!a || ! aCoeff || !aFactorsIndexes){
        fprintf(stderr,"OUT OF MEM AT a Malloc\n");
        free(a);free(aCoeff);
        return NULL;
    }
    u_int a_factors_num=0;
    mpz_init_set_ui(*a,1);
    /// get a ideally value ~ SQRT(2N)/M
    mpz_t ideal_a_value;
    mpz_init_set(ideal_a_value,configuration->N);
    mpz_mul_ui(ideal_a_value,ideal_a_value,2);
    mpz_sqrt(ideal_a_value,ideal_a_value);
    mpz_div_ui(ideal_a_value,ideal_a_value,configuration->M);
    gmp_printf("ideal value of a is about: %Zd\n",ideal_a_value);

    /// get a factors center
    mpz_t a_factors_center; mpz_init_set(a_factors_center,ideal_a_value);
    mpz_root(a_factors_center,a_factors_center,s);
    if (mpz_cmp_ui(a_factors_center,MIN_FACTOR_A_COEFF)<0) {
        gmp_printf("too little a factor center achived by (sqrt(2N)/M)^(1/s) : %Zd... rounding up to :%d\n",a_factors_center,MIN_FACTOR_A_COEFF);
        mpz_set_ui(a_factors_center, MIN_FACTOR_A_COEFF);
    }
    //get factor center index
    u_int64_t a_factor_center=mpz_get_ui(a_factors_center);
    u_int64_t a_factor_center_indx;
    for (a_factor_center_indx=0; a_factor_center_indx < factorBaseSize ; ++a_factor_center_indx) {
        if(a_factor_center-1<factorbase[a_factor_center_indx])
            break;
    }
    //aggregate factor until reached threashold
    mpz_t a_tmp;mpz_init_set(a_tmp,*a);

    // aggregate a factors swinging around a_factor_center_indx until the ideal threshold has been reached
    u_int64_t prime,primeIndx, swing,i;
    for(primeIndx = a_factor_center_indx,swing=0,i=1; mpz_cmp(a_tmp, ideal_a_value) < 0 ;swing++) { //centered aggregation of factors
        prime=factorbase[primeIndx];
        mpz_mul_ui(a_tmp, a_tmp, prime);
        aFactorsIndexes[a_factors_num++]=primeIndx;
        primeIndx=((swing)%2 ? (a_factor_center_indx+ (i++)):(a_factor_center_indx-i));
    }
#ifdef GEN_A_UNDER_THRESHOLD
    mpz_div_ui(a_tmp,a_tmp,prime);              //remove last prime for a gen under threshold ( last prime caused for condition become false)
    aFactorsIndexes[a_factors_num--]=0;
#endif
    mpz_set(*a,a_tmp);                          //finally set computed a coeff
    gmp_printf("generated a coeff with s:%d factors\t a:%Zd\n",s,*a);

    aCoeff->a=a;
    aCoeff->a_factors_num=a_factors_num;
    aCoeff->a_factors_indexes_FB=aFactorsIndexes;
    configuration->a_coefficient=*aCoeff;

    mpz_clears(a_tmp,ideal_a_value,a_factors_center,NULL);
    return aCoeff;
}
void nextPolynomial_b_i(mpz_t* b, unsigned int i, PRECOMPUTES *precomputes){
    //gen b i+1
    ///get v s.t. 2^v|| 2i --> maximal power of 2 dividing 2i
    unsigned int v,i_2v;
//    for (v=1;  !(i & (1 << (v-1))); v++) ;
    for (v=0; 0 == ((2 * i) % (1 << v)); v++) { //starting from 2^0 find first power of 2 not diving 2i -> prev is v value
    } v--;  //for definition of 2^v || 2i
    i_2v= (i % (1<<v)) ? (i/ (1<<v))+1 : (i/ (1<<v));               //integral ceil [i/2^v]
    int signExp=(i_2v%2)?-1:1;

    /// compute next b i+1
    v--;                                        //C indexing  v in [0,s)
    mpz_t* B_v=&(precomputes->polPrecomputedData.B_l[v]);
    mpz_t *Bainv2vp_pntr;
    mpz_t tmp; mpz_init_set_si(tmp,2*signExp);
    mpz_mul(tmp,tmp,*B_v);                                  //tmp = B_v*2*(-1)^(ceil(i/2<<v))
    mpz_add(*b,*b,tmp);
#ifdef VERBOSE
    gmp_printf("\ngenerated next b i:%d with v:%d\t b:%Zd\n",i,v,*b);Z
#endif
    int64_t Bainv2vp;
    int64_t prime;
    /// update sieve jumps
    printf("updating sieving jumps, signExp:%d\n",signExp);
    for (u_int64_t p = 0; p < precomputes->factorbaseSize; ++p) {
        prime=precomputes->factorbase[p];
        //// SKIP A FACTORS IN SIEVE JUMPS UPDATE TODO BETTER WAY WITHOUT A COEFF PASSING EVERYWHERE --> copy in precomputations?
        if(!mpz_cmp_ui(precomputes->a_inv_mod_p[p],0))
            continue;
        Bainv2vp_pntr = &(precomputes->polPrecomputedData.B_ainv_2Bj_p[v + p * precomputes->polPrecomputedData.s]);
        Bainv2vp = mpz_get_si(*Bainv2vp_pntr);
        precomputes->polPrecomputedData.sol1p[p]-=(Bainv2vp*signExp);
        precomputes->polPrecomputedData.sol1p[p]=MOD(precomputes->polPrecomputedData.sol1p[p],prime);
        precomputes->polPrecomputedData.sol2p[p]-=(Bainv2vp*signExp);
        precomputes->polPrecomputedData.sol2p[p]=MOD(precomputes->polPrecomputedData.sol2p[p],prime);
    }
}
struct Precomputes *preComputations(CONFIGURATION *configuration, struct polynomial *dstPolynomial, A_COEFF *aCoeff) {
    //smart precomputations stored to reduce computational cost for each sieve iteration on each polynomial
    //if aCoeff is NULL a coeff will be computed otherwise that value will be used
    struct Precomputes* precomputations;
    if(!(precomputations=malloc(sizeof(*precomputations)))){
        fprintf(stderr,"malloc fail on precomputes\n");
        return NULL;
    }

    struct Precomputes* result=precomputations;
    DYNAMIC_VECTOR primes=ReadPrimes(PRIMES_32B_PATH, configuration->B);
    if(!primes.pntr){
        fprintf(stderr,"primes read error\n");
        return NULL;
    }
    DYNAMIC_VECTOR factorBase=ReadFactorBase(primes, configuration->N);
    if(!factorBase.pntr){
        fprintf(stderr,"factor base gen error\n");
        result=NULL; goto exit;
    }
    precomputations->primes=primes;
    precomputations->factorbase=factorBase.pntr; precomputations->factorbaseSize=factorBase.vectorSize;


    if (!aCoeff) aCoeff =gen_a_centered(factorBase.pntr,factorBase.vectorSize,4,configuration);
    //SQRT(N) MOD P GEN
    if(genSqrtN_modp(precomputations)==EXIT_FAILURE){
        fprintf(stderr,"PRECOMPUTATIONS ERROR ON SQRT(N) FOR EACH p IN FACTOR BASE\n");
        result=NULL; goto exit;
    }
    printf("sqrt of N mod p in primes generated\n");

    ///B_l GENERATION
    mpz_t* b_first=genB(configuration , precomputations);
    if(!b_first){
        fprintf(stderr,"error in B generation\n");
        result=NULL; goto exit;
    }

    gmp_printf("B_l generated, first b coeff: %Zd\n",*b_first);

    ///a^-1 && 2Bl(a^-1)
    if(gen_a_inv_modp(&(configuration->a_coefficient),precomputations)==EXIT_FAILURE){
        fprintf(stderr,"a^-1 PRECOMPUTATION ERROR");
        result=NULL; goto exit;
    }
    printf("a^-1 mod p in primes generated\n");
    //soln1,2_p GENERATION ---> SIEVING JUMPS
    if(genSievingJumps(precomputations,b_first)==EXIT_FAILURE){
        fprintf(stderr,"error in Sieving jumps precomputations \n");
        result=NULL; goto exit;
    }
    printf("sieving jumps computed\n");

    ///init first polynomial with computed a and b
    mpz_init_set(dstPolynomial->a,*(aCoeff->a));
    configuration->a_coefficient=*aCoeff;               //extra shallow copy of a coefficient
    mpz_init_set(dstPolynomial->b,*b_first);
    dstPolynomial->N=&configuration->N;
    exit:
        if(result==NULL){       //null set equally to same failure happened
            free(primes.pntr);free(factorBase.pntr);free(precomputations);
        }
        return result;
}

#ifdef TEST_POLYNOMIALS
int main(){
#else
void main_(){
#endif
    const char* n_str="100000030925519250968982645360649";
    CONFIGURATION *configuration = initConfiguration(n_str, 0, 0, 0, 0);
    //// get first polynomial:
    struct polynomial pol;
    PRECOMPUTES *precomputes = preComputations(configuration, &pol,NULL);
    printPrecomputations(precomputes, 10);
    gmp_printf("first polynomial:\ta:%Zd\tb:%Zd\n", pol.a, pol.b);
    printSievingJumps(precomputes, 10);
    //// polynomials family generation:
    for (u_int j = 1; j + 1 < 1u << (configuration->a_coefficient.a_factors_num - 1); ++j) {
        nextPolynomial_b_i(&(pol.b), j, precomputes);
        gmp_printf("polynomial:%d\ta=%Zd;\tb=%Zd;\n", j, pol.a, pol.b);
        printSievingJumps(precomputes, 10);
        if (checkSieveJumps(precomputes, &pol) == EXIT_FAILURE) {
            exit(EXIT_FAILURE);
        }
    }
}