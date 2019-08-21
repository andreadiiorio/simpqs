//
// Created by andysnake on 19/08/19.
//

#include "SIMPQS.h"
#include "utils.h"
#include <stdio.h>


#define DIGITNUM_PRINT 5
void printPrecomputations(struct Precomputes* precomputes,int blockPrint){
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
    printf("\n sieve jumps <sol1_p,sol2_p> \n");
    for(__uint64_t i=0;i<sizeFB-blockPrint;i+=blockPrint){
        printf("%lu:\t\t",i);
        for(int j=0;j<blockPrint;j++){
            printf(" <%*lu,%*lu> ",DIGITNUM_PRINT,precomputes->polPrecomputedData.sol1p[i+j],DIGITNUM_PRINT,precomputes->polPrecomputedData.sol2p[i+j]);
        }
        printf("\n");
    }
    //TODO
//    precomputes->polPrecomputedData.B_ainv_2Bj_p
//    precomputes->polPrecomputedData.B_l
}