//
// Created by andysnake on 19/08/19.
//

#ifndef QUADRATIC_SIEVE_UTILS_H
#define QUADRATIC_SIEVE_UTILS_H

//initialize multiple precision int from sourceString in dest
#define strToMpz(dest,sourceStr)\
    if (mpz_init_set_str(dest,sourceStr, 10) == -1) {\
    printf("Cannot load N %s\n", argv[1]);\
    exit(2);\
}

#define LONG_MUL_OVERFLOW_P(a, b) \
   __builtin_mul_overflow_p (a, b, (__typeof__ ((a) * (b))) 0)

void printPrecomputations(struct Precomputes* precomputes,int blockPrint);
#endif //QUADRATIC_SIEVE_UTILS_H
