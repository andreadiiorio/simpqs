/*
 *  based on work of martani Jan 7, 2012
 *  enhanced  for distribuited SIMPQS by andysnake96
 */

#include <gmp.h>
#include <inttypes.h>
#include "worker/sievingSIMPQS.h"
#include "../utils/utils.h"


#ifndef MATRIX_H_
#define MATRIX_H_



typedef struct {
    mpz_t *MATRIX;              // bit rows TODO SAME INDEXING OF (partial) relations
    uint64_t rowsN;             //number of callocced rows  in MATRIX field
    mpz_t *IDENTITY;
	uint64_t colsN;
	uint64_t next_free_row;		/* used to insert rows in the matrix, point to the next free row to be inserted position */

} MATRIX;

/* allocates space for m rows * n columns matrix for MATRIX and IDENTITY */
void init_matrix(MATRIX *matrix, uint64_t rows, uint64_t cols);
int initMatrixFromRows(MATRIX *matrix, uint64_t rowsN,struct ArrayEntry* rowsRaw, uint64_t cols);

void push_row(MATRIX *matrix, mpz_t row);

void print_matrix_matrix(MATRIX *matrix);

void print_matrix_identity(MATRIX *matrix);

/* performs a Gauss elimination on matrix->MATRIX, result (linear dependence) will be in the matrix->IDENTITY */
void gauss_elimination(MATRIX *matrix);

/* does not check for bounds, the caller must */
void get_matrix_row(mpz_t rop, MATRIX *matrix, uint64_t row_index);

void get_identity_row(mpz_t rop, MATRIX *matrix, uint64_t row_index);

int quadraticRelationTry(REPORTS* reports,MATRIX* matrix);

#endif /* MATRIX_H_ */
