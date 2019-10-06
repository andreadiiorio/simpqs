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
    mpz_t *MATRIX;              // exponent vector rows ( in all QS variant values in GF2, sum == xor -> huge HardWare accelaration)
    uint64_t rowsN;             //number of exponent vector rows (same num of full reports)
    mpz_t *IDENTITY;            //Identity square matrix rowsN x rowsN
    /*
     * it hold that: full report i among all reports is identified in Identity matrix with a 1 in col i
     * before gauss step it hold also: full report i -> exp vector in MATRIX at i -> 1 at [i][i] of IDENTITY matrix
     * after  gauss step it hold: IDENTITY row  r identify summed rows in start MATRIX ( same indexing of source reports )
     * to produce final row r in MATRIX, simply summing related rows watching at 1s in IDENTITY[r]
     * if M' is MATRIX after gauss elimination ->  M'[i]=I'[i]*M  TODO CHECK OR DELETE LINE**
     */
	uint64_t colsN;
	uint64_t next_free_row;		/* used to insert rows in the matrix, point to the next free row to be inserted position */

} MATRIX;

void init_matrix(MATRIX *matrix, uint64_t rows, uint64_t cols);
int initMatrixFromRows(MATRIX *matrix, uint64_t rowsN,struct ArrayEntry* rowsRaw, uint64_t cols);

void gauss_elimination_over_GF2(MATRIX *matrix);

/* does not check for bounds, the caller must */

int quadraticRelationTry(REPORTS* reports,MATRIX* matrix);

#endif /* MATRIX_H_ */
