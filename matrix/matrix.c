/*
 *  based on work of martani Jan 7, 2012
 *  enhanced to SIMPQS by andysnake96
 */
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <inttypes.h>
#include <worker/sievingSIMPQS.h>
#include "matrix.h"
#include "../utils/utils.h"

void initMatrixFromRows(MATRIX *matrix, uint64_t rowsN,mpz_t* rows, uint64_t cols) {
    matrix->MATRIX = rows;
    matrix->rowsN = rowsN;
    matrix->colsN = cols;
    matrix->next_free_row = rowsN+1;
}

void init_matrix(MATRIX *matrix, uint64_t rows, uint64_t cols) {
	matrix->MATRIX = (mpz_t *) calloc(rows, sizeof(mpz_t));
	matrix->rowsN = rows;
	matrix->colsN = cols;
	matrix->next_free_row = 0;
}
int setIdentityMatrix(MATRIX* matrix){
    //alloc & set identity matrix inside matrix
    if (!(matrix->IDENTITY= (mpz_t *) calloc(matrix->rowsN, sizeof(*(matrix->MATRIX))))){
        return EXIT_FAILURE;
    }
    for (uint64_t i = 0; i < matrix->colsN; ++i) {
        mpz_setbit(matrix->IDENTITY[i],i);              //set [i][i] in identity matrix
    }
    return EXIT_SUCCESS;
}

int extendMatrixByEntries(MATRIX* matrix,struct ArrayEntry* rowsEntries,uint64_t rowsNum){
    //extend matrix with rowsNum rows
    //realloc matrix rows array
    uint64_t newMatrixRowsN= matrix->rowsN + rowsNum;
    matrix->rowsN=newMatrixRowsN;
    if (!(matrix->MATRIX=realloc(matrix->MATRIX,newMatrixRowsN* sizeof(*(matrix->MATRIX))))){
        return EXIT_FAILURE;
    }
    //copy array values into matrix
    for (uint64_t i=0; i < rowsNum; i++, matrix->next_free_row++) {
        mpz_set(matrix->MATRIX[matrix->next_free_row],(rowsEntries[i]).exp_vector);                             //copy row to matrix
//        mpz_clear((rowsEntries[i]).exp_vector);
        // TODO SMART CHECKS
        // -DICT CUSTOMIZED FOR REDUNDANT ROWs
        // -0s ROW CHECK --> ALREADY FOUNDED QUADRATIC RELATION :==)))
    }
    return EXIT_SUCCESS;
}
int extendMatrixByRows(MATRIX* matrix,mpz_t* rows,uint64_t rowsNum){
    //extend matrix with rowsNum rows
    //realloc matrix rows array
    uint64_t newMatrixRowsN= matrix->rowsN + rowsNum;
    matrix->rowsN=newMatrixRowsN;
    if (!(matrix->MATRIX=realloc(matrix->MATRIX,newMatrixRowsN* sizeof(*(matrix->MATRIX))))){
        return EXIT_FAILURE;
    }
    //copy array values into matrix
    for (uint64_t i=0; i < rowsNum; i++, matrix->next_free_row++) {
        mpz_set(matrix->MATRIX[matrix->next_free_row],rows[i]);                             //copy row to matrix
        mpz_clear(rows[i]);                                             //clear source row pointer

        // TODO SMART CHECKS
        // -DICT CUSTOMIZED FOR REDUNDANT ROWs
        // -0s ROW CHECK --> ALREADY FOUNDED QUADRATIC RELATION :==)))
    }
    return EXIT_SUCCESS;
}

int concatMatrixes(MATRIX* matrixDst, MATRIX* matrix2){
    //concat MatrixDst with rows of matrix2 coping mpz_t rows

    uint64_t newMatrixRowsN= matrixDst->rowsN + matrix2->rowsN;
    matrixDst->rowsN=newMatrixRowsN;
    if (!(matrixDst->MATRIX=realloc(matrixDst->MATRIX, newMatrixRowsN * sizeof(*(matrixDst->MATRIX))))){
        return EXIT_FAILURE;
    }
    //copy array values into matrix1
    for (uint64_t i=0; i < matrix2->rowsN; i++, matrixDst->next_free_row++) {
        mpz_set(matrixDst->MATRIX[matrixDst->next_free_row], matrix2->MATRIX[i]);                             //copy row to matrix1

        // TODO SMART CHECKS
        // -DICT CUSTOMIZED FOR REDUNDANT ROWs
        // -0s ROW CHECK --> ALREADY FOUNDED QUADRATIC RELATION :==)))
    }
    return EXIT_SUCCESS;
}
void push_row(MATRIX *matrix, mpz_t row) {

	mpz_set(matrix->MATRIX[matrix->next_free_row], row);
	mpz_init2(matrix->IDENTITY[matrix->next_free_row], matrix->colsN); /* initializes a n bit vector all set to 0 */
	mpz_setbit(matrix->IDENTITY[matrix->next_free_row], matrix->next_free_row); /* set the next_free_row bit to 1 */
	matrix->next_free_row++;
}

void print_matrix_matrix(MATRIX *matrix) {
	uint64_t i, j;

	printf("\nMATRIX\n");
	for (i = 0; i < matrix->rowsN; i++) {
		printf("[");
		for (j = 0; j < matrix->colsN; j++) {
			printf("%2d", mpz_tstbit(matrix->MATRIX[i], j));
		}
		printf(" ]\n");
	}
	printf("\n");
}

void print_matrix_identity(MATRIX *matrix) {
	uint64_t i, j;

	printf("\nIDENTITY\n");
	for (i = 0; i < matrix->rowsN; i++) {
		printf("[");
		for (j = 0; j < matrix->colsN; j++) {
			printf("%2d", mpz_tstbit(matrix->IDENTITY[i], j));
		}
		printf(" ]\n");
	}
	printf("\n");
}

void free_matrix(MATRIX *matrix) {
	free(matrix->MATRIX);
	free(matrix->IDENTITY);
}

/* performs a Gauss elimination on matrix->MATRIX, M:I---->M':I'
 * result (linear dependence) will be in the matrix->IDENTITY
 *   TODO BECAUSE EACH OP TO MAKE MATRIX IN ROW ECHELON FORM it has been applied on Identity Matrix too
 *      so M'[i]=I'[i]*M  **
 *      target so it's to find linear dependencies (0 vectors in M') and re build in terms of source terms **
 * */
void gauss_elimination(MATRIX *matrix) {
	printf("\nPerforming Gauss elimination..\n");
	mpz_t *m = matrix->MATRIX;
	mpz_t *I = matrix->IDENTITY;
	uint64_t col, row, next_row;
	int64_t next_pivot;
	for (next_row = 0, col = 0; col < MIN(matrix->colsN, matrix->rowsN); col++) /* for all rows*/
	{
		next_pivot = -1;
		for (row = next_row; row < matrix->rowsN; row++) /* search for the next pivot*/
		{
			if (mpz_tstbit(m[row], col)) {
				next_pivot = row; /* row contains the next pivot */
				next_row++;
				break;
			}
		}

		if (next_pivot == -1)
			continue;

		if (next_pivot != next_row - 1) /* current row is not the pivot, switch rows */
		{
		    printf("switch rows...\n");fflush(0);
			mpz_swap(m[next_pivot], m[next_row - 1]);
			mpz_swap(I[next_pivot], I[next_row - 1]);
		}

		for (row = next_row; row < matrix->rowsN; row++) {
			if (mpz_tstbit(m[row], col)) {
				mpz_xor(m[row], m[row], m[next_row - 1]); /* XOR the rows to eliminate the 1 in position (row, next_row-1)*/
				mpz_xor(I[row], I[row], I[next_row - 1]);
			}
		}
	}
}

/* does not check for bounds, the caller must */
void get_matrix_row(mpz_t rop, MATRIX *matrix, uint64_t row_index) {
	mpz_set(rop, matrix->MATRIX[row_index]);
}

void get_identity_row(mpz_t rop, MATRIX *matrix, uint64_t row_index) {
	mpz_set(rop, matrix->IDENTITY[row_index]);
}

int test() {
	MATRIX matrix;

	init_matrix(&matrix, 6, 6);
	mpz_t r1, r2, r3;
	mpz_init(r1);
	mpz_init(r2);
	mpz_init(r3);

	mpz_set_ui(r1, 12);
	mpz_set_ui(r2, 11);
	mpz_set_ui(r3, 23);

	push_row(&matrix, r1);
	push_row(&matrix, r2);
	push_row(&matrix, r3);

	mpz_set_ui(r3, 1);
	push_row(&matrix, r3);

	mpz_set_ui(r3, 1);
	push_row(&matrix, r3);

	mpz_set_ui(r3, 23);
	push_row(&matrix, r3);

	print_matrix_matrix(&matrix);
	print_matrix_identity(&matrix);

	gauss_elimination(&matrix);

	print_matrix_matrix(&matrix);
	print_matrix_identity(&matrix);

	return 0;
}

//#define TEST_MATRIX
#ifdef  TEST_MATRIX
int main(){
#else
void __main(){
#endif
    test();
}