/*
 *  based on work of martani Jan 7, 2012
 *  enhanced for SIMPQS by andysnake96
 */
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <inttypes.h>
#include "matrix.h"
#include <math.h>
#include <unistd.h>
#include <signal.h>

int setIdentityMatrix(MATRIX* matrix){
    /*
     * set the identity matrix associated to exponent vector matrix,
     * each initial row r can be identified by 1 at r column and r row in the identity matrix
     * because of that identity matrix is wider the exponent vector matrix matrix
     */
    if (!(matrix->IDENTITY= (mpz_t *) calloc(matrix->rowsN, sizeof(*(matrix->IDENTITY))))){
        fprintf(stderr,"Out of mem at Identity matrix calloc\n");
        return EXIT_FAILURE;
    }
    for (u_int64_t j = 0; j <matrix->rowsN ; ++j) {
        mpz_init2(matrix->IDENTITY[j],matrix->rowsN);               //allocate rows for matrix identity
        mpz_setbit(matrix->IDENTITY[j],j);                      //set [j][j] in identity matrix
    }
    return EXIT_SUCCESS;
}

#define CLEAR_ARRAY_ENTRY(arrayEntry)\
    mpz_clears(arrayEntry->x,arrayEntry->element,arrayEntry->exp_vector,NULL);\
    if(arrayEntry->largePrime->_mp_d)\
        mpz_clear(arrayEntry->largePrime);

int initMatrixFromRows(MATRIX *matrix, uint64_t rowsN,struct ArrayEntry* rowsRaw, uint64_t cols) {
    //initialize exponent vector matrix from reports aggregated in rowsRaw
    //identity matrix will be setted after the standard one

    matrix->MATRIX = calloc(1, sizeof(*(matrix->MATRIX))*rowsN);
    if(!(matrix->MATRIX)){
        fprintf(stderr,"Out of meme at matrix calloc\n");
        return EXIT_FAILURE;
    }

    for (u_int64_t i = 0; i < rowsN; ++i) {
        mpz_init_set(matrix->MATRIX[i],rowsRaw[i].exp_vector);
    }
    matrix->rowsN = rowsN;
    matrix->colsN = cols;
    matrix->next_free_row = rowsN;
    printf("new matrix rows: %lu x cols:%lu\n",matrix->rowsN,cols);
    return setIdentityMatrix(matrix);
}

void print_matrix_matrix(MATRIX *matrix) {
	uint64_t i, j;

	printf("\nMATRIX\n");
	for (i = 0; i < matrix->rowsN; i++) {
		printf("%4lu- [",i);
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


void checkMatrixReducted(MATRIX* matrix){
//    print_matrix_matrix(matrix);
    //scrols col checking if only element i,i is 1
    int bit,result=EXIT_SUCCESS;
    for (u_int64_t i = 0; i < matrix->colsN; ++i) {
        for (u_int64_t j = 0; j < matrix->rowsN; ++j) {
            bit=mpz_tstbit(matrix->MATRIX[j],i);
            if((j>i && bit!=0)){
                result=EXIT_FAILURE;
                fprintf(stderr,"invalid bit %d at row;col: %lu\t%lu\n",bit,j,i);
            }
        }
    }
    for (u_int64_t i = 0,j=0; i < matrix->colsN && j<matrix->rowsN; ++i,++j) {
        bit=mpz_tstbit(matrix->MATRIX[j],i);
        while(bit==0 && i< matrix->colsN){
            bit=mpz_tstbit(matrix->MATRIX[j],++i);
        }//founded pivot -> all bellow has to 0
        for (u_int64_t k = j+1; k < matrix->rowsN; ++k) {
            if((bit=mpz_tstbit(matrix->MATRIX[k],i))!=0){
                fprintf(stderr,"not null bit under pivot(at %lu , %lu) at row %lu\n",j,i,k);
                result=EXIT_FAILURE;
            }
        }
    }
    if(result==EXIT_FAILURE)
        kill(getpid(),SIGFPE);
} //TODO DEBUG

/* performs a Gauss elimination on matrix->MATRIX, M:I---->M':I'
 * result (linear dependence) will be in the matrix->IDENTITY
 *   TODO BECAUSE EACH OP TO MAKE MATRIX IN ROW ECHELON FORM it has been applied on Identity Matrix too
 *      so M'[i]=I'[i]*M  **
 *      target so it's to find linear dependencies (0 vectors in M') and re build in terms of source terms **
 * */
void gauss_elimination_over_GF2(MATRIX *matrix) {
	printf("\nPerforming Gauss elimination..\n");
	mpz_t *m = matrix->MATRIX;
	mpz_t *I = matrix->IDENTITY;
	uint64_t col, row, next_row;
	int64_t next_pivot;
	for (next_row = 0, col = 0; col < MIN(matrix->colsN, matrix->rowsN); col++)
	//for all cols is expected to found the pivot at next_row
	{
		next_pivot = -1;                              //contains pivot row index or -1 if not founded in current column
		for (row = next_row; row < matrix->rowsN; row++) /* search for the next pivot in the column */
		{
			if (mpz_tstbit(m[row], col)) {
				next_pivot = row; /* row contains the next pivot */
				next_row++;
				break;
			}
		}

		if (next_pivot == -1)   //[0] col
        {
//            fprintf(stderr, "NOT FOUNDED PIVOT at col\t%lu\n",col);
            continue;
        }

		if (next_pivot != (next_row - 1)) /* current row is not the pivot, switch rows */
		{
//		    printf("Sw\n");fflush(0);
			mpz_swap(m[next_pivot], m[next_row - 1]);
			mpz_swap(I[next_pivot], I[next_row - 1]);
		}

		/* XOR the rows to eliminate all the 1 bellow the pivot*/
		for (row = next_row; row < matrix->rowsN; row++) {
			if (mpz_tstbit(m[row], col)) {
				mpz_xor(m[row], m[row], m[next_row - 1]);
				mpz_xor(I[row], I[row], I[next_row - 1]);
			}
		}
	}
    fflush(0);printf("Check if matrix is correctly reduced\n");
	checkMatrixReducted(matrix);
}

int quadraticRelationTry(REPORTS* reports,MATRIX* matrix){
    /* start backward from last row in the matrix  to find 0 vectors in final matrix --> linear depenecies --> square factor combination*/
    checkReports(reports,true);
    int linearRelFoundedN= 0,result=EXIT_SUCCESS;
    mpz_t matrixRow, identityMatrixRow,X, Y_element,tmp,tmp2;
    mpz_inits(matrixRow, identityMatrixRow, X, Y_element, tmp, tmp2, NULL);

    //// search for 0 vectors matrix backward from end
    uint64_t row_index = matrix->rowsN-1;
    mpz_set(matrixRow,matrix->MATRIX[row_index--]);
    while (mpz_cmp_ui(matrixRow, 0) == 0) {
        linearRelFoundedN++;
        mpz_set(matrixRow,matrix->MATRIX[row_index--]);
    }
    printf("\tLinear dependent relations found : %d\nTring GCDs\n", linearRelFoundedN);sleep(1);
    ////// Factor trys
    struct ArrayEntry arrayEntry;mpz_inits(arrayEntry.x,arrayEntry.largePrime,arrayEntry.element,arrayEntry.exp_vector,NULL);  //TODO DEBUG
    fflush(0);
    for (row_index+=2;row_index<matrix->rowsN;row_index++){
        mpz_set_ui(X, 1);
        mpz_set_ui(Y_element, 1);
        mpz_set_ui(tmp,0);                              //TODO CHECK EXP VECT SUM AT 0
        /*I'[i]*M=M'[i] <-- M'[i] is 0 vector <- founded linear dependency <-> BSmooth factors quadratic combination
         * build quadratic relation using founded linear dependency combining previusly founded X::f(x)=BSmoothVal::Y indexed on M
         */
        mpz_set(identityMatrixRow,matrix->IDENTITY[row_index]);
        mpz_init2(matrixRow,matrix->colsN);             //TODO DEBUG MEM CURRPUTION WITH GAUSS
        for (u_int64_t i = 0; i < matrix->rowsN; i++) {
            if (mpz_tstbit(identityMatrixRow, i)) {
                mpz_mul(X, X, reports->bsmoothEntries[i].x);
                mpz_mod(X, X, reports->n);
                mpz_mul(Y_element, Y_element,reports->bsmoothEntries[i].element);
                mpz_xor(matrixRow,matrixRow ,reports->bsmoothEntries[i].exp_vector);
            }
        }
        //TODO DEBUG CHECK
        if(mpz_cmp_ui(matrixRow,0))
            fprintf(stderr,"WTF GAUSS MISSED IN MEM BUG\n");
//        mpz_set(arrayEntry.x,X);mpz_set(arrayEntry.element,Y_element);mpz_init2(arrayEntry.exp_vector,matrix->colsN);
//        CHECK_X_SQURARE_CONGRUENT_Y_MOD_N(&arrayEntry, tmp, tmp2, true);
        fflush(0);
        mpz_sqrtrem(tmp2, tmp, Y_element);
        if(!mpz_cmp_ui(tmp,0)){
            gmp_printf("%lu)\tX=\t%Zd;\tY=%Zd;m()\n",row_index ,X,Y_element);fflush(0);
        }

        else{
            gmp_fprintf(stderr,"%lu)\tX=\t%Zd;\tY=%Zd;m()\n\n",row_index ,X,Y_element);
            for (int j = 0; j < matrix->colsN; j++) {
                fprintf(stderr,"%2d", mpz_tstbit(matrix->MATRIX[row_index], j));
            }
            fprintf(stderr," \n\n");
            for (int j = 0; j < matrix->colsN; j++) {
                fprintf(stderr,"%2d", mpz_tstbit(identityMatrixRow, j));
            }
            fprintf(stderr," \n\n");

            for (int j = 0; j < matrix->colsN; j++) {
                fprintf(stderr,"%2d", mpz_tstbit(matrixRow, j));
            }
            fprintf(stderr," \n\n");

            exit(EXIT_FAILURE);
        }
    
        fflush(0);
        mpz_mod(Y_element, tmp2, reports->n);
        mpz_sub(tmp, X, Y_element);
        mpz_add(tmp2, X, Y_element);
        mpz_gcd(tmp, tmp,reports->n );
        mpz_gcd(tmp2, tmp2,reports->n );
        if ((mpz_cmp(tmp,reports->n ) != 0 && mpz_cmp_ui(tmp, 1) != 0) || (mpz_cmp(tmp2,reports->n ) != 0 && mpz_cmp_ui(tmp2, 1) != 0)){
            mpz_cdiv_q(Y_element,reports->n , tmp);                        //set the other factor
            gmp_printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\tFOUNDED FACTORS: %Zd \t :%Zd\t!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",tmp,tmp2);
            gmp_printf("with x:%Zd \t\t y:%Zd \n",X,Y_element);
            result=EXIT_SUCCESS;
            goto exit;
        }
    }
    printf("useless trys: %d\n",linearRelFoundedN);
    result=EXIT_FAILURE;
    exit:
    mpz_clears(X, Y_element, matrixRow, identityMatrixRow, NULL);
    return result;
}
