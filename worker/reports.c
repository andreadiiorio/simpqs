#include "sievingSIMPQS.h"
#include <stdio.h>
#include <unistd.h>

void print_reports(REPORTS *reports, u_int64_t colsN, bool printMatrix) {
    u_int64_t i, j;
    //////////// standard reports
    printf("num primes in matrix rows: %lu\n",colsN);
    printf("\n printing: %lu reports\n",reports->relationsNum);
    mpz_t* row;
    for (i = 0; i < reports->relationsNum && printMatrix; i++) {
        gmp_printf("%d> %30Zd  ",i,reports->bsmoothEntries[i].element);
        row=&(reports->bsmoothEntries[i].exp_vector);
        for (j = 0; j < colsN; j++) {
            printf("%d", mpz_tstbit(*row,j));
        }
        printf(" \n");
    }
    printf("\n");

    ////////// partial reports
    printf("\n printing: %lu partial reports\n",reports->partialRelationsNum);
    for (u_int k = 0; k < reports->partialRelationsNum && printMatrix; ++k) {
        row=&(reports->largePrimesEntries[k].exp_vector);
        gmp_printf("%d> %30Zd\t",k,reports->largePrimesEntries[k].largePrime);
        for (j = 0; j < colsN; j++) {
            printf("%d", mpz_tstbit(*row,j));
        }
        printf(" \n");
    }
}

//#define NULL_FILTER_RACE_COND
int mergeReports(REPORTS *dstReports, const REPORTS *new_reports) {
    //merge founded reports  into dstReports
    /// realloc array  entries of enough space to hold newly founded (partial) reports
    u_int64_t newPartialReportsN= dstReports->partialRelationsNum + new_reports->partialRelationsNum;
    u_int64_t newReportsN= dstReports->relationsNum + new_reports->relationsNum;
    if (newReportsN && !(dstReports->bsmoothEntries=realloc(dstReports->bsmoothEntries,newReportsN* sizeof(*(new_reports->bsmoothEntries))))) {
        return EXIT_FAILURE;
    }
    if (newPartialReportsN && !(dstReports->largePrimesEntries=realloc(dstReports->largePrimesEntries,newPartialReportsN* sizeof(*(new_reports->largePrimesEntries))))) {
        return EXIT_FAILURE;
    }
    //// shallow copy of founded (partial) reports in dstReports
    for (u_int i = 0; i < new_reports->relationsNum; ++i) {
#ifdef NULL_FILTER_RACE_COND
        if (!(mpz_cmp_ui(new_reports->bsmoothEntries[i].element, 0))) {
            struct ArrayEntry* report=&new_reports->bsmoothEntries[i];  //TODO DEBUG
            fprintf(stderr, "skipping report at :%p\n", report);
            continue;
        }
#endif
        dstReports->bsmoothEntries[dstReports->relationsNum + i] = new_reports->bsmoothEntries[i];
    }
    for (u_int i = 0; i < new_reports->partialRelationsNum; ++i) {
#ifdef NULL_FILTER_RACE_COND
        if (!(mpz_cmp_ui(new_reports->largePrimesEntries[i].largePrime, 0))) {
            fprintf(stderr, "skipping partial report at :%d\n", i);
            continue;
        }
#endif
        dstReports->largePrimesEntries[dstReports->partialRelationsNum + i] = new_reports->largePrimesEntries[i];
    }
    /// update (partial) reports sizes
    dstReports->partialRelationsNum=newPartialReportsN;
    dstReports->relationsNum=newReportsN;
    return EXIT_SUCCESS;
}
int mergeReportsFast(REPORTS *dstReports, const REPORTS *new_reports) { //TODO DEBUG MEMMOVE CORRUPT EVERYTHING
    //merge founded reports  into dstReports
    /// realloc array  entries of enough space to hold newly founded (partial) reports
    u_int64_t newPartialReportsN= dstReports->partialRelationsNum + new_reports->partialRelationsNum;
    u_int64_t newReportsN= dstReports->relationsNum + new_reports->relationsNum;
    if (newReportsN && !(dstReports->bsmoothEntries=realloc(dstReports->bsmoothEntries,++newReportsN* sizeof(*(new_reports->bsmoothEntries))))) {
        return EXIT_FAILURE;
    }
    if (newPartialReportsN && !(dstReports->largePrimesEntries=realloc(dstReports->largePrimesEntries,++newPartialReportsN* sizeof(*(new_reports->largePrimesEntries))))) {
        return EXIT_FAILURE;
    }
    /// hard copy of every report, because of reports are already a copy of source array entries block, there is no need to recopy, just not free until computation is over
    memmove(&(dstReports->bsmoothEntries)+dstReports->relationsNum, &(new_reports->bsmoothEntries), (new_reports->relationsNum )* sizeof(*(new_reports->bsmoothEntries)));
    memmove(&(dstReports->largePrimesEntries)+dstReports->partialRelationsNum, &(new_reports->largePrimesEntries), (new_reports->partialRelationsNum ) * sizeof(*(new_reports->largePrimesEntries)));
    /// update (partial) reports sizes
    dstReports->partialRelationsNum=newPartialReportsN;
    dstReports->relationsNum=newReportsN;
    return EXIT_SUCCESS;
}


#define RESET_MPZ_ON_COPY           //0 RESET before overwrite re initializing mpz fields in the destination of the copy
void arrayEntryCopy(struct ArrayEntry *destEntry, struct ArrayEntry *entry) {
    //perform values copy of entry on destEntry with shallow copy of entire entry and mpz copy by mpz_init_set
    //TODO .elemnt COPY IS UNNECESSARY because of j inside struct already copied

    *destEntry=*entry;                  //SHALLOW COPY OF STRUCT ENTRY
#ifdef RESET_MPZ_ON_COPY            //reset MPZ_T fields because of gmp doc not multiple initializing .... bad lib:(
    memset(&(destEntry->exp_vector),0, sizeof(mpz_t));
    memset(&(destEntry->largePrime),0, sizeof(mpz_t));
    memset(&(destEntry->element),0, sizeof(mpz_t)); //SEE TODO
#endif
    if (&(destEntry->largePrime)!=NULL) {
        mpz_init_set(destEntry->largePrime,entry->largePrime);
    }
    mpz_init_set(destEntry->exp_vector,entry->exp_vector);
    mpz_init_set(destEntry->element,entry->element);    //SEE TODO
}

