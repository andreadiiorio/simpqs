#include <stdio.h>
#include "sievingSIMPQS.h"
#include <unistd.h>
#include <utils/utils.h>
#include <signal.h>

#define ENTRY_MAX_LEN 4096;

#define WRITE_SEPARATOR(fd,sep) \
    fputc(sep,fd);

#define WRITE_WORD_SEPARATOR_raw(fd) \
if(fwrite(&WORD_SEPARATOR, sizeof(WORD_SEPARATOR),1,fd)!=1){\
        fprintf(stderr,"separator write error\n");\
        result=EXIT_FAILURE;goto exit;}

const char WORD_SEPARATOR='\t';
const char LINE_SEPARATOR='\n';
unsigned int dbg_h,dbg_exp=0; u_int64_t tmpFactor;mpz_t tmpL;           //TODO DEBUG
//DYNAMIC_VECTOR* primes_B;
DYNAMIC_VECTOR FB;
CONFIGURATION* Configuration;
//#define NULL_FILTER_RACE_COND //TODO ROBSTNESS DEBUG MACRO


/*
 * DEBUG CHECKS ON REPORT ELEMENTS:
 *  - CHECKED IF X^2==a_fx MOD N
 *  - CHECKED IF BIT_i=1 IN EXP VECTOR OF a_fx CORRESPOND TO A p_i THAT DIVDE a_fx AN ODD NUMBER OF TIMES (notting dbg_exp)
 */
void CHECK_X_SQURARE_CONGRUENT_Y_MOD_N(struct ArrayEntry* arrayEntry, mpz_t tmp, mpz_t tmp2, bool largePrimeAggregatedEntryCheck) {
    dbg_h = EXIT_SUCCESS;
    mpz_pow_ui(tmp, arrayEntry->x, 2);
    mpz_mod(tmp, tmp, Configuration->N);
    mpz_mod(tmp2,arrayEntry->element, Configuration->N);
    if (mpz_cmp(tmp, tmp2)) {
        gmp_fprintf(stderr, "Mismatch X(j)^2:\t%Zd !== af(j):\t %Zd\n ", tmp, tmp2);
        dbg_h = EXIT_FAILURE;
    }
    if (mpz_tstbit(arrayEntry->exp_vector, 0) && mpz_cmp_ui(arrayEntry->element, 0) >= 0) {    //test sign bit
        fprintf(stderr, "missed sign bit\n");
    }
    for (unsigned long i = 1,bit=0; i < FB.vectorSize; ++i) {              //test odd power bits
        bit=mpz_tstbit(arrayEntry->exp_vector, i);
        dbg_exp = 0;
        mpz_set_ui(tmp,((u_int64_t *) FB.pntr)[i - 1]);
        while (mpz_divisible_p(arrayEntry->element, tmp)) {
            dbg_exp++;
            mpz_mul_ui(tmp,tmp,((u_int64_t *) FB.pntr)[i - 1]);
        }
        if ((dbg_exp%2)!=bit) {
            fprintf(stderr, "bit missetted at i:%d exp->%d \n", i,dbg_exp);
            dbg_h = EXIT_FAILURE;
        }
    }
    if(largePrimeAggregatedEntryCheck && arrayEntry->largePrime->_mp_size) {
        mpz_set(tmp,arrayEntry->largePrime);
        dbg_exp=0;
        while(mpz_divisible_p(arrayEntry->element,tmp)) {
            dbg_exp++;
            mpz_mul(tmp,tmp,arrayEntry->largePrime);
        }
        if(dbg_exp%2!=0) {
            gmp_fprintf(stderr, "not even exponent with aggregated large prime:\t%Zd\t^ %d\n",arrayEntry->largePrime, dbg_exp);
//            dbg_h = EXIT_FAILURE;
        }
    }
    if (dbg_h == EXIT_FAILURE) {
        gmp_fprintf(stderr,"\nY=%Zd\n", arrayEntry->element);
        kill(0,SIGFPE);
    }
}
int checkReports(REPORTS *reports, bool aggregatedLargePrimeCheck) {
    int result=EXIT_SUCCESS;
    mpz_t tmp,tmp2;mpz_inits(tmp,tmp2,NULL);
    struct ArrayEntry* baseArrayToCheck=reports->bsmoothEntries;
    int sizeArray=reports->relationsNum;
    //// check BSMooth Entries
    for (int j = 0; j < sizeArray; ++j) {
        fflush(0);
        CHECK_X_SQURARE_CONGRUENT_Y_MOD_N(baseArrayToCheck+j, tmp, tmp2,aggregatedLargePrimeCheck);
    }
    mpz_clears(tmp,tmp2,NULL);
    return result;
}

int mergeReports(REPORTS *dstReports, const REPORTS *new_reports) {
    //merge founded reports  into dstReports
    /// realloc array  entries of enough space to hold newly founded (partial) reports
    u_int64_t newPartialReportsN = dstReports->partialRelationsNum, newReportsN = dstReports->relationsNum;
    newPartialReportsN = dstReports->partialRelationsNum + new_reports->partialRelationsNum;
    newReportsN = dstReports->relationsNum + new_reports->relationsNum;
    if (newReportsN && !(dstReports->bsmoothEntries = realloc(dstReports->bsmoothEntries, newReportsN *
                                                                                          sizeof(*(new_reports->bsmoothEntries))))) {
        return EXIT_FAILURE;
    }
    if (newPartialReportsN && !(dstReports->largePrimesEntries = realloc(dstReports->largePrimesEntries,
                                                                         newPartialReportsN *
                                                                         sizeof(*(new_reports->largePrimesEntries))))) {
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
    dstReports->partialRelationsNum = newPartialReportsN;
    dstReports->relationsNum = newReportsN;
#ifdef DEBUG_CHECK
    return checkReports(dstReports, false);
#endif
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
static int compareLargePrimeReports(const void* p1,const void* p2){
    struct ArrayEntry* largePrimeEntry1;
    struct ArrayEntry* largePrimeEntry2;
    largePrimeEntry1=((struct ArrayEntry*) p1);
    largePrimeEntry2=((struct ArrayEntry*) p2);
    return mpz_cmp(largePrimeEntry1->largePrime,largePrimeEntry2->largePrime);          //return int according to comparison of large primes
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
    memset(&(destEntry->x),0, sizeof(mpz_t)); //SEE TODO
#endif
    mpz_init_set(destEntry->largePrime,entry->largePrime);
    mpz_init_set(destEntry->exp_vector,entry->exp_vector);
    mpz_init_set(destEntry->element,entry->element);    //SEE TODO
    mpz_init_set(destEntry->x,entry->x);    //SEE TODO
}
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


/*
 * REPORTS_FILE_ORGANIZATION
 * N\ta\tb\tNumReports\tNumPartialReports\n
 * j\tx\tBSmoothEntryElement\texp_vector\n
 * .......
 * \n\n
 * LargePrime\tj\tx\tBSmoothEntryElement\texp_vector\n
 * ....
 */
int saveReports(REPORTS *reports, u_int64_t colsN, bool printReports,bool polynomialFamily ,struct polynomial *polynomial) {
    // serialize founded reports to raw file using gmp mpz serialization function
    // reports file structure is defined above under name REPORTS_FILE_ORGANIZATION
    //if polynomial family flag is set b will skipped in filename indicating serialized reports are for the whole polynomial family
    //if printReports is true reports will be printed in stdout without fflush per line
    const int formattedOutputBufLen=colsN+2+2*ENTRY_MAX_LEN;
    char* outBuf=calloc(1,sizeof(*outBuf)*formattedOutputBufLen);
    if(!outBuf){
        fprintf(stderr,"Out of mem in report saving\n");
    }

    //////////// writing file name
    int writtenInBuf,writtenInFile,result=EXIT_SUCCESS;
    //repot filename identified by N & polynomial
    if(!polynomialFamily) {
        if ((writtenInBuf = gmp_snprintf(outBuf, formattedOutputBufLen, "%s%Zd_%Zd_%Zd_%lu_%lu%s",
                                        REPORTS_FILENAME_PREFIX, *(polynomial->N), polynomial->a.a, polynomial->b,
                                         reports->relationsNum, reports->partialRelationsNum,
                                         REPORTS_FILENAME_SUFFIX)) >= formattedOutputBufLen) {
            fprintf(stderr, "snprintf hasn't founded enough space to output to buf\n");
            free(outBuf);
            return EXIT_FAILURE;
        }
    }else { // if polynomialFamily is true reports are the aggregation of all polynomial family reports-> skip b in name to distinguish
        if ((writtenInBuf = gmp_snprintf(outBuf, formattedOutputBufLen, "%s%Zd_%Zd_%lu_%lu%s",
                                         REPORTS_FILENAME_PREFIX, *(polynomial->N), polynomial->a.a,
                                         reports->relationsNum, reports->partialRelationsNum,
                                         REPORTS_POLYNOMIAL_FAMILY_FILENAME_SUFFIX)) >= formattedOutputBufLen) {
            fprintf(stderr, "snprintf hasn't founded enough space to output to buf\n");
            free(outBuf);
            return EXIT_FAILURE;
        }

    }
//    fprintf(stderr,"serializing %s\n",outBuf);
    FILE*reportFp=fopen(outBuf, "w");
    
    if(!reportFp){
        perror("failed to create reports file");
        result=EXIT_FAILURE;goto exit;
    }

    //////////// writing file header
    if(!(writtenInFile=mpz_out_raw(reportFp, *(polynomial->N)))){                                      //N
        fprintf(stderr,"header write error\n");
        result=EXIT_FAILURE;goto exit;
    }
    WRITE_SEPARATOR(reportFp,'\t')
    if(!(writtenInFile=mpz_out_raw(reportFp, (polynomial->a.a)))){                                       //a
        fprintf(stderr,"header write error \t a \n");
        result=EXIT_FAILURE;goto exit;
    }
    WRITE_SEPARATOR(reportFp,'\t')
    if(!(writtenInFile=mpz_out_raw(reportFp, (polynomial->b)))){                                       //b
        fprintf(stderr,"header write error\t b\n");
        result=EXIT_FAILURE;goto exit;
    }
    WRITE_SEPARATOR(reportFp,'\t')
    if(fwrite(&reports->relationsNum, sizeof(reports->relationsNum), 1, reportFp) != 1){              //report #
        fprintf(stderr,"separator write error\n");
        result=EXIT_FAILURE;goto exit;
    }
    WRITE_SEPARATOR(reportFp,'\t')
    if(fwrite(&reports->partialRelationsNum, sizeof(reports->partialRelationsNum), 1, reportFp) != 1){ //partial reports #
        fprintf(stderr,"separator write error\n");
        result=EXIT_FAILURE;goto exit;
    }
    WRITE_SEPARATOR(reportFp,'\n')


    //////////// standard reports
    mpz_t* row;
    for (u_int64_t i= 0; i < reports->relationsNum ; i++) {
        printReports ? gmp_printf("%d> %30Zd  ", i, reports->bsmoothEntries[i].element) : 0;
        row=&(reports->bsmoothEntries[i].exp_vector);
        for (u_int64_t j = 0; j < colsN && printReports; j++) {
            printf("%d", mpz_tstbit(*row,j));
        }
        printReports ? printf(" \n") : 0;
        //row raw file write
        if(fwrite(&(reports->bsmoothEntries[i].j), sizeof(reports->bsmoothEntries->j), 1, reportFp) != 1){   //j
            fprintf(stderr,"j write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\t')

        if(!(writtenInFile=mpz_out_raw(reportFp, (reports->bsmoothEntries[i].x)))){                          //x
            fprintf(stderr,"x write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\t')
        if(!(writtenInFile=mpz_out_raw(reportFp, (reports->bsmoothEntries[i].element)))){                  //element
            fprintf(stderr,"bsmooth element write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\t')
        if(!(writtenInFile=mpz_out_raw(reportFp, (*row)))){                                                //row
            fprintf(stderr," row write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\n')
    }
    //WRITE REPORTS SEPARATOR
    WRITE_SEPARATOR(reportFp,'\n')
    WRITE_SEPARATOR(reportFp,'\n')

    ////////// partial reports
    for (u_int64_t k = 0; k < reports->partialRelationsNum ; ++k) {
        row=&(reports->largePrimesEntries[k].exp_vector);
        printReports ?gmp_printf("%d> %30Zd\t",k,reports->largePrimesEntries[k].largePrime):0;
        for (u_int64_t j = 0; j < colsN && printReports; j++) {
            printf("%d", mpz_tstbit(*row,j));
        }
        printReports ?printf("\n"):0;
        //row raw file write
        if(!(writtenInFile=mpz_out_raw(reportFp, (reports->largePrimesEntries[k].largePrime)))){          //large prime
            fprintf(stderr,"large prime  write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\t')
        if(fwrite(&(reports->largePrimesEntries[k].j), sizeof(reports->bsmoothEntries->j), 1, reportFp) != 1){   //j
            fprintf(stderr,"j write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\t')

        if(!(writtenInFile=mpz_out_raw(reportFp, (reports->largePrimesEntries[k].x)))){                          //x
            fprintf(stderr,"x write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\t')
        if(!(writtenInFile=mpz_out_raw(reportFp, (reports->largePrimesEntries[k].element)))){              //element
            fprintf(stderr,"large prime element write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\t')
        if(!(writtenInFile=mpz_out_raw(reportFp, (*row)))){                                                //row
            fprintf(stderr," row write error\n");
            result=EXIT_FAILURE;goto exit;
        }
        WRITE_SEPARATOR(reportFp,'\n')
    }

    exit:
    free(outBuf);
    fclose(reportFp);
    return result;
}

void _deleteLocalReports(bool deletePolynomialReports) {
    char *cmd;
    cmd = "find -iname \"reports_*\" | xargs -d \"\\n\" rm";
    if(deletePolynomialReports)
        cmd = "find -iname \"reports_*reportslist\" | xargs -d \"\\n\" rm";
    FILE *fp = popen(cmd, "r");
    if (fp == NULL) {
        printf("Failed to run command find\n");
        exit(EXIT_FAILURE);
    }
    pclose(fp);
}

char** findReportsLocally(unsigned int numReports,const char* reportSuffix) {
    // find all reports file saved locally by siever processes;
    // returned an array of strings containg paths of reports file or NULL if an error occurred
    // if less then numReports has been found NULL will be setted first empty entry of the array and error message printed

    const char* FIND_REPORTS_BASH_CMD_LINUX="find \"$(pwd -P)\" -iname \"reports_*";
    const size_t MAX_FIND_CMD_SIZE=(strlen(FIND_REPORTS_BASH_CMD_LINUX)+5+strlen(REPORTS_POLYNOMIAL_FAMILY_FILENAME_SUFFIX));
    char* findCmdBuf=malloc(sizeof(*findCmdBuf)*MAX_FIND_CMD_SIZE);
    if(!findCmdBuf) {
        fprintf(stderr, "Out of mem at find cmd buf allocate :(\n");
        return NULL;
    }
    snprintf(findCmdBuf,MAX_FIND_CMD_SIZE,"%s\"%s",FIND_REPORTS_BASH_CMD_LINUX,reportSuffix);
    printf("searching reports file with cmd str:\t%s\n",findCmdBuf);
    const int MAX_PATH_SIZE = 2048;
//    const char* FIND_REPORTS_BASH_CMD_LINUX="find -iname 'reports_*'";
    char **paths = calloc(numReports,  sizeof(*paths));
    if (!paths) {
        fprintf(stderr,"Out of Mem at path** calloc\n");
        return NULL;
    }
    /* Open the command  find for reading. */
    FILE *findPipe = popen(findCmdBuf, "r");
    if (findPipe == NULL) {
        printf("Failed to run command find\n");
        return NULL;
    }
    int commandOutTerminated=0;
    for (unsigned int i = 0; i < numReports && !commandOutTerminated; ++i) {
        //for each expected output file get file path  by executed command find line by line
        char **destPathBuf = &(paths[i]);
        if(!(*destPathBuf = malloc(sizeof(*(*destPathBuf)) * MAX_PATH_SIZE))){
            fprintf(stderr,"Out of mem at path allocation :%d\n",i);
            goto exit_failure;
        }
        commandOutTerminated = (fgets(*destPathBuf, MAX_PATH_SIZE - 1, findPipe)) == NULL;
        //// remove trailing newline
        for (int j=0;j<MAX_PATH_SIZE;j++){
            if ((*destPathBuf)[j] == '\n') {
                (*destPathBuf)[j] = 0;
                break;
            }
        }
        if (commandOutTerminated){
            free(*destPathBuf);
            *destPathBuf=NULL;
            printf("terminated before expected\n");
        }

        printf("founded reports file:%s \n",*destPathBuf);
    }
    pclose(findPipe);
    return paths;
    exit_failure:
    if(paths){
        for (unsigned int i = 0; i < numReports; ++i) {
            free(paths[i]);
        }
    }
    free(paths);free(findCmdBuf);
    return NULL;
}
REPORTS *loadReports(char *filePath) {
    FILE* reportsFp=fopen(filePath, "r");
    if(!reportsFp){
        perror(filePath);
        return NULL;
    }
    mpz_t N;mpz_init(N);
    struct polynomial polynomial=(struct polynomial){.N=&N};
    int nread,result=EXIT_SUCCESS;unsigned char sep;
    REPORTS* reports=calloc(1,sizeof(*reports));
    if(!reports){
        fprintf(stderr,"Out of MEM at reports malloc\n");
        result=EXIT_FAILURE;goto exit;
    }
    ////////// reading file header
    if(!(nread =mpz_inp_raw(N,reportsFp))){                                 //N
        fprintf(stderr,"invalid read of N in header \n");
        result=EXIT_FAILURE;goto exit;
    }
    mpz_init_set(reports->n,N);                                         //set n in reports for convenience on later restart
    sep=fgetc(reportsFp);
    if(!(nread =mpz_inp_raw(polynomial.a.a,reportsFp))){                      //a
        fprintf(stderr,"invalid read of N in header \n");
        result=EXIT_FAILURE;goto exit;
    }
    sep=fgetc(reportsFp);
    if(!(nread =mpz_inp_raw(polynomial.b,reportsFp))){                      //b
        fprintf(stderr,"invalid read of N in header \n");
        result=EXIT_FAILURE;goto exit;
    }
    sep=fgetc(reportsFp);

    u_int64_t reportsN,partialReportsN;
    if(fread(&reportsN, sizeof(reportsN), 1, reportsFp) != 1){          //Num reports
        fprintf(stderr,"invalid read of number of reports\n");
        result=EXIT_FAILURE;goto exit;
    }
    sep=fgetc(reportsFp);

    if(fread(&partialReportsN, sizeof(partialReportsN),1,reportsFp)!=1){    //Num partial reports
        fprintf(stderr,"invalid read of number of partial reports\n");
        result=EXIT_FAILURE;goto exit;
    }
    sep=fgetc(reportsFp);               //newline
    gmp_printf("de serializing N:%Zd\ta: %Zd\tb:%Zd\t reportsN:%lu\tpartialReportsN:%lu\n",N,polynomial.a.a,polynomial.b,reportsN,partialReportsN);fflush(0);

    ///// allocate reports Num
    if(!(reports->bsmoothEntries=malloc(sizeof(*(reports->bsmoothEntries))*reportsN))){
        fprintf(stderr,"Out of mem at reports allocation\n");
        result=EXIT_FAILURE;goto exit;
    }
    if(!(reports->largePrimesEntries=malloc(sizeof(*(reports->largePrimesEntries))*partialReportsN))){
        fprintf(stderr,"Out of mem at partial reports allocation\n");
        result=EXIT_FAILURE;goto exit;
    }
    reports->relationsNum=reportsN;
    reports->partialRelationsNum=partialReportsN;

    ///////////// standard reports
    for (u_int64_t i = 0; i < reportsN; ++i) {
        mpz_inits(reports->bsmoothEntries[i].element,reports->bsmoothEntries[i].exp_vector,reports->bsmoothEntries[i].x,NULL);
        //row raw file  read
        if(fread(&(reports->bsmoothEntries[i].j), sizeof(reports->bsmoothEntries->j),1,reportsFp)!=1){    //j
            fprintf(stderr,"invalid read of number of partial reports\n");
            result=EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);               //newline
        if (!(nread= mpz_inp_raw(reports->bsmoothEntries[i].x,reportsFp))) {                                //x
            fprintf(stderr, "x element read error\n");
            result = EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);
        if (!(nread= mpz_inp_raw(reports->bsmoothEntries[i].element,reportsFp))) {                          //element
            fprintf(stderr, "bsmooth element read error\n");
            result = EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);

        if (!(nread= mpz_inp_raw(reports->bsmoothEntries[i].exp_vector,reportsFp))) {                       //row
            fprintf(stderr, "bsmooth element read error\n");
            result = EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);                    //newline
        memset(reports->bsmoothEntries[i].largePrime,0, sizeof(*(reports->largePrimesEntries[i].largePrime) )); //mark as UNSETTED(0) large prime inside standard entry

    }
    ///// remove separator of reports
    sep=fgetc(reportsFp);                      //newline
    sep=fgetc(reportsFp);                   //newline

    ///// get Large prime entries
    for (u_int64_t i = 0; i < partialReportsN; ++i) {
        mpz_inits(reports->largePrimesEntries[i].largePrime,reports->largePrimesEntries[i].element,reports->largePrimesEntries[i].x,reports->largePrimesEntries[i].exp_vector,NULL);
        //row raw file  read
        if (!(nread= mpz_inp_raw(reports->largePrimesEntries[i].largePrime,reportsFp))) {               //large prime
            fprintf(stderr, "large prime read err\n");
            result = EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);
        if(fread(&(reports->largePrimesEntries[i].j), sizeof(reports->bsmoothEntries->j),1,reportsFp)!=1){    //j
            fprintf(stderr,"j read err \n");
            result=EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);               //
        if (!(nread= mpz_inp_raw(reports->largePrimesEntries[i].x,reportsFp))) {                            //x
            fprintf(stderr, "large prime x read error\n");
            result = EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);
        if (!(nread= mpz_inp_raw(reports->largePrimesEntries[i].element,reportsFp))) {                          //element
            fprintf(stderr, "large prime bsmooth element read error\n");
            result = EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);
        if (!(nread= mpz_inp_raw(reports->largePrimesEntries[i].exp_vector,reportsFp))) {               //row
            fprintf(stderr, "large prime row element read error\n");
            result = EXIT_FAILURE;goto exit;
        }
        sep=fgetc(reportsFp);                   //newline

    }

    exit:
        fclose(reportsFp);
        if(result==EXIT_FAILURE) {
            free(reports->largePrimesEntries);
            free(reports->bsmoothEntries);
            free(reports);
            return NULL;
        }
#ifdef DEBUG_CHECK
    checkReports(reports, false);
#endif
        return reports;
}

void sortPartialReports(REPORTS* reports){
    qsort(reports->largePrimesEntries,reports->partialRelationsNum, sizeof(*(reports->largePrimesEntries)),compareLargePrimeReports);
//    for (u_int64_t i = 0; i < reports->partialRelationsNum; ++i) {
//        gmp_printf("%4lu\t-- %9Zd\t\t",i,reports->largePrimesEntries[i].largePrime);
//        if(i%3==0)
//            printf("\n");
//    }
//    fflush(0);
}
int  pairLargePrimeEntries(REPORTS *report, unsigned int startAddr, unsigned int matchFounded) {
    /*
     * Pair large prime partial reports from startAddr to startAddr+matchFounded
     * will be coupled first partial report in the series with the others
     * TODO CHECK ASK?
     * large prime will be divided out from reports element ( P^2==0 in the exp vector mod 2)
     * x and resulting element  will be coupled with a multiplication
     * exp vectors will be ORed
     *      because if e1=x^2*y (0,1) and e2=x*y^2 (1,0) => e1*e2=x^3*y^3 (1,1 in the exponent vector)
     */
    unsigned int newSizeReports=report->relationsNum + matchFounded,newReportIndx=report->relationsNum;
    REALLOC_WRAP((newSizeReports), report->relationsNum, report->bsmoothEntries, (matchFounded))
            return EXIT_FAILURE;
        }}
    struct ArrayEntry* partialReport0=report->largePrimesEntries+startAddr,*partialReportk,*newReport=report->bsmoothEntries+newReportIndx;
    for (unsigned int k = 1; k <=matchFounded; ++k) {
        //// couple partial report 0 with partial report k>0
        partialReportk=report->largePrimesEntries+startAddr+k;
        arrayEntryCopy(newReport,partialReport0);
        //TODO !!!!!!!!!!!!!!!!!!!! remove large prime from elements because  P^2 =0  in exp vector
//        mpz_divexact(partialReportk->element,partialReportk->element,partialReportk->largePrime);
//        mpz_divexact(newReport->element,newReport->element,newReport->largePrime);
        // couple partial reports multipling x and element
        mpz_mul(newReport->x,newReport->x,partialReportk->x);
        mpz_mul(newReport->element,newReport->element,partialReportk->element);
        mpz_xor(newReport->exp_vector,newReport->exp_vector,partialReportk->exp_vector);
#ifdef   DEBUG_CHECK
        gmp_printf("coupled large prime entries in \t X=%Zd;\tY=%Zd\n",newReport->x,newReport->element);fflush(0);
        if (mpz_cmp(partialReport0->largePrime,partialReportk->largePrime)) {
            gmp_fprintf(stderr, "Not matching large primes %Zd\t%Zd\n", partialReport0->largePrime,partialReportk->largePrime);
            exit(EXIT_FAILURE);
        }
        mpz_init_set(tmpL,partialReport0->largePrime);
        dbg_exp=0;
        while (mpz_divisible_p(newReport->element,tmpL)) {
            mpz_mul(tmpL, tmpL,partialReport0->largePrime);
            dbg_exp++;
        }
        if(dbg_exp%2 || dbg_exp!=2) {
            fprintf(stderr, "Invalid Large prime exponent :%d\n", dbg_exp);
            exit(EXIT_FAILURE);
        }
#endif
        newReport++;
    }
    return EXIT_SUCCESS;
}
int pairPartialReports(REPORTS* reports) {
    //pair partialReports (that has been sorted by large prime
    sortPartialReports(reports);
    int matchFoundedAll = 0, i = 0;
    unsigned int matchStartAddr = 0;
    struct ArrayEntry *prevPartialReport = reports->largePrimesEntries;
    struct ArrayEntry *partialReport = (reports->largePrimesEntries + 1);
    for (unsigned int j = 1; j < reports->partialRelationsNum; partialReport = &(reports->largePrimesEntries[++j])) {
        if (!mpz_cmp(partialReport->largePrime, prevPartialReport->largePrime)) {  //founded a match
            if((i++)==0)
                matchStartAddr=j-1;
        }
        else if (i > 0) {                                                         //end of match series
            if(pairLargePrimeEntries(reports, matchStartAddr,  i)==EXIT_FAILURE)
                return EXIT_FAILURE;
            i = 0;
            matchStartAddr = j+1;
            matchFoundedAll++;
        }
        prevPartialReport = partialReport;
    }
    printf("founded %d\tlarge prime matches series\n",matchFoundedAll);
#ifdef DEBUG_CHECK
    checkReports(reports,true);
#endif
    return EXIT_SUCCESS;
}