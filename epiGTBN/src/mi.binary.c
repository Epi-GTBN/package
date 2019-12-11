#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include "include/matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAX 10000

int data[MAX][3*MAX]; // store binary data
int count; // should be the same value as nsamples
int CaseCount; // total case number (not total control number)
int SNPnumber; // should be nvars-1

// use binary AND to count appearance
// SNP >= 1, factor = 0/1/2
int CompareSNP(int SNP1, int SNP2, int factor){
    // Rprintf("CompareSNP called...\n");
    // Rprintf("in CompareSNP count = %d, CaseCount = %d\n", count, CaseCount);
    int result[3 * MAX] = {0};
    int i = 0, k = 0, common = 0;
    // Rprintf("factor * count = %d\n", factor * count);
    // Rprintf("(factor + 1) * count = %d\n", (factor + 1) * count);
    for(i = factor * count; i < (factor + 1) * count; i++){
        result[k] = data[SNP1 - 1][i] & data[SNP2 - 1][i];
        if(result[k]) common++;
        k++;
    }
    
    // Rprintf("the k is %d\n", k);
    
    // for(i = 0; i < k; i++){
    //     Rprintf("%d", result[i]);
    // }
    // Rprintf("\n");
    
    // Rprintf("in CompareSNP return (the common number) is %d\n", common);
    return common;
}


// use binary AND to count appearance
// SNP >= 1, factor = 0/1/2, type = 1(case)/0(control)
int comparePartSet(int SNP1, int factor1, int SNP2, int factor2, int type){
    // Rprintf("comparePartSet called...\n");
    // Rprintf("in comparePartSet count = %d, CaseCount = %d\n", count, CaseCount);
    int i, sum, common = 0;
    // case
    if(type){
        i = 0;
        sum = CaseCount;
    }else{ // control
        i = CaseCount;
        sum = count;
    }

    // Rprintf("factor1 * count = %d\n", factor1 * count);
    // Rprintf("factor2 * count = %d\n", factor2 * count);
    for(; i < sum; i++){
        int flag = data[SNP1 - 1][(factor1 * count) + i] & data[SNP2 - 1][(factor2 * count) + i];
        if(flag) common++;
    }

    // Rprintf("in comparePartSet return (the common number) is %d\n", common);
    return common;
}

double H_one() {
    // Rprintf("H_one called...\n");
    return - 0.5 * log(0.5) - 0.5 * log(0.5);
}

double H_two(int SNP1, int SNP2, int nsamples) {
    // Rprintf("H_two called...\n");
    // Rprintf("nsamples = %d\n", nsamples);
    // Rprintf("in H_two count = %d, CaseCount = %d\n", count, CaseCount);
    // p(0, 0)
    double p0_0 = CompareSNP(SNP1, SNP2, 0) / (double)nsamples;
    // Rprintf("CompareSNP(SNP1, SNP2, 0) = %d\n", CompareSNP(SNP1, SNP2, 0));
    double log_p0_0 = p0_0 == 0 ? 0 : log(p0_0);
    // p(1, 1)
    double p1_1 = CompareSNP(SNP1, SNP2, 1) / (double)nsamples;
    // Rprintf("CompareSNP(SNP1, SNP2, 1) = %d\n", CompareSNP(SNP1, SNP2, 1));
    double log_p1_1 = p1_1 == 0 ? 0 : log(p1_1);
    // p(2, 2)
    double p2_2 = CompareSNP(SNP1, SNP2, 2) / (double)nsamples;
    // Rprintf("CompareSNP(SNP1, SNP2, 2) = %d\n", CompareSNP(SNP1, SNP2, 2));
    double log_p2_2 = p2_2 == 0 ? 0 : log(p2_2);
    // p(0, 1)
    double p0_1 = (comparePartSet(SNP1, 0, SNP2, 1, 0) + comparePartSet(SNP1, 0, SNP2, 1, 1)) / (double)nsamples;
    // Rprintf("comparePartSet(SNP1, 0, SNP2, 1, 0) + comparePartSet(SNP1, 0, SNP2, 1, 1) = %d\n", comparePartSet(SNP1, 0, SNP2, 1, 0) + comparePartSet(SNP1, 0, SNP2, 1, 1));
    double log_p0_1 = p0_1 == 0 ? 0 : log(p0_1);
    // p(0, 2)
    double p0_2 = (comparePartSet(SNP1, 0, SNP2, 2, 0) + comparePartSet(SNP1, 0, SNP2, 2, 1)) / (double)nsamples;
    // printf("comparePartSet(SNP1, 0, SNP2, 2, 0) + comparePartSet(SNP1, 0, SNP2, 2, 1) = %d\n", comparePartSet(SNP1, 0, SNP2, 2, 0) + comparePartSet(SNP1, 0, SNP2, 2, 1));
    double log_p0_2 = p0_2 == 0 ? 0 : log(p0_2);
    // p(1, 0)
    double p1_0 = (comparePartSet(SNP1, 1, SNP2, 0, 0) + comparePartSet(SNP1, 1, SNP2, 0, 1)) / (double)nsamples;
    // printf("comparePartSet(SNP1, 1, SNP2, 0, 0) + comparePartSet(SNP1, 1, SNP2, 0, 1) = %d\n", comparePartSet(SNP1, 1, SNP2, 0, 0) + comparePartSet(SNP1, 1, SNP2, 0, 1));
    double log_p1_0 = p1_0 == 0 ? 0 : log(p1_0);
    // p(1, 2)
    double p1_2 = (comparePartSet(SNP1, 1, SNP2, 2, 0) + comparePartSet(SNP1, 1, SNP2, 2, 1)) / (double)nsamples;
    // printf("comparePartSet(SNP1, 1, SNP2, 2, 0) + comparePartSet(SNP1, 1, SNP2, 2, 1) = %d\n", comparePartSet(SNP1, 1, SNP2, 2, 0) + comparePartSet(SNP1, 1, SNP2, 2, 1));
    double log_p1_2 = p1_2 == 0 ? 0 : log(p1_2);
    // p(2, 0)
    double p2_0 = (comparePartSet(SNP1, 2, SNP2, 0, 0) + comparePartSet(SNP1, 2, SNP2, 0, 1)) / (double)nsamples;
    // printf("comparePartSet(SNP1, 2, SNP2, 0, 0) + comparePartSet(SNP1, 2, SNP2, 0, 1) = %d\n", comparePartSet(SNP1, 2, SNP2, 0, 0) + comparePartSet(SNP1, 2, SNP2, 0, 1));
    double log_p2_0 = p2_0 == 0 ? 0 : log(p2_0);
    // p(2, 1)
    double p2_1 = (comparePartSet(SNP1, 2, SNP2, 1, 0) + comparePartSet(SNP1, 2, SNP2, 1, 1)) / (double)nsamples;
    // printf("comparePartSet(SNP1, 2, SNP2, 1, 0) + comparePartSet(SNP1, 2, SNP2, 1, 1) = %d\n", comparePartSet(SNP1, 2, SNP2, 1, 0) + comparePartSet(SNP1, 2, SNP2, 1, 1));
    double log_p2_1 = p2_1 == 0 ? 0 : log(p2_1);
    // Rprintf("log(p0_0) = %lf\n", log_p0_0);
    // Rprintf("log(p0_1) = %lf\n", log_p0_1);
    // Rprintf("log(p0_2) = %lf\n", log_p0_2);
    // Rprintf("log(p1_0) = %lf\n", log_p1_0);
    // Rprintf("log(p1_1) = %lf\n", log_p1_1);
    // Rprintf("log(p1_2) = %lf\n", log_p1_2);
    // Rprintf("log(p2_0) = %lf\n", log_p2_0);
    // Rprintf("log(p2_1) = %lf\n", log_p2_1);
    // Rprintf("log(p2_2) = %lf\n", log_p2_2);
    return (- (p0_0 * log_p0_0) - (p0_1 * log_p0_1) - (p0_2 * log_p0_2) - (p1_0 * log_p1_0) - (p1_1 * log_p1_1) - (p1_2 * log_p1_2) - (p2_0 * log_p2_0) - (p2_1 * log_p2_1) - (p2_2 * log_p2_2));
}

double H_three(int SNP1, int SNP2, int nsamples) {
    // Rprintf("H_three called...\n");
    // Rprintf("nsamples = %d\n", nsamples);
    // Rprintf("in H_three count = %d, CaseCount = %d\n", count, CaseCount);
    // p(0, 0, 0) Class, SNP1, SNP2
    double p0_0_0 = comparePartSet(SNP1, 0, SNP2, 0, 0) / (double)nsamples;
    // Rprintf("comparePartSet(SNP1, 0, SNP2, 0, 0) = %d\n", comparePartSet(SNP1, 0, SNP2, 0, 0));
    double log_p0_0_0 = p0_0_0 == 0 ? 0 : log(p0_0_0);
    // p(0, 1, 1)
    double p0_1_1 = comparePartSet(SNP1, 1, SNP2, 1, 0) / (double)nsamples;
    // Rprintf("comparePartSet(SNP1, 1, SNP2, 1, 0) = %d\n", comparePartSet(SNP1, 1, SNP2, 1, 0));
    double log_p0_1_1 = p0_1_1 == 0 ? 0 : log(p0_1_1);
    // p(0, 2, 2)
    double p0_2_2 = comparePartSet(SNP1, 2, SNP2, 2, 0) / (double)nsamples;
    // Rprintf("comparePartSet(SNP1, 2, SNP2, 2, 0) = %d\n", comparePartSet(SNP1, 2, SNP2, 2, 0));
    double log_p0_2_2 = p0_2_2 == 0 ? 0 : log(p0_2_2);
    // p(0, 0, 1)
    double p0_0_1 = comparePartSet(SNP1, 0, SNP2, 1, 0) / (double)nsamples;
    // printf("comparePartSet(SNP1, 0, SNP2, 1, 0) = %d\n", comparePartSet(SNP1, 0, SNP2, 1, 0));
    double log_p0_0_1 = p0_0_1 == 0 ? 0 : log(p0_0_1);
    // p(0, 0, 2)
    double p0_0_2 = comparePartSet(SNP1, 0, SNP2, 2, 0) / (double)nsamples;
    // printf("comparePartSet(SNP1, 0, SNP2, 2, 0) = %d\n", comparePartSet(SNP1, 0, SNP2, 2, 0));
    double log_p0_0_2 = p0_0_2 == 0 ? 0 : log(p0_0_2);
    // p(0, 1, 0)
    double p0_1_0 = comparePartSet(SNP1, 1, SNP2, 0, 0) / (double)nsamples;
    // printf("comparePartSet(SNP1, 1, SNP2, 0, 0) = %d\n", comparePartSet(SNP1, 1, SNP2, 0, 0));
    double log_p0_1_0 = p0_1_0 == 0 ? 0 : log(p0_1_0);
    // p(0, 1, 2)
    double p0_1_2 = comparePartSet(SNP1, 1, SNP2, 2, 0) / (double)nsamples;
    // printf("comparePartSet(SNP1, 1, SNP2, 2, 0) = %d\n", comparePartSet(SNP1, 1, SNP2, 2, 0));
    double log_p0_1_2 = p0_1_2 == 0 ? 0 : log(p0_1_2);
    // p(0, 2, 0)
    double p0_2_0 = comparePartSet(SNP1, 2, SNP2, 0, 0) / (double)nsamples;
    // printf("comparePartSet(SNP1, 2, SNP2, 0, 0) = %d\n", comparePartSet(SNP1, 2, SNP2, 0, 0));
    double log_p0_2_0 = p0_2_0 == 0 ? 0 : log(p0_2_0);
    // p(0, 2, 1)
    double p0_2_1 = comparePartSet(SNP1, 2, SNP2, 1, 0) / (double)nsamples;
    // printf("comparePartSet(SNP1, 2, SNP2, 1, 0) = %d\n", comparePartSet(SNP1, 2, SNP2, 1, 0));
    double log_p0_2_1 = p0_2_1 == 0 ? 0 : log(p0_2_1);
    // p(1, 0, 0) Class, SNP1, SNP2
    double p1_0_0 = comparePartSet(SNP1, 0, SNP2, 0, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 0, SNP2, 0, 1) = %d\n", comparePartSet(SNP1, 0, SNP2, 0, 1));
    double log_p1_0_0 = p1_0_0 == 0 ? 0 : log(p1_0_0);
    // p(1, 1, 1)
    double p1_1_1 = comparePartSet(SNP1, 1, SNP2, 1, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 1, SNP2, 1, 1) = %d\n", comparePartSet(SNP1, 1, SNP2, 1, 1));
    double log_p1_1_1 = p1_1_1 == 0 ? 0 : log(p1_1_1);
    // p(1, 2, 2)
    double p1_2_2 = comparePartSet(SNP1, 2, SNP2, 2, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 2, SNP2, 2, 1) = %d\n", comparePartSet(SNP1, 2, SNP2, 2, 1));
    double log_p1_2_2 = p1_2_2 == 0 ? 0 : log(p1_2_2);
    // p(1, 0, 1)
    double p1_0_1 = comparePartSet(SNP1, 0, SNP2, 1, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 0, SNP2, 1, 1) = %d\n", comparePartSet(SNP1, 0, SNP2, 1, 1));
    double log_p1_0_1 = p1_0_1 == 0 ? 0 : log(p1_0_1);
    // p(1, 0, 2)
    double p1_0_2 = comparePartSet(SNP1, 0, SNP2, 2, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 0, SNP2, 2, 1) = %d\n", comparePartSet(SNP1, 0, SNP2, 2, 1));
    double log_p1_0_2 = p1_0_2 == 0 ? 0 : log(p1_0_2);
    // p(1, 1, 0)
    double p1_1_0 = comparePartSet(SNP1, 1, SNP2, 0, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 1, SNP2, 0, 1) = %d\n", comparePartSet(SNP1, 1, SNP2, 0, 1));
    double log_p1_1_0 = p1_1_0 == 0 ? 0 : log(p1_1_0);
    // p(1, 1, 2)
    double p1_1_2 = comparePartSet(SNP1, 1, SNP2, 2, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 1, SNP2, 2, 1) = %d\n", comparePartSet(SNP1, 1, SNP2, 2, 1));
    double log_p1_1_2 = p1_1_2 == 0 ? 0 : log(p1_1_2);
    // p(1, 2, 0)
    double p1_2_0 = comparePartSet(SNP1, 2, SNP2, 0, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 2, SNP2, 0, 1) = %d\n", comparePartSet(SNP1, 2, SNP2, 0, 1));
    double log_p1_2_0 = p1_2_0 == 0 ? 0 : log(p1_2_0);
    // p(1, 2, 1)
    double p1_2_1 = comparePartSet(SNP1, 2, SNP2, 1, 1) / (double)nsamples;
    // printf("comparePartSet(SNP1, 2, SNP2, 1, 1) = %d\n", comparePartSet(SNP1, 2, SNP2, 1, 1));
    double log_p1_2_1 = p1_2_1 == 0 ? 0 : log(p1_2_1);
    return (- (p0_0_0 * log_p0_0_0) - (p0_0_1 * log_p0_0_1) - (p0_0_2 * log_p0_0_2) - (p0_1_0 * log_p0_1_0) - (p0_1_1 * log_p0_1_1) - (p0_1_2 * log_p0_1_2) - (p0_2_0 * log_p0_2_0) - (p0_2_1 * log_p0_2_1) - (p0_2_2 * log_p0_2_2) - (p1_0_0 * log_p1_0_0) - (p1_0_1 * log_p1_0_1) - (p1_0_2 * log_p1_0_2) - (p1_1_0 * log_p1_1_0) - (p1_1_1 * log_p1_1_1) - (p1_1_2 * log_p1_1_2) - (p1_2_0 * log_p1_2_0) - (p1_2_1 * log_p1_2_1) - (p1_2_2 * log_p1_2_2));
}

// dev2 branch, has to convert R numeric matrix to C 0/1 array in every arc MI calc, time-consuming.
// extract convertion part to a new function `num2binary` in dev3 branch
// double mi2(const int *d, int varNum1, int varNum2, int varNum3, int nsamples, int nvars) {
//     // Rprintf("mi2 called...\n");
//     // Rprintf("nsamples, nvars = %d, %d\n", nsamples, nvars);
//     int i, j; 
//     count = 0; 
//     CaseCount = 0;
    
//     // row Class, first row
//     for(i = 0; i < nsamples; i++){
//         // Rprintf("d[CMC(varNum1-1, i, nvars)]: %d\n", d[CMC(varNum1-1, i, nvars)]);
//         if(d[CMC(varNum1-1, i, nvars)] >= 0 && d[CMC(varNum1-1, i, nvars)] <= 1){
//             // Rprintf("in class _if_ block\n");
//             count++;
//             if(d[CMC(varNum1-1, i, nvars)] == 1) CaseCount++;
//         }
//     }
//     // Rprintf("in mi2 count = %d, CaseCount = %d\n", count, CaseCount);
//     //ControlCount = count - CaseCount;
//     //printf("the CaseCount is %d\n", CaseCount);
//     SNPnumber = 0;
//     //printf("count is %d\n", count);
    
//     // process the rest of the vars
//     // i starts from 1, skip row Class
//     for(i = 1; i < nvars; i++) {
//         int k = 0;
//         for(j = 0; j < nsamples; j++) {
//             if(d[CMC(i, j, nvars)] >= 0 && d[CMC(i, j, nvars)] <= 2) {
//                 // Rprintf("data[%d][%d] = 1\n", SNPnumber, (d[CMC(i, j, nvars)] - 0)*count + k);
//                 data[SNPnumber][(d[CMC(i, j, nvars)] - 0)*count + k] = 1;
//                 k++;
//             }
//         }
//         SNPnumber++;
//     }
    
//     // Rprintf("data matrix presented below...\n");
//     // Rprintf("-----------------------------------------------------------\n");
//     // j = 0;
//     // for(i = 0; i < SNPnumber; i++){
//     //     for(j = 0; j < 3 * count; j++){
//     //         Rprintf("%d-", data[i][j]);
//     //     }
//     //     Rprintf("\n");
//     // }

//     // Rprintf("H_one is %lf\n", H_one());
//     // Rprintf("H_two is %lf\n", H_two(varNum2, varNum3, nsamples));
//     // Rprintf("H_three is %lf\n", H_three(varNum2, varNum3, nsamples));
//     // Rprintf("return is %lf\n", H_one() + H_two(varNum2, varNum3, nsamples) - H_three(varNum2, varNum3, nsamples));
//     return H_one() + H_two(varNum2, varNum3, nsamples) - H_three(varNum2, varNum3, nsamples);
// }

void num2binary(const int *d, int varNum1, int nsamples, int nvars) {
    // Rprintf("dev3@branch: num2binary called...\n");
    int i, j; 
    count = 0; 
    CaseCount = 0;
    SNPnumber = 0;
    memset(data, 0, sizeof(data));
    // row Class, first row
    for(i = 0; i < nsamples; i++){
        // Rprintf("d[CMC(varNum1-1, i, nvars)]: %d\n", d[CMC(varNum1-1, i, nvars)]);
        if(d[CMC(varNum1-1, i, nvars)] >= 0 && d[CMC(varNum1-1, i, nvars)] <= 1){
            // Rprintf("in class _if_ block\n");
            count++;
            if(d[CMC(varNum1-1, i, nvars)] == 1) CaseCount++;
        }
    }
    // Rprintf("in mi2 count = %d, CaseCount = %d\n", count, CaseCount);
    //ControlCount = count - CaseCount;
    //printf("the CaseCount is %d\n", CaseCount);
    //printf("count is %d\n", count);
    
    // process the rest of the vars
    // i starts from 1, skip row Class
    for(i = 1; i < nvars; i++) {
        int k = 0;
        for(j = 0; j < nsamples; j++) {
            if(d[CMC(i, j, nvars)] >= 0 && d[CMC(i, j, nvars)] <= 2) {
                // Rprintf("data[%d][%d] = 1\n", SNPnumber, (d[CMC(i, j, nvars)] - 0)*count + k);
                data[SNPnumber][(d[CMC(i, j, nvars)] - 0)*count + k] = 1;
                k++;
            }
        }
        SNPnumber++;
    }
}

double mi2(int varNum2, int varNum3, int nsamples) {
    // Rprintf("dev3@branch: mi2 called...\n");
    return H_one() + H_two(varNum2, varNum3, nsamples) - H_three(varNum2, varNum3, nsamples);
}

// dev2 branch, has to convert R numeric matrix to C 0/1 array in every arc MI calc, time-consuming.
// extract convertion part to a new function `num2binary` in dev3 branch
// SEXP mi2R (SEXP Rdata, SEXP RvarNum1, SEXP RvarNum2, SEXP RvarNum3, SEXP Rncols, SEXP Rnrows)
// {
//     const int *rdata;
//     const int *varNum1, *varNum2, *varNum3;
//     const int *nrows, *ncols;
//     SEXP res;

//     PROTECT(Rdata = AS_INTEGER(Rdata));
//     PROTECT(RvarNum1 = AS_INTEGER(RvarNum1));
//     PROTECT(RvarNum2 = AS_INTEGER(RvarNum2));
//     PROTECT(RvarNum3 = AS_INTEGER(RvarNum3));
//     PROTECT(Rnrows= AS_INTEGER(Rnrows));
//     PROTECT(Rncols= AS_INTEGER(Rncols));
//     PROTECT(res = NEW_NUMERIC(1));

//     rdata = INTEGER_POINTER(Rdata);
//     varNum1 = INTEGER_POINTER(RvarNum1);
//     varNum2 = INTEGER_POINTER(RvarNum2);
//     varNum3 = INTEGER_POINTER(RvarNum3);
//     nrows= INTEGER_POINTER(Rnrows);
//     ncols= INTEGER_POINTER(Rncols);
    
//     REAL(res)[0] = mi2(rdata, *varNum1, *varNum2, *varNum3, *ncols, *nrows);
//     UNPROTECT(7);
//     return res;
// }

SEXP mi2R (SEXP RvarNum2, SEXP RvarNum3, SEXP Rncols)
{
    // dev3@branch: backend for arc MI calc
    const int *varNum2, *varNum3;
    const int *ncols;
    SEXP res;

    PROTECT(RvarNum2 = AS_INTEGER(RvarNum2));
    PROTECT(RvarNum3 = AS_INTEGER(RvarNum3));
    PROTECT(Rncols= AS_INTEGER(Rncols));
    PROTECT(res = NEW_NUMERIC(1));

    varNum2 = INTEGER_POINTER(RvarNum2);
    varNum3 = INTEGER_POINTER(RvarNum3);
    ncols= INTEGER_POINTER(Rncols);
    
    REAL(res)[0] = mi2(*varNum2, *varNum3, *ncols);
    UNPROTECT(4);
    return res;
}

void num2binaryR (SEXP Rdata, SEXP RvarNum1, SEXP Rncols, SEXP Rnrows)
{
    // dev3@branch: backend for R numeric matrix to C 0/1 array transformation
    const int *rdata;
    const int *nrows, *ncols;
    const int *varNum1;

    PROTECT(Rdata = AS_INTEGER(Rdata));
    PROTECT(RvarNum1 = AS_INTEGER(RvarNum1));
    PROTECT(Rnrows= AS_INTEGER(Rnrows));
    PROTECT(Rncols= AS_INTEGER(Rncols));

    rdata = INTEGER_POINTER(Rdata);
    varNum1 = INTEGER_POINTER(RvarNum1);
    nrows= INTEGER_POINTER(Rnrows);
    ncols= INTEGER_POINTER(Rncols);
    
    num2binary(rdata, *varNum1, *ncols, *nrows);
    UNPROTECT(4);
}
