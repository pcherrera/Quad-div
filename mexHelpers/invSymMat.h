#ifndef _INV_SYM_MAT_H_
#define _INV_SYM_MAT_H_

/* THe below defined function calculates inverse of a symmetric positive function */

/* DEFINE LAPACK routine names on LINUX machines! */
#if !defined(_WIN32) && !defined(_WIN64)
#define dpotrf dpotrf_
#define dpotri dpotri_
#endif

#include "mex.h"
#include "blas.h"

/* Function declarations */
void dpotrf( char*, mwSize*, double*, mwSize*, mwSignedIndex* );
void dpotri( char*, mwSize*, double*, mwSize*, mwSignedIndex* );

/* function that computes inverse, note that input data is overwritten */
inline void invSymMat(mwSize n, double * M) {
    char *up = "U"; /* Upper triangle matrix */
    mwSignedIndex info; /* OUTPUT FLAG */
    mwIndex i, j;
    
    dpotrf( up, &n, M, &n, &info ); /* Cholesky decomposition */
    dpotri( up, &n, M, &n, &info ); /* Inverse using the before computed Chol.decomp. */
    
    
    /* fill lower triangle */
    for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
            M[n*i+j]=M[j*n+i];
        }
    }
}

#endif