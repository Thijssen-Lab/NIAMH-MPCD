#ifndef CTOOLS_H
#define CTOOLS_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are nonspecific routines which
   I find useful when coding.
*/

void zerowarning( double test, double zero,double mult );
void wait4u( );

// Core routines for NaN debugging
int isNaN(double x);
int isNaNs(double *x, int n);

// Routines to set various forms of arrays to be zero'd
void zerovec( double VEC[],int dimension );
void zeromat( int dim1, int dim2, double **MAT);
void zerovec_v(int count, int dim, ...); // variadic
void zeromat_v(int count, int dim1, int dim2, ...); // variadic

#endif
