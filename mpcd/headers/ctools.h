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

int isNaN(double x);
int isNaNs(double *x, int n);

void zerovec( double VEC[],int dimension );
void zeromat( int dim1, int dim2, double **MAT);

#endif
