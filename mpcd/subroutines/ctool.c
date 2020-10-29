# include<stdio.h>
# include<math.h>

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ HELPFUL ROUTINES ************ */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void wait4u() {
/*
    Waits for the user to press enter
*/
	printf( "\nPress enter to continue\n" );
	while ( 1 ) {
		if( '\n' == getchar() )
		break;
	}
}
void zerowarning( double test, double zero,double mult ) {
/*
    Prints a warning message if the value is
    within the precision cutoff for zero
*/
	if( fabs( test ) <= zero*mult ) printf( "Warning: Value %e is within %lf times magnitude of zero check value of %e.\n",test,mult,zero );
}
