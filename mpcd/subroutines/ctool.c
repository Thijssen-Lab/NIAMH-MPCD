///
/// @file
///

# include<stdio.h>
# include<math.h>

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ HELPFUL ROUTINES ************ */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
///
/// @brief		Waits for the user to press enter.
///		
void wait4u() {
	printf( "\nPress enter to continue\n" );
	while ( 1 ) {
		if( '\n' == getchar() )
		break;
	}
}

///
/// @brief			    Prints a warning message if the value is within the precision cutoff for zero. 
///
/// @param test			The value that will be checked.
/// @param zero 		The check value for zero.
/// @param mult 		The coefficient that changes the magnitude of `zero`.
///
void zerowarning( double test, double zero,double mult ) {
	if( fabs( test ) <= zero*mult ) printf( "Warning: Value %e is within %lf times magnitude of zero check value of %e.\n",test,mult,zero );
}
