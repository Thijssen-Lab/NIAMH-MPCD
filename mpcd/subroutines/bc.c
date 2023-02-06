///
/// @file
///
/// @brief This file applies boundary conditions (bc) to mpcd particles.
///
/// This file applies boundary conditions (bc) to mpcd particles crossing through a boundary.
///

# include<stdio.h>
# include<math.h>
# include<time.h>
# include<string.h>
# include<stdlib.h>

# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/SRDclss.h"
# include "../headers/rand.h"
# include "../headers/pout.h"
# include "../headers/ctools.h"
# include "../headers/mtools.h"
# include "../headers/pout.h"
# include "../headers/mpc.h"
# include "../headers/therm.h"
# include "../headers/bc.h"
# include "../headers/init.h"
# include "../headers/lc.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ********** PARTICLES PASSING BCs ********* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief TODO
///
/// This routine calculates the W parameter which determines whether the particle
/// has crossed through the boundary or not.
///	TODO mention what the W values mean.
///
/// @param WALL The boundary.
/// @param P The individual mpcd particle.
/// @see calcWavyW()
/// @see non4foldSymmCalcW()
/// @return Returns W.
///
double calcW( bc WALL,particleMPC P ) {
/*
   This function calculates W which is used to
   determine if boundary conditions should be
   applied to a particleMPC.
*/
	double terms, W=0.0;
	int i;

	if( feq(WALL.ROTSYMM[0],4.0) && feq(WALL.ROTSYMM[1],4.0) ) {
		for( i=0; i<DIM; i++ ) {
			terms = WALL.A[i] * ( P.Q[i]-WALL.Q[i] );
			if( WALL.ABS ) terms=fabs(terms);
			terms = smrtPow( terms,WALL.P[i] );
			W += terms;
		}
		terms = WALL.R;
		if( WALL.ABS ) terms=fabs(terms);
		terms = smrtPow( terms,WALL.P[3] );
		W -= terms;
		//Check if need wavy wall complications
		if( !feq(WALL.B[0],0.0) ) W += calcWavyW(WALL,P.Q,W);
		//Check if invert wall
		if( WALL.INV ) W *= -1.0;
	}
	else {
		W = non4foldSymmCalcW( WALL,P.Q,DIM );
	}
	return W;
}

///
/// @brief TODO.
///
/// TODO.
///
/// @param WALL The boundary.
/// @param POS[] The position of the particle.
/// @param W[] The W value.
/// @see calcW()
/// @return Returns TODO.
///
double calcWavyW( bc WALL,double POS[], double W ) {
/*
   This function calculates additions to W for wavy
   walls in the calcW function.
*/
	double W1=0.0,W2=0.0;
	int i,flag;
	//Planes
	flag=0;
	for( i=0; i<DIM; i++ ) if( !feq(WALL.P[i],1.0) ) flag+=1;
	if(!flag){
		if ( !feq(WALL.B[1],0.0) ) {
			W1 = (WALL.A[1]*POS[0]-WALL.A[0]*POS[1]) / sqrt( WALL.A[0]*WALL.A[0] + WALL.A[1]*WALL.A[1] );
		}
		if( DIM>2 && !feq(WALL.B[2],0.0) ){
			W2 = WALL.A[0]*WALL.A[2]*POS[0] + WALL.A[1]*WALL.A[2]*POS[1] - (WALL.A[0]*WALL.A[0]+WALL.A[1]*WALL.A[1])*POS[2];
			W2 /= sqrt( pow(WALL.A[0],4) + pow(WALL.A[1],4) + pow(WALL.A[0],2)*pow(WALL.A[2],2) + 2.0*pow(WALL.A[0],2)*pow(WALL.A[1],2) + pow(WALL.A[1],2)*pow(WALL.A[2],2) );
		}
	}
	else {
		//Ellipsoidal
		flag = 0;
		for( i=0; i<DIM; i++ ) if( !feq(WALL.P[i],2.0) || feq(WALL.A[i],0.0) ) flag+=1;
		if(!flag){
			if ( !feq(WALL.B[1],0.0) ) 	W1 = atan2( (WALL.A[1] * ( POS[1]-WALL.Q[1] ) ),(WALL.A[0] * ( POS[0]-WALL.Q[0] ) ) );
			if( DIM>2 && !feq(WALL.B[2],0.0) ){
				W2 = atan2( sqrt( pow(WALL.A[0] * ( POS[0]-WALL.Q[0] ),2) + pow(WALL.A[1] * ( POS[1]-WALL.Q[1] ),2) ) , ( WALL.A[2] * ( POS[2]-WALL.Q[2] ) ) );
			}
		}
		else if(DIM>2 && flag==1){
			//Cylinders
			flag = -1;
			for( i=0; i<DIM; i++) if ( feq(WALL.A[i],0.0) ) flag=i;
			if (flag!=-1) {
				int j, k, l;
				if (flag==0) j = 0, k = 1, l = 2;
				else if (flag==1) j = 1, k = 2, l = 0;
				else j = 2, k = 0, l = 1;
				if ( !feq(WALL.B[1],0.0) ) W1 = atan2( (WALL.A[l] * ( POS[l]-WALL.Q[l] ) ),(WALL.A[k] * ( POS[k]-WALL.Q[k] ) ) );
				if ( !feq(WALL.B[2],0.0) ) W2 = atan2( sqrt( W + pow( WALL.R,WALL.P[3] ) ), sqrt( POS[j]) );
			}
		}
	}

	return WALL.B[0] * cos( WALL.B[1]*W1 ) * cos( WALL.B[2]*W2 );
}

///
/// @brief TODO.
///
/// TODO.
///
/// @param movingWall TODO.
/// @param stillWall TODO.
/// @param flagCentre TODO.
/// @see BC_BCcollision()
/// @return Returns W.
///
double calcW_BC( bc movingWall,bc stillWall,int flagCentre ) {
/*
   This function calculates W which is used to
   determine if boundary conditions should be
   applied to a moving BC.
	 intFlag determines if we are considering the centre of the movingWall or the surface
*/
	double terms, modR, W=0.0;

	if( feq(stillWall.ROTSYMM[0],4.0) && feq(stillWall.ROTSYMM[1],4.0) ) {
		W = surf_func(stillWall, movingWall.Q, DIM);
		if( !flagCentre ) {
			terms = stillWall.R;
			if( stillWall.ABS ) terms=fabs(terms);
			terms = smrtPow( terms,stillWall.P[3] );
			if( stillWall.INV ) terms *= -1.0; //because W already inverted
			W += terms; //adds back R term
			//subtracts the correct R terms
			modR = stillWall.R + movingWall.R;
			if( stillWall.INV ) modR *= -1.0;
			W -= modR;
		}
	}
	else {
		printf( "Error:\tNon 4-fold symmetry not yet programmed for BC-BC interactions\n" );
		exit( 1 );
	}

	return W;
}

///
/// @brief Calculates the particle crosstime.
///
/// This routine calculates the particle crosstime.
/// Since this particle is found inside the boundary, the particle must have
/// crossed the boundary before the end of the streaming step.
/// This method interpolates back the path taken by the particle and finds this
/// time when the particle crossed the boundary.
/// TODO : the above might not be right. It might be forward interpolation, but do we even need this?
///
/// @param p The individual mpcd particle.
/// @param WALL The boundary.
/// @param tc_pos TODO.
/// @param tc_neg TODO.
/// @param t_step The time step interval.
/// @see secant_time()
///
void crosstime( particleMPC p,bc WALL,double *tc_pos, double *tc_neg,double t_step ) {
/*
    Calculate when the particleMPC crosses the bc
    by solving the trajectory equation

*/
	double a=0.0,b=0.0,c=0.0;
	int i;

	// Planar Wall
	if( WALL.PLANAR || ( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) && feq(WALL.P[3],1.0) ) ) {
		*tc_pos = WALL.R;
		for( i=0; i<DIM; i++) *tc_pos += WALL.A[i]*(WALL.Q[i]-p.Q[i]);
		*tc_pos /= (WALL.A[0]*p.V[0] + WALL.A[1]*p.V[1] + WALL.A[2]*p.V[2]);
		//There is only one time.
		*tc_neg = *tc_pos;
	}
	// Ellipsoid
	//else if( feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0) ) {
	else if ((DIM == 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0)) || (DIM > 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0))){
	//else if( feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0)) {
		for( i=0; i<DIM; i++ ) {
			a += WALL.A[i]*WALL.A[i]*p.V[i]*p.V[i];
			b += WALL.A[i]*WALL.A[i]*p.V[i]*(p.Q[i]-WALL.Q[i]);
			c += WALL.A[i]*WALL.A[i]*(p.Q[i]*p.Q[i]-2.0*p.Q[i]*WALL.Q[i]+WALL.Q[i]*WALL.Q[i]);
		}
		b *= 2.0;
		c -= smrtPow(WALL.R,WALL.P[3]);
		//Use the quadratic formula
		*tc_neg = - b - sqrt(b*b-4.0*a*c);
		*tc_neg /= (2.0*a);
		*tc_pos = - b + sqrt(b*b-4.0*a*c);
		*tc_pos /= (2.0*a);
	}
	else {
		//Must use secant method to determine cross times
		*tc_pos = secant_time( p,WALL,t_step );
		*tc_neg = *tc_pos;
	}
}

///
/// @brief Calculates the particle crosstime.
///
/// This routine calculates the particle crosstime.
/// Since this particle is found inside the boundary, the particle must have
/// crossed the boundary before the end of the streaming step.
/// This method interpolates back the path taken by the particle and finds this
/// time when the particle crossed the boundary.
///
/// @param p The individual mpcd particle.
/// @param WALL The boundary.
/// @param tc_pos TODO.
/// @param tc_neg TODO.
/// @param t_step The time step interval.
/// @see chooseBC()
/// @see secant_time()
///
void crosstimeReverse( particleMPC p,bc WALL,double *tc_pos, double *tc_neg,double t_step ) {
/*
    Calculate when the particleMPC crossed the bc
    by solving the trajectory equation using negative velocity

*/
	double a=0.0,b=0.0,c=0.0;
	int i;

	// Planar Wall
	if ( ( WALL.PLANAR ) || ( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) && feq(WALL.P[3],1.0) && feq(WALL.B[0],0.0) ) ){
		*tc_pos = WALL.R;
		for( i=0; i<DIM; i++) *tc_pos += WALL.A[i]*(WALL.Q[i]-p.Q[i]);
		*tc_pos /= (WALL.A[0]*(-p.V[0]) + WALL.A[1]*(-p.V[1]) + WALL.A[2]*(-p.V[2]));
		//There is only one time.
		*tc_neg = *tc_pos;
	}
	// Ellipsoid
	//else if( feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0) ) {
	else if ((DIM == 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.B[0],0.0) ) || (DIM > 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0) && feq(WALL.B[0],0.0) )){
		for( i=0; i<DIM; i++ ) {
			a += WALL.A[i]*WALL.A[i]*(-p.V[i])*(-p.V[i]);
			b += WALL.A[i]*WALL.A[i]*(-p.V[i])*(p.Q[i]-WALL.Q[i]);
			c += WALL.A[i]*WALL.A[i]*(p.Q[i]*p.Q[i]-2.*p.Q[i]*WALL.Q[i]+WALL.Q[i]*WALL.Q[i]);
		}
		b *= 2.0;
		c -= smrtPow(WALL.R,WALL.P[3]);
		//Use the quadratic formula
		*tc_neg = - b - sqrt(b*b-4.*a*c);
		*tc_neg /= (2.*a);
		*tc_pos = - b + sqrt(b*b-4.*a*c);
		*tc_pos /= (2.*a);
	}
	else {
		//Must use secant method to determine cross times
		*tc_pos = secant_time( p,WALL,t_step );
		//If secant method fails rewind the full timestep
		if(*tc_pos<0 || *tc_pos > t_step) {
				*tc_pos = t_step;
		}
		*tc_neg = *tc_pos;
	}
}

///
/// @brief Calculates W for a boundary identified as a plane.
///
/// TODO
///
/// @param WALL The boundary.
/// @param P The individual mpcd particle.
/// @return Returns W.
///
double calcW_PLANE( bc WALL,particleMPC P ) {
/*
   This function calculates W like calcW BUT assumes it is a plane
   i.e. faster
*/
	double W=0.;
	int i;


	for( i=0; i<DIM; i++ ) W += WALL.A[i] * ( P.Q[i]-WALL.Q[i] );
	W -= WALL.R;
	if( WALL.INV ) W *= -1.;
	return W;
}

///
/// @brief TODO
///
/// TODO
///
/// @param p The individual mpcd particle.
/// @param WALL The boundary.
/// @param t_step The time step interval.
/// @see crosstimeReverse()
/// @return TODO
///
double secant_time( particleMPC p,bc WALL,double t_step ) {
	double Qi[DIM],QiM1[DIM];
	double ti,tiM1,root;
	double fi,fiM1;
	int i,iter=0;


	//Rewind the particle back to it's old position
	for( i=0;i<DIM;i++ ) {
		QiM1[i] = p.Q[i];
		Qi[i] = trans(t_step,-p.V[i],p.Q[i]);
	}

	tiM1 = 0.0;
	ti = t_step;
	root = ti;

	//Secant Loop
	do {
		iter++;
		// Calculate the surface function for the particles' positions at these times
		fi = surf_func( WALL,Qi,DIM );
		fiM1 = surf_func( WALL,QiM1,DIM );

		if ( !feq(fi, fiM1) ) root = ti - fi * ( ti - tiM1 )/ ( fi - fiM1 );

		//updates values so always on opposite sides of 0
		if (root<TOL/100.) {
			ti = tiM1;
			tiM1 = root;
		}
		else {
			tiM1 = ti;
			ti = root;
		}

		// Calculate the particles' positions at these times
		for( i=0;i<DIM;i++ ) {
			QiM1[i] = trans(tiM1,-p.V[i],p.Q[i]);
			Qi[i] = trans(ti,-p.V[i],p.Q[i]);
		}
	} while( iter<10 && fabs( ti-tiM1 ) > TOL/100. &&  fabs( fi-fiM1 ) > TOL/100. );
	return root;
}

///
/// @brief Shifts the boundary according to periodic boundary conditions.
///
/// Shift the position of boundary according to periodic boundary conditions (PBC).
/// This follows the steps:
/// - Finds whether to shift (if the particle - wall separation is closer after a full PBC shift).
/// - Calculates the shift as the relevant length of the system (x,y or z direction).
/// - Applies the shift and updates the wall position.
/// @param shift The amount that the wall shifted.
/// @param WALL The boundary.
/// @param pp The individual particle.
/// @see shiftbackBC()
/// @note The shift value is returned through the `shift` variable so that the original
/// boundary position can be restored later using the shiftbackBC() method.
///
void shiftBC( double *shift,bc *WALL,particleMPC *pp ) {

	int k;

	// Zero
	for( k=0; k<DIM; k++ ) shift[k] = 0.0;
	//Don't shift planar surfaces
	//if( fneq(WALL->P[0],1.0) && fneq(WALL->P[1],1.0) && fneq(WALL->P[2],1.0) ) {
	if( fneq(WALL->P[0],1.0) && fneq(WALL->P[1],1.0) ) {
		for( k=0; k<DIM; k++ ) {
			// Determine if should shift
			if( pp->Q[k] - WALL->Q[k] >= 0.5*(double)XYZ[k] ) shift[k] = (double)XYZ[k];
			else if( pp->Q[k] - WALL->Q[k] <= -0.5*(double)XYZ[k] ) shift[k] = -1.0*(double)XYZ[k];
			// Shift BCs
			WALL->Q[k] += shift[k];
		}
	}
}

///
/// @brief Returns the boundary back to original position.
///
/// The boundary has been shifted by the periodic boundary conditions previously in shiftBC().
/// This routine uses the shift value applied previously to return the wall back to the original position.
/// @param shift The amount that the wall shifted.
/// @param WALL The boundary.
/// @see shiftBC()
///
void shiftbackBC( double *shift,bc *WALL ) {
/*
     Shifts the BC back.
*/
	int i;
	for( i=0; i<DIM; i++ ) WALL->Q[i] -= shift[i];
}

///
/// @brief Rotates the boundary.
///
/// If the boundary has a listed non-zero orientation, then it must be rotated to match that orientation.
/// This is performed by rotating the surrounding fluid particles (including
/// positions, velocities and orientations) about the boundary itself.
/// This routine does the rotation and rotation back by having a sign passed to it.
/// @param WALL The boundary.
/// @param pp The individual particle.
/// @param sign The sign multiplier for rotation foward and backwards.
/// @param LC The flag for the fluid being liquid crystalline.
/// @see rotateBC()
/// @see rotatebackBC()
/// @note The `sign` is negative (so the rotation angle will be negative) for the
/// forward particle rotation, since it's the particles rather than the boundary being rotated.
/// Correspondingly the `sign` is positive for the rotation back.
///
void MPC_BCrotation( bc *WALL,particleMPC *pp, double sign, int LC ) {

	int i;
	double rotM[_3D][_3D];		//The rotation matrix
	double oldQ[_3D] = {0.0,0.0,0.0};				//Particle position relative to the centre of the BC
	double newQ[_3D] = {0.0,0.0,0.0};				//Particle position relative to the centre of the BC
	double ax[_3D] = {1.0,0.0,0.0};	//x-axis
	double ay[_3D] = {0.0,1.0,0.0};	//y-axis
	double az[_3D] = {0.0,0.0,1.0};	//z-axis

	#ifdef DBG
		if( DBUG==DBGMPCBC || DBUG==DBGBCMPC) {
			printf( "Rotate BC\n" );
		}
	#endif
	for( i=0; i<DIM; i++ ) oldQ[i] = pp->Q[i]-WALL->Q[i];
	if( DIM>_2D ) setRotMatrix3D( rotM,sign*WALL->O[0],sign*WALL->O[1],sign*WALL->O[2] );
	else setRotMatrix2D( rotM,sign*WALL->O[2] );
	//Rotate the position into place
	dotprodMatVec( rotM,oldQ,newQ,DIM );
	for( i=0; i<DIM; i++ ) pp->Q[i] = WALL->Q[i] + newQ[i];
	//Rotate the velocity vector about each axis
	//Z-axis
	rodriguesRotation( pp->V,az,sign*WALL->O[2] );
	//X- & Y-axes
	if( DIM>_2D ) {
		rodriguesRotation( pp->V,ay,sign*WALL->O[1] );
		rodriguesRotation( pp->V,ax,sign*WALL->O[0] );
	}
	//Rotate the LC orientation/director about each axis
	if(LC) {
		//Z-axis
		rodriguesRotation( pp->U,az,sign*WALL->O[2] );
		//X- & Y-axes
		if( DIM>_2D ) {
			rodriguesRotation( pp->U,ay,sign*WALL->O[1] );
			rodriguesRotation( pp->U,ax,sign*WALL->O[0] );
		}
	}
}

///
/// @brief Rotates the boundary according to the boundary orientation.
///
/// Applies the forward boundary rotation by instead rotating the particles about the boundary (in the reverse direction).
/// @param WALL The boundary.
/// @param pp The individual particle.
/// @param LC The flag for the fluid being liquid crystalline.
/// @see MPC_BCrotation()
/// @see rotatebackBC()
///
void rotateBC( bc *WALL,particleMPC *pp, int LC ) {

	// NOTICE: The current implementation is very wasteful. Every ***particle***
	// is rotated about the centre of each BC.
	// While this is simplest, there are very many particles.

	if(WALL->REORIENT) MPC_BCrotation( WALL,pp,-1.0,LC );
}

///
/// @brief Returns the boundary back to original alignment.
///
/// Applies the backwards boundary rotation (or undoes rotations from rotateBC())
/// by rotating the particles back to original position, velocities and orientations.
/// @param WALL The boundary.
/// @param pp The individual particle.
/// @param LC The flag for the fluid being liquid crystalline.
/// @see MPC_BCrotation()
/// @see rotateBC()
///
void rotatebackBC( bc *WALL,particleMPC *pp, int LC ) {
	if(WALL->REORIENT) MPC_BCrotation( WALL,pp,1.0,LC );
}

///
/// @brief Calculates the magnitude of the impulse
///
/// Calculates the impulse magnitude for the collision of body1 and body2.
/// @param n Unit vector direction of impulse.
/// @param V Velocity for body1.
/// @param U Velocity for body2.
/// @param Pv Centre of mass of body1.
/// @param Pu Centre of mass of body2.
/// @param Wv Angular velocity of body1.
/// @param Wu Angular velocity of body2.
/// @param invMv Total inverse mass of body1.
/// @param invMu Total inverse mass of body2.
/// @param Iv Inverse moment of inertia tensor of body1.
/// @param Iu Inverse moment of inertia tensor of body2.
/// @param r Point of contact.
/// @param E Coefficient of restitution. E=0 for inelastic, E=1 for elastic.
/// @return Impulse magnitude
///
double impulse( double n[],double V[],double U[],double Pv[], double Pu[],double Wv[],double Wu[],double invMv,double invMu,double Iv[][_3D],double Iu[][_3D],double r[],double E ) {

	double J, denom;		//The conservation constant and its denominator
	double Rv[_3D],Ru[_3D];		//Centre wrt point of contact
	double Rxn_v[_3D],Rxn_u[_3D];	//The cross product of vector to contact point and normal
	double irn[_3D];		//The inverse moment of inertia tensor dotted into the cross product
	int i;

	for(i=0;i<_3D;i++) {
		Rxn_v[i]=0.;
		Rxn_u[i]=0.;
		Rv[i]=0.;
		Ru[i]=0.;
		irn[i]=0.;
	}

	//Calculate centre wrt point of contact
	for(i=0;i<DIM;i++) Rv[i] = r[i]-Pv[i];
	for(i=0;i<DIM;i++) Ru[i] = r[i]-Pu[i];

	J = dotprod( V,n,DIM ) - dotprod( U,n,DIM );
	denom = invMv + invMu;

	crossprod( Rv,n,Rxn_v );
// 	dotprodmat( Rxn_v,Iv,irn,_3D );
	dotprodMatVec( Iv,Rxn_v,irn,_3D );

	J += dotprod( Rxn_v,Wv,_3D );
	denom += dotprod( Rxn_v,irn,DIM );

	crossprod( Ru,n,Rxn_u );
// 	dotprodmat( Rxn_u,Iu,irn,_3D );
	dotprodMatVec( Iu,Rxn_u,irn,_3D );

	J -= dotprod( Rxn_u,Wu,_3D );
	denom += dotprod( Rxn_u,irn,_3D );

	J *= -1.*(1. + E)/denom;
	return J;
}

///
/// @brief Applies boundary conditions to particle velocity.
///
/// Applies boundary conditions (specific to `WALL`) to particle velocity. TODO brief description how.
/// @param pp The individual particle.
/// @param WALL The boundary.
/// @param n The normal vector to the boundary
/// @param SP TODO
/// @param KBT TODO
/// @see MPC_BCcollision()
/// @note TODO give examples of this routine
///
void velBC( particleMPC *pp,bc *WALL,double n[_3D],spec *SP,double KBT ) {
/*
    This subroutine applies the BC transformation to velocity.
		The BC is assumed to be in a rest frame (slow moving compared to the MPCD particle)
*/
	double V[_3D],VN[_3D],VT[_3D],VR[_3D],R[_3D],zip[_3D];
	double IIpart[_3D][_3D],IIwall[_3D][_3D];
	double IMpart,IMwall;		//Inverse mass
	double J=1.;			//Impulse
	double rand;			//Random number
	int i,j;

	//Zero everything
	for( i=0; i<_3D; i++ ) {
		V[i] = 0.;
		VN[i] = 0.;
		VT[i] = 0.;
		VR[i] = 0.;
		R[i] = 0.;
		zip[i] = 0.;
	}

	//Calculate rotational velocity of BC at the contact point
	for( i=0; i<DIM; i++ ) R[i] = pp->Q[i] - WALL->Q[i];
	crossprod( WALL->L,R,VR );
	//Begin Determining the direction of the impulse
	for( i=0; i<DIM; i++ ) V[i] = pp->V[i] - WALL->V[i] - VR[i];
	// 	Calculate normal and tangential components of velocity
	proj( V,n,VN,DIM );
	tang( V,VN,VT,DIM );

	//Find the inverted moment of innertia and mass
	//Point particles do not have angular momentum
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) IIpart[i][j] = 0.0;
	IMpart = 1.0 / (SP + (pp->SPID))->MASS;
	//BC
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) IIwall[i][j] = 0.0;
	IMwall = 0.0;
	if ( WALL->DSPLC ) {
		invert3x3(IIwall,WALL->I);
		IMwall = 1.0 / WALL->MASS;
	}
	else if ( WALL->DSPLC > 1 ) {
		printf( "Error: WALL.DSPLC must be 0 or 1.\n" );
		exit( 1 );
	}
	//The energy/momentum/angular momentum conserving impulse method
	if( WALL->COLL_TYPE == BC_IMP ) {
		vel_trans( WALL,VN,VT,n );
		//Combine normal and tangential components
		for( i=0; i<DIM; i++ ) V[i] = VN[i] + VT[i];
		//Add the cartesion shift
		for( i=0; i<DIM; i++ ) V[i] += WALL->DVxyz[i];
		//The impulse's direction is the difference between the two (final minus initial)
		for( i=0; i<DIM; i++ ) V[i] -= ( pp->V[i] - WALL->V[i] - VR[i] );
		//Normalize V. We only want the unit vector. Conservation will determine its magnitude
		norm( V,DIM );

		J = impulse( V,pp->V,WALL->V,pp->Q,WALL->Q,zip,WALL->L,IMpart,IMwall,IIpart,IIwall,pp->Q,WALL->E );
	}
	//The rule method such as bounceback or reflection or periodic which does NOT necesarily conserve momentum
	else if( WALL->COLL_TYPE == BC_SURF_RULES ) {

		vel_trans( WALL,VN,VT,n );
		//Combine normal and tangential components
		for( i=0; i<DIM; i++ ) V[i] = VN[i] + VT[i];
		//Add the cartesion shift
		for( i=0; i<DIM; i++ ) V[i] += WALL->DVxyz[i];
		//Move the velocity out of the particle's rest frame and back into the lab frame
		for( i=0; i<DIM; i++ ) V[i] += WALL->V[i] + VR[i];
		//Make V the change in momentum instead of the velocity by subtracting the old velocity and times mass
		for( i=0; i<DIM; i++ ) {
			V[i] -= pp->V[i];
			V[i] *= (SP + (pp->SPID))->MASS;
		}
		//Set J=1 because the total change in momentum (i.e. impulse) is in V.
		J = 1.0;
	}
	else if( WALL->COLL_TYPE == BC_THERMO_SURF ) {
		vel_trans( WALL,VN,VT,n );
		norm( VN,DIM );
		norm( VT,DIM );
		// Apply Thermo boundary conditions
		rand = genrand_gaussMB(WALL->KBT,(SP + (pp->SPID))->MASS);
		for( i=0; i<DIM; i++ ) VT[i] *= rand;
		rand = genrand_rayleigh( sqrt(WALL->KBT/(SP + (pp->SPID))->MASS) );
		for( i=0; i<DIM; i++ ) VN[i] *= rand;
		//Combine normal and tangential components
		for( i=0; i<DIM; i++ ) V[i] = VN[i] + VT[i];
		//Add the cartesion shift
		for( i=0; i<DIM; i++ ) V[i] += WALL->DVxyz[i];

		/************************************************************************/
		/************************************************************************/
		//Needs a thermostat to be on!
		/************************************************************************/
		/************************************************************************/
		//Move the velocity out of the particle's rest frame and back into the lab frame
		for( i=0; i<DIM; i++ ) V[i] += WALL->V[i] + VR[i];
		//Don't because BC not moving during collision - it is paused

		//Make V the change in momentum instead of the velocity by subtracting the old velocity and times mass
		for( i=0; i<DIM; i++ ) {
			V[i] -= pp->V[i];
			V[i] *= (SP + (pp->SPID))->MASS;
		}
		//Set J=1 because the total change in momentum (i.e. impulse) is in V.
		J = 1.0;
	}
	else if( WALL->COLL_TYPE == BC_HALF_RULES ) {
		printf("Error: Rule-based solvent/BC collisions at t/2 not yet operational.\n");
		exit(1);
	}
	else if( WALL->COLL_TYPE == BC_THERMO_HALF ){
		printf("Error: Probabilistic solvent/BC collisions at t/2 not yet operational.\n");
		exit(1);
	}
	else {
		printf("Error: Solvent/BC Collision type unrecognized.\n");
		exit(1);
	}

	//Use the impulse to set the velocity of particleMPC
	for( i=0; i<DIM; i++ ) pp->V[i] += V[i] * J * IMpart;
	// #ifdef DBG
	// 	if( DBUG == DBGMPCBC ) {
	// 		printf( "\tParticle dV: [%lf,%lf,%lf]\n",V[0]*J*IMpart,V[1]*J*IMpart,V[2]*J*IMpart );
	// 	}
	// #endif

	if( WALL->DSPLC ) {
		//Set the velocity of BC
		for( i=0; i<DIM; i++ ) WALL->dV[i] -= V[i] * J * IMwall;
		#ifdef DBG
			if( DBUG == DBGMPCBC ) {
				printf( "\tWall dV: [%lf,%lf,%lf]\n",-V[0]*J*IMwall,V[1]*J*IMwall,V[2]*J*IMwall );
			}
		#endif
		//Set the angular velocity of BC
		if(DIM>1) {
			//Recall R is vector from contact point to centre of mass
			//Use VT as crossprod of separation with direction of imp
			//NOTE: Zero for reflective
			crossprod( R,V,VT );
			//Use VN as dotprod of mom inertia tens with VT
			dotprodMatVec( IIwall,VT,VN,_3D );
			for( i=0; i<_3D; i++) WALL->dL[i] -= VN[i] * J;
			#ifdef DBG
				if( DBUG == DBGMPCBC ) {
					printf( "\tWall dL: [%lf,%lf,%lf]\n",-VN[0]*J,VN[1]*J,VN[2]*J );
				}
			#endif
		}
		else for( i=0; i<_3D; i++) WALL->dL[i] = 0.0; // fallback for 1D
	}
}
///
/// @brief Applies boundary conditions to particle position.
///
/// Applies boundary conditions (specific to `WALL`) to particle position.
/// @param pp The individual particle.
/// @param WALL The boundary.
/// @param n Normal vector to the boundary.
/// @see MPC_BCcollision()
/// @note This routine is called for the application of periodic boundary conditions.
///
void posBC( particleMPC *pp,bc WALL,double n[] ) {
/*
    This subroutine applies the BC transformation to position.
    IT DOES NOT STREAM! THAT IS DONE IN SEPARATE ROUTINE.
*/
	double PN[DIM],PT[DIM];
	int i;

	proj( pp->Q,n,PN,DIM );		//Calculate normal component of the position
	tang( pp->Q,PN,PT,DIM );		//Calculate tangential component of position

	//Transform the position
	for( i=0; i<DIM; i++ ) {
		PN[i] += WALL.DN*n[i]*( 1.-2.*TOL );	//Transform normal component. Multiply by 1-2TOL to ensure inside box incase of numerical error
		PT[i] += WALL.DT;				//Transform tangential component
	}
	//Combine normal and tangential components
	for( i=0; i<DIM; i++ ) pp->Q[i] = PN[i] + PT[i];
}
///
/// @brief Applies boundary conditions to the position of another boundary.
///
/// In the case of a collision between two boundaries, this routine applies
/// the boundary condition transformation (specific to `WALL1`) to a second boundary's (`WALL2`) position.
/// @param WALL1 The first boundary.
/// @param WALL2 The second boundary (that has position updated from interaction with `WALL1`).
/// @param n Normal vector to the boundary (of `WALL1`).
/// @see BC_BCcollision()
///
void BCBCpos( bc *WALL1 ,bc *WALL2,double n[] ) {
	double PN[DIM],PT[DIM];
	int i;

	proj( WALL1->Q,n,PN,DIM );		//Calculate normal component of the position
	tang( WALL1->Q,PN,PT,DIM );		//Calculate tangential component of position

	//Transform the position
	for( i=0; i<DIM; i++ ) {
		PN[i] += WALL2->DN * n[i];	//Transform normal component
		PT[i] += WALL2->DT;		//Transform tangential component
	}
	//Combine normal and tangential components
	for( i=0; i<DIM; i++ ) WALL1->Q[i] = PN[i] + PT[i];
}
// void BCBCvel( double cp[_3D],bc *W1,bc *W2,double n[_3D],double KBT ) {
// /*
//     This subroutine applies the BC transformation of W2 to velocity of W1.
// */
// 	double V[_3D],VN[_3D],VT[_3D],VR1[_3D],VR2[_3D],R[_3D];
// 	double II1[_3D][_3D],II2[_3D][_3D];	//Inverse moment of inertia
// 	double IM1,IM2;			//Inverse mass
// 	double J=1.;			//Impulse
// 	double rand;			//Random number
// 	int i,j;
//
// 	//Zero everything
// 	for( i=0; i<_3D; i++ ) {
// 		V[i] = 0.;
// 		VN[i] = 0.;
// 		VT[i] = 0.;
// 		VR1[i] = 0.;
// 		VR2[i] = 0.;
// 		R[i] = 0.;
// 	}
//
// 	//Rotational velocities
// 	for( i=0; i<DIM; i++ ) R[i] = cp[i] - W1->Q[i];
// 	crossprod( W1->L,R,VR1 );
// 	for( i=0; i<DIM; i++ ) R[i] = cp[i] - W2->Q[i];
// 	crossprod( W2->L,R,VR2 );
//
// 	//Calculated the total velocity at the contact point
// 	for( i=0; i<DIM; i++ ) V[i] = W1->V[i] - W2->V[i] + VR1[i] - VR2[i];
//
// 	//Calculate normal and tangential components of velocity
// 	proj( V,n,VN,DIM );
// 	tang( V,VN,VT,DIM );
//
// 	//Calculate inverse mass and moments of inertia
// 	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) {
// 		II1[i][j] = 0.;
// 		II2[i][j] = 0.;
// 	}
// 	IM1 = 0.;
// 	IM2 = 0.;
// 	if ( W1->DSPLC ) {
// 		invert3x3(II1,W1->I);
// 		IM1 = 1. / W1->M;
// 	}
// 	if ( W2->DSPLC ) {
// 		invert3x3(II2,W2->I);
// 		IM2 = 1. / W2->M;
// 	}
//
// 	//The energy/momentum/angular momentum conserving impulse method
// 	if( W2->COLL_TYPE == BC_IMP ) {
// 		vel_trans( W2,VN,VT,n );
// 		//Combine normal and tangential components
// 		for( i=0; i<DIM; i++ ) V[i] = VN[i] + VT[i];
// 		//The impulse's direction is the difference between the two (final minus initial)
// 		for( i=0; i<DIM; i++ ) V[i] -= (W1->V[i] - W2->V[i] + VR1[i] - VR2[i]);
// 		//Normalize V. We only want the unit vector. Conservation will determine its magnitude
// 		norm( V,DIM );
//
// 		J = impulse( V,W1->V,W2->V,W1->Q,W2->Q,W1->L,W2->L,IM1,IM2,II1,II2,cp,W2->E );
// 	}
// 	//The rule method such as bounceback or reflection or periodic which does NOT necesarily conserve momentum
// 	else if( W2->COLL_TYPE == BC_SURF_RULES ) {
// 		vel_trans( W2,VN,VT,n );
// 		//Combine normal and tangential components
// 		for( i=0; i<DIM; i++ ) V[i] = VN[i] + VT[i];
// 		//Move the velocity out of the particle's rest frame and back into the lab frame
// 		for( i=0; i<DIM; i++ ) V[i] += W2->V[i] + VR2[i];
// 		//Make V the change in momentum instead of the velocity by subtracting the old velocity and times mass
// 		for( i=0; i<DIM; i++ ) {
// 			V[i] -= W1->V[i];
// 			V[i] *= W1->M;
// 		}
// 		//Set J=1 because the total change in momentum (i.e. impulse) is in V.
// 		J = 1.;
// 	}
// 	else if( W2->COLL_TYPE == BC_THERMO_SURF ) {
// 		vel_trans( W2,VN,VT,n );
// 		norm( VN,DIM );
// 		norm( VT,DIM );
// 		// Apply Thermo boundary conditions
// 		rand = genrand_gauss();
// 		for( i=0; i<DIM; i++ ) VT[i] *= rand;
// 		rand = genrand_rayleigh( sqrt(KBT/W1->M ) );
// 		for( i=0; i<DIM; i++ ) VN[i] *= rand;
// 		//Combine normal and tangential components
// 		for( i=0; i<DIM; i++ ) V[i] = VN[i] + VT[i];
//
// 		//Move the velocity out of the particle's rest frame and back into the lab frame
// 		for( i=0; i<DIM; i++ ) V[i] += W2->V[i] + VR2[i];
// 		//Make V the change in momentum instead of the velocity by subtracting the old velocity and times mass
// 		for( i=0; i<DIM; i++ ) {
// 			V[i] -= W1->V[i];
// 			V[i] *= W1->M;
// 		}
// 		//Set J=1 because the total change in momentum (i.e. impulse) is in V.
// 		J = 1.;
// 	}
// 	else if( W2->COLL_TYPE == BC_HALF_RULES ) {
// 		printf("Error: Rule-based solvent/BC collisions at t/2 not yet operational.\n");
// 		exit(1);
// 	}
// 	else if( W2->COLL_TYPE == BC_THERMO_HALF ){
// 		printf("Error: Probabilistic solvent/BC collisions at t/2 not yet operational.\n");
// 		exit(1);
// 	}
// 	else {
// 		printf("Error: Solvent/BC Collision type unrecognized.\n");
// 		exit(1);
// 	}
//
// 	//Use the impulse to set the velocity of W1
// 	for( i=0; i<DIM; i++ ) W1->V[i] += V[i] * J * IM1;
// 	for( i=0; i<DIM; i++ ) VN[i] = cp[i] - W1->Q[i];
// 	crossprod( VN,V,VT );
// // 	dotprodmat( VT,II1,VN,_3D );
// 	dotprodMatVec( II1,VT,VN,_3D );
// 	for( i=0; i<_3D; i++) W1->L[i] += VN[i] * J;
//
// 	//Use the impulse to set the velocity of W2
// 	if( W2->DSPLC ) {
// 		for( i=0; i<DIM; i++ ) W2->V[i] -= V[i] * J * IM2;
// 		for( i=0; i<DIM; i++ ) VN[i] = cp[i] - W2->Q[i];
// 		crossprod( VN,V,VT );
// // 		dotprodmat( VT,II2,VN,_3D );
// 		dotprodMatVec( II2,VT,VN,_3D );
// 		for( i=0; i<_3D; i++) W2->L[i] -= VN[i] * J;
// 	}
// }

///
/// @brief Transforms the particle's normal and tangential velocity components.
///
/// Applies boundary conditions to the particle's normal velocity component and tangential component.
/// These normal and tangential components were identified in velBC().
/// The applied boundary conditions include transformations to the magnitude of the velocity,
/// and flips to the orientation of the velocity at the boundary.
/// @param WALL The boundary.
/// @param VN Normal component of the velocity.
/// @param VT Tangential component of the velocity.
/// @param norm Normal vector (e.g. to the boundary).
/// @see velBC()
///
void vel_trans( bc *WALL,double VN[],double VT[],double norm[] ) {

	int i;
	for( i=0; i<DIM; i++ ) {
		VN[i] *= WALL->MVN;
		// Adding a velocity component (magnitude is DVN) in the normal direction to the boundary
		VN[i] += WALL->DVN*norm[i];
		VT[i] *= WALL->MVT;
		// Since VT has no component normal to the boundary, we add (magnitude of DVT) to all remaining components.
		VT[i] += WALL->DVT;
	}
}
double chooseT( double tstep,double tp,double tn,int p,int *flag ) {
/*
    This routine choose which of the
    calculated cross times is the proper
    one for the collision event
*/
	double tc=tstep;		//chosen time
	double zero = -TOL;
	double step = tstep+TOL;

	if( tp>step && tn>step ) {
		*flag = 1;
		tc=0.;
	}
	else if( tp<zero && tn<zero ){
		*flag = 1;
		tc=0.;
	}
	else if( tp<zero ) tc = tn;
	else if( tn<zero ) tc = tp;
	else if( tn<tp ) tc = tn;
	else tc = tp;

	//Check the chosen time
	if( tc<zero ) {
		*flag = 1;
		tc=0.;
	}
	if( tc>step ) {
		*flag = 1;
		tc=0.;
	}
	return tc;
}
void chooseP( bc WALL,particleMPC *pp,double *chosenW,int *chosenP ) {
/*
   Check if any of the particles are now inside the BC.
   Find the particle
 */
 	int i;
	double tempW = 1.0;
	double shift[DIM];
	*chosenW = 1.0;

	for( i=0; i<GPOP; i++ ) {
		//Shift the BC due to any periodic BCs
		shiftBC( shift,&WALL,(pp+i) );
		rotateBC( &WALL,(pp+i),0 );
		tempW = calcW( WALL,*(pp+i) );
		if( tempW < -TOL ) {
			//Particle within BC
			*chosenW = tempW;
			*chosenP = i;
			break;
		}
		//Shift the BC back to it's real position
		rotatebackBC( &WALL,(pp+i),0 );
		shiftbackBC( shift,&WALL );
	}
}
void BC_BCcollision( bc *movingWall,bc *stillWall,double t_step,int *flag ) {
/*
   Check if movingWall violates stillWall (still wall is not necesarily still - but the moving wall is moving at this part of the code).
   If so apply collision
   CURRENTLY ONLY DEALS WITH SPHERES
*/
	int i;
	double W;
	double n[_3D];		//Normal
	double V[_3D],VN[_3D],VT[_3D],VR[_3D],R[_3D];

	//Zero everything
	for( i=0; i<_3D; i++ ) {
		V[i] = 0.;
		VN[i] = 0.;
		VT[i] = 0.;
		VR[i] = 0.;
		R[i] = 0.;
	}

	//See if movingWall violates stillWall
	W = calcW_BC( *movingWall,*stillWall,0 );

	if( W<=0.0 ) {
		#ifdef DBG
			if( DBUG == DBGBCBC ) {
				printf( "\tBC-BC collision identified.\n" );
				printf( "\tW=%e\n",W );
				bccoord( *movingWall );
			}
		#endif

		/* ****************************************** */
		/* ******** VELOCITY TRANSFORMATION ********* */
		/* ****************************************** */
		// Calculate the normal
		// NOTE this probably won't work if p!=1,2
		normal( n,*stillWall,movingWall->Q,DIM );
		norm( n,DIM );
		//Calculate rotational velocity of BC
		for( i=0; i<DIM; i++ ) R[i] = movingWall->Q[i] - stillWall->Q[i];
		crossprod( stillWall->L,R,VR );
		//Determining the relative velocities of the collision points on each surface
		for( i=0; i<DIM; i++ ) V[i] = movingWall->V[i] - stillWall->V[i] - VR[i];
		#ifdef DBG
			if( DBUG == DBGBCBC ) {
				printf( "\tRelative collision velocity:" );
				pvec( V,DIM );
			}
		#endif
		//Calculate normal and tangential components of velocity
		proj( V,n,VN,DIM );
		tang( V,VN,VT,DIM );
		#ifdef DBG
			if( DBUG == DBGBCBC ) printf( "\tCalculate tangential vel.\n" );
		#endif
		//No-slip - particle must bounce back or reflect depending on stillWall
		for( i=0; i<DIM; i++ ) {
			VT[i] *= stillWall->MVT;
			VT[i] += stillWall->DVT;
		}
		//Bounce off wall
		#ifdef DBG
			if( DBUG == DBGBCBC ) printf( "\tCalculate normal vel.\n" );
		#endif
		//If the mobile boundary is passing a PBC then it's velocity should not change.
		if( feq(stillWall->MVN,1.0) && feq(stillWall->MVT,1.0) ) {
			#ifdef DBG
				if( DBUG == DBGBCBC ) printf( "\t\tMoving BC passing PBC.\n" );
			#endif
			//Just let them move it not collision-type interaction
			for( i=0; i<DIM; i++ ) VN[i] += stillWall->DVN;
		}
		//Otherwise the mobile boundary is colliding with the still wall
		else{
			#ifdef DBG
				if( DBUG == DBGBCBC ) printf( "\t\tMoving BC will rebound.\n" );
			#endif
			//Only apply the velocity transoformation if the BCs are moving towards each other --- the should be (but just in case)
			if( dotprod( n,VN,_3D )<0. ) {
				#ifdef DBG
					if( DBUG == DBGBCBC ) printf( "\t\tMoving BC approaching still BC ---> transform normal vel.\n" );
				#endif
				//The particle is moving toward the wall
				for( i=0; i<DIM; i++ ) {
					VN[i] *= stillWall->MVN;
					VN[i] += stillWall->DVN;
				}
			}
			//If they're moving apart then just let them
			else {
				#ifdef DBG
					if( DBUG == DBGBCBC ) printf( "\t\tMoving BC already receding from still BC ---> do nothing.\n" );
				#endif
				//The particle is moving away from the wall
			}
			//Flag that it did move
			*flag+=1;
		}
		#ifdef DBG
			if( DBUG == DBGBCBC ) {
				printf( "\t\tTransformed vel:\n" );
				printf( "\t\t\tNormal: " );
				pvec( VN,DIM );
				printf( "\t\t\tTangential: " );
				pvec( VT,DIM );
			}
		#endif

		#ifdef DBG
			if( DBUG == DBGBCBC ) printf( "\tAssign new vel.\n" );
		#endif
		//Combine normal and tangential components
		for( i=0; i<DIM; i++ ) V[i] = VN[i] + VT[i];
		//Add the cartesion shift
		for( i=0; i<DIM; i++ ) V[i] += stillWall->DVxyz[i];
		//Set velocity
		for( i=0; i<DIM; i++ ) movingWall->V[i] = V[i] + stillWall->V[i] + VR[i];
		//Elastically change the velocity of the "still wall"
		if( stillWall->DSPLC==1 ) {
			for( i=0; i<DIM; i++ ) stillWall->V[i] = -1.0*movingWall->V[i] * (float)(movingWall->MASS) / (float)(stillWall->MASS);
		}
		#ifdef DBG
			if( DBUG == DBGBCBC ) {
				printf( "\tNew vel: " );
				pvec( movingWall->V,DIM );
			}
		#endif

		// #ifdef DBG
		// 	if( DBUG == DBGBCBC ) printf( "\tAssign new ang vel.\n" );
		// #endif
		// //Set angular velocity
		// norm( R,DIM );
		// // crossprod( VT,R,movingWall->L );
		// // for( i=0; i<_3D; i++ ) movingWall->dL[i] += (VT[i] + VR[i]) * smrtPow(movingWall->R,1./movingWall->P);
		// // printf( "BCBC ang vel: " );
		// // pvec( movingWall->dL,_3D );
		// #ifdef DBG
		// 	if( DBUG == DBGBCBC ) {
		// 		printf( "\tCurrent ang vel: " );
		// 		pvec( movingWall->dL,_3D );
		// 	}
		// #endif

		/* ****************************************** */
		/* ******** POSITION TRANSFORMATION ********* */
		/* ****************************************** */
		#ifdef DBG
			if( DBUG == DBGBCBC ) printf( "\tAssign new pos.\n" );
		#endif
		// //Put particle back where it was if it's a bounce-back wall
		// if( feq(stillWall->MVN,-1.0) && feq(stillWall->MVT,-1.0) ) {
		// 	#ifdef DBG
		// 		if( DBUG == DBGBCBC ) printf( "\tBounce-back so replace\n" );
		// 	#endif
		// 	for( i=0; i<DIM; i++ ) movingWall->Q[i] = movingWall->Q_old[i];
		// }
		// else{
		// 	//Only apply the positional BCs if the CENTRE of the particle has passed the other BC
		// 	W = calcW_BC( *movingWall,*stillWall,1 );
		// 	if( W<=0.0 ) BCBCpos( movingWall,stillWall,n );
		// }
		// PBC
		if( feq(stillWall->MVN,1.0) && feq(stillWall->MVT,1.0) ) {
			//Only apply the positional BCs if the CENTRE of the particle has passed the other BC
			W = calcW_BC( *movingWall,*stillWall,1 );
			if( W<=0.0 ) BCBCpos( movingWall,stillWall,n );
		}
		else {
			//Put particle back where it was if it's not a PBC (all BC-BC are bounce-back)
			#ifdef DBG
				if( DBUG == DBGBCBC ) printf( "\tBounce-back so replace\n" );
			#endif
			for( i=0; i<DIM; i++ ) movingWall->Q[i] = movingWall->Q_old[i];
		}
		#ifdef DBG
			if( DBUG == DBGBCBC ) {
				printf( "\tNew pos: " );
				pvec( movingWall->Q,DIM );
			}
		#endif
		// #ifdef DBG
		// 	if( DBUG == DBGBCBC ) wait4u();
		// #endif
	}
}
void BC_MPCcollision(bc WALL[],int BCcurrent,particleMPC *pp,spec *pSP,double KBT,double GRAV[],double t_step,simptr simMD,int MDmode,int LC,int *bcCNT,int *reCNT,int *rethermCNT) {
/*
   This subroutine takes care of all the collision occuring
   because a BC is translating through SRD particles.
	 It does so by giving any collided particles negative the BC's velocity, rewinding them and playing them from there.
   It shouldn't have to worry about collisions with other BCs which were handled in BC_BCcollision()
*/
	int i;
	int chosenP=GPOP+1;					//Particle to go with t_min
	int flag = 1;								//flag for if should keep looping. Loop while 1.
	double pV[DIM],bcV[DIM],bcL[_3D];		//temporary velocities
	double W;

	#ifdef DBG
		if( DBUG == DBGBCMPC ) {
			printf( "BC %d\n",BCcurrent );
			bccoord( WALL[BCcurrent] );
		}
	#endif
	while( flag ) {
		/* ****************************************** */
		/* ************ BCs ON PARTICLES ************ */
		/* ****************************************** */
		//We must check if any of the particles are now inside the BC,
		//Find the one that the collision should have occured with and the time for that collision
		#ifdef DBG
			if( DBUG == DBGBCMPC ) printf( "Choose a particle: " );
		#endif
		// chooseP( WALL[BCcurrent],pp,&time,&W,&chosenP,t_step,GRAV );
		chooseP( WALL[BCcurrent],pp,&W,&chosenP );
		#ifdef DBG
			if( DBUG == DBGBCMPC ) printf( "BC=%d; particle=%d; W=%e\n",BCcurrent,chosenP,W );
		#endif
		//If the particle doesn't collide with any MPCD particles exit
		if( chosenP>=GPOP ) flag = 0;
		//If no particles were inside then we are done.
		if( W > -TOL ) flag = 0;
		else{
			//A particle got inside a BC since the BC translated
			#ifdef DBG
				if( DBUG == DBGBCMPC ) {
					printf( "%d\n",chosenP );
					printf( " W: %e\n", W );
					pvec( (pp+chosenP)->Q,DIM );
					pvec( (pp+chosenP)->V,DIM );
					printf( " distance: %lf\n", distpoints( WALL[BCcurrent].Q,(pp+chosenP)->Q,DIM ) );
				}
			#endif
			//Save the particle and BC's velocities
			for( i=0; i<DIM; i++ ) pV[i]=(pp+chosenP)->V[i];
			for( i=0; i<DIM; i++ ) bcV[i]=WALL[BCcurrent].V[i];
			for( i=0; i<_3D; i++ ) bcL[i]=WALL[BCcurrent].L[i];
			//Give the particle negative the BC's velocity
			for( i=0; i<DIM; i++ ) (pp+chosenP)->V[i]=-bcV[i];
			//The change in reference frame means BC vel=0 (at the surface)
			for( i=0; i<DIM; i++ ) WALL[BCcurrent].V[i]=-pV[i];
			for( i=0; i<_3D; i++ ) WALL[BCcurrent].L[i]=0.0;
			#ifdef DBG
				if( DBUG == DBGBCMPC ) {
					printf( "Change the reference frame:\nBC %d\n",BCcurrent );
					pvec( WALL[BCcurrent].Q,DIM );
					pvec( WALL[BCcurrent].V,DIM );
					printf( "MPC Particle %d\n",chosenP );
					pvec( (pp+chosenP)->Q,DIM );
					pvec( (pp+chosenP)->V,DIM );
					printf( "Full rewind:" );
					rewind_BC( &WALL[BCcurrent],t_step );
					printf( " distance: %lf\n", distpoints( WALL[BCcurrent].Q,(pp+chosenP)->Q,DIM ) );
					stream_BC( &WALL[BCcurrent],t_step );
				}
			#endif
			//With the BC's negative velocity work out all the BCs for this particle
			MPC_BCcollision( pp,chosenP,WALL,pSP,KBT,t_step,LC,bcCNT,reCNT,rethermCNT,0 );
			#ifdef DBG
				if( DBUG == DBGBCMPC ) {
					printf( "New particle position:\n" );
					pvec( (pp+chosenP)->Q,DIM );
					pvec( (pp+chosenP)->V,DIM );
					printf( " distance: %lf\n", distpoints( WALL[BCcurrent].Q,(pp+chosenP)->Q,DIM ) );
				}
			#endif
			//Give the particle back it's old velocity (plus the velocity it got from these collisions)
			#ifdef DBG
				if( DBUG == DBGBCMPC ) 	printf( "Change the reference frame back\n" );
			#endif
			for( i=0; i<DIM; i++ ) (pp+chosenP)->V[i] += pV[i] + bcV[i];
			for( i=0; i<DIM; i++ ) WALL[BCcurrent].V[i] += pV[i] + bcV[i];
			for( i=0; i<_3D; i++ ) WALL[BCcurrent].L[i] = bcL[i];
			//Done. Reset particle number test
			chosenP=GPOP+1;
			// wait4u();
		}
	}
}
void chooseBC( bc WALL[],int currentP,particleMPC *pp,double *t_minColl,double *chosenW,int *chosenBC,double t_left ) {
/*
   We must check if the particle is inside any of the BCs
   This subroutine finds the BC and the time of collision
*/

	int i,flag;
	double t1,t2,tc;
	double tempW, shift[DIM];

	*t_minColl = t_left;
	*chosenW=1.0;
	tempW = 1.0;
	flag=0;

	for( i=0; i<NBC; i++ ) {
		// Planar BCs
		if( WALL[i].PLANAR ) {
			tempW=calcW_PLANE( WALL[i],*(pp+currentP) );
			if( tempW < 0.0 ) {
				//printf("ChoosePlanar\n");
				//Calculate crosstime
				crosstimeReverse( *(pp+currentP),WALL[i],&t1,&t2,t_left );
				tc=t_left-t1;
				if( tc < *t_minColl ) {
					*t_minColl = tc;
					*chosenBC = i;
					*chosenW = tempW;
				}
			}
		}
		//Planar wavy BCs
		else if ( !feq(WALL[i].B[0],0.0) && feq(WALL[i].P[0],1.0) && feq(WALL[i].P[1],1.0) && feq(WALL[i].P[2],1.0) ) {
			tempW=calcW( WALL[i],*(pp+currentP) );
			if( tempW < 0.0 ) {
				//Calculate crosstime
				crosstimeReverse( *(pp+currentP),WALL[i],&t1,&t2,t_left );
				t1=t_left-t1;
				t2=t_left-t2;
				tc = chooseT( t_left,t1,t2,currentP,&flag );
				if( flag ) {
					printf( "Error: Cross time unacceptable: %lf.\n",tc );
				}
				if( tc < *t_minColl ) {
					*t_minColl = tc;
					*chosenBC = i;
					*chosenW = tempW;
				}
			}
		}
		//Non-planar BCs
		else {
			//Shift BC due to periodic BCs
			shiftBC( shift,&WALL[i],(pp+currentP) );
			rotateBC( &WALL[i],(pp+currentP),0 );
			tempW=calcW( WALL[i],*(pp+currentP) );
			if( tempW < 0.0 ) {
				//Calculate crosstime
				crosstimeReverse( *(pp+currentP),WALL[i],&t1,&t2,t_left );
				t1=t_left-t1;
				t2=t_left-t2;
				// #ifdef DBG
				// 	if( DBUG==DBGMPCBC || DBUG==DBGBCMPC) {
				// 		printf( "t1=%e, t2=%e\n",t1,t2 );
				// 	}
				// #endif
				tc = chooseT( t_left,t1,t2,currentP,&flag );
				if( flag ) {
					printf( "Error: Cross time unacceptable: %lf.\n",tc );
					//exit( 1 );
				}
				if( tc < *t_minColl ) {
					*t_minColl = tc;
					*chosenBC = i;
					*chosenW = tempW;
				}
			}
			//Shift BC back
			rotatebackBC( &WALL[i],(pp+currentP),0 );
			shiftbackBC( shift,&WALL[i] );
		}
	}
}
void MPC_BCcollision( particleMPC *pp,int currentP,bc WALL[],spec *pSP,double KBT,double t_step,int LC,int *bcCNT,int *reCNT,int *rethermCNT, int flagMPConBC ) {
	/*
	   This routine does the collision of MPCD particles with a BC
		 Normal use is flagMPConBC=1 but this routine is also used for BC on MPCD when flagMPConBC=0
	*/
	double t_left = t_step;	//time left to move for
	double t_coll=0.0;			//time to go from current position to the collision point
	int chosenBC=0;					//BC that the collision occurs with
	int flag = 1;						//flag for if should keep looping. Loop while 1.
	int i,cnt=0,cnt2=0,flagCNT=0;
	double n[_3D] = {0.0,0.0,0.0};	//Normal to the surface
	double W=1.0,Wn=1.0,Wo=1.0;
	double shift[DIM];

	while( flag ) {
		//We must check if the particle is inside any of the BCs
		// printf( "\nChoosing BC...\n" );
		chooseBC( WALL,currentP,pp,&t_coll,&W,&chosenBC,t_left );
		//If no particles were inside then we are done.
		if( W > -TOL ) flag = 0;
		//Otherwise, COLLISON
		else {
			// #ifdef DBG
			// 	if( DBUG==DBGMPCBC || DBUG==DBGBCMPC) {
			// 		printf( "Chosen BC=%d, time left=%lf cross time=%lf\n",chosenBC,t_left,t_coll );
			// 	}
			// #endif
			cnt++;
			//Shift the BC due to any periodic BCs
			shiftBC( shift,&WALL[chosenBC],(pp+currentP) );
			rotateBC( &WALL[chosenBC],(pp+currentP),LC );
			if( t_coll>-TOL ) {
				// printf( "Enter collision W=%lf t_coll=%lf...\n",W,t_coll );
				//We have the BC to collide with and the time at which it collided
				//Rewind the particle back to when it should have collided with the surface
				rewind_P( (pp+currentP),t_left-t_coll );
				//Now the BC and the particle are touching. Let them collide.
				//Find normal
				normal( n,WALL[chosenBC],(pp+currentP)->Q,DIM );
				//Normalize the normal vector
				norm( n,DIM );
				//Apply the BC to velocity
				velBC( (pp+currentP),&WALL[chosenBC],n,pSP,KBT );
				//Apply the BC to position
				posBC( (pp+currentP),WALL[chosenBC],n );
				//Apply the BC to orientation
				if( LC!=ISOF ) oriBC( (pp+currentP),pSP,(WALL+chosenBC),n );
				//Update the time to stream for
				t_left-=t_coll;
				Wo = calcW( WALL[chosenBC],*(pp+currentP) );
				if( flagMPConBC )	{
					// printf( "Stream particle...\n" );
					//Let the MPC particle stream the rest of the way
					//ONLY for MPC on BC. When used for BC on MPC this step must be skipped (already at surface)
					stream_P( (pp+currentP),t_left );
					#ifdef DBG
						if( DBUG==DBGBCMAX ) {
							printf( "\tParticle %d and BC %d\n",currentP,chosenBC );
							printf( "\tPost streaming particle coordinates:\n" );
							pcoord( *(pp+currentP) );
							printf( "\tPost streaming BC coordinates:\n" );
							bccoord( WALL[chosenBC] );
							printf( "\tMake sure streamed out: W_beforeStream=%e",Wo );
							Wn = calcW( WALL[chosenBC],*(pp+currentP) );
							printf( " W_new=%e\n",Wn );
							printf( "\tDistance: %lf\n", distpoints( WALL[chosenBC].Q,(pp+currentP)->Q,DIM ) );
							printf( "\tNumber of collision attempts %d\n", cnt );
							if( Wn<-TOL ) {
								printf("Fuck me too\n");
								wait4u();
							}
						}
						// if( DBUG==DBGMPCBC  || DBUG==DBGBCMPC ) printf( "\tNumber of collision attempts %d\n", cnt );
					#endif
				}
				else{
						//When BC on MPC don't stream. Leave on surface BUT do make sure obeys PBCs
						// Ultimate PBC conditions for particles AT surface of moving BC
						// printf( "Set on surface...\n" );

						// THIS SEEMS TO BE THE PROBLEM!!! PARTICLE SET ON SURFACE BUT WHEN ROTATE BACK, NOT PERFECT AND FAILS???
						//NOTICE: This doesn't seem to be setting it anywhere at all!!!

						for( i=0; i<DIM; i++ ) {
							if( (pp+currentP)->Q[i] < 0.0 ) (pp+currentP)->Q[i] += XYZ[i];
							else if( (pp+currentP)->Q[i] > XYZ[i] ) (pp+currentP)->Q[i] -= XYZ[i];
						}
				}
				//Make sure streamed out
				W = calcW( WALL[chosenBC],*(pp+currentP) );
				//If didn't then place on surface
				if( W<-TOL ) rewind_P( (pp+currentP),t_left );
			}
			else {
				// printf( "DONOT enter collision...\n" );
				rewind_P( (pp+currentP),(t_left-t_coll) );
			}
			//Shift the BC back to it's real position
			// printf( "Shift and rotate BC back...\n" );
			rotatebackBC( &WALL[chosenBC],(pp+currentP),LC );
			shiftbackBC( shift,&WALL[chosenBC] );

			if( cnt>NBOUNCE ) {
				flagCNT=1;
				if( cnt<2*NBOUNCE && cnt2==0 ) {
					*reCNT += 1;
					//Try to let the particle rewind out a bit extra
					if( t_left>TOL ) t_left*=1.1;
					//If effectively zero, try resetting
					else {
						t_left=0.9*t_step;
						t_coll=0.1*t_step;
					}
					#ifdef DBG
						if( DBUG==DBGMPCBC || DBUG==DBGBCMPC ) printf( "\tParticle %d and BC %d stuck. cnt=%d cnt2=%d t_left=%lf\n",currentP,chosenBC,cnt,cnt2,t_left );
					#endif
				}
				else {
					*rethermCNT += 1;
					//Try a new velocity
					#ifdef DBG
						if( DBUG==DBGMPCBC  || DBUG==DBGBCMPC ) printf( "\tRe-thermalize particle %d\n",currentP );
					#endif
					push( (pp+currentP)->V,KBT,pSP[(pp+currentP)->SPID].VDIST, pSP[(pp+currentP)->SPID].MASS,NULL );
					t_left=t_step;
					t_coll=0.0;
					cnt2++;
					cnt=0;
				}
			}
			//Return to the top and try to move again for the rest of the time.
		}
	}
	*bcCNT += flagCNT;
	#ifdef DBG
		if( DBUG==DBGMPCBC  || DBUG==DBGBCMPC ) printf( "\tNumber of collision attempts %d\n", cnt );
	#endif
}
double *normal( double *n,bc WALL,double *point,int dimension ) {
	//2D!?
/*
   Find the normal to the surface at this point (particleMPC is
   presently ON surface). We take the gradient of
   ( a(x-h) )^p + (b(y-k))^p + (c(z-l))^p - r =0
   since the gradient is equal to the normal. For powers of
   1 and 2 we take shortcuts and have programmed in the specific
   solution, for higher powers we use a more general solution
*/
	int i;

	if( feq(WALL.ROTSYMM[0],4.0) && feq(WALL.ROTSYMM[1],4.0) ) {
		// Planar
		if( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) ) for( i=0; i<dimension; i++ ) n[i] = WALL.A[i];
		// Ellipsoid
		else if ((DIM == 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0)) || (DIM > 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0))) for( i=0; i<dimension; i++ ){
			n[i] = 2.0*WALL.A[i]*WALL.A[i]*(point[i]-WALL.Q[i]);
		}
		// squircle 4 and 6
		else if ((DIM == 2 && feq(WALL.P[0],4.0) && feq(WALL.P[1],4.0)) || (DIM > 2 && feq(WALL.P[0],4.0) && feq(WALL.P[1],4.0) && feq(WALL.P[2],4.0))) for( i=0; i<dimension; i++ ){
			n[i] = 4.0*WALL.A[i]*WALL.A[i]*WALL.A[i]*WALL.A[i]*(point[i]-WALL.Q[i])*(point[i]-WALL.Q[i])*(point[i]-WALL.Q[i]);
		}
		else if ((DIM == 2 && feq(WALL.P[0],6.0) && feq(WALL.P[1],6.0)) || (DIM > 2 && feq(WALL.P[0],6.0) && feq(WALL.P[1],6.0) && feq(WALL.P[2],6.0))) for( i=0; i<dimension; i++ ){
			n[i] = 4.0*WALL.A[i]*WALL.A[i]*WALL.A[i]*WALL.A[i]*WALL.A[i]*WALL.A[i];
			n[i] *= (point[i]-WALL.Q[i])*(point[i]-WALL.Q[i])*(point[i]-WALL.Q[i])*(point[i]-WALL.Q[i])*(point[i]-WALL.Q[i]);
		}
		// All others
		else {
			if( WALL.ABS ) for( i=0; i<dimension; i++ ) {
				n[i] = WALL.A[i]*WALL.A[i]*WALL.P[i]*(point[i]-WALL.Q[i]) * smrtPow( fabs(WALL.A[i]*(point[i]-WALL.Q[i])) , WALL.P[i]-2.0 );
			}
			else for( i=0; i<dimension; i++ ) {
				n[i] = WALL.P[i]*smrtPow( WALL.A[i],WALL.P[i] )*smrtPow( point[i]-WALL.Q[i],WALL.P[i]-1.0 );
			}
		}

		// adds correction for wavy boundaries
		if ( !feq(WALL.B[0],0.0) ) {
			normalWavy(n,WALL,point,dimension);
		}

	}
	else normalNon4foldSymm( n,WALL,point,dimension );

	return n;
}

double *normalWavy( double *n,bc WALL,double *point,int dimension ) {
	double W1=0.0,W2=0.0, div;
	double dw1[3], dw2[3];
	int i;
	for( i=0; i<3; i++ ) {
		dw1[i]=0.0;
		dw2[i]=0.0;
	}
	// Planar
	if( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) ) {
		if ( !feq(WALL.B[1],0.0) ) {
			W1 = (WALL.A[1]*point[0]-WALL.A[0]*point[1]) / sqrt( WALL.A[0]*WALL.A[0] + WALL.A[1]*WALL.A[1] );
			dw1[0] = WALL.A[1]/sqrt(WALL.A[0]*WALL.A[0]+WALL.A[1]*WALL.A[1]);
			dw1[1] = -WALL.A[0]/sqrt(WALL.A[0]*WALL.A[0]+WALL.A[1]*WALL.A[1]);
			dw1[2] = 0;
		}
		if( dimension>2 && !feq(WALL.B[2],0.0) ){
			div = sqrt( pow(WALL.A[0],4) + pow(WALL.A[1],4) + pow(WALL.A[0],2)*pow(WALL.A[2],2) + 2.0*pow(WALL.A[0],2)*pow(WALL.A[1],2) + pow(WALL.A[1],2)*pow(WALL.A[2],2) );
			W2 = ( WALL.A[0]*WALL.A[2]*point[0] + WALL.A[1]*WALL.A[2]*point[1] - (WALL.A[0]*WALL.A[0]+WALL.A[1]*WALL.A[1])*point[2] ) /div;
			dw2[0] = (WALL.A[0]*WALL.A[2]) / div;
			dw2[1] = (WALL.A[1]*WALL.A[2]) / div;
			dw2[2] = -(WALL.A[0]*WALL.A[0] + WALL.A[1]*WALL.A[1]) / div;
		}

	}
	else {
		//Ellipsoidal
		int flag = 0;
		for( i=0; i<dimension; i++ ) if( !feq(WALL.P[i],2.0) || feq(WALL.A[i],0.0) ) flag+=1;
		if(!flag){
			div = pow(WALL.A[0]*(point[0]-WALL.Q[0]),2) + pow(WALL.A[1]*(point[1]-WALL.Q[1]),2);
			W1 = atan2(WALL.A[1]*(-WALL.Q[1]+point[1]),WALL.A[0]*(-WALL.Q[0]+point[0]));
			if ( !feq(WALL.B[1],0.0) ) {
				dw1[0] = WALL.A[0]*WALL.A[1]*(WALL.Q[1]-point[1])/div;
				dw1[1] = WALL.A[0]*WALL.A[1]*(point[0]-WALL.Q[0])/div;
				dw1[2] = 0;
			}
			if( dimension>2 && !feq(WALL.B[2],0.0) ){
				dw2[0] = WALL.A[2]*WALL.A[0]*(point[0]-WALL.Q[0])*(point[2]-WALL.Q[2])/(sqrt(div)*(div+pow(WALL.A[2]*(point[2]-WALL.Q[2]),2)));
				dw2[1] = WALL.A[2]*WALL.A[1]*(point[1]-WALL.Q[1])*(point[2]-WALL.Q[2])/(sqrt(div)*(div+pow(WALL.A[2]*(point[2]-WALL.Q[2]),2)));
				dw2[2] = -WALL.A[2]*sqrt(div)/( div+pow(WALL.A[2]*(point[2]-WALL.Q[2]),2) );
			}
		}
		else if(dimension>2 && flag==1){
			//Cylinders
			flag = -1;
			for( i=0; i<DIM; i++) if ( feq(WALL.A[i],0.0) ) flag=i;
			if (flag!=-1) {
				int j, k, l;
				if (flag==0) j = 0, k = 1, l = 2;
				else if (flag==1) j = 1, k = 2, l = 0;
				else j = 2, k = 0, l = 1;

				div = pow(WALL.A[k]*(point[k]-WALL.Q[k]),2) + pow(WALL.A[l]*(point[l]-WALL.Q[l]),2);

				if ( !feq(WALL.B[1],0.0) ) {
					dw1[k] = WALL.A[k]*WALL.A[l]*(WALL.Q[l]-point[l])/div;
					dw1[l] = WALL.A[k]*WALL.A[l]*(point[k]-WALL.Q[k])/div;
					dw1[j] = 0;
				}
				if ( !feq(WALL.B[2],0.0) ) {
					dw2[k] = WALL.A[k]*WALL.A[k]*(point[k]-WALL.Q[k])*point[j]/( sqrt(div)*(div+point[k]) );
					dw2[l] = WALL.A[l]*WALL.A[l]*(point[l]-WALL.Q[l])*point[j]/( sqrt(div)*(div+point[l]) );;
					dw2[j] = -sqrt(div)/(div+point[j]);
				}
			}
		}
	}
	for( i=0; i<3; i++ ) n[i] -= WALL.B[0]*( WALL.B[1]*dw1[i]*sin(WALL.B[1]*W1)*cos(WALL.B[2]*W2) + WALL.B[2]*dw2[i]*cos(WALL.B[1]*W1)*sin(WALL.B[2]*W2) );
	return n;
}

double *normalNon4foldSymm( double *n,bc WALL,double *point,int dimension ) {
/*
   Find the normal to the surface at this point (particleMPC is
   presently ON surface). We take the gradient of
   the more complex non 4-fold symmetry form of the BC surfaces
*/
	double r,phi,theta;
	double cosT,sinT,cosP,sinP;
	double dr,dphi,dtheta;
	double m,l;
	int i;

	m = WALL.ROTSYMM[0];
	l = WALL.ROTSYMM[1];
	r=0.0;
	for( i=0; i<dimension; i++ ) r += ( point[i]-WALL.Q[i] )*( point[i]-WALL.Q[i] );
	r=sqrt(r);
	phi=atan2( point[1]-WALL.Q[1] , point[0]-WALL.Q[0] );
	cosP=cos(0.25*m*phi);
	sinP=sin(0.25*m*phi);
	if( dimension>_2D ) {
		theta=acos( (point[2]-WALL.Q[2])/r );
		cosT=cos(0.25*l*theta);
		sinT=sin(0.25*l*theta);
	}
	else{
		theta=0.0;
		cosT=0.0;
		sinT=1.0;
	}

	// Calculate the derivatives in spherical/polar coordinates
	if( WALL.ABS ) {
		dr = (WALL.P[3]*WALL.R*WALL.R/r/r/r) * smrtPow( fabs(WALL.R/r) , WALL.P[3]-2.0 );
		dphi = (-0.25*m*WALL.P[0]/WALL.A[0]/WALL.A[0])*smrtPow(sinT,2)*sinP*cosP*smrtPow( fabs(cosP*sinT/WALL.A[0]),WALL.P[0]-2.0 );
		dphi += (-0.25*m*WALL.P[1]/WALL.A[1]/WALL.A[1])*smrtPow(sinT,2)*sinP*cosP*smrtPow( fabs(sinP*sinT/WALL.A[1]),WALL.P[1]-2.0 );
		if( dimension>_2D ) {
			dtheta = (0.25*l*WALL.P[0]/WALL.A[0]/WALL.A[0])*sinT*cosT*smrtPow(cosP,2)*smrtPow( fabs(cosP*sinT/WALL.A[0]),WALL.P[0]-2.0 );
			dtheta += (0.25*l*WALL.P[1]/WALL.A[1]/WALL.A[1])*sinT*cosT*smrtPow(sinP,2)*smrtPow( fabs(sinP*sinT/WALL.A[1]),WALL.P[1]-2.0 );
			dtheta -= (0.25*l*WALL.P[2]/WALL.A[2]/WALL.A[2])*sinT*cosT*smrtPow( fabs(cosT/WALL.A[2]),WALL.P[2]-2.0 );
		}
		else dtheta=0.0;
	}
	else {
		dr = (WALL.P[3]/r) * smrtPow( WALL.R/r , WALL.P[3] );
		dphi = (-0.25*m*WALL.P[0]/WALL.A[0]) * smrtPow( sinT,WALL.P[0] ) * sinP * smrtPow( cosP/WALL.A[0],WALL.P[0]-1.0 );
		dphi += (0.25*m*WALL.P[1]/WALL.A[1]) * smrtPow( sinT,WALL.P[1] ) * cosP * smrtPow( sinP/WALL.A[1],WALL.P[1]-1.0 );
		if( dimension>_2D ) {
			dtheta = (0.25*l*WALL.P[0]/WALL.A[0]) * smrtPow( cosP,WALL.P[0] ) * cosT * smrtPow( sinT/WALL.A[0],WALL.P[0]-1.0 );
			dtheta += (0.25*l*WALL.P[1]/WALL.A[1]) * smrtPow( sinP,WALL.P[1] ) * cosT * smrtPow( sinT/WALL.A[1],WALL.P[1]-1.0 );
			dtheta -= (0.25*l*WALL.P[2]/WALL.A[2]) * sinT * smrtPow( cosT/WALL.A[2],WALL.P[2]-1.0 );
		}
		else dtheta=0.0;
	}

	//Cartesian derivatives
	if( dimension>_2D ) {
		n[0] = cos(phi)*sin(theta)*dr - (sin(phi)/r/sin(theta))*dphi + (cos(phi)*cos(theta)/r)*dtheta;
		n[1] = sin(phi)*sin(theta)*dr + (cos(phi)/r/sin(theta))*dphi + (sin(phi)*cos(theta)/r)*dtheta;
		n[2] = cos(theta)*dr - (sin(theta)/r)*dtheta;
	}
	else{
		n[0] = cos(phi)*dr - (sin(phi)/r)*dphi;
		n[1] = sin(phi)*dr + (cos(phi)/r)*dphi;
	}

	return n;
}
void rudimentaryPBC( particleMPC *pp,int axis ) {
	/*
	   This routine does a "hard" PBC along a given axis
	*/
	double thisXYZ = (double)XYZ[axis];

	if( pp->Q[axis] < 0.0 ) pp->Q[axis] += thisXYZ;
	else if( pp->Q[axis] > thisXYZ ) pp->Q[axis] -= thisXYZ;
}
void rudimentaryBBBC( particleMPC *pp,int axis ) {
	/*
	   This routine does a "hard-code" bounce-back wall along a given axis (impermiable, no-slip walls above and below)
	*/
	int d;
	double thisXYZ = (double)XYZ[axis];

	if( pp->Q[axis] < 0.0 ) {
		pp->Q[axis] *= -1.0;
		for( d=0; d<DIM; d++ ) pp->V[d] *= -1.0;
	}
	else if( pp->Q[axis] > thisXYZ ) {
		pp->Q[axis] = 2.0*thisXYZ - pp->Q[axis];
		for( d=0; d<DIM; d++ ) pp->V[d] *= -1.0;
	}
}
void rudimentaryPBC_box( particleMPC *pp ) {
	/*
	   This routine does a "hard" PBC along a given axis
	*/
	int d;
	for( d=0; d<DIM; d++ ) rudimentaryPBC( pp,d );
}
void rudimentaryBBBC_box( particleMPC *pp ) {
	/*
	   This routine does a "hard" PBC along a given axis
	*/
	int d;
	for( d=0; d<DIM; d++ ) rudimentaryBBBC( pp,d );
}
void rudimentaryChannel_x( particleMPC *pp ) {
	/*
	   This routine does a "hard-code" BBBC along x and PBC along y and z (if 3D)
	*/
	rudimentaryBBBC( pp,0 );
	rudimentaryPBC( pp,1 );
	if( DIM>=_3D ) rudimentaryPBC( pp,2 );
}
void rudimentaryChannel_y( particleMPC *pp ) {
	/*
	   This routine does a "hard-code" BBBC along y and PBC along x and z (if 3D)
	*/
	rudimentaryBBBC( pp,1 );
	rudimentaryPBC( pp,0 );
	if( DIM>=_3D ) rudimentaryPBC( pp,2 );
}
void rudimentaryChannel_z( particleMPC *pp ) {
	/*
	   This routine does a "hard-code" BBBC along z (if 3D) and PBC along x and y
	*/
	if( DIM>=_3D ) rudimentaryBBBC( pp,2 );
	rudimentaryPBC( pp,0 );
	rudimentaryPBC( pp,1 );
}
