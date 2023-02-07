///
/// @file
/// @brief Math functions applied in MPCD
///
/// Different math methods are collecetd there for easy access in any place of the code
# include<math.h>
# include<stdio.h>
# include<stdlib.h>

# include "../headers/definitions.h"
# include "../headers/SRDclss.h"
# include "../headers/pout.h"
# include "../headers/ctools.h"
# include "../headers/bc.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************** MATH ROUTINES ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
///
/// @brief A "smart" Pow method that will only call C-math pow if necessary (non-natural y).
///
/// First checks if y is int and below the smart y limit if not then do just C-math pow
///
/// @param x any real number
/// @param y any real number
///
double smrtPow(double x, double y){

	const int yLim = 10; // an arbitrary limit for smart y usage

	//check if y is int and below the smart y limit
	if ((y - (int)y == 0) && (y <= yLim)){
		int i;
		double result = 1;
		for (i = 0; i < y; i++) result *= x; // dumb power

		return result;
	} else return pow(x, y); // otherwise just do C-math pow
}
///
/// @brief Check if two floats are equal within set TOL
///
/// TOL is real number that is defined in definitions.h
///
///@param x any real number
///@param y any real number
///
int feq(double x,double y) {
		return fabs(x-y)<=TOL;
}
///
/// @brief Check if two floats are NOT equal within set TOL
///
/// TOL is real number that is defined in definitions.h
///
///@param x any real number
///@param y any real number
///
int fneq(double x,double y) {
		return fabs(x-y)>=TOL;
}
///
/// @brief Evaluates the Levi Civita index counter
///
/// @param i first index
/// @param j second index
/// @param k third index
///
int levicivita( int i,int j,int k ) {

	signed int result;
	if( (i==1 && j==2 && k==3) || (i==3 && j==1 && k==2) || (i==2 && j==3 && k==1) ) {
		result = 1;
	}
	else if( (i==3 && j==2 && k==1) || (i==1 && j==3 && k==2) || (i==2 && j==1 && k==3) ) {
		result = -1;
	}
	else result = 0;
	return result;
}
///
/// @brief Takes the dot product of two vectors and returns a scalar
///
/// @param x first vector
/// @param y second vector
/// @param dimension dimensionality of the vectors
///
double dotprod( double x[], double y[],int dimension ) {
	int i;
	double result = 0.;
	for( i=0; i<dimension; i++ ) {
		result += x[i] * y[i];
	}
	return result;
}
///
/// @brief Takes the dot product of a matrix to a vector and returns a vector
///
/// @param M matrix
/// @param v vector
/// @param result output vector
/// @param dimension dimensionality of the operands
///
void dotprodMatVec( double M[][3],double v[],double result[],int dimension ) {
	int i,j;
	for( i=0; i<dimension; i++ ) result[i] = 0.;
	for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) result[i] += M[i][j]*v[j];
}
///
/// @brief Takes the dot product of a vector to a matrix and returns a vector
///
/// @param v vector
/// @param M matrix
/// @param result output vector
/// @param dimension dimensionality of the operands
///
void dotprodVecMat( double v[], double M[][3],double result[],int dimension ) {
	int i,j;
	for( i=0; i<dimension; i++ ) result[i] = 0.;
	for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) result[i] += v[j] * M[j][i];
}
///
/// @brief Takes the dot product of a matrix to another matrix and returns a matrix
///
/// @param A first matrix
/// @param B second matrix
/// @param result output matrix
/// @param dimension dimensionality of the operands
///
void dotprodMatMat( double A[][3],double B[][3],double result[][3],int dimension ) {
	int i,j,k;
	for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) result[i][j] = 0.;
	for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) for( k=0; k<dimension; k++ ) result[i][j] += A[i][k]*B[k][j];
}
///
/// @brief Takes the cross product of two 3D vectors and sets it as the third
///	Must be 3D since in 2D, result is in 3rd dimension
///
/// @param x first vector
/// @param y second vector
/// @param result output vector
///
void crossprod( double x[3], double y[3], double result[3] ) {
	int i;
	for( i=0; i<_3D; i++ ) result[i]=0.; // init
	// manually compute cross product terms
	result[0] = x[1]*y[2] - x[2]*y[1];
	result[1] = x[2]*y[0] - x[0]*y[2];
	result[2] = x[0]*y[1] - x[1]*y[0];
}
///
/// @brief Old version of the cross product operation.
///	 This version was found to be slow (gprof said this and it's calls to
///	 levicivita took >35% runtime!!!!)
///
/// @param x first vector
/// @param y second vector
/// @param result output vector
///
void oldcrossprod( double x[3], double y[3], double result[3] ) {
	int i,j,k;
	signed int eps;
	for( i=0; i<_3D; i++ ) result[i]=0.;
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) {
				eps = levicivita( i+1,j+1,k+1 );
				result[i] += ((double) eps) * x[j] * y[k];
	}
}
///
/// @brief Finds the outer product of two vectors, returns a matrix
///
/// @param x first vector
/// @param y second vector
/// @param result output matrix
/// @param dimension dimensionality of the operands
///
void outerprod( double x[], double y[], double result[][_3D],int dimension ) {
	int i,j;
	for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) result[i][j] = 0.;
	for( i=0; i<dimension; i++ ) for( j=0; j<dimension; j++ ) result[i][j] = x[i]*y[j];
}
///
/// @brief Finds the magnitude of the vector, returns scalar value
///
/// @param x input vector
/// @param dimension dimensionality of the vector
///
double length( double x[],int dimension ) {
	int i;
	double result = 0.;
	for( i=0; i<dimension; i++ ) result += x[i] * x[i];
	result = sqrt( result );
	return result;
}
///
/// @brief Normalizes the input vector
///
/// @param x input vector
/// @param dimension dimensionality of the vector
///
void norm( double x[],int dimension ) {
	int i;
	double l = 0.;
	l = length( x,dimension );
	if( fneq(l,0.0) ) for ( i=0; i<dimension; i++ ) x[i] = x[i] / l;
}
///
/// @brief Normalizes the input vector and returns as separate one
///
/// @param xin input vector
/// @param xout output vector
/// @param dimension dimensionality of the vectors
///
void normCopy( double xin[],double xout[],int dimension ) {
	int i;
	double l = 0.;
	l = length( xin,dimension );
	if( fneq(l,0.0) ) for ( i=0; i<dimension; i++ ) xout[i] = xin[i] / l;
}
///
/// @brief Finds the normal vector (n) to a plane defined by x and y
///  
///	 It assumes 3D because even in 2D, the result must be in 3rd dimension
///
/// @param x first input vector
/// @param y second input vector
/// @param n dimensionality of the vectors
///
void normalplane( double x[3], double y[3], double n[3] ) {
	crossprod( x,y,n );
	norm( n,_3D );
}
///
/// @brief Gives the vector projection of v onto n (which is most often the normal of a plane-normal compontent)
///
/// @param v input vector to project
/// @param n input vector which is used for projection to
/// @param VN output vector projection
/// @param dimension dimensionality of the vectors
///
void proj( double v[],double n[],double VN[],int dimension ) {
	int i;
	double x;
	x = dotprod( v,n,dimension );
	for( i=0; i<dimension; i++ ) VN[i] = x*n[i];
}
///
/// @brief Gives the tangential component of the vector
///
/// @param v input vector
/// @param VN input normal component of the vector
/// @param VT output tangential component of the vector
/// @param dimension dimensionality of the vectors
///
void tang( double v[],double VN[],double VT[],int dimension ) {
	int i;
	for( i=0; i<dimension; i++ ) VT[i] = v[i] - VN[i];
}
///
/// @brief Returns the cosign of the angle between two vectors
///
/// @param v1 first input vector
/// @param v2 second input vector
/// @param dimension dimensionality of the vectors
///
double cosang( double v1[],double v2[],int dimension ) {
	double cosa;
	cosa = dotprod( v1,v2,dimension );
	cosa /= length( v1,dimension );
	cosa /= length( v2,dimension );
	return cosa;
}
///
/// @brief arctan that returns a signed angle
/// 
/// @param y first input scalar
/// @param x second input scalar
///
double atan2( double y,double x ) {
	double at=0.0;
	if( x>0.0 ) at=atan(y/x);
	else if( x<0.0 && y>=0.0 ) at=atan(y/x)+pi;
	else if( x<0.0 && y<0.0 ) at=atan(y/x)-pi;
	else if( feq(x,0.0) && y>0.0 ) at=0.5*pi;
	else if( feq(x,0.0) && y<0.0 ) at=-0.5*pi;

	//Map from (-pi,pi] to [0,2pi)
	//if( at<0.0 )at+=2.0*pi;
	return at;
}
///
/// @brief Finds the Unsigned angle between two vectors, return scalar
///
/// @param v1 first input vector
/// @param v2 second input vector
/// @param dimension dimensionality of the vectors
///
double absAngle( double v1[], double v2[], int dimension ) {
	return acos( cosang(v1,v2,dimension) );
}
///
///	@brief Finds the SIGNED angle between two vectors
///	If not in 3D, need to set to be 3D
///	// |A·B| = |A| |B| COS(θ)
///	// |A×B| = |A| |B| SIN(θ)
///	return Math.Atan2(Cross(A,B), Dot(A,B));
///
/// @param v1 first input vector
/// @param v2 second input vector
/// @param dimension dimensionality of the vectors
///
double signedAngle( double v1[], double v2[], int dimension ) {
	double A[_3D],B[_3D],cross[_3D];
	double s,c;
	int i;

	for( i=0; i<_3D; i++ ) A[i]=0.0;
	for( i=0; i<_3D; i++ ) B[i]=0.0;
	for( i=0; i<dimension; i++ ) A[i]=v1[i];
	for( i=0; i<dimension; i++ ) B[i]=v2[i];

	crossprod( A,B,cross );
	s = length(cross,_3D);
	c = dotprod(A,B,dimension);
	return atan2(s,c);
}
///
///	@brief Calculates the distance between two points
///
/// @param P1 first input point
/// @param P2 second input point
/// @param dimension dimensionality of the points
///
double distpoints( double P1[_3D],double P2[_3D],int dimension ) {
	double dist = 0.;
	int i;

	for( i=0; i<dimension; i++ ) dist += (P2[i]-P1[i]) * (P2[i]-P1[i]);
	return sqrt(dist);
}
///
///	@brief Calculates the distance from a point to a surface
///
/// There is still a question regarding planar surfaces
///
/// @param WALL input boundary
/// @param P input point
///
double distsurf( bc WALL,double P[_3D] ) {
	double len = 0.;
	double dist = 0.;
	int i;
	//I'm not sure if this works for non-planar surfaces!!!
	for( i=0; i<_3D; i++ ) dist += WALL.A[i] * P[i];
 	dist -= WALL.R;
	for( i=0; i<_3D; i++ ) len += WALL.A[i] * WALL.A[i];
	len = sqrt( len );
	dist = fabs( dist ) / len;
	return dist;
}
///
///	@brief Calculates the distance from a point to a plane
///
/// @param WALL input boundary
/// @param x input point x coordinate
/// @param y input point y coordinate
/// @param z input point z coordinate
///
double distplane( bc WALL,double x, double y, double z ) {
	double len = 0.;
	double dist = 0.;
	int i;
	dist = WALL.A[0] * x;
	dist += WALL.A[1] * y;
	dist += WALL.A[2] * z;
	dist -= WALL.R;
	for( i=0; i<_3D; i++ ) len += WALL.A[i] * WALL.A[i];
	len = sqrt( len );
	dist = fabs( dist ) / len;
	return dist;
}
///
///	@brief Calculates Pythagorean theorem
///
/// @param x input real number
/// @param y input real number
///
double pythag( double x, double y ) {
	return sqrt( x*x + y*y );
}
///
///	@brief Magnitude of x times sign of y
///
/// @param x input real number
/// @param y input real number
///
double SIGN( double x,double y ) {
	if( y>0. ) return fabs(x);
	else return -1.*fabs(x);
}
///
///	@brief Numerically calculates the moment of inertia of a bc structure.
///		WARNING!!! CURRENTLY JUST APPROXIMATES EVERYTHING AS THE CLOSEST ELLIPSOID!!!
///
/// @param body input boudary condition
/// @param XYZ parameters of the control volume (not in use in this function)
/// @param dimension dimensionality of the structure
///
void latticeEstMomInert( bc *body,int XYZ[],int dimension ) {
	int i,j;
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) body->I[i][j] = 0.0;

	if( dimension==_3D ) {
			body->I[0][0] = body->MASS*( body->AINV[1]*body->AINV[1] + body->AINV[2]*body->AINV[2])*body->R*body->R/5.0;
			body->I[1][1] = body->MASS*( body->AINV[0]*body->AINV[0] + body->AINV[2]*body->AINV[2])*body->R*body->R/5.0;
			body->I[2][2] = body->MASS*( body->AINV[0]*body->AINV[0] + body->AINV[1]*body->AINV[1])*body->R*body->R/5.0;
	}
	else if( dimension==_2D ) {
			body->I[0][0] = body->MASS*( body->AINV[1]*body->AINV[1] )*body->R*body->R/5.0;
			body->I[1][1] = body->MASS*( body->AINV[0]*body->AINV[0] )*body->R*body->R/5.0;
			body->I[2][2] = body->MASS*( body->AINV[0]*body->AINV[0] + body->AINV[1]*body->AINV[1])*body->R*body->R/5.0;
	}
	else {
		printf( "Warning: Moment of inertia tensor zero because dimensionality not 3 or 2D." );
	}
}
///
///	@brief Calculates the moment of inertia of a bc structure.
/// Formula for it is different depending on the boundary conditions type
///
/// @param body input boudary condition
/// @param XYZ parameters of the control volume
/// @param dimension dimensionality of the inputs
///
void mominert( bc *body,int XYZ[],int dimension ) {
	int i,j;
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) body->I[i][j] = 0.0;

	if( dimension==_3D ) {
		//Planes (set to zero)
		if( feq(body->P[0],1.0) && feq(body->P[1],1.0) && feq(body->P[2],1.0) ) {
			for( i=0; i<_3D; i++ ) body->I[i][i] = 0.0;
		}
		//Sphere & ellipsoids
		if( feq(body->P[0],2.0) && feq(body->P[1],2.0) && feq(body->P[2],2.0) ) {
			//Sphere
			if( body->A[0]==1.0 && body->A[1]==1.0 && body->A[2]==1.0 ) {
				for( i=0; i<_3D; i++ ) body->I[i][i] = 2.0*body->MASS*body->R*body->R/5.0;
			}
			//Ellipsoids
			else {
				body->I[0][0] = body->MASS*( body->AINV[1]*body->AINV[1] + body->AINV[2]*body->AINV[2])*body->R*body->R/5.0;
				body->I[1][1] = body->MASS*( body->AINV[0]*body->AINV[0] + body->AINV[2]*body->AINV[2])*body->R*body->R/5.0;
				body->I[2][2] = body->MASS*( body->AINV[0]*body->AINV[0] + body->AINV[1]*body->AINV[1])*body->R*body->R/5.0;
			}
		}
		//Cylinders
		// Cylinder along z
		else if( feq(body->P[0],2.0) && feq(body->P[1],2.0)  && feq(body->A[2],0.0) ) {
			body->I[0][0] = body->MASS*( 3.0*body->R*body->R + (double)XYZ[2] )/12.0;
			body->I[1][1] = body->MASS*( 3.0*body->R*body->R + (double)XYZ[2] )/12.0;
			body->I[2][2] = body->MASS*body->R*body->R/2.0;
		}
		// Cylinder along y
		else if( feq(body->P[0],2.0) && feq(body->P[2],2.0)  && feq(body->A[1],0.0) ) {
			body->I[0][0] = body->MASS*( 3.0*body->R*body->R + (double)XYZ[1] )/12.0;
			body->I[1][1] = body->MASS*body->R*body->R/2.0;
			body->I[2][2] = body->MASS*( 3.0*body->R*body->R + (double)XYZ[1] )/12.0;
		}
		// Cylinder along x
		else if( feq(body->P[1],2.0) && feq(body->P[2],2.0)  && feq(body->A[0],0.0) ) {
			body->I[0][0] = body->MASS*body->R*body->R/2.0;
			body->I[1][1] = body->MASS*( 3.0*body->R*body->R + (double)XYZ[0] )/12.0;
			body->I[2][2] = body->MASS*( 3.0*body->R*body->R + (double)XYZ[0] )/12.0;
		}
		//All others
		else latticeEstMomInert( body,XYZ,dimension );
	}
	else if( dimension==_2D ) {
		//Planes (set to zero)
		if( feq(body->P[0],1.0) && feq(body->P[1],1.0) ) for( i=0; i<_3D; i++ ) body->I[i][i] = 0.0;
		//Circle & ellipses
		if(
			feq(body->P[0],2.0) && feq(body->P[1],2.0) ) {
			//Disc/Circle
			if( feq(body->A[0],1.0) && feq(body->A[1],1.0) ) {
				body->I[2][2] = body->MASS*body->R*body->R/2.0;
				body->I[1][1] = body->MASS*body->R*body->R/4.0;
				body->I[0][0] = body->MASS*body->R*body->R/4.0;
			}
			//Ellipses
			else {
				body->I[0][0] = body->MASS*( body->AINV[1]*body->AINV[1] )*body->R*body->R/5.0;
				body->I[1][1] = body->MASS*( body->AINV[0]*body->AINV[0] )*body->R*body->R/5.0;
				body->I[2][2] = body->MASS*( body->AINV[0]*body->AINV[0] + body->AINV[1]*body->AINV[1])*body->R*body->R/5.0;
			}
		}
		//All others
		else latticeEstMomInert( body,XYZ,dimension );
	}
	else {
		printf( "Warning: Moment of inertia tensor zero because dimensionality not 3 or 2D." );
	}
}
///
///	@brief This routine numerically estimates the volume of the BC object.
///		WARNING!!! CURRENTLY JUST APPROXIMATES EVERYTHING AS THE CLOSEST ELLIPSOID!!!
///
/// @param body input boudary condition
/// @param XYZ parameters of the control volume
/// @param dimension dimensionality of the inputs
///
double latticeEstVol( bc *body,int XYZ[],int dimension ) {
	double vol = 0.0;
	if( dimension==_3D ) vol = 4.0*pi*( body->AINV[0]*body->AINV[1]*body->AINV[2]*smrtPow(body->R,3) )/3.0;
	else if( dimension==_2D ) vol = pi*body->AINV[0]*body->AINV[1]*smrtPow(body->R,2);
	else {
		printf( "Warning: Volume zero because dimensionality not 3 or 2D." );
	}
	return vol;
}
///
///	@brief This routine returns the volume of the BC object
///   for dimension=3 and the area for dimension=2
///   i.e. it returns the dimension-dimensional volume
///
/// @param body input boudary condition
/// @param XYZ parameters of the control volume
/// @param dimension dimensionality of the inputs
///
void dim_vol( bc *body,int XYZ[],int dimension ) {
	body->VOL = 0.0;
	if( dimension==_3D ) {
		//Sphere & ellipsoids
		if( feq(body->P[0],2.0) && feq(body->P[1],2.0) && feq(body->P[2],2.0) ) {
			//Sphere
			if( feq(body->A[0],1.0) && feq(body->A[1],1.0) && feq(body->A[2],1.0) ) body->VOL = 4.0*pi*smrtPow( body->R,3.0 )/3.0;
			//Ellipsoids
			else body->VOL = 4.0*pi*( body->AINV[0]*body->AINV[1]*body->AINV[2]*smrtPow(body->R,3) )/3.0;
		}
		//All others
		else body->VOL = latticeEstVol( body,XYZ,dimension );
	}
	else if( dimension==_2D ) {
		//Circle & ellipses
		if( feq(body->P[0],2.0) && feq(body->P[1],2.0) ) {
			//Circle
			if( feq(body->A[0],1.0) && feq(body->A[1],1.0) ) body->VOL = pi * smrtPow( body->R,body->P[3] );
			//Ellipses
			else body->VOL = pi*body->AINV[0]*body->AINV[1]*smrtPow(body->R,2);
		}
		//All others
		else body->VOL = latticeEstVol( body,XYZ,dimension );
	}
	else {
		printf( "Warning: Volume zero because dimensionality not 3 or 2D." );
	}
}
///
///	@brief The random initializers may give a net momentum
///   to the system. We do not allow this by doing a
///   Galilean transformation to rest frame.
///
///  The function finds the total net momentum of the whole system and subtract it from the velocities of all the stuff
///  
/// @param pp reference to the particles in the system
/// @param WALL reference to the boundaries in the system
/// @param simMD reference to the simulation parameters
/// @param KBT Thermal energy
/// @param VEL The average speed of the particles in the system
/// @param POP Total number of particles in the system
/// @param NBC Total number of of boundaries present in the system
/// @param MDmode The MD coupling mode
/// @param dimension Dimensions of the system
///
void galileantrans( particleMPC *pp,bc WALL[],simptr simMD,spec SP[],double KBT,double VEL[],int POP,int NBC,int MDmode,int dimension ) {
	int i,j;
	double NET[_3D];		//Net momentum
	double M,totM=0.0;			//Mass and total mass

	for( i=0; i<_3D; i++ ) NET[i] = 0.;
	// Sum up the net momentum
	for( i=0; i<POP; i++ ) {
		M = SP[(pp+i)->SPID].MASS;
		totM += M;
		for( j=0; j<dimension; j++ ) NET[j] += M*(pp+i)->V[j];
	}
	for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
		M = WALL[i].MASS;
		totM += M;
		for( j=0; j<dimension; j++ ) NET[j] += M*WALL[i].V[j];
	}
	if( MDmode == MDinMPC ) for( i=0; i<(simMD->atom.n); i++ ){
		M = (double) (simMD->atom.items+i)->mass;
		totM += M;
		NET[0] += M*(simMD->atom.items+i)->vx;
		if(dimension>=_2D) NET[1] += M*(simMD->atom.items+i)->vy;
		if(dimension>=_3D) NET[2] += M*(simMD->atom.items+i)->vz;
	}

	//Average the net momentum
	for( i=0; i<dimension; i++ ) NET[i] /= totM;
	//Take into account that the user may have GIVEN the system an average velocity VEL
	for( i=0; i<dimension; i++ ) NET[i] -= VEL[i];

	// Subtract off the net momentum
	for( i=0; i<POP; i++ ) for( j=0; j<dimension; j++ ) (pp+i)->V[j] -= NET[j];
	for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) for( j=0; j<dimension; j++ ) WALL[i].V[j] -= NET[j];
	if( MDmode == MDinMPC ) for( i=0; i<(simMD->atom.n); i++ ){
		(simMD->atom.items+i)->vx -= NET[0];
		if(dimension>=_2D) (simMD->atom.items+i)->vy -= NET[1];
		if(dimension>=_3D) (simMD->atom.items+i)->vz -= NET[2];
	}
}
///
///	@brief Zeros the components of the positions and
///   velocities of objects that are greater
///   dimension than the simulation (just paranoid)
///  
/// @param pp reference to the particles in the system
/// @param WALL reference to the boundaries in the system
/// @param simMD reference to the simulation parameters
/// @param GPOP Total number of particles in the system
/// @param NBC Total number of of boundaries present in the system
/// @param MDmode The MD coupling mode
/// @param dimension Dimensions of the system
///
void zeroExtraDims( particleMPC *pp,bc WALL[],simptr simMD,int GPOP,int NBC,int MDmode,int dimension ) {
	int i;
	if( dimension<_3D ) {
		for( i=0; i<GPOP; i++ ) {
			pp[i].Q[2] = 0.;
			pp[i].V[2] = 0.;
		}
		for( i=0; i<NBC; i++ ) {
			WALL[i].Q[2] = 0.;
			WALL[i].V[2] = 0.;
		}
		if( MDmode == MDinMPC ) for( i=0; i<(simMD->atom.n); i++ ) {
			(simMD->atom.items+i)->rz = 0.;
			(simMD->atom.items+i)->vz = 0.;
		}
	}
	if( dimension<_2D ) {
		for( i=0; i<GPOP; i++ ) {
			pp[i].Q[1] = 0.;
			pp[i].V[1] = 0.;
		}
		for( i=0; i<NBC; i++ ) {
			WALL[i].Q[1] = 0.;
			WALL[i].V[1] = 0.;
		}
		if( MDmode == MDinMPC ) for( i=0; i<(simMD->atom.n); i++ ) {
			(simMD->atom.items+i)->ry = 0.;
			(simMD->atom.items+i)->vy = 0.;
		}
	}
}
///
///	@brief Generic histogram binning algorithm
///
///  Create a histogram from the values in the input
///  
/// @param values actual values of the velocity in every MPCD cell
/// @param hist histogram of the velocity in each of the 3D component
/// @param minRange Minimum range of the histogram
/// @param maxRange Maximum range of the histogram
/// @param POP Total volume of the system
///
void histbin( double values[],int hist[BINS],double minRange,double maxRange,int POP ) {
	int i,bin,binsM1;
	double invDenom;

	binsM1=BINS-1;
	invDenom=1./( maxRange-minRange );
	for( i=0; i<POP; i++ ) {
		bin=binsM1*( values[i]-minRange )*invDenom;
		if(bin<0) printf( "\t%d\n",bin );
		else if (bin>=BINS) {
			printf( "\twell shit\n" );
			printf( "\tmin=%lf,max=%lf,\tvalue=%lf\n",minRange,maxRange,values[i] );
			printf( "\t%d\n",bin );
		}
		// printf( "\t%d\n",bin );
		else hist[bin]++;
	}
}
///
///	@brief This routine is the parallel axis theorem.
///    It takes a inertia tensor I about the centre
///    of mass and calculates the I about a
///    displaced by R
///
///  
/// @param I input an inertia tensor
/// @param R input displacement
/// @param M input mass
/// @param dimension Dimensions of the input values
///
void parallelaxis( double I[][_3D],double R[],double M,int dimension ) {
	int i,j,k;

	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) {
		I[i][j] -= R[i] * R[j];
		if( i == j ) for( k=0; k<_3D; k++ ) I[i][j] += R[k] * R[k];
		I[i][j] *= M;
	}
}
///
///	@brief Operate in the frame of reference of the bc.
///    The routine labframe must proceed it.
///		
///    Subtract the velocity of the walls from the velocity of the objects
///  
/// @param V input velocity
/// @param WALL reference to the boundary
/// @param dimension Dimensions of the input values
///
void restframe( double V[],bc WALL,int dimension ) {
	int i;
	for( i=0; i<dimension; i++ ) V[i] -= WALL.V[i];
}
///
///	@brief Operate in the frame of reference of the bc.
///    The routine restframe must preceed it.
///		
///    Add the velocity of the walls to the velocity of the objects
///  
/// @param V input velocity
/// @param WALL reference to the boundary
/// @param dimension Dimensions of the input values
///
void labframe( double V[],bc WALL,int dimension ) {
	int i;
	for( i=0; i<dimension; i++ ) V[i] += WALL.V[i];
}
///
///	@brief Finds the determinant of a 2x2 matrix
///
///
/// @param m input 2x2matrix
///
double det2x2( double m[_2D][_2D] ) {
	return m[0][0]*m[1][1] - m[0][1]*m[1][0];
}
///
///	@brief Finds the determinant of a 3x3 matrix
///
///
/// @param m input 3x3matrix
///
double det3x3( double m[_3D][_3D] ) {
	double c1,c2,c3;

	c1 = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]);
	c2 = m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]);
	c3 = m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
	return( c1 - c2 + c3);
}
///
///	@brief Finds the determinant of a nxn matrix (not more than 3x3)
///    Recursive definition of determinate using expansion by minors
///    Stolen from http://paulbourke.net/miscellaneous/determinant/
///    BUT
///    I HATE passing to a double pointer so I'll just stick to det2x3 and det3x3
///
///
/// @param a reference to the nxn matrix
/// @param n dimensionality of the matrix
///
double determinant( double **a,int n ) {
	int i,j,j1,j2;
	double det = 0.;
	double **m = NULL;

	if( n<1 ) { /* Error */ }
	else if( n==1 ) det = a[0][0]; //Shouldn't get used
	else if( n==2 ) det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	else if( n==3 ) {
		det = 0.;
		for( j1=0; j1<n; j1++ ) {
			m = malloc((n-1)*sizeof(double *));
			for( i=0; i<n-1; i++ ) m[i] = malloc((n-1)*sizeof(double));
			for( i=0; i<n-1; i++ ) for( j=0; j<n-1; j++ ) m[i][j] =0.0;
			for( i=1; i<n; i++ ) {
				j2 = 0;
				for( j=0; j<n; j++ ) {
					if (j == j1) continue;
					m[i-1][j2] = a[i][j];
					j2++;
				}
			}
			det += smrtPow(-1.0,1.0+j1+1.0) * a[0][j1] * determinant( m,n-1 );
			for( i=0; i<n-1; i++ ) free( m[i] );
			free( m );
		}
	}
	else {
		printf( "Error: Determinant only programmed for 3x3 max size.\n" );
		exit(EXIT_FAILURE);
	}

	return det;
}
///
///	@brief Finds the trace of the matrix
///
///
/// @param a reference to the nxn matrix
/// @param n dimensionality of the matrix
///
double trace( double **a,int n ) {
	int i;
	double tr=0.;
	for( i=0; i<n; i++ ) tr+=a[i][i];
	return tr;
}
///
///	@brief Inverts a 2x2 matrix
///
///
/// @param m reference to the 2x2 matrix
///
void invert2x2( double m[_2D][_2D] ) {
	double det;
	double n[_2D][_2D];		//The inverted matrix
	int i,j;

	if( fabs( m[0][0] ) <= TOL && fabs( m[0][1] ) <= TOL && fabs( m[1][0] ) <= TOL && fabs( m[1][1] ) <= TOL ) {
		//Zero matrix: Leave as zero
		return;
	}
	else if( fabs( det2x2( m ) ) <= TOL ) {
	//else if( fabs(determinant( &m[0],2 ) ) <= TOL ) {
		//Warning: Small determinant. Consider matrix singular: leave matrix untouched
		return;
	}
	else{
		det = det2x2( m );

		n[0][0] = m[1][1];
		n[0][1] = -1. * m[0][1];
		n[1][0] = -1. * m[1][0];
		n[1][1] = m[0][0];

		for( i=0; i<_2D; i++ ) for( j=0; j<_2D; j++ ) n[i][j] /= det;
		for( i=0; i<_2D; i++ ) for( j=0; j<_2D; j++ ) m[i][j] = n[i][j];
	}
}
///
///	@brief Returns the i,j cofactor for a 3x3 matrix m
///
///
/// @param m reference to the 3x3 matrix
/// @param i cofactor index
/// @param j cofactor index
///
double cofactor3x3( double m[_3D][_3D],int i,int j ) {
	double a00,a01,a10,a11;
	int c;

	c = ( 2*((i+j)/2) == (i+j) ) ? 1 : -1;

	if( i==0 ) {
		if( j==0 ) {
			a00 = m[1][1]; a01 = m[1][2];
			a10 = m[2][1]; a11 = m[2][2];
		}
		else if( j==1 ) {
			a00 = m[1][0]; a01 = m[1][2];
			a10 = m[2][0]; a11 = m[2][2];
		}
		else {
			a00 = m[1][0]; a01 = m[1][1];
			a10 = m[2][0]; a11 = m[2][1];
		}
	}
	else if( i==1 ) {
		if( j==0 ) {
			a00 = m[0][1]; a01 = m[0][2];
			a10 = m[2][1]; a11 = m[2][2];
		}
		else if( j==1 ) {
			a00 = m[0][0]; a01 = m[0][2];
			a10 = m[2][0]; a11 = m[2][2];
		}
		else {
			a00 = m[0][0]; a01 = m[0][1];
			a10 = m[2][0]; a11 = m[2][1];
		}
	}
	else {
		if( j==0 ) {
			a00 = m[0][1]; a01 = m[0][2];
			a10 = m[1][1]; a11 = m[1][2];
		}
		else if( j==1 ) {
			a00 = m[0][0]; a01 = m[0][2];
			a10 = m[1][0]; a11 = m[1][2];
		}
		else {
			a00 = m[0][0]; a01 = m[0][1];
			a10 = m[1][0]; a11 = m[1][1];
		}
	}
	return c * ( a00*a11 - a01*a10 );
}
///
///	@brief Returns the i,j cofactor for a 3x3 matrix
///
///
/// @param m_inv reference to the inverted 3x3 matrix
/// @param m reference to the 3x3 matrix to invert
///
void invert3x3(double m_inv[_3D][_3D],double m[_3D][_3D]) {
	double det;
	int i,j;

	det = det3x3(m);
// 	det = determinant( &m[0],3 );
	for( i=0; i<_3D; i++ )for( j=0; j<_3D; j++ ) m_inv[j][i] = cofactor3x3( m,i,j ) / det;
}
///
///	@brief Calculates the energy, the linear momentum
///    and the angular momentum of a point particleMPC
///    and an object
///
///
/// @param VA velocity of the particle
/// @param MA mass of the particle
/// @param QA position of the particle
/// @param VB velocity of the boundary
/// @param MB mass of the boundary
/// @param QB position of the boundary
/// @param IB tensor of the boundary
/// @param dimension dimensionality of the input values
///
void conservation( double VA[],int MA,double QA[],double VB[],int MB,double QB[],double WB[],double IB[_3D][_3D],int dimension ) {
	int i,j;
	double E,TE;
	double P[_3D],TP[_3D];
	double L[_3D],TL[_3D];
	double R[_3D];
	//Zero vectors
	for( i=0; i<dimension; i++ ) {
		P[i] = 0.;
		TP[i] = 0.;
		L[i] = 0.;
		TL[i] = 0.;
	}
	//Energy
	//Kinetic particleMPC
	E = 0.;
	for( i=0; i<dimension; i++ ) E += VA[i]*VA[i];
	E *= 0.5 * MA;
	printf( "Kinetic energy of the particleMPC = %lf\n",E );
	TE = E;
	//Kinetic BC
	E = 0.;
	for( i=0; i<dimension; i++ ) E += VB[i]*VB[i];
	E *= 0.5 * MB;
	printf( "Kinetic energy of the object = %lf\n",E );
	TE += E;
	//Rotational BC
	E = 0.;
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) E += WB[i] * IB[i][j] * WB[j];
	E *= 0.5;
	printf( "Rotational energy of the object = %lf\n",E );
	TE += E;

	//Linear Momentum
	//Momentum particleMPC
	for( i=0; i<dimension; i++ ) {
		P[i] = VA[i]*MA;
		TP[i] += P[i];
	}
	printf( "Linear Momentum of the particleMPC =" );
	pvec( P,dimension );
	//Momentum BC
	for( i=0; i<dimension; i++ ) {
		P[i] = VB[i]*MB;
		TP[i] += P[i];
	}
	printf( "Linear Momentum of the object =" );
	pvec( P,dimension );

	//Angular Momentum (about object centre)
	//Momentum point particleMPC
	for( i=0;i<dimension;i++ ) {
		P[i] = VA[i]*MA;
		R[i] = QA[i] - QB[i];
	}
	crossprod( P,R,L );
	for( i=0; i<_3D; i++ ) TL[i] += L[i];
	printf( "Angular Momentum of the particleMPC (about centre of object) =" );
	pvec( L,_3D );
	//Momentum BC
// 	dotprodmat( WB,IB,L,_3D );
	dotprodMatVec( IB,WB,L,_3D );
	for( i=0; i<_3D; i++ ) TL[i] += L[i];
	printf( "Angular Momentum of the object =" );
	pvec( L,_3D );

	//Output
	printf( "\nTotal Energy: %lf\n",TE );
	printf( "Total Linear Momentum:" );
	pvec( TP,_3D );
	printf( "Total Angular Momentum:" );
	pvec( TL,_3D );
}
///
///	@brief This function calculates W which is used to
///	    determine if boundary conditions should be
///	    applied to a particleMPC. It is a more generic form of calcW()
///		Non-4-fold symmetries
///
///
/// @param WALL input boundary to calculate
/// @param POS position of the particle
/// @param dimension dimensionality of the input values
///
double non4foldSymmCalcW( bc WALL,double POS[], int dimension ) {
	double terms, W=0.0;
	int i;
	double r,phi,theta;
	double cosT,sinT,cosP,sinP;

	r=0.0;
	for( i=0; i<dimension; i++ ) r += ( POS[i]-WALL.Q[i] )*( POS[i]-WALL.Q[i] );
	r=sqrt(r);
	phi=atan2( POS[1]-WALL.Q[1],POS[0]-WALL.Q[0] );
	cosP=cos(0.25*WALL.ROTSYMM[0]*phi);
	sinP=sin(0.25*WALL.ROTSYMM[0]*phi);
	if( dimension>_2D ) {
		theta=acos( (POS[2]-WALL.Q[2])/r );
		cosT=cos(0.25*WALL.ROTSYMM[1]*theta);
		sinT=sin(0.25*WALL.ROTSYMM[1]*theta);
	}
	else{
		theta=0.0;
		cosT=0.0;
		sinT=1.0;
	}
	// First term
	terms = WALL.A[0]*cosP*sinT;
	if( WALL.ABS ) terms=fabs(terms);
	terms = smrtPow( terms,WALL.P[0] );
	W += terms;
	// Second term
	terms = WALL.A[1]*sinP*sinT;
	if( WALL.ABS ) terms=fabs(terms);
	terms = smrtPow( terms,WALL.P[1] );
	W += terms;
	// Third term
	terms = WALL.A[2]*cosT;
	if( WALL.ABS ) terms=fabs(terms);
	terms = smrtPow( terms,WALL.P[2] );
	W += terms;
	// Fourth terms
	terms = WALL.R/r;
	if( WALL.ABS ) terms=fabs(terms);
	terms = smrtPow( terms,WALL.P[3] );
	W -= terms;
	if( WALL.INV ) W *= -1.0;

	return W;
}
///
///	@brief This function evaluates the surface
///   function at the position POS - exactly like calcW
///
///
/// @param WALL input boundary to calculate
/// @param POS position of the particle
/// @param dimension dimensionality of the input values
///
double surf_func( bc WALL,double POS[], int dimension ) {
	double terms, W=0.0;
	int i;

	if( feq(WALL.ROTSYMM[0],4.0) && feq(WALL.ROTSYMM[1],4.0) ) {
		for( i=0; i<dimension; i++ ) {
			terms = WALL.A[i] * ( POS[i]-WALL.Q[i] );
			if( WALL.ABS ) terms=fabs(terms);
			terms = smrtPow( terms,WALL.P[i] );
			W += terms;
		}
		terms = WALL.R;
		if( WALL.ABS ) terms=fabs(terms);
		terms = smrtPow( terms,WALL.P[3] );
		W -= terms;
		//Check if need wavy wall complications
		if( !feq(WALL.B[0],0.0) ) W += calcWavyW(WALL,POS,W);
		//Check if invert wall
		if( WALL.INV ) W *= -1.;
	}
	else {
		W = non4foldSymmCalcW( WALL,POS,dimension );
	}
	return W;
}
///
///	@brief Find the two eigenvalues for m for a 2x2 matrix
///
///
/// @param m reference to the 2x2 matrix
/// @param eigval output Eigenvalues
///
void eigenvalues2x2( double **m,double eigval[] ) {
/*
    
*/
	double det=determinant( m,_2D );
	double trace=m[0][0]+m[1][1];
	double sq=sqrt(trace*trace*0.25-det);
	eigval[0]=trace*0.5 + sq;
	eigval[1]=trace*0.5 - sq;
}
void eigenvectors2x2( double **m,double eigval[],double eigvec[][_2D] ) {
/*
    Find the two eigenvectors (normalized) for m for a 2x2 matrix
*/
	if( fneq(m[1][0],0.0) ) {
		//First eigenvalue
		eigvec[0][0]=eigval[0]-m[1][1];
		eigvec[0][1]=m[1][0];
		norm( eigvec[0],_2D );
		//Second eigenvalue
		eigvec[1][0]=eigval[1]-m[1][1];
		eigvec[1][1]=m[1][0];
		norm( eigvec[1],_2D );
	}
	else if( fneq(m[0][1],0.0) ) {
		//First eigenvalue
		eigvec[0][0]=m[0][1];
		eigvec[0][1]=eigval[0]-m[0][0];
		norm( eigvec[0],_2D );
		//Second eigenvalue
		eigvec[1][0]=m[0][1];
		eigvec[1][1]=eigval[1]-m[0][0];
		norm( eigvec[1],_2D );
	}
	else {
		//Matrix is diagonal
		if( feq(m[0][0],m[1][1]) ) {
			//The eigenvalues are the same
			//eigval[0]==eigval[1] == m[0][0]==m[1][1]
			//Therefore the eigenvectors are the same
			eigvec[0][0]=sqrt(2.);
			eigvec[0][1]=eigvec[0][0];
			eigvec[1][0]=eigvec[0][0];
			eigvec[1][1]=eigvec[0][0];
		}
		else if( feq(m[0][0],eigval[0]) ) {
			eigvec[0][0]=1.0;
			eigvec[0][1]=0.0;
			eigvec[1][0]=0.0;
			eigvec[1][1]=1.0;
		}
		else if( feq(m[1][1],eigval[0]) ){
			eigvec[0][0]=0.0;
			eigvec[0][1]=1.0;
			eigvec[0][0]=1.0;
			eigvec[0][1]=0.0;
		}
		else printf("Warning: 2D eigensolver failed.\n");
	}
}
void eigenvalues3x3( double **m,double eigval[] ) {
/*
    Find the three eigenvalues for m for a 3x3 matrix
    MUST BE SYMMETRIC
    http://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
    and Smith, Communications of the ACM 4 (4): 168, 1961.
*/
	int i,j;
	double B[_3D][_3D];
	double q,p1,p2,p,ip,r,phi;

	//Check if symmetric --- hopefully a waste of time and should be removed once confident that it is always symmetric
	for( i=0; i<_3D; i++ ) for( j=i+1; j<_3D; j++ ) if( fneq(m[i][j],m[j][i]) ) {
		printf( "Error: Matrix is not symmetric.\n" );
		ptens( m,_3D );
		exit(EXIT_FAILURE);
	}

	p1=m[0][1]*m[0][1] + m[0][2]*m[0][2] + m[1][2]*m[1][2];
	// Diagonal matrix
	if( feq(p1,0.0) ) {
		for( i=0; i<_3D; i++ ) eigval[i]=m[i][i];
		//Sort
		if( eigval[1]>eigval[0] ) {
			q=eigval[0];
			eigval[0]=eigval[1];
			eigval[1]=q;
		}
		if( eigval[2]>eigval[0] ) {
			q=eigval[0];
			eigval[0]=eigval[2];
			eigval[2]=q;
		}
		if( eigval[2]>eigval[1] ) {
			q=eigval[1];
			eigval[1]=eigval[2];
			eigval[2]=q;
		}
	}
	else {
		p2=0.;
		q = (m[0][0]+m[1][1]+m[2][2])/3.;
		for( i=0; i<_3D; i++ ) p2 +=  (m[i][i]-q) * (m[i][i]-q);
		p2 += 2.*p1;
		p = sqrt( p2/6.);
		ip = 1./p;
		for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) B[i][j] = m[i][j];
		for( i=0; i<_3D; i++ ) B[i][i] -= q;
		for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) B[i][j] *= ip;
		r = det3x3(B) * 0.5;

		if( r <= -1. ) phi = pi / 3.;
		else if( r >= 1. ) phi = 0.;
		else phi = acos(r) / 3.;

		eigval[0] = q + 2.*p*cos( phi );
		eigval[2] = q + 2.*p*cos( phi + (2.*pi/3.) );
		eigval[1] = 3.*q - eigval[0] - eigval[2];	// Cuz trace(A) = eig1 + eig2 + eig3
	}
}
void eigenvectors3x3( double **m,double eigval[],double eigvec[][_3D] ) {
/*
    Find the three eigenvectors (normalized) for m for a SYMMETRIC 3x3 matrix
    Uses the Cayley-Hamilton theorem from (http://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices).
    Don't worry about generalized eigenvector stuff for eigenvalue multiplicities greater than 1.
    Check if diagonal.
*/
	int row,col,i,k;
	double a,b;

	// Check if diagonal
	a=0.;
	for( col=0;col<_3D;col++ ) for( row=0;row<_3D;row++ ) if(col!=row) a+=m[row][col];
	a*=a;
	
	// TODO (but not really needed). This was a short cut, but the (largest) eigenvalue and corresponding eigenvector were not matched.
	//if(a<=TOL) for( col=0;col<_3D;col++ ) for( row=0;row<_3D;row++ ) eigvec[col][row]=m[row][col];
	//else {
		//Cayley-Hamilton gives eigenvector of k to be ANY column (as long as it's not zero). So pick the first non-zero column
	
	k=0;		//eigval[0]
	for( col=0;col<_3D;col++ ) {
		for( row=0;row<_3D;row++ ) eigvec[k][row]=0.;
		for( row=0;row<_3D;row++ ) for( i=0;i<_3D;i++ ) {
			a=m[row][i];
			if( i==row ) a-=eigval[1];
			b=m[i][col];
			if( i==col ) b-=eigval[2];
			eigvec[k][row]+=a*b;
		}
		// //Make sure didn't get a zero value --- If did then look for new solution; else stop
		if( !( feq(eigvec[k][0],0.0) && feq(eigvec[k][1],0.0) && feq(eigvec[k][2],0.0) ) ) break;
	}
	k=1;		//eigval[1]
	for( col=0;col<_3D;col++ ) {
		for( row=0;row<_3D;row++ ) eigvec[k][row]=0.;
		for( row=0;row<_3D;row++ ) for( i=0;i<_3D;i++ ) {
			a=m[row][i];
			if( i==row ) a-=eigval[0];
			b=m[i][col];
			if( i==col ) b-=eigval[2];
			eigvec[k][row]+=a*b;
		}
		//Make sure didn't get a zero value --- If did then look for new solution; else stop
		if( !( feq(eigvec[k][0],0.0) && feq(eigvec[k][1],0.0) && feq(eigvec[k][2],0.0) ) ) break;
	}
	k=2;		//eigval[2]
	for( col=0;col<_3D;col++ ) {
		for( row=0;row<_3D;row++ ) eigvec[k][row]=0.;
		for( row=0;row<_3D;row++ ) for( i=0;i<_3D;i++ ) {
			a=m[row][i];
			if( i==row ) a -= eigval[0];
			b=m[i][col];
			if( i==col ) b-=eigval[1];
			eigvec[k][row]+=a*b;
		}
		//Make sure didn't get a zero value --- If did then look for new solution; else stop
		if( !( feq(eigvec[k][0],0.0) && feq(eigvec[k][1],0.0) && feq(eigvec[k][2],0.0) ) ) break;
	}
	//}
	//Normalize
	for( k=0;k<_3D;k++ ) norm( eigvec[k],_3D );
}
void solveEigensystem( double **m,int dimension,double eigval[] ) {
/*
    Finds the eigenvalues and vectors of the real, symmetric matrix m by analytical methods
    The matrix m is lost.
    It becomes the eigenvectors: the kth column of m returns the normalized eigenvector corresponding to eigval[k].
*/
	int i,j;

	if( dimension==_2D ) {
		double eigvec[dimension][dimension];
		eigenvalues2x2( m,eigval );
		eigenvectors2x2( m,eigval,eigvec );
		for( i=0;i<dimension;i++ ) for( j=0;j<dimension;j++ ) m[i][j]=eigvec[i][j];
	}
	else if( dimension==_3D ) {
		double eigvec[dimension][dimension];
		eigenvalues3x3( m,eigval );
		eigenvectors3x3( m,eigval,eigvec );
		for( i=0;i<dimension;i++ ) for( j=0;j<dimension;j++ ) m[i][j]=eigvec[i][j];
	}
	else if( dimension==_1D ) {
		// In 1D, the idea of eigenvectors and values is a bit silly
		// But for programmatic reasons, to run in 1D we just set the vector to be it's only possible axis and value to be the value
		eigval[0]=m[0][0];
		m[0][0]=1.0;
	}
	else {
		printf( "Error: Solving the eigensystem for dimensions greater than 3 is not coded (DIM=%d).\n",dimension );
		exit(EXIT_FAILURE);
	}
}
double centredDeriv( double xM1,double xP1,double dt ) {
/*
    Find the derivative of x by a centred derivative
*/
	double deriv=0.5*(xP1-xM1)/dt;
	return deriv;
}
double forwardDeriv( double x0,double xP1,double dt ) {
/*
    Find the derivative of x by a centred derivative
*/
	double deriv=(xP1-x0)/dt;
	return deriv;
}
double backwardDeriv( double x0,double xM1,double dt ) {
/*
    Find the derivative of x by a centred derivative
*/
	double deriv=(x0-xM1)/dt;
	return deriv;
}
double simps( double F[],double dx,int n ) {
/*
    Find the integral of a discrete function with equal steps in x
*/
	int i,halfN;
	double t1=0.0,t2=0.0;
	halfN=n/2;
	for( i=1; i<halfN-1; i++ ) t1+=F[2*i];
	t1*=2.;
	for( i=1; i<halfN; i++ ) t2+=F[2*i-1];
	t2*=4.;
	return dx*(F[0] + t1 + t2 + F[n-1])/3.;
}

double stdNum( cell ***CL,int GPOP,int XYZ[3],int XYZ_P1[3] ) {
/*
    Find the standard deviation of the number of particles in each cell
    NOTICE: previously s1 was average number per cell but this was erroneous. Should just be sum
*/
	int a,b,c;
	double std,NC,s2,s1;

	NC=(double) XYZ[0]*XYZ[1]*XYZ[2];
	//Average number density (cell size always a=1)
	// s1=((double) GPOP)/NC;
	// Sum
	s1=(double) GPOP;
	s2=0.;

	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		s2+= CL[a][b][c].POP*CL[a][b][c].POP;
	}
	std=sqrt( (NC*s2-s1*s1)/(NC*(NC-1.)) );
	return std;
}
void rodriguesRotation( double vec[],double rotAx[],double theta ) {
	/*
	 This routine rotates one vector (vec) about an axis of rotation vector (rotAx) by an angle theta and writes over the vector
	 For some reason the rotation appears to shrink vec's magnitude slightly. Therefore rescale.
	 NOTICE: rotAx MUST be a UNIT vector so immediately normalized
	 NOTICE: that EVEN if this is 2D it MUST be 3D because the cross product cp will be in the perpendicular direction
	 */
	int i;
	double cp[_3D],dp=0.0;
	double old,new;

	norm( rotAx,_3D );
	old=length( vec,_3D );
	crossprod( rotAx,vec,cp );
	dp=dotprod( rotAx,vec,_3D );
	for( i=0; i<_3D; i++ ) vec[i] = vec[i]*cos(theta) + cp[i]*sin(theta) + rotAx[i]*dp*(1.0-cos(theta));
	//Unfortunately, the rotation seems to shrink vec
	new=length( vec,_3D );
	for( i=0; i<_3D; i++ ) vec[i]*=(old/new);
}
void setRotMatrix3D( double M[][3],double angx,double angy,double angz ) {
/*
   Just sets a rotation matrix based on angles about the cartesian axes
*/
	double cosx,sinx,cosy,siny,cosz,sinz;
	cosx=cos(angx);
	cosy=cos(angy);
	cosz=cos(angz);
	sinx=sin(angx);
	siny=sin(angy);
	sinz=sin(angz);

	M[0][0] = cosy*cosz;
	M[0][1] = cosx*sinz + sinx*siny*cosz;
	M[0][2] = sinx*sinz - cosx*siny*cosz;
	M[1][0] = -cosy*sinz;
	M[1][1] = cosx*cosz - sinx*siny*sinz;
	M[1][2] = sinx*cosz + cosx*siny*sinz;
	M[2][0] = siny;
	M[2][1] = -sinx*cosy;
	M[2][2] = cosx*cosy;
}
void setRotMatrix2D( double M[][3],double angz ) {
/*
   Just sets a rotation matrix based on angles about the cartesian z axis (uses as input a 3x3 array regardless of dimensionality)
*/
		M[0][0] = cos(angz);
	M[0][1] = -sin(angz);
	M[1][0] = -M[0][1];
	M[1][1] = M[0][0];
}
void skewSymmetricCrossProductMatrix( double *v,double result[][3] ) {
/*
    Generate the skew-symmetric cross-product matrix needed for the rotation in findRotationMatrix()
*/
	result[0][0]=0.;
	result[1][1]=0.;
	result[2][2]=0.;
	result[0][1]=-v[2];
	result[1][0]=v[2];
	result[0][2]=v[1];
	result[2][0]=-v[1];
	result[1][2]=-v[0];
	result[2][1]=v[0];
}
void rotationMatrix( double rotMat[][3],double vx[][3],double c,double s ) {
/*
    Find the rotation matrix for the subroutine findRotationMatrix()
*/
	double unity[_3D][_3D],vx2[_3D][_3D];
	int i,j;

	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) {
		unity[i][j]=0.;
		rotMat[i][j]=0.;
	}
	for( i=0; i<_3D; i++ ) unity[i][i]=1.;
	dotprodMatMat( vx,vx,vx2,_3D );
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) rotMat[i][j] = unity[i][j]+vx[i][j]+vx2[i][j]*(1.-c)/s/s;
}
void findRotationMatrix( double rotMat[][3],double *original,double *final ) {
/*
    Find the rotation matrix necesary to rotation the original vector parallel to the final vector
*/
	double a[_3D],b[_3D],v[_3D];
	double vx[_3D][_3D];
	double s,c;

	normCopy(original,a,_3D);
	normCopy(final,b,_3D);
	crossprod( a, b, v );
	s = length( v,_3D );
	c = dotprod( a,b,_3D );
	skewSymmetricCrossProductMatrix( v,vx );
	rotationMatrix( rotMat,vx,c,s );
}
void dirdirCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
/*
    Find the spacial autocorrelation function of director
*/
	int a,b,c,d;
	int aa,bb,cc;
	int cnt[maxXYZ];

	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] = 0.;
		cnt[d] = 0;
	}

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
			d=round( sqrt( (double)d ) );
			if( d>=maxXYZ ) printf( "maxXYZ=%d\td=%d\n",maxXYZ,d );
			cnt[d] += 1;
			avCorr[d] += fabs( dotprod( CL[a][b][c].DIR, CL[aa][bb][cc].DIR, dimension ) );
		}
	}
	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
		avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1); ;
	}
	//Transformation to make the correlation function go from unity to zero
	//corr0=1 and corrINF=2/3 ideally
	//I think now that I should leave this for post analysis!!!
	//corrINF=avCorr[maxXYZ/2];
	//for( d=0; d<maxXYZ; d++ ) avCorr[d] = (avCorr[d]-corrINF)/(corr0-corrINF);
}
void densdensCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
/*
    Find the spacial autocorrelation function of director
*/
	int a,b,c,d;
	int aa,bb,cc;
	int cnt[maxXYZ];

	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] = 0.;
		cnt[d] = 0;
	}

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
			d=round( sqrt( (double)d ) );
			if( d>=maxXYZ ) printf( "maxXYZ=%d\td=%d\n",maxXYZ,d );
			cnt[d] += 1;
			avCorr[d] += (double) (CL[a][b][c].POP*CL[aa][bb][cc].POP);
		}
	}
	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
		avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1);
	}
}
void orderorderCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
/*
    Find the spacial autocorrelation function of scalar order parameter
*/
	int a,b,c,d;
	int aa,bb,cc;
	int cnt[maxXYZ];

	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] = 0.;
		cnt[d] = 0;
	}

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
			d=round( sqrt( (double)d ) );
			if( d>=maxXYZ ) printf( "maxXYZ=%d\td=%d\n",maxXYZ,d );
			cnt[d] += 1;
			avCorr[d] += CL[a][b][c].S * CL[aa][bb][cc].S;
		}
	}
	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
		avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1);
	}
}
// void phiphiCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
// /*
//     Find the spacial autocorrelation function of binary phase (phi)
// */
// 	int a,b,c,d;
// 	int aa,bb,cc;
// 	int cnt[maxXYZ];
//
// 	for( d=0; d<maxXYZ; d++ ) {
// 		avCorr[d] = 0.;
// 		cnt[d] = 0;
// 	}
//
// 	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
// 		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
// 			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
// 			d=round( sqrt( (double)d ) );
// 			if( d>=maxXYZ ) printf( "maxXYZ=%d\td=%d\n",maxXYZ,d );
// 			cnt[d] += 1;
// 			avCorr[d] += CL[a][b][c].PHI * CL[aa][bb][cc].PHI;
// 		}
// 	}
// 	for( d=0; d<maxXYZ; d++ ) {
// 		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
// 		avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1);
// 	}
// }
void velvelNormedCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
/*
    Find the spacial autocorrelation function of velocity
*/
	int a,b,c,d;
	int aa,bb,cc;
	int cnt[maxXYZ];

	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] = 0.;
		cnt[d] = 0;
	}

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
			d=(int) sqrt( (double)d );
			if( d>=maxXYZ ) printf( "Warning: maxXYZ=%d\td=%d\n",maxXYZ,d );
			cnt[d] += 1;
			avCorr[d] += dotprod( CL[a][b][c].VCM, CL[aa][bb][cc].VCM, dimension );
		}
	}
	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
		avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1);
	}
}
void vortvortNormedCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
/*
    Find the spacial autocorrelation function of velocity
*/
	int a,b,c,d;
	int aa,bb,cc;
	double w1[_3D],w2[_3D];
	int cnt[maxXYZ];

	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] = 0.;
		cnt[d] = 0;
	}

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
	  w1[0]=(CL[a][b][c].E[2][1] - CL[a][b][c].E[1][2]);
		w1[1]=(CL[a][b][c].E[0][2] - CL[a][b][c].E[2][0]);
		w1[2]=(CL[a][b][c].E[1][0] - CL[a][b][c].E[0][1]);
		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
			d=round( sqrt( (double)d ) );
			if( d>=maxXYZ ) printf( "maxXYZ=%d\td=%d\n",maxXYZ,d );
			w2[0]=(CL[aa][bb][cc].E[2][1] - CL[aa][bb][cc].E[1][2]);
			w2[1]=(CL[aa][bb][cc].E[0][2] - CL[aa][bb][cc].E[2][0]);
			w2[2]=(CL[aa][bb][cc].E[1][0] - CL[aa][bb][cc].E[0][1]);
			cnt[d] += 1;
			avCorr[d] += dotprod( w1,w2, dimension );
		}
	}
	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
		avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1);
	}
}
void velvelCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
/*
    Find the spacial autocorrelation function of velocity
*/
	int a,b,c,d;
	int aa,bb,cc;
	int cnt[maxXYZ];

	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] = 0.;
		cnt[d] = 0;
	}

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
			d=round( sqrt( (double)d ) );
			if( d>=maxXYZ ) printf( "Warning: maxXYZ=%d\td=%d\n",maxXYZ,d );
			cnt[d] += 1;
			avCorr[d] += dotprod( CL[a][b][c].VCM, CL[aa][bb][cc].VCM, dimension );
		}
	}
	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
		// Velocity is NOT normalized.
		// avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1);
	}
}
void vortvortCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension ) {
/*
    Find the spacial autocorrelation function of velocity
*/
	int a,b,c,d;
	int aa,bb,cc;
	double w1[_3D],w2[_3D];
	int cnt[maxXYZ];

	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] = 0.;
		cnt[d] = 0;
	}

	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) {
	  w1[0]=(CL[a][b][c].E[2][1] - CL[a][b][c].E[1][2]);
		w1[1]=(CL[a][b][c].E[0][2] - CL[a][b][c].E[2][0]);
		w1[2]=(CL[a][b][c].E[1][0] - CL[a][b][c].E[0][1]);
		for( aa=0; aa<XYZ[0]; aa++ ) for( bb=0; bb<XYZ[1]; bb++ ) for( cc=0; cc<XYZ[2]; cc++ ) {
			d = (aa-a)*(aa-a) + (bb-b)*(bb-b) + (cc-c)*(cc-c);
			d=round( sqrt( (double)d ) );
			if( d>=maxXYZ ) printf( "maxXYZ=%d\td=%d\n",maxXYZ,d );
			w2[0]=(CL[aa][bb][cc].E[2][1] - CL[aa][bb][cc].E[1][2]);
			w2[1]=(CL[aa][bb][cc].E[0][2] - CL[aa][bb][cc].E[2][0]);
			w2[2]=(CL[aa][bb][cc].E[1][0] - CL[aa][bb][cc].E[0][1]);
			cnt[d] += 1;
			avCorr[d] += dotprod( w1,w2, dimension );
		}
	}
	for( d=0; d<maxXYZ; d++ ) {
		avCorr[d] /= (double) (cnt[d]?(cnt[d]):1);
		// Vorticity is NOT normalized.
		// avCorr[d] /= (double) (avCorr[0]?(avCorr[0]):1);
	}
}
void normCorr( double *corr,int maxXYZ ) {
/*
    Normalize an unnormalized correlation function
*/
	int i;
	double corr0=0.;

	corr0=corr[0];
	for( i=0; i<maxXYZ; i++ ) corr[i] = corr[i]/corr0;
}
void FT_spherical( double *f,double *F,double *rad,int n,int dimension ) {
/*
    Transform a spherically symmetric function into its Fourier transform
*/
	int r,k;
	double integrand[n],waveNum;
	double dr,pi2;

	dr=(rad[n-1]-rad[0])/((double)n);
	pi2=2.0*pi;
	if( dimension==_3D ) {
		for( k=0; k<n; k++ ) {
			waveNum=pi2/rad[k];
			for( r=0; r<n; r++ ) integrand[r]=f[r]*rad[r]*sin( waveNum*rad[r] );
			F[k] = 2.0*pi2*simps( integrand,dr,n )/waveNum;
		}
	}
	else if( dimension==_2D ) {
		for( k=0; k<n; k++ ) {
			waveNum=pi2/rad[k];
			// for( r=0; r<n; r++ ) integrand[r]=f[r]*sin( waveNum*rad[r] );
			// F[k] = 2.0*simps( integrand,dr,n )/waveNum;
			for( r=0; r<n; r++ ) integrand[r]=f[r]*rad[r]*j0( waveNum*rad[r] );
			F[k] = pi2*simps( integrand,dr,n );
		}
	}
	else {
		printf( "Error: Fourier transform failed." );
		exit(EXIT_FAILURE);
	}
}
void FTspectrum( double *corr,double *spect,int maxXYZ,int dimension ) {
/*
    Transform a correlation function into a spectrum
*/
	int i;
	double rad[maxXYZ],waveNum,pi2;

	pi2=2.0*pi;
	for( i=0; i<maxXYZ; i++ ) rad[i] = (double) i;
	FT_spherical( corr,spect,rad,maxXYZ,dimension );
	if( dimension==_3D ) for( i=1; i<maxXYZ; i++ ) {
		waveNum=pi2/rad[i];
		spect[i] *= waveNum*waveNum*pi2;
	}
	else if( dimension==_2D ) for( i=1; i<maxXYZ; i++ ) {
		waveNum=pi2/rad[i];
		spect[i] *= waveNum*pi;
	}
	else {
		printf( "Error: Spectrum failed." );
		exit(EXIT_FAILURE);
	}
}
int checkNAN_vec( double vec[],int dimension ) {
	/*
	    Check that no component of a vector is NAN or INF
	*/
	int d=0,flag=0;
	for( d=0; d<dimension; d++ ) {
		if(isnan( vec[d] )) flag=1;
		else if(isinf( vec[d] )) flag=1;
	}
	return flag;
}
void checkNAN_Q( cell ***CL,int XYZ_P1[3],int pauseFlag,int dimension ) {
	/*
	    Check that no position values are NANs or INFs
	*/
	int i,j,k,flag,cnt=0;
	particleMPC *cp;	//Pointer to current item in list
	//Search each cell for particleMPCs that have left the cell
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		if( CL[i][j][k].pp != NULL ) {
			cp = CL[i][j][k].pp;
			while(cp != NULL) {
				flag = checkNAN_vec( cp->Q,dimension );
				if( flag != 0 ) {
					printf( "\tWarning: Bad position found. Particle in cell [%d,%d,%d] with POP=%d. Q=",i,j,k,CL[i][j][k].POP );
					pvec( cp->Q,dimension );
					cnt++;
				}
				//Increment link in list
				cp = cp->next;
			}
		}
	}
	if( cnt>0 && pauseFlag ) wait4u();
}
void checkNAN_V( cell ***CL,int XYZ_P1[3],int pauseFlag,int dimension ) {
	/*
	    Check that no velocity values are NANs or INFs
	*/
	int i,j,k,flag,cnt=0;
	particleMPC *cp;	//Pointer to current item in list

	//Search each cell for particleMPCs that have left the cell
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		if( CL[i][j][k].pp != NULL ) {
			cp = CL[i][j][k].pp;
			while(cp != NULL) {
				flag = checkNAN_vec( cp->V,dimension );
				if( flag != 0 ) {
					printf( "\tWarning: Bad velocity found. Particle in cell [%d,%d,%d] with POP=%d. V=",i,j,k,CL[i][j][k].POP );
					pvec( cp->V,dimension );
					cnt++;
				}
				//Increment link in list
				cp = cp->next;
			}
		}
	}
	if( cnt>0 && pauseFlag ) wait4u();
}
