# include <math.h>
# include <sys/time.h>
# include <stdio.h>
# include <unistd.h>
# include <stdint.h>

# include "../headers/definitions.h"
# include "../headers/SRDclss.h"
# include "../headers/mtools.h"
# include "../headers/pout.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************* MERSENNE TWISTER *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
static unsigned long mt[NN];	//Mersenne twister variable: the array for the state vector
static int mti=NN+1;		//Mersenne twister variable: mti==NN+1 means mt[NN] is not initialized

unsigned long MT_RandomSeedSRD (unsigned long seed)
{
// STOLEN FROM FRED!

    // Get a seed from time*pid if seed=0
    // 	if (!seed) seed = time(0)*getpid();
    // if (!seed) seed = time(0);
    struct timeval tv;
    gettimeofday(&tv, NULL); // Get the time to use the microseconds as an "random" seed
    if (!seed) seed = tv.tv_usec+getpid();

    // Initialize mersenne twister array
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<NN; mti++) {
        mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
    }

    return (seed);
}

void MT_init_genrand(unsigned long s){
/*
   Mersenne twister
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
   Initializes mt[NN] with a seed
*/
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<NN; mti++) {
        mt[mti] =
                (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}
void MT_init_by_array(unsigned long init_key[], int key_length){
/*
   Mersenne twister
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
   Iinitialize by an array with array-length
   init_key is the array for initializing keys
   key_length is its length
*/
    int i, j, k;

    struct timeval tv;
    gettimeofday(&tv, NULL); // Get the time to use the microseconds as an "random" seed
    MT_init_genrand(tv.tv_usec);
    //init_genrand(19650218UL);

    i=1; j=0;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
                + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
                - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}
unsigned long MT_genrand_int32(void){
/*
   Mersenne twister
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
   Generates a random number on [0,0xffffffff]-interval
*/
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    struct timeval tv;

    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    if (mti >= NN) { /* generate NN words at one time */
        int kk;

        if (mti == NN+1) {  /* if init_genrand() has not been called, */
            gettimeofday(&tv, NULL); // Get the time to use the microseconds as an "random" seed
            MT_init_genrand(tv.tv_usec);
            //init_genrand(5489UL); /* a default initial seed is used */
        }
        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
    y = mt[mti++];
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** xoshiro128++ ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/*  Taken from: https://prng.di.unimi.it/xoshiro128plusplus.c

    Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

    To the extent possible under law, the author has dedicated all copyright
    and related and neighboring rights to this software to the public domain
    worldwide. This software is distributed without any warranty.

    See <http://creativecommons.org/publicdomain/zero/1.0/>. */

static unsigned long X_state[4]; // RNG state
int X_seeded = 0; // flag to show whether this has been seeded or not.

static inline unsigned long X_rotl(const long int x, int k) {
    /* Rotate x left by k bits */
    return (x << k) | (x >> (32 - k));
}

void X_init_genrand(unsigned long s) {
    /*
     * Initialize the RNG state with the seed using SplitMix64
     *
     * Note that unlike MT: We need to prepare 4 pseudo-random values to initialise the RNG state
     * Hence, we initialise a temporary instance of SplitMix64, using the existing seed, to generate the initial state
    */
    int i; // counting variable

    /* SplitMix64 code taken from: https://github.com/svaarala/duktape/blob/master/misc/splitmix64.c
     * Written in 2015 by Sebastiano Vigna (vigna@acm.org)
        To the extent possible under law, the author has dedicated all copyright
        and related and neighboring rights to this software to the public domain
        worldwide. This software is distributed without any warranty.
        See <http://creativecommons.org/publicdomain/zero/1.0/>. */

    unsigned long sm_state = s; // splitmix64 state

    for (i = 0; i < 4; i++) {
        // generate the next value in splitmix
        unsigned long z = (sm_state += UINT64_C(0x9E3779B97F4A7C15));
        z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);

        X_state[i] = z ^ (z >> 31); // output value from splitmix
    }

    X_seeded = 1; // mark as seeded
}

unsigned long X_RandomSeedSRD (unsigned long seed) {
    /*
     * Perform the same seeding as in the MT
    */
    struct timeval tv;
    gettimeofday(&tv, NULL); // Get the time to use the microseconds as an "random" seed
    if (!seed) seed = tv.tv_usec+getpid();

    if (X_seeded == 0) {
        X_init_genrand(seed);
    }

    return seed;
}

unsigned long X_genrand_int32(void) {
    /*
     * Performs the next() step from the xoshiro128++ algorithm, generating a random integer
    */
    if (X_seeded == 0) { // ensure seed is properly set, if not then seed with a random value
        X_RandomSeedSRD(0);
    }

    const unsigned long result = X_rotl(X_state[0] + X_state[3], 7) + X_state[0];
    const unsigned long t = X_state[1] << 9;

    X_state[2] ^= X_state[0];
    X_state[3] ^= X_state[1];
    X_state[1] ^= X_state[2];
    X_state[0] ^= X_state[3];

    X_state[2] ^= t;

    X_state[3] = X_rotl(X_state[3], 11);

    return result;
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ interface methods *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

unsigned long RandomSeedSRD (unsigned long seed)
{
    /*
     * Check if a random seed is required. If so then generates one pseudo-randomly.
     * Then initialise the random number generators using a pseudo-random seed
     */
#ifdef RNG_MERSENNE
    return MT_RandomSeedSRD(seed);
#else
    return X_RandomSeedSRD(seed);
#endif
}

void init_genrand(unsigned long s){
    /*
     * Initialise the random number generators WITHOUT creating a random seed.
     */
#ifdef RNG_MERSENNE
    MT_init_genrand(s);
#else
    X_init_genrand(s);
#endif
}

unsigned long genrand_int32(void){
    /*
     * Base RNG method. Returns a random unsigned long.
     */
#ifdef RNG_MERSENNE
    return MT_genrand_int32();
#else
    return X_genrand_int32();
#endif
}

long genrand_int31(void){
/*
   Mersenne twister
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
   Generates a random number on [0,0x7fffffff]-interval
*/
    return (long)(genrand_int32()>>1);
}
double genrand_real(void){
/*
   Maps the generated random number onto a [0,1) interval
   divisor controls the precision of the random number. The higher the divisor, the more precise the double.
*/
    const double divisor = 4294967296.0; // 2^32 by default
    return (genrand_int32() % (unsigned long) divisor) * (1.0/divisor); // modulo ensures this is always bounded
}
double genrand_pmOne(void){
	/*
	 Randomly generates +1.0 or -1.0 with 50% probability
	*/
	double rand = genrand_real();
	if(rand<=0.5) return -1.0;
	else return 1.0;
}
void genrand_coneNP( double vec[],double theta,int dimension ) {
	/*
	 Generate a random, uniformly distributed normalized vector for a direction within cone around "north pole"/x-axis
	*/
	double z,phi,ct;

	ct=cos(theta);
	if( dimension==_3D ) {
		// Sample z on [cos(theta),1] and phi on [0,2*pi]
		z=genrand_real()*(1.0-ct)+ct;
		phi=genrand_real()*2.0*pi;
		vec[0]=z;
		vec[1]=sqrt(1.0-z*z)*cos(phi);
		vec[2]=sqrt(1.0-z*z)*sin(phi);
	}
	else if( dimension==_2D ) {
		// Sample angle homogeneously
		phi=theta*(1.-2.*genrand_real());
		vec[0]=cos(phi);
		vec[1]=sin(phi);
	}
    else if( dimension==_1D ) {
		// A bit of a funny definition of a cone in 1D.
        // Checks if cos(phi) 2D would be parallel or antiparallel
		phi=theta*0.5*pi*(1.-2.*genrand_real());
        if(phi>1.0 || phi<-1.0) vec[0]=-1.0;
        else vec[0]=1.0;
	}
	else printf("Warning: genrand_coneNP() only programmed for DIM={3,2,1}, not DIM=%d\n",dimension);
}
void genrand_cone( double axis[],double vecOut[],double theta,int dimension ) {
	/*
	 Generate a random, uniformly distributed normalized vector for a direction within cone around axis
	 First generates around the x-axis usng genrand_coneNP() then rotates
	 Even for dimension=2 the vectors must be 3D
	*/
	double rotAx[_3D],randVec[_3D]={0.0},xaxis[_3D]={0.0};
	double angle;
	int i;
	xaxis[0]=1.0;

	//Generate random vec about z-axis/north pole
	genrand_coneNP( randVec,theta,dimension );
	//Rotate north pole to align with the axis direction
	//Find the axis that the x-axis must be rotated about
	crossprod( axis,xaxis,rotAx );
	angle = absAngle( axis,xaxis,_3D );
	//Even if 2D, this must be 3D cuz then rotAX will be in 3rd dimension.
	//This rotation is why 3D vectors were required
	rodriguesRotation( randVec,rotAx,angle );
	for( i=0; i<dimension; i++ ) vecOut[i]=randVec[i];
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************** GAUSSIAN DIST ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
float genrand_gauss( void ) {
/*
   Box-Muller transformation to turn a uniform random number
   between 0-1 into a gaussian distribution of mean 0 and
   StDev of 1. From
   http://www.taygeta.com/random/gaussian.html
*/
	float x1,x2,w,y1;
// 	float y2;
	do{
		x1 = 2. * genrand_real() - 1.;
		x2 = 2. * genrand_real() - 1.;
		w = x1*x1+x2*x2;
	} while( w >= 1. );
	w = sqrt( (-2. * log( w )) / w );
	y1 = x1*w;
	return y1;
}
double genrand_gaussMB(double KBT,double M) {
/*
  Scale the normal distribution into a Gaussian distribution
  of variance KBT/M, which produces the Maxwell-Boltzmann distribution
*/
	double sigma=sqrt(KBT/M);
	return sigma*genrand_gauss();
}
double genrand_gaussGen(double mu,double sigma) {
/*
  Scale the normal distribution into a Gaussian distribution
  with a mean mu and standard deviation sigma
*/
	return sigma*genrand_gauss()+mu;
}
double genrand_exp(double lambda) {
/*
  An exponetial distribution using Box-Muller. Here lambda = 1/rate ie the mean
*/
	return -lambda*log( genrand_real() );
}
int genrand_poisson(double lambda) {
/*
  Poisson distributed random numbers by Junhao, based on Knuth
	From
	https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
*/
	double L,eSTEP,STEP=50.0;	//STEP is chosen for double precision
	double r,myExp=M_E,p=1.0;
	int k=0;

	// For "small" values of lambda
	if( lambda<STEP ) {
		L=exp(-lambda);
		do{
			k+=1;
			r = genrand_real();
			p*=r;
		}while( p>L );
	}
	// For "large" values of lambda the exponent is not well known
	else if( lambda<10*STEP ) {
		L=lambda;
		eSTEP=exp(STEP);
		do{
			k+=1;
			r = genrand_real();
			p*=r;
			if( p<myExp && L>0.0 ) {
				if( L>STEP ){
					p*=eSTEP;
					L-=STEP;
				}
				else{
					p*=L;
					L=-1.0;
				}
			}
		}while( p>1.0 );
	}
	// For large values of lambda just use a normal distribution approximation
	else {
		L=genrand_gaussGen( lambda,sqrt(lambda) );
		k=1+(int)L;
	}
	return k-1;
}
float genrand_rayleigh( float std ) {
/*
   Inverse transformation to turn a uniform random number
   between 0-1 into a xe^{-x^2} distribution. This used to
   require a numeric solution (Newton's method - iterated until
   percent diff between iterations < PRCNT).
   BUT there is actually an easy transformation.
*/
	float r,lnr;	// Input uniform random number
	float x;	// Distributed random number

	// Uniform Random number
	do{
		r = genrand_real();
	} while ( feq(r,0.0) );
	lnr = log( r );
	x = std * sqrt( -2.*lnr );

	return x;
}
void genrand_sphere( double vec[],int dimension ) {
	/*
	 Generate a random, normalized vector uniformly distributed on a sphere
	*/
// 	int i;
// 	for( i=0; i<dimension; i++ ) vec[i]=genrand_gauss();
// 	norm( vec,dimension );

	double u,theta,t1;

	if( dimension==_3D ) {
		u=genrand_pmOne()*genrand_real();
		theta=genrand_real()*pi;
		t1=sqrt(1.-u*u);
		vec[0] = t1*cos(theta);
		vec[1] = t1*sin(theta);
		vec[2] = u;
	}
	else if( dimension==_2D ) {
		theta=genrand_real()*pi;
		vec[0]=cos(theta);
		vec[1]=sin(theta);
	}
	else if( dimension==_1D ) vec[0]=genrand_pmOne();
}
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************* RANDOM VECTORS ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
double *ranshift( double *v,int doShift,int dimension ) {
/*
   Produces a vectors with uniformly distributed
   random orientations, where no component is
   greater than unity.
	 Each component of the random shift is [0,1)
*/
	double x,y,z;
	if( doShift ) {
		if( dimension >= _3D ) z = genrand_real();
		else z = 0.0;
		if( dimension >= _2D ) y = genrand_real();
		else y=0.0;
		x = genrand_real();
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}
	else {
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
	}
	return v;
}
double *ranvec3D( double *v ) {
/*
   Produces unit vectors in 3 dimensions with
   uniformly distributed random orientations, using
   a standard rejection method [rap95].
*/
	double x,y,s,r;
	do{
		r = genrand_real();
		x = 2.*r - 1.;
		r = genrand_real();
		y = 2.*r - 1.;
		s = x*x + y*y;
	}while( s > 1. );
	v[2] = 1. - 2.*s;
	s = 2.*sqrt( 1.-s );
	v[0] = s * x;
	v[1] = s * y;
	return v;
}
double *ranvec2D( double *v ) {
/*
   Produces unit vectors in 2 dimensions with
   uniformly distributed random orientations, by
   randomly generating x and solving for y by
   the Pythagorean theorem
*/
	double x,y;
	x = genrand_real();
	y = sqrt( 1.-x*x );
	v[0] = x;
	v[1] = y;
	return v;
}
double *ranvec( double *v,int dimension ) {
/*
   Produces a random vector in 2 or 3D
*/
	if( dimension==_3D ) ranvec3D( v );
	if( dimension==_2D ) ranvec2D( v );
	return v;
}
int rand_particle( int POP ) {
/*
   Randomly picks one of the MPC particles
*/
	return (int)( genrand_real()*(double)POP );
}
