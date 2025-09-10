#ifndef SWIMMERS_H
#define SWIMMERS_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

void readswimmers( char fpath[],specSwimmer *specS,swimmer **sw );
void setswimmers( specSwimmer *SS,swimmer *swimmers,bc WALL[],int stepsMD,double dt );
void swimout( FILE *fout,swimmer swimmers[],double t );
void swimoriout( FILE *fout,swimmer swimmers[],double t );
void runtumbleout( FILE *fout,swimmer *sw );
void swimmerPBC_dr (double *dr );
void swimmerOri( double n[],swimmer *sw );
double swimmerWCA( double r,double eps );
double swimmerFENE( double r,double k );
double swimmerHookean( double r,double k );
double swimmerSpring6( double r,double k );
void swimmerVerlet_nonInteracting( specSwimmer SS,swimmer *s,double dt,int springType,bc WALL[],int i);
void swimmerVerlet_all( specSwimmer SS,swimmer swimmers[],double dt,int springType, bc WALL[]);
void smonoDist( double r[],double *dr,smono m1, smono m2 );
double smonoForceMag_sameSwimmer( double dr,specSwimmer SS,swimmer *s,int springType );
void smonoForce_sameSwimmer( double a[],specSwimmer SS,swimmer *s,int springType );
double smonoForceMag_differentSwimmers( double dr,specSwimmer SS );
void smonoForce_differentSwimmers( double a[],specSwimmer SS,smono s1,smono s2 );
void integrateSwimmers( specSwimmer SS,swimmer swimmers[],bc WALL[],int stepsMD,double timeStep,double MAG[],int springType );
void swimmerDipole( specSwimmer SS,swimmer swimmers[],cell ***CL,spec SP[],double timeStep,particleMPC *SRDparticles,bc WALL[],simptr simMD,int noHI2 );
void swimmerForceDipole( specSwimmer SS,swimmer *sw,cell ***CL,spec SP[],double timeStep,int noHI2 );
void swimmerRotletDipole( specSwimmer SS,swimmer *sw,cell ***CL,spec SP[],double timeStep,int noHI2 );
void swimmerMagTorque( specSwimmer SS,swimmer *sw,double dt,double MAG[] );
void allSwimmersMagTorque( specSwimmer SS,swimmer swimmers[],double timeStep,int stepsMD,double MAG[] );
void runTumbleDynamics( specSwimmer *SS,swimmer swimmers[],bc WALL[],int stepsMD,double MAG[],double dt,int RTOUT,FILE *fruntumble );

void bininSwimmers( specSwimmer SS,swimmer swimmers[],cell ***CL );
void binSwimmers( cell ***CL,int shifted );
void addlinkSwimmer( cell *CL,smono *s );
void removelinkSwimmer( smono *current,cell *CL );

//Swimmer-BC interactions kept in mdbc.c

#endif
