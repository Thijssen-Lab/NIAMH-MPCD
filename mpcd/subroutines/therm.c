# include<stdio.h>
# include<math.h>
# include<stdlib.h>

# include "../headers/SRDclss.h"
# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/mpc.h"
# include "../headers/rand.h"
# include "../headers/mtools.h"
# include "../headers/therm.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** THERMO ***************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
double thermostat( double KBT,double KBTNOW,double t,double RELAX,int TSTECH ) {
/*
   This routine calculates the velocity scaling factor,g
*/
	double g;
	//No thermostat used
	if( TSTECH == NOTHERM ) g = 1.;
	//Velocity rescaling
	else if( TSTECH==VSC ) g = sqrt( KBT/KBTNOW );
	//Berendsen rescaling
	else if( TSTECH==BEREND ) g = sqrt( 1.+(KBT/KBTNOW-1.)*t/RELAX );
	else if( TSTECH==MAXV ) {
		g = sqrt( 1.+(KBT/KBTNOW-1.)*t/RELAX );	//Same as BEREND
// 		g = sqrt( 1.+(KBT/KBTNOW-1.)*t/RELAX );
// 		g = 1.+(KBT/KBTNOW-1.)*t/RELAX;
// 		g = KBT/KBTNOW;
// 		g = 1.0;
// 		if( g>10.0 ) printf( "%lf\n",g );
	}
	else{
			printf( "Error:\tThermostat unacceptable.\n" );
			exit( 1 );
		}
	return g;
}
void scaleT( double KBT,double KBTNOW,double t,double RELAX,double VEL[],double VELNOW[],int TSTECH,spec SP[],int LC,bc *WALL,particleMPC *p,cell ***CL ) {
/*
   This routine calls the thermostat and scales the velocities
*/
	int i,j,k;
	double TSC;	//The velocity scaling factor (commonly lambda)

	//Local thermostat
	if( TSTECH==HEYES ) {
		bin( CL,SP,WALL,KBT,LC,0 );
		//bin( CL );
		for( i=0;i<XYZ_P1[0];i++ ) for( j=0;j<XYZ_P1[1];j++ ) for( k=0;k<XYZ_P1[2];k++ ) heyes_cell( CL[i][j][k],KBT,KBTNOW,RELAX );
	}
// 	Global thermostats
	else if( TSTECH!=NOTHERM ) {
		if( TSTECH==MAXV ) {
			//Relaxes to a given velocity (actually it relaxes to a given kinetic energy rather than thermal energy (or shifted thermal energy))
			TSC = thermostat( KBT,KBTNOW,t,RELAX,TSTECH );
			for( j=0; j<DIM; j++ ) {
				if( fneq(VEL[j],0.0) ) {
// 					TSC = thermostat( VEL[j],VELNOW[j],t,RELAX,TSTECH );
					for( i=0; i<GPOP; i++ ) (p+i)->V[j] = VEL[j] + TSC * ((p+i)->V[j] -VEL[j]);
					for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
						(WALL+i)->V[j] *= TSC;
						(WALL+i)->L[j] *= TSC;
					}
				}
			}
		}
		else {
			//Relaxes to zero
			TSC = thermostat( KBT,KBTNOW,t,RELAX,TSTECH );
			//Scale the velocities
			for( i=0; i<GPOP; i++ ) for( j=0; j<DIM; j++ ) (p+i)->V[j] *= TSC;
			for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) for( j=0; j<DIM; j++ ) {
				(WALL+i)->V[j] *= TSC;
				(WALL+i)->L[j] *= TSC;
			}
		}
	}
}
void heyes_cell( cell CL,double KBT,double KBTNOW,double RELAX ) {
	int d;
	double sc_fctr;	//Scaling factor
	double prob;		//Probability of rescaling
	particleMPC *pMPC;	//Temporary pointer to MPC particles

	sc_fctr= 1.+RELAX*genrand_real( );
	if( genrand_real( ) < 0.5 ) sc_fctr = 1./sc_fctr;
	prob = 0;

	//Create probability of rescaling
	if( CL.pp != NULL ) {
		pMPC = CL.pp;
		while(pMPC != NULL) {
			for( d=0;d<DIM;d++ ) prob += (pMPC->V[d]-CL.VCM[d])*(pMPC->V[d]-CL.VCM[d]);
			//Increment link in list
			pMPC = pMPC->next;
		}
		prob *= -0.5*(sc_fctr*sc_fctr-1.)*CL.MASS;
		prob /= KBT;
		prob = exp( prob );
		prob *= pow( sc_fctr, DIM*(CL.POP - 1) );
		if( prob>1. ) prob = 1.;
	}

	//Randomly check against this probability of acceptance. Scale if passes
	if( genrand_real() < prob ) if( CL.pp != NULL ) {
		pMPC = CL.pp;
		while(pMPC != NULL) {
			for( d=0;d<DIM;d++ ) pMPC->V[d] *= sqrt( KBT/KBTNOW );
			//Increment link in list
			pMPC = pMPC->next;
		}
	}
}
double TEMP( particleMPC *pp,spec SP[],bc WALL[],double VEL[] ) {
/*
   This routine sums the energy to calculate the temperature
   m*av( v^2 )/2 = dim *KBT / 2
*/
	double KBT = 0.;
	double temp = 0.;
	int i,j,k,c = 0;
	//Fluid particleMPC contribution
	for( i=0; i<GPOP; i++ ) {
		temp = ((pp+i)->V[0]-VEL[0])*((pp+i)->V[0]-VEL[0]);
		temp += ((pp+i)->V[1]-VEL[1])*((pp+i)->V[1]-VEL[1]);
		temp += ((pp+i)->V[2]-VEL[2])*((pp+i)->V[2]-VEL[2]);
		KBT +=  SP[(pp+i)->SPID].MASS * temp;
// 		KBT += SP[(pp+i)->SPID].M * ( (pp+i)->V[0]*(pp+i)->V[0]+(pp+i)->V[1]*(pp+i)->V[1]+(pp+i)->V[2]*(pp+i)->V[2] );
	}
	for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
		//BC object's kinetic contribution
		temp = (WALL[i].V[0]-VEL[0])*(WALL[i].V[0]-VEL[0]);
		temp += (WALL[i].V[1]-VEL[1])*(WALL[i].V[1]-VEL[1]);
		temp += (WALL[i].V[2]-VEL[2])*(WALL[i].V[2]-VEL[2]);
		KBT += WALL[i].MASS * temp;
// 		KBT += WALL[i].M * (WALL[i].V[0]*WALL[i].V[0] + WALL[i].V[1]*WALL[i].V[1] + WALL[i].V[2]*WALL[i].V[2]);
		//BC object's rotational contribution
		for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) KBT += WALL[i].L[j] * WALL[i].I[j][k] * WALL[i].L[k];
		c++;
	}
	//Equipartition function ( sum(mv^2/2) = DIM*N*kBT/2 )
	KBT /= ( (double)((GPOP+c)*DIM) );
	return KBT;
}
double calcE( particleMPC *p,spec *pSP,bc WALL[] ) {
	double E,TE = 0.;
	int i,j,k;
	for( i=0; i<GPOP; i++ ) {
		E = 0.;
		for( j=0; j<DIM; j++ ) E += (p+i)->V[j] * (p+i)->V[j];
		E *= 0.5 * (pSP+(p+i)->SPID)->MASS;
		TE += E;
	}
	for( i=0; i<NBC; i++ ) if( WALL[i].DSPLC ) {
		//Kinetic
		E = 0.;
		for( j=0; j<DIM; j++ ) E += WALL[i].V[j]*WALL[i].V[j];
		E *= 0.5 * WALL[i].MASS;
		TE += E;
		//Rotational
		E = 0.;
		for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) E += WALL[i].L[j] * WALL[i].I[j][k] * WALL[i].L[j];
		E *= 0.5;
		TE += E;
	}
	return TE;
}
double calcE_MPC( particleMPC *p,spec *pSP ) {
/*
   This routine calculates the kinetic energy of a single particle
*/
	double E = 0.;
	int j;

	for( j=0; j<DIM; j++ ) E += p->V[j] * p->V[j];
	E *= 0.5 * (pSP+p->SPID)->MASS;
	return E;
}
double calcE_BC( bc *WALL ) {
/*
   This routine calculates the energy of a single BC solute
*/
	double E=0., TE=0.;
	int j,k;

	if( WALL->DSPLC ) {
		//Translational
		for( j=0; j<DIM; j++ ) E += WALL->V[j] * WALL->V[j];
		E *= 0.5 * WALL->MASS;
		TE += E;
		//Rotational
		E = 0.;
		for( j=0; j<_3D; j++ ) for( k=0; k<_3D; k++ ) E += WALL->L[j] * WALL->I[j][k] * WALL->L[k];
		E *= 0.5;
		TE += E;
	}

	return TE;
}
double calcE_LC( cell ***CL,int LC,double MFPOT ) {
/*
   This routine sums the potential orientational energy of all nematic particles
   Potential energy so is really negative
*/
	double wmf=0.;
	double S,un,DIR[_3D],u[_3D];
	particleMPC *tmpc;	//Temporary particleMPC
	int a,b,c,i;
	//double invdim=1./((double)DIM);

	if( LC ) for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) if( CL[a][b][c].POP > 1 ) {
		S = CL[a][b][c].S;
		for( i=0; i<DIM; i++ ) DIR[i] = CL[a][b][c].DIR[i];
		tmpc = CL[a][b][c].pp;
		while( tmpc != NULL ) {
			for( i=0; i<DIM; i++ ) u[i] = tmpc->U[i];
			un = dotprod( u,DIR,DIM );
			wmf += S*un*un;
			//wmf += (1.-S)*invdim;		//Don't include constant (wrt u.n) term
			tmpc = tmpc->next;
		}
	}
	wmf*=MFPOT;
	return wmf;
}
void avVel( cell ***CL,double AVVEL[] ) {
/*
   This routine finds the average global scalar order parameter.
*/
	int a,b,c,d;
	for( d=0; d<DIM; d++ ) AVVEL[d]=0.;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) for( d=0; d<DIM; d++ ) AVVEL[d] += CL[a][b][c].VCM[d];
	for( d=0; d<DIM; d++ ) AVVEL[d] /= (double)(XYZ[0]*XYZ[1]*XYZ[2]);
}
double avEnstrophy( cell ***CL ) {
/*
   This routine finds the global average enstrophy E (mean squared vorticity)
*/
	int a,b,c;
	double w[_3D];
	double E;

	E=0.;
	for( a=0; a<XYZ_P1[0]; a++ ) for( b=0; b<XYZ_P1[1]; b++ ) for( c=0; c<XYZ_P1[2]; c++ ) {
		w[0]=(CL[a][b][c].E[2][1] - CL[a][b][c].E[1][2]);
		w[1]=(CL[a][b][c].E[0][2] - CL[a][b][c].E[2][0]);
		w[2]=(CL[a][b][c].E[1][0] - CL[a][b][c].E[0][1]);
		E += dotprod( w,w, _3D );
	}
	E /= (XYZ[0]*XYZ[1]*XYZ[2]);
	E *= 0.5;
	return E;
}
