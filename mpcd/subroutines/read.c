# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include<math.h>

# include "../headers/SRDclss.h"
# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/pout.h"
# include "../headers/mtools.h"
# include "../headers/cJson.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ******* ROUTINES FOR READING FILES ******* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void checkRead( int flag,char failure[],char file[]) {
	/*
	   If a read from file fileName fails
		 then write to fsynopsis and stop the program
	*/
	if(flag<0) {
		printf( "\nError: Reached end of file before read %s from %s.\n",failure,file );
		exit( 1 );
	}
	else if(flag==0) {
		printf( "\nError: Failed to read %s from %s.\n",failure,file );
		exit( 1 );
	}
}
void readin( char fpath[],inputList *in,spec **SP,particleMPC **pSRD,cell ****CL,int *MDmode ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the values to what
   is in the input file input.inp
*/
	FILE *finput;
	int i,j,MS,read;
	double MF;
	char STR[100],inSTR[100];

	strcpy( inSTR,fpath );
	strcat( inSTR,"input.inp" );
	finput = fopen( inSTR, "r" );
	if( !finput ) {					// file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",inSTR );
		exit( 1 );
	}

	//Read the dimensionality of the system
	read=fscanf( finput,"%d %s",&DIM,STR );
	checkRead( read,"dimensionality",inSTR);
	//Read system dimensions
	for( i=0; i<_3D; i++ ) {
		read=fscanf( finput,"%i %s",&XYZ[i],STR );
		checkRead( read,"system size",inSTR);
	}
	for(i=0; i<_3D; i++ ) XYZ_P1[i] = XYZ[i]+1;
	//Read the thermal energy KBT
	read=fscanf( finput,"%lf %s",&(in->KBT),STR );
	checkRead( read,"thermal energy",inSTR);
	//Read if do initial transformation to rest frame
	read=fscanf( finput,"%d %s",&(in->RFRAME),STR );
	checkRead( read,"Galilean transform",inSTR);
	//Read how often transformation to rest frame
	read=fscanf( finput,"%d %s",&(in->zeroNetMom),STR );
	checkRead( read,"rest frame",inSTR);
	//Read if do random shift of cells
	read=fscanf( finput,"%d %s",&(in->GALINV),STR );
	checkRead( read,"random shift of cells",inSTR);
	//Read thermostat mode/technique
	read=fscanf( finput,"%d %s",&(in->TSTECH),STR );
	checkRead( read,"thermostat mode",inSTR);
	//Read collision operator technique/mode
	read=fscanf( finput,"%d %s",&(in->RTECH),STR );
	checkRead( read,"collision mode",inSTR);
	//Read liquid crystal (0 if not, 1 if LC)
	read=fscanf( finput,"%d %s",&(in->LC),STR );
	checkRead( read,"liquid crystal mode",inSTR);
	//Read the thermal relaxation time scale
	read=fscanf( finput,"%lf %s",&(in->TAU),STR );
	checkRead( read,"relaxation time",inSTR);
	//Read the rotation angle used
	read=fscanf( finput,"%lf %s",&(in->RA),STR );
	checkRead( read,"rotation angle",inSTR);
	//Read the Langevin thermostat friction coefficient
	read=fscanf( finput,"%lf %s",&(in->FRICCO),STR );
	checkRead( read,"friction coefficient",inSTR);
	//Read the liquid crystal mean-field potential
	read=fscanf( finput,"%lf %s",&(in->MFPOT),STR );
	checkRead( read,"LC mean-field potential",inSTR);
	//Read the constant external acceleration
	for( i=0; i<_3D; i++ ) {
		read=fscanf( finput,"%lf %s",&(in->GRAV[i]),STR );
		checkRead( read,"acceleration",inSTR);
	}
	//Read the constant external magnetic field
	for( i=0; i<_3D; i++ ) {
		read=fscanf( finput,"%lf %s",&(in->MAG[i]),STR );
		checkRead( read,"magnetic field",inSTR);
	}
	//Read the time (or number of loops) of the warmup
	read=fscanf( finput,"%d %s",&(in->warmupSteps),STR );
	checkRead( read,"warmup time",inSTR);
	//Read the total time (or number of loops) of the simulation
	read=fscanf( finput,"%d %s",&(in->simSteps),STR );
	checkRead( read,"total time",inSTR);
	//Read the time step of one iteration
	read=fscanf( finput,"%lf %s",&(in->dt),STR );
	checkRead( read,"time step",inSTR);
	//Read the random seed (0 if read from time)
	read=fscanf( finput,"%ld %s",&(in->seed),STR );
	checkRead( read,"random seed",inSTR);
	//Read the MD coupling mode
	read=fscanf( finput,"%d %s",MDmode,STR );
	checkRead( read,"MD coupling",inSTR);
	//Read the number of MD steps per SRD step
	read=fscanf( finput,"%d %s",&(in->stepsMD),STR );
	checkRead( read,"MD time steps per SRD step",inSTR);
	//Read the species properties
	//Read the number of species
	read=fscanf( finput,"%d %s",&NSPECI,STR );
	checkRead( read,"number of species",inSTR);
	//Allocate the needed amount of memory for the species SP
	(*SP) = (spec*) malloc( NSPECI * sizeof( spec ) );
	for( i=0; i<NSPECI; i++ ) {
		//Read the species' mass
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"mass",inSTR);
		(*SP+i)->MASS = MF;
		//Read the species' population
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"population",inSTR);
		(*SP+i)->POP = MS;
		//Read the species' position distribution function
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"positional distribution",inSTR);
		(*SP+i)->QDIST = MS;
		//Read the species' velocity distribution function
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"velocity distribution",inSTR);
		(*SP+i)->VDIST = MS;
		//Read the species' orienation distribution function
		read=fscanf( finput,"%i %s",&MS,STR );
		checkRead( read,"orientational distribution",inSTR);
		(*SP+i)->ODIST = MS;
		//Read the binary fluid interaction matrix for this species with all other species
		for( j=0; j<NSPECI; j++ ) {
			//Read the species' interaction matrix with other species
			read=fscanf( finput,"%lf %s",&MF,STR );
			checkRead( read,"interaction matrix",inSTR);
			(*SP+i)->M[j] = MF;
		}
		//Read the species' rotational friction coefficient
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"rotational friction",inSTR);
		(*SP+i)->RFC = MF;
		//Read the species' tumbling parameter
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"effective rod length",inSTR);
		(*SP+i)->LEN = MF;
		//Read the species' tumbling parameter
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"tumbling parameter",inSTR);
		(*SP+i)->TUMBLE = MF;
		//Read the species' susceptibility to shear
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"shear susceptibility",inSTR);
		(*SP+i)->CHIHI = MF;
		//Read the species' magnetic susceptibility anisotropy
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"magnetic susceptibility",inSTR);
		(*SP+i)->CHIA = MF;
		//Read the species' activity
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"activity",inSTR);
		(*SP+i)->ACT = MF;
		//Read the species' damping friction coefficient
		read=fscanf( finput,"%lf %s",&MF,STR );
		checkRead( read,"damping friction",inSTR);
		(*SP+i)->DAMP = MF;
	}
	//Total Number of particleMPCs
	GPOP = 0;
	for( i=0; i<NSPECI; i++ ) GPOP += (*SP+i)->POP;
	(*pSRD) = (particleMPC*) malloc( GPOP * sizeof( particleMPC ) );

	//Allocate memory for the cells
	//Allocate rows (x first)
	*CL = (cell***) malloc( XYZ_P1[0] * sizeof( cell** ) );
	//For each x-element allocate the y columns
	for( i=0; i<XYZ_P1[0]; i++ ) {
		(*CL)[i] = (cell**) malloc( XYZ_P1[1] * sizeof( cell* ) );
		//For each y-element allocate the z columns
		for( j=0; j<XYZ_P1[1]; j++ ) {
			(*CL)[i][j] = (cell*) malloc( XYZ_P1[2] * sizeof( cell ) );
		}
	}

	fclose( finput );
}
void readpc( char fpath[],outputFlagsList *out ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the flags which tells
  the program what data to output
*/
	FILE *finput;
	char STR[100],inSTR[100];
	int read;

	strcpy( inSTR,fpath );
	strcat( inSTR,"printcom.inp" );
	finput = fopen( inSTR,"r" );
	if( !finput ) {				// file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",inSTR );
		exit( 1 );
	}
	//Read debug mode
	read=fscanf( finput,"%d %s",&DBUG,STR );
	checkRead( read,"debug mode",inSTR);
	//Read how often the detailed species' trajectories are printed
	read=fscanf( finput,"%d %s",&(out->TRAJOUT),STR );
	checkRead( read,"detailed species' trajectories",inSTR);
	//Read which (number of) species whose detailed trajectories are printed
	read=fscanf( finput,"%d %s",&(out->printSP),STR );
	checkRead( read,"number of species whose detailed trajectories are printed",inSTR);
	//Read how often the coarse data is outputted
	read=fscanf( finput,"%d %s",&(out->COAROUT),STR );
	checkRead( read,"coarse data",inSTR);
	//Read how often the flow field is outputted
	read=fscanf( finput,"%d %s",&(out->FLOWOUT),STR );
	checkRead( read,"flow field",inSTR);
	//Read how often the total average MPCD velocity is outputted
	read=fscanf( finput,"%d %s",&(out->AVVELOUT),STR );
	checkRead( read,"total average MPCD velocity",inSTR);
	//Read how often the local director and scalar order parameter fields are outputted
	read=fscanf( finput,"%d %s",&(out->ORDEROUT),STR );
	checkRead( read,"director and scalar order parameter fields",inSTR);
	//Read how often the local order parameter tensor field is outputted
	read=fscanf( finput,"%d %s",&(out->QTENSOUT),STR );
	checkRead( read,"order parameter tensor field",inSTR);
	//Read how often the local order parameter tensor field in reciprical space is outputted
	read=fscanf( finput,"%d %s",&(out->QKOUT),STR );
	checkRead( read,"order parameter tensor field in reciprical space",inSTR);
	//Read how often the orientational energy field is outputted
	read=fscanf( finput,"%d %s",&(out->ENFIELDOUT),STR );
	checkRead( read,"orientational energy field",inSTR);
	//Read how often the colour/conc/phi field is outputted
	read=fscanf( finput,"%d %s",&(out->SPOUT),STR );
	checkRead( read,"colour/conc field",inSTR);
	//Read how often the pressure field is outputted
	read=fscanf( finput,"%d %s",&(out->PRESOUT),STR );
	checkRead( read,"pressure field",inSTR);
	//Read how often the orientational energy from neighbours is outputted
	read=fscanf( finput,"%d %s",&(out->ENNEIGHBOURS),STR );
	checkRead( read,"orientational energy from neighbours",inSTR);
	//Read how often the total average scalar order parameter is outputted
	read=fscanf( finput,"%d %s",&(out->AVSOUT),STR );
	checkRead( read,"average scalar order parameter",inSTR);
	//Read how often the standard deviation of the number per cell is outputted
	read=fscanf( finput,"%d %s",&(out->DENSOUT),STR );
	checkRead( read,"standard deviation of the number per cell",inSTR);
	//Read how often the total average enstrophy is outputted
	read=fscanf( finput,"%d %s",&(out->ENSTROPHYOUT),STR );
	checkRead( read,"total average enstrophy",inSTR);
	//Read how often the velocity distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTVELOUT),STR );
	checkRead( read,"velocity distribution",inSTR);
	//Read how often the speed distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTSPEEDOUT),STR );
	checkRead( read,"speed distribution",inSTR);
	//Read how often the vorticity distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTVORTOUT),STR );
	checkRead( read,"vorticity distribution",inSTR);
	//Read how often the enstrophy distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTENSTROUT),STR );
	checkRead( read,"enstrophy distribution",inSTR);
	//Read how often the director distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTDIROUT),STR );
	checkRead( read,"director distribution",inSTR);
	//Read how often the scalar order parameter distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTSOUT),STR );
	checkRead( read,"scalar order parameter distribution",inSTR);
	//Read how often the number per cell distribution should be outputted
	read=fscanf( finput,"%d %s",&(out->HISTNOUT),STR );
	checkRead( read,"number per cell distribution",inSTR);
	//Read how often the solid trajectory data is outputted
	read=fscanf( finput,"%d %s",&(out->SOLOUT),STR );
	checkRead( read,"solid trajectory",inSTR);
	//Read how often the topological charge field data is outputted
	read=fscanf( finput,"%d %s",&(out->DEFECTOUT),STR );
	checkRead( read,"topological charge field",inSTR);
	//Read how often the energy data is outputted
	read=fscanf( finput,"%d %s",&(out->ENOUT),STR );
	checkRead( read,"energy",inSTR);
	//Read how often the velocity-velocity correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CVVOUT),STR );
	checkRead( read,"velocity-velocity correlation",inSTR);
	//Read how often the director-director correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CNNOUT),STR );
	checkRead( read,"director-director correlation",inSTR);
	//Read how often the vorticity-vorticity correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CWWOUT),STR );
	checkRead( read,"vorticity-vorticity correlation",inSTR);
	//Read how often the density-density correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CDDOUT),STR );
	checkRead( read,"density-density correlation",inSTR);
	//Read how often the order-order correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CSSOUT),STR );
	checkRead( read,"order-order correlation",inSTR);
	//Read how often the phase-phase (binary fluid) correlation data is outputted
	read=fscanf( finput,"%d %s",&(out->CPPOUT),STR );
	checkRead( read,"phase-phase (binary fluid) correlation",inSTR);
	//Read how often the energy spectrum data is outputted
	read=fscanf( finput,"%d %s",&(out->ENERGYSPECTOUT),STR );
	checkRead( read,"energy spectrum",inSTR);
	//Read how often the enstrophy spectrum data is outputted
	read=fscanf( finput,"%d %s",&(out->ENSTROPHYSPECTOUT),STR );
	checkRead( read,"enstrophy spectrum",inSTR);
	//Read how often the Binder cumulant is outputted
	read=fscanf( finput,"%d %s",&(out->BINDER),STR );
	checkRead( read,"Binder cumulant",inSTR);
	//Read the Binder cumulant bin size
	read=fscanf( finput,"%d %s",&(out->BINDERBIN),STR );
	checkRead( read,"Binder cumulant bin size",inSTR);
	//Read how often the swimmers' positions are outputted
	read=fscanf( finput,"%d %s",&(out->SWOUT),STR );
	checkRead( read,"swimmers' positions",inSTR);
	//Read how often the swimmers' orientations are outputted
	read=fscanf( finput,"%d %s",&(out->SWORIOUT),STR );
	checkRead( read,"swimmers' orientations",inSTR);
	//Read if the swimmers' run/tumble statistics are recorded
	read=fscanf( finput,"%d %s",&(out->RTOUT),STR );
	checkRead( read,"swimmers' run/tumble",inSTR);
	//Read if a synopsis printed
	read=fscanf( finput,"%d %s",&(out->SYNOUT),STR );
	checkRead( read,"synopsis",inSTR);
	//Read if checkpointing is on
	read=fscanf( finput,"%d %s",&(out->CHCKPNT),STR );
	checkRead( read,"checkpointing",inSTR);

	fclose( finput );
}
void bcin( FILE *fbc,bc *WALL,char fname[] ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the values to what
   is in the input file bc.inp
 */
	int x,read;
	double l;
	char LABEL[160];

	read=fscanf( fbc,"%s",LABEL );
	checkRead( read,"bc",fname);

	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->COLL_TYPE = x;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->PHANTOM = x;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->E = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->Q[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->Q[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->Q[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->V[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->V[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->V[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->O[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->O[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->O[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->L[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->L[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->L[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->G[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->G[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->G[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->AINV[0] = l;
	if( fneq(WALL->AINV[0],0.0) ) WALL->A[0] = 1.0/WALL->AINV[0];
	else WALL->A[0] = 0.0;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->AINV[1] = l;
	if( fneq(WALL->AINV[1],0.0) ) WALL->A[1] = 1.0/WALL->AINV[1];
	else WALL->A[1] = 0.0;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->AINV[2] = l;
	if( fneq(WALL->AINV[2],0.0) ) WALL->A[2] = 1.0/WALL->AINV[2];
	else WALL->A[2] = 0.0;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->ROTSYMM[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->ROTSYMM[1] = l;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->ABS = x;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[2] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->P[3] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->R = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DT = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVT = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVxyz[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVxyz[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DVxyz[2] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MVN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MVT = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUN = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUT = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUxyz[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUxyz[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MUxyz[2] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DUxyz[0] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DUxyz[1] = l;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->DUxyz[2] = l;

	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->KBT = l;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->DSPLC = x;
	read=fscanf( fbc,"%d %s",&x,LABEL );
	checkRead( read,"bc",fname);
	WALL->INV = x;
	read=fscanf( fbc,"%lf %s",&l,LABEL );
	checkRead( read,"bc",fname);
	WALL->MASS = l;

	// Set the planar flag
	if( feq(WALL->P[0],1.0) && feq(WALL->P[1],1.0) && feq(WALL->P[2],1.0) ) {
		// Left or right wall
		if( fneq(WALL->A[0],0.0) && feq(WALL->A[1],0.0) && feq(WALL->A[2],0.0) ) WALL->PLANAR = 1;
		// Top or bottom wall
		else if( feq(WALL->A[0],0.0) && fneq(WALL->A[1],0.0) && feq(WALL->A[2],0.0) ) WALL->PLANAR = 1;
		// Far or near wall
		else if( feq(WALL->A[0],0.0) && feq(WALL->A[1],0.0) && fneq(WALL->A[2],0.0) ) WALL->PLANAR = 1;
		else WALL->PLANAR = 0;
	}
	else WALL->PLANAR = 0;
}
void setPBC( bc *WALL ) {
/*
   Determines if a boundary is a planar periodic boundary
*/
	int i;
	//Check if any axis is a planar PBC
	for( i=0; i<_3D; i++ ) {
		// Is the BC a plane?
		if( feq(WALL->P[i],1.0) ) {
			// Does it point in the correct direction?
			if( feq( fabs(WALL->AINV[i]),1.0 ) ) {
				// Does it keep the velocity pointed in the correct direction?
				if( WALL->MVN>0.0 && WALL->MVT>0.0 ) {
					//Does it modify the position?
					if( fneq(WALL->DN,0.0) || fneq(WALL->DT,0.0) ) {
						//Then it is a periodic boundary condition
						XYZPBC[i]=1;
					}
				}
			}
		}
	}
}

void readbc( char fpath[],bc **WALL ) {
/*
   By calling bcin this function sets each of
   the boundaries
*/
	FILE *fbc;
	int i,read;
	char STR[100],inSTR[100];

	strcpy( inSTR,fpath );
	strcat( inSTR,"bc.inp" );
	fbc = fopen(inSTR, "r");
	if( !fbc ) { // file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",inSTR );
		exit( 1 );
	}

	read=fscanf( fbc,"%s",STR );
	checkRead( read,"bc",inSTR);
	//Read the number of boundary conditions being used
	read=fscanf( fbc,"%d %s",&NBC,STR );
	checkRead( read,"bc",inSTR);
	//Allocate the needed amount of memory for the boundary conditions WALLs
	(*WALL) = (bc*) malloc( NBC * sizeof( bc ) );
	//Read the parameters for each of the BCs
	for( i=0; i<NBC; i++ ) bcin( fbc,(*WALL+i),inSTR );
	//Determine if any BCs are periodic boundaries
	for( i=0; i<_3D; i++ ) XYZPBC[i]=0;
	for( i=0; i<NBC; i++ ) setPBC( (*WALL+i) );
	fclose( fbc );
}
void readchckpnt( char fpath[],inputList *in,spec **SP,particleMPC **pSRD,cell ****CL,int *MDmode,bc **WALL,outputFlagsList *out,int *runtime,int *warmtime,kinTheory *theory,double *AVVEL, double *AVS,double avDIR[_3D],double *S4,double *stdN,double *KBTNOW,double AVV[_3D],double AVNOW[_3D] ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the values to what
   is in the input file input.inp
*/
	FILE *finput;
	int i,j;
	char STR[100];

	strcpy( STR,fpath );
	strcat( STR,"checkpoint.dat" );
	finput = fopen( STR, "r" );
	if( !finput ) {					// file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",STR );
		exit( 1 );
	}

	if(fscanf( finput,"%d %d %d %d %lf %lf",&DIM,&XYZ[0],&XYZ[1],&XYZ[2],&(in->KBT),KBTNOW ));
	else printf("Warning: Failed to read dimensionality, size or temperature.\n");
	for(i=0; i<_3D; i++ ) XYZ_P1[i] = XYZ[i]+1;
	if(fscanf( finput,"%d %d %d %d %d",&(in->RFRAME),&(in->zeroNetMom),&(in->TSTECH),&(in->RTECH),&(in->LC) ));
	else printf("Warning: Failed to Galilean transform, rest frame, thermostat mode, collision mode or liquid crystal mode.\n");
	if(fscanf( finput,"%lf %lf %lf %lf",&(in->TAU),&(in->RA),&(in->FRICCO),&(in->MFPOT) ));		//Read the thermal relaxation time scale
	else printf("Warning: Failed to read relaxation time, rotation angle, friction coefficient or mean-field potential.\n");

	if(fscanf( finput,"%lf %lf %lf",&(in->GRAV[0]),&(in->GRAV[1]),&(in->GRAV[2]) ));	//Read the constant external acceleration
	else printf("Warning: Failed to read acceleration.\n");
	if(fscanf( finput,"%lf %lf %lf",&(in->MAG[0]),&(in->MAG[1]),&(in->MAG[2]) ));	//Read the constant external magnetic field
	else printf("Warning: Failed to read magnetic field.\n");

	if(fscanf( finput,"%d %d %lf",&(in->warmupSteps),&(in->simSteps),&(in->dt) ));		//Read time
	else printf("Warning: Failed to read time.\n");
	if(fscanf( finput,"%ld",&(in->seed) ));	//Read the random seed (0 if read from time)
	else printf("Warning: Failed to read random seed.\n");

	if(fscanf( finput,"%d %d",MDmode,&(in->stepsMD) ));	//Read the MD coupling mode
	else printf("Warning: Failed to read MD coupling.\n");
	if(fscanf( finput,"%d %d",&GPOP,&NSPECI ));	//Read the number of MPC particles
	else printf("Warning: Failed to read total number of particles or number of species.\n");

	if(fscanf( finput,"%d %d %lf %lf %d %d",runtime,warmtime,&(in->C),&(in->S),&(in->GRAV_FLAG),&(in->MAG_FLAG) ));//Read program variables
	else printf("Warning: Failed to read various program variables.\n");
	if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&(theory->MFP), &(theory->VISC), &(theory->THERMD), &(theory->SDIFF), &(theory->sumM), AVVEL, AVS, &avDIR[0], &avDIR[1], &avDIR[2], S4, stdN, &nDNST, &mDNST ));//Read program variables
	else printf("Warning: Failed to read various program variables.\n");
	if(fscanf( finput,"%lf %lf %lf %lf %lf %lf",&AVV[0], &AVV[1], &AVV[2], &AVNOW[0], &AVNOW[1], &AVNOW[2] ));//Read program variables
	else printf("Warning: Failed to read average velocities.\n");

	//Read output
	if(fscanf( finput,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",&DBUG, &(out->TRAJOUT), &(out->printSP), &(out->COAROUT), &(out->FLOWOUT), &(out->AVVELOUT), &(out->ORDEROUT), &(out->QTENSOUT), &(out->QKOUT), &(out->AVSOUT), &(out->SOLOUT), &(out->ENOUT), &(out->ENFIELDOUT), &(out->ENNEIGHBOURS), &(out->ENSTROPHYOUT), &(out->DENSOUT), &(out->CVVOUT), &(out->CNNOUT), &(out->CWWOUT), &(out->CDDOUT), &(out->CSSOUT), &(out->CPPOUT), &(out->BINDER), &(out->BINDERBIN), &(out->SYNOUT), &(out->CHCKPNT) ));
	else printf("Warning: Failed to read output.\n");
	if(fscanf( finput,"%d %d %d %d %d %d %d",&(out->HISTVELOUT), &(out->HISTSPEEDOUT), &(out->HISTVORTOUT), &(out->HISTENSTROUT), &(out->HISTDIROUT), &(out->HISTSOUT), &(out->HISTNOUT) ));
	else printf("Warning: Failed to read histogram output.\n");
	if(fscanf( finput,"%d %d %d",&(out->ENERGYSPECTOUT), &(out->ENSTROPHYSPECTOUT), &(out->DEFECTOUT) ));
	else printf("Warning: Failed to read histogram output.\n");

	//Allocate the needed amount of memory for the species SP
	(*SP) = (spec*) malloc( NSPECI * sizeof( spec ) );
	for( i=0; i<NSPECI; i++ ) {
		if(fscanf( finput,"%lf %i %i %i %i %lf %lf %lf %lf %lf %lf %lf",&((*SP+i)->MASS), &((*SP+i)->POP), &((*SP+i)->QDIST), &((*SP+i)->VDIST), &((*SP+i)->ODIST), &((*SP+i)->RFC), &((*SP+i)->LEN), &((*SP+i)->TUMBLE), &((*SP+i)->CHIHI), &((*SP+i)->CHIA), &((*SP+i)->ACT), &((*SP+i)->DAMP) ));	//Read the species' mass
		else printf("Warning: Failed to read species %i.\n",i);
		for( j=0; j<NSPECI; j++ ) {
			//Read the species' interaction matrix with other species
			if(fscanf( finput,"%lf ",&((*SP+i)->M[j]) ));	//Read the species' interactions
			else printf("Warning: Failed to read species %d interaction with %d.\n",i,j);
		}
	}
	//Check total number of particleMPCs
	j = 0;
	for( i=0; i<NSPECI; i++ ) j += (*SP+i)->POP;
	if(GPOP!=j) {
		printf("Warning: GPOP does not match sum of species populations.\n");
		GPOP=j;
	}
	(*pSRD) = (particleMPC*) malloc( GPOP * sizeof( particleMPC ) );

	//Allocate memory for the cells
	//Allocate rows (x first)
	*CL = (cell***) malloc( XYZ_P1[0] * sizeof( cell** ) );
	//For each x-element allocate the y columns
	for( i=0; i<XYZ_P1[0]; i++ ) {
		(*CL)[i] = (cell**) malloc( XYZ_P1[1] * sizeof( cell* ) );
		//For each y-element allocate the z columns
		for( j=0; j<XYZ_P1[1]; j++ ) {
			(*CL)[i][j] = (cell*) malloc( XYZ_P1[2] * sizeof( cell ) );
		}
	}
	//Read the BCs
	if(fscanf( finput,"%d",&NBC ));		//Read the number of BCs
	else printf("Warning: Failed to read number of BCs.\n");
	for( i=0; i<NBC; i++ ) {
		if(fscanf( finput,"%d %d %lf %lf %lf %lf %lf %lf %lf",&((*WALL+i)->COLL_TYPE), &((*WALL+i)->PHANTOM), &((*WALL+i)->E), &((*WALL+i)->Q[0]), &((*WALL+i)->Q[1]), &((*WALL+i)->Q[2]), &((*WALL+i)->V[0]), &((*WALL+i)->V[1]), &((*WALL+i)->V[2]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &((*WALL+i)->L[0]), &((*WALL+i)->L[1]), &((*WALL+i)->L[2]), &((*WALL+i)->G[0]), &((*WALL+i)->G[1]), &((*WALL+i)->G[2]), &((*WALL+i)->A[0]), &((*WALL+i)->A[1]), &((*WALL+i)->A[2]),&((*WALL+i)->P[0]),&((*WALL+i)->P[1]),&((*WALL+i)->P[2]),&((*WALL+i)->P[3]), &((*WALL+i)->R) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &((*WALL+i)->DN), &((*WALL+i)->DT), &((*WALL+i)->DVN), &((*WALL+i)->DVT), &((*WALL+i)->DVxyz[0]), &((*WALL+i)->DVxyz[1]), &((*WALL+i)->DVxyz[2]), &((*WALL+i)->MVN), &((*WALL+i)->MVT), &((*WALL+i)->MUN), &((*WALL+i)->MUT), &((*WALL+i)->MUxyz[0]), &((*WALL+i)->MUxyz[1]), &((*WALL+i)->MUxyz[2]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%lf %lf %lf %lf %d %d %lf", &((*WALL+i)->DUxyz[0]), &((*WALL+i)->DUxyz[1]), &((*WALL+i)->DUxyz[2]), &((*WALL+i)->KBT), &((*WALL+i)->DSPLC), &((*WALL+i)->INV), &((*WALL+i)->MASS) ));
		else printf("Warning: Failed to read BC %d.\n",i);
		if(fscanf( finput,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &((*WALL+i)->PLANAR), &((*WALL+i)->W), &((*WALL+i)->VOL), &((*WALL+i)->Q_old[0]), &((*WALL+i)->Q_old[1]), &((*WALL+i)->Q_old[2]), &((*WALL+i)->I[0][0]), &((*WALL+i)->I[0][1]), &((*WALL+i)->I[0][2]), &((*WALL+i)->I[1][0]), &((*WALL+i)->I[1][1]), &((*WALL+i)->I[1][2]), &((*WALL+i)->I[2][0]), &((*WALL+i)->I[2][1]), &((*WALL+i)->I[2][2]) ));
		else printf("Warning: Failed to read BC %d.\n",i);
	}

	//Read the MPCD particles
	for( i=0; i<GPOP; i++ ) {
		if(fscanf( finput,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&((*pSRD+i)->S_flag), &((*pSRD+i)->SPID), &((*pSRD+i)->q), &((*pSRD+i)->Q[0]), &((*pSRD+i)->Q[1]), &((*pSRD+i)->Q[2]), &((*pSRD+i)->V[0]), &((*pSRD+i)->V[1]), &((*pSRD+i)->V[2]), &((*pSRD+i)->U[0]), &((*pSRD+i)->U[1]), &((*pSRD+i)->U[2]), &((*pSRD+i)->T[0]), &((*pSRD+i)->T[1]), &((*pSRD+i)->T[2]) ));
		else printf("Warning: Failed to read MPCD particle %d.\n",i);
	}

	fclose( finput );
}

void readarg( int argc, char* argv[], char ip[],char op[], int *inMode ) {
	int arg;
	int strln;
	// Default
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf("Read arguments\n");
	#endif
	strcpy(ip,"mpcd/data/");
	strcpy(op,"mpcd/data/");
	for(arg=1; arg<argc; arg++) {
		// Check for a dash
		if(argv[arg][0]=='-') {
			switch (argv[arg][1]) {
				case 'i':
					*inMode = 0;
					arg++;
					sprintf(ip,"%s",argv[arg]);
					break;
				case 'L': // legacy 
					if (argv[arg][2] == 'i'){ // arg 'Li' for legacy input
						*inMode = 1;
						arg++;
						sprintf(ip,"%s",argv[arg]);
						// Make sure that the directory ends with a "/"
						strln = strlen(ip);
						if( ip[strln-1]!='/' ) strcat( ip,"/" );
						break;
					}
				case 'o':
					arg++;
					sprintf(op,"%s",argv[arg]);
					// Make sure that the directory ends with a "/"
					strln = strlen(op);
					if( op[strln-1]!='/' ) strcat( op,"/" );
					break;
				case 'v':
					printVersionSummary( );
				case 'h':
					printf("\nMPCv090\n");
					printf("For the Polymer Physics Research Group\n");
					printf("At the University of Ottawa\n");
					printf("by Tyler Shendruk\n\nUsage:\n");
					printf("\t-i\t[path to input files]\t\tdefault='mpcd/data/'\n\t-o\t[path to output data files]\tdefault='mpcd/data/'\n\t-v\tprint version summary\n\t-h\t[this help menu]\n\n");
					exit(EXIT_SUCCESS);
					break;
				default:
					printf("Error:\n** Invalid arument: %s\n** Try -h\n",argv[arg]);
					exit(EXIT_FAILURE);
			}
		}
	}
}

/*
	Checks if a given BC, as a cJSON object, is valid and contains all the 
	 	necessary parameters it needs. 
	Returns 1 if valid, 0 if not.
*/
int checkBC(cJSON *bc){
	/*
		Need to check if the following json tags exist, do so using cJSON 
			primitives:
		- "aInv"
		- "P"
		- "R"
		- "DN"
		- "MVN"
		- "MVT"
	*/
	char tagList[6][5] = {"aInv", "P", "R", "DN", "MVN", "MVT"};

	cJSON *temp = NULL;
	int i;
	for (i = 0; i < 6; i++){
		temp = cJSON_GetObjectItemCaseSensitive(bc, tagList[i]);
		if (temp == NULL) return 0;
	}

	return 1; // if you're here without returning then all succesful
}

void readJson( char fpath[], inputList *in, spec **SP, particleMPC **pSRD, 
   cell ****CL, int *MDMode, outputFlagsList *out, bc **WALL, 
   specSwimmer *specS, swimmer **sw){
/*
	Main method for reading in Json, parsing it, and populating ALL inputs
*/
	int i, j; // counting variables

	char* fileStr = NULL;
	if(getFileStr(fpath, &fileStr) != 0){ // read, return on error
		exit(EXIT_FAILURE);
	} 
	printf("==== Read JSON =====\n%s", fileStr);
	printf("\n\n");

	//now can actually parse the json
	cJSON *jObj = cJSON_Parse(fileStr); // create the json object 
	if (jObj == NULL) // error checking
	{
		const char *error_ptr = cJSON_GetErrorPtr();
		if (error_ptr != NULL)
		{
			fprintf(stderr, "Json read error. \nError before: %s\n", error_ptr);
		}
		exit(EXIT_FAILURE);
	}

	// set up input validation routines
	linkedList *jsonTagList = NULL;
	initLL(&jsonTagList);
	linkedList *arrayList = NULL;
	initLL(&arrayList);

	////////////////////////////////////////////////////////////////////////////
	// Perform parsing here. 
	// This is done in the order declared in docs/InputGuide.md
	////////////////////////////////////////////////////////////////////////////
	
	// 1. Old input.inp ////////////////////////////////////////////////////////
	// scroll up to void readin() to see better descriptions & definitions for these

	// dimensionality and domain bounds array
	cJSON *arrDomain = NULL;
	getCJsonArray(jObj, &arrDomain, "domain", jsonTagList, arrayList, 0);
	if(arrDomain != NULL){ // if the can be found in the json
		DIM = cJSON_GetArraySize(arrDomain);
		if(DIM != 2 && DIM != 3){ // check dimensionality is valid
			printf("Error: Dimensionality must be 2 or 3.\n");
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < _3D; i++) {
			if (i == 2 && DIM == 2) XYZ[i] = 1; // if 2D, set z to 1
			else XYZ[i] = cJSON_GetArrayItem(arrDomain, i)->valueint; // get the value
		}		
	} else { // if array cannot be found then fallback to default
		DIM = 2;
		XYZ[0] = 30;
		XYZ[1] = 30;
		XYZ[2] = 1;
	}

	// first set of primitives
	in->KBT = getJObjDou(jObj, "kbt", 1, jsonTagList); // kbt
	in->dt = getJObjDou(jObj, "dt", 0.1, jsonTagList); // dt
	in->simSteps = getJObjInt(jObj, "simSteps", 2000, jsonTagList); // simSteps
	in->warmupSteps = getJObjInt(jObj, "warmUp", 0, jsonTagList); // warmupSteps
	in->RFRAME = getJObjInt(jObj, "rFrame", 1, jsonTagList); // RFRAME
	in->zeroNetMom = getJObjInt(jObj, "zeroNetMom", 0, jsonTagList); // zeroNetMom
	in->GALINV = getJObjInt(jObj, "galInv", 2, jsonTagList); // GALINV
	in->TSTECH = getJObjInt(jObj, "tsTech", 0, jsonTagList); // TSTECH
	in->RTECH = getJObjInt(jObj, "rTech", 2, jsonTagList); // RTECH
	in->LC = getJObjInt(jObj, "lc", 0, jsonTagList); // LC
	in->TAU = getJObjDou(jObj, "tau", 0.5, jsonTagList); // TAU
	in->RA = getJObjDou(jObj, "rotAng", 1.570796, jsonTagList); // rotAng
	in->FRICCO = getJObjDou(jObj, "fricCoef", 1.0, jsonTagList); // fricCo
	in->MFPOT = getJObjDou(jObj, "mfpot", 10.0, jsonTagList); // mfpPot
	
	// grav array
	cJSON *arrGrav = NULL;
	getCJsonArray(jObj, &arrGrav, "grav", jsonTagList, arrayList, 0);
	if (arrGrav != NULL) { // if grav has been found then....
		if (cJSON_GetArraySize(arrGrav) != _3D) { // check dimensionality is valid
			printf("Error: Grav must be a 3D array.\n");
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < _3D; i++) { // get the value
			in->GRAV[i] = cJSON_GetArrayItem(arrGrav, i)->valuedouble; 
		}	
	} else { // if no grav specified then fallback
		in->GRAV[0] = 0;
		in->GRAV[1] = 0;
		in->GRAV[2] = 0;
	}

	// mag array
	cJSON *arrMag = NULL;
	getCJsonArray(jObj, &arrMag, "mag", jsonTagList, arrayList, 0);
	if (arrGrav != NULL) { // if grav has been found then....
		if (cJSON_GetArraySize(arrMag) != _3D) { // check dimensionality is valid
			printf("Error: Mag must be a 3D array.\n");
			exit(EXIT_FAILURE);
		}

		for (i = 0; i < _3D; i++) { // get the value
			in->MAG[i] = cJSON_GetArrayItem(arrMag, i)->valuedouble; 
		}	
	} else { // if no grav specified then fallback
		in->MAG[0] = 0;
		in->MAG[1] = 0;
		in->MAG[2] = 0;
	}

	// second set of primitives
	in->seed = getJObjInt(jObj, "seed", 0, jsonTagList); // seed
	MDmode = getJObjInt(jObj, "mdMode", 0, jsonTagList); // mdMode
	in->stepsMD = getJObjInt(jObj, "stepsMD", 20, jsonTagList); // stepsMD

	// 2. Species //////////////////////////////////////////////////////////////
	// scroll up to void readin() to see better descriptions & definitions for these

	cJSON *arrSpecies = NULL;
	getCJsonArray(jObj, &arrSpecies, "species", jsonTagList, arrayList, 1);
	if(arrSpecies != NULL){ // if this can be found in the json
		NSPECI = cJSON_GetArraySize(arrSpecies); // get the number of species
		
		//Allocate the needed amount of memory for the species SP
		(*SP) = (spec*) malloc( NSPECI * sizeof( spec ) );
		for (i = 0; i < NSPECI; i++) { // loop through the species
			cJSON *objElem = cJSON_GetArrayItem(arrSpecies, i); // get the species object

			// now get first set of primitives
			(*SP+i)->MASS = getJObjDou(objElem, "mass", 1.0, jsonTagList); // mass
			(*SP+i)->POP = getJObjInt(objElem, "pop", 18000, jsonTagList); // pop
			(*SP+i)->QDIST = getJObjInt(objElem, "qDist", 0, jsonTagList); // qDist
			(*SP+i)->VDIST = getJObjInt(objElem, "vDist", 0, jsonTagList); // vDist
			(*SP+i)->ODIST = getJObjInt(objElem, "oDist", 2, jsonTagList); // oDist

			//Read the binary fluid interaction matrix for this species with all other species
			cJSON *arrBFM = NULL;
			getCJsonArray(jObj, &arrBFM, "interMatr", jsonTagList, arrayList, 0);
			if (arrBFM != NULL) { // if grav has been found then....
				if (cJSON_GetArraySize(arrBFM) != NSPECI) { // check dimensionality is valid
					printf("Error: Interaction matrices must have columns of length equal to the number of species.\n");
					exit(EXIT_FAILURE);
				}

				for (j = 0; j < NSPECI; j++) { // get the value
					(*SP+i)->M[j] = cJSON_GetArrayItem(arrBFM, j)->valuedouble; 
				}	
			} else { // if no grav specified then fallback
				for (j = 0; j < NSPECI; j++) { // get the value
					(*SP+i)->M[j] = 0; 
				}	
			}

			// get second set of primitives
			(*SP+i)->RFC = getJObjDou(objElem, "rfc", 0.01, jsonTagList); // rfCoef
			(*SP+i)->LEN = getJObjDou(objElem, "len", 0.007, jsonTagList); // len
			(*SP+i)->TUMBLE = getJObjDou(objElem, "tumble", 2.0, jsonTagList); // tumble
			(*SP+i)->CHIHI = getJObjDou(objElem, "shearSusc", 0.5, jsonTagList); // chiHi
			(*SP+i)->CHIA = getJObjDou(objElem, "magnSusc", 0.001, jsonTagList); // chiA
			(*SP+i)->ACT = getJObjDou(objElem, "act", 0.05, jsonTagList); // act
			(*SP+i)->DAMP = getJObjDou(objElem, "damp", 0.0, jsonTagList); // damp
		}
	} else { // if nothing found in the JSON then fallback to the default
		// setting up a single species with default parameters
		//		note this is just copied from the above code w lines changed
		NSPECI = 1;

		(*SP) = (spec*) malloc( NSPECI * sizeof( spec ) );
		for (i = 0; i < NSPECI; i++) { // loop through the species
			// now get first set of primitives
			(*SP+i)->MASS = 1; // mass
			(*SP+i)->POP = 18000; // pop
			(*SP+i)->QDIST = 0; // qDist
			(*SP+i)->VDIST = 0; // vDist
			(*SP+i)->ODIST = 2; // oDist
			for (j = 0; j < NSPECI; j++) { // interaction matrix
				(*SP+i)->M[j] = 0; 
			}	
			(*SP+i)->RFC = 0.01; // rfCoef
			(*SP+i)->LEN = 0.007; // len
			(*SP+i)->TUMBLE = 2; // tumble
			(*SP+i)->CHIHI = 0.5; // chiHi
			(*SP+i)->CHIA = 0.001; // chiA
			(*SP+i)->ACT = 0.05; // act
			(*SP+i)->DAMP = 0; // damp
	}
	
	// 3. Printcom /////////////////////////////////////////////////////////////
	// scroll up to void readpc() to see better descriptions & definitions for these

	DBUG = getJObjInt(jObj, "debugOut", 3, jsonTagList); // dbug
	out->TRAJOUT = getJObjInt(jObj, "trajOut", 0, jsonTagList); // trajOut
	out->printSP = getJObjInt(jObj, "trajSpecOut", 0, jsonTagList); // printSP
	out->COAROUT = getJObjInt(jObj, "coarseOut", 0, jsonTagList); // coarOut
	out->FLOWOUT = getJObjInt(jObj, "flowOut", 0, jsonTagList); // flowOut
	out->AVVELOUT = getJObjInt(jObj, "avVelOut", 0, jsonTagList); // avVelOut
	out->ORDEROUT = getJObjInt(jObj, "dirSOut", 0, jsonTagList); // orderOut
	out->QTENSOUT = getJObjInt(jObj, "qTensOut", 0, jsonTagList); // qTensOut
	out->QKOUT = getJObjInt(jObj, "qkOut", 0, jsonTagList); // qKOut
	out->ENFIELDOUT = getJObjInt(jObj, "oriEnOut", 0, jsonTagList); // enFieldOut
	out->SPOUT = getJObjInt(jObj, "colourOut", 0, jsonTagList); // spOut
	out->PRESOUT = getJObjInt(jObj, "pressureOut", 0, jsonTagList); // presOut
	out->ENNEIGHBOURS = getJObjInt(jObj, "neighbourEnOut", 0, jsonTagList); // enNeighbours
	out->AVSOUT = getJObjInt(jObj, "avSOut", 0, jsonTagList); // avSOut
	out->DENSOUT = getJObjInt(jObj, "densSDOut", 0, jsonTagList); // densOut
	out->ENSTROPHYOUT = getJObjInt(jObj, "enstrophyOut", 0, jsonTagList); // enStrophOut
	out->HISTVELOUT = getJObjInt(jObj, "histVelOut", 0, jsonTagList); // histVelOut
	out->HISTSPEEDOUT = getJObjInt(jObj, "histSpeedOut", 0, jsonTagList); // histSpeedOut
	out->HISTVORTOUT = getJObjInt(jObj, "histVortOut", 0, jsonTagList); // histVortOut
	out->HISTENSTROUT = getJObjInt(jObj, "histEnsOut", 0, jsonTagList); // histEnstrophyOut
	out->HISTDIROUT = getJObjInt(jObj, "histDirOut", 0, jsonTagList); // histDirOut
	out->HISTSOUT = getJObjInt(jObj, "histSOut", 0, jsonTagList); // histSOut
	out->HISTNOUT = getJObjInt(jObj, "histNOut", 0, jsonTagList); // histNOut
	out->SOLOUT = getJObjInt(jObj, "solidTrajOut", 0, jsonTagList); // solOut
	out->DEFECTOUT = getJObjInt(jObj, "topoFieldOut", 0, jsonTagList); // defectOut
	out->ENOUT = getJObjInt(jObj, "energyOut", 0, jsonTagList); // enOut
	out->CVVOUT = getJObjInt(jObj, "velCorrOut", 0, jsonTagList); // cvvOut
	out->CNNOUT = getJObjInt(jObj, "dirCorrOut", 0, jsonTagList); // cnnOut
	out->CWWOUT = getJObjInt(jObj, "vortCorrOut", 0, jsonTagList); // cwwOut
	out->CDDOUT = getJObjInt(jObj, "densCorrOut", 0, jsonTagList); // cddOut
	out->CSSOUT = getJObjInt(jObj, "orderCorrOut", 0, jsonTagList); // cssOut
	out->CPPOUT = getJObjInt(jObj, "phaseCorrOut", 0, jsonTagList); // cppOut
	out->ENERGYSPECTOUT = getJObjInt(jObj, "energySpecOut", 0, jsonTagList); // energySpectOut
	out->ENSTROPHYSPECTOUT = getJObjInt(jObj, "enstrophySpecOut", 0, jsonTagList); // enstrophySpectOut
	out->BINDER = getJObjInt(jObj, "binderOut", 0, jsonTagList); // binderOut
	out->BINDERBIN = getJObjInt(jObj, "binderBin", 0, jsonTagList); // binderBinOut
	out->SWOUT = getJObjInt(jObj, "swimQOut", 0, jsonTagList); // swOut
	out->SWORIOUT = getJObjInt(jObj, "swimOOut", 0, jsonTagList); // swOriOut
	out->RTOUT = getJObjInt(jObj, "swimROut", 0, jsonTagList); // swVelOut
	out->SYNOUT = getJObjInt(jObj, "synopsisOut", 1, jsonTagList); // swSynOut
	out->CHCKPNT = getJObjInt(jObj, "checkpointOut", 0, jsonTagList); // chkpntOut

	// 3. Boundaries ///////////////////////////////////////////////////////////
	// scroll up to void bcin() to see better descriptions & definitions for these

	cJSON *arrBC = NULL;
	getCJsonArray(jObj, &arrBC, "BC", jsonTagList, arrayList, 1);
	if(arrBC != NULL){ // if this can be found in the json
		NBC = cJSON_GetArraySize(arrBC); // get the number of BCs
		
		//Allocate the needed amount of memory for the BCs
		(*WALL) = (bc*) malloc( NBC * sizeof( bc ) );
		for (i = 0; i < NBC; i++) { // loop through the BCs
			cJSON *objElem = cJSON_GetArrayItem(arrBC, i); // get the BC object
			bc *currWall = *(WALL + i); // get the pointer to the BC we want to write to

			// Do a check if the necessary parameters are present 
			if (!checkBC(objElem)){
				printf("Error: BC %d is missing essential parameters\n", i);
				exit(EXIT_FAILURE);
			}

			///TODO: need to parse these
			currWall->COLL_TYPE = x;
			currWall->PHANTOM = x;
			currWall->E = l;

			currWall->Q[0] = l;
			currWall->Q[1] = l;
			currWall->Q[2] = l;

			currWall->V[0] = l;
			currWall->V[1] = l;
			currWall->V[2] = l;

			currWall->O[0] = l;
			currWall->O[1] = l;
			currWall->O[2] = l;

			currWall->L[0] = l;
			currWall->L[1] = l;
			currWall->L[2] = l;

			currWall->G[0] = l;
			currWall->G[1] = l;
			currWall->G[2] = l;

			currWall->AINV[0] = l;
			if( fneq(currWall->AINV[0],0.0) ) currWall->A[0] = 1.0/currWall->AINV[0];
			else currWall->A[0] = 0.0;
			currWall->AINV[1] = l;
			if( fneq(currWall->AINV[1],0.0) ) currWall->A[1] = 1.0/currWall->AINV[1];
			else currWall->A[1] = 0.0;
			currWall->AINV[2] = l;
			if( fneq(currWall->AINV[2],0.0) ) currWall->A[2] = 1.0/currWall->AINV[2];
			else currWall->A[2] = 0.0;

			currWall->ROTSYMM[0] = l;
			currWall->ROTSYMM[1] = l;
			currWall->ABS = x;

			currWall->P[0] = l;
			currWall->P[1] = l;
			currWall->P[2] = l;
			currWall->P[3] = l;
			currWall->R = l;

			currWall->DN = l;
			currWall->DT = l;
			currWall->DVN = l;
			currWall->DVT = l;
			currWall->DVxyz[0] = l;
			currWall->DVxyz[1] = l;
			currWall->DVxyz[2] = l;
			currWall->MVN = l;
			currWall->MVT = l;

			currWall->MUN = l;
			currWall->MUT = l;
			currWall->MUxyz[0] = l;
			currWall->MUxyz[1] = l;
			currWall->MUxyz[2] = l;
			currWall->DUxyz[0] = l;
			currWall->DUxyz[1] = l;
			currWall->DUxyz[2] = l;

			currWall->KBT = l;
			currWall->DSPLC = x;
			currWall->INV = x;
			currWall->MASS = l;
		}
	} else { // otherwise default
		// allocate memory for BCs as necessary
		///TODO: also take into account the dimension

		//set up PBCs on the xy plane based on the domain dimensions
		///TODO

		if(DIM == _3D){ // if in 3d then also do the z axis
			///TODO	
		}
	}

	/// TODO: :))))))))))))

	// input verification step
   	verifyJson(jObj, jsonTagList, arrayList);

	// clear memory
	free(fileStr);

	exit(EXIT_SUCCESS); ///TODO: temp while testing, remove eventually
}