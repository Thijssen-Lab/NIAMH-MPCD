#ifndef GLOBALS_H
#define GLOBALS_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ GLOBAL VARIABLES ************ */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

int DBUG;	//Debugging flag
int DIM;	//Dimensions
int MDmode;	//The MD coupling mode

int GPOP;	//Total number of particles in the system
int NSPECI;	//NSPECI is the number of species present in the system
int NBC;	//NBC is the number of boundaries present in the system
int NS;		//NS is the number of swimmers in the system
double nDNST;	//The particle number density of the fluid
double mDNST;	//The mass density of the fluid

int XYZ[3];		//x,y and z dimensions of the control volume
int XYZ_P1[3];	//x,y and z dimensions of the control volume plus 1
int XYZPBC[3];	//Flags whether x,y or z dimensions are wrapped with periodic BCs
int maxXYZ;		//Maximum dimension

#endif
