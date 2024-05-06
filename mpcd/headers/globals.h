///
/// @file
///
/// @brief Contains global non-constant variables used throughout the code
///
/// Contains all true global variables used throughout the code. Most of these relate to the control volume size, and
/// size of lists.
///

#ifndef GLOBALS_H
#define GLOBALS_H

# include "definitions.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ GLOBAL VARIABLES ************ */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/// @brief The debugging/ verbosity level of the simulation.
int DBUG;
/// @brief The dimension of the simulation. Must be 1, 2, or 3.
int DIM;
/// @brief The MD mode of the simulation. Must be 0, 1, or 2.
///
/// @see noMD
/// @see MDinMPC
/// @see MPCinMD
int MDmode;

/// @brief The total number of MPCD particles in the simulation.
int GPOP;
/// @brief The total number of species in the simulation.
int NSPECI;
/// @brief The total number of boundaries in the simulation.
int NBC;
/// @brief The total number of swimmers in the simulation.
int NS;
/// @brief The volume accessible to the simulation fluid. Determined by Monte Carlo.
double VOL;
/// @brief The particle number density of the simulation fluid. Found using the volume VOL determined by Monte Carlo.
double nDNST;
/// @brief The mass density of the simulation fluid. Found using the volume VOL determined by Monte Carlo. 
double mDNST;

/// @brief The x, y, z dimensions of the control volume.
int XYZ[3];
/// @brief The x, y, z dimensions of the control volume plus 1.
int XYZ_P1[3];
/// @brief Flags as to whether the x, y, or z dimensions are wrapped in periodic boundary conditions.
int XYZPBC[3];
/// @brief The maximum length that can fit in the control volume (the diagonal length).
int maxXYZ;

/// @brief The name of the input file for the MD simulation.
char* mdInputFile;

#endif
