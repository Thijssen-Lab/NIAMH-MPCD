# Guide to MPCD JSON Input

This is just meant to be a quick and dirty guide to the new MPCD JSON input system. It is not meant to be a complete reference. For basics on the JSON file format take a look at the tutorial series starting on [this page](https://www.w3schools.com/js/js_json_intro.asp).

**Note:** The MPCD JSON input system has no bindings for MD. Please use legacy input for MD sims for now.

Everything listed below should be in a single JSON file. 
Every parameter in the simulation can be fed in with a corresponding name/value pair. 
The name, or "tag", uniquely represents each parameter, is case sensitive, but can be given in any order within the json file. 
If a tag is not specified, then that parameter will assume the default value in the table below --- This means that, with the exception of some parameters used when declaring a BC, you may only need to specify several parameters in your input file.

The table below lists all tags, their types, their default values, and their description.
Tags will generally be named the same as they are named in the code (in the various input structs), unless there is a good reason to change it for legibility or clarity.

Some tags will take an array of custom objects (for example, species and boundary conditions). 
These are each defined in their own table following the main one.

## Main Input Tag Table

All tags suffixed with "out", unless otherwise specified, take a value representing how frequently (in timesteps) you want that quantity to be dumped to file.

Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`domain`        | array(int)    | [30, 30]      | The system size, given as [X, Y] or [X, Y, Z]. Array dimensions correspond to DIM. If run in 1D or D>3, sim won't run.
`kbt`           | double        | 1             | The thermal energy of the system. Should ALWAYS be 1. 
`dt`            | double        | 0.1           | The size of the simulation timestep
`simSteps`      | int           | 2000          | How many total timesteps to run the simulation for
`warmUp`        | int           | 0             | How much warmup time to run before starting the simulation. Note that during warmup no output is written to file
`rFrame`        | int           | 1             | Takes of net non-zero average momentum post initialisation to minimise drift. 1 = on, 0 = off
`zeroNetMom`    | int           | 0             | This substracts any excess momentum every time step (generally have this off)
`galInv`        | int           | 1             | Whether to enable the random shift of particles to counter gallilean invariance
`tsTech`        | int           | 0             | Which thermostat technique to use. Note some collsion operators have thermostats built in. See definitions.h for list
`rTech`         | int           | 2             | Which collision operator to use. See definitions.h for list
`lc`            | int           | 0             | Liquid crystal mode. 2 = global S, 1 = local S, 0 = off
`tau`           | double        | 0.5           | Thermal relaxation time scale
`rotAng`        | double        | 1.570796      | This is the angle used in the original SRD collision operator
`fricCoef`      | double        | 1.0           | Friction coefficient for langevin thermostat
`mfpot`         | double        | 10            | Liquid crystal mean field potential
`grav`          | array(double) | [0,0,0]       | Constant acceleration due to external force. MUST be 3D
`mag`           | array(double) | [0,0,0]       | Constant external magnetic field. MUST be 3D
`seed`          | int           | 0             | Seed for random number generator. 0 for pseudorandom seed. Set to -1 to load a checkpoint.
`mdMode`        | int           | 0             | Enable the MD sim coupling. 1 = on, 0 = off
`stepsMD`       | int           | 20            | MD time steps per MPCD time step
`species`       | array(species)| 1 default spec  | An array of species objects. See the species table for species tags.
---             | ---           | ---           | ---
`debugOut`      | int           | 3             | Debug (verbosity) level. See definitions.h for list.
`trajOut`       | int           | 0             | Detailed species trajectories
`trajSpecOut`   | int           | 0             | Which number of species whose detailed trajectories to output
`coarseOut`     | int           | 0             | Coarse grain data
`flowOut`       | int           | 0             | Flow field
`avVelOut`      | int           | 0             | Total average MPCD velocity
`dirSOut`       | int           | 0             | Director and scalar order parameter fields
`qTensOut`      | int           | 0             | Q tensor field
`qkTensOut`      | int           | 0             | Reciprocal Q tensor field
`oriEnOut`      | int           | 0             | Orientational energy field
`colourOut`     | int           | 0             | Colour/ phi/ species-type field
`pressureOut`   | int           | 0             | Pressure field
`neighbourEnOut`| int           | 0             | Orientational energy from neighbours
`avSOut`        | int           | 0             | Total average scalar order parameter
`densSDOut`     | int           | 0             | SD of the number per cell
`enstrophyOut`  | int           | 0             | Total average enstrophy
`histVelOut`    | int           | 0             | Velocity distribution
`histSpeedOut`  | int           | 0             | Speed distribution
`histVortOut`   | int           | 0             | Vorticity distribution
`histEnsOut`    | int           | 0             | Enstrophy distribution
`histDirOut`    | int           | 0             | Director distribution
`histSOut`      | int           | 0             | Scalar order parameter distribution
`histNOut`      | int           | 0             | Number per cell distribution
`solidTrajOut`  | int           | 0             | Solid BC trajectories
`topoFieldOut`  | int           | 0             | Topological charge field
`energyOut`     | int           | 0             | System energy
`velCorrOut`    | int           | 0             | Velocity autocorrelation
`dirCorrOut`    | int           | 0             | Director autocorrelation
`vortCorrOut`   | int           | 0             | Vorticity autocorrelation
`densCorrOut`   | int           | 0             | Densty autocorrelation
`orderCorrOut`  | int           | 0             | Scalar order parameter autocorrelation
`phaseCorrOut`  | int           | 0             | Binary fluid correlation
`energySpecOut` | int           | 0             | Energy spectra
`enstrophySpecOut`| int           | 0             | Enstrophy spectra
`binderOut`     | int           | 0             | Binder cumulant
`binderBin`     | int           | 0             | Binder cumulant bin size
`swimQOut`      | int           | 0             | Swimmer positions
`swimOOut`      | int           | 0             | Swimmer orientations
`swimROut`      | int           | 0             | Swimmer run/ tumble
`synopsisOut`   | int           | 1             | Synopsis output. Highly recommended to be on. 1 = on, 0 = off
`checkpointOut` | int           | 0             | Simulation checkpointing
---             | ---           | ---           | ---
`BC`            | array(BC)     | PBCs around domain | The array of boundary objects. See the BC table for BC tags.
---             | ---           | ---           | ---
`typeSwim`      | int           | 2             | Swimmer type. Default is a dumbell swimmer with excluded volume interactions
`nSwim`         | int           | 0             | Swimmer population
`qDistSwim`     | int           | 0             | Initial distribution of swimmer positions
`oDistSwim`     | int           | 0             | Initial distribution of swimmer orientations
`headMSwim`     | int           | 20            | Mass of head monomer
`midMSwim`      | int           | 20            | Mass of middle monomer
`hspIdSwim`     | int           | 1             | Multiphase fluid particle type of the head monomer in the swimmer
`mspIdSwim`     | int           | 1             | Multiphase fluid particle type of the middle monomer in the swimmer
`fsSwim`        | double        | 20            | Magnitude of swimmer propulsion force
`dsSwim`        | double        | 1             | Dipole strength
`tsSwim`        | double        | 0             | Magnitude of swimmer torque
`sizeShrinkSwim`| double        | 0.1           | How much do LJ sigma ro shrink when tumbling
`springShrinkSwim`| double        | 0.1           | How much the spring constant is shrunk when tumbling
`kSwim`         | double        | 30            | Spring constant
`roSwim`        | double        | 4             | Spring seperation
`sigSwim`       | double        | 4             | Diameter approx sigma
`epsSwim`       | double        | 1             | Interaction energy
`runTSwim`      | double        | 0             | Average run time in units of MPCD timesteps dt
`tumTSwim`      | double        | 0             | Average tumble time in units of MPCD timesteps dt
`shrTSwim`      | double        | 2             | Set time to shrink/ extend in units of MPCD timesteps dt
`magMomSwim`    | double        | 1             | Magnetic moment/ strength
`fixDistSwim`   | double        | 0             | The fixed distance from the wall for DUMBELL_NEARWALL mode


## Species Tag Table
Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`mass`          | double        | 1             | Mass of this species of particles  
`pop`           | int           | 18000         | Number of particles of this species
`qDist`         | int           | 0             | Positional distribution function for the species. See definitions.h for a list.
`vDist`         | int           | 0             | Velocity distribution function for the species. See definitions.h for a list.
`oDist`         | int           | 2             | Orientation distribution function for the species. See definitions.h for a list.
`interMatr`     | array(double) | [0,...]       | Interaction matrix for this species, against other species. Must be of the same length as the number of species. Default will autopopulate with 0s
`rfc`           | double        | 0.01          | Nematogen rotational friction coefficient
`len`           | double        | 0.007         | Effective rod length (specifically for solid boundary interactions)
`tumble`        | double        | 2             | Tumbling parameter
`shearSusc`     | double        | 0.5           | Shear susceptibility
`magnSusc`      | double        | 0.001         | Magnetic susceptibility
`act`           | double        | 0.05          | Species activity
`damp`          | double        | 0             | Damping friction to kill hydrodynamics. Between 0 and 1.

## BC Tag Table

Unlike the previous tags, if you declare a new BC object then there are some tags that MUST be declared. 
These are noted as NECESSARY in the default value column.

Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`colType`       | int           | 1             | Which BC type do you want to use. See definitions.h for a list.
`phantom`       | int           | 0             | Use phantom particles if 1, 0 otherwise
`E`             | double        | -1            | Coefficient of restitution
`Q`             | array(double) | [0,0,0]       | Position of the boundary, MUST be 3D
`V`             | array(double) | [0,0,0]       | Velocity of the boundary, MUST be 3D
`O`             | array(double) | [0,0,0]       | Orientation of the boundary about (x,y,z), MUST be 3D
`L`             | array(double) | [0,0,0]       | Angular velocity of the boundary, MUST be 3D
`G`             | array(double) | [0,0,0]       | External acceleration (ie, due to gravity) of the boundary, MUST be 3D
`aInv`          | array(double) | NECESSARY     | Sets the geometry of the surface. Principal semi-axes of the ellipsoid (see SRDClass for explanation). MUST be 3D. If you declare a BC then this MUST be given, or the simulation will not run
`rotSym`        | array(double) | [4,4]         | Sets rotational symmetry of shapes. Must be of form (a,b)
`abs`           | int           | 0             | Flags if each term should be absolute only. 1 = yes, 0 = no
`P`             | array(double) | NECESSARY     | Sets the geometry of the surface. Must be of form (a,b,c,d). See SRDClass for explanation. If you declare a BC then this MUST be given, or the simulation will not run
`R`             | double        | NECESSARY     | Sets the geometry of the surface. See SRDClass for explanation. If you declare a BC then this MUST be given, or the simulation will not run
`DN`            | double        | NECESSARY     | Displacement of the particle in the normal direction. If you declare a BC then this MUST be given, or the simulation will not run
`DT`            | double        | 0             | Displacement of the particle in the tangential direction
`DVN`           | double        | 0             | Add velocity in the normal direction on contact
`DVT`           | double        | 0             | Add velocity in the tangential direction on contact
`DVxyz`         | array(double) | [0,0,0]       | Add velocity in the cartesian directions on contact. Must be 3D
`MVN`           | double        | NECESSARY     | Multiplies the velocity in the normal direction on contact. If you declare a BC then this MUST be given, or the simulation will not run
`MVT`           | double        | NECESSARY     | Multiplies the velocity in the tangential direction on contact. If you declare a BC then this MUST be given, or the simulation will not run
`MUN`           | double        | 1             | Multiplies the orientation in the normal direction on contact
`MUT`           | double        | 1             | Multiplies the orientation in the tangential direction on contact
`MUxyz`         | array(double) | [1,1,1]       | Multiplies the orientation in the cartesian direction on contact. Must be 3D
`DUxyz`         | array(double) | [0,0,0]       | Add the orientation in the cartesian direction on contact. Must be 3D
`kbt`           | double        | 1             | Temperature of the wall
`dsplc`         | int           | 0             | Whether the wall can displace/ is mobile. 0 = no, 1 = yes
`inv`           | int           | 0             | Whether to invert the bc (ie, multiply the A's by -1). 0 = no, 1 = yes
`mass`          | double        | 1             | Mass of the wall in MPCD units. Should be the same density as the fluid if its displaceable

## Example
Below is an example JSON input file that uses all tags but sets them to the defaults.
As a reminder, if you wish to use the default value for a tag, you can leave it out of the JSON file, so this is primarily for illustrative purposes.

```json
{

}
```