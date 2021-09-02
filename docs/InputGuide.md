# Guide to MPCD JSON Input

This is just meant to be a quick and dirty guide to the new MPCD JSON input system. It is not meant to be a complete reference. For basics on the JSON file format take a look at the tutorial series starting on [this page](https://www.w3schools.com/js/js_json_intro.asp).

Everything listed below should be in a single JSON file. 
Every parameter in the simulation can be fed in with a corresponding name/value pair. 
The name, or "tag", uniquely represents each parameter, is case sensitive, but can be given in any order within the json file. 

The table below lists all tags, their types, their default values, and their description.
Tags will generally be named the same as they are named in the code (in the `inputList` struct) unless they are completely nonsensical.
Some tags will take an array of custom objects (for example, species and boundary conditions). 
These are each defined in their own table following the main one.

## Main Input Tag Table

All tags suffixed with "out", unless otherwise specified, take a value representing how frequently (in timesteps) you want that quantity to be dumped to file.

Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`domain`        | array(int)    |               | The system size, given as [X, Y] or [X, Y, Z]. Array dimensions correspond to DIM.
`kbt`           | double        |               | The thermal energy of the system
`dt`            | double        |               | 
`simSteps`      | int           |               | 
`warmUp`        | int           |               | 
`rFrame`        | int           |               | 
`zeroNetMom`    | int           |               | 
`galInv`        | int           |               | 
`tsTech`        | int           |               | 
`rTech`         | int           |               | 
`lc`            | int           |               | 
`tau`           | double        |               | 
`rotAng`        | double        |               | 
`fricCoef`      | double        |               | 
`mfp`           | double        |               | 
`grav`          | array(double) |               | 
`mag`           | array(double) |               | 
`seed`          | int           |               | 
`mdMode`        | int           |               |
`stepsMD`       | int           |               |
`species`       | array(species)| default spec  | An array of species objects. See the species table for species tags.
---             | ---           | ---           | ---
`debugOut`      | int           |               | 
`trajOut`       | int           |               | 
`trajSpecOut`   | int           |               | 
`coarseOut`     | int           |               | 
`flowOut`       | int           |               | 
`avVelOut`      | int           |               | 
`dirSOut`       | int           |               | 
`qTensOut`      | int           |               | 
`oriEnOut`      | int           |               | 
`colourOut`     | int           |               | 
`pressureOut`   | int           |               | 
`neighbourEnOut`| int           |               | 
`avSOut`        | int           |               | 
`densSDOut`     | int           |               | 
`enstrophyOut`  | int           |               | 
`histVelOut`    | int           |               | 
`histSpeedOut`  | int           |               | 
`histVortOut`   | int           |               | 
`histEnsOut`    | int           |               | 
`histDirOut`    | int           |               | 
`histSOut`      | int           |               | 
`histNOut`      | int           |               | 
`solidTrajOut`  | int           |               | 
`topoFieldOut`  | int           |               | 
`energyOut`     | int           |               | 
`velCorrOut`    | int           |               | 
`dirCorrOut`    | int           |               | 
`vortCorrOut`   | int           |               | 
`densCorrOut`   | int           |               | 
`orderCorrOut`  | int           |               | 
`phaseCorrOut`  | int           |               | Binary fluid correlation
`energySpecOut` | int           |               | 
`enstrophySpecOut`| int           |               | 
`binderOut`     | int           |               | 
`binderBin`     | int           |               | 
`swimQOut`      | int           |               | 
`swimOOut`      | int           |               | 
`swimROut`      | int           |               | Swimmer run/ tumble
`synopsisOut`   | int           |               | Synopsis output. Highly recommended to be on. 1 = on, 0 = off.
`checkpointOut` | int           |               | 
---             | ---           | ---           | ---
`BC`            | array(BC)     |               | The array of boundary objects. See the BC table for BC tags.
---             | ---           | ---           | ---
`typeSwim`      | int           |               | 
`nSwim`         | int           |               | 
`qDistSwim`     | int           |               | 
`oDistSwim`     | int           |               | 
`headMSwim`     | int           |               | 
`midMSwim`      | int           |               | 
`hspIdSwim`     | int           |               | 
`mspIdSwim`     | int           |               | 
`fsSwim`        | int           |               | 
`dsSwim`        | int           |               | 
`dsSwim`        | double        |               | 
`sizeShrinkSwim`| double        |               | 
`springShrinkSwim`| double        |               | 
`kSwim`         | double        |               | 
`roSwim`        | double        |               | 
`sigSwim`       | double        |               | 
`epsSwim`       | double        |               | 
`runTSwim`      | double        |               | 
`tumTSwim`      | double        |               | 
`shrTSwim`      | double        |               | 
`magMomSwim`    | double        |               | 
`fixDistSwim`   | double        |               | 


## Species Tag Table
Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`mass`          | double        |               |   
`pop`           | int           |               | Number of particles of this species
`qDist`         | int           |               | Positional distribution function for the species. See definitions.h for a list.
`vDist`         | int           |               | Velocity distribution function for the species. See definitions.h for a list.
`oDist`         | int           |               | Orientation distribution function for the species. See definitions.h for a list.
`interMatr`     | array(double) |               | Interaction matrix for this species, against other species. Must be of the same length as the number of species.
`rfc`           | double        |               | 
`len`           | double        |               | 
`tumble`        | double        |               | 
`shearSusc`     | double        |               | Shear susceptibility
`magnSusc`      | double        |               | Magnetic susceptibility
`act`           | double        |               | Species activity
`damp`          | double        |               | Damping friction to kill hydrodynamics. Between 0 and 1.

## BC Tag Table
Tag             | Type          | Default Value | Description
---             | ---           | ---           | ---
`colType`       | int           |               | 
`phantom`       | int           |               | Use phantom particles if 1, 0 otherwise
`E`             | double        |               | Coefficient of restitution
`Q`             | array(double) |               | Position of the boundary, MUST be 3D
`V`             | array(double) |               | Velocity of the boundary, MUST be 3D
`O`             | array(double) |               | Orientation of the boundary about (x,y,z), MUST be 3D
`L`             | array(double) |               | Angular velocity of the boundary, MUST be 3D
`G`             | array(double) |               | External acceleration (ie, due to gravity) of the boundary, MUST be 3D
`G`             | array(double) |               | External acceleration (ie, due to gravity) of the boundary, MUST be 3D
`aInv`          | array(double) |               | Sets the geometry of the surface. Principal semi-axes of the ellipsoid (see SRDClass for explanation). MUST be 3D. 
`rotSym`        | array(double) |               | Sets rotational symmetry of shapes. Must be of form (a,b)
`abs`           | int           |               | 
`P`             | array(double) |               | Sets the geometry of the surface. Must be of form (a,b,c,d). See SRDClass for explanation
`R`             | double        |               | Sets the geometry of the surface. See SRDClass for explanation
`DN`            | double        |               | 
`DT`            | double        |               | 
`DVN`           | double        |               | 
`DVT`           | double        |               | 
`DVxyz`         | array(double) |               | Must be 3D
`MVN`           | double        |               | 
`MVT`           | double        |               | 
`MUN`           | double        |               | 
`MUT`           | double        |               | 
`MUxyz`         | array(double) |               | Must be 3D
`DUxyz`         | array(double) |               | Must be 3D
`kbt`           | double        |               | 
`dsplc`         | int           |               | 
`inv`           | int           |               | Whether to invert the bc. 0 = no, 1 = yes
`mass`          | double        |               | 