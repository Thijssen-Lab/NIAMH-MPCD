# Guide for MPCD Developers

## Contents
1. [Introduction](#introduction-dev)
2. [Git/ Github](#git-github)
3. [Code Style](#code-style)
   - [Comments](#comments)
   - [Hard-coded Constants](#hard-coded-constants)
   - [Dimensional Constants](#dimensional-constants)
   - [Global Variables](#global-variables)
   - [Warnings](#warnings)
4. [Implementing New Input/ Output Options](#implementing-new-input-output-options)
   - [JSON Input](#json-input)
   - [Output Files](#output-files)
5. [Documentation & Documentation Style Guide](#documentation--documentation-style-guide)
   - [File Example](#file-example)
   - [Method Example](#method-example)
   - [Global Variable Example](#global-variable-example)
   - [Pre-Procesor Example](#pre-procesor-example)
   - [Struct Example](#struct-example)

## Introduction      {#introduction-dev}
"[MPCD] is a scientific code, so rules are meant to be broken" --- Tyler Shendruk, 2023

This is a scientific code, so we expect research to take the priority, however we wish that all users are considerate to each other. 
This is the primary purpose of this document: To set some guidelines for how we can all work together to make the code as clean and easy to use as possible.

## Git/ Github       {#git-github}
MPCD is hosted on Github, and uses Git for version control. 
We put this first and foremost to ensure everyone understands the normal process for contributing to the code:
1. Either make a new branch prefixed by your initials (for Shendruk lab members), or make a fork to your personal Github account (for all others).
2. Make your changes, and commit them to your branch/ fork.
3. Check your code!
   - Ensure your code compiles, runs as expected, and does not break any existing functionality (check via the sample inputs).
   - Ensure your code satisfies all guidelines in this document.
4. Make a pull request (PR) to the main repository on Github, and assign it to a relevant Shendruk lab member for review. 
   - If possible, contact them via Slack or email to inform them.
5. If you are asked to make corrections as part of your PR, correct them and ask for a re-review. Repeat this process until your changes are approved.
6. Once approved, your code will be merged into the main repository.

## Code Style        {#code-style}
MPCD is a code designed and used by physicists. As such, we do not enforce a strict code style --- The physics and research is _always_ the priority.
Despite this, remember that the code is used and developed by man people, so we need to keep it clean for each other. As such, there are some general rules we would like to follow:

### Comments         {#comments}
Comments should be used to explain the code broadly, and should not be used to explain every individual line of code.

Good uses for comments are:
- Explaining a block of code doing a particular task.
- Explaining the implementation of a complex physics algorithm, an equation, or a particularly technical operation.
- Adding a citation to explain the source for an equation/ algorithm.

Comments should also be avoided as a preface to methods, immediately post-declaration. These should instead be included as part of the method's [documentation](#documentation--documentation-style-guide).

### Hard-coded Constants         {#hard-coded-constants}
Hard-coded numerical values should be avoided wherever possible. If you _need_ to use some hard-coded value, check `definitions.h` to see if it (or something similar) is already defined. 
If you still desperately need a new constant, add it as a `#define` pre-processor statement in `definitions.h`.

### Dimensional Constants        {#dimensional-constants}
A particular code-tic that is used throughout the code, is when referring to a dimension (whether it be in an array, checking against a dimension, etc), then you should use the defined `_1D`, `_2D`, or `_3D` pre-processor defines. 
This is to ensure that it is clear we are referring to a dimension when handling data structures in the code. 

For example:
``````
double vec[_3D]; 
``````
makes it clear that we are referring to and intend `vec` to be a 3 dimensional mathematical vector.

### Global Variables       {#global-variables}
Global variables should be avoided wherever possible. 
If following regular programming practices, you should _not_ require any global variable that isn't already defined. 
The list of global variables available throughout the code can be seen in `globals.h`.

### Warnings         {#warnings}
Try to avoid introducing new warnings into the code, when compiled with the `make` command in the root directory.

## Implementing New Input/ Output Options       {#implementing-new-input-output-options}
There may sometimes be a need to implement either new JSON input options for the code, or to add new output files. 
For both of these, there is a specific process that you must follow in order to ensure consistency between files & inputs.

### JSON Input       {#json-input}
Adding new JSON inputs is the most common of these two tasks you may need to do. 
This is documented in detail at the end of `/Docs/InputGuide.md`, and the process will not be repeated here to avoid duplication. 

A summary of the key parts are:
1. Ensure the variable you wish to read data into can be accessed by `readJson()` in `read.c`.
2. Use one of the cJSON read API methods (detailed in both the input guide and `cJson.c`) to read the data into the variable. 
   - Place it, so it is with all the other reads in `readJson()`, rather than at the very beginning or at the very end.
   - Arrays, in particular arrays of custom objects, are more complicated. Copy one of the existing arrays in the code and modify it appropriately.
3. **Update the `InputGuide.md` file!** This is especially important, as it is the only concise guide to all inputs!

### Output Files        {#output-files}
New output files should use the same header as all other output files. 
All non-integer numerical values should be written using scientific notation. 
Columns should be labelled with intuitive variable names and separated by tabs (both header and data). 
- Fields should redundantly write the time step on each position (with no empty line between time steps), and should iterate through z then y then x coordinates. Output should write 3D even for simulations performed in 2D. 
- System-wide measurements should write the time on each line. 
- Histograms should write the time on its own line followed by the bins and associated counts. 

## Documentation & Documentation Style Guide       {#documentation--documentation-style-guide}
The code uses Doxygen to generate API documentation for the code. 
As such, each of the following should be documented:
- All `.c` files, and any `.h` file that is not purely function pre-declarations.
- All methods
- All global variables (these are primarily in `globals.h`)
- All (non-header guard) pre-processor defines (these are primarily in `definitions.h`)
- All structs (these are primarily in `SRDclss.h`)

Examples, and any key points, for each of these are included below.

The key points used throughout are:
- We use **only** `///` for documentation.
- Everything, except global variables/ member variables, are required to have an `@brief`.
- There should be an empty "doxygen line" (ie `///`) at the beginning, end, and between `@file`, `@brief`, and details, and parameters.

### File Example        {#file-example}
``````
///
/// @file
///
/// @brief Prints data output files.
///
/// A collection of functions for constructing and printing different raw data outputs to .dat files. The types of files produced must be specified in input.json.
///
``````
- Empty doxygen line between `@file`, `@brief`, and the detailed explanation.

### Method Example         {#method-example}
``````
        ...
    }
}

///
/// @brief	Determines if the BC must be shifted due to the periodicity of the control volume.
///   
/// It checks if the BC is periodic, then it calculates the shift and shifts the BC.
/// 
/// @param shift	This is how much the boundary must be shifted, gets calculated inside the routine.
/// @param WALL		One of the walls of the BCs that particle is interacting with.
/// @param atom 	The MD particle.
/// @see			shiftBC()
///
void shiftBC_MD( double *shift,bc *WALL,particleMD *atom ) {
    ...
}
``````
- Empty line between previous method, and doxygen comment for next method
- `@param` for each method argument, in the order they appear
- `@return` if the method returns a value
- `@see` if there is a method that is strongly related

### Global Variable Example         {#global-variable-example}
``````
/// @brief The dimension of the simulation. Must be 1, 2, or 3.
int DIM;
``````
- `@brief` on the line preceding the variable declaration.

### Pre-Processor Defines Example         {#pre-procesor-example}
``````
/// @brief Number of bins used for distributions. Best if an odd integer.
# define BINS 101
``````
- `@brief` on the line preceding the define declaration.

### Struct Example         {#struct-example}
``````
    ...
} specSwimmer;

///
/// @brief Struct that represents a single monomer.
///
/// This container structure is used to store all of the parameters associated with a single monomer.
///
typedef struct smono {
	double Q[_3D];				///< Position of the monomer.
	double V[_3D];				///< Velocity of the monomer.
	double A[_3D];				///< Acceleration of the monomer.
	int HorM;					///< Signfies whether monomer is a head, or a middle (0=head; 1=middle). Tail takes the mass of the head and doesn't participate.
	struct smono *next;			///< Pointer to next particle in cell list.
	struct smono *previous;		///< Pointer to previous particle in cell list.
} smono;
``````
- Empty line between previous struct, and doxygen comment for next struct.
- Each member variable should be documented inline with `///<` - the `<` is important!