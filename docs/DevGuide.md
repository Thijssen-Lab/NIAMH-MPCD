# Guide for MPCD Developers

## Contents
1. [Introduction](#introduction)
2. [Code Style](#code-style)
   - [Comments](#comments)
   - [Hard-coded Constants](#hard-coded-constants)
   - [Dimensional Constants](#dimensional-constants)
   - [Global Variables](#global-variables)
3. [Implementing New Input/ Output Options](#implementing-new-input-output-options)
4. [Documentation](#documentation)

TODO

## Introduction
This document is intended to provide a rough guide for editors of the MPCD code. 
Here we specify some rough rules we would like others to follow when contributing to the code.

## Code Style
"[MPCD] is a scientific code, so rules are meant to be broken" --- Tyler Shendruk, 2023

MPCD is a code designed and used by physicists. As such, we do not enforce a strict code style --- The physics and research is _always_ the priority.
Despite this, remember that the code is used and developed by man people, so we need to keep it clean for each other. As such, there are some general rules we would like to follow:

### Comments
Comments should be used to explain the code broadly, and should not be used to explain every individual line of code.

Good uses for comments are:
- Explaining a block of code doing a particular task.
- Explaining the implementation of a complex physics or technical algorithm, or an equation.
- Adding a citation to explain the source for an equation/ algorithm.

Comments should also be avoided as a preface to methods, immediately post-declaration. These should instead be included as part of the method's [documentation](#documentation).

### Hard-coded Constants
Hard-coded numerical values should be avoided wherever possible. If you _need_ to use some hard-coded value, check `definitions.h` to see if it (or something similar) is already defined. 
If you still desperately need a new constant, add it as a `#define` pre-processor statement in `definitions.h`.

### Dimensional Constants
A particular code-tic that is used throughout the code, is when referring to a dimension (whether it be in an array, checking against a dimension, etc), then you should use the defined `_1D`, `_2D`, or `_3D` pre-processor defines. 
This is to ensure that it is clear we are referring to a dimension when handling data structures in the code. 

For example:
```c
double vec[_3D]; 
```
makes it clear that we are referring to and intend `vec` to be a 3 dimensional mathematical vector.

### Global Variables
Global variables should be avoided wherever possible. 
If following regular programming practices, you should _not_ require any global variable that isn't already defined. 
The list of global variables available throughout the code can be seen in `globals.h`.

## Implementing New Input/ Output Options
TODO

## Documentation
TODO

### Documentation Style Guide