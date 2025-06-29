# NIAMH-MPCD
This is the NoIsy Algorithm for Mesoscale Hydrodynamics Multi-Particle Collision Dynamics (NIAMH-MPCD) code developed by the [Shendruk Lab](https://tnshendruk.com). The primary distribution for the simulator is available at [Shendruk-Lab/NIAMH-MPCD](https://github.com/Shendruk-Lab/NIAMH-MPCD).

If you just have obtained the code, and you are looking for documentation on how to use it, you might first want to have a look at the User's Guide at [docs/userguide/html](https://github.com/Shendruk-Lab/NIAMH-MPCD/docs/userguide/html).

## Useful Links
- [Shendruk Lab website](https://tnshendruk.com)
- [Main NAIMH-MPCD repository](https://github.com/Shendruk-Lab/NIAMH-MPCD)
- [Code documentation](https://github.com/Shendruk-Lab/NIAMH-MPCD/docs/DevGuide.md)
- [User's guide](https://github.com/Shendruk-Lab/NIAMH-MPCD/docs/userguide/html)

## Make
After downloading, build the code by running `make` in the main directory. 
Other make options are detailed in the `Makefile`.

Build code documentation using `make docs`, and view them in `docs/doxyout/html/index.html`.

## MPCD
Multi-Particle Collision Dynamics is a particle-based mesoscale simulation technique for complex fluids in low Reynolds number, which fully incorporates thermal fluctuations and hydrodynamic interactions. This code can be used to simulate systems such as nematic liquid crystals, active and passive, with inclusions such as polymers, swimmers and complex mobile boundaries in 2D and 3D. Details on the MPCD approach can be found in the following review articles:

* [Simulation of complex fluids by multi-particle-collision dynamics](https://www.sciencedirect.com/science/article/pii/S0010465505001700)

* [Mesoscopic model for solvent dynamics](https://aip.scitation.org/doi/abs/10.1063/1.478857)

* [Multi-particle collision dynamics algorithm for nematic fluids](https://pubs.rsc.org/en/content/articlehtml/2015/sm/c5sm00839e)

* [Modeling hydrodynamic interactions in soft materials with multiparticle collision dynamics](https://www.sciencedirect.com/science/article/pii/S2211339819300024) 

## Credits
Code started and majority written by Tyler Shendruk ([email](mailto:t.shendruk@ed.ac.uk)), for the Polymer Physics Research Group at the University of Ottawa.

Originally started in October 2008 and known as Stochastic Rotation Dynamics.

The MD simulator used was written by Frédéric Tessier, started in 2005 at the University of Ottawa.

Major contributors to the code include Louise Head, Timofey Kozhukhov, Zahra Valei, and Benjamin Loewe at the University of Edinburgh.

Uses the [cJson library](https://github.com/DaveGamble/cJSON) by Dave Gramble, and the corresponding [parser](https://github.com/T-Kozhukhov/cJson-Parser) by Timofey Kozhukhov under MIT license.