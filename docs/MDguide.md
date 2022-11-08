# Guide to Run MD Coupled To MPCD
This is just meant to be a quick and short guide to MD + MPCD JSON input system.
As you can see in the examples provided (in the `sampleInputs/md2D` and `sampleInputs/md3D` folders, all is required is to include the MD input file in JSON input file:
```json
    "mdIn":             "md.inp",
```
All the parameters related to MD part can be changed in md.inp file. the only **important** point needed to be considered is how the **timesteps** are coupled.
MD particles are coupled to MPCD particles during collision steps. For the sake of **stability** in the MD part, before entering the collision step ( during each streaming steps of MPCD ) MD particles should go through a number of mdsteps. This can be specified by changing **"stepsMD"** in JSON file (50 steps is safe!). Consequently parameter **dt** in file md.inp should be changed as well, dt_MPCD/stepsMD = dt_MD.
The same thing is applied to the frequency of getting data printed out. For instance, if the frequency for the director field output is every 100 streaming steps and we want the same number of data for MD particles, and there are 50 MD steps every MPCD streaming step, then all the numbers in **step counters** part in md.inp file should be 100 X 50 = 5000 . 

The 2D folder is a 2D system with periodic boundary condition. The 3D folder is a cylinder, along the x direction, of liquid crystal with planar anchoring and slip condition at cyliner surface.
