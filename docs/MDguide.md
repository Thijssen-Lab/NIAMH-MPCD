# Guide to Run MD Coupled To MPCD
This is just meant to be a quick and short guide to MD + MPCD JSON input system.
As you can see in the examples provided (in the `sampleInputs/md2D` and `sampleInputs/md3D` folders, all is required is to include the MD input file path in JSON input file:
```json
    "mdIn":             "md.inp",
```
Almost all the parameters related to MD part can be set in md.inp file. The only MD parameter needed to be set in JSON input file is **stepsMD**. Pay attention to how **time steps** of MD and MPCD are coupled, because it is important for **stability** of MD simulation and **output frequency**. 

MD particles are coupled to MPCD particles during collision steps. For the sake of **stability** in the MD part, before entering the collision step ( during each streaming steps of MPCD ) MD particles should go through a number of md steps. This can be specified by changing **stepsMD** in JSON file (50 steps is safe!). Consequently, parameter **dt** in file md.inp should be changed as well, `dt_MPCD/stepsMD = dt_MD`. 

The same thing is applied to the frequency of data printed out. For instance, if the frequency for the director field output is every 100 streaming steps and we want the same number of data for MD particles, and there are 50 MD steps every MPCD streaming step, then all the numbers in **step counters** part in md.inp file should be 100 X 50 = 5000.

Furthermore, in **simulation phases** section, the third number in **nStep** specifies the time step up to which data will be printed out, after this time step there will not be any output from MD. 

The `sampleInputs/md2D` is a 2D system with periodic boundary condition. The `sampleInputs/md3D` is a cylinder, along the x direction, of liquid crystal with planar anchoring and slip condition at cyliner surface.
