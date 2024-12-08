# DOMA
**D**iscrete **O**rdinates **M**ini-**A**pp.

A multigroup discrete ordinates radiation transport solver for 1D/2D/3D orthogonal cartesian grids. This started as my final project for MCSC 6020U at Ontario Tech, and is now being improved for my final projects in NPRE 560 and NPRE 555. Currently has:

- Diamond Differencing and Theta-Weighted Diamond Differencing spatial discretization (with upwinding);
- Naive (non-wavefront) serial transport sweeps;
- Scattering source iteration;
- Gauss-Seidel iteration for converging the in-scattering and fission matrices;
- Arbitrary orthogonal cartesian meshes with an input syntax similar to MOOSE's CartesianMeshGenerator.
