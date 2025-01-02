# DOMA
**D**iscrete **O**rdinates **M**ini-**A**pp.

A multigroup discrete ordinates radiation transport solver for 1D/2D/3D orthogonal cartesian grids. This started as my final project for MCSC 6020U at Ontario Tech, and is now being improved for my final projects in NPRE 555 (soon) and NPRE 560 at UIUC. This toy transport solver currently implements:

- Transient calculations using the transient fixed source method (constant delayed neutron precursor fission source between timesteps);
- Power iteration for criticality eigenvalue (k_eff) calculations;
- Gauss-Seidel iteration for converging the in-scattering and fission matrices;
- Scattering source iteration;
- Naive (non-wavefront) serial transport sweeps;
- (Theta-Weighted) Diamond Differencing spatial discretization (with upwinding);
- Arbitrary cartesian meshes with an input syntax similar to MOOSE's CartesianMeshGenerator.
