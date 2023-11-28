# DOMA
**D**iscrete **O**rdinates **M**ini-**A**pp.

A grey radiation transport solver for 3D orthogonal cartesian grids, my final project for MCSC 6020U at Ontario Tech. Currently has:

- Diamond Differencing and Step Characteristics spatial discretization;
- Naive (non-wavefront) serial transport sweeps;
- Scattering source iteration;
- Arbitrary orthogonal cartesian meshes with an input syntax similar to MOOSE's CartesianMeshGenerator.
