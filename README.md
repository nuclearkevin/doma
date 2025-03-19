# DOMA
**D**iscrete **O**rdinates **M**ini-**A**pp.

A multigroup discrete ordinates radiation transport solver for 1D/2D/3D orthogonal cartesian grids. This started as project for MCSC 6020U at Ontario Tech; it is under continuous improvement for my final projects in NPRE 555 and NPRE 560 at UIUC. This toy transport solver currently implements:

- Transient calculations using the transient fixed source method (constant delayed neutron precursor fission source between timesteps);
- Power iteration for criticality eigenvalue (k_eff) calculations;
- Gauss-Seidel iteration for converging the in-scattering and fission matrices;
- Scattering source iteration;
- Naive (non-wavefront) serial transport sweeps;
- Theta-weighted upwinded diamond differences for spatial discretization;
- Arbitrary cartesian meshes with an input syntax similar to MOOSE's CartesianMeshGenerator.

<p align="center">
  <img src="/images/2D_g1_lin.png" width="200" />
  <img src="/images/2D_g2_lin.png" width="200" />
  <img src="/images/2D_g3_lin.png" width="200" />
</p>

<p align="center">
  <img src="/images/3D_g1_lin.png" width="200" />
  <img src="/images/3D_g2_lin.png" width="200" />
  <img src="/images/3D_g3_lin.png" width="200" />
</p>
