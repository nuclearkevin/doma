# DOMA
**D**iscrete **O**rdinates **M**ini-**A**pp.

A multigroup discrete ordinates radiation transport solver for 1D/2D/3D orthogonal cartesian grids. This started as project for MCSC 6020U at Ontario Tech; it is under continuous improvement for my final projects in NPRE 555 and NPRE 560 at UIUC. This toy transport solver currently implements:

- Transient calculations using the transient fixed source method (constant delayed neutron precursor fission source between timesteps);
- Power iteration for criticality eigenvalue (k_eff) calculations;
- Gauss-Seidel iteration for converging the in-scattering and fission matrices;
- Scattering source iteration;
- Parallel-in-angle naive transport sweeps;
- Theta-weighted upwinded diamond differences for spatial discretization;
- Arbitrary cartesian meshes with an input syntax similar to MOOSE's CartesianMeshGenerator.

## Dependencies

- OpenMP (included by default in GCC 4.2.0). If on Ubuntu and using older GCC compilers: `sudo apt install libomp-dev`.

## Sample problems

<p align="center">
  <img src="/images/2D_g1_lin.png" width="200" />
  <img src="/images/2D_g2_lin.png" width="200" />
  <img src="/images/2D_g3_lin.png" width="200" />
</p>

<p align="center">
Eigenfluxes in a 2D, 3 group subcritical assembly: cases/2D_assembly_eig.xml
</p>

<p align="center">
  <img src="/images/3D_g1_lin.png" width="200" />
  <img src="/images/3D_g2_lin.png" width="200" />
  <img src="/images/3D_g3_lin.png" width="200" />
</p>

<p align="center">
Eigenfluxes in a 3D, 3 group subcritical assembly: cases/3D_assembly_eig.xml
</p>

<p align="center">
  <img src="/images/3D_src_g1_log.png" width="200" />
  <img src="/images/3D_src_g2_log.png" width="200" />
  <img src="/images/3D_src_g3_log.png" width="200" />
</p>

<p align="center">
Source-driven fluxes in a 3D, 3 group subcritical assembly: cases/3D_assembly_fs.xml
</p>

<p align="center">
  <img src="/images/3D_kob_3.png" width="600" />
</p>

<p align="center">
Scalar flux in the fixed-source Kobayashi benchmark (problem 3): cases/3D_kobayashi_3.xml
</p>
