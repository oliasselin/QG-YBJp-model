QG-YBJ+ model
=============

This is a numerical model for the two-way interaction of near-inertial waves with (Lagrangian-mean) balanced eddies. Wave evolution is governed by the YBJ+ equation (Asselin & Young 2019). The traditional quasigeostrophic equation dictates the evolution of potential vorticity, which includes the wave feedback term of Xie & Vanneste (2015). The model is pseudo-spectral in the horizontal and uses second-order finite differences to evaluate vertical and time derivatives. 

Code written by Olivier Asselin


Brief overview of files
=======================

#Essentials


parametersXXX.f90: contains all the parameters determining the simulation.

init.f90:          initialization of all basic arrays, stratification profile, initial condition for eddies and waves.

IO_ncf.f90:        all things netCDF input/output.

lcometXXX          compiling and job launching script

main_waqg.f90:     main program performing the integration



#Under the hood


elliptic.f90:      routines pertaining to inversion of q for psi, and LA for A. 

derivatives.f90:   contains various subroutines computing derivatives and nonlinear terms via the transform method.

fft.f90            all things Fourier transforms

mpi.f90            all things parallelization via MPI



#Deprecated


diagnostic.f90:    contains a bunch of various old diagnostics (obsolete)

files.f90:         initialize all text files needed (obsolete)

special.f90:       contains a couple special functions for diagnostics (obsolete)
