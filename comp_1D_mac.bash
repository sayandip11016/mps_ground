#!/bin/bash
# script to call the Mac gnu fortran 95 compiler for program driver_mps_1D.f90
gfortran -O3 -llapack -lblas double_precision_mod.f90 matrix_mod.f90 vardefs_v5.f90 driver_mps_1D.f90 -o a_1D.o
