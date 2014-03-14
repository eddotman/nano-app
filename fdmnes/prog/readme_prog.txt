Here are the different fdmnes routines:
   main.f
   clemf0.f
   coabs.f
   convolution.f
   dirac.f
   fdm.f
   fprime.f
   general.f
   leture.f
   mat.f
   metric.f
   minim.f
   not_mpi.f
   optic.f
   potential.f
   scf.f
   selec.f
   spgroup.f
   sphere.f
   sub_util.f
   tab_data.f
   tddft.f
   tensor.f

They must be compiled and linked together.

mpif.h is a file which is automatically included during the compilation.

sub_util.f contains different lapack and blas routines.
It can be replaced by a call to these libraries.

When using the MPI compiler for parallelisation, the routine not_mpi.f and the file mpif.h
must be replaced by the call to the corresponding library

Makefile_example is an example of very simple makefile for the Portland-linux compiler.

