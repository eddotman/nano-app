# Make file for the fdmnes package


pgf90 -c -fast -O3 -Mvect=sse main.f clemf0.f coabs.f convolution.f dirac.f fdm.f fprime.f general.f 
pgf90 -c -fast -O3 -Mvect=sse lecture.f mat.f metric.f minim.f not_mpi.f  potential.f optic.f
pgf90 -c -fast -O3 -Mvect=sse scf.f selec.f sphere.f spgroup.f sub_util.f tab_data.f tddft.f tensor.f 
pgf90 -o fdmnes main.o clemf0.o coabs.o convolution.o dirac.o fdm.o fprime.o general.o lecture.o mat.o metric.o minim.o not_mpi.o potential.o scf.o selec.o spgroup.o sphere.o sub_util.o tab_data.o tddft.o tensor.o optic.o
