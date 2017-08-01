To compile this test code:

$: mpif90 -I /usr/lib64/mpi/gcc/openmpi/include pi_test.mpi.f90 -o mcpi

to run on n processor: nproc

$: mpirun -np nproc mcpi
