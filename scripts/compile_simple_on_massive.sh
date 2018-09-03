#!/bin/bash

MASSIVE=m3.massive.org.au
MASSIVE_USERNAME=$1
MASSIVE_ACCOUNT=el85

## STANDARD BUILD: GCC 5, fftw, static, OpenMP

cat <<EOF  > buildsimple
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load fftw/3.3.5-gcc5 cmake/3.5.2 git/2.8.1 gcc/5.4.0
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuild ] && rm -rf tmpbuild
mkdir tmpbuild
cd tmpbuild
cmake  -DUSE_OPENACC=OFF -DUSE_MPI=OFF  -DUSE_CUDA=OFF -DCMAKE_BUILD_TYPE=DEBUG .. -Wdev  --debug-trycompile > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
. add2.bashrc
ctest -V
OMP_NUM_THREADS=4 simple_test_omp
EOF
# scp buildsimple ${MASSIVE_USERNAME}@${MASSIVE}:~
# scp ${MASSIVE_USERNAME}@${MASSIVE} sbatch -A ${MASSIVE_ACCOUNT} ~/buildsimple


## LATEST GCC BUILD: GCC 8, fftw, static, OpenMP

cat <<EOF > buildsimplegcc8
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load  git fftw cmake/3.5.2  gcc/8.1.0
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuild-gcc8 ] && rm -rf tmpbuild-gcc8
mkdir tmpbuild-gcc8
cd tmpbuild-gcc8
cmake  -DUSE_OPENACC=OFF -DUSE_MPI=OFF  -DUSE_CUDA=OFF -DCMAKE_BUILD_TYPE=DEBUG ..  -Wdev  --debug-trycompile > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
. add2.bashrc
ctest -V  > log_check 2> log_check_err
OMP_NUM_THREADS=4 simple_test_omp
EOF
# scp buildsimplegcc8 ${MASSIVE_USERNAME}@${MASSIVE}:~
# scp ${MASSIVE_USERNAME}@${MASSIVE} sbatch -A ${MASSIVE_ACCOUNT} ~/buildsimplegcc8

## INTEL BUILD: IFORT 17, MKL, static, OpenMP

cat <<EOF  > buildsimpleintel
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load git cmake/3.5.2 intel/2017u4
export INTEL_DIR=/usr/local/intel/2017u4/
 . ${INTEL_DIR}/bin/compilervars.sh intel64
 . ${INTEL_DIR}/mkl/bin/mklvars.sh intel64 lp64
#  . /usr/local/intel/2015.0.090/mpirt/bin/intel64/mpivars.sh intel64 lp64
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuild-intel ] && rm -rf tmpbuild-intel
mkdir tmpbuild-intel
cd tmpbuild-intel
#FC=mpifort CC=mpicc CXX=mpicxx 
 FC=ifort CC=icc CXX=icpc LDFLAGS="" cmake -DVERBOSE=ON -DSIMPLE_BUILD_GUI=OFF -DUSE_OPENACC=OFF -DUSE_MPI=OFF -DCMAKE_BUILD_TYPE=RELEASE -Wdev -DINTEL_OMP_OVERRIDE=ON -DUSE_AUTO_PARALLELISE=ON --debug-trycompile .. > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
. add2.bashrc
ctest -V > log_check 2> log_check_err
OMP_NUM_THREADS=4 simple_test_omp
simple_test_openacc
mpirun -np 10 simple_test_mpi
simple_test_openacc
EOF

## CUDA/MPI BUILD: GCC 5, CUDA 8, OpenMPI 1.10.3, OpenMP, FFTW 3.3.5, shared

cat <<EOF > buildsimplecudampi
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load  cuda/8.0.61 fftw/3.3.5-gcc5 cmake/3.5.2 git/2.8.1 gcc/5.4.0 openmpi/1.10.3-gcc5-mlx
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuildall ] && rm -rf tmpbuildall
mkdir tmpbuildall
cd tmpbuild
FC=mpifort CC=mpicc CXX=mpicxx cmake  -DVERBOSE=ON -DSIMPLE_BUILD_GUI=OFF -DUSE_OPENACC=ON -DUSE_MPI=ON -DUSE_CUDA=ON -DCMAKE_BUILD_TYPE=DEBUG -BUILD_SHARED_LIBS=ON .. > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
. add2.bashrc
ctest -V > log_check 2> log_check_err
OMP_NUM_THREADS=4 simple_test_omp
simple_test_openacc
srun simple_test_mpi
module load virtualgl
vglrun simple_test_cuda
EOF
#scp buildsimplecudampi ${MASSIVE_USERNAME}@${MASSIVE}:~
# scp ${MASSIVE_USERNAME}@${MASSIVE} sbatch -A ${MASSIVE_ACCOUNT} ~/buildsimplecudampi

## CUDA/MPI BUILD: GCC 4.9, CUDA 8, OpenMPI 1.10.3, OpenMP, FFTW 3.3.4, shared

cat <<EOF > buildsimple-gcc4.9-mpi
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load  fftw/3.3.4-gcc cmake/3.5.2 git/2.8.1 gcc/4.9.3 openmpi/1.10.3-gcc4-mlx-verbs
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuildgcc4 ] && rm -rf tmpbuildgcc4
mkdir tmpbuildgcc4
cd tmpbuild
FC=mpifort CC=mpicc CXX=mpicxx cmake  -DVERBOSE=ON -DSIMPLE_BUILD_GUI=OFF -DUSE_OPENACC=OFF -DUSE_MPI=ON -DUSE_CUDA=OFF -DCMAKE_BUILD_TYPE=DEBUG -BUILD_SHARED_LIBS=ON .. > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
. add2.bashrc
ctest -V > log_check 2> log_check_err
OMP_NUM_THREADS=4 simple_test_omp
srun ./bin/simple_test_mpi
EOF
#scp buildsimple-gcc49-mpi ${MASSIVE_USERNAME}@${MASSIVE}:~
# scp ${MASSIVE_USERNAME}@${MASSIVE} sbatch -A ${MASSIVE_ACCOUNT} ~/buildsimple-gcc49-mpi
