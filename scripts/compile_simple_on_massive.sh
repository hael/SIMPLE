#!/bin/bash

MASSIVE=m3-login2.massive.org.au
MASSIVE_USERNAME=$1
MASSIVE_ACCOUNT=el85
cat >>EOF
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load  cuda/8.0.61 fftw/3.3.5-gcc5 cmake/3.5.2 git/2.8.1 gcc/5.4.0 openmpi 
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuild ] && rm -rf tmpbuild
mkdir tmpbuild
cd tmpbuild
cmake  -DUSE_OPENACC=ON  -DUSE_CUDA=ON -DCMAKE_BUILD_TYPE=DEBUG .. > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
source add2.bashrc
OMP_NUM_THREADS=4 simple_test_omp
simple_test_openacc
mpirun -np 10 simple_test_mpi
simple_test_cuda
EOF > buildsimple
scp buildsimple ${MASSIVE_USERNAME}@${MASSIVE}:~
scp ${MASSIVE_USERNAME}@${MASSIVE} sbatch -A ${MASSIVE_ACCOUNT} ~/buildsimple

cat >>EOF
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load  cuda/8.0.61 fftw/3.3.5-gcc5 cmake/3.5.2 git/2.8.1 openmpi/1.10.7-intel
export INTEL_DIR=/usr/local/intel/2015.0.090/
 . /usr/local/intel/2015.0.090/bin/compilervars.sh intel64 lp64
 . /usr/local/intel/2015.0.090/mkl/bin/mklvars.sh intel64 lp64
 . /usr/local/intel/2015.0.090/mpirt/bin/intel64/mpivars.sh intel64 lp64
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuild-intel ] && rm -rf tmpbuild-intel
mkdir tmpbuild-intel
cd tmpbuild-intel
 FC=ifort CC=icc CXX=icpc LDFLAGS= cmake -DVERBOSE=ON -DSIMPLE_BUILD_GUI=OFF -DUSE_OPENACC=ON -DUSE_MPI=ON -DCMAKE_BUILD_TYPE=DEBUG -Wdev  --debug-trycompile .. > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
source add2.bashrc
OMP_NUM_THREADS=4 simple_test_omp
simple_test_openacc
mpirun -np 10 simple_test_mpi
simple_test_cuda
EOF > buildsimpleintel



cat >>EOF
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=${MASSIVE_ACCOUNT}
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load  cuda/8.0.61 fftw/3.3.5-gcc5 cmake/3.5.2 git/2.8.1 gcc/5.4.0 openmpi 
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuild ] && rm -rf tmpbuild
mkdir tmpbuild
cd tmpbuild
FC=mpifort CC=mpicc CXX=mpicxx cmake  -DVERBOSE=ON -DSIMPLE_BUILD_GUI=OFF -DUSE_OPENACC=ON -DUSE_MPI=ON -DUSE_CUDA=ON -DCMAKE_BUILD_TYPE=DEBUG .. > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
source add2.bashrc
OMP_NUM_THREADS=4 simple_test_omp
simple_test_openacc
mpirun -np 10 simple_test_mpi
simple_test_cuda
EOF > buildsimple
scp buildsimple ${MASSIVE_USERNAME}@${MASSIVE}:~
scp ${MASSIVE_USERNAME}@${MASSIVE} sbatch -A ${MASSIVE_ACCOUNT} ~/buildsimple
