#!/bin/bash

MASSIVE=m3-login2.massive.org.au
MASSIVE_USERNAME=$1

cat >>EOF
#!/bin/bash
#SBATCH --job-name=CompileSIMPLE
#SBATCH --account=el85
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=12
#SBATCH --qos=shortq

module load fftw/3.3.5-gcc5 cmake/3.5.2 git/2.8.1 gcc/5.4.0 openmpi/1.10.3-gcc5 cuda/8.0.61 relion gctf/1.06_cuda8
cd ~/${SLURM_JOB_ACCOUNT}_scratch/${MASSIVE_USERNAME}/SIMPLE3.0
git pull --rebase
[ -d tmpbuild ] && rm -rf tmpbuild
mkdir tmpbuild
cd tmp build
cmake  -DUSE_OPENACC=ON -DUSE_MPI=ON -DUSE_CUDA=ON -DCMAKE_BUILD_TYPE=RELEASE .. > log_cmake 2> log_cmake_err
make -j12 install > log_make 2> log_make_err
source add2.bashrc

OMP_NUM_THREADS=4 simple_test_omp
simple_test_openacc
mpirun -np 10 simple_test_mpi
simple_test_cuda
EOF > buildsimple
scp buildsimple ${MASSIVE_USERNAME}@${MASSIVE}:~
scp ${MASSIVE_USERNAME}@${MASSIVE} sbatch -A el85 ~/buildsimple
