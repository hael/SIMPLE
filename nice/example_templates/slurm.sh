#!/usr/bin/bash
#SBATCH -c XXXNCPUXXX
#SBATCH --mem=10000

export SIMPLE_DEFAULT_PARTITION_TIME=true
export SIMPLE_PATH=/home/jcaesar/Code/SIMPLE/build
export SIMPLE_EMAIL=''
export SIMPLE_QSYS=slurm

echo $$ > nice.pid
echo "Running on $HOSTNAME"
XXXSIMPLEXXX
