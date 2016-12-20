#!/usr/bin/perl
################################################################################
# Input variables that describes different cluster environments for input to   #
# simple_clusterDistr.pm                                                       #
#                                                                              #
#       Package simple_ClusterSpecs                                            #
#                                                                              #
# Lets the user input the global variables required for distributed simple     #
# execution.                                                                   #
################################################################################
#
package simple_clusterSpecs;
use strict;
use warnings;
our (@ISA, @EXPORT);
use Exporter;
@ISA = ("Exporter");
@EXPORT = qw($ALGNDOC_FBODY %LOCAL_DISTR_ENV %MASSIVE_DISTR_ENV %MASSIVE2_DISTR_ENV %MONARCH_DISTR_ENV %OXFORD_DISTR_ENV %OXFORD2_DISTR_ENV %OXFORD3_DISTR_ENV  $CVL_DISTR_ENV);

# LOCAL MODULE VARIABLES
my$NAME_DISTR       = 'distr_simple';
# my$NAME_SHMEM_DISTR = 'shmemdistr_simple';

# GLOBAL EXPORTED VARIABLES
our$ALGNDOC_FBODY = 'algndoc_';

#####################################################################
# DEFINES LOCAL ENVIROMENT DISTRIBUTED EXECUTION (FOR WORKSTATIONS) #
#####################################################################
our%LOCAL_DISTR_ENV;
$LOCAL_DISTR_ENV{'SUBMITCMD'}='';
$LOCAL_DISTR_ENV{'SCRIPT'}="#!/bin/bash
cd <<<EXECDIR>>> 
<<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

####################################################################
# DEFINES DISTRIBUTED EXECUTION ON THE MASSIVE 1 CLUSTER           #
####################################################################
our%MASSIVE_DISTR_ENV;
$MASSIVE_DISTR_ENV{'SUBMITCMD'}='qsub';
$MASSIVE_DISTR_ENV{'SCRIPT'}="#!/bin/bash
#PBS -N $NAME_DISTR
#PBS -l nodes=1:ppn=<<<NTHR>>>,mem=<<<MEMSTR>>>
#PBS -l walltime=<<<HOURS>>>:<<<MINUTES>>>:0
#PBS -o outfile.\$PBS_JOBID
#PBS -e errfile.\$PBS_JOBID
cd <<<EXECDIR>>> 
<<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

####################################################################
# DEFINES DISTRIBUTED EXECUTION ON THE MASSIVE 2 CLUSTER           #
####################################################################
our%MASSIVE2_DISTR_ENV;
$MASSIVE2_DISTR_ENV{'SUBMITCMD'}='sbatch';
$MASSIVE2_DISTR_ENV{'SCRIPT'}="#!/bin/bash
#SBATCH --mail-user=<<<<EMAIL>>>>
#SBATCH --mail-type=FAIL
#SBATCH --job-name=$NAME_DISTR
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=<<<NTHR>>>
#SBATCH --mem=<<<MEMSTR>>>
#SBATCH --time=0-<<<HOURS>>>:<<<MINUTES>>>:0
#SBATCH --output=outfile.%j
#SBATCH --error=errfile.%j
#SBATCH --partition=cryoem
#SBATCH --qos=vip_m2
cd <<<EXECDIR>>> 
<<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

####################################################################
# DEFINES DISTRIBUTED EXECUTION ON THE MONARCH CLUSTER             #
####################################################################
our%MONARCH_DISTR_ENV;
$MONARCH_DISTR_ENV{'SUBMITCMD'}='sbatch';
$MONARCH_DISTR_ENV{'SCRIPT'}="#!/bin/bash
#SBATCH --mail-user=<<<<EMAIL>>>>
#SBATCH --mail-type=FAIL
#SBATCH --job-name=$NAME_DISTR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<<<NTHR>>>
#SBATCH --mem=<<<MEMSTR>>>
#SBATCH --time=0-<<<HOURS>>>:<<<MINUTES>>>:0
#SBATCH --output=outfile.%j
#SBATCH --error=errfile.%j
#SBATCH --partition=cryo,batch
#SBATCH --account=p06
#SBATCH --qos=cryo_qos
cd <<<EXECDIR>>> 
<<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

####################################################################
# DEFINES DISTRIBUTED EXECUTION ON SUSANS CLUSTER IN OXFORD        #
####################################################################
our%OXFORD_DISTR_ENV;
$OXFORD_DISTR_ENV{'SUBMITCMD'}='qsub';
$OXFORD_DISTR_ENV{'SCRIPT'}="#!/bin/bash
#PBS -N $NAME_DISTR
#PBS -l nodes=1:ppn=<<<NTHR>>>,mem=<<<MEMSTR>>>
#PBS -l walltime=<<<HOURS>>>:<<<MINUTES>>>:0
#PBS -o outfile.\$PBS_JOBID
#PBS -e errfile.\$PBS_JOBID
#PBS -V
#PBS -l naccesspolicy=UNIQUEUSER
cd <<<EXECDIR>>> 
mpirun -np 1 --bind-to-socket --cpus-per-proc <<<NTHR>>> <<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

####################################################################
# DEFINES DISTRIBUTED EXECUTION ON THE OXFORD2 CLUSTER             #
####################################################################
our%OXFORD2_DISTR_ENV;
$OXFORD2_DISTR_ENV{'SUBMITCMD'}='sbatch';
$OXFORD2_DISTR_ENV{'SCRIPT'}="#!/bin/bash
#SBATCH --job-name=$NAME_DISTR
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=<<<NTHR>>>
#SBATCH --mem=<<<MEMSTR>>>
#SBATCH --time=0-<<<HOURS>>>:<<<MINUTES>>>:0
#SBATCH --output=outfile.%j
#SBATCH --error=errfile.%j
cd <<<EXECDIR>>> 
<<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

####################################################################
# DEFINES DISTRIBUTED EXECUTION ON THE OXFORD 3 CLUSTER             #
####################################################################
our%OXFORD3_DISTR_ENV;
$OXFORD3_DISTR_ENV{'SUBMITCMD'}='sbatch';
$OXFORD3_DISTR_ENV{'SCRIPT'}="#!/bin/bash
#SBATCH --job-name=$NAME_DISTR
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=<<<NTHR>>>
#SBATCH --mem=<<<MEMSTR>>>
#SBATCH --time=0-<<<HOURS>>>:<<<MINUTES>>>:0
#SBATCH --output=outfile.%j
#SBATCH --error=errfile.%j
#SBATCH --partition=Intel_Cluster

cd <<<EXECDIR>>> 
<<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

####################################################################
# DEFINES DISTRIBUTED EXECUTION ON CVL                             #
####################################################################
our%CVL_DISTR_ENV;
$CVL_DISTR_ENV{'SUBMITCMD'}='sbatch';
$CVL_DISTR_ENV{'SCRIPT'}="#!/bin/bash
#SBATCH --mail-user=<<<<EMAIL>>>>
#SBATCH --mail-type=FAIL
#SBATCH --job-name=$NAME_DISTR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<<<NTHR>>>
#SBATCH --mem=<<<MEMSTR>>>
#SBATCH --time=0-<<<HOURS>>>:<<<MINUTES>>>:0
#SBATCH --output=outfile.%j
#SBATCH --error=errfile.%j
#SBATCH --partition=compute
#SBATCH --account=cvl
cd <<<EXECDIR>>> 
<<<CMDSTRING>>> fromp=<<<START>>> top=<<<STOP>>> part=<<<PART>>> nthr=<<<NTHR>>> outfile=$ALGNDOC_FBODY<<<PART>>>.txt > OUT<<<PART>>>\nexit\n";

1;