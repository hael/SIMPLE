#!/usr/bin/perl
package simple_clusterDistr;
use lib '/Users/hansolo/src/fortran/lib/simple_wfrederic/Simple_Restruct.projet/scripts';
use lib '/Users/hansolo/src/fortran/lib/simple_wfrederic/Simple_Restruct.projet';
use strict;
use warnings;
our (@ISA, @EXPORT);
use Exporter;
@ISA = ("Exporter");
use simple_user_input;
use simple_clusterSpecs;
@EXPORT = qw($SIMPLESYS $SUBMITCMD $TIME_PER_IMAGE $SIMPLEBIN generate_distr_script generate_shmem_distr_script);

#####################################################################
# USER-DEFINED VARIABLES THAT CONTROL WHICH CLUSTER ENVIRONMENT TO  #
# USE AND HOW TO USE THE RESOURCES                                  #
#####################################################################
our$SIMPLESYS      = 'LOCAL';                   # Name of system
my%DISTR_ENV       = %LOCAL_DISTR_ENV;          # Defines the environment for distributed execution     
my$EMAIL           = 'hans.elmlund@monash.edu'; # e-mail for failure report
my$NTHR            = 4;                         # number of threads (CPUs per core)
my$MEMSTR          = '500';                     # string descriptor for memory
our$TIME_PER_IMAGE = 100;                       # time per image (in seconds)

# SETTINGS FOR COMLIN_CLUSTER ON MASSIVE2
# $NTHR           = 1
# $MEMSTR         = '500';
# $TIME_PER_IMAGE = 5;

# SETTINGS FOR PRIME2/PRIME2_CLUSTER ON MASSIVE2
# $NTHR           = 8
# $MEMSTR         = '32000';
# $TIME_PER_IMAGE = =100;

# SETTINGS FOR PRIME2/PRIME2_CLUSTER ON OXFORD
# $NTHR           = 8
# $MEMSTR         = '9gb';
# $TIME_PER_IMAGE = =100;

#####################################################################
# OTHER VARIABLES TO EXPORT                                         #
#####################################################################
our$SIMPLEBIN = $SIMPLE_PATH.'/bin'; # location of compiled binaries
our$SUBMITCMD = $DISTR_ENV{'SUBMITCMD'};   

#####################################################################
# SUBROUTINES TO EXPORT                                             #
#####################################################################
sub generate_distr_script{
    my$hours     = shift;
    my$minutes   = shift;
    my$execdir   = shift;
    my$cmdstring = shift;
    my$start     = shift;
    my$stop      = shift;
    my$part      = shift;
    my $script   = $DISTR_ENV{'SCRIPT'};
    $script =~ s/\<\<\<EMAIL\>\>\>/$EMAIL/g;
    $script =~ s/\<\<\<NTHR\>\>\>/$NTHR/g;
    $script =~ s/\<\<\<MEMSTR\>\>\>/$MEMSTR/g;
    $script =~ s/\<\<\<HOURS\>\>\>/$hours/g;
    $script =~ s/\<\<\<MINUTES\>\>\>/$minutes/g;
    $script =~ s/\<\<\<EXECDIR\>\>\>/$execdir/g;
    $script =~ s/\<\<\<CMDSTRING\>\>\>/$cmdstring/g;
    $script =~ s/\<\<\<START\>\>\>/$start/g;
    $script =~ s/\<\<\<STOP\>\>\>/$stop/g;
    $script =~ s/\<\<\<PART\>\>\>/$part/g;
    if( $NTHR == 1 ){
        my $str = "#SBATCH --ntasks-per-socket=1\n";
        $script =~ s/$str//;
    }
    return $script;    
}

sub generate_shmem_distr_script{
    my$hours     = shift;
    my$minutes   = shift;
    my$execdir   = shift;
    my$cmdstring = shift;
    my $script   = $DISTR_ENV{'SHMEMSCRIPT'};
    $script =~ s/\<\<\<EMAIL\>\>\>/$EMAIL/g;
    $script =~ s/\<\<\<NTHR\>\>\>/$NTHR/g;
    $script =~ s/\<\<\<MEMSTR\>\>\>/$MEMSTR/g;
    $script =~ s/\<\<\<HOURS\>\>\>/$hours/g;
    $script =~ s/\<\<\<MINUTES\>\>\>/$minutes/g;
    $script =~ s/\<\<\<EXECDIR\>\>\>/$execdir/g;
    $script =~ s/\<\<\<CMDSTRING\>\>\>/$cmdstring/g;
    return $script;    
}

1;