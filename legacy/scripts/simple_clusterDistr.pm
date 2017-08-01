#!/usr/bin/perl
package simple_clusterDistr;
use lib './';
use lib '../';
use strict;
use warnings;
our (@ISA, @EXPORT);
use Exporter;
@ISA = ("Exporter");
use simple_user_input;
use simple_clusterSpecs;
@EXPORT = qw(generate_distr_script); # generate_shmem_distr_script

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
1;
