#!/usr/bin/perl
use strict;
use warnings;
my $file = $ARGV[0] or die "Need detect_atoms output file on the command line\n";
open(my $data, '<', $file) or die "Could not open '$file' $!\n";
while (my $line = <$data>) {
    chomp $line; 
    if( $line =~ /# atoms, final\s+(\d+)/ ){
        print "NATOMS: ", $1, "\n";
    }
    if( $line =~ /VALID_CORR Average:\s+(\d\.\d+)/ ){
        print "VALID_CORR, avg: ", $1, "\n";
    }
    if( $line =~ /VALID_CORR Sigma  :\s+(\d\.\d+)/ ){
        print "VALID_CORR, sig: ", $1, "\n";
    }
}


