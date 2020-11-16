#!/usr/bin/perl

use warnings;
use strict;
use Tie::File;

if( scalar(@ARGV) < 1 ){
    die "Need filetable of Relion *ctffind3.log files as input\n";
}
my @logfiles;
tie @logfiles, 'Tie::File', $ARGV[0] or die "Cannot tie to file: $ARGV[0]\n";
open(CTFHANDLE, ">simple_deftab.txt") or die "Cannot open file: simple_deftab.txt\n";
foreach my$file (@logfiles){
    my@logfile;
    tie @logfile, 'Tie::File', $file or die "Cannot tie to file: $file\n";
    foreach my$line (@logfile){
        if( $line =~ /Final Values/ ){
            my@vals = split(/\s+/, $line);
            my$dfx    = $vals[1]/10000.;
            my$dfy    = $vals[2]/10000.;
            my$angast = $vals[3];
            print CTFHANDLE "dfx=$dfx dfy=$dfy angast=$angast\n";
        }
    }
    untie @logfile;
}
close(CTFHANDLE);

