#!/usr/bin/perl
use warnings;
use strict;
use Tie::File;
if( scalar(@ARGV) < 2 ){
    die "Need (1) filetable of SAM dat files (2) binning fatcor (3) box size in pixels\n";
}
my @samfiles;
tie @samfiles, 'Tie::File', $ARGV[0] or die "Cannot tie to file: $ARGV[0]\n";
my$bin   = $ARGV[1];
my$box   = $ARGV[2];
my$boxo2 = int($box/2);
foreach my$file (@samfiles){
    if( -e $file ){
        my@samfile;
        tie @samfile, 'Tie::File', $file or die "Cannot tie to file: $file\n";
        my$boxfile = $file;
        $boxfile =~ s/dat/box/;
        open(BOXHANDLE, ">$boxfile") or die "Cannot open file: $boxfile\n";
        foreach my$line (@samfile){
            my @entries = split(/\s+/,$line);
            printf BOXHANDLE "%u\t%u\t%u\t%u\t%d\n", $entries[3]+$boxo2, $entries[4]-$boxo2, $box, $box, -3

        }
        close(BOXHANDLE);
        untie @samfile;
    }else{
        print "WARNING! This file does not exist: ", $file; 
    }
}
untie @samfiles;