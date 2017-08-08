#!/usr/bin/perl

use warnings;
use strict;
use Tie::File;

if( scalar(@ARGV) < 2 ){
    die "Need (1) filetable of Relion star files (2) box size in pixels\n";
}
my @starfiles;
tie @starfiles, 'Tie::File', $ARGV[0] or die "Cannot tie to file: $ARGV[0]\n";
my$box   = $ARGV[1];
my$boxo2 = int($box/2);
foreach my$file (@starfiles){
    if( -e $file ){
        my@starfile;
        tie @starfile, 'Tie::File', $file or die "Cannot tie to file: $file\n";
        my$boxfile = $file;
        $boxfile =~ s/star/box/;
        open(BOXHANDLE, ">$boxfile") or die "Cannot open file: $boxfile\n";
        foreach my$line (@starfile){
            if($line =~ /\d+\.\d+\s+\d+\.\d+/){
                my @entries = split(/\s+/,$line);
                printf BOXHANDLE "%u\t%u\t%u\t%u\t%d\n", $entries[1]-$boxo2, $entries[2]-$boxo2, $box, $box, -3
            }
        }
        close(BOXHANDLE);
        untie @starfile;
    }else{
        print "WARNING! This file does not exist: ", $file; 
    }
}
untie @starfiles;