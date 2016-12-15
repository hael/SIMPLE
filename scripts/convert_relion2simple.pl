#!/usr/bin/perl

use warnings;
use strict;

while(my $file=<STDIN>){
    chomp($file);
    open(FOO, "<$file") or die "Cannot open: $_\n";
    while(my $line=<FOO>){
        chomp($line);
        if( $line !~ /_r/ ){
            my@row  = split(/\s+/,$line);
            if( defined($row[9]) ){
                my $dfx = $row[9]/10000.;
                my $dfy = $row[10]/10000.;
                print "dfx=$dfx dfy=$dfy angast=$row[11]\n";
            }        
        }
    }
    close(FOO);
}
