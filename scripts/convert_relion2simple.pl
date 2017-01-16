#!/usr/bin/perl

use warnings;
use strict;

while(my $line=<STDIN>){
    chomp($line);
    if( $line !~ /_r/ ){
        my@row  = split(/\s+/,$line);
        if( defined($row[8]) ){
            my $kv     = $row[8];
            my $dfx    = $row[9]/10000.;
            my $dfy    = $row[10]/10000.;
            my $angast = $row[11];
            my $cs     = $row[12];
            my $fraca  = 0.1;
            print "kv=$kv dfx=$dfx dfy=$dfy angast=$angast cs=$cs fraca=$fraca\n";
        }
    }
}
