#!/usr/bin/perl

use warnings;
use strict;

while(my $line=<STDIN>){
    chomp($line);
    if( $line !~ /_r/ ){
        my@row  = split(/\s+/,$line);
        if( defined($row[5]) ){
            my $kv     = $row[5];
            my $dfx    = $row[6]/10000.;
            my $dfy    = $row[7]/10000.;
            my $angast = $row[8];
            my $cs     = $row[9];
            my $fraca  = $row[13];
            print "kv=$kv dfx=$dfx dfy=$dfy angast=$angast cs=$cs fraca=$fraca\n";
        }
    }
}
