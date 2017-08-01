#!/usr/bin/perl

use strict;
use warnings;

(@ARGV != 2)&&(
die "Give to script:
1) number of objects
2) number of partitions\n");

my$nobjs = $ARGV[0];
my$npart = $ARGV[1];

my$nobjs_per_part = int(int($nobjs)/int($npart));
my$leftover = int($nobjs)-$nobjs_per_part*int($npart);
print "nobjs_per_part: $nobjs_per_part\n";
print "leftover: $leftover\n";
my $stop  = 0;
my $start = 0;
for(my $i=1; $i<=$npart; $i++){
    if( $i == $npart ){
        $start = $stop+1;
        $stop  = $nobjs;
    }else{
        if( $leftover == 0 ){
            $start = $stop+1;
            $stop  = $start+$nobjs_per_part-1;
        }else{
            $stop  = $i*($nobjs_per_part+1);
            $start = $stop-($nobjs_per_part+1)+1;
            $leftover = $leftover-1;
        }
    }
    print "Part: $i Start: $start Stop: $stop\n";
}