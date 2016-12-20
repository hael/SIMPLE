#!/usr/bin/perl

use warnings;
use strict;

my @prgnames;
open( EXEC,    "<", "./production/simple/simple_exec/simple_exec.f90") or die "Cannot open: $!\n";
while(my$line=<EXEC>){
    if( $line =~ /case\(\s*\'(\w+)\'\s*\)/ ){
        push(@prgnames,$1);
    }
}
close( EXEC );
my @sorted_prgnames = sort { lc($a) cmp lc($b) } @prgnames;
foreach my$prgname (@sorted_prgnames){
    print "$prgname\n";
}

