#!/usr/bin/perl
use strict;
use warnings;
my$ncls    = 0;
my$linecnt = 0;
my$frac    = 0.00;
while(<>){
  $linecnt = $linecnt+1;
  if( $_ =~ /frac\=(\d+\.\d+)\s/ ){
    $frac = $frac+$1;
  }
  if( $_ =~ /class\=(\d+)\.\d+\s/ ){
    if( $1 > $ncls ){
      $ncls = $1;
    }
  }
}
print "NR LINES: $linecnt\n";
$frac=$frac/$linecnt;
print "AVG FRAC: $frac\n";
print "NR CLUSTERS: $ncls\n";
my$nrefs=($frac/100.00)*$ncls;
print "AVG NREFS: $nrefs\n";
my$nrefs_tot = $linecnt*$nrefs;
print "TOT NR REFS: $nrefs_tot\n";
