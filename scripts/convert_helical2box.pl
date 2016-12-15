#!/usr/bin/perl
use warnings;
use strict;
if( scalar(@ARGV) < 2 ){
  die "Please input (1) filename of table of helical parameter files and (2) box size\n";
}
my$helitab = $ARGV[0];
my$box     = $ARGV[1];
open(HELI, "<$helitab") or die "Cannot open file: $helitab for reading: $!\n";
while(my$heliparam = <HELI>){
  chomp($heliparam);
  my$boxfile = $heliparam;
  $boxfile =~ s/\..+/\.box/;
  open(BOX, ">$boxfile") or die "Cannot open file: $boxfile for reading: $!\n";
  open(HELIP, "<$heliparam") or die "Cannot open file: $heliparam for reading: $!\n";
  while(my$line = <HELIP>){ 
    if( $line =~ /^#/ ){
    }else{
      my@nrs = split('\s+', $line);
      my$xcoord = rounder($nrs[0])-$box/2;
      my$ycoord = rounder($nrs[1])-$box/2;
      print BOX "$xcoord $ycoord $box $box -3\n";
    }
  }
  close(HELIP);
  close(BOX);
}
close(HELI);

sub rounder{
    my$number     = shift;
    my$number_int = int($number);
    my$diff       = $number-$number_int;
    if( $diff > 0.5 or $diff == 0.5 ){
        return $number_int+1;
    }else{
        return $number_int;
    }
}
