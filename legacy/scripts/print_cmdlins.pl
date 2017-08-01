#!/usr/bin/perl
@bins = glob("simple*");
foreach $bin (@bins){
  system($bin);
}
