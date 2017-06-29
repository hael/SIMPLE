#!/usr/bin/perl
use warnings;
use strict;
# read the simple_params.f90 into an array
my @lines;
my @vars;
my$src_path;
my$build_path;
$src_path = $ENV{'SIMPLE_SOURCE_PATH'};
$build_path = $ENV{'SIMPLE_PATH'};
open(PARAMS, "< ".$src_path."/src/simple_main/simple_params.f90") or die "Cannot open simple_params.f90\n";
@lines = <PARAMS>;
close(PARAMS);
# extract the relevant lines
LINESLOOP:{
  foreach (@lines) {
    chomp($_);
    if( $_ =~ /^\s+contains/ ){
      last LINESLOOP;
    }
    if( $_ =~ /^\s+character/ or ($_ =~ /^\s+integer/ or ($_ =~ /^\s+real/ or $_ =~ /^\s+logical/))){
      $_ =~ s/.+:://;
      my @splitted = split(",",$_);
      push(@vars,@splitted);
    }
  }
}
# remove unwanted stuff
foreach my$i (0..$#vars){
  $vars[$i] =~ s/\=.+//;
  if( $vars[$i] =~ /\(/ or ( $vars[$i] =~ /\)/ or $vars[$i] =~ /^\s+$/ ) ){
    delete $vars[$i];
  }
  if( defined($vars[$i]) ){ $vars[$i] =~ s/\s//g };
}
@vars = grep defined, @vars;
# add wanted stuff (volumes)
foreach my$i (1 .. 20){
  my$str = 'vol'.$i;
  push(@vars,$str);
}
unlink($build_path."/tests/simple_varlist.txt");
open(VLIST, "> ".$build_path."/tests/simple_varlist.txt") or die "Cannot open simple_varlist.txt\n";
foreach (@vars){
  print VLIST $_, "\n";
}
close(VLIST);
