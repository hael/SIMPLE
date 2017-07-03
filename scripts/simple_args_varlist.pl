#!/usr/bin/perl
use warnings;
use strict;
# read the simple_params.f90 into an array
my @lines;
my @vars;

my$src_path;
my$simple_path;
$src_path = $ENV{'SIMPLE_SOURCE_PATH'};
$simple_path = $ENV{'SIMPLE_PATH'};
if ( -e  $src_path and -d $src_path ){
    die 'simple_args_varlist cannot find SIMPLE_PATH or is not set';
}
if ( -e  $src_path and -d $src_path ){
    die 'simple_args_varlist cannot find SIMPLE_SOURCE_PATH';
}

my$varlistfile;
my$simple_argsfile;
if ( -e  $ENV{'SIMPLE_PATH'}.'/lib64' and -d $ENV{'SIMPLE_PATH'}.'/lib64' ){
$varlistfile=$ENV{'SIMPLE_PATH'}.'/lib64/simple/simple_varlist.txt';
$simple_argsfile=$ENV{'SIMPLE_PATH'}.'/lib64/simple/simple_args.f90';
}else {
$varlistfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_varlist.txt';
$simple_argsfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_args.f90';
}
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
unlink($simple_path."/tests/simple_varlist.txt");
open(VLIST, "> ".$simple_path."/tests/simple_varlist.txt") or die "Cannot open simple_varlist.txt\n";
foreach (@vars){
  print VLIST $_, "\n";
}
close(VLIST);
