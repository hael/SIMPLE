#!/usr/bin/perl
use warnings;
use strict;
#use File::Basename;
#use File::Grep qw( fgrep );
use Cwd qw(cwd);
my $dir = cwd;
# read the simple_params.f90 into an array
my @lines;
my @vars;

my$varlistfile;
my$tmp_varlist;
my$Gitversionfile;
my$simple_argsfile;
my$source_dir="";
my$sline="";
print " In simple_args_varlist.pl Current dir: ", $dir, "\n";

if ( -d $ENV{'SIMPLE_PATH'}.'/lib64'){
    $varlistfile=$ENV{'SIMPLE_PATH'}.'/lib64/simple/simple_varlist.txt';
}else {
    $varlistfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_varlist.txt';
}
$tmp_varlist=$ENV{'SIMPLE_PATH'}.'/simple_varlist.tmp';

$Gitversionfile=$ENV{'SIMPLE_PATH'}."/lib/simple/SimpleGitVersion.h";
print " Reading $Gitversionfile";
# $source_dir = fgrep{ /SIMPLE_SOURCE_PATH/ } $Gitversionfile;
open my $fh, '<', $Gitversionfile or die "Could not open file  $Gitversionfile:
+ $!";

while ( my$line = <$fh> ) {
    print $line;
    if ( $line =~ m/SIMPLE_SOURCE_PATH/ ){
        $sline = $line =~ s/.* = \"([^"]*)\".*/$1/gr;
        last;
    }
}
print "SIMPLE SOURCE DIR LINE:", $sline, "\n";
$source_dir = $sline;
close $fh  or die "Could not close file $Gitversionfile: $!";
#$source_dir =`sed 's/SIMPLE_SOURCE_PATH=\"\([^"]*\)\"/\1/'  $Gitversionfile`;

print "SIMPLE SOURCE DIR:", $source_dir, "\n";

open(PARAMS, "<".$source_dir."/src/main/simple_params.f90") or die "Cannot open simple_params.f90\n simple_args_varlist.pl must be called from <simple source directory>/src/main ";
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


unlink(${tmp_varlist});
open(VLIST, "> ".${tmp_varlist}) or die "Cannot open simple_varlist.txt\n";
foreach (@vars){
  print VLIST $_, "\n";
}
close(VLIST);
if ( -e ${varlistfile} && (compare( ${varlistfile}, ${tmp_varlist})==0) ){
    unlink(${tmp_varlist})# do nothing
} else {
    rename  ${tmp_varlist},${varlistfile}
}
