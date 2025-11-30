#!/usr/bin/perl

use warnings;
use strict;
use File::Compare;

# read the simple_parameters.f90 file into an array
my @lines;
my @vars;
my $varlistfile;     my $tmp_varlist;
my $simple_argsfile; my $tmp_argsfile;
if ( -d $ENV{'SIMPLE_PATH'}.'/lib64'){
    $varlistfile=$ENV{'SIMPLE_PATH'}.'/lib64/simple/simple_varlist.txt';
    $simple_argsfile=$ENV{'SIMPLE_PATH'}.'/lib64/simple/simple_args.f90';
}else {
    $varlistfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_varlist.txt';
    $simple_argsfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_args.f90';
}
$tmp_varlist  = $ENV{'SIMPLE_PATH'}.'/simple_varlist.tmp';
$tmp_argsfile = $ENV{'SIMPLE_PATH'}.'/simple_args.tmp';
open(PARAMS, "< simple_parameters.f90") or die "Cannot open simple_parameters.f90\n";
@lines = <PARAMS>;
close(PARAMS);

# extract the relevant lines
LINESLOOP:{
    foreach (@lines) {
        chomp($_);
        if( $_ =~ /^\s+contains/ ){
            last LINESLOOP;
        }
        if( $_ =~ /^\s+character/ or ($_ =~ /^\s+integer/ or ($_ =~ /^\s+real/ or ($_ =~ /^\s+logical/ or $_ =~ /^\s+type\(string\)/)))){
            $_ =~ s/.+:://;
            my @splitted = split(",",$_);
            push(@vars,@splitted);
        }
    }
}

# remove unwanted stuff
foreach my$i (0..$#vars){
    $vars[$i] =~ s/\=.+//;
    $vars[$i] =~ s/\!.+//;
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

# file-handling
unlink(${tmp_varlist});
open(VLIST, "> ".${tmp_varlist}) or die "Cannot open simple_varlist.txt\n";
foreach (@vars){
    print VLIST $_, "\n";
}
close(VLIST);
if ( -e ${varlistfile} && (compare( ${varlistfile}, ${tmp_varlist})==0) ){
    unlink(${tmp_varlist})
} else {
    rename ${tmp_varlist}, ${varlistfile}
}

# generate the code
unlink($tmp_argsfile);
open(MODULE, "> ".$tmp_argsfile) or die "Cannot open simple_args.f90\n";
print MODULE "! for error checking of the SIMPLE command line arguments

module simple_args
include 'simple_lib.f08'
implicit none

public :: args
private
#include \"simple_local_flags.inc\"

integer, parameter :: NARGMAX=1000, MINVARS=100

type args
    private
    type(string) :: args(NARGMAX)
  contains
    procedure :: is_present
end type

interface args
    module procedure constructor
end interface

contains

function constructor( ) result( self )
    type(args) :: self\n";
my$j;
foreach my$i (0 .. $#vars){
  $j = $i+1;
  my$str = "    self%args(".$j.") = '".$vars[$i]."'\n";
  print MODULE $str;
}
$j++;
print MODULE  "    self%args(".$j.") = 'job_memory_per_task'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'projname'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'simple_path'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'time_per_image'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'user_account'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'user_email'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'user_project'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'qsys_partition'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'qsys_qos'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'qsys_reservation'\n";
$j++;
print MODULE  "    self%args(".$j.") = 'job_name'\n";
$j++;
print MODULE  "    self%args(".$j.") = ''\n";
print MODULE  "end function

function is_present( self, arg ) result( yep )
    class(args),      intent(in) :: self
    character(len=*), intent(in) :: arg
    integer :: i
    logical :: yep
    yep = .false.
    do i=1,NARGMAX
        if( self%args(i)%to_char() .eq. arg )then
            yep = .true.
            return
        endif
    end do
end function

end module simple_args\n";
close(MODULE);

# file-handling
if ( -e ${simple_argsfile} && (compare(${simple_argsfile}, ${tmp_argsfile}) == 0)) {
      unlink(${tmp_argsfile}) # do nothing
     }
 else {
    rename  ${tmp_argsfile},${simple_argsfile}
 }
