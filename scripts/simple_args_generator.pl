#!/usr/bin/perl

use warnings;
use strict;
use File::Compare;

# read the simple_params.f90 into an array
my @lines;
my @vars;
my$varlistfile;my$tmp_varlist;
my$simple_argsfile;
my$tmp_argsfile;
if ( -d $ENV{'SIMPLE_PATH'}.'/lib64'){
    $varlistfile=$ENV{'SIMPLE_PATH'}.'/lib64/simple/simple_varlist.txt';
    $simple_argsfile=$ENV{'SIMPLE_PATH'}.'/lib64/simple/simple_args.f90';
}else {
    $varlistfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_varlist.txt';
    $simple_argsfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_args.f90';
}
$tmp_varlist=$ENV{'SIMPLE_PATH'}.'/simple_varlist.tmp';
$tmp_argsfile=$ENV{'SIMPLE_PATH'}.'/simple_args.tmp';
open(PARAMS, "< simple_params.f90") or die "Cannot open simple_params.f90\n";
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

# generate the code
unlink($tmp_argsfile);
open(MODULE, "> ".$tmp_argsfile) or die "Cannot open simple_args.f90\n";
print MODULE "! for error checking of the SIMPLE command line arguments
module simple_args
use simple_defs
implicit none

public :: args, test_args
private

integer, parameter :: NARGMAX=500

type args
    private
    character(len=STDLEN) :: args(NARGMAX)
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
print MODULE  "    self%args(".$j.") = ''\n";
print MODULE  "end function

function is_present( self, arg ) result( yep )
    class(args), intent(in)      :: self
    character(len=*), intent(in) :: arg
    integer :: i
    logical :: yep
    yep = .false.
    do i=1,NARGMAX
        if( self%args(i) .eq. arg )then
            yep = .true.
            return
        endif
    end do
end function

subroutine test_args()
    use simple_filehandling, only: sys_stat, get_fileunit, nlines
    use simple_syscalls, only: exec_cmdline
    type(args) :: as
    character(len=STDLEN) :: vfilename,arg, errarg1, errarg2, errarg3, spath, srcpath
    integer :: funit, n, i
    integer,allocatable :: buff(:)
    integer :: status,length1,length2
    character(len=STDLEN) :: varlist = 'simple/simple_varlist.txt'
    write(*,'(a)') '**info(simple_args_unit_test): testing it all'
    write(*,'(a)') '**info(simple_args_unit_test, part 1): testing for args that should be present'
    as = args()
    funit = get_fileunit()
    write(*,'(a)') '**info(simple_args_unit_test): getting SIMPLE_PATH env variable'
    call getenv(\"SIMPLE_PATH\",spath)
    spath=adjustl(trim(spath))
    print *, 'get_environment_variable found SIMPLE_PATH    ', trim(spath)
    print *, 'simple_defs compile-time variable SIMPLE_PATH:', trim(SIMPLE_PATH)
    write(*,'(a)') '**info(simple_args_unit_test): getting SIMPLE_SOURCE_PATH variable'
    print *, 'simple_defs variable SIMPLE_SOURCE_PATH', trim(SIMPLE_SOURCE_PATH)
    srcpath = SIMPLE_SOURCE_PATH
    srcpath=adjustl(trim(srcpath))

    print *, 'appending varlist '
    vfilename = trim(adjustl(spath)) // '/lib/' // trim(adjustl(varlist))
    vfilename = trim(adjustl(vfilename))
    print *, 'varlist: ', trim(adjustl(vfilename))
    write(*,'(a,a)') '**info(simple_args_unit_test): checking varlist file ',trim(adjustl(vfilename))
    call sys_stat(vfilename,status,buff)
    if(status /= 0)then
      print *,' varlist not in lib/simple/,  checking lib64/simple'
      vfilename = trim(adjustl(spath)) // '/lib64/' // trim(adjustl(varlist))
      vfilename = trim(adjustl(vfilename))
      print *, 'varlist: ', trim(adjustl(vfilename))
      write(*,'(a)') '**info(simple_args_unit_test): checking varlist file'
      call sys_stat(vfilename,status,buff)
      if(status /= 0)then
        print *,' varlist not in lib64/simple/,  calling simple_args_varlist.pl'
        call exec_cmdline(\"cd \"//srcpath//\";simple_args_varlist.pl\")
        call sys_stat(vfilename,status,buff)
        if(status /= 0)then
          print *,' varlist still not in lib/simple/ after calling simple_args_varlist.pl'
        end if
      end if
    end if
    n = nlines(vfilename)
    open(unit=funit, status='old', action='read', file=vfilename)
    do i=1,n
        read(funit,*) arg
        if( as%is_present(arg) )then
            ! alles gut
        else
            write(*,'(a)') 'this argument should be present: ', arg
            stop 'part 1 of the unit test failed'
        endif
    end do
    close(funit)
    errarg1 = 'XXXXXXX'
    errarg2 = 'YYYYYY'
    errarg3 = 'ZZZZZZ'
    write(*,'(a)') '**info(simple_args_unit_test, part 2): testing for args that should NOT be present'
    if( as%is_present(errarg1) .or. as%is_present(errarg2) .or. as%is_present(errarg3) )then
        write(*,'(a)') 'the tested argumnets should NOT be present'
        stop 'part 2 of the unit test failed'
    endif
    write(*,'(a)') 'SIMPLE_ARGS_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
end subroutine

end module simple_args\n";
close(MODULE);

if ( -e ${simple_argsfile} && (compare(${simple_argsfile}, ${tmp_argsfile}) == 0)) {
      unlink(${tmp_argsfile}) # do nothing
     }
 else {
    rename  ${tmp_argsfile},${simple_argsfile}
 }
