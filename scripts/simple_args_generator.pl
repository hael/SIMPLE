#!/usr/bin/perl
use warnings;
use strict;
# read the simple_params.f90 into an array
my @lines;
my @vars;
my$varlistfile;
my$simple_argsfile;
$varlistfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_varlist.txt';
$simple_argsfile=$ENV{'SIMPLE_PATH'}.'/lib/simple/simple_args.f90';
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

unlink(${varlistfile});
open(VLIST, "> ".${varlistfile}) or die "Cannot open simple_varlist.txt\n";
foreach (@vars){
  print VLIST $_, "\n";
}
close(VLIST);
# generate the code
unlink($simple_argsfile);
open(MODULE, "> ".$simple_argsfile) or die "Cannot open simple_args.f90\n";
print MODULE "!==Class simple_args
!
!\> \\brief simple_args is for error checking of the SIMPLE command line arguments.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution
! or modification is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2011-08-18.
!
!==Changes are documented below
!
module simple_args
use simple_defs
implicit none

public :: args, test_args
private
#include \"simple_local_flags.inc\"
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
  my$str = "        self%args(".$j.") = '".$vars[$i]."'\n";
  print MODULE $str;
}
$j++;
print MODULE  "        self%args(".$j.") = ''\n";
print MODULE  "    end function

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

    subroutine test_args(vlist)
        use simple_filehandling, only: get_fileunit, nlines
        character(len=STDLEN), intent(inout), optional :: vlist
        type(args) :: as
        character(len=STDLEN) :: arg, errarg1, errarg2, errarg3, spath, srcpath
        integer :: funit, n, i
        integer,dimension(13) :: buff
        integer :: status,length
        write(*,'(a)') '**info(simple_args_unit_test): testing it all'
        write(*,'(a)') '**info(simple_args_unit_test, part 1): testing for args that should be present'
        as = args()
        funit = get_fileunit()
        write(*,'(a)') '**info(simple_args_unit_test): getting SIMPLE_PATH env variable'
        call get_environment_variable(\"SIMPLE_PATH\",spath,LENGTH=length,STATUS=status) ! 'simple absolute path'
        if( status /= 0 ) ErrorPrint(\"get_environment_variable failed to find SIMPLE_PATH!!!\")
        write(*,'(a)') '**info(simple_args_unit_test): getting SIMPLE_SOURCE_PATH env variable'
        call get_environment_variable(\"SIMPLE_SOURCE_PATH\",srcpath,LENGTH=length,STATUS=status) ! 'simple absolute source path'
        if( status /= 0 ) ErrorPrint(\"get_environment_variable failed to find SIMPLE_SOURCE_PATH!!!\")
        spath=adjustl(trim(spath))
        srcpath=adjustl(trim(spath))
        if(present(vlist))then
            vlist=adjustl(trim(vlist))
        else
            vlist = adjustl(trim(spath)) // '/lib/simple/simple_varlist.txt'
        end if
        write(*,'(a)') '**info(simple_args_unit_test): checking varlist file'
        call stat(vlist , buff, status)
        if(status /= 0)then
            WarnPrint \"simple_args_unit_test: varlist not in <build dir>/lib/simple/,  calling simple_args_varlist\"
            call system(\"simple_args_varlist.pl\",status)
            if(status /= 0) ErrorPrint(\"simple_args_unit_test:  simple_args_varlist.pl failed\")
        end if
        n = nlines(vlist)
        open(unit=funit, status='old', action='read', file=vlist)
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
