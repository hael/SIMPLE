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

public :: args, test_args
private
#include \"simple_local_flags.inc\"

integer, parameter :: NARGMAX=500, MINVARS=100

type args
    private
    character(len=LONGSTRLEN) :: args(NARGMAX)
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
print MODULE  "    self%args(".$j.") = 'qsys_name'\n";
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
    type(args) :: as
    character(len=LONGSTRLEN) :: vfilename,arg, errarg1, errarg2, errarg3, spath, cwd, bpath,cmd
    integer               :: funit, n, i,io_stat, status
    integer, allocatable  :: buff(:)
    character(len=*),parameter :: varlist = 'simple_varlist.txt'
    write(*,'(a)') '**info(simple_args_unit_test): testing it all'
    write(*,'(a)') '**info(simple_args_unit_test, part 1): testing for args that should be present'
    as = args()
    write(*,'(a)') '**info(simple_args_unit_test): getting SIMPLE_PATH env variable'
    status = simple_getenv(\"SIMPLE_PATH\",spath)
    if (status /= 0) THROW_HARD( 'test_args SIMPLE_PATH env variable not found')
    spath = trim(adjustl(spath))
    print *, 'get_environment_variable found SIMPLE_PATH ', trim(adjustl(spath))
    write(*,'(a)') '**info(simple_args_unit_test): getting current working directory'
    call simple_getcwd(cwd)
    cwd = trim(adjustl(cwd))
    print *, 'Current working directory ', trim(adjustl(cwd))
    write(*,'(a)') '**info(simple_args_unit_test): getting directory that contains varlist'
    if( .not. dir_exists(SIMPLE_BUILD_PATH) .and. &
    (.not. dir_exists(trim(adjustl(SIMPLE_BUILD_PATH))//'/lib/simple'))) then
    THROW_HARD( 'test_args SIMPLE_BUILD_PATH  not found')
    else
        bpath= filepath( trim(adjustl(SIMPLE_BUILD_PATH)),'lib/simple')
    endif
    print *, 'Simple varlist directory ', trim(adjustl(bpath))

    print *, 'appending varlist '
    vfilename = filepath(trim(adjustl(bpath)),trim(adjustl(varlist)))
    vfilename = trim(adjustl(vfilename))
    print *, 'varlist: ', trim(adjustl(vfilename))
    write(*,'(a,a)') '**info(simple_args_unit_test): checking varlist file ',trim(adjustl(vfilename))
    if(.not. file_exists(vfilename))then
        print *,' varlist not in lib/simple/,  checking lib64/simple'
        vfilename = filepath(trim(SIMPLE_BUILD_PATH),'lib64/simple',trim(adjustl(varlist)))
        vfilename = trim(adjustl(vfilename))
        print *, 'varlist: ', trim(adjustl(vfilename))
        write(*,'(a)') '**info(simple_args_unit_test): checking varlist file'
        call simple_file_stat(vfilename,status,buff)
        if(status /= 0)then
            print *,' varlist not in lib64/simple/,  calling simple_args_varlist.pl'
            cmd='cd '//filepath(trim(adjustl(spath)),'src/main')//' && simple_args_varlist.pl'
            call execute_command_line(trim(cmd), wait=.true.)
            call simple_file_stat(vfilename,status,buff)
            if(status /= 0)then
                print *,' varlist still not in lib/simple/ after calling simple_args_varlist.pl'
            endif
        endif
    endif
    n = nlines(vfilename)
    if(n<MINVARS) THROW_HARD( 'test_args SIMPLE_BUILD_PATH  not found')

    call fopen(funit, file=trim(adjustl(vfilename)), status='old', action='read', iostat=io_stat)
    if(io_stat /= 0) then
     call fileiochk(\"simple_args::test  Unable to open \"//trim(adjustl(vfilename)),io_stat)
THROW_HARD( 'simple_args::test  Unable to open varlist ')
end if
    do i=1,n
        read(funit,*) arg
        if( as%is_present(arg) )then
            ! alles gut
        else
            write(*,'(a)') 'this argument should be present: ', arg
            THROW_HARD( 'part 1 of the unit test failed')
        endif
    end do
    call fclose(funit, errmsg=\"simple_args::test close\")
    errarg1 = 'XXXXXXX'
    errarg2 = 'YYYYYY'
    errarg3 = 'ZZZZZZ'
    write(*,'(a)') '**info(simple_args_unit_test, part 2): testing for args that should NOT be present'
    if( as%is_present(errarg1) .or. as%is_present(errarg2) .or. as%is_present(errarg3) )then
        write(*,'(a)') 'the tested argumnets should NOT be present'
        THROW_HARD( 'part 2 of the unit test failed')
    endif
    write(*,'(a)') 'SIMPLE_ARGS_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
end subroutine

end module simple_args\n";
close(MODULE);

# file-handling
if ( -e ${simple_argsfile} && (compare(${simple_argsfile}, ${tmp_argsfile}) == 0)) {
      unlink(${tmp_argsfile}) # do nothing
     }
 else {
    rename  ${tmp_argsfile},${simple_argsfile}
 }
