program simple_test_picker_comp
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_picker_utils
use simple_parameters

implicit none
character(len = 200) :: micname 
real, parameter :: SMPD = 1.72, MOLDIAM = 150. ! edit these based on your data
integer :: nthr, nptcls, iostat
real :: moldiams_optimal(2)
character(len=LONGSTRLEN) :: cwd
character(len=LONGSTRLEN) :: boxfile_out
type(parameters), target :: params

!$ nthr = omp_get_max_threads()
!$ call omp_set_num_threads(nthr)
nthr_glob = nthr

call simple_getcwd(cwd)
allocate(CWD_GLOB, source=trim(cwd))

! set up parameters
params_glob => params
params_glob%moldiam = MOLDIAM
params_glob%pcontrast = 'black' ! not sure if you want black or white contrast?

! Testing to make sure the file opens correctly
! PUT THE NAME OF YOUR FILE INSIDE THE TRIM FUNCTION
micname = trim( CWD_GLOB // trim('/SLC22A6_intg_1.mrc'))

open(unit=16,file=micname,iostat=iostat,status='old')
if (iostat .eq. 0) then
    print *, 'File ' // trim(micname) // ' exists'
else
    print *, 'An error occured when trying to open the file'
end if
close(16)

! perform picking
! see simple_picker_utils for details on how this works 
call exec_gaupick(micname=micname, boxfile_out=boxfile_out, smpd=SMPD, nptcls=nptcls)


end program simple_test_picker_comp
