program simple_test_speed
include 'simple_lib.f08'
use simple_ssim
use simple_image, only: image
use simple_cmdline, only: cmdline
implicit none
type(cmdline)     :: cline
type(image)       :: img1,img2
real              :: starttime, stoptime, smpd
logical           :: be_verbose=.false.
character(len=STDLEN) :: oldCWDfolder,curDir,datestr
character(len=STDLEN) :: simple_testbench,timestr,folder, fileref, filecompare
integer :: io_stat, ldim1(3),ldim2(3), nptcls1, nptcls2
integer(8) :: count1, count2, clockrate
real :: ssim_result
smpd = 1.

call date_and_time(date=datestr)
folder = trim('./SIMPLE_TEST_SSIM_'//datestr)
call simple_mkdir( trim(folder) , status=io_stat)
if(io_stat/=0) call simple_stop("simple_mkdir failed")
print *," Changing directory to ", folder
call simple_chdir( trim(folder),  oldCWDfolder , status=io_stat)
if(io_stat/=0) call simple_stop("simple_chdir failed")
call simple_getcwd(curDir)
print *," Current working directory ", curDir
count1=tic()
! if( command_argument_count() < 2 )then
!     write(*,'(a)') 'simple_test_ssim nthr=<number of threads> stk1=<file> stk2=<file> [verbose=<yes|no{no}>]'
!     stop
! endif
! call cline%parse_oldschool
! call cline%checkvar('nthr', 1)
! call cline%check
! be_verbose = .false.
! if( cline%defined('verbose') )then
!     if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
!         be_verbose = .true.
!     endif
! endif
io_stat =  simple_getenv("SIMPLE_TESTBENCH_DATA", simple_testbench)
simple_testbench = trim(adjustl(simple_testbench))
if( io_stat /= 0 ) call simple_stop(' Cannot test ssim without SIMPLE_TESTBENCH_DATA environment variable set')
write (*,*) " Checking testbench data folder :"//trim(simple_testbench)

fileref = trim(simple_testbench)//'/recvol_1wcm/reference/recvol_state01.mrc'
filecompare = trim(simple_testbench)//'/recvol_1wcm/reconstruct3D_test/1_reconstruct3D/recvol_state01.mrc'
if(.not. file_exists (trim(fileref)))&
    call simple_stop(' Cannot test ssim without valid SIMPLE_TESTBENCH_DATA directory. File not found '//&
    &trim(fileref) )
if(.not.  file_exists (trim(filecompare)))then
    call execute_command_line("../../benchmarks/test_recvol")
end if
if(.not.  file_exists (trim(filecompare)))&
    call simple_stop(' Cannot test ssim -- test_recvol failed to create recvol_state01.mrc')

call find_ldim_nptcls(trim(fileref), ldim1, nptcls1)
call find_ldim_nptcls(trim(filecompare), ldim2, nptcls2)
if ( nptcls1 /= nptcls2 .or. any(ldim1 /= ldim2)) &
    call simple_stop ('files mismatch')
print *,' File 1: ', ldim1
print *,' File 2: ', ldim2

call img1%new(ldim1, smpd)
call img2%new(ldim2, smpd)
write(*,'(a,1x,a)') " Reading ", trim(fileref)
call img1%read(trim(fileref))
write(*,'(a,1x,a)') " Reading ", trim(filecompare)
call img2%read(trim(filecompare))
write(*,'(a,1x,f9.2,1x,f9.2)') " Mean of images  ", img1%mean(), img2%mean()
ssim_result = ssim(img1,img2)
write(*,'(a,1x,f9.2)') " Structural-similarity metric between images is ", ssim_result
write(*,'(a,1x,f9.2)') 'time elapsed (s): ', toc(count1)
end program simple_test_speed
