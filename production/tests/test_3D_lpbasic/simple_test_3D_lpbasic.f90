program simple_test_3D_lpbasic

include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_image,              only: image
use simple_parameters,         only: parameters
implicit none

type(parameters)            :: p
type(cmdline)               :: cline
type(image)                 :: even, odd
character(*), parameter     :: fname = '/Users/wietfeldthc/Documents/simpleTestImages/vol_noisy.mrc'
character(100)              :: resChar
integer                     :: ldim(3) = [256, 256, 256], ifoo
real                        :: smpd = 1.0, res

! Load res using commandline
if( command_argument_count() .NE. 1 )then
    write(logfhandle,'(a)') 'Usage: simple_test_3D_lpbasic res=[real]'
    write(logfhandle,'(a)') 'Example: simple_test_3D_lpbasic res=6.0'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    res = 6.0
else
    call get_command_argument(1, resChar)
    read(resChar, *)res
endif
print *, 'Using a resolution of ', res
!call cline%checkvar('lpthresh',    1)
!call cline%check
!call p%new(cline)

! Load the image.
call find_ldim_nptcls(fname, ldim, ifoo, smpd = smpd)
call even%new(ldim, smpd)
call odd%new(ldim, smpd)
print *, 'Beign reading volume'
call even%read(fname)
call odd%copy(even)

! forward FT
print *, 'Begin FFT'
call even%fft()
call odd%fft()

! Apply the lo pass filter.
print *, 'Begin low pass filter'
call even%lp(calc_fourier_index(res, ldim(1), smpd))

! Inverse FT
print *, 'Begin IFFT'
call even%ifft
call odd%ifft

! Write the clean image and the noisey image to compare
print *, 'Begin writing output'
call even%write('/Users/wietfeldthc/Documents/simpleTestImages/vol_lpfilter.mrc')
call odd%write('/Users/wietfeldthc/Documents/simpleTestImages/vol_nofilter.mrc')

print *, 'Test 3D LP Basic Complete!'
end program simple_test_3D_lpbasic