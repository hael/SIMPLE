! This program is to illustrate a low-pass filter on a noisy 3D volume
! as well as the structure of a Fortran program
program simple_test_3D_lpbasic

! Include the following libraries and modules
include 'simple_lib.f08'
use simple_image,              only: image

implicit none  ! Ensures all variables must have an explicit type (all modern Fortran programs should include this statement)

! Declaration of variables
type(image)                 :: denoised
real                        :: smpd=1.0      ! Sampling distance (Angstroms per voxel length)
real                        :: res           ! Resolution of output image
integer                     :: ldim(3), ifoo ! ldim is the dimensions of the 3D volume
character(*), parameter     :: fn_in = '/Users/wietfeldthc/Documents/simpleTestImages/lpForAlannaJun27/vol_noisy.mrc'
character(256)              :: fn_out, resChar

! Load the resolution from the commandline
if( command_argument_count() /= 2 ) then
    write(logfhandle,'(a)') 'Usage: simple_test_3D_lpbasic res outvol.mrc'
    write(logfhandle,'(a)') 'res: resolution of denoised image [real]'
    write(logfhandle,'(a)') 'outvol.mrc: name of output volume with .mrc extension'
    write(logfhandle,'(a)') 'Example: simple_test_3D_lpbasic 6.0 example_out.mrc'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    res = 6.0
    fn_out = 'example_out.mrc'
else
    call get_command_argument(1, resChar)
    read(resChar, *)res ! Converts character to real
    call get_command_argument(2, fn_out)
endif
! '(a,f5.2,a)' means format the output with three components: a string, a floating-point number with
! a length of 5 characters and a precision of 2 decimal points, and another string.
write(logfhandle,'(a,f5.2,a)') 'Using a resolution of ', res, 'Å'

! Load the image
call find_ldim_nptcls(fn_in, ldim, ifoo, smpd = smpd) ! ifoo is not used but must be included in the function call
call denoised%new(ldim, smpd)
print *, 'Reading volume'
call denoised%read(fn_in)

! Denoise the image using a low-pass filter
print *, 'Applying Fast Fourier Transform'
call denoised%fft()
print *, 'Applying low-pass filter'
call denoised%lp(calc_fourier_index(res, ldim(1), smpd))
print *, 'Applying Inverse Fast Fourier Transform'
call denoised%ifft()

! Write the denoised image
print *, 'Writing output'
call denoised%write(fn_out)

! Cleanup allocated memory (not necessary in this case but good practice)
call denoised%kill()

print *, 'Test 3D LP Basic Complete!'
end program simple_test_3D_lpbasic