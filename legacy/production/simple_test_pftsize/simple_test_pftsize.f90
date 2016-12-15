program simple_test_pftsize
use simple_defs     ! singleton
use simple_cmdline, only: cmdline
use simple_params,  only: params
use simple_image,   only: image
use simple_math,    only: round2even
use simple_timing
implicit none
type(params)         :: p
type(image)          :: img
type(cmdline)        :: cline
integer(longer)      :: nrnr, sz_all_pfts, nrots, kfromto(2), sz_one_pft, sz_one_real, sz_one_complex
real                 :: mbytes, gbytes, rval
complex, allocatable :: pft(:,:)
complex              :: cval
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'SIMPLE_TEST_PFTSIZE box=<box size(pixels)> smpd=<sampling dist(in A)>'
    write(*,'(a)') ' nspace=<number of refs> lp=<low-pass limit(in A)> [msk=<mask radius(pixels)>]'
    stop
endif
call cline%parse
call cline%checkvar('box',    1)
call cline%checkvar('smpd',   2)
call cline%checkvar('nspace', 3)   
call cline%checkvar('lp',     4)
call cline%check
call cline%set('prg', 'test_pftsize')
p = params(cline) ! parameters generated
sz_one_real    = sizeof(rval)
write(*,'(a,19x,i3)') '>>> SIZE OF ONE REAL:', sz_one_real
sz_one_complex = sizeof(cval)
write(*,'(a,16x,i3)') '>>> SIZE OF ONE COMPLEX:', sz_one_complex
nrots       = round2even(twopi*real(p%ring2))
write(*,'(a,17x,i4)') '>>> NUMBER OF ROTATIONS:', nrots
write(*,'(a,6x,i4)') '>>> NUMBER OF PROJECTION DIRECTIONS:', p%nspace
call img%new([p%box,p%box,1], p%smpd)
kfromto(1)  = 2
kfromto(2)  = img%get_find(p%lp)
write(*,'(a,16x,i3,1x,a,1x,i3)') '>>> FOURIER INDEX RANGE:', kfromto(1), ':', kfromto(2) 
allocate(pft(nrots,kfromto(1):kfromto(2)))
sz_one_pft  = sizeof(pft)
write(*,'(a,12x,i7)') '>>> SIZE OF ONE PFT IN BYTES:', sz_one_pft
! calculate the number of numbers
nrnr = nrots*(kfromto(2)-kfromto(1)+1)*p%nspace ! we would have a factor two but we match only half the transforms
write(*,'(a,1x,i10)') '>>> THE TOTAL NUMBER OF COMPLEX NUMBERS:', nrnr
write(*,'(a,10x,f9.2)') '>>> DIM COMPLEX SQUARE MATRIX:', dsqrt(dble(nrnr))
sz_all_pfts = p%nspace*sz_one_pft ! we would have a factor two but we match only half the transforms
write(*,'(a,10x,i12)') '>>> SIZE OF ALL PFTS IN BYTES:', sz_all_pfts
write(*,'(a,10x,i12)') '>>> SIZE OF ALL PFTS IN BYTES:', sz_one_complex*nrnr
mbytes      = real(dble(sz_all_pfts)/dble(1e6))
write(*,'(a,6x,f9.2)') '>>> SIZE OF ALL PFTS IN MEGABYTES:', mbytes
gbytes      = real(dble(sz_all_pfts)/dble(1e9))
write(*,'(a,6x,f9.2)') '>>> SIZE OF ALL PFTS IN GIGABYTES:', gbytes
end program
