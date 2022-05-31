program simple_test_3D_opt_filt
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_commander_resolest, only: opt_3D_filter_commander
use simple_image,              only: image
implicit none
type(parameters)              :: p
type(cmdline)                 :: cline, cline_opt_filt
type(image)                   :: even, odd, even_copy, odd_copy, noise, res_map
type(opt_3D_filter_commander) :: xopt_3D_filter
integer                       :: j, nyq, ifoo, smooth_ext, rc
real                          :: res_fsc05, res_fsc0143, ave, sdev, maxv, minv, med
real, allocatable             :: res(:), corrs(:)
character(len=20)             :: filter
character(len=:), allocatable :: cmd
logical                       :: mrc_exists
real, parameter               :: LP_LB_PHASE = 7.
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_3D_opt_filt smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
    write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    inquire(file="1JYX.mrc", exist=mrc_exists)
    if( .not. mrc_exists )then
        write(*, *) 'Downloading the example dataset...'
        cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
        write(*, *) 'Converting .pdb to .mrc...'
        cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
        call execute_command_line(cmd, exitstat=rc)
        cmd = 'rm 1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('vol1'   , '1JYX.mrc')
    call cline%set('mskdiam', 180.)
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('vol1',    3)
call cline%checkvar('mskdiam', 4)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
call even%new(     p%ldim, p%smpd)
call odd%new(      p%ldim, p%smpd)
call noise%new(    p%ldim, p%smpd)
call even_copy%new(p%ldim, p%smpd)
call odd_copy%new( p%ldim, p%smpd)
call even%read(p%vols(1))
call odd%copy      (even)
call even_copy%copy(even)
call odd_copy%copy (odd)
! spherical masking
call even%mask(p%msk, 'soft')
call odd%mask( p%msk, 'soft')
! forward FT
call even%fft()
call odd%fft()
! calculate FSC
res = even%get_res()
nyq = even%get_filtsz()
allocate(corrs(nyq))
call even%fsc(odd, corrs)
call even%ifft
call odd%ifft
!do j=1,nyq
!    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
!end do
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(*, *) 'Comparing clean volume vs clean volume...'
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
! phase randomize from 7A
call even%stats('foreground', ave, sdev, maxv, minv)
! add noise in a small center region of the even
call noise%gauran(0., 15. * sdev)
call noise%mask(0.4 * p%msk, 'soft')
call even%add(noise)
call even%write('vol_noisy.mrc')
call even%fft
call odd%fft
call even%fsc(odd, corrs)
call even%ifft
call odd%ifft
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(*, *) 'Comparing clean volume vs noisy volume...'
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
call even%phase_rand(LP_LB_PHASE)
call even%write('vol_phase_rand.mrc')
call odd%write('vol_clean.mrc')
! calculate FSC
call even%fft
call odd%fft
call even%fsc(odd, corrs)
call even%ifft
call odd%ifft
call even%kill
call odd%kill
call noise%kill
!do j=1,nyq
!    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
!end do
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(*, *) 'Comparing clean volume vs phase-randomized (beyond '//trim(int2str(int((LP_LB_PHASE))))//' A) (noisy) volume...'
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
! calculate optimal filter
write(*, *) 'Filtering noisy volume in progress...'
cline_opt_filt = cline
filter         = 'butterworth'
smooth_ext     = 1
if( cline%defined('smooth_ext') ) smooth_ext = p%smooth_ext
if( cline%defined('filter') )     filter     = p%filter
if( .not. cline%defined('match_filt') ) call cline_opt_filt%set('match_filt', 'no')
call cline_opt_filt%set('vol1'      , 'vol_noisy.mrc')
call cline_opt_filt%set('vol2'      , 'vol_clean.mrc')
call cline_opt_filt%set('filter'    , trim(filter))
call cline_opt_filt%set('mkdir'     , 'no')
call cline_opt_filt%set('smooth_ext', real(smooth_ext))
call p%new(cline_opt_filt)
call xopt_3D_filter%execute(cline_opt_filt)
! comparing the nonuniform result with the original data
call even%new(p%ldim, p%smpd)
call odd%new( p%ldim, p%smpd)
call even%copy(even_copy)
call odd%read('nonuniform_opt_3D_filter_'//trim(filter)//'_ext_'//int2str(smooth_ext)//'_odd.mrc')
! spherical masking
call even%mask(p%msk, 'soft')
call odd%mask( p%msk, 'soft')
! forward FT
call even%fft()
call odd%fft()
! calculate FSC
res = even%get_res()
nyq = even%get_filtsz()
call even%fsc(odd, corrs)
call even%ifft
call odd%ifft
!do j=1,nyq
!    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
!end do
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(*, *) 'Comparing clean volume vs FILTERED noisy volume...'
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
write(*, *) 'Filtering phase-randomized volume in progress...'
cline_opt_filt = cline
filter         = 'butterworth'
smooth_ext     = 0
if( cline%defined('smooth_ext') ) smooth_ext = p%smooth_ext
if( cline%defined('filter') )     filter     = p%filter
if( .not. cline%defined('match_filt') ) call cline_opt_filt%set('match_filt', 'no')
call cline_opt_filt%set('vol1'      , 'vol_phase_rand.mrc')
call cline_opt_filt%set('vol2'      , 'vol_clean.mrc')
call cline_opt_filt%set('filter'    , trim(filter))
call cline_opt_filt%set('mkdir'     , 'no')
call cline_opt_filt%set('smooth_ext', real(smooth_ext))
call p%new(cline_opt_filt)
call xopt_3D_filter%execute(cline_opt_filt)
! comparing the nonuniform result with the original data
call even%new(p%ldim, p%smpd)
call odd%new( p%ldim, p%smpd)
call even%copy(even_copy)
call odd%read('nonuniform_opt_3D_filter_'//trim(filter)//'_ext_'//int2str(smooth_ext)//'_odd.mrc')
! spherical masking
call even%mask(p%msk, 'soft')
call odd%mask( p%msk, 'soft')
! forward FT
call even%fft()
call odd%fft()
! calculate FSC
res = even%get_res()
nyq = even%get_filtsz()
call even%fsc(odd, corrs)
call even%ifft
call odd%ifft
!do j=1,nyq
!    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
!end do
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(*, *) 'Comparing clean volume vs FILTERED phase-randomized volume...'
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
call res_map%new(p%ldim, p%smpd)
call res_map%read('opt_resolution_odd_map_nonuniform_filter_'//trim(filter)//'_ext_'//int2str(smooth_ext)//'.mrc')
call res_map%stats( 'foreground', ave, sdev, maxv, minv, p%msk, med )
write(*, *) 'resolution map: average = ', ave, '; max = ', maxv, '; min = ', minv, '; median = ', med
end program simple_test_3D_opt_filt