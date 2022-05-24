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
integer                       :: j, nyq, ifoo, smooth_ext
real                          :: res_fsc05, res_fsc0143, ave, sdev, maxv, minv, med, lp_lb
real, allocatable             :: res(:), corrs(:)
character(len=20)             :: filter
real, pointer                 :: res_rmat(:, :, :)=>null()
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'simple_test_3D_opt_filt smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
    stop
endif
call cline%parse_oldschool
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
do j=1,nyq
    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
end do
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
! phase randomize from 7A
call even%stats('foreground', ave, sdev, maxv, minv)
! add noise in a small center region of the even
call noise%gauran(0., 15. * sdev)
call noise%mask(0.4 * p%msk, 'soft')
call even%add(noise)
lp_lb = 7.
call even%phase_rand(lp_lb)
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
do j=1,nyq
    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
end do
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
! calculate optimal filter
cline_opt_filt = cline
filter         = 'lp'
smooth_ext     = 1
if( cline%defined('smooth_ext') ) smooth_ext = p%smooth_ext
if( cline%defined('filter') )     filter     = p%filter
call cline_opt_filt%set('vol1',   'vol_phase_rand.mrc')
call cline_opt_filt%set('vol2',   'vol_clean.mrc')
call cline_opt_filt%set('filter', trim(filter))
call cline_opt_filt%set('mkdir',  'no')
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
do j=1,nyq
    write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', corrs(j)
end do
call get_resolution(corrs, res, res_fsc05, res_fsc0143)
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
call res_map%new(p%ldim, p%smpd)
call res_map%read('opt_resolution_odd_map_nonuniform_filter_'//trim(filter)//'_ext_'//int2str(smooth_ext)//'.mrc')
call res_map%stats( 'foreground', ave, sdev, maxv, minv, p%msk, med )
write(*, *) 'resolution map: average = ', ave, '; max = ', maxv, '; min = ', minv, '; median = ', med
end program simple_test_3D_opt_filt