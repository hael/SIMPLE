program simple_test_2D_opt_filt
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_commander_resolest, only: opt_2D_filter_commander
use simple_commander_volops,   only: reproject_commander
use simple_image,              only: image
implicit none
type(parameters)              :: p
type(cmdline)                 :: cline, cline_opt_filt, cline_projection
type(image)                   :: even, odd, noise, res_map
type(opt_2D_filter_commander) :: xopt_2D_filter
type(reproject_commander)     :: xreproject
integer                       :: j, nyq, nptcls, smooth_ext, rc, iptcl
real                          :: res_fsc05, res_fsc0143, ave, sdev, maxv, minv, med, lp_lb
real, allocatable             :: res(:), corrs(:)
character(len=20)             :: filter
character(len=:), allocatable :: cmd
logical                       :: mrc_exists
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_2D_opt_filt smpd=xx nthr=yy stk=stk.mrc, mskdiam=zz'
    write(logfhandle,'(a)') 'Example: projections of https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
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
        write(*, *) 'Projecting 1JYX.mrc...'
        call cline_projection%set('vol1'      , '1JYX.mrc')
        call cline_projection%set('smpd'      , 1.)
        call cline_projection%set('pgrp'      , 'c1')
        call cline_projection%set('mskdiam'   , 180.)
        call cline_projection%set('nspace'    , 50.)
        call cline_projection%set('nthr'      , 16.)
        call xreproject%execute(cline_projection)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('stk'    , 'reprojs.mrcs')
    call cline%set('mskdiam', 180.)
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('stk' ,    3)
call cline%checkvar('mskdiam', 4)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%stk, p%ldim, nptcls)
p%ldim(3) = 1 ! because we operate on stacks
call even%new(     p%ldim, p%smpd)
call odd%new(      p%ldim, p%smpd)
call noise%new(    p%ldim, p%smpd)
do iptcl = 1, p%nptcls
    write(*, *) 'Particle # ', iptcl
    call even%read(p%stk, iptcl)
    call odd%copy      (even)
    ! spherical masking
    call even%mask(p%msk, 'soft')
    call odd%mask( p%msk, 'soft')
    ! forward FT
    call even%fft()
    call odd%fft()
    ! calculate FSC
    res = even%get_res()
    nyq = even%get_filtsz()
    if( .not. allocated(corrs) ) allocate(corrs(nyq))
    call even%fsc(odd, corrs)
    call even%ifft
    call odd%ifft
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs clean particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call even%stats('foreground', ave, sdev, maxv, minv)
    ! add noise in a small center region of the even
    call noise%gauran(0., 5. * sdev)
    call noise%mask(0.4 * p%msk, 'soft')
    call even%add(noise)
    call even%write('vol_noisy.mrc', iptcl)
    call odd%write('vol_clean.mrc', iptcl)
    call even%fft
    call odd%fft
    call even%fsc(odd, corrs)
    call even%ifft
    call odd%ifft
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs noisy particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call odd%zero_and_unflag_ft
    call even%zero_and_unflag_ft
    call sleep(1)
enddo
! calculate optimal filter
write(*, *) 'Filtering noisy stack in progress...'
cline_opt_filt = cline
filter         = 'lp'
smooth_ext     = 1
if( cline%defined('smooth_ext') ) smooth_ext = p%smooth_ext
if( cline%defined('filter') )     filter     = p%filter
if( .not. cline%defined('match_filt') ) call cline_opt_filt%set('match_filt', 'no')
call cline_opt_filt%set('stk'       , 'vol_noisy.mrc')
call cline_opt_filt%set('stk2'      , 'vol_clean.mrc')
call cline_opt_filt%set('filter'    , trim(filter))
call cline_opt_filt%set('mkdir'     , 'no')
call cline_opt_filt%set('smooth_ext', real(smooth_ext))
call p%new(cline_opt_filt)
call xopt_2D_filter%execute(cline_opt_filt)
do iptcl = 1, p%nptcls
    write(*, *) 'Particle # ', iptcl
    ! comparing the nonuniform result with the original data
    call even%read(p%stk2,  iptcl)
    call odd%read('nonuniform_opt_2D_filter_'//trim(filter)//'_ext_'//int2str(smooth_ext)//'_odd.mrc', iptcl)
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
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs FILTERED noisy particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call odd%zero_and_unflag_ft
    call even%zero_and_unflag_ft
    call sleep(1)
enddo
end program simple_test_2D_opt_filt