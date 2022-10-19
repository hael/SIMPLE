program simple_test_uniform_filter_2D
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_commander_volops,   only: reproject_commander
use simple_image,              only: image
implicit none
type(parameters)              :: p
type(cmdline)                 :: cline, cline_projection
type(image)                   :: img_clean, img_noisy, noise, res_map, img_filt
type(reproject_commander)     :: xreproject
integer                       :: k, nyq, nptcls, smooth_ext, rc, iptcl, ext
real                          :: res_fsc05, res_fsc0143, ave, sdev, maxv, minv, med
real, allocatable             :: res(:), corrs(:)
character(len=:), allocatable :: cmd
logical                       :: mrc_exists
real,    parameter            :: LP_LOWRES_PHASE = 7.
integer, parameter            :: NSEARCH = 1000
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
        call cline_projection%set('nspace'    , 6.)
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
call img_noisy%new(     p%ldim, p%smpd)
call img_clean%new(      p%ldim, p%smpd)
call noise%new(    p%ldim, p%smpd)
do iptcl = 1, p%nptcls
    write(*, *) 'Particle # ', iptcl
    call img_noisy%read(p%stk, iptcl)
    call img_clean%copy(img_noisy)
    ! spherical masking
    call img_noisy%mask(p%msk, 'soft')
    call img_clean%mask( p%msk, 'soft')
    ! forward FT
    call img_noisy%fft()
    call img_clean%fft()
    ! calculate FSC
    res = img_noisy%get_res()
    nyq = img_noisy%get_filtsz()
    if( .not. allocated(corrs) ) allocate(corrs(nyq))
    call img_noisy%fsc(img_clean, corrs)
    call img_noisy%ifft
    call img_clean%ifft
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs clean particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call img_noisy%stats('foreground', ave, sdev, maxv, minv)
    ! add noise in a small center region of the img_noisy
    ! call noise%gauran(0., 5. * sdev)
    ! call noise%mask(0.4 * p%msk, 'soft')
    ! call img_noisy%add(noise)
    ! add background noise
    call noise%gauran(0., 1.5 * sdev)
    call noise%mask(1.4 * p%msk, 'soft')
    call img_noisy%add(noise)
    call img_noisy%write('stk_noisy.mrc', iptcl)
    call img_clean%write('stk_clean.mrc', iptcl)
    call img_noisy%fft
    call img_clean%fft
    call img_noisy%fsc(img_clean, corrs)
    call img_noisy%ifft
    call img_clean%ifft
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs noisy particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call img_noisy%phase_rand(LP_LOWRES_PHASE)
    call img_noisy%write('stk_phase_rand.mrc', iptcl)
    call img_clean%phase_rand(LP_LOWRES_PHASE)
    call img_clean%write('stk_clean_phase_rand.mrc', iptcl)
    call img_noisy%fft
    call img_clean%fft
    call img_noisy%fsc(img_clean, corrs)
    call img_noisy%ifft
    call img_clean%ifft
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs phase-randomized (beyond '//trim(int2str(int((LP_LOWRES_PHASE))))//' A) (noisy) particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call img_clean%zero_and_unflag_ft
    call img_noisy%zero_and_unflag_ft
    call sleep(1)
enddo
call img_noisy%kill()
call img_clean%kill()    
write(*, *) 'Filtering noisy stack in progress...'
call img_noisy%new(p%ldim, p%smpd)
call img_clean%new(p%ldim, p%smpd)
call img_filt %new(p%ldim, p%smpd)
do iptcl = 1, p%nptcls
    write(*, *) 'Particle # ', iptcl
    ! comparing the nonuniform result with the original data
    call img_noisy%read('stk_noisy.mrc', iptcl)
    call img_noisy%uniform_filter_2D(img_filt, NSEARCH)
    call img_filt%write('stk_noisy_filt.mrc', iptcl)
    call img_clean%read('stk_clean.mrc', iptcl)
    ! spherical masking
    call img_filt %mask(p%msk, 'soft')
    call img_clean%mask( p%msk, 'soft')
    ! forward FT
    call img_filt%fft()
    call img_clean%fft()
    ! calculate FSC
    res = img_filt%get_res()
    nyq = img_filt%get_filtsz()
    call img_filt%fsc(img_clean, corrs)
    call img_filt%ifft
    call img_clean%ifft
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs FILTERED noisy particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call img_clean%zero_and_unflag_ft
    call img_noisy%zero_and_unflag_ft
    call img_filt %zero_and_unflag_ft
    call sleep(1)
enddo
call img_noisy%kill()
call img_clean%kill()
call img_filt %kill()
write(*, *) 'Filtering phase-randomized stack in progress...'
call img_noisy%new(p%ldim, p%smpd)
call img_clean%new(p%ldim, p%smpd)
call img_filt %new(p%ldim, p%smpd)
do iptcl = 1, p%nptcls
    write(*, *) 'Particle # ', iptcl
    ! comparing the nonuniform result with the original data
    call img_noisy%read('stk_noisy.mrc', iptcl)
    call img_noisy%uniform_filter_2D(img_filt, NSEARCH)
    call img_filt%write('stk_phase_filt.mrc', iptcl)
    call img_clean%read('stk_clean.mrc', iptcl)
    ! spherical masking
    call img_filt %mask(p%msk, 'soft')
    call img_clean%mask( p%msk, 'soft')
    ! forward FT
    call img_filt %fft()
    call img_clean%fft()
    ! calculate FSC
    res = img_filt%get_res()
    nyq = img_filt%get_filtsz()
    call img_filt%fsc(img_clean, corrs)
    call img_filt%ifft
    call img_clean%ifft
    call get_resolution(corrs, res, res_fsc05, res_fsc0143)
    write(*, *) 'Comparing clean particle vs FILTERED phase-randomized particle...'
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
    write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
    call img_clean%zero_and_unflag_ft
    call img_noisy%zero_and_unflag_ft
    call img_filt %zero_and_unflag_ft
    call sleep(1)
enddo
call img_noisy%kill()
call img_clean%kill()
call img_filt %kill()
end program simple_test_uniform_filter_2D