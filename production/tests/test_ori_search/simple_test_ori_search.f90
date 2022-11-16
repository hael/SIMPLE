program simple_test_ori_search
include 'simple_lib.f08'
use simple_image,             only: image
use simple_parameters,        only: parameters
use simple_cmdline,           only: cmdline
use simple_cartft_corrcalc,   only: cartft_corrcalc
use simple_projector
use simple_sym
use simple_ori
implicit none
character(len=:),   allocatable :: cmd
integer,            parameter   :: N_PTCLS = 1, N_SAMPLES = 1000, N_ITERS = 200
type(cmdline)           :: cline
type(parameters)        :: p
integer                 :: ifoo, rc, iptcl, iter, isample, cnt
type(projector)         :: vol_proj
type(sym)               :: pgrpsyms
type(ori)               :: o, o_truth, o_arr(N_ITERS), o_best, o_init
type(image)             :: o_proj
type(cartft_corrcalc)   :: cftcc
logical                 :: mrc_exists, l_match_filt
real                    :: corr, corr_arr(N_ITERS), p_cur, p_best
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_ori_search smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
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
call o_proj%new([p%ldim(1), p%ldim(2), 1], p%smpd)
call vol_proj%new(p%ldim, p%smpd)
call vol_proj%read(p%vols(1))
call vol_proj%fft()
call vol_proj%expand_cmat(1.)
call pgrpsyms%new('c1')
call o_truth%new(.true.)
call o%new(.true.)
l_match_filt = (p%match_filt .eq. 'yes')
iptcl        = 1 ! working on just one particle
call cftcc%new(vol_proj, vol_proj, [1, N_PTCLS], l_match_filt=l_match_filt)
call pgrpsyms%rnd_euler(o_truth)
call vol_proj%fproject(o_truth, o_proj)
call cftcc%set_ptcl(iptcl, o_proj)
call o_proj%write('TEST_PTCL.mrc', 1)
! randomizing the orientation and compare with the truth
print *, '-- SHC on probability --'
call pgrpsyms%rnd_euler(o)
o_init = o
corr   = cftcc%project_and_correlate(iptcl, o)
o_arr(1)    = o
corr_arr(1) = corr
p_best      = 0.
cnt         = 1
do iter = 2, N_ITERS
    do isample = 1, N_SAMPLES
        call pgrpsyms%rnd_euler(o)
        corr  = cftcc%project_and_correlate(iptcl, o)
        p_cur = cftcc%ori_chance( iptcl, o, o_arr, corr_arr, R = 100., n = cnt )
        if( p_cur > p_best )then
            cnt           = cnt + 1
            corr_arr(cnt) = corr
            o_arr(   cnt) = o
            p_best        = p_cur
            o_best        = o
            print *, 'iter = ', iter, '; p_best = ', p_best, '; angle_diff = ', rot_angle(o_best%get_mat(), o_truth%get_mat())
            exit
        endif
    enddo       
enddo
print *, o_best%get_euler()
print *, o_truth%get_euler()
! normal SHC
print *, '-- SHC on corr only --'
o      = o_init
p_best = 0.
do iter = 1, N_ITERS
    do isample = 1, N_SAMPLES
        call pgrpsyms%rnd_euler(o)
        corr = cftcc%project_and_correlate(iptcl, o)
        if( corr > p_best )then
            p_best = corr
            o_best = o
            print *, 'iter = ', iter, '; corr best = ', p_best, '; angle_diff = ', rot_angle(o_best%get_mat(), o_truth%get_mat())
            exit            ! SHC
        endif
    enddo
enddo
print *, o_best%get_euler()
print *, o_truth%get_euler()
call vol_proj%kill
call o%kill
call o_truth%kill
call o_proj%kill
end program simple_test_ori_search
