program simple_test_cont_inplane
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
integer,            parameter   :: N_PTCLS = 1, N_SAMPLES = 1000, N_ITERS_SHC = 150, N_ITERS_PROB = N_ITERS_SHC
type(cmdline)           :: cline
type(parameters)        :: p
integer                 :: ifoo, rc, iptcl, iter, isample, cnt, j
type(projector)         :: vol_proj
type(sym)               :: pgrpsyms
type(ori)               :: o, o_truth, o_arr(N_ITERS_PROB+1), o_best, o_init
type(image)             :: o_proj
type(cartft_corrcalc)   :: cftcc
logical                 :: mrc_exists
real                    :: corr, corr_arr(N_ITERS_PROB+1), p_cur, p_best, corr_best
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_cont_inplane smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
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
iptcl        = 1 ! working on just one particle
call cftcc%new(vol_proj, vol_proj, [1, N_PTCLS])
call pgrpsyms%rnd_euler(o_truth)
call vol_proj%fproject(o_truth, o_proj)
call cftcc%set_ptcl(iptcl, o_proj)
call o_proj%ifft()
call o_proj%write('TEST_PTCL.mrc', 1)
! normal SHC
print *, '-- SHC on corr only --'
o         = o_init
corr_best = 0.
o_best    = o
do iter = 1, N_ITERS_SHC
    o = o_best
    do isample = 1, N_SAMPLES
        call pgrpsyms%rnd_euler(o)
        corr = cftcc%project_and_correlate(iptcl, o)
        if( corr > corr_best )then
            corr_best = corr
            o_best    = o
            print *, 'iter = ', iter, '; corr best = ', corr_best
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
end program simple_test_cont_inplane