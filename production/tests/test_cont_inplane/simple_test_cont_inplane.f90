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
integer,            parameter   :: N_PTCLS = 1, N_SAMPLES = 500, N_ITERS_SHC = 100, INPL_ITERS = 100
type(cmdline)           :: cline
type(parameters)        :: p
integer                 :: ifoo, rc, iptcl, iter, isample, cnt, j, iinpl
type(projector)         :: vol_proj
type(sym)               :: pgrpsyms
type(ori)               :: o, o_truth, o_best, o_init, o_ret_1, o_ret_2, o_ret_3
type(image)             :: o_proj
type(cartft_corrcalc)   :: cftcc
logical                 :: mrc_exists
real                    :: corr, p_cur, p_best, corr_best
real, allocatable       :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
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
p%kfromto(1) = 2
p%kfromto(2) = 40
allocate( sigma2_noise(p%kfromto(1):p%kfromto(2), 1:N_PTCLS), source=1. )
call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
call o_proj%new([p%ldim(1), p%ldim(2), 1], p%smpd)
call vol_proj%new(p%ldim, p%smpd)
call vol_proj%read(p%vols(1))
call vol_proj%fft()
call vol_proj%expand_cmat(1.)
call pgrpsyms%new('c1')
call o_truth%new(.true.)
call o%new(.true.)
iptcl = 1 ! working on just one particle
call cftcc%new(vol_proj, vol_proj, [1, N_PTCLS])
call cftcc%assign_sigma2_noise(sigma2_noise)
call pgrpsyms%rnd_euler(o_truth)
call vol_proj%fproject(o_truth, o_proj)
call cftcc%set_ptcl(iptcl, o_proj)
call o_proj%ifft()
call o_proj%write('TEST_PTCL.mrc', 1)
! normal SHC
print *, '-- SHC on corr only --'
call pgrpsyms%rnd_euler(o)
o_init    = o
corr      = cftcc%project_and_correlate(iptcl, o)
corr_best = corr
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
o_ret_1 = o_best
call o_proj%fft()
call vol_proj%fproject(o_ret_1, o_proj)
call o_proj%ifft()
call o_proj%write('TEST_PTCL_1.mrc', 1)
! SHC on {e1,e2,e3} then e3
print *, '-- SHC e3 on corr only --'
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
            ! greedy w.r.t e3
            do iinpl = 1, INPL_ITERS
                call o%rnd_inpl
                corr = cftcc%project_and_correlate(iptcl, o)
                if( corr > corr_best )then
                    corr_best = corr
                    o_best    = o
                endif
            enddo
            print *, 'iter = ', iter, '; corr best = ', corr_best
            exit            ! SHC
        endif
    enddo
enddo
print *, o_best%get_euler()
print *, o_truth%get_euler()
o_ret_2 = o_best
call o_proj%fft()
call vol_proj%fproject(o_ret_2, o_proj)
call o_proj%ifft()
call o_proj%write('TEST_PTCL_2.mrc', 1)
call o_proj%fft()
call vol_proj%fproject(o_init, o_proj)
call o_proj%ifft()
call o_proj%write('TEST_PTCL_INIT.mrc', 1)
! killing
call o_ret_1%kill
call o_ret_2%kill
call vol_proj%kill
call o%kill
call o_truth%kill
call o_best%kill
call o_proj%kill
end program simple_test_cont_inplane