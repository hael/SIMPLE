program simple_test_project_shift
include 'simple_lib.f08'
use simple_image,             only: image
use simple_parameters,        only: parameters
use simple_cmdline,           only: cmdline
use simple_cartft_corrcalc,   only: cartft_corrcalc
use simple_cftcc_shsrch_grad, only: cftcc_shsrch_grad
use simple_projector
use simple_sym
use simple_ori
implicit none
character(len=:),   allocatable :: cmd
integer,            parameter   :: N_PTCLS = 10
type(cmdline)           :: cline
type(parameters)        :: p
integer                 :: ifoo, rc, i, nevals(2)
type(projector)         :: vol_proj
type(sym)               :: pgrpsyms
type(ori)               :: o, o_cur
type(image)             :: o_proj
type(cartft_corrcalc)   :: cftcc
type(cftcc_shsrch_grad) :: cftcc_shsrch
logical                 :: mrc_exists
real                    :: corr, euls(3), lims(2,2), cxy(3)
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_project_shift smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
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
call o%new(.true.)
call cftcc%new(vol_proj, vol_proj, [1, N_PTCLS], l_match_filt=.false.) ! 2 particles for now
do i = 1, N_PTCLS
    call pgrpsyms%rnd_euler(o)
    call vol_proj%fproject(o, o_proj)
    call cftcc%set_ptcl(i, o_proj)
    call o_proj%ifft()
    ! call o_proj%write('TEST.mrc', i)
enddo
o_cur = o
corr  = cftcc%project_and_correlate(N_PTCLS, o_cur)
print *, 'corr comparing projection against itself: ', corr
euls(1) = o_cur%e1get() + 3.
euls(2) = o_cur%e2get() - 3.
euls(3) = o_cur%e3get() + 3.
call o%set_euler(euls)
corr = cftcc%project_and_correlate(N_PTCLS, o)
print *, 'corr comparing projection against itself at 3 angerr: ', corr
corr = cftcc%project_and_correlate(N_PTCLS, o, [-3., 3.])
print *, 'corr comparing projection against itself at 3 angerr, 3 sherr: ', corr
euls(1) = o_cur%e1get() - 5.
euls(2) = o_cur%e2get() + 5.
euls(3) = o_cur%e3get() - 5.
call o%set_euler(euls)
corr = cftcc%project_and_correlate(N_PTCLS, o)
print *, 'corr comparing projection against itself at 5 angerr: ', corr
corr = cftcc%project_and_correlate(N_PTCLS, o, [-1., 1.])
print *, 'corr comparing projection against itself at 5 angerr, 1 sherr: ', corr
lims(1,1) = -6.
lims(1,2) =  6.
lims(2,1) = -6.
lims(2,2) =  6.
call cftcc_shsrch%new(lims)
call cftcc_shsrch%set_pind(N_PTCLS)
call o_proj%fft()
call o_proj%shift2Dserial([-2., 2.]) 
call cftcc%set_ref(o_proj)
cxy = cftcc_shsrch%minimize(nevals)
print *, 'shsrch corr: ', cxy(1), ' at shift = ', cxy(2:3)
call vol_proj%kill
call o%kill
call o_proj%kill
end program simple_test_project_shift
