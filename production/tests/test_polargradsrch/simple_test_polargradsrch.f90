program simple_test_polargradsrch
include 'simple_lib.f08'
use simple_polarft_corrcalc,    only: polarft_corrcalc
use simple_cmdline,             only: cmdline
use simple_builder,             only: builder, build_glob
use simple_parameters,          only: parameters
use simple_pftcc_orisrch_grad,  only: pftcc_orisrch_grad
use simple_strategy2D3D_common, only: set_bp_range
use simple_ori,                 only: ori
implicit none
type(parameters)         :: p
type(polarft_corrcalc)   :: pftcc
type(cmdline)            :: cline
type(builder)            :: b
type(pftcc_orisrch_grad) :: orisrch
type(ori)                :: o
real                     :: cxy(3), shvec(2), shift_err, ang_err
logical :: found_better
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_srch vol1=xx mskdiam=xx lp=xx'
    write(logfhandle,'(a)') ' smpd=xx>'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('vol1', 1)
call cline%checkvar('mskdiam',  2)
call cline%checkvar('smpd', 3)
call cline%checkvar('lp',   4)
call cline%set('nptcls',1.0)
call cline%set('ctf','no')
call cline%check
call b%init_params_and_build_strategy3D_tbox(cline,p)
call set_bp_range(cline)

ang_err   = 16.
shift_err = 8.
call b%eulspace%get_ori(irnd_uni(p%nspace), o)
print *,'Ground truth:'
call o%print_ori
print *,'Shift= 0.0 0.0'
print *,'---------------------'

call pftcc%new(p%nptcls, [1, p%nptcls], .false.)
call b%vol%read(p%vols(1))
call b%vol%mask(p%msk,'soft')
if( p%gridding.eq.'yes' ) call b%vol%div_w_instrfun(p%interpfun, alpha=p%alpha)
call b%vol%fft()
call b%vol%expand_cmat(p%alpha,norm4proj=.true.)
call b%vol%fproject_polar(1, o, pftcc, iseven=.true., mask=b%l_resmsk)
call pftcc%cp_even_ref2ptcl(1,1)
call pftcc%set_eo(1, .true. )


call o%e1set(o%e1get() + 2.*(ran3()-0.5)*ang_err)
call o%e2set(o%e2get() + 2.*(ran3()-0.5)*ang_err)
call o%e3set(o%e3get() + 2.*(ran3()-0.5)*ang_err)
shvec(1) = 2.*(ran3()-0.5)*shift_err
shvec(2) = 2.*(ran3()-0.5)*shift_err
print *,'Starting point:'
call o%print_ori
print *,'Shift= ',shvec
print *,'---------------------'

call pftcc%shift_ptcl(1,shvec)
call pftcc%memoize_ffts

call orisrch%new
call orisrch%set_particle(1)
cxy = orisrch%minimize(o, ang_err, shift_err, found_better)
print *,'Solution:'
call o%print_ori
print *,'cc   = ',cxy(1),found_better



end program simple_test_polargradsrch
