program simple_test_grad_cartftcc
include 'simple_lib.f08'
use simple_cartft_corrcalc, only: cartft_corrcalc
use simple_eval_cartftcc,   only: eval_cartftcc
use simple_cmdline,         only: cmdline
use simple_builder,         only: builder
use simple_parameters,      only: parameters
use simple_projector_hlev
use simple_timer
use simple_oris
use simple_image
implicit none
type(parameters)         :: p
type(cartft_corrcalc)    :: cftcc
type(cmdline)            :: cline
type(builder)            :: b
integer,     parameter   :: NSPACE=1    ! set to 1 for fast test
real                     :: corrs(NSPACE), corrs2(NSPACE), grad(2, NSPACE)
type(image), allocatable :: imgs(:)
real,        allocatable :: pshifts(:,:)
integer                  :: iref, iptcl, loc(1), cnt, x, y
type(eval_cartftcc)      :: evalcc
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)') 'simple_test_eval_cartftcc lp=xx smpd=yy nthr=zz vol1=vol1.mrc'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('lp',   1)
call cline%checkvar('smpd', 2)
call cline%checkvar('nthr', 3)
call cline%checkvar('vol1', 4)
call cline%set('ctf', 'no')
call cline%set('match_filt', 'no')
call cline%check
call p%new(cline)
p%kfromto(1) = 3
p%kfromto(2) = calc_fourier_index(p%lp, p%box, p%smpd)
p%kstop      = p%kfromto(2)
call b%build_general_tbox(p, cline)
! prep projections
call b%vol%read(p%vols(1))
call b%eulspace%new(NSPACE, is_ptcl=.false.)
call b%eulspace%spiral
p%nptcls = NSPACE
imgs     = reproject(b%vol, b%eulspace)
call b%vol%fft
call b%vol%expand_cmat(KBALPHA) ! necessary for re-projection
! prep correlator
call cftcc%new(p%nptcls, [1, p%nptcls], .false.)
do iref = 1, p%nptcls
    call imgs(iref)%fft
    call cftcc%set_ref(iref, imgs(iref), .true.)
end do
do iptcl = 1,p%nptcls
    call cftcc%set_ptcl(iptcl,imgs(iptcl))
end do
! prep evaluator
call evalcc%new(b%vol, b%vol, NSPACE)
allocate(pshifts(p%nptcls, 2), source=0.)
do iref = 1,p%nptcls
    ! mapping ran3 (0 to 1) to [-5, 5]
    pshifts(iref, 1) = floor(ran3()*10.99) - 5
    pshifts(iref, 2) = floor(ran3()*10.99) - 5
    call evalcc%set_ori(iref, b%eulspace%get_euler(iref), pshifts(iref, :))
end do
cnt = 0
do iptcl = 1,p%nptcls
    call evalcc%project_and_correlate(iptcl, corrs)
    print *, corrs
    call evalcc%project_and_correlate(iptcl, corrs, grad)
    print *, corrs, grad
    loc = maxloc(corrs)
    if( loc(1) == iptcl ) cnt = cnt + 1
    if( .not. loc(1) == iptcl ) write(*, *) pshifts(iptcl, :)
end do
print *, 'initial corr = ', corrs(loc)
! numerical gradient
iptcl = 1
print *, pshifts(iref, 1)
do iref = 1,p%nptcls
    ! dcc/dx
    call evalcc%set_ori(iref, b%eulspace%get_euler(iref), [pshifts(iref, 1) - 0.00001, pshifts(iref, 2)])
    call evalcc%project_and_correlate(iptcl, corrs, grad)
    call evalcc%set_ori(iref, b%eulspace%get_euler(iref), [pshifts(iref, 1) + 0.00001, pshifts(iref, 2)])
    call evalcc%project_and_correlate(iptcl, corrs2, grad)
    print *, 'numerical d_corr/dx = ', (corrs2(1) - corrs(1))/0.00002
    ! dcc/dy
    call evalcc%set_ori(iref, b%eulspace%get_euler(iref), [pshifts(iref, 1), pshifts(iref, 2) - 0.00001])
    call evalcc%project_and_correlate(iptcl, corrs)
    call evalcc%set_ori(iref, b%eulspace%get_euler(iref), [pshifts(iref, 1), pshifts(iref, 2) + 0.00001])
    call evalcc%project_and_correlate(iptcl, corrs2)
    print *, 'numerical d_corr/dy = ', (corrs2(1) - corrs(1))/0.00002
end do
end program simple_test_grad_cartftcc
