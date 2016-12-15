program simple_oasis_ctfit
use simple_cmdline      ! singleton
use simple_defs         ! singleton
use simple_build,       only: build
use simple_params,      only: params
use simple_opt_spec,    only: opt_spec
use simple_jiffys,      only: simple_end, progress, alloc_err
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_math,        only: numgrad, deg2rad, rad2deg
use simple_rnd,         only: ran3, seed_rnd
use simple_stat,        only: moment, pearsn
use gnufor2,            only: plot
implicit none
type(params)              :: p
type(build)               :: b
type(opt_factory)         :: ofac_simplex, ofac_oasis
type(opt_spec)            :: spec_simplex, spec_oasis
class(optimizer), pointer :: opt_simplex=>null(), opt_oasis=>null()
integer, parameter        :: NDIM=4
real                      :: limits(NDIM,2), lowest_cost, rtol, gmin, cave, csdev, cvar, sgmin
logical                   :: cyclic(NDIM), err
integer                   :: i, j, success
real, allocatable         :: costs(:)
if( command_argument_count() < 6 )then
     write(*,'(a)', advance='no') 'SIMPLE_OASIS_CTFIT oritab=<SIMPLE alignment doc> box=<image size(in pixels)>'
     write(*,'(a)', advance='no') ' smpd=<sampling distance(in A)> astigerr=<astigmatism error(in microns)> '
     write(*,'(a)', advance='no') '  snr=<signal2noise ratio> lp=<low-pass limit> [defocus=<defocus(in microns)' 
     write(*,'(a)', advance='no') '{3.0}>] [deferr=<defocus error(in microns){1.0}>] [bfac=<bfactor(in A**2){500}>]' 
     write(*,'(a)', advance='no') ' [kv=<acceleration voltage(in kV){300.}>] [fraca=<frac amp contrast{0.07}>] '
     write(*,'(a)', advance='no') '[cs=<spherical aberration constant(in mm){2.7}>] [nrestarts<nr of opt restarts{1}>]'
     write(*,'(a)') ' [hp=<high-pass limit(in A)>]'
     stop
endif
call parse_cmdline
if( .not. defined_cmd_arg('bfac') )then
    call set_cmdline('bfac', 500.)
endif
if( .not. defined_cmd_arg('ctf') )then
    call set_cmdline('ctf', 'yes')
endif
call cmdcheckvar('oritab',   1)
call cmdcheckvar('box',      2)
call cmdcheckvar('smpd',     3)
call cmdcheckvar('astigerr', 4)
call cmdcheckvar('snr',      5)
call cmdcheckvar('lp',       6)
call cmdcheck                              ! checks args and prints cmdline to cmdline.txt
p = params()                               ! constants & derived constants produced
call b%build_general_tbox(p)               ! general objects built
! set limits and cyclic
limits(1:2,1) = max(p%deferr,p%defocus-2.*p%deferr) ! xydefrange 1
limits(1:2,2) = p%defocus+2.*p%deferr               ! xydefrange 2
cyclic(1:2)   = .false.                             ! defocus not cyclic
limits(3,1)   = 0.                                  ! angastrange 1
limits(3,2)   = twopi                               ! angastrange 2
cyclic(3)     = .true.                              ! angast is cyclic
limits(4,1)   = p%bfac-100.                         ! bfacrange 1
limits(4,2)   = p%bfac+100.                         ! bfacrange 2
cyclic(4)     = .false.                             ! bfac not cyclic
!print *, 'defocus:', p%defocus
!print *, 'bfac:', p%bfac
!print *, 'LIMITS'
!print *, 'lim1:', limits(1,1), limits(1,2), cyclic(1)
!print *, 'lim2:', limits(2,1), limits(2,2), cyclic(2)
!print *, 'lim3:', limits(3,1), limits(3,2), cyclic(3)
!print *, 'lim4:', limits(4,1), limits(4,2), cyclic(4)
! specify optimizer
call spec_simplex%specify('simplex', NDIM, ftol=1e-3, gtol=1e-3, limits=limits, cyclic=cyclic, nrestarts=3)
call spec_oasis%specify('oasis',     NDIM, ftol=1e-3, gtol=1e-3, limits=limits, cyclic=cyclic, nrestarts=1)
call spec_simplex%set_costfun(costfun)  ! set pointer to costfun
call spec_oasis%set_costfun(costfun) ! set pointer to costfun
! generate optimizer object with the factory
call ofac_simplex%new(spec_simplex, opt_simplex)
call ofac_oasis%new(spec_oasis, opt_oasis)
! Initialize
call seed_rnd
success = 0
allocate(costs(b%a%get_noris()))

! call spistk_found%open(, 'replace')
sgmin = 0.
do i=1,b%a%get_noris()
    call progress(i,b%a%get_noris())
    ! make noise-corrupted oracle
    call b%tfun%ctf2img(b%img_ctfsq, b%a%get(i,'dfx'), 'square', dfy=b%a%get(i,'dfy'),&
    angast=deg2rad(b%a%get(i,'angast')), bfac=p%bfac)
    b%img_copy = b%img_ctfsq%ft2img('real')
    call b%img_copy%write('refimgs'//p%ext, i)
    call b%img_ctfsq%add_gauran(p%snr)
    b%img_copy = b%img_ctfsq%ft2img('real')
    call b%img_copy%write('refimgs_noisy'//p%ext, i)
    ! calculate cost value at global minimum
    call b%tfun%ctf2img(b%img, b%a%get(i,'dfx'), 'square', dfy=b%a%get(i,'dfy'),&
    angast=deg2rad(b%a%get(i,'angast')), bfac=p%bfac)
    gmin = -b%img%corr(b%img_ctfsq, lp_dyn=p%lp)
    sgmin = sgmin+gmin
    ! generate starting points
    do j=1,NDIM
        spec_simplex%x(j) = spec_simplex%limits(j,1)+ran3()*(spec_simplex%limits(j,2)-spec_simplex%limits(j,1))
    end do
    spec_simplex%x(1) = (spec_simplex%x(2)+spec_simplex%x(2))/2 ! enforce radial symmetry initially
    spec_simplex%x(2) = spec_simplex%x(1)
    ! simplex minimization first
    call opt_simplex%minimize(spec_simplex, lowest_cost)
    spec_oasis%x = spec_simplex%x
    ! then oasis
    call opt_oasis%minimize(spec_oasis, lowest_cost)
    costs(i) = lowest_cost
    ! calculate relative tolerance
    rtol=2.0*abs(gmin-lowest_cost)/(abs(gmin)+abs(lowest_cost)+TINY) ! relative tolerance
    if(rtol <= spec_oasis%ftol) success = success+1
    call b%tfun%ctf2img(b%img, spec_oasis%x(1), 'square', dfy=spec_oasis%x(2),&
    angast=spec_oasis%x(3), bfac=spec_oasis%x(4))
    b%img_copy = b%img%ft2img('real')
    call b%img_copy%write('foundimgs'//p%ext, i)
end do
call moment(costs, cave, csdev, cvar, err)
write(*,'(a,f4.0,a,f4.2,a,f4.2,a,f4.2)') 'successes: ',(100.*real(success ))/real(b%a%get_noris()),&
' gmin: ', sgmin/real(b%a%get_noris()), ' cave: ', cave, ' csdev: ', csdev
call simple_end('**** SIMPLE_OASIS_CTFIT NORMAL STOP ****')

contains
    
    function costfun( x, D ) result( cost )
        integer, intent(in) :: D
        real, intent(in)    :: x(D)
        real :: corr, cost
        ! from parameter vector make CTF image
        call b%tfun%ctf2img(b%img, x(1), 'square', dfy=x(2),&
        angast=x(3), bfac=x(4))
        ! correlate with noise-corrupted oracle
        corr = b%img%corr(b%img_ctfsq, lp_dyn=p%lp)
        ! cost is negative corr
        cost = -corr
    end function
    
end program 
