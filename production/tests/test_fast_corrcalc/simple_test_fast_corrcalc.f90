program simple_test_fast_corrcalc
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_cmdline,          only: cmdline
use simple_builder,          only: builder
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_polarizer,        only: polarizer
use simple_pftcc_shsrch_grad,only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_timer
implicit none
type(cmdline)          :: cline
type(builder)          :: b
type(parameters)       :: p
type(polarft_corrcalc) :: pftcc
type(polarizer)        :: img_copy, img
type(pftcc_shsrch_grad):: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
integer(timer_int_kind)        :: t
real(timer_int_kind)        :: rt

logical                :: be_verbose=.false.
real,    parameter     :: SHMAG=1.0
integer, parameter     :: N_PTCLS = 1
real,    allocatable   :: corrs(:), corrs2(:), sig(:), sig2(:)
real                   :: shift(2), corrmax, corr, cxy(3), lims(2,2), rmsd, df
real(dp)               :: f, grad(2), f_err, gradx_err, grady_err
real(dp)               :: f_dev, grad_dev(2)
integer                :: xsh, ysh, xbest, ybest, i, irot, nrots
real, allocatable      :: sigma2_noise(:,:)
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_fast_corrcalc stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
    write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('stk',      1)
call cline%checkvar('mskdiam',  2)
call cline%checkvar('smpd',     3)
call cline%check
be_verbose = .false.
if( cline%defined('verbose') )then
    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
        be_verbose = .true.
    endif
endif
! call cline%set('objfun','cc')
call p%new(cline)
p%kfromto(1) = 1
p%kfromto(2) = 60
allocate( sigma2_noise(p%kfromto(1):p%kfromto(2), 1:N_PTCLS), source=1.0 )
call b%build_general_tbox(p, cline)
call pftcc%new(N_PTCLS, [1,N_PTCLS], p%kfromto)
call pftcc%assign_sigma2_noise(sigma2_noise)
nrots = pftcc%get_nrots()
allocate(corrs(nrots),corrs2(nrots))
call img_copy%init_polarizer(pftcc, p%alpha)
call b%img%read(p%stk, 1)
call b%img%norm
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 1, isptcl=.false., iseven=.true., mask=b%l_resmsk)

! call b%img%add_gauran(1.)
! call b%img%norm
! img_copy = b%img
! call img_copy%fft()
! call img_copy%polarize(pftcc, 1, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
pftcc%pfts_ptcls(:,:,1) = pftcc%pfts_refs_even(:,:,1)

! dummy ctf matrix
pftcc%with_ctf = .true.
! pftcc%with_ctf = .false.
df = 2.0
! call pftcc%calc_polar_ctf(1, p%smpd,200.,2.7,0.1, df,df, 0.178)
call pftcc%calc_polar_ctf(1, p%smpd,200.,2.7,0.1, df,df, 0.178)
! pftcc%ctfmats(:,:,1) = 1.0
pftcc%pfts_ptcls(:,:,1) = pftcc%pfts_ptcls(:,:,1) * pftcc%ctfmats(:,:,1)

! memoizations
call pftcc%memoize_sqsum_ptcl(1)
call pftcc%memoize_ffts
call pftcc%memoize_refs
call pftcc%allocate_ptcls_memoization
call pftcc%memoize_ptcls

call pftcc%gencorrs(1,1, corrs)
call pftcc%gencorrs_dev(1,1, corrs2)
! do irot = 1,pftcc%get_nrots()
!     print *,irot,corrs(irot), corrs2(irot), abs(corrs(irot)-corrs2(irot))
! enddo
corrs = abs(corrs-corrs2)
print *,'gencorr abs err min max: ',minval(corrs),maxval(corrs)

! stop

t = tic()
do irot = 1,500
    call pftcc%gencorrs_dev(1,1, corrs2)
enddo
print *,'gencorrs dev t: ',toc(t)

t = tic()
do irot = 1,500
    call pftcc%gencorrs(1,1, corrs)
enddo
print *,'gencorrs t    : ',toc(t)

! stop

! shift
shift = 4.*[ran3()-0.5, ran3()-0.5]
print *,'shift: ',shift
call pftcc%gencorrs(1,1, shift, corrs)
call pftcc%gencorrs_dev(1,1, shift, corrs2)
corrs = abs(corrs-corrs2)
print *,'gencorr shifted abs err min max: ',minval(corrs),maxval(corrs)

call pftcc%gencorrs(1,1, shift, corrs)
call pftcc%gencorrs_dev(1,1, shift, corrs2)
do irot = 1,pftcc%get_nrots()
    ! print *,irot,corrs(irot), corrs2(irot), abs(corrs(irot)-corrs2(irot))
    corrs(irot)  = pftcc%gencorr_for_rot_8(1,1, real(shift,dp),irot)
    corrs2(irot) = pftcc%gencorr_for_rot_8_dev(1,1, real(shift,dp), irot)
    ! print *,irot,corrs(irot), corrs2(irot), abs(corrs(irot)-corrs2(irot))
    corrs(irot) = abs(corrs(irot)-corrs2(irot))
enddo
print *,'gencorr for_rot_8 abs err min max: ',minval(corrs),maxval(corrs)

! stop

t = tic()
do irot = 1,pftcc%get_nrots()
    call pftcc%gencorrs_dev(1,1, shift,corrs2)
enddo
print *,'gencorrs shifted dev t ',toc(t)

t = tic()
do irot = 1,pftcc%get_nrots()
    call pftcc%gencorrs(1,1, shift,corrs)
enddo
print *,'gencorrs shifted     t ',toc(t)

t = tic()
do irot = 1,pftcc%get_nrots()
    corrs(irot) = pftcc%gencorr_for_rot_8_dev(1,1,real(shift,dp),irot)
enddo
print *,'gencorr_for_rot_8 dev t ',toc(t)

t = tic()
do irot = 1,pftcc%get_nrots()
    corrs2(irot) = pftcc%gencorr_for_rot_8(1,1,real(shift,dp),irot)
enddo
print *,'gencorr_for_rot_8    t ',toc(t)

! stop
f_err     = 0.d0
gradx_err = 0.d0
grady_err = 0.d0
do irot = 1,nrots
    call pftcc%gencorr_grad_for_rot_8(1,1, real(shift,dp),irot,f,grad)
    call pftcc%gencorr_grad_for_rot_8_dev(1,1, real(shift,dp),irot,f_dev,grad_dev)
    f_err     = f_err + abs(f-f_dev)
    gradx_err = gradx_err + abs(grad(1)-grad_dev(1))
    grady_err = grady_err + abs(grad(2)-grad_dev(2))
    ! print *,irot,f,f_dev,grad,grad_dev
enddo
print *,'gencorr_grad_for_rot_8 f abs err avg: ', f_err / real(nrots,dp)
print *,'gencorr_grad_for_rot_8 x abs err avg: ', gradx_err / real(nrots,dp)
print *,'gencorr_grad_for_rot_8 y abs err avg: ', grady_err / real(nrots,dp)

t = tic()
do irot = 1,nrots
    call pftcc%gencorr_grad_for_rot_8(1,1, real(shift,dp),irot,f,grad)
enddo
print *,'gencorr_grad_for_rot_8     t ',toc(t)

t = tic()
do irot = 1,pftcc%get_nrots()
    call pftcc%gencorr_grad_for_rot_8_dev(1,1, real(shift,dp),irot,f_dev,grad_dev)
enddo
print *,'gencorr_grad_for_rot_8 dev t ',toc(t)

! sigma contrib
allocate(sig(pftcc%kfromto(1):pftcc%kfromto(2)),sig2(pftcc%kfromto(1):pftcc%kfromto(2)))

f_err  = 0.d0
do irot = 1,nrots
    call pftcc%gencorr_sigma_contrib(1,1, shift,irot,sig)
    call pftcc%gencorr_sigma_contrib_dev(1,1, shift,irot,sig2)
    f_err = f_err + sum(abs(sig-sig2))
    ! print *,irot,' --- ',sum(abs(sig-sig2))
enddo
print *,'gencorr_sigma_contrib abs err avg: ', f_err / real(nrots,dp)

t = tic()
do irot = 1,pftcc%get_nrots()
    call pftcc%gencorr_sigma_contrib(1,1, shift,irot,sig)
enddo
print *,'gencorr_sigma_contrib     t ',toc(t)

t = tic()
do irot = 1,pftcc%get_nrots()
    call pftcc%gencorr_sigma_contrib_dev(1,1, shift,irot,sig2)
enddo
print *,'gencorr_sigma_contrib dev t ',toc(t)

! vs. not fourier accelerated version
call pftcc%gencorrs(1,1, shift, corrs)
call pftcc%gencorrs_dev(1,1, shift, corrs2)
do irot = 1,pftcc%get_nrots()
    corr         = pftcc%calc_corr_rot_shift(1,1,shift,irot)
    corrs(irot)  = abs(corrs(irot)-corr)
    corrs2(irot) = abs(corrs2(irot)-corr)
    ! print *,irot,corrs(irot), corrs2(irot), pftcc%calc_corr_rot_shift(1,1,shift,irot)
enddo
print *,'gencorr cc     w vs. w/o err min max: ',minval(corrs),maxval(corrs)
print *,'gencorr cc dev w vs. w/o err min max: ',minval(corrs2),maxval(corrs2)

! cleanup
call pftcc%kill_memoized_ptcls
call pftcc%kill_memoized_refs
call pftcc%kill
end program simple_test_fast_corrcalc    