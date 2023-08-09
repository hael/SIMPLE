program simple_test_fast_corrcalc
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_cmdline,          only: cmdline
use simple_builder,          only: builder
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_polarizer,        only: polarizer
use simple_pftcc_shsrch_grad,only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
implicit none
type(cmdline)          :: cline
type(builder)          :: b
type(parameters)       :: p
type(polarft_corrcalc) :: pftcc
type(polarizer)        :: img_copy, img
type(pftcc_shsrch_grad):: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
logical                :: be_verbose=.false.
real,    parameter     :: SHMAG=1.0
integer, parameter     :: N_PTCLS = 1
real,    allocatable   :: corrs(:), sig(:)
real                   :: shift(2), corrmax, corr, cxy(3), lims(2,2), rmsd, df
integer                :: i, irot, nrots
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
call cline%checkvar('iptcl',    4)
call cline%check
be_verbose = .false.
if( cline%defined('verbose') )then
    if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
        be_verbose = .true.
    endif
endif
call cline%set('objfun','cc')
call p%new(cline)
p%kfromto(1) = 1
p%kfromto(2) = 60
allocate( sigma2_noise(p%kfromto(1):p%kfromto(2), 1:N_PTCLS), source=1.0 )
call b%build_general_tbox(p, cline)
call pftcc%new(N_PTCLS, [1,N_PTCLS], p%kfromto)
call pftcc%assign_sigma2_noise(sigma2_noise)
nrots = pftcc%get_nrots()
allocate(corrs(nrots))
call img_copy%init_polarizer(pftcc, p%alpha)
call b%img%read(p%stk, p%iptcl)
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
call pftcc%calc_polar_ctf(1, p%smpd,200.,2.7,0.1, df,df+0.1, 0.178)
! pftcc%ctfmats(:,:,1) = 1.0
pftcc%pfts_ptcls(:,:,1) = pftcc%pfts_ptcls(:,:,1) * pftcc%ctfmats(:,:,1)

! memoizations
call pftcc%memoize_sqsum_ptcl(1)
call pftcc%memoize_refs
call pftcc%memoize_ptcls


! shift
shift = 0.
! shift = 4.*[ran3()-0.5, ran3()-0.5]
print *,'shift: ',shift

! vs. not fourier accelerated version
call pftcc%gencorrs(1,1, shift, corrs)
do irot = 1,pftcc%get_nrots()
    corr         = pftcc%calc_corr_rot_shift(1,1,shift,irot)
    corrs(irot)  = abs(corrs(irot)-corr)
    ! print *,irot,corrs(irot), pftcc%calc_corr_rot_shift(1,1,shift,irot)
enddo
print *,'gencorr cc     w vs. w/o err min max: ',minval(corrs),maxval(corrs)

! cleanup
call pftcc%kill
end program simple_test_fast_corrcalc    