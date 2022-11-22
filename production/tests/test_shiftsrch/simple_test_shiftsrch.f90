program simple_test_shiftsrch
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
type(polarizer)        :: img_copy
type(pftcc_shsrch_grad):: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
logical                :: be_verbose=.false.
real,    parameter     :: SHMAG=1.0
integer, parameter     :: N_PTCLS = 9
real,    allocatable   :: corrs(:), norm_const(:, :)
real                   :: corrmax, corr, cxy(3), lims(2,2)
integer                :: xsh, ysh, xbest, ybest, i, irot
real, allocatable      :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)',advance='no') 'simple_test_shiftsrch stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
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
call p%new(cline)
p%kfromto(1) = 2
p%kfromto(2) = 40
allocate( sigma2_noise(p%kfromto(1):p%kfromto(2), 1:N_PTCLS), source=1. )
call b%build_general_tbox(p, cline)
call pftcc%new(N_PTCLS, [1,N_PTCLS], p%kfromto, .false.)
call pftcc%assign_sigma2_noise(sigma2_noise)
allocate(corrs(pftcc%get_nrots()), norm_const(pftcc%get_nrots(), 2))
call img_copy%init_polarizer(pftcc, p%alpha)
call b%img%read(p%stk, 1)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 1, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 1, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(1, [SHMAG,0.,0.]) ! left
call img_copy%ifft()
call img_copy%write('shifted.mrc', 1)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 2, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 2, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(2, [0.,SHMAG,0.]) ! down
call img_copy%ifft()
call img_copy%write('shifted.mrc', 2)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 3, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 3, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(3, [-SHMAG,0.,0.]) ! right
call img_copy%ifft()
call img_copy%write('shifted.mrc', 3)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 4, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 4, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(4, [0.,SHMAG,0.]) ! up
call img_copy%ifft()
call img_copy%write('shifted.mrc', 4)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 5, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 5, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(5, [SHMAG,SHMAG,0.]) ! left + down
call img_copy%ifft()
call img_copy%write('shifted.mrc', 5)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 6, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 6, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(6, [-SHMAG,-SHMAG,0.]) ! right + up
call img_copy%ifft()
call img_copy%write('shifted.mrc', 6)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 7, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 7, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(7, [-SHMAG,SHMAG,0.]) ! right + down
call img_copy%ifft()
call img_copy%write('shifted.mrc', 7)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 8, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 8, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(8, [SHMAG,-SHMAG,0.]) ! left + up
call img_copy%ifft()
call img_copy%write('shifted.mrc', 8)
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 9, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, 9, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(9, [0.,0.,0.]) ! no shift
call img_copy%ifft()
call img_copy%write('shifted.mrc', 9)
call pftcc%memoize_ffts
lims(1,1) = -6.
lims(1,2) =  6.
lims(2,1) = -6.
lims(2,2) =  6.
call grad_shsrch_obj%new(lims, opt_angle=.true.)
call grad_shsrch_obj%set_indices(8, 8)
irot = 0
cxy  = grad_shsrch_obj%minimize(irot)
print *, cxy(1), cxy(2:3), irot
! do i=1,N_PTCLS
!     corrmax = -1.
!     do xsh=-2,2
!         do ysh=-2,2
!             call pftcc%gencorrs(i, i, real([xsh,ysh]), corrs)
!             corr  = maxval(corrs)

!             print *, 'corr: ', corr, xsh, ysh

!             if( corr > corrmax )then
!                 corrmax = corr
!                 xbest   = xsh
!                 ybest   = ysh
!             endif
!         end do
!     end do
!     print *, xbest, ybest
! end do
end program simple_test_shiftsrch
