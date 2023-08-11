program simple_test_polar2cartesian
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_regularizer,      only: regularizer
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
type(regularizer)      :: reg_obj
type(polarizer)        :: img_copy, img
type(pftcc_shsrch_grad):: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient

logical                :: be_verbose=.false.
real,    parameter     :: SHMAG=1.0
integer, parameter     :: N_PTCLS = 1
real,    allocatable   :: corrs(:)
real                   :: shift(2), corr, cxy(3), lims(2,2), rmsd, df
complex, allocatable   :: cmat(:,:)
complex(dp), allocatable :: pft(:,:), pft_bak(:,:)
integer                :: xsh, ysh, xbest, ybest, i, irot, nrots, box
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
call cline%set('objfun','cc')
call cline%set('ctf','no')
call p%new(cline)
p%kfromto(1) = 1
p%kfromto(2) = 80
allocate( sigma2_noise(p%kfromto(1):p%kfromto(2), 1:N_PTCLS), source=1.0 )
call b%build_general_tbox(p, cline)
call pftcc%new(N_PTCLS, [1,N_PTCLS], p%kfromto)
call reg_obj%new(pftcc)
nrots = pftcc%get_nrots()
allocate(corrs(nrots))
call img_copy%init_polarizer(pftcc, p%alpha)

call b%img%read(p%stk, 1)
call b%img%norm
call b%img%mask(p%msk,'soft',backgr=0.)

! reference
call b%img%write('reference.mrc')
img_copy = b%img
call img_copy%fft()
call img_copy%polarize(pftcc, 1, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call pftcc%polar2cartesian(1,.true.,cmat,box)
call img%new([box,box,1],1.0)
call img%set_cmat(cmat)
call img%shift_phorig()
call img%ifft
call img%write('reference_polar.mrc')

! particle <= scrambled reference
call b%img%rtsq(25.,20.0,10.0, img_copy)
call img_copy%write('particle_shifted_rotated.mrc')
call img_copy%fft()
call img_copy%polarize(pftcc, 1, isptcl=.true., iseven=.true., mask=b%l_resmsk)
call pftcc%memoize_ptcls
call pftcc%polar2cartesian(1,.false.,cmat,box)
call img%new([box,box,1],1.0)
call img%set_cmat(cmat)
call img%shift_phorig()
call img%ifft
call img%write('particle_shifted_rotated_polar.mrc')

call pftcc%gencorrs(1,1,corrs)
irot = maxloc(corrs,dim=1)
lims(1,:) = [-40.,40.]
lims(2,:) = [-40.,40.]
call grad_shsrch_obj%new(lims, opt_angle=.true., coarse_init=.true.)
call grad_shsrch_obj%set_indices(1,1)
cxy  = grad_shsrch_obj%minimize(irot)
print *, cxy(1), cxy(2:3), 360.-pftcc%get_rot(irot)

! stashing ptcl pft
pft_bak = cmplx(pftcc%pfts_ptcls(:,:,1),kind=dp)
allocate(pft(nrots/2,pftcc%nk))

! 0: -shift then irot
pftcc%pfts_ptcls(:,:,1) = cmplx(pft_bak,kind=sp)
call pftcc%shift_ptcl(1,-cxy(2:3))
call reg_obj%rotate_polar(cmplx(pftcc%pfts_ptcls(:,:,1),kind=dp), pft, irot)
pftcc%pfts_ptcls(:,:,1) = cmplx(pft,kind=sp)
call pftcc%polar2cartesian(1,.false.,cmat,box)
call img%new([box,box,1],1.0)
call img%set_cmat(cmat)
call img%shift_phorig()
call img%ifft
call img%write('particle_aligned_0_polar.mrc')

! 1: +shift then irot
pftcc%pfts_ptcls(:,:,1) = cmplx(pft_bak,kind=sp)
call pftcc%shift_ptcl(1,cxy(2:3))
call reg_obj%rotate_polar(cmplx(pftcc%pfts_ptcls(:,:,1),kind=dp), pft, irot)
pftcc%pfts_ptcls(:,:,1) = cmplx(pft,kind=sp)
call pftcc%polar2cartesian(1,.false.,cmat,box)
call img%new([box,box,1],1.0)
call img%set_cmat(cmat)
call img%shift_phorig()
call img%ifft
call img%write('particle_aligned_1_polar.mrc')


! reverse rotation
irot =  nrots+1 - (irot-1)
if( irot > nrots ) irot = irot - nrots

! 2: -shift then 360-irot
pftcc%pfts_ptcls(:,:,1) = cmplx(pft_bak,kind=sp)
call pftcc%shift_ptcl(1,-cxy(2:3))
call reg_obj%rotate_polar(cmplx(pftcc%pfts_ptcls(:,:,1),kind=dp), pft, irot)
pftcc%pfts_ptcls(:,:,1) = cmplx(pft,kind=sp)
call pftcc%polar2cartesian(1,.false.,cmat,box)
call img%new([box,box,1],1.0)
call img%set_cmat(cmat)
call img%shift_phorig()
call img%ifft
call img%write('particle_aligned_2_polar.mrc')

! 3: +shift then 360-irot
pftcc%pfts_ptcls(:,:,1) = cmplx(pft_bak,kind=sp)
call pftcc%shift_ptcl(1,cxy(2:3))
call reg_obj%rotate_polar(cmplx(pftcc%pfts_ptcls(:,:,1),kind=dp), pft, irot)
pftcc%pfts_ptcls(:,:,1) = cmplx(pft,kind=sp)
call pftcc%polar2cartesian(1,.false.,cmat,box)
call img%new([box,box,1],1.0)
call img%set_cmat(cmat)
call img%shift_phorig()
call img%ifft
call img%write('particle_aligned_3_polar.mrc')


call pftcc%kill
end program simple_test_polar2cartesian