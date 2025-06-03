program simple_test_line_sim
include 'simple_lib.f08'
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_cmdline,           only: cmdline
use simple_builder,           only: builder, build_glob
use simple_image,             only: image
use simple_parameters,        only: parameters, params_glob
use simple_polarizer,         only: polarizer
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_pftcc_shsrch_fm,   only: pftcc_shsrch_fm
use simple_commander_volops,  only: reproject_commander
use simple_optimizer,         only: optimizer
use simple_opt_factory,       only: opt_factory
use simple_opt_spec,          only: opt_spec
use simple_image
implicit none
type(cmdline)                 :: cline, cline_projection
type(builder)                 :: b
type(parameters)              :: p
type(polarft_corrcalc)        :: pftcc
type(polarizer)               :: img_copy, polar_img
type(pftcc_shsrch_grad)       :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
type(pftcc_shsrch_fm)         :: fm_correlator
type(reproject_commander)     :: xreproject
character(len=:), allocatable :: cmd
type(image)                   :: ref_img, ptcl_img, img1, img2, img1_copy, img2_copy, img
logical                :: be_verbose=.false.
real,    parameter     :: SHMAG=1.0
integer, parameter     :: N_PTCLS = 10
real,    allocatable   :: corrs(:), norm_const(:, :)
complex, allocatable   :: cmat(:,:)
real                   :: corrmax, corr, cxy(3), lims(2,2), shvec(2), cur_sh(2), grad(2)
real(dp)               :: curval, maxval
integer                :: xsh, ysh, xbest, ybest, i, j, irot, rc, line_irot, cur_irot, nrots, prev_irot
real, allocatable      :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
logical                :: mrc_exists
class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
integer,          parameter :: NDIM=2, NRESTARTS=1
type(opt_factory) :: ofac                           ! the optimization factory object
type(opt_spec)    :: spec                           ! the optimizer specification object
type(image)       :: ccimgs(2)
character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
real              :: lowest_cost, ang, offset(2), cc
integer           :: box
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)',advance='no') 'ERROR! Usage: simple_test_line_sim stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
    write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
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
        write(*, *) 'Projecting 1JYX.mrc...'
        call cline_projection%set('vol1'      , '1JYX.mrc')
        call cline_projection%set('smpd'      , 1.)
        call cline_projection%set('pgrp'      , 'c1')
        call cline_projection%set('mskdiam'   , 180.)
        call cline_projection%set('nspace'    , 6.)
        call cline_projection%set('nthr'      , 16.)
        call xreproject%execute(cline_projection)
        call cline%set('stk'    , 'reprojs.mrcs')
        call cline%set('smpd'   , 1.)
        call cline%set('nthr'   , 16.)
        call cline%set('stk'    , 'reprojs.mrcs')
        call cline%set('mskdiam', 180.)
    endif
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
allocate( sigma2_noise(p%kfromto(1):p%kfromto(2), 1:N_PTCLS), source=1. )
call b%build_general_tbox(p, cline)
call pftcc%new(N_PTCLS, [1,N_PTCLS], p%kfromto)
call pftcc%assign_sigma2_noise(sigma2_noise)
allocate(corrs(pftcc%get_nrots()), norm_const(pftcc%get_nrots(), 2))
call img_copy%new([p%box_crop,p%box_crop,1],p%smpd_crop)
call img_copy%init_polarizer(pftcc, p%alpha)
call b%img%read(p%stk, 1)
call b%img%norm
call b%img%fft
call b%img%clip_inplace([p%box_crop,p%box_crop,1])
call img_copy%polarize(pftcc, b%img, 1, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 1, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(1, [SHMAG,0.,0.]) ! left
call img_copy%polarize(pftcc, b%img, 2, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 2, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(2, [0.,SHMAG,0.]) ! down
call img_copy%polarize(pftcc, b%img, 3, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 3, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(3, [-SHMAG,0.,0.]) ! right
call img_copy%polarize(pftcc, b%img, 4, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 4, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(4, [0.,SHMAG,0.]) ! up
call img_copy%polarize(pftcc, b%img, 5, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 5, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(5, [SHMAG,SHMAG,0.]) ! left + down
call img_copy%polarize(pftcc, b%img, 6, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 6, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(6, [-SHMAG,-SHMAG,0.]) ! right + up
call img_copy%polarize(pftcc, b%img, 7, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 7, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(7, [-SHMAG,SHMAG,0.]) ! right + down
call img_copy%polarize(pftcc, b%img, 8, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 8, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(8, [SHMAG,-SHMAG,0.]) ! left + up
call img_copy%polarize(pftcc, b%img, 9, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 9, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(9, [0.,0.,0.]) ! no shift
call b%img%read(p%stk, 2)
call b%img%norm
call b%img%fft
call b%img%clip_inplace([p%box_crop,p%box_crop,1])
call img_copy%polarize(pftcc, b%img, 10, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call img_copy%polarize(pftcc, b%img, 10, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
call pftcc%shift_ptcl(10, [0.,0.,0.]) ! no shift
call pftcc%set_with_ctf(.false.)
call pftcc%memoize_refs
do i = 1, N_PTCLS
    call pftcc%memoize_sqsum_ptcl(i)
enddo
call pftcc%memoize_ptcls
lims(1,1) = -6.
lims(1,2) =  6.
lims(2,1) = -6.
lims(2,2) =  6.
call grad_shsrch_obj%new(lims, opt_angle=.true.)
call grad_shsrch_obj%set_indices(1, 5)
irot = 1
cxy  = grad_shsrch_obj%minimize(irot)
print *, 'irot = ', irot
print *, cxy(1), cxy(2:3)
!
call pftcc%polar2cartesian(9,.true.,cmat,box)
call ref_img%new([box,box,1],1.0)
call ref_img%set_cmat(cmat)
call ref_img%shift_phorig()
call ref_img%ifft
call ref_img%write('ref_ptcl.mrc',1)
!
call pftcc%rotate_iptcl(9, irot=100, sh=[0.5,-0.5])
call pftcc%polar2cartesian(9,.false.,cmat,box)
call ptcl_img%new([box,box,1],1.0)
call ptcl_img%set_cmat(cmat)
call ptcl_img%shift_phorig()
call ptcl_img%ifft
call ptcl_img%write('ref_ptcl.mrc',2)
!
call pftcc%memoize_sqsum_ptcl(1)
call pftcc%memoize_ptcls
call grad_shsrch_obj%set_indices(9, 10)
irot = 1
cxy  = 0.
cxy  = grad_shsrch_obj%minimize(irot)
print *, 'irot = ', irot
print *, 'sh   = ', cxy(2:3)
if( irot == 0 ) cxy(2:3) = 0.
print *, 'best line sim    (full    inplane/shift ) = ', pftcc%bestline_sim(9, irot, cxy(2:3), 10)
print *, 'best line sim fm (full    inplane/offset) = ', pftcc%bestline_sim_fm(9, 10, irot, cxy(2:3))
print *, 'best line sim mag(full    inplane/offset) = ', pftcc%bestline_sim_mag(9, 10, irot, cxy(2:3))
! debugging, i.e.
call pftcc%polar2cartesian(9,.true.,cmat,box)
call ref_img%new([box,box,1],1.0)
call ref_img%set_cmat(cmat)
call ref_img%shift_phorig()
call ref_img%ifft
call ref_img%write('rotated_ref.mrc',1)
call pftcc%polar2cartesian(10,.false.,cmat,box)
call ptcl_img%new([box,box,1],1.0)
call ptcl_img%set_cmat(cmat)
call ptcl_img%shift_phorig()
call ptcl_img%ifft
call ptcl_img%write('rotated_ref.mrc',2)
! using fm correlator for offsets inplane rotation before using bestline_sim
call fm_correlator%new(p%trs,1.,opt_angle=.false.)
call ref_img%kill
call ref_img%new([p%box_crop,p%box_crop,1],p%smpd_crop)
call ref_img%read(p%stk, 1)
call ref_img%norm
call ref_img%fft
call ref_img%clip_inplace([p%box_crop,p%box_crop,1])
call ptcl_img%kill
call ptcl_img%new([p%box_crop,p%box_crop,1],p%smpd_crop)
call ptcl_img%read(p%stk, 2)
call ptcl_img%norm
call ptcl_img%fft
call ptcl_img%clip_inplace([p%box_crop,p%box_crop,1])
call ccimgs(1)%new( ref_img%get_ldim(), p%smpd_crop, wthreads=.false.)
call ccimgs(2)%new(ptcl_img%get_ldim(), p%smpd_crop, wthreads=.false.)
call fm_correlator%calc_phasecorr(9, 10, ref_img, ptcl_img, ccimgs(1), ccimgs(2), cc, rotang=ang, shift=offset)
print *, 'cc = ', cc
irot = pftcc%get_roind(ang)
print *, 'fm irot   = ', irot
print *, 'fm offset = ', offset
print *, 'best line sim    (full fm inplane/offset) = ', pftcc%bestline_sim(9, irot, offset, 10)
print *, 'best line sim fm                          = ', pftcc%bestline_sim_fm(9, 10)
print *, 'best line sim fm (full fm inplane/offset) = ', pftcc%bestline_sim_fm(9, 10, irot, offset)
print *, 'best line sim mag                         = ', pftcc%bestline_sim_mag(9, 10)
print *, 'best line sim mag(full fm inplane/offset) = ', pftcc%bestline_sim_mag(9, 10, irot, offset)
end program simple_test_line_sim
