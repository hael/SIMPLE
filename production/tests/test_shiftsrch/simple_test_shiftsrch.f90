program simple_test_shiftsrch
include 'simple_lib.f08'
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_cmdline,           only: cmdline
use simple_builder,           only: builder
use simple_image,             only: image
use simple_parameters,        only: parameters, params_glob
use simple_polarizer,         only: polarizer
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_commander_volops,  only: reproject_commander
implicit none
type(cmdline)                 :: cline, cline_projection
type(builder)                 :: b
type(parameters)              :: p
type(polarft_corrcalc)        :: pftcc
type(polarizer)               :: img_copy
type(pftcc_shsrch_grad)       :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
type(reproject_commander)     :: xreproject
character(len=:), allocatable :: cmd
logical                :: be_verbose=.false.
real,    parameter     :: SHMAG=1.0
integer, parameter     :: N_PTCLS = 9
real,    allocatable   :: corrs(:), norm_const(:, :)
real                   :: corrmax, corr, cxy(3), lims(2,2), sh(2)
integer                :: xsh, ysh, xbest, ybest, i, irot, rc
real, allocatable      :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
logical                :: mrc_exists
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)',advance='no') 'ERROR! Usage: simple_test_shiftsrch stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
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
p%kfromto(1) = 2
p%kfromto(2) = 40
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
call pftcc%gencorr_sigma_contrib(5,5,[SHMAG,SHMAG],1,sigma2_noise(:,5))
call pftcc%assign_sigma2_noise(sigma2_noise)
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
call pftcc%set_with_ctf(.false.)
call b%img%ifft
call b%img%read(p%stk, 5)
call b%img%norm
call b%img%fft
call b%img%clip_inplace([p%box_crop,p%box_crop,1])
call img_copy%fft
call img_copy%polarize(pftcc, b%img, 1, isptcl=.false., iseven=.true., mask=b%l_resmsk)
call pftcc%memoize_refs
do i = 1, N_PTCLS
    call pftcc%memoize_sqsum_ptcl(i)
enddo
call pftcc%memoize_ptcls
lims(1,1) = -6.
lims(1,2) =  6.
lims(2,1) = -6.
lims(2,2) =  6.
call grad_shsrch_obj%new(lims, opt_angle=.false., multirefs=.true.)
call grad_shsrch_obj%set_indices([5, 5], [.3, .7], 5)
irot = 1
cxy  = grad_shsrch_obj%minimize(irot)
print *, cxy(1), cxy(2:3), irot
params_glob%nstates = 2
do i=5,5
    call pftcc%gencorrs([i, i], [.5, .5], i, corrs)
    print *, 'corr: ', maxval(corrs)
    corrmax = 0.
    do xsh=-2,2
        do ysh=-2,2
            call pftcc%gencorrs([i, i], [.5, .5], i, real([xsh,ysh]), corrs)
            corr  = maxval(corrs)

            print *, 'corr: ', corr, xsh, ysh

            if( corr > corrmax )then
                corrmax = corr
                xbest   = xsh
                ybest   = ysh
            endif
        enddo
    enddo
    print *, xbest, ybest, corrmax
enddo
call pftcc%calc_shift(5, 2, sh)
print *, 'calculated shift = ', sh
end program simple_test_shiftsrch
