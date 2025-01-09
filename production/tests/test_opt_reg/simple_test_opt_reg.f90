program simple_test_opt_reg
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_cmdline,           only: cmdline
use simple_builder,           only: builder
use simple_image,             only: image
use simple_parameters,        only: parameters
use simple_optimizer,         only: optimizer
use simple_opt_factory,       only: opt_factory
use simple_opt_spec,          only: opt_spec
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_strategy2D3D_common
use simple_simulator
use simple_ctf
use simple_ori
use simple_classaverager
use simple_euclid_sigma2
implicit none
type(cmdline)            :: cline
type(builder)            :: b
type(parameters)         :: p
type(polarft_corrcalc)   :: pftcc
type(pftcc_shsrch_grad)  :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
type(ctf)                :: tfun
type(ori)                :: o
type(oris)               :: os
type(ctfparams)          :: ctfparms
type(euclid_sigma2)      :: eucl
type(opt_factory)        :: ofac                           ! the optimization factory object
type(opt_spec)           :: spec                           ! the optimizer specification object
class(optimizer),pointer :: opt_ptr=>null()                ! the generic optimizer object
type(image), allocatable :: match_imgs(:), ptcl_match_imgs(:)
integer,     allocatable :: pinds(:)
real,        allocatable :: truth_sh(:,:), lims(:,:), sigma2_group(:,:,:)
real,        pointer     :: rmat_cavg_even(:,:,:), rmat_tmp(:,:,:)
real,        parameter   :: SHMAG = 3.0, SNR = 0.05, BFAC = 10.
integer,     parameter   :: N_PTCLS = 100, NRESTARTS = 10, SH_ITERS = 5
integer                  :: iptcl, nptcls2update, ithr, ndim, iter, irot, ne, no
real                     :: lowest_cost, cxy(3), lims_init(2,2), lims_sh(2,2)
character(len=8)         :: str_opts                       ! string descriptors for the NOPTS optimizers
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)',advance='no') 'ERROR! required arguments: '
endif
call cline%parse_oldschool
call cline%checkvar('stk',      1)
call cline%checkvar('mskdiam',  2)
call cline%checkvar('smpd',     3)
call cline%checkvar('lp',       4)
call cline%check
call cline%set('oritype','ptcl2D')
if( .not.cline%defined('objfun') ) call cline%set('objfun', 'euclid')
call cline%set('ml_reg', 'no')
call cline%set('ncls',   1.)
call cline%set('kv',     300)
call cline%set('cs',     2.7)
call cline%set('fraca',  0.1)
call cline%set('nptcls',  N_PTCLS)

! general input
call b%init_params_and_build_strategy2D_tbox(cline, p)
call b%spproj%projinfo%new(1,is_ptcl=.false.)
call b%spproj%projinfo%set(1,'projname', 'test')
call b%spproj%os_ptcl2D%kill
call b%spproj%os_ptcl3D%kill
p%fromp          = 1
p%top            = p%nptcls
p%frcs           = trim(FRCS_FILE)
ctfparms%smpd    = p%smpd
ctfparms%kv      = p%kv
ctfparms%cs      = p%cs
if( trim(p%ctf) .eq. 'no' )then
ctfparms%ctfflag = CTFFLAG_NO
endif
ctfparms%fraca   = p%fraca
tfun             = ctf(p%smpd, p%kv, p%cs, p%fraca)

! generate particles
call b%img%read(p%stk, p%iptcl)
call prepimgbatch(N_PTCLS)
call os%new(p%nptcls,is_ptcl=.true.)
if( trim(p%ctf) .eq. 'no' )then
else
    call os%rnd_ctf(p%kv, p%cs, p%fraca, 2.5, 1.5, 0.001)
endif
allocate(truth_sh(p%fromp:p%top,2))
do iptcl = p%fromp,p%top
    call os%set(iptcl,'state',1.)
    call os%set(iptcl,'w',    1.)
    call os%set(iptcl,'class',1.)
    if( mod(iptcl, 2) == 0 )then
        b%imgbatch(iptcl) = b%imgbatch(iptcl-1)
        truth_sh(iptcl,:) = truth_sh(iptcl-1,:)
        call os%set(iptcl,'x', truth_sh(iptcl,1))
        call os%set(iptcl,'y', truth_sh(iptcl,2))
    else
        call b%img%read(p%stk, 1)
        truth_sh(iptcl,:) = [gasdev( 0., SHMAG), gasdev( 0., SHMAG)]
        call os%set(iptcl,'x', truth_sh(iptcl,1))
        call os%set(iptcl,'y', truth_sh(iptcl,2))
        call os%get_ori(iptcl, o)
        call b%img%pad(b%img_pad)
        call b%img_pad%fft
        call b%img_pad%shift2Dserial(truth_sh(iptcl,:) )
        call simimg(b%img_pad, o, tfun, p%ctf, SNR, bfac=BFAC)
        call b%img_pad%clip(b%imgbatch(iptcl))
    endif
    call b%imgbatch(iptcl)%write('particles.mrc',iptcl)
enddo
do iptcl = p%fromp,p%top
    call os%set(iptcl,'x', 0.)
    call os%set(iptcl,'y', 0.)
enddo
call b%spproj%add_single_stk('particles.mrc', ctfparms, os)
call b%spproj_field%partition_eo
call b%spproj_field%sample4update_all([p%fromp,p%top],nptcls2update, pinds, .true.)

! pftcc
call pftcc%new(p%nptcls, [1,p%nptcls], p%kfromto)
call eucl%new('dummy.dat', p%box)
call eucl%allocate_ptcls
allocate(match_imgs(p%ncls),ptcl_match_imgs(nthr_glob))
call pftcc%reallocate_ptcls(p%nptcls, pinds)
do ithr = 1,nthr_glob
    call ptcl_match_imgs(ithr)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
enddo

! references
call restore_read_polarize_cavgs(0)

! particles
!$omp parallel do default(shared) private(iptcl,ithr)&
!$omp schedule(static) proc_bind(close)
do iptcl = 1,p%nptcls
    ithr  = omp_get_thread_num() + 1
    call prepimg4align(iptcl, b%imgbatch(iptcl), ptcl_match_imgs(ithr))
    call b%imgbatch(iptcl)%ifft
    call b%img_crop_polarizer%polarize(pftcc, ptcl_match_imgs(ithr), iptcl, .true., mod(iptcl,2)==0)
    call pftcc%set_eo(iptcl, mod(iptcl,2)==0 )
end do
!$omp end parallel do
call pftcc%create_polar_absctfmats(b%spproj, 'ptcl2D')

! prep for shift search
lims_sh(:,1)    = -p%trs
lims_sh(:,2)    =  p%trs
lims_init(:,1)  = -SHC_INPL_TRSHWDTH
lims_init(:,2)  =  SHC_INPL_TRSHWDTH
call grad_shsrch_obj%new(lims_sh, lims_init=lims_init, maxits=p%maxits_sh, opt_angle=.true., coarse_init=.false.)

! prep for vol gradient descent
ndim      = p%box_crop**2
allocate(lims(ndim,2))
str_opts  = 'lbfgsb'
lims(:,1) = -100.
lims(:,2) =  100.
call spec%specify(str_opts, ndim, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
call spec%set_costfun(costfct)                                      ! set pointer to costfun
call spec%set_gcostfun(gradfct)                                     ! set pointer to gradient of costfun
call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory

! initial sigma2
allocate( sigma2_group(2,1,1:fdim(p%box)-1), source=0. )
do iter = 1, SH_ITERS
    ne = 0
    no = 0
    do iptcl = p%fromp,p%top
        call b%spproj_field%get_ori(iptcl, o)
        call eucl%calc_sigma2(pftcc, iptcl, o, 'class')
        if( o%get_eo() == 0 )then
            ne = ne+1
            sigma2_group(1,1,:) = sigma2_group(1,1,:) + eucl%sigma2_part(:,iptcl)
        else
            no = no+1
            sigma2_group(2,1,:) = sigma2_group(2,1,:) + eucl%sigma2_part(:,iptcl)
        endif
    enddo
    sigma2_group(1,:,:) = sigma2_group(1,:,:) / real(ne)
    sigma2_group(2,:,:) = sigma2_group(2,:,:) / real(no)
    call write_groups_starfile(sigma2_star_from_iter(iter-1), sigma2_group, 1)
    call eucl%read_groups(b%spproj_field)
    do iptcl = p%fromp,p%top
        call pftcc%memoize_sqsum_ptcl(iptcl)
    enddo
    call pftcc%memoize_ptcls

    ! shift search
    do iptcl = p%fromp,p%top
        call grad_shsrch_obj%set_indices(1, iptcl)
        irot = 1 ! zero angle
        cxy  = grad_shsrch_obj%minimize(irot=irot)
        if( iptcl == 3 .or. iptcl == 4 )then
            print *,'iptcl = ',   iptcl
            print *,'irot  ',     irot
            print *,'score ',     cxy(1)
            print *,'cur shift ', b%spproj_field%get_2Dshift(iptcl)
            print *,'opt shift ', cxy(2:3)
            print *,'truth ',     truth_sh(iptcl,:)
        endif
        call b%spproj_field%set_shift(iptcl, cxy(2:3)) !!
    enddo

    call restore_read_polarize_cavgs(iter)

    if( iter < SH_ITERS )then
        call cavgs_even(1)%get_rmat_ptr(rmat_cavg_even)
        spec%x = reshape(rmat_cavg_even(1:p%box_crop,1:p%box_crop,1), (/ndim/))
        call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
        rmat_cavg_even(1:p%box_crop,1:p%box_crop,1) = reshape(spec%x, (/p%box_crop, p%box_crop/))
        p%refs_even = trim(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//'_even'//p%ext
        p%refs_odd  = trim(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//'_odd' //p%ext
        call cavger_write(trim(p%refs_even), 'even')
        call b%clsfrcs%read(FRCS_FILE)
        call cavger_read(trim(p%refs_odd),   'odd' )
        call b%img_crop_polarizer%init_polarizer(pftcc, p%alpha)
        call match_imgs(1)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
        call prep2Dref(cavgs_even(1), match_imgs(1), 1, iseven=.true., center=.false.)
        call b%img_crop_polarizer%polarize(pftcc, match_imgs(1), 1, isptcl=.false., iseven=.true.)
        call prep2Dref(cavgs_odd(1), match_imgs(1), 1, iseven=.false., center=.false.)
        call b%img_crop_polarizer%polarize(pftcc, match_imgs(1), 1, isptcl=.false., iseven=.false.)
        call pftcc%memoize_refs
    endif
enddo

! last one is truth
do iptcl = p%fromp,p%top
    call b%spproj_field%set_shift(iptcl, truth_sh(iptcl,:))
enddo
call restore_read_polarize_cavgs(iter)

call opt_ptr%kill
deallocate(opt_ptr)

contains

    subroutine restore_read_polarize_cavgs( iter )
        integer, intent(in) :: iter
        p%which_iter = iter
        call cavger_kill()
        call cavger_new
        call cavger_transf_oridat( b%spproj )
        call cavger_assemble_sums( .false. )
        call cavger_merge_eos_and_norm
        call cavger_calc_and_write_frcs_and_eoavg(p%frcs, p%which_iter)
        p%refs      = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//p%ext
        p%refs_even = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//'_even'//p%ext
        p%refs_odd  = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//'_odd'//p%ext
        call cavger_write(trim(p%refs),      'merged')
        call cavger_write(trim(p%refs_even), 'even'  )
        call cavger_write(trim(p%refs_odd),  'odd'   )
        call b%clsfrcs%read(FRCS_FILE)
        call cavger_read(trim(p%refs_even), 'even' )
        call cavger_read(trim(p%refs_odd),  'odd' )
        call b%img_crop_polarizer%init_polarizer(pftcc, p%alpha)
        call match_imgs(1)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
        call prep2Dref(cavgs_even(1), match_imgs(1), 1, iseven=.true., center=.false.)
        call b%img_crop_polarizer%polarize(pftcc, match_imgs(1), 1, isptcl=.false., iseven=.true.)
        call prep2Dref(cavgs_odd(1), match_imgs(1), 1, iseven=.false., center=.false.)
        call b%img_crop_polarizer%polarize(pftcc, match_imgs(1), 1, isptcl=.false., iseven=.false.)
        call pftcc%memoize_refs
    end subroutine restore_read_polarize_cavgs

    function costfct( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        real :: x_norm(d), y_norm(d)
        x_norm = x / sqrt(sum(x**2))
        r      = 0.
        do iptcl = p%fromp,p%top
            ithr = omp_get_thread_num() + 1
            call prepimg4align(iptcl, b%imgbatch(iptcl), ptcl_match_imgs(ithr))
            call b%imgbatch(iptcl)%ifft
            call ptcl_match_imgs(ithr)%ifft
            call ptcl_match_imgs(ithr)%get_rmat_ptr(rmat_tmp)
            y_norm = reshape(rmat_tmp(1:p%box_crop,1:p%box_crop,1), (/ndim/))
            y_norm = y_norm / sqrt(sum(y_norm**2))
            r      = r + acos(sum(x_norm * y_norm))
        end do
        r = r / real(p%nptcls)
        print *, r
    end function

    subroutine gradfct( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        real :: abs_x, x_norm(d), y_norm(d), xy
        abs_x  = sqrt(sum(x**2))
        x_norm = x / abs_x
        grad   = 0.
        do iptcl = p%fromp,p%top
            ithr = omp_get_thread_num() + 1
            call prepimg4align(iptcl, b%imgbatch(iptcl), ptcl_match_imgs(ithr))
            call b%imgbatch(iptcl)%ifft
            call ptcl_match_imgs(ithr)%ifft
            call ptcl_match_imgs(ithr)%get_rmat_ptr(rmat_tmp)
            y_norm = reshape(rmat_tmp(1:p%box_crop,1:p%box_crop,1), (/ndim/))
            y_norm = y_norm / sqrt(sum(y_norm**2))
            xy     = sum(x_norm * y_norm)
            grad   = grad - (y_norm * abs_x - xy * x) / abs_x**2 / sqrt(1. - xy**2) / real(p%nptcls)
        end do
    end subroutine

end program simple_test_opt_reg