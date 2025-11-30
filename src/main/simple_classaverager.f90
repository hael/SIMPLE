module simple_classaverager
include 'simple_lib.f08'
!$ use omp_lib
use simple_builder,           only: build_glob
use simple_ctf,               only: ctf
use simple_discrete_stack_io, only: dstack_io
use simple_euclid_sigma2,     only: euclid_sigma2, eucl_sigma2_glob
use simple_image,             only: image, image_ptr
use simple_parameters,        only: params_glob
use simple_stack_io,          only: stack_io
implicit none

public :: cavger_new, cavger_transf_oridat, cavger_gen2Dclassdoc, cavger_assemble_sums,&
cavger_merge_eos_and_norm, cavger_calc_and_write_frcs_and_eoavg, cavger_write, cavger_read, cavger_read_all,&
cavger_readwrite_partial_sums, cavger_assemble_sums_from_parts, cavger_kill, cavgs_even, cavgs_odd, cavgs_merged,&
cavger_read_euclid_sigma2, transform_ptcls
private
#include "simple_local_flags.inc"

type ptcl_record
    type(ctf) :: tfun              !< transfer function
    real      :: pw        = 0.0   !< particle weight
    real      :: dfx       = 0.0   !< defocus in x (microns)
    real      :: dfy       = 0.0   !< defocus in y (microns)
    real      :: angast    = 0.0   !< angle of astigmatism (in degrees)
    real      :: phshift   = 0.0   !< additional phase shift from the Volta
    real      :: e3                !< in-plane rotations
    real      :: shift(2)          !< rotational origin shift
    integer   :: pind       = 0    !< particle index
    integer   :: eo         = -1   !< even is 0, odd is 1, default is -1
    integer   :: class             !< class assignment
    integer   :: ind_in_stk = 0    !< index in stack
end type ptcl_record

integer                          :: ctfflag                   !< ctf flag <yes=1|no=0|flip=2>
integer                          :: istart      = 0, iend = 0 !< particle index range
integer                          :: partsz      = 0           !< size of partition
integer                          :: ncls        = 0           !< # classes
integer                          :: ldim(3)        = [0,0,0] !< logical dimension of image
integer                          :: ldim_crop(3)   = [0,0,0] !< logical dimension of cropped image
integer                          :: ldim_pd(3)     = [0,0,0] !< logical dimension of image, padded
integer                          :: ldim_croppd(3) = [0,0,0] !< logical dimension of cropped image, padded
real                             :: smpd       = 0.          !< sampling distance
real                             :: smpd_crop  = 0.          !< cropped sampling distance
type(ptcl_record),   allocatable :: precs(:)                 !< particle records
type(image), target, allocatable :: cavgs_even(:)            !< class averages
type(image), target, allocatable :: cavgs_odd(:)             !< -"-
type(image),         allocatable :: cavgs_even_wfilt(:)      !< class averages wiener filtered
type(image),         allocatable :: cavgs_odd_wfilt(:)       !< -"-
type(image), target, allocatable :: cavgs_merged(:)          !< -"-
type(image),         allocatable :: cavgs_merged_wfilt(:)    !< -"-
type(image),         allocatable :: ctfsqsums_even(:)        !< CTF**2 sums for Wiener normalisation
type(image),         allocatable :: ctfsqsums_odd(:)         !< -"-
type(image),         allocatable :: ctfsqsums_even_wfilt(:)  !< CTF**2 sums for Wiener normalisation
type(image),         allocatable :: ctfsqsums_odd_wfilt(:)   !< -"-
type(image),         allocatable :: ctfsqsums_merged(:)      !< -"-
type(image),         allocatable :: ctfsqsums_merged_wfilt(:)!< -"-
type(image),         allocatable :: cavgs_even_part(:)           !< partial class averages
type(image),         allocatable :: cavgs_odd_part(:)            !< -"-
type(image),         allocatable :: cavgs_even_wfilt_bak(:)      !< partial class averages wiener filtered
type(image),         allocatable :: cavgs_odd_wfilt_bak(:)       !< -"-
type(image),         allocatable :: ctfsqsums_even_bak(:)        !< CTF**2 sums for Wiener normalisation
type(image),         allocatable :: ctfsqsums_odd_bak(:)         !< -"-
type(image),         allocatable :: ctfsqsums_even_wfilt_bak(:)  !< CTF**2 sums for Wiener normalisation
type(image),         allocatable :: ctfsqsums_odd_wfilt_bak(:)   !< -"-
type(euclid_sigma2)              :: eucl_sigma
integer,             allocatable :: prev_eo_pops(:,:)
logical,             allocatable :: pptcl_mask(:)
logical                          :: l_stream           = .false.  !< flag for cluster2D_stream
logical                          :: phaseplate         = .false.  !< Volta phaseplate images or not
logical                          :: l_ml_reg           = .false.  !< Maximum-Likelihood regularization
logical                          :: l_alloc_read_cavgs = .true.   !< whether to allocate sums and read partial sums

integer(timer_int_kind) :: t_class_loop,t_batch_loop, t_gridding, t_init, t_tot
real(timer_int_kind)    :: rt_class_loop,rt_batch_loop, rt_gridding, rt_init, rt_tot
type(string)            :: benchfname

contains

    subroutine cavger_new( pinds, alloccavgs )
        integer, optional, intent(in) :: pinds(:)
        logical, optional, intent(in) :: alloccavgs
        l_alloc_read_cavgs = .true.
        if( present(alloccavgs) ) l_alloc_read_cavgs = alloccavgs
        call cavger_kill(dealloccavgs=l_alloc_read_cavgs)
        allocate(pptcl_mask(params_glob%fromp:params_glob%top), source=.true.)
        if( present(pinds) )then
            pptcl_mask           = .false.
            pptcl_mask(pinds(:)) = .true.
        endif
        ncls          = params_glob%ncls
        ! work out range and partsz
        if( params_glob%l_distr_exec )then
            istart    = params_glob%fromp
            iend      = params_glob%top
        else
            istart    = 1
            iend      = params_glob%nptcls
        endif
        partsz        = size(pptcl_mask)
        ! CTF logics
        ctfflag       = build_glob%spproj%get_ctfflag_type('ptcl2D',iptcl=params_glob%fromp)
        ! set CTF-related flag
        phaseplate    = build_glob%spproj%has_phaseplate('ptcl2D')
        l_stream      = (trim(params_glob%stream) .eq. 'yes') .and. params_glob%l_wiener_part
        ! smpd
        smpd          = params_glob%smpd
        smpd_crop     = params_glob%smpd_crop
        ! set ldims
        ldim          = [params_glob%box,  params_glob%box,  1]
        ldim_crop     = [params_glob%box_crop,  params_glob%box_crop,  1]
        ldim_croppd   = [params_glob%box_croppd,params_glob%box_croppd,1]
        ldim_pd       = [params_glob%boxpd,params_glob%boxpd,1]
        ! ML-regularization
        l_ml_reg      = params_glob%l_ml_reg
        ! build arrays
        allocate(precs(partsz), prev_eo_pops(ncls,2))
        prev_eo_pops = 0
        call alloc_cavgs_sums
    end subroutine cavger_new

    ! setters/getters

    !>  \brief  transfers orientation data to the instance
    subroutine cavger_transf_oridat( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        class(oris), pointer :: spproj_field
        type(ctfparams)      :: ctfvars(nthr_glob)
        integer              :: i, icls, cnt, iptcl, ithr, stkind
        call spproj%ptr2oritype(params_glob%oritype, spproj_field)
        ! build index map
        cnt = 0
        do iptcl=istart,iend
            cnt = cnt + 1
            ! exclusion
            precs(cnt)%pind = 0
            if( spproj_field%get_state(iptcl) == 0  ) cycle
            if( spproj_field%get(iptcl,'w') < SMALL ) cycle
            if( .not.pptcl_mask(iptcl)              ) cycle
            precs(cnt)%pind = iptcl
        enddo
        ! fetch data from project
        !$omp parallel do default(shared) private(cnt,iptcl,ithr,stkind) schedule(static) proc_bind(close)
        do cnt = 1,partsz
            iptcl              = precs(cnt)%pind
            if( iptcl == 0 ) cycle
            ithr               = omp_get_thread_num() + 1
            precs(cnt)%eo      = spproj_field%get_eo(iptcl)
            precs(cnt)%pw      = spproj_field%get(iptcl,'w')
            ctfvars(ithr)      = spproj%get_ctfparams(params_glob%oritype,iptcl)
            precs(cnt)%tfun    = ctf(params_glob%smpd_crop, ctfvars(ithr)%kv, ctfvars(ithr)%cs, ctfvars(ithr)%fraca)
            precs(cnt)%dfx     = ctfvars(ithr)%dfx
            precs(cnt)%dfy     = ctfvars(ithr)%dfy
            precs(cnt)%angast  = ctfvars(ithr)%angast
            precs(cnt)%phshift = 0.
            if( phaseplate ) precs(cnt)%phshift = ctfvars(ithr)%phshift
            precs(cnt)%class   = spproj_field%get_class(iptcl)
            precs(cnt)%e3      = spproj_field%e3get(iptcl)
            precs(cnt)%shift   = spproj_field%get_2Dshift(iptcl)
            call spproj%map_ptcl_ind2stk_ind(params_glob%oritype, iptcl, stkind, precs(cnt)%ind_in_stk)
        end do
        !$omp end parallel do
        prev_eo_pops = 0
        if( l_stream .and. spproj%os_cls2D%get_noris() == ncls )then
            do i = 1,ncls
                icls = spproj%os_cls2D%get_class(i)
                if( .not.spproj%os_cls2D%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(icls,1) = spproj%os_cls2D%get_int(i,'prev_pop_even')
                prev_eo_pops(icls,2) = spproj%os_cls2D%get_int(i,'prev_pop_odd')
            enddo
        endif
    end subroutine cavger_transf_oridat

    subroutine cavger_read_euclid_sigma2
        type(string) :: fname
        if( l_ml_reg )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            call eucl_sigma%read_part(  build_glob%spproj_field)
            call eucl_sigma%read_groups(build_glob%spproj_field)
        end if
    end subroutine cavger_read_euclid_sigma2

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population, average correlation and weight
    subroutine cavger_gen2Dclassdoc( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), target, intent(inout) :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(params_glob%ncls)
        real(dp) :: corrs(params_glob%ncls), ws(params_glob%ncls)
        real     :: frc05, frc0143, rstate, w
        integer  :: i, iptcl, icls, pop, nptcls
        select case(trim(params_glob%oritype))
            case('ptcl2D')
                ptcl_field => spproj%os_ptcl2D
                cls_field  => spproj%os_cls2D
            case('ptcl3D')
                ptcl_field => spproj%os_ptcl3D
                cls_field  => spproj%os_cls3D
            case DEFAULT
                THROW_HARD('Unsupported ORITYPE: '//trim(params_glob%oritype))
        end select
        nptcls = ptcl_field%get_noris()
        pops   = 0
        corrs  = 0.d0
        ws     = 0.d0
        !$omp parallel do default(shared) private(iptcl,rstate,icls,w) schedule(static)&
        !$omp proc_bind(close) reduction(+:pops,corrs,ws)
        do iptcl=1,nptcls
            rstate = ptcl_field%get(iptcl,'state')
            if( rstate < 0.5 ) cycle
            w = ptcl_field%get(iptcl,'w')
            if( w < SMALL ) cycle
            icls = ptcl_field%get_class(iptcl)
            if( icls<1 .or. icls>params_glob%ncls )cycle
            pops(icls)  = pops(icls)  + 1
            corrs(icls) = corrs(icls) + real(ptcl_field%get(iptcl,'corr'),dp)
            ws(icls)    = ws(icls)    + real(w,dp)
        enddo
        !$omp end parallel do
        if( l_stream  .and. cls_field%get_noris()==ncls .and. params_glob%update_frac<.99 )then
            do i = 1,ncls
                icls = cls_field%get_class(i)
                if( .not.cls_field%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(icls,1) = cls_field%get_int(i,'prev_pop_even')
                prev_eo_pops(icls,2) = cls_field%get_int(i,'prev_pop_odd')
                pop = sum(prev_eo_pops(icls,:))
                if( pop == 0 ) cycle
                corrs(icls) = corrs(icls) + real(pop) * cls_field%get(i,'corr')
                ws(icls)    = ws(icls)    + real(pop) * cls_field%get(i,'w')
                pops(icls)  = pops(icls)  +      pop
            enddo
        endif
        where(pops>1)
            corrs = corrs / real(pops)
            ws    = ws / real(pops)
        elsewhere
            corrs = -1.
            ws    = 0.
        end where
        call cls_field%new(params_glob%ncls, is_ptcl=.false.)
        do icls=1,params_glob%ncls
            pop = pops(icls)
            call build_glob%clsfrcs%estimate_res(icls, frc05, frc0143)
            call cls_field%set(icls, 'class',     icls)
            call cls_field%set(icls, 'pop',       pop)
            call cls_field%set(icls, 'res',       frc0143)
            call cls_field%set(icls, 'corr',      corrs(icls))
            call cls_field%set(icls, 'w',         ws(icls))
            if( pop > 0 )then
                call cls_field%set(icls, 'state', 1.0) ! needs to be default val if no selection has been done
            else
                call cls_field%set(icls, 'state', 0.0) ! exclusion
            endif
        end do
        call score_classes(cls_field)
    end subroutine cavger_gen2Dclassdoc

    !>  \brief  is for allocation of the sums array to size used for reading
    subroutine alloc_cavgs_sums
        integer :: icls
        if( l_alloc_read_cavgs )then
            allocate(cavgs_even_part(ncls),cavgs_odd_part(ncls))
            if( l_stream ) allocate(cavgs_even_wfilt_bak(ncls),cavgs_odd_wfilt_bak(ncls))
            allocate(cavgs_even(ncls), cavgs_odd(ncls), cavgs_merged(ncls),&
            &ctfsqsums_even(ncls),ctfsqsums_odd(ncls), ctfsqsums_merged(ncls))
            if( l_stream )then
                allocate(cavgs_even_wfilt(ncls), cavgs_odd_wfilt(ncls), ctfsqsums_merged_wfilt(ncls),&
                &ctfsqsums_even_wfilt(ncls), ctfsqsums_odd_wfilt(ncls), cavgs_merged_wfilt(ncls))
            endif
            !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
            do icls=1,ncls
                call cavgs_even(icls)%new(  ldim_croppd, smpd_crop,wthreads=.false.)
                call cavgs_odd(icls)%new(   ldim_croppd, smpd_crop,wthreads=.false.)
                call cavgs_merged(icls)%new(ldim_croppd, smpd_crop,wthreads=.false.)
                call ctfsqsums_even(icls)%new(  ldim_croppd, smpd_crop,wthreads=.false.)
                call ctfsqsums_odd(icls)%new(   ldim_croppd, smpd_crop,wthreads=.false.)
                call ctfsqsums_merged(icls)%new(ldim_croppd, smpd_crop,wthreads=.false.)
                if( l_stream )then
                    call cavgs_merged_wfilt(icls)%new(ldim_croppd, smpd_crop,wthreads=.false.)
                    call cavgs_even_wfilt(icls)%new(  ldim_croppd, smpd_crop,wthreads=.false.)
                    call cavgs_odd_wfilt(icls)%new(   ldim_croppd, smpd_crop,wthreads=.false.)
                    call ctfsqsums_even_wfilt(icls)%new(  ldim_croppd, smpd_crop,wthreads=.false.)
                    call ctfsqsums_odd_wfilt(icls)%new(   ldim_croppd, smpd_crop,wthreads=.false.)
                    call ctfsqsums_merged_wfilt(icls)%new(ldim_croppd, smpd_crop,wthreads=.false.)
                endif
            end do
            !$omp end parallel do
        endif
    end subroutine alloc_cavgs_sums

    !>  \brief  is for initialization of the sums to size used in restoration
    subroutine init_cavgs_sums
        integer :: icls
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)%new(    ldim_croppd, smpd_crop,wthreads=.false.)
            call cavgs_odd(icls)%new(     ldim_croppd, smpd_crop,wthreads=.false.)
            call cavgs_merged(icls)%new(  ldim_croppd, smpd_crop,wthreads=.false.)
            call cavgs_even(icls)%zero_and_flag_ft
            call cavgs_odd(icls)%zero_and_flag_ft
            call cavgs_merged(icls)%zero_and_flag_ft
            call ctfsqsums_even(icls)%zero_and_flag_ft
            call ctfsqsums_odd(icls)%zero_and_flag_ft
            call ctfsqsums_merged(icls)%zero_and_flag_ft
            if( l_stream )then
                call cavgs_merged_wfilt(icls)%new(ldim_croppd,smpd_crop,wthreads=.false.)
                call cavgs_even_wfilt(icls)%zero_and_flag_ft
                call cavgs_odd_wfilt(icls)%zero_and_flag_ft
                call ctfsqsums_even_wfilt(icls)%zero_and_flag_ft
                call ctfsqsums_odd_wfilt(icls)%zero_and_flag_ft
                call ctfsqsums_merged_wfilt(icls)%zero_and_flag_ft
            endif
        end do
        !$omp end parallel do
    end subroutine init_cavgs_sums

    !>  \brief  is for calculating class population
    integer function class_pop( class )
        integer, intent(in) :: class
        class_pop = sum(eo_class_pop(class))
    end function class_pop

    !>  \brief  is for calculating even/odd class population
    function eo_class_pop( class ) result( pops )
        integer, intent(in) :: class
        integer :: pops(2), iprec
        pops = 0
        do iprec=1,partsz
            if( precs(iprec)%pind > 0 .and. precs(iprec)%class .eq. class )then
                if( precs(iprec)%eo == 1 )then
                    pops(2) = pops(2) + 1
                else
                    pops(1) = pops(1) + 1
                endif
            endif
        end do
    end function eo_class_pop

    ! calculators

    !>  \brief  is for assembling the sums in distributed/non-distributed mode
    !           using gridding interpolation in Fourier space
    subroutine cavger_assemble_sums( do_frac_update )
        logical, intent(in)        :: do_frac_update
        integer, parameter         :: READBUFFSZ = 1024
        complex, parameter         :: zero = cmplx(0.,0.)
        type(kbinterpol)           :: kbwin
        type(dstack_io)            :: dstkio_r
        type(image_ptr)            :: pcmat(nthr_glob), prhomat(nthr_glob)
        type(image), allocatable   :: cgrid_imgs(:), read_imgs(:), cgrid_imgs_crop(:)
        complex,     allocatable   :: cmats(:,:,:)
        real,        allocatable   :: rhos(:,:,:), tvals(:,:,:)
        type(string) :: stk_fname
        complex :: fcompl, fcompll
        real    :: loc(2), mat(2,2), dist(2), add_phshift, tval, kw, maxspafreqsq, reg_scale, crop_scale
        integer :: batch_iprecs(READBUFFSZ), fdims_crop(3), logi_lims_crop(3,2)
        integer :: phys(2), win_corner(2), cyc_lims_cropR(2,2),cyc_lims_crop(3,2), sigma2_kfromto(2)
        integer :: iprec, i, sh, nyq_crop, ind_in_stk, iprec_glob, nptcls_eff, radfirstpeak
        integer :: wdim, h, k, l, m, ll, mm, icls, iptcl, interp_shlim, batchind
        integer :: first_stkind, fromp, top, istk, nptcls_in_stk, nstks, last_stkind
        integer :: ibatch, nbatches, istart, iend, ithr, nptcls_in_batch, first_pind, last_pind
        if( .not. params_glob%l_distr_exec ) write(logfhandle,'(a)') '>>> ASSEMBLING CLASS SUMS'
        ! init cavgs
        if( l_alloc_read_cavgs )then
            call init_cavgs_sums
            if( do_frac_update )then
                call cavger_readwrite_partial_sums('read')
                call cavger_apply_weights(1.0-params_glob%update_frac)
            endif
        else
            if( do_frac_update )then
                call cavger_apply_weights(1.0-params_glob%update_frac)
            else
                call init_cavgs_sums
            endif
        endif
        kbwin  = kbinterpol(KBWINSZ, params_glob%alpha)
        wdim   = kbwin%get_wdim()
        ! Number stacks
        first_pind = params_glob%fromp
        call build_glob%spproj%map_ptcl_ind2stk_ind(params_glob%oritype, first_pind, first_stkind, ind_in_stk)
        last_pind = 0
        do i = partsz,1,-1
            if( precs(i)%pind > 0 )then
                last_pind = precs(i)%pind
                exit
            endif
        enddo
        call build_glob%spproj%map_ptcl_ind2stk_ind(params_glob%oritype, last_pind, last_stkind,  ind_in_stk)
        nstks     = last_stkind - first_stkind + 1
        reg_scale = 1.0
        if( l_ml_reg )then
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
            reg_scale         = real(ldim(1)) / real(ldim_pd(1))
        endif
        crop_scale = real(ldim_croppd(1)) / real(ldim_pd(1))
        ! Objects allocations
        allocate(read_imgs(READBUFFSZ), cgrid_imgs(READBUFFSZ), cgrid_imgs_crop(READBUFFSZ))
        !$omp parallel do schedule(static) default(shared) proc_bind(close) private(i)
        do i = 1,READBUFFSZ
            call read_imgs(i)%new(ldim, params_glob%smpd, wthreads=.false.)
            call cgrid_imgs(i)%new(ldim_pd, params_glob%smpd, wthreads=.false.)
            call cgrid_imgs_crop(i)%new(ldim_croppd, params_glob%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
        logi_lims_crop = cgrid_imgs_crop(1)%loop_lims(2)
        cyc_lims_crop  = cgrid_imgs_crop(1)%loop_lims(3)
        nyq_crop       = cgrid_imgs_crop(1)%get_lfny(1)
        fdims_crop     = cgrid_imgs_crop(1)%get_array_shape()
        cyc_lims_cropR(:,1) = cyc_lims_crop(1,:)
        cyc_lims_cropR(:,2) = cyc_lims_crop(2,:)
        allocate(tvals(fdims_crop(1),fdims_crop(2),READBUFFSZ),cmats(fdims_crop(1),fdims_crop(2),READBUFFSZ),&
        &rhos(fdims_crop(1),fdims_crop(2),READBUFFSZ))
        interp_shlim = nyq_crop + 1
        call dstkio_r%new(smpd, ldim(1))
        ! Main loop
        iprec_glob = 0 ! global record index
        do istk = first_stkind,last_stkind
            ! Particles range in stack
            fromp = build_glob%spproj%os_stk%get_fromp(istk)
            top   = build_glob%spproj%os_stk%get_top(istk)
            nptcls_in_stk = top - fromp + 1 ! # of particles in stack
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, max(params_glob%fromp,fromp), stk_fname, ind_in_stk)
            ! batch loop
            nbatches = ceiling(real(nptcls_in_stk)/real(READBUFFSZ)) ! will be 1 most of the tme
            do ibatch = 1,nbatches
                batch_iprecs = 0                                     ! records in batch, if zero skip
                istart = (ibatch - 1)              * READBUFFSZ + 1  ! first index in current batch, will be 1      most of the time
                iend   = min(nptcls_in_stk, istart + READBUFFSZ - 1) ! last  index in current batch, will be nptcls_in_stk most of the time
                nptcls_in_batch = iend-istart+1
                batchind   = 0
                nptcls_eff = 0                                       ! # particles to process in batch
                do i = istart,iend
                    iptcl    = fromp + i - 1                         ! global particle index
                    batchind = batchind + 1                          ! index in batch
                    if( iptcl < params_glob%fromp ) cycle            ! taking care of limits
                    if( iptcl > params_glob%top )   cycle
                    iprec_glob = iprec_glob + 1                      ! global particle record
                    batch_iprecs(batchind) = iprec_glob              ! particle record in batch
                    if( precs(iprec_glob)%pind == 0 ) cycle
                    nptcls_eff = nptcls_eff + 1
                    call dstkio_r%read(stk_fname, precs(iprec_glob)%ind_in_stk, read_imgs(batchind)) ! read
                enddo
                if( nptcls_eff == 0 ) cycle
                ! Interpolation loop
                !$omp parallel default(shared) proc_bind(close)&
                !$omp private(i,ithr,icls,iprec,win_corner,add_phshift,mat,h,k,l,m,ll,mm,dist,loc,sh,phys,kw,tval,fcompl,fcompll)
                !$omp do schedule(static)
                do i = 1,nptcls_in_batch
                    iprec = batch_iprecs(i)
                    if( iprec == 0 ) cycle
                    if( precs(iprec)%pind == 0 ) cycle
                    cmats(:,:,i) = zero
                    rhos(:,:,i)  = 0.0
                    tvals(:,:,i) = 0.0
                    ! normalize & pad & FFT
                    call read_imgs(i)%norm_noise_pad_fft(build_glob%lmsk, cgrid_imgs(i))
                    ! Fourier cropping
                    call cgrid_imgs(i)%clip(cgrid_imgs_crop(i))
                    ! shift
                    call cgrid_imgs_crop(i)%shift2Dserial(-precs(iprec)%shift*crop_scale)
                    ! apply CTF to image, stores CTF values
                    add_phshift = 0.
                    if( phaseplate ) add_phshift = precs(iprec)%phshift
                    if( l_stream )then
                        call precs(iprec)%tfun%eval_and_apply(cgrid_imgs_crop(i), ctfflag, logi_lims_crop, fdims_crop(1:2), tvals(:,:,i), &
                        & precs(iprec)%dfx, precs(iprec)%dfy, precs(iprec)%angast, add_phshift, .false.)
                    else
                        call precs(iprec)%tfun%eval_and_apply(cgrid_imgs_crop(i), ctfflag, logi_lims_crop, fdims_crop(1:2), tvals(:,:,i), &
                        & precs(iprec)%dfx, precs(iprec)%dfy, precs(iprec)%angast, add_phshift, .not.params_glob%l_wiener_part)
                    endif
                    ! ML-regularization
                    if( l_ml_reg )then
                        do h=logi_lims_crop(1,1),logi_lims_crop(1,2)
                            do k=logi_lims_crop(2,1),logi_lims_crop(2,2)
                                sh = nint(hyp(h,k))                             ! shell in padded image
                                if( sh > interp_shlim )cycle
                                sh = nint(reg_scale*sqrt(real(h*h)+real(k*k)))  ! shell at original scale
                                if( sh > sigma2_kfromto(2) ) cycle
                                phys = cgrid_imgs_crop(i)%comp_addr_phys(h,k)
                                if( sh >= sigma2_kfromto(1) )then
                                    call cgrid_imgs_crop(i)%mul_cmat_at(  phys(1),phys(2),1, 1./eucl_sigma2_glob%sigma2_noise(sh,precs(iprec)%pind))
                                    tvals(phys(1),phys(2),i) = tvals(phys(1),phys(2),i) / sqrt(eucl_sigma2_glob%sigma2_noise(sh,precs(iprec)%pind))
                                else
                                    call cgrid_imgs_crop(i)%mul_cmat_at(  phys(1),phys(2),1, 1./eucl_sigma2_glob%sigma2_noise(1,precs(iprec)%pind))
                                    tvals(phys(1),phys(2),i) = tvals(phys(1),phys(2),i) / sqrt(eucl_sigma2_glob%sigma2_noise(1,precs(iprec)%pind))
                                endif
                            enddo
                        enddo
                    endif
                    ! Rotation matrix
                    call rotmat2d(-precs(iprec)%e3, mat)
                    ! Interpolation
                    do h = logi_lims_crop(1,1),logi_lims_crop(1,2)
                        do k = logi_lims_crop(2,1),logi_lims_crop(2,2)
                            sh = nint(hyp(h,k))
                            if( sh > interp_shlim )cycle
                            ! Rotation
                            loc        = matmul(real([h,k]),mat)
                            win_corner = floor(loc) ! bottom left corner
                            dist       = loc - real(win_corner)
                            ! Bi-linear interpolation
                            l     = cyci_1d(cyc_lims_cropR(:,1), win_corner(1))
                            ll    = cyci_1d(cyc_lims_cropR(:,1), win_corner(1)+1)
                            m     = cyci_1d(cyc_lims_cropR(:,2), win_corner(2))
                            mm    = cyci_1d(cyc_lims_cropR(:,2), win_corner(2)+1)
                            ! l, bottom left corner
                            phys   = cgrid_imgs_crop(i)%comp_addr_phys(l,m)
                            kw     = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
                            fcompl = kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                            tval   = kw * tvals(phys(1),phys(2),i)
                            ! l, bottom right corner
                            phys   = cgrid_imgs_crop(i)%comp_addr_phys(l,mm)
                            kw     = (1.-dist(1))*dist(2)
                            fcompl = fcompl + kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                            tval   = tval   + kw * tvals(phys(1),phys(2),i)
                            if( l < 0 ) fcompl = conjg(fcompl) ! conjugation when required!
                            ! ll, upper left corner
                            phys    = cgrid_imgs_crop(i)%comp_addr_phys(ll,m)
                            kw      = dist(1)*(1.-dist(2))
                            fcompll =         kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                            tval    = tval  + kw * tvals(phys(1),phys(2),i)
                            ! ll, upper right corner
                            phys    = cgrid_imgs_crop(i)%comp_addr_phys(ll,mm)
                            kw      = dist(1)*dist(2)
                            fcompll = fcompll + kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                            tval    = tval    + kw * tvals(phys(1),phys(2),i)
                            if( ll < 0 ) fcompll = conjg(fcompll) ! conjugation when required!
                            ! update with interpolated values
                            phys = cgrid_imgs_crop(i)%comp_addr_phys(h,k)
                            cmats(phys(1),phys(2),i) = fcompl + fcompll
                            rhos(phys(1),phys(2),i)  = tval*tval
                        end do
                    end do
                enddo
                !$omp end do
                ! Sum over classes
                !$omp do schedule(static)
                do icls = 1,ncls
                    ithr = omp_get_thread_num() + 1
                    do i = 1,nptcls_in_batch
                        iprec = batch_iprecs(i)
                        if( iprec == 0 ) cycle
                        if( precs(iprec)%pind == 0 ) cycle
                        if( precs(iprec)%class == icls )then
                            select case(precs(iprec)%eo)
                                case(0,-1)
                                    call cavgs_even(icls)%get_cmat_ptr(pcmat(ithr)%cmat)
                                    call ctfsqsums_even(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
                                case(1)
                                    call cavgs_odd(icls)%get_cmat_ptr(pcmat(ithr)%cmat)
                                    call ctfsqsums_odd(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
                            end select
                            pcmat(ithr)%cmat(:,:,1)   = pcmat(ithr)%cmat(:,:,1)   + precs(iprec)%pw * cmats(:,:,i)
                            prhomat(ithr)%cmat(:,:,1) = prhomat(ithr)%cmat(:,:,1) + precs(iprec)%pw * cmplx(rhos(:,:,i),0.0)
                        endif
                    enddo
                enddo
                !$omp end do
                !$omp end parallel
                if( l_stream .and. ctfflag /= CTFFLAG_NO )then
                    !$omp parallel default(shared) proc_bind(close)&
                    !$omp private(i,icls,iprec,win_corner,add_phshift,mat,ithr,h,k,l,m,ll,mm,dist,loc,sh,phys,kw,tval,fcompl,fcompll,maxspafreqsq,radfirstpeak)
                    !$omp do schedule(static)
                    do i = 1,nptcls_in_batch
                        iprec = batch_iprecs(i)
                        if( iprec == 0 ) cycle
                        if( precs(iprec)%pind == 0 ) cycle
                        ! apply CTF to image, stores CTF values
                        add_phshift = 0.
                        if( phaseplate ) add_phshift = precs(iprec)%phshift
                        call precs(iprec)%tfun%eval_and_apply_before_first_zero(cgrid_imgs_crop(i), ctfflag, logi_lims_crop, fdims_crop(1:2), tvals(:,:,i), &
                        & precs(iprec)%dfx, precs(iprec)%dfy, precs(iprec)%angast, add_phshift, maxspafreqsq)
                        radfirstpeak = ceiling( sqrt(maxspafreqsq) * real(ldim_croppd(1)) )
                        radfirstpeak = min( max(radfirstpeak,0), interp_shlim)
                        ! ML-regularization
                        if( l_ml_reg )then
                            do h = logi_lims_crop(1,1),radfirstpeak
                                do k = -radfirstpeak,radfirstpeak
                                    sh = nint(hyp(h,k))
                                    if( sh > radfirstpeak )cycle
                                    if( sh > interp_shlim )cycle
                                    sh = nint(reg_scale*sqrt(real(h*h)+real(k*k)))
                                    if( sh > sigma2_kfromto(2) ) cycle
                                    phys = cgrid_imgs(i)%comp_addr_phys(h,k)
                                    if( sh >= sigma2_kfromto(1) )then
                                        tvals(phys(1),phys(2),i) = tvals(phys(1),phys(2),i) / sqrt(eucl_sigma2_glob%sigma2_noise(sh,precs(iprec)%pind))
                                    else
                                        tvals(phys(1),phys(2),i) = tvals(phys(1),phys(2),i) / sqrt(eucl_sigma2_glob%sigma2_noise(1,precs(iprec)%pind))
                                    endif
                                enddo
                            enddo
                        endif
                        ! Rotation matrix
                        call rotmat2d(-precs(iprec)%e3, mat)
                        ! Interpolation
                        do h = logi_lims_crop(1,1),radfirstpeak
                            do k = -radfirstpeak,radfirstpeak
                                sh = nint(hyp(h,k))
                                if( sh > radfirstpeak )cycle
                                if( sh > interp_shlim )cycle
                                ! Rotation
                                loc        = matmul(real([h,k]),mat)
                                win_corner = floor(loc) ! bottom left corner
                                dist       = loc - real(win_corner)
                                ! Bi-linear interpolation
                                l     = cyci_1d(cyc_lims_cropR(:,1), win_corner(1))
                                ll    = cyci_1d(cyc_lims_cropR(:,1), win_corner(1)+1)
                                m     = cyci_1d(cyc_lims_cropR(:,2), win_corner(2))
                                mm    = cyci_1d(cyc_lims_cropR(:,2), win_corner(2)+1)
                                ! l, bottom left corner
                                phys   = cgrid_imgs_crop(i)%comp_addr_phys(l,m)
                                kw     = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
                                fcompl = kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                                tval   = kw * tvals(phys(1),phys(2),i)
                                ! l, bottom right corner
                                phys   = cgrid_imgs_crop(i)%comp_addr_phys(l,mm)
                                kw     = (1.-dist(1))*dist(2)
                                fcompl = fcompl + kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                                tval   = tval   + kw * tvals(phys(1),phys(2),i)
                                if( l < 0 ) fcompl = conjg(fcompl) ! conjugaison when required!
                                ! ll, upper left corner
                                phys    = cgrid_imgs_crop(i)%comp_addr_phys(ll,m)
                                kw      = dist(1)*(1.-dist(2))
                                fcompll =         kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                                tval    = tval  + kw * tvals(phys(1),phys(2),i)
                                ! ll, upper right corner
                                phys    = cgrid_imgs_crop(i)%comp_addr_phys(ll,mm)
                                kw      = dist(1)*dist(2)
                                fcompll = fcompll + kw * cgrid_imgs_crop(i)%get_cmat_at(phys(1), phys(2),1)
                                tval    = tval    + kw * tvals(phys(1),phys(2),i)
                                if( ll < 0 ) fcompll = conjg(fcompll) ! conjugaison when required!
                                ! update with interpolated values
                                phys = cgrid_imgs_crop(i)%comp_addr_phys(h,k)
                                cmats(phys(1),phys(2),i) = fcompl + fcompll
                                rhos(phys(1),phys(2),i)  = tval*tval
                            end do
                        end do
                    enddo
                    !$omp end do
                    ! Sum over classes
                    !$omp do schedule(static)
                    do icls = 1,ncls
                        ithr = omp_get_thread_num() + 1
                        do i = 1,nptcls_in_batch
                            iprec = batch_iprecs(i)
                            if( iprec == 0 ) cycle
                            if( precs(iprec)%pind == 0 ) cycle
                            if( precs(iprec)%class == icls )then
                                select case(precs(iprec)%eo)
                                    case(0,-1)
                                        call cavgs_even_wfilt(icls)%get_cmat_ptr(pcmat(ithr)%cmat)
                                        call ctfsqsums_even_wfilt(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
                                    case(1)
                                        call cavgs_odd_wfilt(icls)%get_cmat_ptr(pcmat(ithr)%cmat)
                                        call ctfsqsums_odd_wfilt(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
                                end select
                                pcmat(ithr)%cmat(:,:,1)   = pcmat(ithr)%cmat(:,:,1)   + precs(iprec)%pw * cmats(:,:,i)
                                prhomat(ithr)%cmat(:,:,1) = prhomat(ithr)%cmat(:,:,1) + precs(iprec)%pw * cmplx(rhos(:,:,i),0.0)
                            endif
                        enddo
                    enddo
                    !$omp end do
                    !$omp end parallel
                endif
            enddo ! end read batches loop
        enddo
        ! copy to part objects
        !$omp parallel do schedule(static) default(shared) proc_bind(close) private(icls)
        do icls = 1,ncls
            call cavgs_even_part(icls)%copy(cavgs_even(icls))
            call cavgs_odd_part(icls)%copy(cavgs_odd(icls))
            if( l_stream )then
                call cavgs_even_wfilt_bak(icls)%copy(cavgs_even_wfilt(icls))
                call cavgs_odd_wfilt_bak(icls)%copy(cavgs_odd_wfilt(icls))
            endif
        enddo
        !$omp end parallel do
        ! Cleanup
        call dstkio_r%kill
        do i = 1,READBUFFSZ
            call read_imgs(i)%kill
            call cgrid_imgs(i)%kill
            call cgrid_imgs_crop(i)%kill
        enddo
        deallocate(read_imgs,cgrid_imgs,cgrid_imgs_crop)
    end subroutine cavger_assemble_sums

    !>  \brief  merges the even/odd pairs and normalises the sums
    subroutine cavger_merge_eos_and_norm
        type(image) :: gridcorrection_img
        integer     :: icls, eo_pop(2), pop
        call cavger_prep_gridding_correction(gridcorrection_img)
        if( l_ml_reg )then
            ! Fourier components & CTF2 need to be stashed
            allocate(ctfsqsums_even_bak(ncls),ctfsqsums_odd_bak(ncls))
            if( l_stream ) allocate(ctfsqsums_even_wfilt_bak(ncls),ctfsqsums_odd_wfilt_bak(ncls))
        endif
        !$omp parallel do default(shared) private(icls,eo_pop,pop) schedule(static) proc_bind(close)
        do icls=1,ncls
            eo_pop = prev_eo_pops(icls,:) + eo_class_pop(icls)
            pop    = sum(eo_pop)
            if(pop == 0)then
                call cavgs_merged(icls)%zero_and_unflag_ft
                call cavgs_even(icls)%zero_and_unflag_ft
                call cavgs_odd(icls)%zero_and_unflag_ft
                call ctfsqsums_merged(icls)%zero_and_flag_ft
                if( l_stream )then
                    call cavgs_merged_wfilt(icls)%zero_and_unflag_ft
                    call cavgs_even_wfilt(icls)%zero_and_unflag_ft
                    call cavgs_odd_wfilt(icls)%zero_and_unflag_ft
                    call ctfsqsums_merged_wfilt(icls)%zero_and_flag_ft
                    call cavgs_even_wfilt_bak(icls)%zero_and_flag_ft
                    call cavgs_odd_wfilt_bak(icls)%zero_and_flag_ft
                endif
            else
                call cavgs_merged(icls)%zero_and_flag_ft
                call cavgs_merged(icls)%add(cavgs_even(icls))
                call cavgs_merged(icls)%add(cavgs_odd(icls))
                call ctfsqsums_merged(icls)%zero_and_flag_ft
                call ctfsqsums_merged(icls)%add(ctfsqsums_even(icls))
                call ctfsqsums_merged(icls)%add(ctfsqsums_odd(icls))
                if( l_stream )then
                    call cavgs_merged_wfilt(icls)%zero_and_flag_ft
                    call cavgs_merged_wfilt(icls)%add(cavgs_even_wfilt(icls))
                    call cavgs_merged_wfilt(icls)%add(cavgs_odd_wfilt(icls))
                    call ctfsqsums_merged_wfilt(icls)%zero_and_flag_ft
                    call ctfsqsums_merged_wfilt(icls)%add(ctfsqsums_even_wfilt(icls))
                    call ctfsqsums_merged_wfilt(icls)%add(ctfsqsums_odd_wfilt(icls))
                endif
                if( l_ml_reg )then
                    call ctfsqsums_even_bak(icls)%copy(ctfsqsums_even(icls))
                    call ctfsqsums_odd_bak(icls)%copy(ctfsqsums_odd(icls))
                    if( l_stream )then
                        call ctfsqsums_even_wfilt_bak(icls)%copy(ctfsqsums_even_wfilt(icls))
                        call ctfsqsums_odd_wfilt_bak(icls)%copy(ctfsqsums_odd_wfilt(icls))
                    endif
                endif
                ! (w*CTF)**2 density correction
                if(eo_pop(1) > 1)then
                    call cavgs_even(icls)%ctf_dens_correct(ctfsqsums_even(icls))
                    if( l_stream ) call cavgs_even_wfilt(icls)%ctf_dens_correct(ctfsqsums_even_wfilt(icls))
                endif
                if(eo_pop(2) > 1)then
                    call cavgs_odd(icls)%ctf_dens_correct(ctfsqsums_odd(icls))
                    if( l_stream ) call cavgs_odd_wfilt(icls)%ctf_dens_correct(ctfsqsums_odd_wfilt(icls))
                endif
                if(pop > 1)then
                    call cavgs_merged(icls)%ctf_dens_correct(ctfsqsums_merged(icls))
                    if( l_stream ) call cavgs_merged_wfilt(icls)%ctf_dens_correct(ctfsqsums_merged_wfilt(icls))
                endif
                call cavgs_even(icls)%ifft()
                call cavgs_odd(icls)%ifft()
                call cavgs_merged(icls)%ifft()
                if( l_stream )then
                    call cavgs_even_wfilt(icls)%ifft()
                    call cavgs_odd_wfilt(icls)%ifft()
                    call cavgs_merged_wfilt(icls)%ifft()
                endif
            endif
            call cavgs_even(icls)%clip_inplace(ldim_crop)
            call cavgs_odd(icls)%clip_inplace(ldim_crop)
            call cavgs_merged(icls)%clip_inplace(ldim_crop)
            if( l_stream )then
                call cavgs_even_wfilt(icls)%clip_inplace(ldim_crop)
                call cavgs_odd_wfilt(icls)%clip_inplace(ldim_crop)
                call cavgs_merged_wfilt(icls)%clip_inplace(ldim_crop)
            endif
            ! gridding correction
            call cavgs_even(icls)%div(gridcorrection_img)
            call cavgs_odd(icls)%div(gridcorrection_img)
            call cavgs_merged(icls)%div(gridcorrection_img)
            if( l_stream )then
                call cavgs_even_wfilt(icls)%div(gridcorrection_img)
                call cavgs_odd_wfilt(icls)%div(gridcorrection_img)
                call cavgs_merged_wfilt(icls)%div(gridcorrection_img)
            endif
        end do
        !$omp end parallel do
        call gridcorrection_img%kill
    end subroutine cavger_merge_eos_and_norm

     !>  \brief  calculates Fourier ring correlations
    subroutine cavger_calc_and_write_frcs_and_eoavg( fname, which_iter )
        class(string), intent(in) :: fname
        integer,       intent(in) :: which_iter
        type(image), allocatable  :: even_imgs(:), odd_imgs(:)
        real,        allocatable  :: frc(:)
        integer :: eo_pop(2), icls, find, pop, filtsz_crop
        filtsz_crop = cavgs_even(1)%get_filtsz()
        allocate(even_imgs(ncls), odd_imgs(ncls), frc(filtsz_crop))
        do icls=1,ncls
            call even_imgs(icls)%copy(cavgs_even(icls))
            call odd_imgs(icls)%copy(cavgs_odd(icls))
        end do
        if( l_ml_reg )then
            !$omp parallel do default(shared) private(icls,frc,find,pop,eo_pop) schedule(static) proc_bind(close)
            do icls=1,ncls
                eo_pop = prev_eo_pops(icls,:) + eo_class_pop(icls)
                pop    = sum(eo_pop)
                if( pop > 0 )then
                    call even_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                    call odd_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                    call even_imgs(icls)%fft()
                    call odd_imgs(icls)%fft()
                    call even_imgs(icls)%fsc(odd_imgs(icls), frc)
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                    find = build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1)
                    ! add noise term to denominator
                    call add_invtausq2rho(ctfsqsums_even(icls), frc)
                    call add_invtausq2rho(ctfsqsums_odd(icls), frc)
                    if( eo_pop(1) < 3 ) call ctfsqsums_even(icls)%add(1.)
                    if( eo_pop(2) < 3 ) call ctfsqsums_odd(icls)%add(1.)
                    if( l_stream )then
                        call add_invtausq2rho(ctfsqsums_even_wfilt(icls), frc)
                        call add_invtausq2rho(ctfsqsums_odd_wfilt(icls), frc)
                        if( eo_pop(1) < 3 ) call ctfsqsums_even_wfilt(icls)%add(1.)
                        if( eo_pop(2) < 3 ) call ctfsqsums_odd_wfilt(icls)%add(1.)
                    endif
                    ! re-generate sums
                    call cavgs_even(icls)%copy(cavgs_even_part(icls))
                    call cavgs_odd(icls)%copy(cavgs_odd_part(icls))
                    call cavgs_merged(icls)%copy(cavgs_even(icls))
                    call cavgs_merged(icls)%add(cavgs_odd(icls))
                    call ctfsqsums_merged(icls)%copy(ctfsqsums_even(icls))
                    call ctfsqsums_merged(icls)%add(ctfsqsums_odd(icls))
                    if( l_stream )then
                        call cavgs_even_wfilt(icls)%copy(cavgs_even_wfilt_bak(icls))
                        call cavgs_odd_wfilt(icls)%copy(cavgs_odd_wfilt_bak(icls))
                        call cavgs_merged_wfilt(icls)%copy(cavgs_even_wfilt(icls))
                        call cavgs_merged_wfilt(icls)%add(cavgs_odd_wfilt(icls))
                        call ctfsqsums_merged_wfilt(icls)%copy(ctfsqsums_even_wfilt(icls))
                        call ctfsqsums_merged_wfilt(icls)%add(ctfsqsums_odd_wfilt(icls))
                    endif
                    ! regularization
                    call cavgs_even(icls)%ctf_dens_correct(ctfsqsums_even(icls))
                    call cavgs_odd(icls)%ctf_dens_correct(ctfsqsums_odd(icls))
                    call cavgs_merged(icls)%ctf_dens_correct(ctfsqsums_merged(icls))
                    if( l_stream )then
                        call cavgs_even_wfilt(icls)%ctf_dens_correct(ctfsqsums_even_wfilt(icls))
                        call cavgs_odd_wfilt(icls)%ctf_dens_correct(ctfsqsums_odd_wfilt(icls))
                        call cavgs_merged_wfilt(icls)%ctf_dens_correct(ctfsqsums_merged_wfilt(icls))
                    endif
                    ! clip to original dimension
                    call cavgs_merged(icls)%ifft()
                    call cavgs_even(icls)%ifft()
                    call cavgs_odd(icls)%ifft()
                    call cavgs_even(icls)%clip_inplace(ldim_crop)
                    call cavgs_odd(icls)%clip_inplace(ldim_crop)
                    call cavgs_merged(icls)%clip_inplace(ldim_crop)
                    if( l_stream )then
                        call cavgs_merged_wfilt(icls)%ifft
                        call cavgs_even_wfilt(icls)%ifft()
                        call cavgs_odd_wfilt(icls)%ifft()
                        call cavgs_even_wfilt(icls)%clip_inplace(ldim_crop)
                        call cavgs_odd_wfilt(icls)%clip_inplace(ldim_crop)
                        call cavgs_merged_wfilt(icls)%clip_inplace(ldim_crop)
                    endif
                    ! average low-resolution info between eo pairs to keep things in register
                    call cavgs_merged(icls)%fft()
                    call cavgs_even(icls)%fft()
                    call cavgs_odd(icls)%fft()
                    call cavgs_even(icls)%insert_lowres_serial(cavgs_merged(icls), find)
                    call cavgs_odd(icls)%insert_lowres_serial(cavgs_merged(icls), find)
                    call cavgs_merged(icls)%ifft()
                    call cavgs_even(icls)%ifft()
                    call cavgs_odd(icls)%ifft()
                    if( l_stream )then
                        call cavgs_merged_wfilt(icls)%fft()
                        call cavgs_even_wfilt(icls)%fft()
                        call cavgs_odd_wfilt(icls)%fft()
                        call cavgs_even_wfilt(icls)%insert_lowres_serial(cavgs_merged_wfilt(icls), find)
                        call cavgs_odd_wfilt(icls)%insert_lowres_serial(cavgs_merged_wfilt(icls), find)
                        call cavgs_merged_wfilt(icls)%ifft
                        call cavgs_even_wfilt(icls)%ifft()
                        call cavgs_odd_wfilt(icls)%ifft()
                    endif
                else
                    call cavgs_even(icls)%clip_inplace(ldim_crop)
                    call cavgs_odd(icls)%clip_inplace(ldim_crop)
                    call cavgs_merged(icls)%clip_inplace(ldim_crop)
                    if( l_stream )then
                        call cavgs_even_wfilt(icls)%clip_inplace(ldim_crop)
                        call cavgs_odd_wfilt(icls)%clip_inplace(ldim_crop)
                        call cavgs_merged_wfilt(icls)%clip_inplace(ldim_crop)
                    endif
                endif
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(icls,frc,find) schedule(static) proc_bind(close)
            do icls=1,ncls
                call even_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                call odd_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                call even_imgs(icls)%fft()
                call odd_imgs(icls)%fft()
                call even_imgs(icls)%fsc(odd_imgs(icls), frc)
                call build_glob%clsfrcs%set_frc(icls, frc, 1)
                ! average low-resolution info between eo pairs to keep things in register
                find = build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1)
                call cavgs_merged(icls)%fft()
                call cavgs_even(icls)%fft()
                call cavgs_odd(icls)%fft()
                call cavgs_even(icls)%insert_lowres_serial(cavgs_merged(icls), find)
                call cavgs_odd(icls)%insert_lowres_serial(cavgs_merged(icls), find)
                call cavgs_merged(icls)%ifft()
                call cavgs_even(icls)%ifft()
                call cavgs_odd(icls)%ifft()
                if( l_stream )then
                    call cavgs_merged_wfilt(icls)%fft()
                    call cavgs_even_wfilt(icls)%fft()
                    call cavgs_odd_wfilt(icls)%fft()
                    call cavgs_even_wfilt(icls)%insert_lowres_serial(cavgs_merged_wfilt(icls), find)
                    call cavgs_odd_wfilt(icls)%insert_lowres_serial(cavgs_merged_wfilt(icls), find)
                    call cavgs_merged_wfilt(icls)%ifft
                    call cavgs_even_wfilt(icls)%ifft()
                    call cavgs_odd_wfilt(icls)%ifft()
                endif
            end do
            !$omp end parallel do
        endif
        ! write FRCs
        call build_glob%clsfrcs%write(fname)
        ! destruct
        do icls=1,ncls
            call even_imgs(icls)%kill
            call odd_imgs(icls)%kill
        end do
        deallocate(even_imgs, odd_imgs, frc)
    end subroutine cavger_calc_and_write_frcs_and_eoavg

    ! private function to add noise term to denominator
    subroutine add_invtausq2rho( ctfsqsum, frc )
        class(image),          intent(inout) :: ctfsqsum
        real,     allocatable, intent(in)    :: frc(:)
        real,     allocatable :: sig2(:), tau2(:), ssnr(:)
        integer,  allocatable :: cnt(:)
        real(dp), allocatable :: rsum(:)
        complex,      pointer :: pctfsqsum(:,:,:)
        real    :: cc, scale, pad_factor, invtau2, fudge
        integer :: flims(3,2), phys(2), h, k, sh, sz, reslim_ind
        fudge = params_glob%tau
        call ctfsqsum%get_cmat_ptr(pctfsqsum)
        flims = ctfsqsum%loop_lims(3)
        sz = size(frc)
        allocate(ssnr(0:sz), rsum(0:sz), cnt(0:sz), tau2(0:sz), sig2(0:sz))
        rsum = 0.d0
        cnt  = 0
        ssnr = 0.0
        tau2 = 0.0
        sig2 = 0.0
        scale = real(ldim_crop(1)) / real(ldim_croppd(1))
        pad_factor = 1.0 / scale**2
        ! SSNR
        do k = 1,sz
            cc      = max(0.001,frc(k))
            cc      = min(0.999,cc)
            ssnr(k) = fudge * cc / (1.-cc)
        enddo
        ! Noise
        !$omp parallel do collapse(2) default(shared) schedule(static)&
        !$omp private(h,k,phys,sh) proc_bind(close) reduction(+:cnt,rsum)
        do h = flims(1,1),flims(1,2)
            do k = flims(2,1),flims(2,2)
                sh = nint(scale * sqrt(real(h*h + k*k)))
                if( sh > sz ) cycle
                phys     = ctfsqsum%comp_addr_phys(h,k)
                cnt(sh)  = cnt(sh) + 1
                rsum(sh) = rsum(sh) + real(pctfsqsum(phys(1),phys(2),1),dp)
            enddo
        enddo
        !$omp end parallel do
        rsum = rsum * pad_factor
        where( rsum > 1.d-10 )
            sig2 = real(real(cnt,dp) / rsum)
        else where
            sig2 = 0.0
        end where
        ! Signal
        tau2 = ssnr * sig2
        ! add Tau2 inverse to denominator
        ! because signal assumed infinite at very low resolution there is no addition
        reslim_ind = max(6, calc_fourier_index(params_glob%hp, params_glob%box, params_glob%smpd))
        !$omp parallel do collapse(2) default(shared) schedule(static)&
        !$omp private(h,k,phys,sh,invtau2) proc_bind(close)
        do h = flims(1,1),flims(1,2)
            do k = flims(2,1),flims(2,2)
                sh = nint(scale*sqrt(real(h*h + k*k)))
                if( (sh < reslim_ind) .or. (sh > sz) ) cycle
                phys = ctfsqsum%comp_addr_phys(h, k)
                if( tau2(sh) > TINY)then
                    invtau2 = 1.0/(pad_factor*fudge*tau2(sh))
                else
                    invtau2 = min(1.e3, 1.e3 * real(pctfsqsum(phys(1),phys(2),1)))
                endif
                pctfsqsum(phys(1),phys(2),1) = pctfsqsum(phys(1),phys(2),1) + cmplx(invtau2,0.)
            enddo
        enddo
        !$omp end parallel do
    end subroutine add_invtausq2rho

    ! calculate classes intensity histogram distance to average histogram
    subroutine score_classes( os )
        use simple_histogram, only:histogram
        class(oris), intent(inout) :: os
        integer,   parameter :: NHISTBINS = 128
        type(image)          :: tmpimg
        type(histogram)      :: hists(ncls), hist_avg
        logical, allocatable :: lmsk(:,:,:), cls_mask(:)
        real    :: minmax(2), mean, tvd, overall_min, overall_max, std
        real    :: innerrad, mskrad
        integer :: icls,n
        cls_mask = nint(os%get_all('pop')) > 0
        n        = count(cls_mask)
        if( n == 0 ) return
        mskrad = real(ldim_crop(1)/2-1)
        call tmpimg%disc(ldim_crop, smpd_crop, mskrad, lmsk)
        overall_min =  huge(0.)
        overall_max = -999999.
        !$omp parallel private(icls,minmax,mean,std,tvd) proc_bind(close) default(shared)
        !$omp do reduction(min:overall_min) reduction(max:overall_max)
        do icls = 1,ncls
            if( cls_mask(icls) )then
                call cavgs_merged(icls)%stats(mean, std, minmax(2), minmax(1), tmpimg)
                overall_min = min(minmax(1),overall_min)
                overall_max = max(minmax(2),overall_max)
                call os%set(icls, 'min',  minmax(1))
                call os%set(icls, 'max',  minmax(2))
                call os%set(icls, 'mean', mean)
                call os%set(icls, 'var',  std*std)
            endif
        enddo
        !$omp end do
        !$omp do
        do icls = 1,ncls
            if( cls_mask(icls) )then
                call hists(icls)%new(cavgs_merged(icls), NHISTBINS,&
                &minmax=[overall_min, overall_max], radius=mskrad)
            endif
        enddo
        !$omp end do
        !$omp single
        call hist_avg%new(hists(findloc(cls_mask, .true., dim=1)))
        do icls = 1,ncls
            if( cls_mask(icls) ) call hist_avg%add(hists(icls))
        enddo
        call hist_avg%div(real(n))
        !$omp end single
        !$omp do
        do icls = 1,ncls
            if( cls_mask(icls) )then
                tvd = hist_avg%tvd(hists(icls))
                call os%set(icls, 'tvd', tvd)
                call hists(icls)%kill
            endif
        enddo
        !$omp end do
        !$omp end parallel
        innerrad = params_glob%msk_crop
        call tmpimg%ring(ldim_crop, smpd_crop, min(innerrad+COSMSKHALFWIDTH, mskrad), innerrad)
        !$omp parallel do private(icls,minmax,mean,std,tvd) proc_bind(close) default(shared)
        do icls = 1,ncls
            if( cls_mask(icls) )then
                call cavgs_merged(icls)%stats(mean, std, minmax(2), minmax(1), tmpimg)
                call os%set(icls, 'mskvar', std*std)
            endif
        enddo
        !$omp end parallel do
        call hist_avg%kill
        call tmpimg%kill
    end subroutine score_classes

    ! I/O

    !>  \brief  writes class averages to disk
    subroutine cavger_write( fname, which )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        type(string) :: fname_wfilt, fname_here
        integer :: icls
        fname_here  = fname
        fname_wfilt = add2fbody(fname_here, params_glob%ext, WFILT_SUFFIX)
        select case(which)
            case('even')
                if( l_stream ) fname_wfilt = add2fbody(fname_here, '_even'//params_glob%ext%to_char(), WFILT_SUFFIX)
                do icls=1,ncls
                    call cavgs_even(icls)%write(fname_here, icls)
                    if( l_stream ) call cavgs_even_wfilt(icls)%write(fname_wfilt, icls)
                end do
            case('odd')
                if( l_stream ) fname_wfilt = add2fbody(fname_here, '_odd'//params_glob%ext%to_char(), WFILT_SUFFIX)
                do icls=1,ncls
                    call cavgs_odd(icls)%write(fname_here, icls)
                    if( l_stream ) call cavgs_odd_wfilt(icls)%write(fname_wfilt, icls)
                end do
            case('merged')
                if( l_stream ) fname_wfilt = add2fbody(fname_here, params_glob%ext%to_char(), WFILT_SUFFIX)
                do icls=1,ncls
                    call cavgs_merged(icls)%write(fname_here, icls)
                    if( l_stream ) call cavgs_merged_wfilt(icls)%write(fname_wfilt, icls)
                end do
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
        call update_stats

    contains

        subroutine update_stats
            integer :: icls, cnt
            real    :: stats(4),minv,maxv,meanv,stdevv
            logical :: l_err
            stats(1)   = huge(stats(1))
            stats(2)   = -stats(1)
            stats(3:4) = 0.
            cnt        = 1
            do icls = 1,ncls
                select case(which)
                    case('even')
                        call cavgs_even(icls)%stats(meanv, stdevv, maxv, minv, errout=l_err)
                    case('odd')
                        call cavgs_odd(icls)%stats(meanv, stdevv, maxv, minv, errout=l_err)
                    case('merged')
                        call cavgs_merged(icls)%stats(meanv, stdevv, maxv, minv, errout=l_err)
                end select
                if( .not.l_err )then
                    cnt = cnt + 1
                    stats(1) = min(stats(1),minv)
                    stats(2) = max(stats(2),maxv)
                    stats(3) = stats(3) + meanv
                    stats(4) = stats(4) + stdevv**2.
                endif
            enddo
            if( cnt > 1 )then
                ! updates header, size, stack & removes box file
                stats(3) = stats(3) / real(cnt)
                stats(4) = sqrt(stats(4) / real(cnt))
                select case(which)
                    case('even')
                        call cavgs_even(1)%update_header_stats(fname,stats)
                    case('odd')
                        call cavgs_odd(1)%update_header_stats(fname,stats)
                    case('merged')
                        call cavgs_merged(1)%update_header_stats(fname,stats)
                end select
            endif
        end subroutine update_stats

    end subroutine cavger_write

    subroutine cavger_read_all()
        if( .not. file_exists(params_glob%refs) ) THROW_HARD('references (REFS) does not exist in cwd')
        call cavger_read(params_glob%refs, 'merged')
        if( file_exists(params_glob%refs_even) )then
            call cavger_read(params_glob%refs_even, 'even')
        else
            call cavger_read(params_glob%refs, 'even')
        endif
        if( file_exists(params_glob%refs_odd) )then
            call cavger_read(params_glob%refs_odd, 'odd')
        else
            call cavger_read(params_glob%refs, 'odd')
        endif
    end subroutine cavger_read_all

    !>  \brief  reads class averages from disk
    subroutine cavger_read( fname, which )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        class(image), pointer :: cavgs(:)
        type(stack_io)        :: stkio_r
        integer :: ldim_read(3), icls
        select case(trim(which))
            case('even')
                cavgs => cavgs_even
            case('odd')
                cavgs => cavgs_odd
            case('merged')
                cavgs => cavgs_merged
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
        ! read
        call stkio_r%open(fname, smpd_crop, 'read', bufsz=ncls)
        ldim_read = stkio_r%get_ldim()
        do icls = 1,ncls
            call cavgs(icls)%new(ldim_read,smpd_crop,wthreads=.false.)
            call stkio_r%read(icls, cavgs(icls))
        end do
        call stkio_r%close
        ! scale
        if( any(ldim_read /= ldim_crop) )then
            if( ldim_read(1) > ldim_crop(1) )then
                ! Cropping is not covered
                THROW_HARD('Incompatible cavgs dimensions! ; cavger_read')
            else if( ldim_read(1) < ldim_crop(1) )then
                ! Fourier padding
                !$omp parallel do proc_bind(close) schedule(static) default(shared) private(icls)
                do icls = 1,ncls
                    call cavgs(icls)%fft
                    call cavgs(icls)%pad_inplace(ldim_crop)
                    call cavgs(icls)%ifft
                end do
                !$omp end parallel do
            endif
        endif
        nullify(cavgs)
    end subroutine cavger_read

    !>  \brief  writes partial class averages to disk (distributed execution)
    subroutine cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        integer        :: icls, ldim_here(3)
        type(string)   :: cae, cao, cte, cto, caewf, caowf, ctewf, ctowf
        type(stack_io) :: stkio(4)
        logical        :: is_ft
        cae   = 'cavgs_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        cao   = 'cavgs_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        cte   = 'ctfsqsums_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        cto   = 'ctfsqsums_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        caewf = 'cavgs_even_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        caowf = 'cavgs_odd_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        ctewf = 'ctfsqsums_even_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        ctowf = 'ctfsqsums_odd_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        select case(trim(which))
            case('read')
                call stkio(1)%open(cae, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                call stkio(2)%open(cao, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                call stkio(3)%open(cte, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                call stkio(4)%open(cto, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                do icls=1,ncls
                    call stkio(1)%read(icls, cavgs_even(icls))
                    call stkio(2)%read(icls, cavgs_odd(icls))
                    call stkio(3)%read(icls, ctfsqsums_even(icls))
                    call stkio(4)%read(icls, ctfsqsums_odd(icls))
                end do
                if( l_stream )then
                    if( .not.file_exists(caewf) )then
                        ! fallback, we use the partial CTF classes instead
                        do icls=1,ncls
                            call cavgs_even_wfilt(icls)%set(cavgs_even(icls))
                            call cavgs_odd_wfilt(icls)%set(cavgs_odd(icls))
                            call ctfsqsums_even_wfilt(icls)%set(ctfsqsums_odd(icls))
                            call ctfsqsums_odd_wfilt(icls)%set(ctfsqsums_odd(icls))
                        end do
                    else
                        call stkio(1)%open(caewf, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                        call stkio(2)%open(caowf, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                        call stkio(3)%open(ctewf, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                        call stkio(4)%open(ctowf, smpd_crop, 'read', bufsz=ncls, is_ft=.true.)
                        do icls=1,ncls
                            call stkio(1)%read(icls, cavgs_even_wfilt(icls))
                            call stkio(2)%read(icls, cavgs_odd_wfilt(icls))
                            call stkio(3)%read(icls, ctfsqsums_even_wfilt(icls))
                            call stkio(4)%read(icls, ctfsqsums_odd_wfilt(icls))
                        end do
                    endif
                endif
                call stkio(1)%close
                call stkio(2)%close
                call stkio(3)%close
                call stkio(4)%close
            case('write')
                is_ft = cavgs_even_part(1)%is_ft()
                ldim_here  = cavgs_even_part(1)%get_ldim()
                call stkio(1)%open(cae, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(2)%open(cao, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(3)%open(cte, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(4)%open(cto, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                do icls=1,ncls
                    call stkio(1)%write(icls, cavgs_even_part(icls))
                    call stkio(2)%write(icls, cavgs_odd_part(icls))
                    call stkio(3)%write(icls, ctfsqsums_even(icls))
                    call stkio(4)%write(icls, ctfsqsums_odd(icls))
                end do
                if( l_stream )then
                    call stkio(1)%open(caewf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    call stkio(2)%open(caowf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    call stkio(3)%open(ctewf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    call stkio(4)%open(ctowf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    do icls=1,ncls
                        call stkio(1)%write(icls, cavgs_even_wfilt_bak(icls))
                        call stkio(2)%write(icls, cavgs_odd_wfilt_bak(icls))
                        call stkio(3)%write(icls, ctfsqsums_even_wfilt(icls))
                        call stkio(4)%write(icls, ctfsqsums_odd_wfilt(icls))
                    end do
                endif
                call stkio(1)%close
                call stkio(2)%close
                call stkio(3)%close
                call stkio(4)%close
            case DEFAULT
                THROW_HARD('unknown which flag; only read & write supported; cavger_readwrite_partial_sums')
        end select
        call cae%kill
        call cao%kill
        call cte%kill
        call cto%kill
        call caewf%kill
        call caowf%kill
        call ctewf%kill
        call ctowf%kill
    end subroutine cavger_readwrite_partial_sums

    subroutine cavger_apply_weights( w )
        real, intent(in) :: w
        integer :: icls
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)%mul(w)
            call ctfsqsums_even(icls)%mul(w)
            call cavgs_odd(icls)%mul(w)
            call ctfsqsums_odd(icls)%mul(w)
            if( l_stream )then
                call cavgs_even_wfilt(icls)%mul(w)
                call ctfsqsums_even_wfilt(icls)%mul(w)
                call cavgs_odd_wfilt(icls)%mul(w)
                call ctfsqsums_odd_wfilt(icls)%mul(w)
            endif
        end do
        !$omp end parallel do
    end subroutine cavger_apply_weights

    !>  \brief  re-generates the object after distributed execution
    subroutine cavger_assemble_sums_from_parts
        use simple_imgfile, only: imgfile
        complex(kind=c_float_complex), pointer :: cmat_ptr1(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr2(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr3(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr4(:,:,:) => null()
        real(kind=c_float),            pointer :: rmat_ptr(:,:,:)  => null()  !< image pixels/voxels (in data)
        integer(timer_int_kind)       ::  t_init,  t_io,  t_workshare_sum,  t_set_sums,  t_merge_eos_and_norm,  t_tot
        real(timer_int_kind)          :: rt_init, rt_io, rt_workshare_sum, rt_set_sums, rt_merge_eos_and_norm, rt_tot
        complex, allocatable :: csums(:,:,:,:)
        type(string)  :: cae, cao, cte, cto, benchfname
        type(image)   :: imgs4read(4)
        type(imgfile) :: ioimg(4)
        integer       :: ipart, icls, array_shape(3), ldim_here(3), fnr, i
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call init_cavgs_sums
        ! construct image objs for read/sum
        ldim_here    = ldim_croppd
        ldim_here(3) = ncls
        call imgs4read(1)%new(ldim_here, smpd_crop)
        call imgs4read(2)%new(ldim_here, smpd_crop)
        call imgs4read(3)%new(ldim_here, smpd_crop)
        call imgs4read(4)%new(ldim_here, smpd_crop)
        call imgs4read(1)%set_ft(.true.)
        call imgs4read(2)%set_ft(.true.)
        call imgs4read(3)%set_ft(.true.)
        call imgs4read(4)%set_ft(.true.)
        call imgs4read(1)%get_cmat_ptr(cmat_ptr1)
        call imgs4read(2)%get_cmat_ptr(cmat_ptr2)
        call imgs4read(3)%get_cmat_ptr(cmat_ptr3)
        call imgs4read(4)%get_cmat_ptr(cmat_ptr4)
        ! construct complex matrices for parallel summation
        array_shape = imgs4read(1)%get_array_shape()
        allocate(csums(array_shape(1),array_shape(2),array_shape(3),4), source=cmplx(0.,0.))
        if( L_BENCH_GLOB )then
            ! end of init
            rt_init = toc(t_init)
            ! initialise incremental timers before loop
            rt_io            = 0.
            rt_workshare_sum = 0.
        endif
        do ipart=1,params_glob%nparts
            if( L_BENCH_GLOB ) t_io = tic()
            ! look for files
            cae = 'cavgs_even_part'    //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            cao = 'cavgs_odd_part'     //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            cte = 'ctfsqsums_even_part'//int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            cto = 'ctfsqsums_odd_part' //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            call ioimg(1)%open(cae, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
            call ioimg(2)%open(cao, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
            call ioimg(3)%open(cte, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
            call ioimg(4)%open(cto, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
            !$omp parallel do default(shared) private(i,rmat_ptr) schedule(static) num_threads(4)
            do i = 1, 4
                call imgs4read(i)%get_rmat_ptr(rmat_ptr)
                call ioimg(i)%rSlices(1,ldim_here(1),rmat_ptr,is_mrc=.true.)
            end do
            !$omp end parallel do
            call ioimg(1)%close
            call ioimg(2)%close
            call ioimg(3)%close
            call ioimg(4)%close
            call cae%kill
            call cao%kill
            call cte%kill
            call cto%kill
            if( L_BENCH_GLOB )then
                rt_io = rt_io + toc(t_io)
                t_workshare_sum = tic()
            endif
            ! parallel summation
            !$omp parallel workshare proc_bind(close)
            csums(:,:,:,1) = csums(:,:,:,1) + cmat_ptr1(:,:,:)
            csums(:,:,:,2) = csums(:,:,:,2) + cmat_ptr2(:,:,:)
            csums(:,:,:,3) = csums(:,:,:,3) + cmat_ptr3(:,:,:)
            csums(:,:,:,4) = csums(:,:,:,4) + cmat_ptr4(:,:,:)
            !$omp end parallel workshare
            if( L_BENCH_GLOB ) rt_workshare_sum = rt_workshare_sum + toc(t_workshare_sum)
        end do
        if( L_BENCH_GLOB ) t_set_sums = tic()
        ! update image objects in parallel
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)    %set_cmat(csums(:,:,icls,1))
            call cavgs_odd(icls)     %set_cmat(csums(:,:,icls,2))
            call ctfsqsums_even(icls)%set_cmat(csums(:,:,icls,3))
            call ctfsqsums_odd(icls) %set_cmat(csums(:,:,icls,4))
        end do
        !$omp end parallel do
        if( l_stream )then
            ! reset sum
            !$omp parallel workshare proc_bind(close)
            csums = cmplx(0.0,0.0)
            !$omp end parallel workshare
            do ipart=1,params_glob%nparts
                if( L_BENCH_GLOB ) t_io = tic()
                ! look for files
                cae = 'cavgs_even_wfilt_part'    //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
                cao = 'cavgs_odd_wfilt_part'     //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
                cte = 'ctfsqsums_even_wfilt_part'//int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
                cto = 'ctfsqsums_odd_wfilt_part' //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
                call ioimg(1)%open(cae, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
                call ioimg(2)%open(cao, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
                call ioimg(3)%open(cte, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
                call ioimg(4)%open(cto, ldim_here, smpd, formatchar='M', readhead=.false., rwaction='read')
                !$omp parallel do default(shared) private(i,rmat_ptr) schedule(static) num_threads(4)
                do i = 1, 4
                    call imgs4read(i)%get_rmat_ptr(rmat_ptr)
                    call ioimg(i)%rSlices(1,ldim_here(1),rmat_ptr,is_mrc=.true.)
                end do
                !$omp end parallel do
                call ioimg(1)%close
                call ioimg(2)%close
                call ioimg(3)%close
                call ioimg(4)%close
                call cae%kill
                call cao%kill
                call cte%kill
                call cto%kill
                if( L_BENCH_GLOB )then
                    rt_io = rt_io + toc(t_io)
                    t_workshare_sum = tic()
                endif
                ! parallel summation
                !$omp parallel workshare proc_bind(close)
                csums(:,:,:,1) = csums(:,:,:,1) + cmat_ptr1(:,:,:)
                csums(:,:,:,2) = csums(:,:,:,2) + cmat_ptr2(:,:,:)
                csums(:,:,:,3) = csums(:,:,:,3) + cmat_ptr3(:,:,:)
                csums(:,:,:,4) = csums(:,:,:,4) + cmat_ptr4(:,:,:)
                !$omp end parallel workshare
                if( L_BENCH_GLOB ) rt_workshare_sum = rt_workshare_sum + toc(t_workshare_sum)
            end do
            if( L_BENCH_GLOB ) t_set_sums = tic()
            ! update image objects in parallel
            !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
            do icls=1,ncls
                call cavgs_even_wfilt(icls)    %set_cmat(csums(:,:,icls,1))
                call cavgs_odd_wfilt(icls)     %set_cmat(csums(:,:,icls,2))
                call ctfsqsums_even_wfilt(icls)%set_cmat(csums(:,:,icls,3))
                call ctfsqsums_odd_wfilt(icls) %set_cmat(csums(:,:,icls,4))
            end do
            !$omp end parallel do
        endif
        if( L_BENCH_GLOB ) rt_set_sums = rt_set_sums + toc(t_set_sums)
        ! destruct
        call imgs4read(1)%kill
        call imgs4read(2)%kill
        call imgs4read(3)%kill
        call imgs4read(4)%kill
        deallocate(csums)
        ! merge eo-pairs and normalize
        if( L_BENCH_GLOB ) t_merge_eos_and_norm = tic()
        call cavger_merge_eos_and_norm()
        if( L_BENCH_GLOB )then
            rt_merge_eos_and_norm = toc(t_merge_eos_and_norm)
            rt_tot                = toc(t_tot)
            benchfname = 'CAVGASSEMBLE_BENCH.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation       : ', rt_init
            write(fnr,'(a,1x,f9.2)') 'I/O                  : ', rt_io
            write(fnr,'(a,1x,f9.2)') 'workshare sum        : ', rt_workshare_sum
            write(fnr,'(a,1x,f9.2)') 'set sums             : ', rt_set_sums
            write(fnr,'(a,1x,f9.2)') 'merge eo-pairs & norm: ', rt_merge_eos_and_norm
            write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation        : ', (rt_init/rt_tot)               * 100.
            write(fnr,'(a,1x,f9.2)') 'I/O                   : ', (rt_io/rt_tot)                 * 100.
            write(fnr,'(a,1x,f9.2)') 'workshare sum         : ', (rt_workshare_sum/rt_tot)      * 100.
            write(fnr,'(a,1x,f9.2)') 'set sums              : ', (rt_set_sums/rt_tot)           * 100.
            write(fnr,'(a,1x,f9.2)') 'merge eo-pairs & norm : ', (rt_merge_eos_and_norm/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for       : ',&
            &((rt_init+rt_io+rt_workshare_sum+rt_set_sums+rt_merge_eos_and_norm)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine cavger_assemble_sums_from_parts

    !>  \brief  corrects for Fourier domain bilinear interpolation
    subroutine cavger_prep_gridding_correction( img )
        class(image), intent(inout) :: img
        real    :: center(3),dist(2),pid,sinc,pad_sc
        integer :: i,j
        call img%new(ldim_crop,smpd_crop)
        center = real(ldim_crop/2 + 1)
        pad_sc = 1. / real(ldim_croppd(1))
        do i = 1, ldim_crop(1)
            dist(1) = pad_sc * (real(i) - center(1))
            do j = 1, ldim_crop(2)
                dist(2) = pad_sc * (real(j) - center(2))
                pid     = PI * sqrt(sum(dist**2.))
                if( pid < TINY )then
                    sinc = 1.
                else
                    sinc = sin(pid) / pid
                endif
                call img%set([i,j,1], sinc*sinc)
            enddo
        enddo
    end subroutine cavger_prep_gridding_correction

    ! destructor

    !>  \brief  is a destructor
    subroutine cavger_kill( dealloccavgs )
        logical, optional, intent(in) :: dealloccavgs
        if( present(dealloccavgs) )then
            if( dealloccavgs ) call dealloc_cavgs_sums
        else
            call dealloc_cavgs_sums
        endif
        if( allocated(pptcl_mask) ) deallocate(pptcl_mask, prev_eo_pops,precs)
    end subroutine cavger_kill

    subroutine dealloc_cavgs_sums
        integer :: icls
        if( allocated(cavgs_merged) )then
            do icls=1,ncls
                call cavgs_even(icls)%kill
                call cavgs_odd(icls)%kill
                call cavgs_merged(icls)%kill
                call ctfsqsums_even(icls)%kill
                call ctfsqsums_odd(icls)%kill
                call ctfsqsums_merged(icls)%kill
                if( l_stream )then
                    call cavgs_even_wfilt(icls)%kill
                    call cavgs_odd_wfilt(icls)%kill
                    call cavgs_merged_wfilt(icls)%kill
                    call ctfsqsums_even_wfilt(icls)%kill
                    call ctfsqsums_odd_wfilt(icls)%kill
                    call ctfsqsums_merged_wfilt(icls)%kill
                endif
            end do
            deallocate(cavgs_even, cavgs_odd, cavgs_merged,&
            &ctfsqsums_even,ctfsqsums_odd, ctfsqsums_merged)
        endif
        if( allocated(cavgs_even_part) )then
            do icls=1,ncls
                call cavgs_even_part(icls)%kill
                call cavgs_odd_part(icls)%kill
            enddo
            deallocate(cavgs_even_part,cavgs_odd_part)
        endif
        if( allocated(cavgs_even_wfilt_bak) )then
            do icls=1,ncls
                call cavgs_even_wfilt_bak(icls)%kill
                call cavgs_odd_wfilt_bak(icls)%kill
            enddo
            deallocate(cavgs_even_wfilt_bak,cavgs_odd_wfilt_bak)
        endif
        if( allocated(ctfsqsums_even_bak) )then
            do icls=1,ncls
                call ctfsqsums_even_bak(icls)%kill
                call ctfsqsums_odd_bak(icls)%kill
                if( l_stream )then
                    call ctfsqsums_even_wfilt_bak(icls)%kill
                    call ctfsqsums_odd_wfilt_bak(icls)%kill
                endif
            enddo
            deallocate(ctfsqsums_even_bak,ctfsqsums_odd_bak)
            if( l_stream ) deallocate(ctfsqsums_even_wfilt_bak,ctfsqsums_odd_wfilt_bak)
        endif
        istart = 0
        iend   = 0
        partsz = 0
        ncls   = 0
    end subroutine dealloc_cavgs_sums

    ! PUBLIC UTILITIES

    subroutine transform_ptcls( spproj, oritype, icls, timgs, pinds, phflip, cavg, imgs_ori, just_transf)
        use simple_sp_project,          only: sp_project
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimgbatch
        class(sp_project),                  intent(inout) :: spproj
        character(len=*),                   intent(in)    :: oritype
        integer,                            intent(in)    :: icls
        type(image),           allocatable, intent(inout) :: timgs(:)
        integer,               allocatable, intent(inout) :: pinds(:)
        logical,     optional,              intent(in)    :: phflip
        type(image), optional,              intent(inout) :: cavg
        type(image), optional, allocatable, intent(inout) :: imgs_ori(:)
        logical,     optional,              intent(in)    :: just_transf
        class(oris),  pointer :: pos
        type(image)           :: img(nthr_glob), timg(nthr_glob)
        type(ctfparams)       :: ctfparms
        type(ctf)             :: tfun
        type(string)          :: str
        complex :: fcompl, fcompll
        real    :: mat(2,2), shift(2), loc(2), dist(2), e3, kw
        integer :: logi_lims(3,2),cyc_lims(3,2),cyc_limsR(2,2),phys(2),win_corner(2)
        integer :: i,iptcl, l,ll,m,mm, pop, h,k, ithr
        logical :: l_phflip, l_imgs, l_just_transf
        l_imgs        = .false.
        l_just_transf = .false.
        if( present(imgs_ori) )then
            if( present(just_transf) )then
                l_imgs        = .false.
                l_just_transf = just_transf
            else
                l_imgs        = .true.
            endif
        endif
        if( allocated(timgs) )then
            do i = 1,size(timgs)
                call timgs(i)%kill
            enddo
            deallocate(timgs)
        endif
        if( l_imgs )then
            if( allocated(imgs_ori) )then
                do i = 1,size(imgs_ori)
                    call imgs_ori(i)%kill
                enddo
                deallocate(imgs_ori)
            endif
        endif
        if(present(cavg)) call cavg%kill
        select case(trim(oritype))
            case('ptcl2D')
                str = 'class'
            case('ptcl3D')
                str = 'proj'
            case DEFAULT
                THROW_HARD('ORITYPE not supported!')
        end select
        call spproj%ptr2oritype( oritype, pos )
        if( allocated(pinds) ) deallocate(pinds)
        call pos%get_pinds(icls, str%to_char(), pinds)
        if( .not.(allocated(pinds)) ) return
        pop = size(pinds)
        if( pop == 0 ) return
        l_phflip = .false.
        if( present(phflip) ) l_phflip = phflip
        if( l_phflip )then
            select case( spproj%get_ctfflag_type(oritype, pinds(1)) )
                case(CTFFLAG_NO)
                    THROW_HARD('NO CTF INFORMATION COULD BE FOUND')
                case(CTFFLAG_FLIP)
                    THROW_WARN('Images have already been phase-flipped, phase flipping is deactivated')
                    l_phflip = .false.
                case(CTFFLAG_YES)
                    ! all good
                case DEFAULT
                    THROW_HARD('UNSUPPORTED CTF FLAG')
                end select
        endif
        allocate(timgs(pop))
        do i = 1,size(timgs)
            call timgs(i)%new([params_glob%box,params_glob%box,1],params_glob%smpd, wthreads=.false.)
        enddo
        if( l_imgs )then
            allocate(imgs_ori(pop))
            do i = 1,size(imgs_ori)
                call imgs_ori(i)%new([params_glob%box,params_glob%box,1],params_glob%smpd, wthreads=.false.)
            enddo
        endif
        ! temporary objects
        call prepimgbatch(pop)
        do ithr = 1, nthr_glob
            call  img(ithr)%new([params_glob%boxpd,params_glob%boxpd,1],params_glob%smpd, wthreads=.false.)
            call timg(ithr)%new([params_glob%boxpd,params_glob%boxpd,1],params_glob%smpd, wthreads=.false.)
        end do
        logi_lims      = img(1)%loop_lims(2)
        cyc_lims       = img(1)%loop_lims(3)
        cyc_limsR(:,1) = cyc_lims(1,:)
        cyc_limsR(:,2) = cyc_lims(2,:)
        call discrete_read_imgbatch(pop, pinds(:), [1,pop])
        !$omp parallel do private(i,ithr,iptcl,shift,e3,ctfparms,tfun,mat,h,k,loc,win_corner,dist,l,ll,m,mm,phys,kw,fcompl,fcompll) &
        !$omp default(shared) schedule(static) proc_bind(close)
        do i = 1,pop
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(i)
            shift = pos%get_2Dshift(iptcl)
            e3    = pos%e3get(iptcl)
            call img(ithr)%zero_and_flag_ft
            call timg(ithr)%zero_and_flag_ft
            ! normalisation
            if( l_just_transf )then
                call imgs_ori(i)%pad_fft(img(ithr))
            else
                call build_glob%imgbatch(i)%norm_noise_pad_fft(build_glob%lmsk,img(ithr))
            endif
            if( l_imgs )then
                call img(ithr)%ifft
                call img(ithr)%clip(imgs_ori(i))
                call img(ithr)%fft
            endif
            ! optional phase-flipping
            if( l_phflip )then
                ctfparms = spproj%get_ctfparams(oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call tfun%apply_serial(img(ithr), 'flip', ctfparms)
            endif
            ! shift
            call img(ithr)%shift2Dserial(-shift)
            ! particle rotation
            call rotmat2d(-e3, mat)
            do h = logi_lims(1,1),logi_lims(1,2)
                do k = logi_lims(2,1),logi_lims(2,2)
                    ! Rotation
                    loc        = matmul(real([h,k]),mat)
                    win_corner = floor(loc) ! bottom left corner
                    dist       = loc - real(win_corner)
                    ! Bi-linear interpolation
                    l      = cyci_1d(cyc_limsR(:,1), win_corner(1))
                    ll     = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
                    m      = cyci_1d(cyc_limsR(:,2), win_corner(2))
                    mm     = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
                    ! l, bottom left corner
                    phys   = img(ithr)%comp_addr_phys(l,m)
                    kw     = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
                    fcompl = kw * img(ithr)%get_cmat_at(phys(1), phys(2),1)
                    ! l, bottom right corner
                    phys   = img(ithr)%comp_addr_phys(l,mm)
                    kw     = (1.-dist(1))*dist(2)
                    fcompl = fcompl + kw * img(ithr)%get_cmat_at(phys(1), phys(2),1)
                    if( l < 0 ) fcompl = conjg(fcompl) ! conjugation when required!
                    ! ll, upper left corner
                    phys    = img(ithr)%comp_addr_phys(ll,m)
                    kw      = dist(1)*(1.-dist(2))
                    fcompll = kw * img(ithr)%get_cmat_at(phys(1), phys(2),1)
                    ! ll, upper right corner
                    phys    = img(ithr)%comp_addr_phys(ll,mm)
                    kw      = dist(1)*dist(2)
                    fcompll = fcompll + kw * img(ithr)%get_cmat_at(phys(1), phys(2),1)
                    if( ll < 0 ) fcompll = conjg(fcompll) ! conjugation when required!
                    ! update with interpolated values
                    phys = img(ithr)%comp_addr_phys(h,k)
                    call timg(ithr)%set_cmat_at(phys(1),phys(2),1, fcompl + fcompll)
                end do
            end do
            call timg(ithr)%ifft
            call timg(ithr)%clip(timgs(i))
        enddo
        !$omp end parallel do
        if( present(cavg) )then
            call cavg%copy(timgs(1))
            do i =2,pop,1
                call cavg%add(timgs(i))
            enddo
            call cavg%div(real(pop))
        endif
        do ithr = 1, nthr_glob
            call img(ithr)%kill
            call timg(ithr)%kill
        end do
        nullify(pos)
    end subroutine transform_ptcls

end module simple_classaverager
