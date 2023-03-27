module simple_classaverager
include 'simple_lib.f08'
!$ use omp_lib
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
use simple_ctf,        only: ctf
use simple_image,      only: image, image_ptr
use simple_stack_io,   only: stack_io
use simple_euclid_sigma2
use simple_fsc
implicit none

public :: cavger_new, cavger_transf_oridat, cavger_gen2Dclassdoc, cavger_assemble_sums,&
cavger_merge_eos_and_norm, cavger_calc_and_write_frcs_and_eoavg, cavger_write, cavger_read,&
cavger_readwrite_partial_sums, cavger_assemble_sums_from_parts, cavger_kill, cavgs_even, cavgs_odd, cavgs_merged,&
cavger_read_euclid_sigma2
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

integer                        :: ctfflag                   !< ctf flag <yes=1|no=0|flip=2>
integer                        :: istart      = 0, iend = 0 !< particle index range
integer                        :: partsz      = 0           !< size of partition
integer                        :: ncls        = 0           !< # classes
integer                        :: ldim(3)        = [0,0,0] !< logical dimension of image
integer                        :: ldim_crop(3)   = [0,0,0] !< logical dimension of cropped image
integer                        :: ldim_pd(3)     = [0,0,0] !< logical dimension of image, padded
integer                        :: ldim_croppd(3) = [0,0,0] !< logical dimension of cropped image, padded
real                           :: smpd       = 0.          !< sampling distance
real                           :: smpd_crop  = 0.          !< cropped sampling distance
type(ptcl_record), allocatable :: precs(:)                 !< particle records
type(image),       allocatable :: cavgs_even(:)            !< class averages
type(image),       allocatable :: cavgs_odd(:)             !< -"-
type(image),       allocatable :: cavgs_even_wfilt(:)      !< class averages wiener filtered
type(image),       allocatable :: cavgs_odd_wfilt(:)       !< -"-
type(image),       allocatable :: cavgs_merged(:)          !< -"-
type(image),       allocatable :: cavgs_merged_wfilt(:)    !< -"-
type(image),       allocatable :: ctfsqsums_even(:)        !< CTF**2 sums for Wiener normalisation
type(image),       allocatable :: ctfsqsums_odd(:)         !< -"-
type(image),       allocatable :: ctfsqsums_even_wfilt(:)  !< CTF**2 sums for Wiener normalisation
type(image),       allocatable :: ctfsqsums_odd_wfilt(:)   !< -"-
type(image),       allocatable :: ctfsqsums_merged(:)      !< -"-
type(image),       allocatable :: ctfsqsums_merged_wfilt(:)!< -"-
type(euclid_sigma2)            :: eucl_sigma
type(image),       allocatable :: cavgs_even_bak(:)            !< class averages
type(image),       allocatable :: cavgs_odd_bak(:)             !< -"-
type(image),       allocatable :: cavgs_even_wfilt_bak(:)      !< class averages wiener filtered
type(image),       allocatable :: cavgs_odd_wfilt_bak(:)       !< -"-
type(image),       allocatable :: ctfsqsums_even_bak(:)        !< CTF**2 sums for Wiener normalisation
type(image),       allocatable :: ctfsqsums_odd_bak(:)         !< -"-
type(image),       allocatable :: ctfsqsums_even_wfilt_bak(:)  !< CTF**2 sums for Wiener normalisation
type(image),       allocatable :: ctfsqsums_odd_wfilt_bak(:)   !< -"-
integer,           allocatable :: prev_eo_pops(:,:)
logical,           allocatable :: pptcl_mask(:)
logical                        :: l_stream      = .false.  !< flag for cluster2D_stream
logical                        :: phaseplate    = .false.  !< Volta phaseplate images or not
logical                        :: l_ml_reg      = .false.
logical                        :: exists        = .false.  !< to flag instance existence

integer(timer_int_kind) :: t_class_loop,t_batch_loop, t_gridding, t_init, t_tot
real(timer_int_kind)    :: rt_class_loop,rt_batch_loop, rt_gridding, rt_init, rt_tot
character(len=STDLEN)   :: benchfname

contains

    subroutine cavger_new( ptcl_mask )
        logical, optional, intent(in) :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer :: icls
        ! destruct possibly pre-existing instance
        call cavger_kill
        if( present(ptcl_mask) )then
            allocate(pptcl_mask(params_glob%fromp:params_glob%top), source=ptcl_mask)
        else
            allocate(pptcl_mask(params_glob%fromp:params_glob%top), source=.true.)
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
        allocate(precs(partsz), cavgs_even(ncls), cavgs_odd(ncls),&
        &cavgs_merged(ncls), ctfsqsums_even(ncls),&
        &ctfsqsums_odd(ncls), ctfsqsums_merged(ncls), prev_eo_pops(ncls,2))
        prev_eo_pops = 0
        if( l_stream )then
            allocate(cavgs_even_wfilt(ncls), cavgs_odd_wfilt(ncls), ctfsqsums_merged_wfilt(ncls),&
            &ctfsqsums_even_wfilt(ncls), ctfsqsums_odd_wfilt(ncls), cavgs_merged_wfilt(ncls))
        endif
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
            call cavgs_odd(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
            call cavgs_merged(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
            call ctfsqsums_even(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
            call ctfsqsums_odd(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
            call ctfsqsums_merged(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
            if( l_stream )then
                call cavgs_merged_wfilt(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
                call cavgs_even_wfilt(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
                call cavgs_odd_wfilt(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
                call ctfsqsums_even_wfilt(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
                call ctfsqsums_odd_wfilt(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
                call ctfsqsums_merged_wfilt(icls)%new(ldim_croppd,params_glob%smpd_crop,wthreads=.false.)
            endif
        end do
        !$omp end parallel do
        ! flag existence
        exists = .true.
    end subroutine cavger_new

    ! setters/getters

    !>  \brief  transfers orientation data to the instance
    subroutine cavger_transf_oridat( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        type(ctfparams)   :: ctfvars(nthr_glob)
        integer           :: i, icls, cnt, iptcl, ithr, stkind
        ! build index map
        cnt = 0
        do iptcl=istart,iend
            cnt = cnt + 1
            ! exclusion
            precs(cnt)%pind = 0
            if( spproj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
            if( spproj%os_ptcl2D%get(iptcl,'w') < TINY ) cycle
            if( .not.pptcl_mask(iptcl)                 ) cycle
            precs(cnt)%pind = iptcl
        enddo
        ! fetch data from project
        !$omp parallel do default(shared) private(cnt,iptcl,ithr,stkind) schedule(static) proc_bind(close)
        do cnt = 1,partsz
            iptcl              = precs(cnt)%pind
            if( iptcl == 0 ) cycle
            ithr               = omp_get_thread_num() + 1
            precs(cnt)%eo      = nint(spproj%os_ptcl2D%get(iptcl,'eo'))
            precs(cnt)%pw      = spproj%os_ptcl2D%get(iptcl,'w')
            ctfvars(ithr)      = spproj%get_ctfparams('ptcl2D',iptcl)
            precs(cnt)%tfun    = ctf(params_glob%smpd_crop, ctfvars(ithr)%kv, ctfvars(ithr)%cs, ctfvars(ithr)%fraca)
            precs(cnt)%dfx     = ctfvars(ithr)%dfx
            precs(cnt)%dfy     = ctfvars(ithr)%dfy
            precs(cnt)%angast  = ctfvars(ithr)%angast
            precs(cnt)%phshift = 0.
            if( phaseplate ) precs(cnt)%phshift = ctfvars(ithr)%phshift
            precs(cnt)%class   = spproj%os_ptcl2D%get_class(iptcl)
            precs(cnt)%e3      = spproj%os_ptcl2D%e3get(iptcl)
            precs(cnt)%shift   = spproj%os_ptcl2D%get_2Dshift(iptcl)
            call spproj%map_ptcl_ind2stk_ind(params_glob%oritype, iptcl, stkind, precs(cnt)%ind_in_stk)
        end do
        !$omp end parallel do
        prev_eo_pops = 0
        if( l_stream .and. spproj%os_cls2D%get_noris() == ncls )then
            do i = 1,ncls
                icls = spproj%os_cls2D%get_class(i)
                if( .not.spproj%os_cls2D%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(icls,1) = nint(spproj%os_cls2D%get(i,'prev_pop_even'))
                prev_eo_pops(icls,2) = nint(spproj%os_cls2D%get(i,'prev_pop_odd'))
            enddo
        endif
    end subroutine cavger_transf_oridat

    subroutine cavger_read_euclid_sigma2
        character(len=STDLEN) :: fname
        if( l_ml_reg )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            call eucl_sigma%read_part(  build_glob%spproj_field, pptcl_mask)
            call eucl_sigma%read_groups(build_glob%spproj_field, pptcl_mask)
        end if
    end subroutine cavger_read_euclid_sigma2

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         poulation, average correlation and weight
    subroutine cavger_gen2Dclassdoc( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        integer  :: pops(params_glob%ncls)
        real(dp) :: corrs(params_glob%ncls), ws(params_glob%ncls), specscores(params_glob%ncls)
        real     :: frc05, frc0143, rstate
        integer  :: i, iptcl, icls, pop, nptcls
        nptcls     = spproj%os_ptcl2D%get_noris()
        pops       = 0
        corrs      = 0.d0
        ws         = 0.d0
        specscores = 0.d0
        !$omp parallel do default(shared) private(iptcl,rstate,icls) schedule(static)&
        !$omp proc_bind(close) reduction(+:pops,corrs,ws,specscores)
        do iptcl=1,nptcls
            rstate = spproj%os_ptcl2D%get(iptcl,'state')
            if( rstate < 0.5 )cycle
            icls = nint(spproj%os_ptcl2D%get(iptcl,'class'))
            if( icls<1 .or. icls>params_glob%ncls )cycle
            pops(icls)       = pops(icls)      + 1
            corrs(icls)      = corrs(icls)     + spproj%os_ptcl2D%get(iptcl,'corr')
            ws(icls)         = ws(icls)        + spproj%os_ptcl2D%get(iptcl,'w')
            specscores(icls) = specscores(icls)+ spproj%os_ptcl2D%get(iptcl,'specscore')
        enddo
        !$omp end parallel do
        if( l_stream  .and. spproj%os_cls2D%get_noris()==ncls .and. params_glob%update_frac<.99 )then
            do i = 1,ncls
                icls = nint(spproj%os_cls2D%get(i,'class'))
                if( .not.spproj%os_cls2D%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(icls,1) = nint(spproj%os_cls2D%get(i,'prev_pop_even'))
                prev_eo_pops(icls,2) = nint(spproj%os_cls2D%get(i,'prev_pop_odd'))
                pop = sum(prev_eo_pops(icls,:))
                if( pop == 0 ) cycle
                corrs(icls)      = corrs(icls)      + real(pop) * spproj%os_cls2D%get(i,'corr')
                ws(icls)         = ws(icls)         + real(pop) * spproj%os_cls2D%get(i,'w')
                specscores(icls) = specscores(icls) + real(pop) * spproj%os_cls2D%get(i,'specscore')
                pops(icls)       = pops(icls) + pop
            enddo
        endif
        where(pops>1)
            corrs      = corrs / real(pops)
            ws         = ws / real(pops)
            specscores = specscores / real(pops)
        elsewhere
            corrs      = -1.
            ws         = 0.
            specscores = 0.
        end where
        call spproj%os_cls2D%new(params_glob%ncls, is_ptcl=.false.)
        do icls=1,params_glob%ncls
            pop = pops(icls)
            call build_glob%clsfrcs%estimate_res(icls, frc05, frc0143)
            call spproj%os_cls2D%set(icls, 'class',     real(icls))
            call spproj%os_cls2D%set(icls, 'pop',       real(pop))
            call spproj%os_cls2D%set(icls, 'res',       frc0143)
            call spproj%os_cls2D%set(icls, 'corr',      real(corrs(icls)))
            call spproj%os_cls2D%set(icls, 'w',         real(ws(icls)))
            call spproj%os_cls2D%set(icls, 'specscore', real(specscores(icls)))
            if( pop > 0 )then
                call spproj%os_cls2D%set(icls, 'state', 1.0) ! needs to be default val if no selection has been done
            else
                call spproj%os_cls2D%set(icls, 'state', 0.0) ! exclusion
            endif
        end do
    end subroutine cavger_gen2Dclassdoc

    !>  \brief  is for initialization of the sums
    subroutine init_cavgs_sums
        integer :: icls
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)%new(ldim_croppd,smpd_crop,wthreads=.false.)
            call cavgs_odd(icls)%new(ldim_croppd,smpd_crop,wthreads=.false.)
            call cavgs_merged(icls)%new(ldim_croppd,smpd_crop,wthreads=.false.)
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
        logical, intent(in)           :: do_frac_update
        integer, parameter            :: READBUFFSZ = 1024
        complex, parameter            :: zero = cmplx(0.,0.)
        type(kbinterpol)              :: kbwin
        type(stack_io)                :: stkio_r
        type(image_ptr)               :: pcmat(nthr_glob), prhomat(nthr_glob)
        type(image),      allocatable :: cgrid_imgs(:), read_imgs(:), cgrid_imgs_crop(:)
        character(len=:), allocatable :: stk_fname
        complex,          allocatable :: cmats(:,:,:)
        real,             allocatable :: rhos(:,:,:), tvals(:,:,:)
        complex :: fcompl, fcompll
        real    :: loc(2), mat(2,2), dist(2), add_phshift, tval, kw, maxspafreqsq, reg_scale, crop_scale
        integer :: batch_iprecs(READBUFFSZ), fdims_crop(3), logi_lims_crop(3,2)
        integer :: phys(2), win_corner(2), cyc_lims_cropR(2,2),cyc_lims_crop(3,2), sigma2_kfromto(2)
        integer :: iprec, i, sh, iwinsz, nyq_crop, ind_in_stk, iprec_glob, nptcls_eff, radfirstpeak
        integer :: wdim, h, k, l, m, ll, mm, icls, iptcl, interp_shlim, batchind
        integer :: first_stkind, fromp, top, istk, nptcls_in_stk, nstks, last_stkind
        integer :: ibatch, nbatches, istart, iend, ithr, nptcls_in_batch, first_pind, last_pind
        if( .not. params_glob%l_distr_exec ) write(logfhandle,'(a)') '>>> ASSEMBLING CLASS SUMS'
        ! init cavgs
        call init_cavgs_sums
        if( do_frac_update )then
            call cavger_readwrite_partial_sums( 'read' )
            call cavger_apply_weights( 1. - params_glob%update_frac )
        endif
        kbwin  = kbinterpol(KBWINSZ, params_glob%alpha)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
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
        !$omp parallel default(shared) proc_bind(close) private(i)
        !$omp do schedule(static)
        do i = 1,READBUFFSZ
            call read_imgs(i)%new(ldim, params_glob%smpd, wthreads=.false.)
            call cgrid_imgs(i)%new(ldim_pd, params_glob%smpd, wthreads=.false.)
            call cgrid_imgs_crop(i)%new(ldim_croppd, params_glob%smpd_crop, wthreads=.false.)
        enddo
        !$omp end do
        !$omp end parallel
        logi_lims_crop = cgrid_imgs_crop(1)%loop_lims(2)
        cyc_lims_crop  = cgrid_imgs_crop(1)%loop_lims(3)
        nyq_crop       = cgrid_imgs_crop(1)%get_lfny(1)
        fdims_crop     = cgrid_imgs_crop(1)%get_array_shape()
        cyc_lims_cropR(:,1) = cyc_lims_crop(1,:)
        cyc_lims_cropR(:,2) = cyc_lims_crop(2,:)
        allocate(tvals(fdims_crop(1),fdims_crop(2),READBUFFSZ),cmats(fdims_crop(1),fdims_crop(2),READBUFFSZ),&
        &rhos(fdims_crop(1),fdims_crop(2),READBUFFSZ))
        interp_shlim = nyq_crop + 1
        ! Main loop
        iprec_glob = 0 ! global record index
        do istk = first_stkind,last_stkind
            ! Particles range in stack
            fromp  = nint(build_glob%spproj%os_stk%get(istk,'fromp'))
            top    = nint(build_glob%spproj%os_stk%get(istk,'top'))
            nptcls_in_stk = top - fromp + 1 ! # of particles in stack
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, max(params_glob%fromp,fromp), stk_fname, ind_in_stk)
            ! open buffer
            call stkio_r%open(stk_fname, smpd, 'read', bufsz=min(nptcls_in_stk,READBUFFSZ))
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
                    call stkio_r%read(precs(iprec_glob)%ind_in_stk, read_imgs(batchind)) ! read
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
            ! close stack
            call stkio_r%close
        enddo
        ! Cleanup
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
            allocate(cavgs_even_bak(ncls),cavgs_odd_bak(ncls),&
            &ctfsqsums_even_bak(ncls),ctfsqsums_odd_bak(ncls))
            if( l_stream )then
                allocate(cavgs_even_wfilt_bak(ncls),cavgs_odd_wfilt_bak(ncls),&
                &ctfsqsums_even_wfilt_bak(ncls),ctfsqsums_odd_wfilt_bak(ncls))
            endif
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
                    call cavgs_even_bak(icls)%copy(cavgs_even(icls))
                    call cavgs_odd_bak(icls)%copy(cavgs_odd(icls))
                    call ctfsqsums_even_bak(icls)%copy(ctfsqsums_even(icls))
                    call ctfsqsums_odd_bak(icls)%copy(ctfsqsums_odd(icls))
                    if( l_stream )then
                        call cavgs_even_wfilt_bak(icls)%copy(cavgs_even_wfilt(icls))
                        call cavgs_odd_wfilt_bak(icls)%copy(cavgs_odd_wfilt(icls))
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
        use simple_masker, only: automask2D
        character(len=*), intent(in) :: fname
        integer,          intent(in) :: which_iter
        type(image), allocatable     :: even_imgs(:), odd_imgs(:)
        real,        allocatable     :: frc(:)
        integer :: eo_pop(2), icls, find, find_plate, pop, filtsz_crop
        filtsz_crop = cavgs_even(1)%get_filtsz()
        allocate(even_imgs(ncls), odd_imgs(ncls), frc(filtsz_crop))
        do icls=1,ncls
            call even_imgs(icls)%copy(cavgs_even(icls))
            call odd_imgs(icls)%copy(cavgs_odd(icls))
        end do
        if( l_ml_reg )then
            !$omp parallel do default(shared) private(icls,frc,find,find_plate,pop,eo_pop) schedule(static) proc_bind(close)
            do icls=1,ncls
                eo_pop = prev_eo_pops(icls,:) + eo_class_pop(icls)
                pop    = sum(eo_pop)
                if( pop > 0 )then
                    call even_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                    call odd_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                    call even_imgs(icls)%fft()
                    call odd_imgs(icls)%fft()
                    call even_imgs(icls)%fsc(odd_imgs(icls), frc)
                    find_plate = 0
                    if( phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                    find = build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1)
                    find = max(find, find_plate)
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
                    call cavgs_even(icls)%copy(cavgs_even_bak(icls))
                    call cavgs_odd(icls)%copy(cavgs_odd_bak(icls))
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
            !$omp parallel do default(shared) private(icls,frc,find,find_plate) schedule(static) proc_bind(close)
            do icls=1,ncls
                call even_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                call odd_imgs(icls)%mask(params_glob%msk_crop, 'soft')
                call even_imgs(icls)%fft()
                call odd_imgs(icls)%fft()
                call even_imgs(icls)%fsc(odd_imgs(icls), frc)
                find_plate = 0
                if( phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
                call build_glob%clsfrcs%set_frc(icls, frc, 1)
                ! average low-resolution info between eo pairs to keep things in register
                find = build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1)
                find = max(find, find_plate)
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

    ! private function to add noise term to denomnator
    subroutine add_invtausq2rho( ctfsqsum, frc )
        class(image),          intent(inout) :: ctfsqsum
        real,     allocatable, intent(in)    :: frc(:)
        real,     allocatable :: sig2(:), tau2(:), ssnr(:)
        integer,  allocatable :: cnt(:)
        real(dp), allocatable :: rsum(:)
        complex,      pointer :: pctfsqsum(:,:,:)
        real,     parameter   :: fudge = 1.0
        real    :: cc, scale, pad_factor, invtau2
        integer :: flims(3,2), phys(2), h, k, sh, sz, reslim_ind
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

    ! I/O

    !>  \brief  writes class averages to disk
    subroutine cavger_write( fname, which )
        character(len=*),  intent(in) :: fname, which
        character(len=:), allocatable :: fname_wfilt
        integer :: icls
        fname_wfilt = add2fbody(fname, params_glob%ext, trim(WFILT_SUFFIX))
        select case(which)
            case('even')
                if( l_stream ) fname_wfilt = add2fbody(fname, '_even'//trim(params_glob%ext), trim(WFILT_SUFFIX))
                do icls=1,ncls
                    call cavgs_even(icls)%write(fname, icls)
                    if( l_stream ) call cavgs_even_wfilt(icls)%write(fname_wfilt, icls)
                end do
            case('odd')
                if( l_stream ) fname_wfilt = add2fbody(fname, '_odd'//trim(params_glob%ext), trim(WFILT_SUFFIX))
                do icls=1,ncls
                    call cavgs_odd(icls)%write(fname, icls)
                    if( l_stream ) call cavgs_odd_wfilt(icls)%write(fname_wfilt, icls)
                end do
            case('merged')
                if( l_stream ) fname_wfilt = add2fbody(fname, params_glob%ext, trim(WFILT_SUFFIX))
                do icls=1,ncls
                    call cavgs_merged(icls)%write(fname, icls)
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

    !>  \brief  reads class averages from disk
    subroutine cavger_read( fname, which )
        character(len=*),  intent(in) :: fname, which
        type(stack_io) :: stkio_r
        integer        :: icls
        select case(which)
            case('even')
                call stkio_r%open(trim(fname), smpd_crop, 'read', bufsz=ncls)
                do icls=1,ncls
                    call cavgs_even(icls)%new(ldim_crop,smpd_crop,wthreads=.false.)
                    call stkio_r%read(icls, cavgs_even(icls))
                end do
                call stkio_r%close
            case('odd')
                call stkio_r%open(trim(fname), smpd_crop, 'read', bufsz=ncls)
                do icls=1,ncls
                    call cavgs_odd(icls)%new(ldim_crop,smpd_crop,wthreads=.false.)
                    call stkio_r%read(icls, cavgs_odd(icls))
                end do
                call stkio_r%close
            case('merged')
                call stkio_r%open(trim(fname), smpd_crop, 'read', bufsz=ncls)
                do icls=1,ncls
                    call cavgs_merged(icls)%new(ldim_crop,smpd_crop,wthreads=.false.)
                    call stkio_r%read(icls, cavgs_merged(icls))
                end do
                call stkio_r%close
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine cavger_read

    !>  \brief  writes partial class averages to disk (distributed execution)
    subroutine cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        integer                       :: icls, ldim_here(3)
        character(len=:), allocatable :: cae, cao, cte, cto, caewf, caowf, ctewf, ctowf
        type(stack_io)                :: stkio(4)
        logical                       :: is_ft
        allocate(cae, source='cavgs_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(cao, source='cavgs_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(cte, source='ctfsqsums_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(caewf, source='cavgs_even_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(caowf, source='cavgs_odd_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(ctewf, source='ctfsqsums_even_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(ctowf, source='ctfsqsums_odd_wfilt_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
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
                is_ft = cavgs_even(1)%is_ft()
                ldim_here  = cavgs_even(1)%get_ldim()
                call stkio(1)%open(cae, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(2)%open(cao, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(3)%open(cte, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(4)%open(cto, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                do icls=1,ncls
                    call stkio(1)%write(icls, cavgs_even(icls))
                    call stkio(2)%write(icls, cavgs_odd(icls))
                    call stkio(3)%write(icls, ctfsqsums_even(icls))
                    call stkio(4)%write(icls, ctfsqsums_odd(icls))
                end do
                if( l_stream )then
                    call stkio(1)%open(caewf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    call stkio(2)%open(caowf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    call stkio(3)%open(ctewf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    call stkio(4)%open(ctowf, smpd_crop, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                    do icls=1,ncls
                        call stkio(1)%write(icls, cavgs_even_wfilt(icls))
                        call stkio(2)%write(icls, cavgs_odd_wfilt(icls))
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
        deallocate(cae, cao, cte, cto)
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
        complex(kind=c_float_complex), pointer :: cmat_ptr1(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr2(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr3(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr4(:,:,:) => null()
        integer(timer_int_kind)       ::  t_init,  t_io,  t_workshare_sum,  t_set_sums,  t_merge_eos_and_norm,  t_tot
        real(timer_int_kind)          :: rt_init, rt_io, rt_workshare_sum, rt_set_sums, rt_merge_eos_and_norm, rt_tot
        complex,          allocatable :: csums(:,:,:,:)
        character(len=:), allocatable :: cae, cao, cte, cto
        character(len=STDLEN)         :: benchfname
        type(image) :: imgs4read(4)
        integer     :: ipart, icls, array_shape(3), ldim_here(3), fnr
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
            allocate(cae, source='cavgs_even_part'    //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            allocate(cao, source='cavgs_odd_part'     //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            allocate(cte, source='ctfsqsums_even_part'//int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            allocate(cto, source='ctfsqsums_odd_part' //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            ! serial read
            call imgs4read(1)%read(cae)
            call imgs4read(2)%read(cao)
            call imgs4read(3)%read(cte)
            call imgs4read(4)%read(cto)
            deallocate(cae, cao, cte, cto)
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
                allocate(cae, source='cavgs_even_wfilt_part'    //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
                allocate(cao, source='cavgs_odd_wfilt_part'     //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
                allocate(cte, source='ctfsqsums_even_wfilt_part'//int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
                allocate(cto, source='ctfsqsums_odd_wfilt_part' //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
                ! serial read
                call imgs4read(1)%read(cae)
                call imgs4read(2)%read(cao)
                call imgs4read(3)%read(cte)
                call imgs4read(4)%read(cto)
                deallocate(cae, cao, cte, cto)
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
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
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
    subroutine cavger_kill
        integer :: icls
        if( exists )then
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
            deallocate( cavgs_even, cavgs_odd, cavgs_merged, ctfsqsums_even,&
            &ctfsqsums_odd, ctfsqsums_merged, pptcl_mask, prev_eo_pops)
            if( allocated(cavgs_even_bak) )then
                do icls=1,ncls
                    call cavgs_even_bak(icls)%kill
                    call cavgs_odd_bak(icls)%kill
                    call ctfsqsums_even_bak(icls)%kill
                    call ctfsqsums_odd_bak(icls)%kill
                    if( l_stream )then
                        call cavgs_even_wfilt_bak(icls)%kill
                        call cavgs_odd_wfilt_bak(icls)%kill
                        call ctfsqsums_even_wfilt_bak(icls)%kill
                        call ctfsqsums_odd_wfilt_bak(icls)%kill
                    endif
                enddo
                deallocate(cavgs_even_bak,cavgs_odd_bak, ctfsqsums_even_bak,ctfsqsums_odd_bak)
                if( l_stream ) deallocate(cavgs_even_wfilt_bak,cavgs_odd_wfilt_bak,ctfsqsums_even_wfilt_bak,ctfsqsums_odd_wfilt_bak)
            endif
            deallocate(precs)
            istart = 0
            iend   = 0
            partsz = 0
            ncls   = 0
            exists = .false.
        endif
    end subroutine cavger_kill

end module simple_classaverager
