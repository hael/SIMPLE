!@descr: Routines to perform the classes restoration and processing
submodule (simple_new_classaverager) simple_new_classaverager_restore
use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
implicit none
#include "simple_local_flags.inc"

contains

    !>  \brief  Constructor
    module subroutine cavger_new( pinds, alloccavgs )
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
        ! smpd
        smpd          = params_glob%smpd
        smpd_crop     = params_glob%smpd_crop
        ! set ldims
        ldim          = [params_glob%box,  params_glob%box,  1]
        ldim_crop     = [params_glob%box_crop,  params_glob%box_crop,  1]
        ldim_croppd   = [params_glob%box_croppd,params_glob%box_croppd,1]
        ldim_pd       = [params_glob%boxpd,params_glob%boxpd,1]
        ! build arrays
        allocate(precs(partsz))
        if( l_alloc_read_cavgs ) call cavgs%new_set(ldim_crop(1:2), ncls)
    end subroutine cavger_new

    ! setters/getters

    !>  \brief  transfers metadata to the instance
    module subroutine cavger_transf_oridat( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        class(oris), pointer :: spproj_field
        type(ctfparams)      :: ctfvars(nthr_glob)
        integer              :: cnt, iptcl, ithr, stkind
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
            precs(cnt)%class   = spproj_field%get_class(iptcl)
            precs(cnt)%e3      = spproj_field%e3get(iptcl)
            precs(cnt)%shift   = spproj_field%get_2Dshift(iptcl)
            call spproj%map_ptcl_ind2stk_ind(params_glob%oritype, iptcl, stkind, precs(cnt)%ind_in_stk)
        end do
        !$omp end parallel do
    end subroutine cavger_transf_oridat

    !>  \brief  for loading sigma2
    module subroutine cavger_read_euclid_sigma2
        type(string) :: fname
        if( params_glob%l_ml_reg )then
            fname = SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'
            call eucl_sigma%new(fname, params_glob%box)
            call eucl_sigma%read_part(  build_glob%spproj_field)
            call eucl_sigma%read_groups(build_glob%spproj_field)
        end if
    end subroutine cavger_read_euclid_sigma2

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population, average correlation and weight
    module subroutine cavger_gen2Dclassdoc( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), target, intent(inout) :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(params_glob%ncls)
        real(dp) :: corrs(params_glob%ncls), ws(params_glob%ncls)
        real     :: frc05, frc0143, rstate, w
        integer  :: iptcl, icls, pop, nptcls
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
    end subroutine cavger_gen2Dclassdoc

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
    !           using convolution interpolation in Fourier space
    module subroutine cavger_assemble_sums_conv( do_frac_update )
        use simple_strategy2D3D_common, only: prepimgbatch
        use simple_math_ft,             only: upsample_sigma2
        logical, intent(in)      :: do_frac_update
        real,          parameter :: KB2        = KBALPHA**2
        integer,       parameter :: READBUFFSZ = 1024
        type(image), allocatable :: tmp_pad_imgs(:)
        type(kbinterpol)         :: kbwin
        type(dstack_io)          :: dstkio_r
        type(string)             :: stk_fname
        complex(kind=c_float_complex), allocatable :: cmats(:,:,:), interp_cmats(:,:,:)
        complex(kind=c_float_complex)              :: fcomp, fcompl
        real,                          allocatable :: tvals(:,:,:), kbw(:,:,:), interp_rhos(:,:,:)
        real,                          allocatable :: sigma2(:,:), sigma2pd(:,:), sqrt_sigma2pd(:,:)
        real    :: loc(2), mat(2,2), tval, croppd_scale, w
        integer :: batch_iprecs(READBUFFSZ), fdims_croppd(3), win(2,2), flims_crop(3,2), phys(2)
        integer :: cyc_lims_croppdR(2,2),cyc_lims_croppd(3,2), sigma2_kfromto(2), cshape_crop(2)
        integer :: iprec, iptcl, i, sh, ind_in_stk, iprec_glob, nptcls_eff, iwinsz
        integer :: wdim, h,k,hh,kk, l, m, icls, nyq_crop, nyq_croppd, batchind, physh, physk
        integer :: first_stkind, fromp, top, istk, nptcls_in_stk, nstks, last_stkind
        integer :: ibatch, nbatches, istart, iend, ithr, nptcls_in_batch, first_pind, last_pind
        logical :: l_conjg
        if( .not. params_glob%l_distr_exec ) write(logfhandle,'(a)') '>>> ASSEMBLING CLASS SUMS'
        ! Init cavgs
        if( l_alloc_read_cavgs )then
            call cavgs%zero_set(.true.)
            if( do_frac_update )then
                call cavger_readwrite_partial_sums('read')
                call apply_weights2cavgs(1.0-params_glob%update_frac)
            endif
        else
            if( do_frac_update ) call apply_weights2cavgs(1.0-params_glob%update_frac)
            call cavgs%zero_set(.true.)
        endif
        ! Interpolation variables
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        allocate(kbw(wdim,wdim,nthr_glob),source=0.)
        ! Indexing the number stacks
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
        ! Memoization, dimensions & limits
        call memoize_ft_maps(ldim_croppd(1:2), params_glob%smpd_crop)
        croppd_scale          = real(ldim_croppd(1)) / real(ldim_pd(1))
        cyc_lims_croppd       = ft_map_lims_nr            ! croppd limits
        cyc_lims_croppdR(:,1) = cyc_lims_croppd(1,:)      ! croppd transposed limits
        cyc_lims_croppdR(:,2) = cyc_lims_croppd(2,:)
        fdims_croppd          = ft_map_get_farray_shape() ! croppd dimensions
        nyq_croppd            = ft_map_get_lfny()         ! croppd Nyquist limit
        cshape_crop = cavgs%even%cshape                   ! crop shape of complex physical matrix
        flims_crop  = cavgs%even%flims                    ! crop Fourier limits of complex matrix
        nyq_crop    = cavgs%even%fit%get_lfny(1)          ! crop Nyquist limit
        ! Work images
        call prepimgbatch(READBUFFSZ)
        call alloc_imgarr(nthr_glob, ldim_pd, smpd, tmp_pad_imgs)
        ! sigma2 prep for regularization
        if( params_glob%l_ml_reg )then
            allocate(sigma2(1:nyq_crop,nthr_glob),sigma2pd(0:nyq_croppd,nthr_glob),&
                    &sqrt_sigma2pd(0:nyq_croppd,nthr_glob),source=0.0)
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
        endif
        ! Work arrays
        allocate(tvals(fdims_croppd(1),fdims_croppd(2),nthr_glob),&
                &cmats(fdims_croppd(1),fdims_croppd(2),nthr_glob),&
                &interp_cmats(cshape_crop(1),cshape_crop(2),READBUFFSZ),&
                &interp_rhos(cshape_crop(1),cshape_crop(2),READBUFFSZ))
        ! Fast particle reading
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
            nbatches = ceiling(real(nptcls_in_stk)/real(READBUFFSZ))
            do ibatch = 1,nbatches
                batch_iprecs = 0                                     ! records in batch, if zero skip
                istart = (ibatch - 1)              * READBUFFSZ + 1  ! first index in current batch
                iend   = min(nptcls_in_stk, istart + READBUFFSZ - 1) ! last  index in current batch
                nptcls_in_batch = iend-istart+1
                batchind   = 0
                nptcls_eff = 0                                       ! # particles to process in batch
                ! Read
                do i = istart,iend
                    iptcl    = fromp + i - 1                         ! global particle index
                    batchind = batchind + 1                          ! index in batch
                    if( iptcl < params_glob%fromp ) cycle            ! taking care of limits
                    if( iptcl > params_glob%top )   cycle
                    iprec_glob = iprec_glob + 1                      ! global particle record
                    batch_iprecs(batchind) = iprec_glob              ! particle record in batch
                    if( precs(iprec_glob)%pind == 0 ) cycle
                    nptcls_eff = nptcls_eff + 1
                    call dstkio_r%read(stk_fname, precs(iprec_glob)%ind_in_stk, build_glob%imgbatch(batchind))
                enddo
                if( nptcls_eff == 0 ) cycle
                ! Interpolation loop
                !$omp parallel default(shared) proc_bind(close)&
                !$omp private(i,ithr,icls,iprec,iptcl,win,mat,h,k,hh,kk,l,m,loc,sh,physh,physk,tval)&
                !$omp private(fcomp,fcompl,l_conjg,phys,w)
                !$omp do schedule(static)
                do i = 1,nptcls_in_batch
                    iprec = batch_iprecs(i)
                    if( iprec == 0 ) cycle
                    iptcl = precs(iprec)%pind
                    if( iptcl == 0 ) cycle
                    ithr  = omp_get_thread_num() + 1
                    cmats(:,:,ithr) = CMPLX_ZERO
                    tvals(:,:,ithr) = 0.0
                    ! prep image: noise normalization, edge tappering, padding, shift, copy into cmats
                    call build_glob%imgbatch(i)%norm_noise_taper_edge_pad_fft_shift_2mat(build_glob%lmsk, tmp_pad_imgs(ithr),&
                        &-precs(iprec)%shift*croppd_scale, ldim_croppd, ft_map_lims, cmats(:,:,ithr))
                    ! apply CTF to padded image, stores CTF values
                    call precs(iprec)%tfun%eval_and_apply(ctfflag, precs(iprec)%dfx, precs(iprec)%dfy,&
                        &precs(iprec)%angast, fdims_croppd(1:2), tvals(:,:,ithr), cmats(:,:,ithr))
                    ! upsample sigma2 & multipliy CTF.image & CTF2
                    if( params_glob%l_ml_reg )then
                        sigma2(sigma2_kfromto(1):nyq_crop,ithr) = eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):nyq_crop,iptcl)
                        call upsample_sigma2(sigma2_kfromto(1), nyq_crop, sigma2(:,ithr), nyq_croppd, sigma2pd(:,ithr))
                        sqrt_sigma2pd(:,ithr) = sqrt(sigma2pd(:,ithr))
                        do h = ft_map_lims(1,1),ft_map_lims(1,2)
                            do k = ft_map_lims(2,1),ft_map_lims(2,2)
                                sh = nint(hyp(h,k))
                                if( sh > nyq_croppd )cycle
                                physh = ft_map_phys_addrh(h,k)
                                physk = ft_map_phys_addrk(h,k)
                                cmats(physh,physk,ithr) = cmats(physh,physk,ithr) /      sigma2pd(sh,ithr)
                                tvals(physh,physk,ithr) = tvals(physh,physk,ithr) / sqrt_sigma2pd(sh,ithr)
                            enddo
                        enddo
                    end if
                    ! Rotation matrix
                    call rotmat2d(-precs(iprec)%e3, mat)
                    ! scale the matrix to map to the padded image
                    mat = KBALPHA * mat
                    ! Kaiser-Bessel Interpolation
                    interp_cmats(:,:,i) = CMPLX_ZERO
                    interp_rhos(:,:,i)  = 0.0
                    ! loop over cropped original image limits
                    do h = flims_crop(1,1), flims_crop(1,2)
                        do k = flims_crop(2,1), flims_crop(2,2)
                            sh = nint(hyp(real(h),real(k)))
                            if( sh > nyq_crop )cycle
                            ! rotation in padded coordinates
                            loc = matmul(real([h,k]),mat)
                            ! interpolation padded window
                            win(1,:) = nint(loc)
                            win(2,:) = win(1,:) + iwinsz
                            win(1,:) = win(1,:) - iwinsz
                            ! interpolation kernel
                            call kbwin%apod_mat_2d(loc, iwinsz, wdim, kbw(:,:,ithr))
                            ! interpolation
                            fcomp = CMPLX_ZERO
                            tval  = 0.0
                            do l = 1, wdim
                                hh      = win(1,1)+ l-1
                                l_conjg = hh < 0
                                hh      = cyci_1d(cyc_lims_croppdR(:,1), hh)
                                fcompl  = CMPLX_ZERO
                                do m = 1, wdim
                                    kk     = win(1,2) +  m-1
                                    kk     = cyci_1d(cyc_lims_croppdR(:,2), kk)
                                    physh  = ft_map_phys_addrh(hh,kk)
                                    physk  = ft_map_phys_addrk(hh,kk)
                                    w      = kbw(l,m,ithr)
                                    fcompl = fcompl + w * cmats(physh, physk, ithr)
                                    tval   = tval   + w * tvals(physh, physk, ithr)
                                enddo
                                fcomp = fcomp + merge(conjg(fcompl), fcompl, l_conjg)
                            end do
                            ! physical address with original dimension
                            phys = cavgs%even%fit%comp_addr_phys(h,k)
                            ! weighting of complex value by particle & padding correction
                            interp_cmats(phys(1), phys(2), i) = precs(iprec)%pw * fcomp * KB2
                            ! weighting of CTF^2 by particle
                            interp_rhos( phys(1), phys(2), i) = precs(iprec)%pw * tval*tval
                        end do
                    end do
                enddo
                !$omp end do
                ! Accumulate sums
                !$omp do schedule(static)
                do icls = 1,ncls
                    do i = 1,nptcls_in_batch
                        iprec = batch_iprecs(i)
                        if( iprec == 0 ) cycle
                        if( precs(iprec)%pind == 0 ) cycle
                        if( precs(iprec)%class == icls )then
                            select case(precs(iprec)%eo)
                            case(0,-1)
                                cavgs%even%cmat(:,:,icls)  = cavgs%even%cmat(:,:,icls)  + interp_cmats(:,:,i)
                                cavgs%even%ctfsq(:,:,icls) = cavgs%even%ctfsq(:,:,icls) + interp_rhos(:,:,i)
                            case(1)
                                cavgs%odd%cmat(:,:,icls)  = cavgs%odd%cmat(:,:,icls)  + interp_cmats(:,:,i)
                                cavgs%odd%ctfsq(:,:,icls) = cavgs%odd%ctfsq(:,:,icls) + interp_rhos(:,:,i)
                            end select
                        endif
                    enddo
                enddo
                !$omp end do
                !$omp end parallel
            enddo ! end read batches loop
        enddo
        ! Cleanup
        call dstkio_r%kill
        call forget_ft_maps
        call dealloc_imgarr(tmp_pad_imgs)
        deallocate(cmats,interp_cmats,tvals,kbw,interp_rhos)
    end subroutine cavger_assemble_sums_conv

    !>  \brief  is for assembling the sums in distributed/non-distributed mode
    !           using interpolation in Fourier space
    module subroutine cavger_assemble_sums( do_frac_update )
        use simple_strategy2D3D_common, only: prepimgbatch
        use simple_math_ft,             only: upsample_sigma2
        logical, intent(in)      :: do_frac_update
        real,          parameter :: KB2        = KBALPHA**2
        integer,       parameter :: READBUFFSZ = 1024
        type(image), allocatable :: tmp_pad_imgs(:)
        type(kbinterpol)         :: kbwin
        type(dstack_io)          :: dstkio_r
        type(string)             :: stk_fname
        complex(kind=c_float_complex), allocatable :: cmats(:,:,:), interp_cmats(:,:,:)
        complex(kind=c_float_complex)              :: fcomp
        real,                          allocatable :: tvals(:,:,:), kbw(:,:,:), interp_rhos(:,:,:)
        real,                          allocatable :: sigma2(:,:), sigma2pd(:,:), sqrt_sigma2pd(:,:)
        real    :: loc(2), mat(2,2), tvalsq, croppd_scale, w
        integer :: batch_iprecs(READBUFFSZ), fdims_croppd(3), win(2,2), flims_crop(3,2), phys(2)
        integer :: cyc_lims_croppdR(2,2),cyc_lims_croppd(3,2), sigma2_kfromto(2), cshape_crop(2)
        integer :: cyc_lims_cropR(2,2),cyc_lims_crop(3,2)
        integer :: iprec, iptcl, i, sh, ind_in_stk, iprec_glob, nptcls_eff, iwinsz
        integer :: wdim, h,k,hh,kk,hp,kp, l, m, icls, nyq_crop, nyq_croppd, batchind, physh, physk
        integer :: first_stkind, fromp, top, istk, nptcls_in_stk, nstks, last_stkind
        integer :: ibatch, nbatches, istart, iend, ithr, nptcls_in_batch, first_pind, last_pind
        logical :: l_conjg
        if( .not. params_glob%l_distr_exec ) write(logfhandle,'(a)') '>>> ASSEMBLING CLASS SUMS'
        ! Init cavgs
        if( l_alloc_read_cavgs )then
            call cavgs%zero_set(.true.)
            if( do_frac_update )then
                call cavger_readwrite_partial_sums('read')
                call apply_weights2cavgs(1.0-params_glob%update_frac)
            endif
        else
            if( do_frac_update ) call apply_weights2cavgs(1.0-params_glob%update_frac)
            call cavgs%zero_set(.true.)
        endif
        ! Interpolation variables
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        allocate(kbw(wdim,wdim,nthr_glob),source=0.)
        ! Indexing the number stacks
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
        ! Work images
        call prepimgbatch(READBUFFSZ)
        call alloc_imgarr(nthr_glob, ldim_pd, smpd, tmp_pad_imgs)
        ! Memoization for padded image
        call memoize_ft_maps(ldim_croppd(1:2), params_glob%smpd_crop)
        ! Dimensions & limits
        croppd_scale          = real(ldim_croppd(1)) / real(ldim_pd(1))
        cyc_lims_croppd       = ft_map_lims_nr            ! croppd limits
        cyc_lims_croppdR(:,1) = cyc_lims_croppd(1,:)      ! croppd transposed limits
        cyc_lims_croppdR(:,2) = cyc_lims_croppd(2,:)
        fdims_croppd          = ft_map_get_farray_shape() ! croppd dimensions
        nyq_croppd            = ft_map_get_lfny()         ! croppd Nyquist limit
        cshape_crop    = cavgs%even%cshape                   ! crop shape of complex physical matrix
        flims_crop     = cavgs%even%flims                    ! crop Fourier limits of complex matrix
        nyq_crop       = cavgs%even%fit%get_lfny(1)          ! crop Nyquist limit
        cyc_lims_crop  = cavgs%even%fit%loop_lims(3)
        cyc_lims_cropR = transpose(cyc_lims_crop(1:2,:))
        ! sigma2 prep for regularization
        if( params_glob%l_ml_reg )then
            allocate(sigma2(1:nyq_crop,nthr_glob),sigma2pd(0:nyq_croppd,nthr_glob),&
                    &sqrt_sigma2pd(0:nyq_croppd,nthr_glob),source=0.0)
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
        endif
        ! Work arrays
        allocate(tvals(fdims_croppd(1),fdims_croppd(2),nthr_glob),&
                &cmats(fdims_croppd(1),fdims_croppd(2),nthr_glob),&
                &interp_cmats(cshape_crop(1),cshape_crop(2),READBUFFSZ),&
                &interp_rhos(cshape_crop(1),cshape_crop(2),READBUFFSZ))
        ! Fast particle reading
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
            nbatches = ceiling(real(nptcls_in_stk)/real(READBUFFSZ))
            do ibatch = 1,nbatches
                batch_iprecs = 0                                     ! records in batch, if zero skip
                istart = (ibatch - 1)              * READBUFFSZ + 1  ! first index in current batch
                iend   = min(nptcls_in_stk, istart + READBUFFSZ - 1) ! last  index in current batch
                nptcls_in_batch = iend-istart+1
                batchind   = 0
                nptcls_eff = 0                                       ! # particles to process in batch
                ! Read
                do i = istart,iend
                    iptcl    = fromp + i - 1                         ! global particle index
                    batchind = batchind + 1                          ! index in batch
                    if( iptcl < params_glob%fromp ) cycle            ! taking care of limits
                    if( iptcl > params_glob%top )   cycle
                    iprec_glob = iprec_glob + 1                      ! global particle record
                    batch_iprecs(batchind) = iprec_glob              ! particle record in batch
                    if( precs(iprec_glob)%pind == 0 ) cycle
                    nptcls_eff = nptcls_eff + 1
                    call dstkio_r%read(stk_fname, precs(iprec_glob)%ind_in_stk, build_glob%imgbatch(batchind))
                enddo
                if( nptcls_eff == 0 ) cycle
                ! Interpolation loop
                !$omp parallel default(shared) proc_bind(close)&
                !$omp private(i,ithr,icls,iprec,iptcl,win,mat,h,k,hh,kk,l,m,loc,sh,physh,physk)&
                !$omp private(fcomp,phys,w,hp,kp,tvalsq,l_conjg)
                !$omp do schedule(static)
                do i = 1,nptcls_in_batch
                    iprec = batch_iprecs(i)
                    if( iprec == 0 ) cycle
                    iptcl = precs(iprec)%pind
                    if( iptcl == 0 ) cycle
                    ithr  = omp_get_thread_num() + 1
                    cmats(:,:,ithr) = CMPLX_ZERO
                    tvals(:,:,ithr) = 0.0
                    ! prep image: noise normalization, edge tappering, padding, shift, copy into cmats
                    call build_glob%imgbatch(i)%norm_noise_taper_edge_pad_fft_shift_2mat(build_glob%lmsk, tmp_pad_imgs(ithr),&
                        &-precs(iprec)%shift*croppd_scale, ldim_croppd, ft_map_lims, cmats(:,:,ithr))
                    ! apply CTF to padded image, stores CTF values
                    call precs(iprec)%tfun%eval_and_apply(ctfflag, precs(iprec)%dfx, precs(iprec)%dfy,&
                        &precs(iprec)%angast, fdims_croppd(1:2), tvals(:,:,ithr), cmats(:,:,ithr))
                    ! upsample sigma2 & multipliy CTF.image & CTF2
                    if( params_glob%l_ml_reg )then
                        sigma2(sigma2_kfromto(1):nyq_crop,ithr) = eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):nyq_crop,iptcl)
                        call upsample_sigma2(sigma2_kfromto(1), nyq_crop, sigma2(:,ithr), nyq_croppd, sigma2pd(:,ithr))
                        sqrt_sigma2pd(:,ithr) = sqrt(sigma2pd(:,ithr))
                        do h = ft_map_lims(1,1),ft_map_lims(1,2)
                            do k = ft_map_lims(2,1),ft_map_lims(2,2)
                                sh = nint(hyp(h,k))
                                if( sh > nyq_croppd )cycle
                                physh = ft_map_phys_addrh(h,k)
                                physk = ft_map_phys_addrk(h,k)
                                cmats(physh,physk,ithr) = cmats(physh,physk,ithr) /      sigma2pd(sh,ithr)
                                tvals(physh,physk,ithr) = tvals(physh,physk,ithr) / sqrt_sigma2pd(sh,ithr)
                            enddo
                        enddo
                    end if
                    ! Rotation matrix
                    call rotmat2d(precs(iprec)%e3, mat)
                    ! Kaiser-Bessel Interpolation
                    interp_cmats(:,:,i) = CMPLX_ZERO
                    interp_rhos(:,:,i)  = 0.0
                    ! loop over cropped original image limits
                    do h = flims_crop(1,1), flims_crop(1,2)
                        hp = h * STRIDE_GRID_PAD_FAC        ! padded coordinate
                        do k = flims_crop(2,1), flims_crop(2,2)
                            sh = nint(hyp(real(h),real(k)))
                            if( sh > nyq_crop )cycle
                            kp = k * STRIDE_GRID_PAD_FAC    ! padded coordinate
                            ! rotation on original lattice
                            loc = matmul(real([h,k]),mat)
                            ! interpolation window limits on original lattice
                            win(1,:) = nint(loc)
                            win(2,:) = win(1,:) + iwinsz
                            win(1,:) = win(1,:) - iwinsz
                            ! kernel window
                            call kbwin%apod_mat_2d(loc, iwinsz, wdim, kbw(:,:,ithr))
                            ! Physical address in padded image
                            physh  = ft_map_phys_addrh(hp,kp)
                            physk  = ft_map_phys_addrk(hp,kp)
                            ! padded image component with particle weight & padding correction
                            fcomp  = precs(iprec)%pw * KB2 * cmats(physh, physk, ithr)
                            fcomp  = merge(conjg(fcomp), fcomp, hp<0)
                            ! padded CTF2 with particle weight
                            tvalsq = precs(iprec)%pw * tvals(physh, physk, ithr)**2
                            ! Splat update
                            do l = 1, wdim
                                hh      = win(1,1)+ l-1
                                l_conjg = hh<0
                                hh      = cyci_1d(cyc_lims_cropR(:,1), hh)
                                do m = 1, wdim
                                    kk   = win(1,2) + m-1
                                    kk   = cyci_1d(cyc_lims_cropR(:,2), kk)
                                    phys = cavgs%even%fit%comp_addr_phys(hh,kk)
                                    w    = kbw(l,m,ithr)
                                    interp_cmats(phys(1),phys(2),i) = interp_cmats(phys(1),phys(2),i) + w * merge(conjg(fcomp), fcomp, l_conjg)
                                    interp_rhos( phys(1),phys(2),i) = interp_rhos( phys(1),phys(2),i) + w * tvalsq
                                enddo
                            end do
                        end do
                    end do
                enddo
                !$omp end do
                ! Accumulate sums
                !$omp do schedule(static)
                do icls = 1,ncls
                    do i = 1,nptcls_in_batch
                        iprec = batch_iprecs(i)
                        if( iprec == 0 ) cycle
                        if( precs(iprec)%pind == 0 ) cycle
                        if( precs(iprec)%class == icls )then
                            select case(precs(iprec)%eo)
                            case(0,-1)
                                cavgs%even%cmat(:,:,icls)  = cavgs%even%cmat(:,:,icls)  + interp_cmats(:,:,i)
                                cavgs%even%ctfsq(:,:,icls) = cavgs%even%ctfsq(:,:,icls) + interp_rhos(:,:,i)
                            case(1)
                                cavgs%odd%cmat(:,:,icls)  = cavgs%odd%cmat(:,:,icls)  + interp_cmats(:,:,i)
                                cavgs%odd%ctfsq(:,:,icls) = cavgs%odd%ctfsq(:,:,icls) + interp_rhos(:,:,i)
                            end select
                        endif
                    enddo
                enddo
                !$omp end do
                !$omp end parallel
            enddo ! end read batches loop
        enddo
        ! Cleanup
        call dstkio_r%kill
        call forget_ft_maps
        call dealloc_imgarr(tmp_pad_imgs)
        deallocate(cmats,interp_cmats,tvals,kbw,interp_rhos)
    end subroutine cavger_assemble_sums

    !>  \brief  merges the even/odd pairs and normalises the sums, merge low resolution
    !    frequencies, calculates & writes FRCs and optionally applies regularization
    module subroutine cavger_restore_cavgs( frcs_fname )
        use simple_gridding, only: prep2D_inv_instrfun4mul
        class(string), intent(in) :: frcs_fname
        real, allocatable :: frc(:)
        type(cavgs_set)    :: cavgs_bak
        type(stack)       :: even_tmp, odd_tmp
        type(image)       :: gridcorr_img
        integer           :: eo_pop(2), icls, ithr, find, pop, filtsz_crop
        ! temporary objects for frc calculation & regularization
        filtsz_crop = fdim(ldim_crop(1))-1
        allocate(frc(filtsz_crop),source=0.0)
        call cavgs_bak%new_set(ldim_crop, ncls)
        call even_tmp%new_stack(ldim_crop, nthr_glob, alloc_ctfsq=.false.)
        call odd_tmp%new_stack( ldim_crop, nthr_glob, alloc_ctfsq=.false.)
        if( params_glob%l_ml_reg ) call cavgs_bak%new_set(ldim_crop, ncls)
        call memoize_ft_maps(ldim_crop(1:2), smpd_crop)
        gridcorr_img = prep2D_inv_instrfun4mul(ldim_crop, ldim_croppd, smpd_crop)
        ! Main loop
        !$omp parallel do default(shared) private(icls,ithr,eo_pop,pop,frc)&
        !$omp schedule(static) proc_bind(close)
        do icls = 1,ncls
            eo_pop = eo_class_pop(icls)
            pop    = sum(eo_pop)
            frc    = 0.0
            if(pop == 0)then
                ! empty class
                call cavgs%even%zero_slice(icls, .false.)
                call cavgs%odd%zero_slice(icls, .false.)
                call cavgs%merged%zero_slice(icls, .false.)
            else
                ithr = omp_get_thread_num() + 1
                ! even + odd
                cavgs%merged%slices(icls)%ft = .true.
                cavgs%merged%cmat(:,:,icls)  = cavgs%even%cmat(:,:,icls)  + cavgs%odd%cmat(:,:,icls)
                cavgs%merged%ctfsq(:,:,icls) = cavgs%even%ctfsq(:,:,icls) + cavgs%odd%ctfsq(:,:,icls)
                ! backup current classes
                if( params_glob%l_ml_reg ) call cavgs_bak%copy_fast(cavgs, icls, .true.)
                ! CTF2 density correction
                if( eo_pop(1) > 1 ) call cavgs%even%ctf_dens_correct(icls)
                if( eo_pop(2) > 1 ) call cavgs%odd%ctf_dens_correct(icls)
                if( pop       > 1 ) call cavgs%merged%ctf_dens_correct(icls)
                ! iFT
                call cavgs%even%ifft(icls)
                call cavgs%odd%ifft(icls)
                call cavgs%merged%ifft(icls)
                ! FRC calculation
                even_tmp%rmat(:,:,ithr)  = cavgs%even%rmat(:,:,icls)
                odd_tmp%rmat(:,:,ithr)   = cavgs%odd%rmat(:,:,icls)
                even_tmp%slices(ithr)%ft = .false.
                odd_tmp%slices(ithr)%ft  = .false.
                call even_tmp%softmask(ithr)
                call odd_tmp%softmask(ithr)
                call even_tmp%fft(ithr)
                call odd_tmp%fft(ithr)
                call even_tmp%frc(odd_tmp, ithr, frc)
                ! ML-regularization: add inverse of noise power to ctfsq & normalize again
                if( params_glob%l_ml_reg )then
                    ! add noise power term to denominator
                    call cavgs_bak%even%add_invnoisepower2rho(icls, filtsz_crop, frc)
                    call cavgs_bak%odd%add_invnoisepower2rho(icls, filtsz_crop, frc)
                    if( eo_pop(1) < 3 ) cavgs_bak%even%ctfsq(:,:,icls) = cavgs_bak%even%ctfsq(:,:,icls) + 1.0
                    if( eo_pop(2) < 3 ) cavgs_bak%odd%ctfsq(:,:,icls)  = cavgs_bak%odd%ctfsq(:,:,icls)  + 1.0
                    cavgs_bak%merged%ctfsq(:,:,icls) = cavgs_bak%even%ctfsq(:,:,icls) + cavgs_bak%odd%ctfsq(:,:,icls)
                    ! re-normalize cavg
                    call cavgs_bak%even%ctf_dens_correct(icls)
                    call cavgs_bak%odd%ctf_dens_correct(icls)
                    call cavgs_bak%merged%ctf_dens_correct(icls)
                    ! transfer back cavgs & iFT
                    call cavgs%copy_fast(cavgs_bak, icls, .true.)
                    call cavgs%even%ifft(icls)
                    call cavgs%odd%ifft(icls)
                    call cavgs%merged%ifft(icls)
                endif
                ! average low-resolution info between eo pairs
                find = build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1)
                call cavgs%merged%fft(icls)
                call cavgs%even%fft(icls)
                call cavgs%odd%fft(icls)
                call cavgs%even%insert_lowres_serial(cavgs%merged, icls, find)
                call cavgs%odd%insert_lowres_serial(cavgs%merged, icls, find)
                ! gridding correction
                call cavgs%merged%ifft(icls)
                call cavgs%even%ifft(icls)
                call cavgs%odd%ifft(icls)
                call gridcorr_img%mul_rmat(cavgs%even%rmat(:ldim_crop(1),:ldim_crop(2),icls:icls))
                call gridcorr_img%mul_rmat(cavgs%odd%rmat(:ldim_crop(1),:ldim_crop(2),icls:icls))
                call gridcorr_img%mul_rmat(cavgs%merged%rmat(:ldim_crop(1),:ldim_crop(2),icls:icls))
            endif
            ! store FRC
            call build_glob%clsfrcs%set_frc(icls, frc, 1)
        end do
        !$omp end parallel do
        ! write FRCs
        call build_glob%clsfrcs%write(frcs_fname)
        ! cleanup
        call gridcorr_img%kill
        call forget_ft_maps
        call even_tmp%kill_stack
        call odd_tmp%kill_stack
        call cavgs_bak%kill_set
        deallocate(frc)
    end subroutine cavger_restore_cavgs

    ! I/O

    module subroutine cavger_write_eo( fname_e, fname_o )
        class(string), intent(in) :: fname_e, fname_o
        call cavgs%even%write(fname_e, .false.)
        call cavgs%odd%write(fname_o, .false.)
    end subroutine cavger_write_eo

    module subroutine cavger_write_all( fname, fname_e, fname_o )
        class(string), intent(in) :: fname, fname_e, fname_o
        call cavgs%merged%write(fname, .false.)
        call cavger_write_eo( fname_e, fname_o )
    end subroutine cavger_write_all

    module subroutine cavger_read_all()
        if( .not. file_exists(params_glob%refs) ) THROW_HARD('references (REFS) does not exist in cwd')
        call read_cavgs(params_glob%refs, 'merged')
        if( file_exists(params_glob%refs_even) )then
            call read_cavgs(params_glob%refs_even, 'even')
        else
            call read_cavgs(params_glob%refs, 'even')
        endif
        if( file_exists(params_glob%refs_odd) )then
            call read_cavgs(params_glob%refs_odd, 'odd')
        else
            call read_cavgs(params_glob%refs, 'odd')
        endif
    end subroutine cavger_read_all

    !>  \brief  submodule utility for reading class averages (image type)
    subroutine read_cavgs( fname, which )
        use simple_imgarr_utils, only: read_stk_into_imgarr
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        class(image), pointer :: cavgs(:)
        integer               :: ldim_read(3), icls
        ! to use a local subroutine instead of pointers
        ! read
        select case(trim(which))
            case('even')
                cavgs_even = read_stk_into_imgarr(fname)
                cavgs => cavgs_even
            case('odd')
                cavgs_odd = read_stk_into_imgarr(fname)
                cavgs => cavgs_odd
            case('merged')
                cavgs_merged = read_stk_into_imgarr(fname)
                cavgs => cavgs_merged
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
        ldim_read = cavgs(1)%get_ldim()
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
    end subroutine read_cavgs

    !>  \brief  writes partial class averages to disk (distributed execution)
    module subroutine cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        type(string)   :: cae, cao, cte, cto
        cae   = 'cavgs_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        cao   = 'cavgs_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        cte   = 'ctfsqsums_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        cto   = 'ctfsqsums_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext%to_char()
        select case(trim(which))
            case('read')
                call cavgs%even%read_cmat(cae)
                call cavgs%odd%read_cmat(cao)
                call cavgs%even%read_ctfsq(cte)
                call cavgs%odd%read_ctfsq(cto)
            case('write')
                call cavgs%even%write(cae,.true.)
                call cavgs%odd%write(cao,.true.)
                call cavgs%even%write_ctfsq(cte)
                call cavgs%odd%write_ctfsq(cto)
            case DEFAULT
                THROW_HARD('unknown which flag; only read & write supported; cavger_readwrite_partial_sums')
        end select
        call cae%kill
        call cao%kill
        call cte%kill
        call cto%kill
    end subroutine cavger_readwrite_partial_sums

    module subroutine apply_weights2cavgs( w )
        real, intent(in) :: w
        !$omp parallel workshare proc_bind(close)
        cavgs%even%cmat  = w * cavgs%even%cmat
        cavgs%even%ctfsq = w * cavgs%even%ctfsq
        cavgs%odd%cmat   = w * cavgs%odd%cmat
        cavgs%odd%ctfsq  = w * cavgs%odd%ctfsq
        !$omp end parallel workshare
    end subroutine apply_weights2cavgs

    !>  \brief  generates the cavgs parts after distributed execution
    module subroutine cavger_assemble_sums_from_parts
        integer(timer_int_kind) ::  t_init,  t_io,  t_sum, t_merge_eos_and_norm,  t_tot
        real(timer_int_kind)    :: rt_init, rt_io, rt_sum, rt_merge_eos_and_norm, rt_tot
        type(string) :: cae, cao, cte, cto, benchfname
        type(stack)  :: cavgs4reade, cavgs4reado
        integer      :: ipart, fnr
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call cavgs%zero_set(.true.)
        call cavgs4reade%new_stack(ldim_crop(1:2), ncls)
        call cavgs4reado%new_stack(ldim_crop(1:2), ncls)
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            rt_io   = 0.
            rt_sum  = 0.
        endif
        do ipart=1,params_glob%nparts
            if( L_BENCH_GLOB ) t_io = tic()
            ! filenames
            cae = 'cavgs_even_part'    //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            cao = 'cavgs_odd_part'     //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            cte = 'ctfsqsums_even_part'//int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            cto = 'ctfsqsums_odd_part' //int2str_pad(ipart,params_glob%numlen)//params_glob%ext%to_char()
            ! read arrays
            call cavgs4reade%read_cmat(cae)
            call cavgs4reade%read_ctfsq(cte)
            call cavgs4reado%read_cmat(cao)
            call cavgs4reado%read_ctfsq(cto)
            if( L_BENCH_GLOB )then
                rt_io = rt_io + toc(t_io)
                t_sum = tic()
            endif
            ! sum arrays
            !$omp parallel workshare proc_bind(close)
            cavgs%even%cmat  = cavgs%even%cmat  + cavgs4reade%cmat
            cavgs%even%ctfsq = cavgs%even%ctfsq + cavgs4reade%ctfsq
            cavgs%odd%cmat   = cavgs%odd%cmat   + cavgs4reado%cmat
            cavgs%odd%ctfsq  = cavgs%odd%ctfsq  + cavgs4reado%ctfsq
            !$omp end parallel workshare
            if( L_BENCH_GLOB ) rt_sum = rt_sum + toc(t_sum)
        enddo
        call cae%kill
        call cao%kill
        call cte%kill
        call cto%kill
        call cavgs4reade%kill_stack
        call cavgs4reado%kill_stack
        ! Restoration of e/o/merged classes
        if( L_BENCH_GLOB ) t_merge_eos_and_norm = tic()
        call cavger_restore_cavgs(params_glob%frcs)
        ! Benchmarck
        if( L_BENCH_GLOB )then
            rt_merge_eos_and_norm = toc(t_merge_eos_and_norm)
            rt_tot                = toc(t_tot)
            benchfname = 'CAVGASSEMBLE_BENCH.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation       : ', rt_init
            write(fnr,'(a,1x,f9.2)') 'I/O                  : ', rt_io
            write(fnr,'(a,1x,f9.2)') 'workshare sum        : ', rt_sum
            write(fnr,'(a,1x,f9.2)') 'merge eo-pairs & norm: ', rt_merge_eos_and_norm
            write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation        : ', (rt_init/rt_tot)               * 100.
            write(fnr,'(a,1x,f9.2)') 'I/O                   : ', (rt_io/rt_tot)                 * 100.
            write(fnr,'(a,1x,f9.2)') 'workshare sum         : ', (rt_sum/rt_tot)                * 100.
            write(fnr,'(a,1x,f9.2)') 'merge eo-pairs & norm : ', (rt_merge_eos_and_norm/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for       : ',&
            &((rt_init+rt_io+rt_sum+rt_merge_eos_and_norm)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine cavger_assemble_sums_from_parts

    ! DESTRUCTOR

    !>  \brief  is a destructor
    module subroutine cavger_kill( dealloccavgs )
        logical, optional, intent(in) :: dealloccavgs
        if( present(dealloccavgs) )then
            if( dealloccavgs ) call dealloc_cavgs
        else
            call dealloc_cavgs
        endif
        if( allocated(pptcl_mask) ) deallocate(pptcl_mask)
        if( allocated(precs)      ) deallocate(precs)
    end subroutine cavger_kill

    !>  \brief submodule private destructor utility
    subroutine dealloc_cavgs
        call dealloc_imgarr(cavgs_even)
        call dealloc_imgarr(cavgs_odd)
        call dealloc_imgarr(cavgs_merged)
        call cavgs%kill_set
        istart  = 0; iend = 0
        partsz  = 0
        ncls    = 0
        ldim    = 0; ldim_crop   = 0
        ldim_pd = 0; ldim_croppd = 0
        smpd    = 0.;smpd_crop   = 0.
        l_alloc_read_cavgs = .true.
    end subroutine dealloc_cavgs

    ! PUBLIC UTILITIES

    module subroutine transform_ptcls( spproj, oritype, icls, timgs, pinds, phflip, cavg, imgs_ori)
        use simple_sp_project,          only: sp_project
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimgbatch
        use simple_memoize_ft_maps
        class(sp_project),                  intent(inout) :: spproj
        character(len=*),                   intent(in)    :: oritype
        integer,                            intent(in)    :: icls
        type(image),           allocatable, intent(inout) :: timgs(:)
        integer,               allocatable, intent(inout) :: pinds(:)
        logical,     optional,              intent(in)    :: phflip
        type(image), optional,              intent(inout) :: cavg
        type(image), optional, allocatable, intent(inout) :: imgs_ori(:)
        class(oris),  pointer :: pos
        type(kbinterpol)      :: kbwin
        type(image)           :: img(nthr_glob), timg(nthr_glob)
        type(ctfparams)       :: ctfparms
        type(ctf)             :: tfun
        type(string)          :: str
        real,     allocatable :: kbw(:,:)
        complex :: fcomp, fcompl
        real    :: mat(2,2), shift(2), loc(2), e3
        integer :: logi_lims(3,2),cyc_lims(3,2),cyc_limsR(2,2),win(2,2)
        integer :: i,iptcl, l,m, pop, h,k, hh,kk, ithr, iwinsz, wdim, physh,physk
        logical :: l_phflip, l_imgs, l_conjg
        l_imgs = .false.
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
        ! interpolation variables
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        allocate(kbw(wdim,wdim),source=0.)
        ! temporary objects
        call prepimgbatch(pop)
        do ithr = 1, nthr_glob
            call  img(ithr)%new([params_glob%boxpd,params_glob%boxpd,1],params_glob%smpd, wthreads=.false.)
            call timg(ithr)%new([params_glob%boxpd,params_glob%boxpd,1],params_glob%smpd, wthreads=.false.)
        end do
        call memoize_ft_maps(img(1)%get_ldim(), img(1)%get_smpd())
        logi_lims      = img(1)%loop_lims(2)
        cyc_lims       = img(1)%loop_lims(3)
        cyc_limsR(:,1) = cyc_lims(1,:)
        cyc_limsR(:,2) = cyc_lims(2,:)
        call discrete_read_imgbatch(pop, pinds(:), [1,pop])
        !$omp parallel do private(i,ithr,iptcl,shift,e3,ctfparms,tfun,mat,h,k,hh,kk,loc,win,l,m,physh,physk,kbw,fcomp,fcompl,l_conjg) &
        !$omp default(shared) schedule(static) proc_bind(close)
        do i = 1,pop
            ithr  = omp_get_thread_num() + 1
            iptcl = pinds(i)
            shift = pos%get_2Dshift(iptcl)
            e3    = pos%e3get(iptcl)
            call img(ithr)%zero_and_flag_ft
            call timg(ithr)%zero_and_flag_ft
            ! normalisation
            call build_glob%imgbatch(i)%norm_noise_taper_edge_pad_fft(build_glob%lmsk,img(ithr))
            if( l_imgs )then
                call img(ithr)%ifft
                call img(ithr)%clip(imgs_ori(i))
                call img(ithr)%fft
            endif
            ! optional phase-flipping
            if( l_phflip )then
                ctfparms = spproj%get_ctfparams(oritype, iptcl)
                tfun     = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                call img(ithr)%apply_ctf(tfun, 'flip', ctfparms)
            endif
            ! shift
            call img(ithr)%shift2Dserial(-shift)
            ! rotation matrix
            call rotmat2d(-e3, mat)
            do h = logi_lims(1,1),logi_lims(1,2)
                do k = logi_lims(2,1),logi_lims(2,2)
                    ! rotation
                    loc = matmul(real([h,k]),mat)
                    ! interpolation window
                    win(1,:) = nint(loc)
                    win(2,:) = win(1,:) + iwinsz
                    win(1,:) = win(1,:) - iwinsz
                    ! interpolation kernel
                    call kbwin%apod_mat_2d(loc, iwinsz, wdim, kbw)
                    ! intepolation
                    fcomp = CMPLX_ZERO
                    do l = 1,wdim
                        hh      = win(1,1)+l-1
                        l_conjg = hh < 0
                        hh      = cyci_1d(cyc_limsR(:,1), hh)
                        fcompl  = CMPLX_ZERO
                        do m = 1,wdim
                            kk     = win(1,2)+m-1
                            kk     = cyci_1d(cyc_limsR(:,2), kk)
                            physh  = ft_map_phys_addrh(hh,kk)
                            physk  = ft_map_phys_addrk(hh,kk)
                            fcompl = fcompl + kbw(l,m) * img(ithr)%get_cmat_at(physh,physk,1)
                        enddo
                        fcomp = fcomp + merge(conjg(fcompl), fcompl, l_conjg)
                    end do
                    physh = ft_map_phys_addrh(h,k)
                    physk = ft_map_phys_addrk(h,k)
                    call timg(ithr)%set_cmat_at(physh, physk, 1, fcomp)
                enddo
            enddo
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
        call forget_ft_maps
        do ithr = 1, nthr_glob
            call img(ithr)%kill
            call timg(ithr)%kill
        end do
        nullify(pos)
    end subroutine transform_ptcls

end submodule simple_new_classaverager_restore
