!@descr: Routines to perform the classes restoration and processing
submodule (simple_classaverager) simple_classaverager_restore
use simple_imgarr_utils, only: alloc_imgarr, dealloc_imgarr
implicit none
#include "simple_local_flags.inc"

contains

    !>  \brief  Constructor
    module subroutine cavger_new( params, build, alloccavgs )
        class(parameters), target, intent(inout) :: params
        class(builder),    target, intent(inout) :: build
        logical, optional,         intent(in)    :: alloccavgs
        p_ptr => params
        b_ptr => build
        l_alloc_read_cavgs = .true.
        if( present(alloccavgs) ) l_alloc_read_cavgs = alloccavgs
        call cavger_kill(dealloccavgs=l_alloc_read_cavgs)
        ncls          = p_ptr%ncls
        ! CTF logics
        ctfflag       = b_ptr%spproj%get_ctfflag_type('ptcl2D',iptcl=p_ptr%fromp)
        ! smpd
        smpd          = p_ptr%smpd
        smpd_crop     = p_ptr%smpd_crop
        ! set ldims
        ldim          = [p_ptr%box,       p_ptr%box,       1]
        ldim_crop     = [p_ptr%box_crop,  p_ptr%box_crop,  1]
        ldim_croppd   = [p_ptr%box_croppd,p_ptr%box_croppd,1]
        ldim_pd       = [p_ptr%boxpd,     p_ptr%boxpd,     1]
        ! instantiate class averages
        if( l_alloc_read_cavgs ) call cavgs%new_set(ldim_crop(1:2), ncls)
        ! populations
        allocate(eo_pops(2,ncls),source=0)
    end subroutine cavger_new

    ! setters/getters

    !>  \brief  transfers metadata to the instance for a subset of particles
    module subroutine cavger_transf_oridat( nptcls, pinds )
        integer,  intent(in) :: nptcls
        integer,  intent(in) :: pinds(nptcls)
        class(oris), pointer :: spproj_field
        type(ctfparams)      :: ctfparms(nthr_glob)
        integer              :: i, iptcl, ithr, stkind
        ! indices
        precs(1:nptcls)%pind = pinds(:)
        if( nptcls < size(precs) ) precs(nptcls+1:)%pind = 0
        ! fetch data from project
        call b_ptr%spproj%ptr2oritype(p_ptr%oritype, spproj_field)
        !$omp parallel do default(shared) private(i,iptcl,ithr,stkind) schedule(static) proc_bind(close)
        do i = 1,nptcls
            iptcl = pinds(i)
            if( iptcl == 0 ) cycle
            precs(i)%pw = spproj_field%get(iptcl,'w')
            if( (spproj_field%get_state(iptcl)==0).or.(precs(i)%pw<SMALL) )then
                precs(i)%pind = 0
            else
                ithr            = omp_get_thread_num() + 1
                precs(i)%pind   = iptcl
                precs(i)%eo     = spproj_field%get_eo(iptcl)
                ctfparms(ithr)  = b_ptr%spproj%get_ctfparams(p_ptr%oritype,iptcl)
                precs(i)%tfun   = ctf(p_ptr%smpd_crop, ctfparms(ithr)%kv, ctfparms(ithr)%cs, ctfparms(ithr)%fraca)
                precs(i)%dfx    = ctfparms(ithr)%dfx
                precs(i)%dfy    = ctfparms(ithr)%dfy
                precs(i)%angast = ctfparms(ithr)%angast
                precs(i)%class  = spproj_field%get_class(iptcl)
                precs(i)%e3     = spproj_field%e3get(iptcl)
                precs(i)%shift  = spproj_field%get_2Dshift(iptcl)
                call b_ptr%spproj%map_ptcl_ind2stk_ind(p_ptr%oritype, iptcl, stkind, precs(i)%ind_in_stk)
            endif
        end do
        !$omp end parallel do
        nullify(spproj_field)
    end subroutine cavger_transf_oridat

    !>  \brief  for loading sigma2
    module subroutine cavger_read_euclid_sigma2
        type(string) :: fname
        if( p_ptr%l_ml_reg )then
            fname = SIGMA2_FBODY//int2str_pad(p_ptr%part,p_ptr%numlen)//'.dat'
            call b_ptr%esig%new(p_ptr, b_ptr%pftc, fname, p_ptr%box)
            call b_ptr%esig%read_part(  b_ptr%spproj_field)
            call b_ptr%esig%read_groups(b_ptr%spproj_field)
        end if
    end subroutine cavger_read_euclid_sigma2

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population, average correlation and weight
    module subroutine cavger_gen2Dclassdoc()
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(p_ptr%ncls)
        real(dp) :: corrs(p_ptr%ncls), ws(p_ptr%ncls)
        real     :: frc05, frc0143, rstate, w
        integer  :: iptcl, icls, pop, nptcls
        select case(trim(p_ptr%oritype))
            case('ptcl2D')
                ptcl_field => b_ptr%spproj%os_ptcl2D
                cls_field  => b_ptr%spproj%os_cls2D
            case('ptcl3D')
                ptcl_field => b_ptr%spproj%os_ptcl3D
                cls_field  => b_ptr%spproj%os_cls3D
            case DEFAULT
                THROW_HARD('Unsupported ORITYPE: '//trim(p_ptr%oritype))
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
            if( icls<1 .or. icls>p_ptr%ncls )cycle
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
        call cls_field%new(p_ptr%ncls, is_ptcl=.false.)
        do icls=1,p_ptr%ncls
            pop = pops(icls)
            call b_ptr%clsfrcs%estimate_res(icls, frc05, frc0143)
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

    ! Calculators

    ! Initialize objects for on-the-fly classes update
    module subroutine cavger_init_online( maxbatchsz, do_frac_update )
        integer, intent(in) :: maxbatchsz
        logical, intent(in) :: do_frac_update
        integer :: fdims_croppd(3), cshape_crop(2), wdim
        ! Zero sums or set to previous with weight
        if( l_alloc_read_cavgs )then
            call cavgs%zero_set(.true.)
            if( do_frac_update )then
                call cavger_readwrite_partial_sums('read')
                call apply_weights2cavgs(1.0-p_ptr%update_frac)
            endif
        else
            if( do_frac_update ) call apply_weights2cavgs(1.0-p_ptr%update_frac)
            call cavgs%zero_set(.true.)
        endif
        ! Interpolation variables
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        wdim   = kbwin%get_wdim()
        allocate(kbw(wdim,wdim,nthr_glob),source=0.)
        ! Work images
        call alloc_imgarr(nthr_glob, ldim_pd, smpd, tmp_pad_imgs)
        ! particle records
        allocate(precs(maxbatchsz))
        precs(:)%pind = 0
        ! populations
        eo_pops(:,:) = 0
        ! Memoization for cropped image, is overwritten just below
        call memoize_ft_maps(ldim_crop(1:2), p_ptr%smpd_crop)
        phys_addrh_crop = ft_map_phys_addrh
        phys_addrk_crop = ft_map_phys_addrk
        ! Memoization for cropped padded image, will be overwritten during search
        call memoize_ft_maps(ldim_croppd(1:2), p_ptr%smpd_crop)
        ! Work arrays
        cshape_crop  = cavgs%even%cshape
        fdims_croppd = ft_map_get_farray_shape()
        allocate(tvals(fdims_croppd(1),fdims_croppd(2),nthr_glob),&
                &cmats(fdims_croppd(1),fdims_croppd(2),nthr_glob),&
                &interp_cmats(cshape_crop(1),cshape_crop(2),maxbatchsz),&
                &interp_rhos(cshape_crop(1),cshape_crop(2),maxbatchsz))
    end subroutine cavger_init_online

    ! Deallocate objects  on-the-fly classes update
    module subroutine cavger_dealloc_online()
        if( allocated(tmp_pad_imgs))then
            call forget_ft_maps
            call dealloc_imgarr(tmp_pad_imgs)
            deallocate(cmats,interp_cmats,tvals,kbw,interp_rhos,precs,&
                &phys_addrh_crop,phys_addrk_crop)
        endif
    end subroutine cavger_dealloc_online

    !>  \brief  is for updating classes sums in distributed/non-distributed mode
    !           with provided images using interpolation in Fourier space
    module subroutine cavger_update_sums( nptcls, ptcl_imgs )
        integer,      intent(in)    :: nptcls
        class(image), intent(inout) :: ptcl_imgs(nptcls)
        real, parameter :: KB2 = KBALPHA**2
        type(fplane_type) :: fplanes(nthr_glob)
        complex :: fcomp
        real    :: loc(2), mat(2,2), tvalsq, croppd_scale, w
        integer :: win(2,2), flims_crop(3,2), phys(2)
        integer :: cyc_lims_cropR(2,2), cyc_lims_crop(3,2), sigma2_kfromto(2)
        integer :: h,k,hh,kk,hp,kp,l,m, icls, nyq_crop
        integer :: iptcl, i, sh, iwinsz, wdim, ithr
        logical :: l_conjg
        ! Interpolation parameters
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        wdim   = kbwin%get_wdim()
        ! Memoization for full padded image (tmp_pad_imgs is allocated with ldim_pd)
        call memoize_ft_maps(ldim_pd(1:2), p_ptr%smpd)
        ! Dimensions & limits
        croppd_scale  = real(ldim_croppd(1)) / real(ldim_pd(1))
        flims_crop    = cavgs%even%flims
        nyq_crop      = cavgs%even%fit%get_lfny(1)
        cyc_lims_crop = cavgs%even%fit%loop_lims(3)
        cyc_lims_cropR = transpose(cyc_lims_crop(1:2,:))
        sigma2_kfromto = [1, nyq_crop]
        if( p_ptr%l_ml_reg ) then
            sigma2_kfromto(1) = lbound(b_ptr%esig%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(b_ptr%esig%sigma2_noise,1)
        end if
        !$omp parallel default(shared) proc_bind(close) &
        !$omp private(i,ithr,iptcl,win,mat,h,k,hh,kk,l,m,loc,sh,hp,kp,phys,w,tvalsq,l_conjg,fcomp)
        !$omp do schedule(static)
        do i = 1, nptcls
            iptcl = precs(i)%pind
            if( iptcl == 0 ) cycle
            ithr = omp_get_thread_num() + 1
            ! prep image: noise normalization, edge tapering, padding, fftshift, FFT
            call ptcl_imgs(i)%norm_noise_taper_edge_pad_fft( b_ptr%lmsk, tmp_pad_imgs(ithr) )
            ! Get and cache CTF parameters
            precs(i)%ctfparams         = b_ptr%spproj%get_ctfparams(p_ptr%oritype, iptcl)
            precs(i)%ctfparams%ctfflag = ctfflag
            precs(i)%ctfparams%dfx     = precs(i)%dfx
            precs(i)%ctfparams%dfy     = precs(i)%dfy
            precs(i)%ctfparams%angast  = precs(i)%angast
            ! Build Fourier plane with shift, CTF and optional ML regularization
            ! gen_fplane4rec applies the minus sign internally.
            if( p_ptr%l_ml_reg ) then
                call tmp_pad_imgs(ithr)%gen_fplane4rec( sigma2_kfromto, p_ptr%smpd_crop, &
                    precs(i)%ctfparams, precs(i)%shift, fplanes(ithr), &
                    b_ptr%esig%sigma2_noise(sigma2_kfromto(1):sigma2_kfromto(2), iptcl) )
            else
                call tmp_pad_imgs(ithr)%gen_fplane4rec( sigma2_kfromto, p_ptr%smpd_crop, &
                    precs(i)%ctfparams, precs(i)%shift, fplanes(ithr) )
            endif
            ! Rotation matrix
            call rotmat2d( precs(i)%e3, mat )
            ! Kaiser-Bessel interpolation
            interp_cmats(:,:,i) = CMPLX_ZERO
            interp_rhos(:,:,i)  = 0.0
            ! loop over cropped original image limits
            do h = flims_crop(1,1), flims_crop(1,2)
                hp = h * OSMPL_PAD_FAC
                do k = flims_crop(2,1), flims_crop(2,2)
                    sh = nint(hyp(real(h),real(k)))
                    if( sh > nyq_crop ) cycle
                    kp = k * OSMPL_PAD_FAC
                    ! rotation on original lattice
                    loc = matmul(real([h,k]), mat)
                    ! interpolation window limits on original lattice
                    win(1,:) = nint(loc)
                    win(2,:) = win(1,:) + iwinsz
                    win(1,:) = win(1,:) - iwinsz
                    ! kernel window
                    call kbwin%apod_mat_2d(loc, iwinsz, wdim, kbw(:,:,ithr))
                    ! Read from the generated Fourier plane.
                    ! gen_fplane4rec stores only k<=0, so use Friedel symmetry for kp>0.
                    if( kp <= 0 ) then
                        fcomp  = precs(i)%pw * KB2 * fplanes(ithr)%cmplx_plane(hp, kp)
                        tvalsq = precs(i)%pw *       fplanes(ithr)%ctfsq_plane(hp, kp)
                    else
                        fcomp  = precs(i)%pw * KB2 * conjg(fplanes(ithr)%cmplx_plane(-hp, -kp))
                        tvalsq = precs(i)%pw *       fplanes(ithr)%ctfsq_plane(-hp, -kp)
                    endif
                    ! Splat update
                    do l = 1, wdim
                        hh      = win(1,1) + l - 1
                        l_conjg = hh < 0
                        hh      = cyci_1d(cyc_lims_cropR(:,1), hh)
                        do m = 1, wdim
                            kk      = win(1,2) + m - 1
                            kk      = cyci_1d(cyc_lims_cropR(:,2), kk)
                            phys(1) = phys_addrh_crop(hh,kk)
                            phys(2) = phys_addrk_crop(hh,kk)
                            w       = kbw(l,m,ithr)
                            interp_cmats(phys(1),phys(2),i) = interp_cmats(phys(1),phys(2),i) + &
                                w * merge(conjg(fcomp), fcomp, l_conjg)
                            interp_rhos(phys(1),phys(2),i)  = interp_rhos(phys(1),phys(2),i)  + &
                                w * tvalsq
                        enddo
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        ! Accumulate sums
        !$omp do schedule(static)
        do icls = 1, ncls
            do i = 1, nptcls
                if( precs(i)%pind == 0 ) cycle
                if( precs(i)%class == icls ) then
                    select case(precs(i)%eo)
                    case(0,-1)
                        cavgs%even%cmat(:,:,icls)  = cavgs%even%cmat(:,:,icls)  + interp_cmats(:,:,i)
                        cavgs%even%ctfsq(:,:,icls) = cavgs%even%ctfsq(:,:,icls) + interp_rhos(:,:,i)
                        eo_pops(1,icls) = eo_pops(1,icls) + 1
                    case(1)
                        cavgs%odd%cmat(:,:,icls)  = cavgs%odd%cmat(:,:,icls)  + interp_cmats(:,:,i)
                        cavgs%odd%ctfsq(:,:,icls) = cavgs%odd%ctfsq(:,:,icls) + interp_rhos(:,:,i)
                        eo_pops(2,icls) = eo_pops(2,icls) + 1
                    end select
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine cavger_update_sums

    !>  \brief  is for generating class averages offline
    module subroutine cavger_assemble_sums( do_frac_update )
        use simple_strategy2D3D_common, only: prepimgbatch, discrete_read_imgbatch, killimgbatch
        logical,  intent(in) :: do_frac_update
        class(oris), pointer :: spproj_field
        type(string)         :: stk_fname
        integer :: pinds(READBUFFSZ)
        integer :: iptcl, i, ind_in_stk, nptcls_eff, batchind, ibatch_start, ibatch_end
        integer :: first_stkind, fromp, istk, nptcls_in_stk, last_stkind, ibatch, nbatches
        ! fetch data from project
        call b_ptr%spproj%ptr2oritype(p_ptr%oritype, spproj_field)
        ! Initialize temporary arrays
        call cavger_init_online(READBUFFSZ, do_frac_update)
        ! Prep for image reading
        call prepimgbatch(p_ptr, b_ptr, READBUFFSZ)
        ! Stack & batch loops
        call b_ptr%spproj%map_ptcl_ind2stk_ind(p_ptr%oritype, p_ptr%fromp, first_stkind, ind_in_stk)
        call b_ptr%spproj%map_ptcl_ind2stk_ind(p_ptr%oritype, p_ptr%top,   last_stkind,  ind_in_stk)
        do istk = first_stkind,last_stkind
            ! Particles range in stack
            fromp         = b_ptr%spproj%os_stk%get_fromp(istk)
            nptcls_in_stk = b_ptr%spproj%os_stk%get_top(istk) - fromp + 1   ! # of particles in stack
            call b_ptr%spproj%get_stkname_and_ind(p_ptr%oritype, max(p_ptr%fromp,fromp), stk_fname, ind_in_stk)
            ! batch loop
            nbatches = ceiling(real(nptcls_in_stk)/real(READBUFFSZ))
            do ibatch = 1,nbatches
                ibatch_start = (ibatch - 1) * READBUFFSZ + 1                    ! first index in current batch
                ibatch_end   = min(nptcls_in_stk, ibatch_start + READBUFFSZ - 1)! last  index in current batch
                ! identify valid particles
                pinds(:)   = 0                                              ! Global valid particle indices in this batch
                batchind   = 0                                              ! ptcl index in batch
                do i = ibatch_start,ibatch_end
                    iptcl = fromp + i - 1                                   ! Global particle index
                    if( (iptcl < p_ptr%fromp).or.(iptcl > p_ptr%top) ) cycle
                    if( (spproj_field%get_state(iptcl)==0)  )          cycle
                    batchind        = batchind + 1                          ! index in batch
                    pinds(batchind) = iptcl
                enddo
                nptcls_eff = batchind                                       ! # valid particles in batch
                if( nptcls_eff == 0 ) cycle
                ! Transfer orientation parameters
                call cavger_transf_oridat( nptcls_eff, pinds(1:nptcls_eff) )
                ! Read images
                call discrete_read_imgbatch(p_ptr, b_ptr, nptcls_eff, pinds(1:nptcls_eff), [1,nptcls_eff])
                ! Interpolate images and update class sums
                call cavger_update_sums(nptcls_eff, b_ptr%imgbatch(1:nptcls_eff))
            enddo   ! batch loop
        enddo       ! stack loop
        ! cleanup
        call stk_fname%kill
        nullify(spproj_field)
        call cavger_dealloc_online
        call killimgbatch(b_ptr)
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
        if( p_ptr%l_ml_reg ) call cavgs_bak%new_set(ldim_crop, ncls)
        call memoize_ft_maps(ldim_crop(1:2), smpd_crop)
        gridcorr_img = prep2D_inv_instrfun4mul(ldim_crop, ldim_croppd, smpd_crop)
        ! Making sure that the public images are allocated with make_cavgs & shared memory
        if( (.not.l_distr_worker_glob).and.(.not.allocated(cavgs_merged)) )then
            call alloc_imgarr(ncls, ldim_crop, smpd_crop, cavgs_even)
            call alloc_imgarr(ncls, ldim_crop, smpd_crop, cavgs_odd)
            call alloc_imgarr(ncls, ldim_crop, smpd_crop, cavgs_merged)
        endif
        ! Main loop
        !$omp parallel do default(shared) private(icls,ithr,eo_pop,pop,frc)&
        !$omp schedule(static) proc_bind(close)
        do icls = 1,ncls
            eo_pop = eo_pops(:,icls)
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
                if( p_ptr%l_ml_reg ) call cavgs_bak%copy_fast(cavgs, icls, .true.)
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
                if( p_ptr%l_ml_reg )then
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
                find = b_ptr%clsfrcs%estimate_find_for_eoavg(icls, 1)
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
            call b_ptr%clsfrcs%set_frc(icls, frc, 1)
            ! Transfer cavg from stack object to image object used in alignment
            ! only in shared memory execution
            if( .not.l_distr_worker_glob )then
                call cavgs_even(icls)%set_rmat(  cavgs%even%rmat(:,:,icls:icls), .false.)
                call cavgs_odd(icls)%set_rmat(   cavgs%odd%rmat(:,:,icls:icls), .false.)
                call cavgs_merged(icls)%set_rmat(cavgs%merged%rmat(:,:,icls:icls), .false.)
            endif
        end do
        !$omp end parallel do
        ! write FRCs
        call b_ptr%clsfrcs%write(frcs_fname)
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
        call cavger_write_merged( fname)
        call cavger_write_eo( fname_e, fname_o )
    end subroutine cavger_write_all

    module subroutine cavger_write_merged( fname)
        class(string), intent(in) :: fname
        call cavgs%merged%write(fname, .false.)
    end subroutine cavger_write_merged

    module subroutine cavger_read_all()
        if( .not. file_exists(p_ptr%refs) ) THROW_HARD('references (REFS) does not exist in cwd')
        call read_cavgs(p_ptr%refs, 'merged')
        if( file_exists(p_ptr%refs_even) )then
            call read_cavgs(p_ptr%refs_even, 'even')
        else
            call read_cavgs(p_ptr%refs, 'even')
        endif
        if( file_exists(p_ptr%refs_odd) )then
            call read_cavgs(p_ptr%refs_odd, 'odd')
        else
            call read_cavgs(p_ptr%refs, 'odd')
        endif
    end subroutine cavger_read_all

    !>  \brief  submodule utility for reading class averages (image type)
    subroutine read_cavgs( fname, which )
        use simple_imgarr_utils, only: read_stk_into_imgarr
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
        class(image), pointer :: pcavgs(:)
        integer               :: ldim_read(3), icls
        ! read
        select case(trim(which))
            case('even')
                cavgs_even = read_stk_into_imgarr(fname)
                pcavgs => cavgs_even
            case('odd')
                cavgs_odd = read_stk_into_imgarr(fname)
                pcavgs => cavgs_odd
            case('merged')
                cavgs_merged = read_stk_into_imgarr(fname)
                pcavgs => cavgs_merged
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
        ldim_read = pcavgs(1)%get_ldim()
        ! scale
        if( any(ldim_read /= ldim_crop) )then
            if( ldim_read(1) > ldim_crop(1) )then
                ! Cropping is not covered
                THROW_HARD('Incompatible cavgs dimensions! ; cavger_read')
            else if( ldim_read(1) < ldim_crop(1) )then
                ! Fourier padding
                !$omp parallel do proc_bind(close) schedule(static) default(shared) private(icls)
                do icls = 1,ncls
                    call pcavgs(icls)%fft
                    call pcavgs(icls)%pad_inplace(ldim_crop)
                    call pcavgs(icls)%ifft
                end do
                !$omp end parallel do
            endif
        endif
        nullify(pcavgs)
    end subroutine read_cavgs

    !>  \brief  writes partial class averages to disk (distributed execution)
    module subroutine cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        type(string)   :: cae, cao, cte, cto
        cae   = 'cavgs_even_part'//int2str_pad(p_ptr%part,p_ptr%numlen)//MRC_EXT
        cao   = 'cavgs_odd_part'//int2str_pad(p_ptr%part,p_ptr%numlen)//MRC_EXT
        cte   = 'ctfsqsums_even_part'//int2str_pad(p_ptr%part,p_ptr%numlen)//MRC_EXT
        cto   = 'ctfsqsums_odd_part'//int2str_pad(p_ptr%part,p_ptr%numlen)//MRC_EXT
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

    !>  \brief  pad partial & ctf squared arrays
    module subroutine cavger_pad_partial_sums( old_box, new_box, n, nparts, numlen )
        integer, intent(in) :: old_box, new_box, n, nparts, numlen
        type(string)        :: ca, ct, str
        type(stack)         :: old, new
        integer :: ipart
        call old%new_stack([old_box, old_box], n, .true.)
        call new%new_stack([new_box, new_box], n, .true.)
        do ipart = 1,nparts
            str = int2str_pad(ipart, numlen)
            ca = string('cavgs_even_part')//str//MRC_EXT
            ct = string('ctfsqsums_even_part')//str//MRC_EXT
            call old%read_cmat(ca)
            call old%read_ctfsq(ct)
            call old%pad(new)
            call new%write(ca,.true.)
            call new%write_ctfsq(ct)
            ca = string('cavgs_odd_part')//str//MRC_EXT
            ct = string('ctfsqsums_odd_part')//str//MRC_EXT
            call old%read_cmat(ca)
            call old%read_ctfsq(ct)
            call old%pad(new)
            call new%write(ca,.true.)
            call new%write_ctfsq(ct)
        enddo
        call old%kill_stack
        call new%kill_stack
        call ca%kill; call ct%kill
    end subroutine cavger_pad_partial_sums

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
        integer      :: iptcl, eo, icls, ipart, fnr
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        ! Calculate populations for restore_cavgs
        !$omp parallel do default(shared) private(iptcl,eo,icls)&
        !$omp proc_bind(close) reduction(+:eo_pops)
        do iptcl = p_ptr%fromp, p_ptr%top
            if( b_ptr%spproj_field%get_state(iptcl) == 0 ) cycle
            if( b_ptr%spproj_field%get(iptcl,'w') < SMALL )cycle
            eo   = b_ptr%spproj_field%get_eo(iptcl) + 1
            icls = b_ptr%spproj_field%get_class(iptcl)
            eo_pops(eo, icls) = eo_pops(eo, icls) + 1
        enddo
        !$omp end parallel do
        ! Assemble class contributions
        call cavgs%zero_set(.true.)
        call cavgs4reade%new_stack(ldim_crop(1:2), ncls)
        call cavgs4reado%new_stack(ldim_crop(1:2), ncls)
        if( L_BENCH_GLOB )then
            rt_init = toc(t_init)
            rt_io   = 0.
            rt_sum  = 0.
        endif
        do ipart=1,p_ptr%nparts
            if( L_BENCH_GLOB ) t_io = tic()
            ! filenames
            cae = 'cavgs_even_part'    //int2str_pad(ipart,p_ptr%numlen)//MRC_EXT
            cao = 'cavgs_odd_part'     //int2str_pad(ipart,p_ptr%numlen)//MRC_EXT
            cte = 'ctfsqsums_even_part'//int2str_pad(ipart,p_ptr%numlen)//MRC_EXT
            cto = 'ctfsqsums_odd_part' //int2str_pad(ipart,p_ptr%numlen)//MRC_EXT
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
        call cavger_restore_cavgs(p_ptr%frcs)
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
        if( allocated(eo_pops) ) deallocate(eo_pops)
    end subroutine cavger_kill

    !>  \brief submodule private destructor utility
    subroutine dealloc_cavgs
        call dealloc_imgarr(cavgs_even)
        call dealloc_imgarr(cavgs_odd)
        call dealloc_imgarr(cavgs_merged)
        call cavgs%kill_set
        ncls    = 0
        ldim    = 0; ldim_crop   = 0
        ldim_pd = 0; ldim_croppd = 0
        smpd    = 0.;smpd_crop   = 0.
        l_alloc_read_cavgs = .true.
    end subroutine dealloc_cavgs

    ! PUBLIC UTILITIES

    module subroutine transform_ptcls( params, build, spproj, oritype, icls, timgs, pinds, phflip, cavg, imgs_ori)
        use simple_sp_project,          only: sp_project
        use simple_strategy2D3D_common, only: discrete_read_imgbatch, prepimgbatch
        use simple_memoize_ft_maps
        class(parameters),                  intent(in)    :: params
        class(builder),                     target, intent(inout) :: build
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
        b_ptr => build
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
            call timgs(i)%new([p_ptr%box,p_ptr%box,1],p_ptr%smpd, wthreads=.false.)
        enddo
        if( l_imgs )then
            allocate(imgs_ori(pop))
            do i = 1,size(imgs_ori)
                call imgs_ori(i)%new([p_ptr%box,p_ptr%box,1],p_ptr%smpd, wthreads=.false.)
            enddo
        endif
        ! interpolation variables
        kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        allocate(kbw(wdim,wdim),source=0.)
        ! temporary objects
        call prepimgbatch(p_ptr, b_ptr, pop)
        do ithr = 1, nthr_glob
            call  img(ithr)%new([p_ptr%boxpd,p_ptr%boxpd,1],p_ptr%smpd, wthreads=.false.)
            call timg(ithr)%new([p_ptr%boxpd,p_ptr%boxpd,1],p_ptr%smpd, wthreads=.false.)
        end do
        call memoize_ft_maps(img(1)%get_ldim(), img(1)%get_smpd())
        logi_lims      = img(1)%loop_lims(2)
        cyc_lims       = img(1)%loop_lims(3)
        cyc_limsR(:,1) = cyc_lims(1,:)
        cyc_limsR(:,2) = cyc_lims(2,:)
        call discrete_read_imgbatch(p_ptr, b_ptr, pop, pinds(:), [1,pop])
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
            call b_ptr%imgbatch(i)%norm_noise_taper_edge_pad_fft(b_ptr%lmsk,img(ithr))
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

end submodule simple_classaverager_restore
