!@descr: for applying CTF to images
submodule (simple_image) simple_image_ctf
use simple_memoize_ft_maps
use simple_math_ctf,      only: ft_map_ctf_kernel
use simple_parameters,    only: parameters
use simple_euclid_sigma2, only: eucl_sigma2_glob
implicit none
#include "simple_local_flags.inc"

contains
 
    !>  \brief  is for generating an image of CTF
    module subroutine ctf2img( self, tfun, dfx_in, dfy_in, angast_in )
        class(image), intent(inout) :: self
        class(ctf),   intent(inout) :: tfun
        real,         intent(in)    :: dfx_in, dfy_in, angast_in
        type(ctfvars) :: ctfvals ! CTF derived variables
        integer :: h,k,physh,physk
        real    :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs, tval
        ! flag image as FT
        call self%set_ft(.true.)
        ! init
        call tfun%init(dfx_in, dfy_in, angast_in) ! conversions
        ctfvals         = tfun%get_ctfvars()
        wl              = ctfvals%wl
        half_wl2_cs     = 0.5 * wl * wl * ctfvals%cs
        sum_df          = ctfvals%dfx + ctfvals%dfy
        diff_df         = ctfvals%dfx - ctfvals%dfy
        angast          = ctfvals%angast
        amp_contr_const = ctfvals%amp_contr_const
        do h = ft_map_lims(1,1), ft_map_lims(1,2)
            do k = ft_map_lims(2,1), ft_map_lims(2,2)
                ! calculate CTF
                tval = ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)
                ! set cmat
                physh = ft_map_phys_addrh(h,k)
                physk = ft_map_phys_addrk(h,k)
                call self%set_cmat_at(physh,physk,1, cmplx(tval,0.))
            end do
        end do
    end subroutine ctf2img

    !>  \brief  is for applying CTF to an image
    module subroutine apply_ctf_wpad( self, tfun, dfx, mode, dfy, angast, bfac )
        class(image),     intent(inout) :: self   !< instance
        class(ctf),       intent(inout) :: tfun   !< CTF object
        real,             intent(in)    :: dfx    !< defocus x-axis
        character(len=*), intent(in)    :: mode   !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: dfy    !< defocus y-axis
        real, optional,   intent(in)    :: angast !< angle of astigmatism
        real, optional,   intent(in)    :: bfac   !< bfactor
        integer         :: ldim(3), ldim_pd(3)
        type(ctfparams) :: ctfparms 
        type(image)     :: img_pd
        ldim = self%get_ldim()
        if( self%is_3d() )then
            write(logfhandle,*) 'ldim: ', ldim
            THROW_HARD('Only 4 2D images; apply_ctf_wpad')
        endif
        ctfparms%smpd   = self%smpd
        ctfparms%dfx    = dfx
        ! defaults
        ctfparms%dfy    = ctfparms%dfx
        ctfparms%angast = 0.
        ! optionals
        if(present(dfy))    ctfparms%dfy    = dfy
        if(present(angast)) ctfparms%angast = angast
        if( self%is_ft() )then
            call self%apply_ctf(tfun, mode, ctfparms)
        else
            ldim_pd(1:2) = 2*ldim(1:2)
            ldim_pd(3)   = 1
            call img_pd%new(ldim_pd, self%smpd)
            call self%pad_mirr(img_pd)
            call img_pd%fft()
            call img_pd%apply_ctf(tfun, mode, ctfparms)
            call img_pd%ifft()
            call img_pd%clip(self)
            call img_pd%kill()
        endif
        if( present(bfac) ) call self%apply_bfac(bfac)
    end subroutine apply_ctf_wpad

    !>  \brief  is for optimised serial application of CTF
    !!          modes: abs, ctf, flip, flipneg, neg, square
    module subroutine apply_ctf( self, tfun, mode, ctfparms )
        class(image),     intent(inout) :: self     !< instance
        class(ctf),       intent(inout) :: tfun     !< CTF object
        character(len=*), intent(in)    :: mode     !< abs, ctf, flip, flipneg, neg, square
        type(ctfparams),  intent(in)    :: ctfparms !< CTF parameters
        type(ctfvars) :: ctfvals
        real, parameter :: ZERO = 0., ONE = 1.0, MIN_SQUARE = 0.001
        integer :: h,k,physh,physk
        real    :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs, tval, t
        logical :: is_abs,is_ctf,is_flip,is_flipneg,is_neg,is_square
        ! Convert mode string to logical flags (SIMD-compatible)
        is_abs     = (mode == 'abs')
        is_ctf     = (mode == 'ctf')
        is_flip    = (mode == 'flip')
        is_flipneg = (mode == 'flipneg')
        is_neg     = (mode == 'neg')
        is_square  = (mode == 'square')
        ! initialize
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast) ! conversions
        ctfvals         = tfun%get_ctfvars()
        wl              = ctfvals%wl
        half_wl2_cs     = 0.5 * wl * wl * ctfvals%cs
        sum_df          = ctfvals%dfx + ctfvals%dfy
        diff_df         = ctfvals%dfx - ctfvals%dfy
        angast          = ctfvals%angast
        amp_contr_const = ctfvals%amp_contr_const
        do h = ft_map_lims(1,1), ft_map_lims(1,2)
            do k = ft_map_lims(2,1), ft_map_lims(2,2)
                ! calculate CTF
                tval = ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)
                ! SIMD-compatible mode handling using masked operations
                t = merge(       abs(tval),                     ZERO, is_abs)     + &
                    merge(           tval,                      ZERO, is_ctf)     + &
                    merge( sign(ONE, tval),                     ZERO, is_flip)    + &
                    merge(-sign(ONE, tval),                     ZERO, is_flipneg) + &
                    merge(          -tval,                      ZERO, is_neg)     + &
                    merge(min(ONE, max(tval*tval, MIN_SQUARE)), ZERO, is_square)
                ! Multiply (this is the key operation)
                physh = ft_map_phys_addrh(h,k)
                physk = ft_map_phys_addrk(h,k)
                self%cmat(physh,physk,1) = t * self%cmat(physh,physk,1)
            end do
        end do
    end subroutine apply_ctf

    module subroutine gen_fplane4rec( self, smpd_crop, ctfparms, shift, l_ml_reg, iptcl, fplane )
        use simple_math_ft, only: upsample_sigma2
        class(image),      intent(inout) :: self
        real,              intent(in)    :: smpd_crop
        class(ctfparams),  intent(in)    :: ctfparms
        real,              intent(in)    :: shift(2)
        logical,           intent(in)    :: l_ml_reg
        integer,           intent(in)    :: iptcl
        type(fplane_type), intent(out)   :: fplane
        type(ctf)                :: tfun
        type(ctfvars)            :: ctfvals
        real, allocatable        :: sigma2_noise_tmp(:), sigma2_noise(:) !< Noise power spectrum for ML regularization
        complex(c_float_complex) :: c, w1, w2, ph0, ph_h, ph_k
        real(dp)                 :: pshift(2)
        type(ftiter) :: fiterator
        ! CTF kernel scalars (precomputed)
        real    :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs, ker, tval, tvalsq
        integer :: physh, physk, sigma2_kfromto(2), h, k, shell, hmin, hmax, kmin, kmax,  frlims_crop(3,2), sigma_nyq
        integer :: box_croppd, box_crop
        logical :: l_ctf, l_flip
        ! Shell LUT: shell = nint(sqrt(r2)) via lookup table
        integer, allocatable :: shell_lut(:)
        integer :: max_r2, r2, abs_hmax, abs_kmax
        ! shift is with respect to the original image dimension
        fplane%shconst = self%get_shconst()
        ! -----------------------
        ! setup the Fourier plane (PADDED/CROPPED)
        ! -----------------------
        box_croppd = self%ldim(1)
        box_crop   = box_croppd / OSMPL_PAD_FAC
        ! cropped Fourier limits & dimensions need to be congruent with how the reconstroctor_eo objects are set up
        call fiterator%new([box_croppd, box_croppd, 1], smpd_crop)
        ! use redundant logical limits (mode 3) so negative h exists for Friedel access (-h,-k)
        fplane%frlims = fiterator%loop_lims(3)
        fplane%nyq    = fiterator%get_lfny(1)
        call fiterator%new([box_crop, box_crop, 1], smpd_crop)
        sigma_nyq     = fiterator%get_lfny(1)
        ! matrices
        if( allocated(fplane%cmplx_plane) ) deallocate(fplane%cmplx_plane)
        if( allocated(fplane%ctfsq_plane) ) deallocate(fplane%ctfsq_plane)
        ! allocate only k<=0 due to Friedel symmetry; the rest will be filled in by conjugation
        hmin = fplane%frlims(1,1); hmax = fplane%frlims(1,2)
        kmin = fplane%frlims(2,1); kmax = fplane%frlims(2,2)
        allocate(fplane%cmplx_plane(hmin:hmax, kmin:0), &
        fplane%ctfsq_plane(hmin:hmax, kmin:0))
        fplane%cmplx_plane = cmplx(0.,0.)
        fplane%ctfsq_plane = 0.
        ! -----------------------
        ! CTF flags
        ! -----------------------
        l_ctf  = (ctfparms%ctfflag /= CTFFLAG_NO)
        l_flip = .false.
        if (l_ctf) then
            l_flip = (ctfparms%ctfflag == CTFFLAG_FLIP)
            tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
            call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
            ! --- optimized kernel scalars ---
            ctfvals         = tfun%get_ctfvars()
            wl              = ctfvals%wl
            half_wl2_cs     = 0.5 * wl * wl * ctfvals%cs
            sum_df          = ctfvals%dfx + ctfvals%dfy
            diff_df         = ctfvals%dfx - ctfvals%dfy
            angast          = ctfvals%angast
            amp_contr_const = ctfvals%amp_contr_const
        end if
        ! -----------------------
        ! ML regularization noise spectrum
        ! -----------------------
        if (l_ml_reg) then
            allocate(sigma2_noise_tmp(1:sigma_nyq), sigma2_noise(0:fplane%nyq), source=0.0)
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_noise_tmp(sigma2_kfromto(1):sigma_nyq) = eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):sigma_nyq, iptcl)
            call upsample_sigma2(sigma2_kfromto(1), sigma_nyq, sigma2_noise_tmp, fplane%nyq, sigma2_noise)
        end if
        ! -----------------------
        ! Shift phase recurrence
        ! -----------------------
        pshift = real(-shift * fplane%shconst(1:2), dp)
        w1     = cmplx( real(cos(pshift(1)), c_float), real(sin(pshift(1)), c_float), kind=c_float_complex )
        w2     = cmplx( real(cos(pshift(2)), c_float), real(sin(pshift(2)), c_float), kind=c_float_complex )
        ph0    = cmplx( real(cos(real(fplane%frlims(1,1),dp)*pshift(1)), c_float), &
                        real(sin(real(fplane%frlims(1,1),dp)*pshift(1)), c_float), kind=c_float_complex )
        ph_k   = cmplx( real(cos(real(fplane%frlims(2,1),dp)*pshift(2)), c_float), &
                        real(sin(real(fplane%frlims(2,1),dp)*pshift(2)), c_float), kind=c_float_complex )
        ! -----------------------
        ! Precompute shell LUT to avoid sqrt in inner loops
        ! r2 = h*h + k*k, shell = nint(sqrt(r2))
        ! -----------------------
        hmin = fplane%frlims(1,1)
        hmax = fplane%frlims(1,2)
        kmin = fplane%frlims(2,1)
        kmax = fplane%frlims(2,2)
        abs_hmax = max(abs(hmin), abs(hmax))
        abs_kmax = max(abs(kmin), abs(kmax))
        max_r2   = abs_hmax*abs_hmax + abs_kmax*abs_kmax
        allocate(shell_lut(0:max_r2))
        do r2 = 0, max_r2
            shell_lut(r2) = nint(sqrt(real(r2)))
        end do
        ! ============================================================
        ! Fill k in [kmin .. 0] explicitly (ONLY STORED REGION)
        ! ============================================================
        do k = kmin, 0
            ph_h = ph0
            do h = hmin, hmax
                r2    = h*h + k*k
                shell = shell_lut(r2)
                if (shell > fplane%nyq) then
                    c      = cmplx(0.0_c_float, 0.0_c_float, kind=c_float_complex)
                    tvalsq = 0.0
                else
                    ! Retrieve Fourier component & apply shift phase
                    physh = ft_map_phys_addrh(h,k)
                    physk = ft_map_phys_addrk(h,k)
                    c     = merge(conjg(self%cmat(physh,physk,1)), self%cmat(physh,physk,1), h < 0) * (ph_k * ph_h)
                    ! CTF (optimized kernel; no phase-plate support)
                    if (l_ctf) then
                        ker = ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)
                        if (l_flip) then
                            ! Consistent with norm_noise_fft_clip_shift_ctf_flip:
                            ! flip sign only; ctfsq is 1 (then ML weighting may scale it)
                            tval   = sign(1.0, ker)
                            tvalsq = 1.0
                            c      = tval * c
                        else
                            tval   = ker
                            tvalsq = tval * tval
                            c      = tval * c
                        end if
                    else
                        tvalsq = 1.0
                    end if
                    ! sigma2 weighting (unchanged semantics)
                    if (l_ml_reg) then
                        c      = c      / sigma2_noise(shell)
                        tvalsq = tvalsq / sigma2_noise(shell)
                    end if
                end if
                fplane%cmplx_plane(h,k) = c
                fplane%ctfsq_plane(h,k) = tvalsq
                ph_h = ph_h * w1
            end do
            ph_k = ph_k * w2
        end do
        if (allocated(shell_lut)       ) deallocate(shell_lut)
        if (allocated(sigma2_noise_tmp)) deallocate(sigma2_noise_tmp)
        if (allocated(sigma2_noise)    ) deallocate(sigma2_noise)
    end subroutine gen_fplane4rec

    ! Calculate Joe's Ice Fraction Score for images, not intended for micrographs
    module subroutine calc_ice_frac( self, tfun, ctfparms, score )
        class(image),     intent(in)    :: self
        class(ctf),       intent(inout) :: tfun
        class(ctfparams), intent(in)    :: ctfparms
        real,             intent(out)   :: score
        real, allocatable :: res(:)
        real    :: phshift, start_freq, end_freq, ice_avg, band_avg, mag_max, mag, g
        integer :: lims(3,2), box, sh, h,k, start_find, end_find, ice_maxind, hmax, kmax, cnt
        score = 0.
        if( abs(self%smpd-ctfparms%smpd) > 1.d-4) THROW_HARD('Inconsistent SMPD; calc_ice_frac')
        if( self%smpd > (ICE_BAND1/2.) ) return
        if( .not.self%is_ft() ) THROW_HARD('Image input must be in the Fourier domain!; calc_ice_frac')
        box = self%get_box()
        res = get_resarr(box, self%smpd)
        lims = self%loop_lims(2)
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        phshift    = 0.
        start_freq = sqrt(tfun%SpaFreqSqAtNthZero(1, phshift, deg2rad(ctfparms%angast)))
        end_freq   = sqrt(tfun%SpaFreqSqAtNthZero(2, phshift, deg2rad(ctfparms%angast)))
        ice_maxind = get_find_at_res(res, ICE_BAND1)
        start_find = max(1,     ice_maxind - 3)
        end_find   = min(box/2, ice_maxind + 3)
        hmax     = -1
        kmax     = -1
        mag_max  = -1.
        band_avg = 0.
        cnt      = 0
        do k = lims(2,1),lims(2,2)
            do h = lims(1,1),lims(1,2)
                sh = nint(hyp(h,k))
                g  = real(sh) / real(box)
                if( g > start_freq .and. g < end_freq )then
                    band_avg = band_avg + csq_fast(self%get_fcomp2D(h,k))
                    cnt      = cnt + 1
                else if( sh >= start_find .and. sh < end_find )then
                    mag = csq_fast(self%get_fcomp2D(h,k))
                    if( mag > mag_max )then
                        hmax = h
                        kmax = k
                        mag_max = mag
                    endif
                endif
            end do
        end do
        if( cnt < 1 ) return
        band_avg = band_avg / real(cnt)
        ice_avg  = 0.
        do k = kmax-1,kmax+1
            do h = hmax-1,hmax+1
                ice_avg = ice_avg + csq_fast(self%get_fcomp2D(h,k))
            enddo
        enddo
        ice_avg = ice_avg / 9.
        score   = ice_avg / band_avg
        if( score > 0.5 ) score = (ice_avg - 0.5*band_avg) / band_avg
    end subroutine calc_ice_frac

    ! KEEP SERIAL
    module subroutine ctf_dens_correct( self_sum, self_rho )
        class(image), intent(inout) :: self_sum
        class(image), intent(inout) :: self_rho
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh    = nint(hyp(h,k,l))
                    phys  = self_sum%comp_addr_phys(h,k,l)
                    denom = real(self_rho%cmat(phys(1),phys(2),phys(3)))
                    if(sh <= nyq .and. abs(denom) > 1.e-10 )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3)) / denom
                    else
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct

    module subroutine ctf_dens_correct_wiener( self_sum, self_rho, ssnr )
        class(image), intent(inout) :: self_sum
        class(image), intent(in)    :: self_rho
        real,         intent(in)    :: ssnr(:)
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(h,k,l))
                    phys = self_sum%comp_addr_phys(h,k,l)
                    if(sh > nyq )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh==0 )then
                            denom = real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        else
                            denom = ssnr(sh)*real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = ssnr(sh)*self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        endif
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct_wiener

end submodule simple_image_ctf
