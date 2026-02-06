!@descr: for applying CTF to images
submodule (simple_image) simple_image_ctf
use simple_memoize_ft_maps
use simple_math_ctf,      only: ft_map_ctf_kernel
use simple_parameters,    only: params_glob
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

    module subroutine gen_fplane4rec( self, ctfparms, shift, iptcl, fplane )
        class(image),      intent(inout) :: self
        class(ctfparams),  intent(in)    :: ctfparms
        real,              intent(in)    :: shift(2)
        integer,           intent(in)    :: iptcl
        type(fplane_type), intent(out)   :: fplane
        type(ctf)                :: tfun
        type(ctfvars) :: ctfvals
        real, allocatable        :: sigma2_noise(:) !< Noise power spectrum for ML regularization
        complex(c_float_complex) :: c, w1, w2, ph0, ph_h, ph_k
        real(dp)                 :: pshift(2)
        type(ftiter) :: fiterator
        ! CTF kernel scalars (precomputed)
        real    :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs, ker, tval, tvalsq
        integer :: physh, physk, sigma2_kfromto(2), h, k, shell, hmin, hmax, kmin, kmax 

        integer :: frlims_crop(3,2),nyq_crop

        logical :: l_ctf, l_flip
        ! Shell LUT: shell = nint(sqrt(r2)) via lookup table
        integer, allocatable :: shell_lut(:)
        integer :: max_r2, r2, abs_hmax, abs_kmax

        ! shift is with respect to the original image dimension
        fplane%shconst = self%get_shconst()
        
        ! -----------------------
        ! setup the Fourier plane
        ! -----------------------

        ! Fourier limits & dimensions should not be the same as the image because this bugs out
        
        ! fplane%frlims  = self%fit%loop_lims(3)
        ! fplane%nyq     = self%fit%get_lfny(1)
        ! print *, 'frclims1_fplane: ', fplane%frlims(1,1), fplane%frlims(1,2)
        ! print *, 'frclims2_fplane: ', fplane%frlims(1,1), fplane%frlims(1,2)
        ! print *, 'frclims3_fplane: ', fplane%frlims(1,1), fplane%frlims(1,2)
        ! print *, 'nyq_fplane:      ', fplane%nyq

        ! cropped Fourier limits & dimensions need to be congruent with how the reconstroctor_eo objects are set up
        ! I guess this will have to change with the new padding policy, but for now it runs
        call fiterator%new([params_glob%box_crop, params_glob%box_crop, 1], params_glob%smpd_crop)
        fplane%frlims = fiterator%loop_lims(3)
        fplane%nyq    = fiterator%get_lfny(1)

        ! matrices
        if( allocated(fplane%cmplx_plane) ) deallocate(fplane%cmplx_plane)
        if( allocated(fplane%ctfsq_plane) ) deallocate(fplane%ctfsq_plane)
        allocate(fplane%cmplx_plane(fplane%frlims(1,1):fplane%frlims(1,2),fplane%frlims(2,1):fplane%frlims(2,2)),&
        &fplane%ctfsq_plane(fplane%frlims(1,1):fplane%frlims(1,2),fplane%frlims(2,1):fplane%frlims(2,2)))
        fplane%cmplx_plane  = cmplx(0.,0.)
        fplane%ctfsq_plane  = 0.
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
        if (params_glob%l_ml_reg) then
            allocate(sigma2_noise(1:fplane%nyq), source=0.0)
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
            sigma2_noise(sigma2_kfromto(1):fplane%nyq) = &
                eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):fplane%nyq, iptcl)
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
        ! Fill k in [kmin .. 0] explicitly
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
                    if (params_glob%l_ml_reg) then
                        if (shell >= sigma2_kfromto(1)) then
                            c      = c      / sigma2_noise(shell)
                            tvalsq = tvalsq / sigma2_noise(shell)
                        else
                            c      = c      / sigma2_noise(sigma2_kfromto(1))
                            tvalsq = tvalsq / sigma2_noise(sigma2_kfromto(1))
                        end if
                    end if
                end if
                fplane%cmplx_plane(h,k) = c
                fplane%ctfsq_plane(h,k) = tvalsq
                ph_h = ph_h * w1
            end do
            ph_k = ph_k * w2
        end do
        ! ============================================================
        ! Fill k in [1 .. kmax] by Friedel symmetry only
        ! No sqrt/shell needed: if (-h,-k) was out of nyq, it was already zeroed.
        ! ============================================================
        do k = 1, kmax
            do h = hmin, hmax
                fplane%cmplx_plane(h,k) = conjg(fplane%cmplx_plane(-h,-k))
                fplane%ctfsq_plane(h,k) =       fplane%ctfsq_plane(-h,-k)
            end do
        end do
        if (allocated(shell_lut))    deallocate(shell_lut)
        if (allocated(sigma2_noise)) deallocate(sigma2_noise)
    end subroutine gen_fplane4rec

    !> Fused routine:
    !!   (1) taper edges in real space (in-place on self)
    !!   (2) compute edge_mean on tapered border
    !!   (3) noise-normalize (masked) in real space (mean/var on tapered self)
    !!   (4) pad + fftshift while copying into a work image, filling pad with background
    !!       that is consistent with normalization: bg_fill = (edge_mean-mean)*invstd
    !!   (5) FFT (r2c) on work
    !!   (6) generate reconstruction fplane (shift + CTF-kernel + ML weighting) directly
    !!       from work%cmat, filling k>0 via Friedel symmetry (no extra sqrt)
    !!
    !! Notes:
    !!  - Requires a pre-planned padded FFT work image (work) with valid plan_fwd.
    !!    This needs to be a heap for parallel adressing in OpenMP regions (big enouhg to justify)
    !!  - Uses the optimized CTF kernel (ft_map_ctf_kernel)
    !!  - FLIP semantics match norm_noise_fft_clip_shift_ctf_flip:
    !!      * FLIP: multiply by sign(1,kernel), and ctfsq=1
    !!      * non-FLIP: multiply by kernel, and ctfsq=kernel^2
    module subroutine norm_noise_taper_edge_pad_fft_gen_fplane4rec( self, lmsk, work, ctfparms, shift, iptcl, fplane )
        class(image),      intent(inout) :: self
        logical,           intent(in)    :: lmsk(self%ldim(1), self%ldim(2), self%ldim(3))
        class(image),      intent(inout) :: work          !< padded FFT workspace (rmat+cmat+plan_fwd)
        class(ctfparams),  intent(in)    :: ctfparms
        real,              intent(in)    :: shift(2)
        integer,           intent(in)    :: iptcl
        type(fplane_type), intent(out)   :: fplane
        ! ------------------------------------------------------------
        ! Part A: taper + edge_mean + noise norm + pad/fftshift into work + FFT
        ! (largely matches norm_noise_taper_edge_pad_fft)
        ! ------------------------------------------------------------
        integer       :: n1, n2, n1o, n2o, h1o, h2o
        integer       :: i, j, npix, starts(3), x0, y0, xo, yo, winsz
        ! Normalization stats in DP (scalars only)
        real(dp)      :: mean_dp, var_dp, sum_dp, sum_sq_dp, rnpix, invstd_dp
        real(c_float) :: mean_sp, invstd_sp, x_sp
        logical       :: do_norm
        ! taper vars (square, 2D, n3=1)
        integer       :: n, wtap, wavg, wbg, j1, j2, border_count
        real(dp)      :: bg_fill_dp, border_sum_dp, edge_mean_dp
        real(c_float) :: alpha_sp, inv_wtap_sp, inv_n_sp
        ! Edge profiles remain SP (cheap+stable)
        real(c_float) :: x_start(self%ldim(2)), x_stop(self%ldim(2))
        real(c_float) :: x_smooth_start(self%ldim(2)), x_smooth_stop(self%ldim(2))
        real(c_float) :: y_start(self%ldim(1)), y_stop(self%ldim(1))
        real(c_float) :: y_smooth_start(self%ldim(1)), y_smooth_stop(self%ldim(1))
        real(c_float) :: edge_avg(self%ldim(1))
        real(c_float_complex) :: c
        ! ------------------------------------------------------------
        ! Part B: build fplane from work%cmat with shift + CTF + ML reg
        ! ------------------------------------------------------------
        type(ctf)                :: tfun
        type(ctfvars)            :: ctfvals
        real, allocatable        :: sigma2_noise(:)
        integer                  :: sigma2_kfromto(2)
        ! phase recurrence for Fourier-domain shift
        complex(c_float_complex) :: w1, w2, ph0, ph_h, ph_k
        real(dp)                 :: pshift(2)
        ! CTF kernel scalars (precomputed)
        real    :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs,  ker, tval, tvalsq
        integer :: h, k, hp, physh, physk, shell, hmin, hmax, kmin, kmax
        logical :: l_ctf, l_flip
        ! shell LUT to avoid sqrt in hot loop
        integer, allocatable :: shell_lut(:)
        integer :: max_r2, r2, abs_hmax, abs_kmax
        ! ============================================================
        ! A0) dimensions
        ! ============================================================
        winsz = nint(COSMSKHALFWIDTH)
        n1    = self%ldim(1)
        n2    = self%ldim(2)
        n     = n1                         ! assume square: n1==n2
        n1o   = work%ldim(1)
        n2o   = work%ldim(2)
        h1o   = n1o/2
        h2o   = n2o/2
        ! ============================================================
        ! A1) taper (in-place on self) + edge_mean
        ! ============================================================
        wtap = min(max(winsz, 0), n/2)
        wbg  = max(1, winsz/2)
        wbg  = min(wbg, n/2)
        if (wtap > 0) then
            wavg = min(max(1, wtap/8), n)
            inv_n_sp = 1.0_c_float / real(wavg, c_float)
            ! ---- Pass 1: X edges ----
            do j = 1, n
                x_start(j) = sum(self%rmat(1:wavg,     j, 1)) * inv_n_sp
                x_stop (j) = sum(self%rmat(n-wavg+1:n, j, 1)) * inv_n_sp
            end do
            edge_avg(:) = 0.5_c_float * (x_start(:) + x_stop(:))
            x_start(:)  = x_start(:) - edge_avg(:)
            x_stop(:)   = x_stop(:)  - edge_avg(:)
            do j = 1, n
                j1 = max(1, j-1)
                j2 = min(n, j+1)
                inv_n_sp = 1.0_c_float / real(j2 - j1 + 1, c_float)
                x_smooth_start(j) = sum(x_start(j1:j2)) * inv_n_sp
                x_smooth_stop (j) = sum(x_stop (j1:j2)) * inv_n_sp
            end do
            inv_wtap_sp = 1.0_c_float / real(wtap, c_float)
            do i = 1, wtap
                alpha_sp = real(wtap - i + 1, c_float) * inv_wtap_sp
                self%rmat(i, 1:n, 1) = self%rmat(i, 1:n, 1) - x_smooth_start(:) * alpha_sp
            end do
            do i = n - wtap + 1, n
                alpha_sp = real(wtap + i - n, c_float) * inv_wtap_sp
                self%rmat(i, 1:n, 1) = self%rmat(i, 1:n, 1) - x_smooth_stop(:) * alpha_sp
            end do
            ! ---- Pass 2: Y edges ----
            inv_n_sp = 1.0_c_float / real(wavg, c_float)
            do i = 1, n
                y_start(i) = sum(self%rmat(i, 1:wavg,     1)) * inv_n_sp
                y_stop (i) = sum(self%rmat(i, n-wavg+1:n, 1)) * inv_n_sp
            end do
            edge_avg(:) = 0.5_c_float * (y_start(:) + y_stop(:))
            y_start(:)  = y_start(:) - edge_avg(:)
            y_stop(:)   = y_stop(:)  - edge_avg(:)
            do i = 1, n
                j1 = max(1, i-1)
                j2 = min(n, i+1)
                inv_n_sp = 1.0_c_float / real(j2 - j1 + 1, c_float)
                y_smooth_start(i) = sum(y_start(j1:j2)) * inv_n_sp
                y_smooth_stop (i) = sum(y_stop (j1:j2)) * inv_n_sp
            end do
            do j = 1, wtap
                alpha_sp = real(wtap - j + 1, c_float) * inv_wtap_sp
                self%rmat(1:n, j, 1) = self%rmat(1:n, j, 1) - y_smooth_start(:) * alpha_sp
            end do
            do j = n - wtap + 1, n
                alpha_sp = real(wtap + j - n, c_float) * inv_wtap_sp
                self%rmat(1:n, j, 1) = self%rmat(1:n, j, 1) - y_smooth_stop(:) * alpha_sp
            end do
        end if
        ! edge_mean on tapered border (DP scalar accumulation)
        border_sum_dp = 0.0_dp
        border_sum_dp = border_sum_dp + sum( real(self%rmat(1:wbg,       1:n,        1), dp) )
        border_sum_dp = border_sum_dp + sum( real(self%rmat(n-wbg+1:n,   1:n,        1), dp) )
        border_sum_dp = border_sum_dp + sum( real(self%rmat(wbg+1:n-wbg, 1:wbg,      1), dp) )
        border_sum_dp = border_sum_dp + sum( real(self%rmat(wbg+1:n-wbg, n-wbg+1:n,  1), dp) )
        border_count  = 4*wbg*n - 4*wbg*wbg
        edge_mean_dp  = border_sum_dp / real(border_count, dp)
        ! ============================================================
        ! A2) noise normalization stats on tapered self
        ! ============================================================
        sum_dp    = 0.0_dp
        sum_sq_dp = 0.0_dp
        npix      = 0
        do i = 1, n1
            do j = 1, n2
                if (.not. lmsk(i,j,1)) then
                    npix = npix + 1
                    x_sp = self%rmat(i,j,1)
                    sum_dp    = sum_dp    + real(x_sp, dp)
                    sum_sq_dp = sum_sq_dp + real(x_sp, dp) * real(x_sp, dp)
                end if
            end do
        end do
        mean_dp   = 0.0_dp
        var_dp    = 0.0_dp
        invstd_dp = 1.0_dp
        do_norm   = .false.
        if (npix > 1) then
            rnpix   = real(npix, dp)
            mean_dp = sum_dp / rnpix
            var_dp  = (sum_sq_dp - sum_dp * sum_dp / rnpix) / real(npix - 1, dp)
            if (is_a_number(real(var_dp, kind=kind(1.0))) .and. var_dp > 0.0_dp) then
                invstd_dp = 1.0_dp / sqrt(var_dp)
                do_norm = (abs(real(mean_dp, kind=kind(1.0))) > TINY .or. invstd_dp /= 1.0_dp)
            end if
        end if
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = real(invstd_dp, c_float)
        ! ============================================================
        ! A3) pad + fftshift while copying into work, fill pad with bg
        !     bg is consistent with normalization if do_norm
        ! ============================================================
        starts = (work%ldim - self%ldim) / 2 + 1
        if (do_norm) then
            bg_fill_dp = (edge_mean_dp - mean_dp) * invstd_dp
        else
            bg_fill_dp = edge_mean_dp
        end if
        work%ft   = .false.
        work%rmat = real(bg_fill_dp, c_float)
        if (do_norm) then
            do j = 1, n2
                y0 = starts(2) + j - 1
                yo = modulo((y0 - 1) + h2o, n2o) + 1
                do i = 1, n1
                    x0 = starts(1) + i - 1
                    xo = modulo((x0 - 1) + h1o, n1o) + 1
                    work%rmat(xo, yo, 1) = (self%rmat(i,j,1) - mean_sp) * invstd_sp
                end do
            end do
        else
            do j = 1, n2
                y0 = starts(2) + j - 1
                yo = modulo((y0 - 1) + h2o, n2o) + 1
                do i = 1, n1
                    x0 = starts(1) + i - 1
                    xo = modulo((x0 - 1) + h1o, n1o) + 1
                    work%rmat(xo, yo, 1) = self%rmat(i,j,1)
                end do
            end do
        end if
        ! ============================================================
        ! A4) FFT (FFTW r2c) + scale (keep consistent with your pad_fft)
        ! ============================================================
        call fftwf_execute_dft_r2c(work%plan_fwd, work%rmat, work%cmat)
        work%cmat = work%cmat * (1.0_c_float / real(n1o*n2o, c_float))
        work%ft   = .true.
        ! ------------------------------------------------------------
        ! Part B: build fplane_type from work%cmat
        ! ------------------------------------------------------------
        ! ============================================================
        ! B0) fplane metadata + allocate planes
        ! Use work's Fourier limits and Nyquist (consistent with padded FFT grid)
        ! ============================================================
        fplane%shconst = work%get_shconst()
        fplane%frlims  = work%fit%loop_lims(3)
        fplane%nyq     = work%fit%get_lfny(1)
        if (allocated(fplane%cmplx_plane)) deallocate(fplane%cmplx_plane)
        if (allocated(fplane%ctfsq_plane)) deallocate(fplane%ctfsq_plane)
        allocate( fplane%cmplx_plane(fplane%frlims(1,1):fplane%frlims(1,2), fplane%frlims(2,1):fplane%frlims(2,2)), &
                fplane%ctfsq_plane(fplane%frlims(1,1):fplane%frlims(1,2), fplane%frlims(2,1):fplane%frlims(2,2)) )
        fplane%cmplx_plane = cmplx(0.0_c_float, 0.0_c_float, kind=c_float_complex)
        fplane%ctfsq_plane = 0.0
        hmin = fplane%frlims(1,1)
        hmax = fplane%frlims(1,2)
        kmin = fplane%frlims(2,1)
        kmax = fplane%frlims(2,2)
        ! ============================================================
        ! B1) CTF kernel scalars (no phase-plate support)
        ! ============================================================
        l_ctf  = (ctfparms%ctfflag /= CTFFLAG_NO)
        l_flip = .false.
        if (l_ctf) then
            l_flip = (ctfparms%ctfflag == CTFFLAG_FLIP)
            tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
            call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
            ctfvals         = tfun%get_ctfvars()
            wl              = ctfvals%wl
            half_wl2_cs     = 0.5 * wl * wl * ctfvals%cs
            sum_df          = ctfvals%dfx + ctfvals%dfy
            diff_df         = ctfvals%dfx - ctfvals%dfy
            angast          = ctfvals%angast
            amp_contr_const = ctfvals%amp_contr_const
        end if
        ! ============================================================
        ! B2) ML regularization noise spectrum
        ! ============================================================
        if (params_glob%l_ml_reg) then
            allocate(sigma2_noise(1:fplane%nyq), source=0.0)
            sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise, 1)
            sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise, 1)
            sigma2_noise(sigma2_kfromto(1):fplane%nyq) = &
                eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):fplane%nyq, iptcl)
        end if
        ! ============================================================
        ! B3) Shift phase recurrence on the padded FFT grid
        ! ============================================================
        pshift = real(-shift * fplane%shconst(1:2), dp)
        w1 = cmplx( real(cos(pshift(1)), c_float), real(sin(pshift(1)), c_float), kind=c_float_complex )
        w2 = cmplx( real(cos(pshift(2)), c_float), real(sin(pshift(2)), c_float), kind=c_float_complex )
        ! starting phase for h: exp(i*hmin*pshift(1))
        ph0 = cmplx( real(cos(real(hmin,dp)*pshift(1)), c_float), &
                    real(sin(real(hmin,dp)*pshift(1)), c_float), kind=c_float_complex )
        ! starting phase for k: exp(i*kmin*pshift(2))
        ph_k = cmplx( real(cos(real(kmin,dp)*pshift(2)), c_float), &
                    real(sin(real(kmin,dp)*pshift(2)), c_float), kind=c_float_complex )
        ! ============================================================
        ! B4) Precompute shell LUT to avoid sqrt in hot loop
        ! ============================================================
        abs_hmax = max(abs(hmin), abs(hmax))
        abs_kmax = max(abs(kmin), abs(kmax))
        max_r2   = abs_hmax*abs_hmax + abs_kmax*abs_kmax
        allocate(shell_lut(0:max_r2))
        do r2 = 0, max_r2
            shell_lut(r2) = nint(sqrt(real(r2)))
        end do
        ! ============================================================
        ! B5) Fill k in [kmin .. 0] explicitly, k>0 by Friedel symmetry
        ! ============================================================
        do k = kmin, 0
            ph_h = ph0
            do h = hmin, hmax
                r2    = h*h + k*k
                shell = shell_lut(r2)
                if (shell > fplane%nyq) then
                    fplane%cmplx_plane(h,k) = cmplx(0.0_c_float, 0.0_c_float, kind=c_float_complex)
                    fplane%ctfsq_plane(h,k) = 0.0
                else
                    physh = ft_map_phys_addrh(h,k)
                    physk = ft_map_phys_addrk(h,k)
                    ! Pull coefficient from padded FFT work buffer; apply Friedel for h<0
                    c = merge(conjg(work%cmat(physh,physk,1)), work%cmat(physh,physk,1), h < 0) * (ph_k * ph_h)
                    ! CTF kernel (no phase-plate support)
                    if (l_ctf) then
                        ker = ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)
                        if (l_flip) then
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
                    ! ML weighting
                    if (params_glob%l_ml_reg) then
                        if (shell >= sigma2_kfromto(1)) then
                            c      = c      / sigma2_noise(shell)
                            tvalsq = tvalsq / sigma2_noise(shell)
                        else
                            c      = c      / sigma2_noise(sigma2_kfromto(1))
                            tvalsq = tvalsq / sigma2_noise(sigma2_kfromto(1))
                        end if
                    end if
                    fplane%cmplx_plane(h,k) = c
                    fplane%ctfsq_plane(h,k) = tvalsq
                end if
                ph_h = ph_h * w1
            end do
            ph_k = ph_k * w2
        end do
        do k = 1, kmax
            do h = hmin, hmax
                fplane%cmplx_plane(h,k) = conjg(fplane%cmplx_plane(-h,-k))
                fplane%ctfsq_plane(h,k) =       fplane%ctfsq_plane(-h,-k)
            end do
        end do
        if (allocated(shell_lut))    deallocate(shell_lut)
        if (allocated(sigma2_noise)) deallocate(sigma2_noise)
    end subroutine norm_noise_taper_edge_pad_fft_gen_fplane4rec

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
