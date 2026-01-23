submodule (simple_image) simple_image_fft
implicit none
#include "simple_local_flags.inc"

contains

    !===========================
    ! Core FFT routines
    !===========================

    module subroutine fwd_ft(self)
        class(image), intent(inout) :: self
        if( self%ft ) return
        call self%shift_phorig
        call fftwf_execute_dft_r2c(self%plan_fwd,self%rmat,self%cmat)
        ! now scale the values so that a ifft() of the output yields the
        ! original image back following FFTW
        self%cmat = self%cmat/real(product(self%ldim))
        self%ft = .true.
    end subroutine fwd_ft

    module subroutine bwd_ft( self )
        class(image), intent(inout) :: self
        if( self%ft )then
            call fftwf_execute_dft_c2r(self%plan_bwd,self%cmat,self%rmat)
            self%ft = .false.
            call self%shift_phorig
        endif
    end subroutine bwd_ft

    module subroutine fft_noshift( self )
        class(image), intent(inout) :: self
        if( self%ft ) return
        call fftwf_execute_dft_r2c(self%plan_fwd,self%rmat,self%cmat)
        self%cmat = self%cmat/real(product(self%ldim))
        self%ft = .true.
    end subroutine fft_noshift

    !===========================
    ! FT / image conversion
    !===========================

    module subroutine img2ft( self, img )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: img
        integer :: h,k,l,lims(3,2),logi(3),phys(3)
        integer :: xcnt,ycnt,zcnt
        if( .not.(self.eqdims.img) )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            write(logfhandle,*) 'img%ldim:  ', img%ldim
            THROW_HARD('non-equal dims; img2ft')
        endif
        call img%zero_and_flag_ft
        xcnt = 0
        ycnt = 0
        zcnt = 0
        lims = self%loop_lims(3)
        do h=lims(1,1),lims(1,2)
            xcnt = xcnt + 1
            if( xcnt > self%ldim(1) ) cycle
            ycnt = 0
            do k=lims(2,1),lims(2,2)
                ycnt = ycnt + 1
                if( ycnt > self%ldim(2) ) cycle
                zcnt = 0
                do l=lims(3,1),lims(3,2)
                    zcnt = zcnt + 1
                    if( zcnt > self%ldim(3) ) cycle
                    logi = [h,k,l]
                    phys = self%comp_addr_phys(logi)
                    call img%set_fcomp(logi, phys, cmplx(self%rmat(xcnt,ycnt,zcnt),0.))
                end do
            end do
        end do
    end subroutine img2ft

    module subroutine ft2img( self, which, img )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: which
        class(image),     intent(inout) :: img
        integer :: h,mh,k,mk,l,ml,lims(3,2),inds(3),phys(3)
        integer :: which_flag
        logical :: didft
        complex :: comp
        if( .not.(self.eqdims.img) )then
            write(logfhandle,*) 'self%ldim: ', self%ldim
            write(logfhandle,*) 'img%ldim:  ', img%ldim
            THROW_HARD('non-equal dims; ft2img')
        endif
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        select case(which)
            case ('real')
                which_flag = 0
            case('power')
                which_flag = 1
            case('sqrt')
                which_flag = 2
            case ('log')
                which_flag = 3
            case('phase')
                which_flag = 4
            case DEFAULT
                THROW_HARD('unsupported mode: '//trim(which)//'; ft2img')
        end select
        call img%zero_and_unflag_ft
        lims = self%loop_lims(3)
        mh   = abs(lims(1,1))
        mk   = abs(lims(2,1))
        ml   = abs(lims(3,1))
        if( .not.self%wthreads .and. self%is_2d() )then
            do k=lims(2,1),lims(2,2)
                inds(2) = min(max(1,k+mk+1),self%ldim(2))
                do h=lims(1,1),lims(1,2)
                    inds(1) = min(max(1,h+mh+1),self%ldim(1))
                    comp    = self%get_fcomp2D(h,k)
                    select case(which_flag)
                        case(0)
                            img%rmat(inds(1),inds(2),1) = real(comp)
                        case(1)
                            img%rmat(inds(1),inds(2),1) = csq(comp)
                        case(2)
                            img%rmat(inds(1),inds(2),1) = sqrt(csq(comp))
                        case(3)
                            img%rmat(inds(1),inds(2),1) = log(csq(comp))
                        case(4)
                            img%rmat(inds(1),inds(2),1) = phase_angle(comp)
                    end select
                end do
            end do
        else
            !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,comp,inds)&
            !$omp schedule(static) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        phys = self%comp_addr_phys([h,k,l])
                        comp = self%get_fcomp([h,k,l],phys)
                        inds(1) = min(max(1,h+mh+1),self%ldim(1))
                        inds(2) = min(max(1,k+mk+1),self%ldim(2))
                        inds(3) = min(max(1,l+ml+1),self%ldim(3))
                        select case(which_flag)
                        case (0)
                            call img%set(inds,real(comp))
                        case(1)
                            call img%set(inds,csq(comp))
                        case(2)
                            call img%set(inds,sqrt(csq(comp)))
                        case (3)
                            call img%set(inds,log(csq(comp)))
                        case(4)
                            call img%set(inds,phase_angle(comp))
                        end select
                    end do
                end do
            end do
            !$omp end parallel do
        endif
        if( didft ) call self%ifft()
    end subroutine ft2img

    !===========================
    ! Padding / normalization
    !===========================

    subroutine pad_fft( self, self_out )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: self_out
        integer       :: n1, n2, n1o, n2o, h1o, h2o
        integer       :: i, j
        integer       :: starts(3), stops(3)
        integer       :: x0, y0, xo, yo
        real(c_float) :: scale_cmat
        ! n3 is always 1 here
        n1  = self%ldim(1)
        n2  = self%ldim(2)
        n1o = self_out%ldim(1)
        n2o = self_out%ldim(2)
        h1o = n1o/2
        h2o = n2o/2
        starts = (self_out%ldim - self%ldim) / 2 + 1
        stops  = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.0_c_float
        ! ============================================================
        ! PAD + FFTSHIFT (fused): write each source pixel directly to
        ! its fftshifted output index.
        ! ============================================================
        do j = 1, n2
            y0 = starts(2) + j - 1
            yo = modulo((y0 - 1) + h2o, n2o) + 1
            do i = 1, n1
                x0 = starts(1) + i - 1
                xo = modulo((x0 - 1) + h1o, n1o) + 1
                self_out%rmat(xo, yo, 1) = self%rmat(i,j,1)
            end do
        end do
        ! ============================================================
        ! FFT (FFTW r2c) + scale with reciprocal (OUTPUT size)
        ! ============================================================
        call fftwf_execute_dft_r2c(self_out%plan_fwd, self_out%rmat, self_out%cmat)
        scale_cmat    = 1.0_c_float / real(n1o*n2o, c_float)
        self_out%cmat = self_out%cmat * scale_cmat
        self_out%ft   = .true.
    end subroutine pad_fft

    module subroutine norm_noise_fft( self, lmsk )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        integer       :: n1, n2, h1, h2, i, j, ii, jj, npix
        real(dp)      :: sum_dp, sum_sq_dp, mean_dp, var_dp, rnpix, xdp
        real(c_float) :: mean_sp, invstd_sp, scale_cmat, rswap
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        h1 = n1/2
        h2 = n2/2
        ! ---- two-pass background mean/variance on unmasked pixels ----
        sum_dp    = 0.0_dp
        sum_sq_dp = 0.0_dp
        npix      = 0
        do j = 1, n2
            do i = 1, n1
                if (.not. lmsk(i,j,1)) then
                    xdp = real(self%rmat(i,j,1), dp)
                    sum_dp    = sum_dp    + xdp
                    sum_sq_dp = sum_sq_dp + xdp*xdp
                    npix      = npix + 1
                end if
            end do
        end do
        mean_dp = 0.0_dp
        var_dp  = 0.0_dp
        if (npix > 1) then
            rnpix   = real(npix, dp)
            mean_dp = sum_dp / rnpix
            var_dp  = (sum_sq_dp - sum_dp*sum_dp/rnpix) / real(npix-1,dp)
        end if
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = 1.0_c_float
        if (var_dp > 0.0_dp) then
            invstd_sp = real(1.0_dp/sqrt(var_dp), c_float)
        end if
        ! ---- fftshift + normalization fused ----
        if (abs(real(mean_dp, kind=kind(1.0))) > TINY .or. invstd_sp /= 1.0_c_float) then
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    rswap = (self%rmat(i, j, 1) - mean_sp) * invstd_sp
                    self%rmat(i, j, 1) = (self%rmat(ii, jj, 1) - mean_sp) * invstd_sp
                    self%rmat(ii, jj, 1) = rswap

                    rswap = (self%rmat(i, jj, 1) - mean_sp) * invstd_sp
                    self%rmat(i, jj, 1) = (self%rmat(ii, j, 1) - mean_sp) * invstd_sp
                    self%rmat(ii, j, 1) = rswap
                end do
            end do
        else
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    rswap = self%rmat(i, j, 1)
                    self%rmat(i, j, 1)   = self%rmat(ii, jj, 1)
                    self%rmat(ii, jj, 1) = rswap

                    rswap = self%rmat(i, jj, 1)
                    self%rmat(i, jj, 1)  = self%rmat(ii, j, 1)
                    self%rmat(ii, j, 1)  = rswap
                end do
            end do
        end if
        ! ---- FFT + scale ----
        call fftwf_execute_dft_r2c(self%plan_fwd, self%rmat, self%cmat)
        scale_cmat = 1.0_c_float / real(n1*n2, c_float)
        self%cmat  = self%cmat * scale_cmat
        self%ft    = .true.
    end subroutine norm_noise_fft

    subroutine norm_noise_pad_fft( self, lmsk, self_out )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        integer       :: n1, n2, n1o, n2o, h1o, h2o, i, j, npix, starts(3), stops(3), x0, y0, xo, yo
        real(dp)      :: mean_dp, var_dp, sum_dp, sum_sq_dp, rnpix
        real(c_float) :: mean_sp, invstd_sp, scale_cmat, x_sp
        logical       :: do_norm
        ! n3 is always 1 here
        n1  = self%ldim(1)
        n2  = self%ldim(2)
        n1o = self_out%ldim(1)
        n2o = self_out%ldim(2)
        h1o = n1o/2
        h2o = n2o/2
        ! ============================================================
        ! NOISE NORMALIZATION: Optimized two-pass with simple sums
        ! ============================================================
        sum_dp    = 0.0_dp
        sum_sq_dp = 0.0_dp
        npix      = 0
        ! First pass: accumulate sums (column-major for Fortran)
        do i = 1, n1
            do j = 1, n2
                if (.not. lmsk(i,j,1)) then
                    npix = npix + 1
                    x_sp = self%rmat(i,j,1)
                    sum_dp = sum_dp + real(x_sp, dp)
                    sum_sq_dp = sum_sq_dp + real(x_sp, dp) * real(x_sp, dp)
                end if
            end do
        end do
        ! Compute mean and variance
        var_dp = 0.0_dp
        mean_dp = 0.0_dp
        if (npix > 1) then
            rnpix = real(npix, dp)
            mean_dp = sum_dp / rnpix
            var_dp = (sum_sq_dp - sum_dp * sum_dp / rnpix) / real(npix - 1, dp)
        end if
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = 1.0_c_float
        do_norm = .false.
        if (is_a_number(real(var_dp, kind=kind(1.0))) .and. var_dp > 0.0_dp) then
            invstd_sp = real(1.0_dp / sqrt(var_dp), c_float)
            do_norm = (abs(real(mean_dp, kind=kind(1.0))) > TINY .or. invstd_sp /= 1.0_c_float)
        end if
        ! ============================================================
        ! PAD + FFTSHIFT (fused) + apply normalization while copying
        !
        ! Instead of:
        !   (1) center copy into self_out
        !   (2) fftshift swap whole self_out
        !
        ! We write each source pixel directly to its fftshifted output index:
        !   x0 = starts(1)+i-1;  xo = ((x0-1 + h1o) mod n1o) + 1
        !   y0 = starts(2)+j-1;  yo = ((y0-1 + h2o) mod n2o) + 1
        ! ============================================================
        starts = (self_out%ldim - self%ldim) / 2 + 1
        stops  = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.0_c_float
        if (do_norm) then
            do j = 1, n2
                y0 = starts(2) + j - 1
                yo = modulo((y0 - 1) + h2o, n2o) + 1
                do i = 1, n1
                    x0 = starts(1) + i - 1
                    xo = modulo((x0 - 1) + h1o, n1o) + 1
                    self_out%rmat(xo, yo, 1) = (self%rmat(i,j,1) - mean_sp) * invstd_sp
                end do
            end do
        else
            do j = 1, n2
                y0 = starts(2) + j - 1
                yo = modulo((y0 - 1) + h2o, n2o) + 1
                do i = 1, n1
                    x0 = starts(1) + i - 1
                    xo = modulo((x0 - 1) + h1o, n1o) + 1
                    self_out%rmat(xo, yo, 1) = self%rmat(i,j,1)
                end do
            end do
        end if
        ! ============================================================
        ! FFT (FFTW r2c) + scale with reciprocal (OUTPUT size)
        ! ============================================================
        call fftwf_execute_dft_r2c(self_out%plan_fwd, self_out%rmat, self_out%cmat)
        scale_cmat    = 1.0_c_float / real(n1o*n2o, c_float)
        self_out%cmat = self_out%cmat * scale_cmat
        self_out%ft   = .true.
    end subroutine norm_noise_pad_fft

    subroutine norm_noise_pad_fft_clip_shift( self, lmsk, self_out, self_out2, shvec )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out, self_out2
        real,         intent(in)    :: shvec(2)
        real,         parameter     :: SHTHRESH = 0.001
        complex(c_float_complex) :: w1, w2, ph0, ph_h, ph_k, phase
        real(dp)      :: sh(2)
        real(dp)      :: mean_dp, var_dp, sum_dp, sum_sq_dp, rnpix
        real(c_float) :: mean_sp, invstd_sp, scale_cmat, x_sp, ratio
        integer       :: lims(3,2), starts(3), stops(3)
        integer       :: n1,n2, n1o,n2o, h1o,h2o, i,j, npix, x0,y0, xo,yo, h,k, kpi,kpo, hp
        logical       :: do_norm, k_neg
        ! n3 is always 1 here
        n1  = self%ldim(1)
        n2  = self%ldim(2)
        n1o = self_out%ldim(1)
        n2o = self_out%ldim(2)
        h1o = n1o/2
        h2o = n2o/2
        ! ============================================================
        ! NOISE NORMALIZATION: Optimized two-pass with simple sums
        ! ============================================================
        sum_dp    = 0.0_dp
        sum_sq_dp = 0.0_dp
        npix      = 0
        ! First pass: accumulate sums (column-major for Fortran)
        do i = 1, n1
            do j = 1, n2
                if (.not. lmsk(i,j,1)) then
                    npix = npix + 1
                    x_sp = self%rmat(i,j,1)
                    sum_dp = sum_dp + real(x_sp, dp)
                    sum_sq_dp = sum_sq_dp + real(x_sp, dp) * real(x_sp, dp)
                end if
            end do
        end do
        ! Compute mean and variance
        var_dp = 0.0_dp
        mean_dp = 0.0_dp
        if (npix > 1) then
            rnpix = real(npix, dp)
            mean_dp = sum_dp / rnpix
            var_dp = (sum_sq_dp - sum_dp * sum_dp / rnpix) / real(npix - 1, dp)
        end if
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = 1.0_c_float
        do_norm = .false.
        if (is_a_number(real(var_dp, kind=kind(1.0))) .and. var_dp > 0.0_dp) then
            invstd_sp = real(1.0_dp / sqrt(var_dp), c_float)
            do_norm = (abs(real(mean_dp, kind=kind(1.0))) > TINY .or. invstd_sp /= 1.0_c_float)
        end if
        ! ============================================================
        ! PAD + FFTSHIFT (fused) + apply normalization while copying
        !
        ! Instead of:
        !   (1) center copy into self_out
        !   (2) fftshift swap whole self_out
        !
        ! We write each source pixel directly to its fftshifted output index:
        !   x0 = starts(1)+i-1;  xo = ((x0-1 + h1o) mod n1o) + 1
        !   y0 = starts(2)+j-1;  yo = ((y0-1 + h2o) mod n2o) + 1
        ! ============================================================
        starts = (self_out%ldim - self%ldim) / 2 + 1
        stops  = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.0_c_float
        if (do_norm) then
            do j = 1, n2
                y0 = starts(2) + j - 1
                yo = modulo((y0 - 1) + h2o, n2o) + 1
                do i = 1, n1
                    x0 = starts(1) + i - 1
                    xo = modulo((x0 - 1) + h1o, n1o) + 1
                    self_out%rmat(xo, yo, 1) = (self%rmat(i,j,1) - mean_sp) * invstd_sp
                end do
            end do
        else
            do j = 1, n2
                y0 = starts(2) + j - 1
                yo = modulo((y0 - 1) + h2o, n2o) + 1
                do i = 1, n1
                    x0 = starts(1) + i - 1
                    xo = modulo((x0 - 1) + h1o, n1o) + 1
                    self_out%rmat(xo, yo, 1) = self%rmat(i,j,1)
                end do
            end do
        end if
        ! ============================================================
        ! FFT (FFTW r2c) + scale with reciprocal (OUTPUT size)
        ! ============================================================
        call fftwf_execute_dft_r2c(self_out%plan_fwd, self_out%rmat, self_out%cmat)
        scale_cmat  = 1.0_c_float / real(n1o*n2o, c_float) ! the scaling occurs below
        self_out%ft = .true.
        ! ============================================================
        ! CLIP & FFTW scaling + optional SHIFT in Fourier space
        ! ============================================================
        lims  = self_out2%fit%loop_lims(2)
        ratio = real(n1o, c_float) / real(self_out2%ldim(1), c_float)
        if (abs(shvec(1)) > SHTHRESH .or. abs(shvec(2)) > SHTHRESH) then
            sh = real(shvec * self_out2%shconst(1:2), dp)
            ! phase increments exp(i*sh1), exp(i*sh2)
            w1 = cmplx( real(cos(sh(1)), c_float), real(sin(sh(1)), c_float), kind=c_float_complex )
            w2 = cmplx( real(cos(sh(2)), c_float), real(sin(sh(2)), c_float), kind=c_float_complex )
            ! starting phase for h: exp(i*hmin*sh1)
            ph0 = cmplx( real(cos(real(lims(1,1),dp)*sh(1)), c_float), &
                        real(sin(real(lims(1,1),dp)*sh(1)), c_float), kind=c_float_complex )
            ! starting phase for k: exp(i*kmin*sh2)
            ph_k = cmplx( real(cos(real(lims(2,1),dp)*sh(2)), c_float), &
                        real(sin(real(lims(2,1),dp)*sh(2)), c_float), kind=c_float_complex )
            do k = lims(2,1), lims(2,2)
                k_neg =  k < 0
                kpi   = merge(k + 1 + n2o              , k + 1, k_neg)
                kpo   = merge(k + 1 + self_out2%ldim(2), k + 1, k_neg)
                ph_h  = ph0
                do h = lims(1,1), lims(1,2)
                    hp = h + 1
                    phase = ph_k * ph_h
                    self_out2%cmat(hp, kpo, 1) = scale_cmat * self_out%cmat(hp, kpi, 1) * phase
                    ph_h = ph_h * w1
                end do
                ph_k = ph_k * w2
            end do
        else
            ! CLIP only
            do k = lims(2,1), lims(2,2)
                k_neg =  k < 0
                kpi   = merge(k + 1 + n2o              , k + 1, k_neg)
                kpo   = merge(k + 1 + self_out2%ldim(2), k + 1, k_neg)
                do h = lims(1,1), lims(1,2)
                    hp = h + 1
                    self_out2%cmat(hp, kpo, 1) = scale_cmat * self_out%cmat(hp, kpi, 1)
                end do
            end do
        endif
        call self_out2%set_smpd(self%smpd * ratio)
        self_out2%ft = .true.
    end subroutine norm_noise_pad_fft_clip_shift

    ! The big players in this subroutine are the normalization (17%) and the fft (70%)  
    module subroutine norm_noise_fft_clip_shift( self, lmsk, self_out, shvec )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        real,         intent(in)    :: shvec(2)
        real,         parameter     :: SHTHRESH = 0.001
        ! FFT / phase
        complex(c_float_complex) :: w1, w2, ph0, ph_h, ph_k, phase
        real(dp)      :: sh(2)
        ! dimensions / indices
        integer       :: n1, n2, h1, h2, i, j, ii, jj, h, k, hp, lims(3,2), kpi, kpo, npix
        ! normalization
        real(dp)      :: sum_dp, sum_sq_dp, mean_dp, var_dp, rnpix, xdp
        real(c_float) :: mean_sp, invstd_sp, scale_cmat, ratio, rswap
        ! logical flag
        logical :: k_neg
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        h1 = n1/2
        h2 = n2/2
        ! ============================================================
        ! NORM_NOISE (two-pass)
        ! ============================================================
        sum_dp    = 0.0_dp
        sum_sq_dp = 0.0_dp
        npix      = 0
        do j = 1, n2
            do i = 1, n1
                if (.not. lmsk(i,j,1)) then
                    xdp = real(self%rmat(i,j,1), dp)
                    sum_dp    = sum_dp    + xdp
                    sum_sq_dp = sum_sq_dp + xdp*xdp
                    npix      = npix + 1
                end if
            end do
        end do
        mean_dp = 0.0_dp
        var_dp  = 0.0_dp
        if (npix > 1) then
            rnpix  = real(npix, dp)
            mean_dp = sum_dp / rnpix
            var_dp  = (sum_sq_dp - sum_dp*sum_dp/rnpix) / real(npix-1,dp)
        end if
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = 1.0_c_float
        if (var_dp > 0.0_dp) then
            invstd_sp = real(1.0_dp/sqrt(var_dp), c_float)
        end if
        ! ============================================================
        ! SHIFT TO PHASE ORIGIN (fftshift) + apply normalization fused
        ! ============================================================
        if (abs(real(mean_dp, kind=kind(1.0))) > TINY .or. invstd_sp /= 1.0_c_float) then
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    ! (1) swap (i,j) <-> (ii,jj), normalize on the fly
                    rswap = (self%rmat(i, j, 1) - mean_sp) * invstd_sp
                    self%rmat(i, j, 1) = (self%rmat(ii, jj, 1) - mean_sp) * invstd_sp
                    self%rmat(ii, jj, 1) = rswap
                    ! (2) swap (i,jj) <-> (ii,j), normalize on the fly
                    rswap = (self%rmat(i, jj, 1) - mean_sp) * invstd_sp
                    self%rmat(i, jj, 1) = (self%rmat(ii, j, 1) - mean_sp) * invstd_sp
                    self%rmat(ii, j, 1) = rswap
                end do
            end do
        else
            ! mean ~ 0 and invstd ~ 1: just swap
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    rswap = self%rmat(i, j, 1)
                    self%rmat(i, j, 1)   = self%rmat(ii, jj, 1)
                    self%rmat(ii, jj, 1) = rswap
                    rswap = self%rmat(i, jj, 1)
                    self%rmat(i, jj, 1)  = self%rmat(ii, j, 1)
                    self%rmat(ii, j, 1)  = rswap
                end do
            end do
        end if
        ! ============================================================
        ! FFT (FFTW r2c) + scale with reciprocal
        ! ============================================================
        call fftwf_execute_dft_r2c(self%plan_fwd, self%rmat, self%cmat)
        scale_cmat = 1.0_c_float / real(n1*n2, c_float)
        self%cmat  = self%cmat * scale_cmat
        self%ft    = .true.
        ! ============================================================
        ! CLIP (+ optional SHIFT in Fourier space) using phase recurrence
        ! ============================================================
        lims  = self_out%fit%loop_lims(2)
        ratio = real(n1, c_float) / real(self_out%ldim(1), c_float)
        if (abs(shvec(1)) > SHTHRESH .or. abs(shvec(2)) > SHTHRESH) then
            sh = real(shvec * self_out%shconst(1:2), dp)
            ! phase increments exp(i*sh1), exp(i*sh2)
            w1 = cmplx( real(cos(sh(1)), c_float), real(sin(sh(1)), c_float), kind=c_float_complex )
            w2 = cmplx( real(cos(sh(2)), c_float), real(sin(sh(2)), c_float), kind=c_float_complex )
            ! starting phase for h: exp(i*hmin*sh1)
            ph0 = cmplx( real(cos(real(lims(1,1),dp)*sh(1)), c_float), &
                        real(sin(real(lims(1,1),dp)*sh(1)), c_float), kind=c_float_complex )
            ! starting phase for k: exp(i*kmin*sh2)
            ph_k = cmplx( real(cos(real(lims(2,1),dp)*sh(2)), c_float), &
                        real(sin(real(lims(2,1),dp)*sh(2)), c_float), kind=c_float_complex )
            do k = lims(2,1), lims(2,2)
                k_neg =  k < 0
                kpi   = merge(k + 1 + n2              , k + 1, k_neg)
                kpo   = merge(k + 1 + self_out%ldim(2), k + 1, k_neg)
                ph_h  = ph0
                do h = lims(1,1), lims(1,2)
                    hp = h + 1
                    phase = ph_k * ph_h
                    self_out%cmat(hp, kpo, 1) = self%cmat(hp, kpi, 1) * phase
                    ph_h = ph_h * w1
                end do
                ph_k = ph_k * w2
            end do
        else
            ! CLIP only
            do k = lims(2,1), lims(2,2)
                k_neg =  k < 0
                kpi   = merge(k + 1 + n2              , k + 1, k_neg)
                kpo   = merge(k + 1 + self_out%ldim(2), k + 1, k_neg)
                do h = lims(1,1), lims(1,2)
                    hp = h + 1
                    self_out%cmat(hp, kpo, 1) = self%cmat(hp, kpi, 1)
                end do
            end do
        endif
        call self_out%set_smpd(self%smpd * ratio)
        self_out%ft = .true.
    end subroutine norm_noise_fft_clip_shift

    module subroutine norm_noise_divwinstrfun_fft( self, lmsk, instrfun )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(in)    :: instrfun
        integer       :: n1, n2, h1, h2, i, j, ii, jj, npix
        real(dp)      :: sum_dp, sum_sq_dp, mean_dp, var_dp, rnpix, xdp
        real(c_float) :: mean_sp, invstd_sp, scale_cmat
        real(c_float) :: v11, v22, v12, v21
        real(c_float) :: u11, u22, u12, u21
        real(c_float) :: rswap, invden
        real(c_float), parameter :: EPS_DEN = 1.e-6_c_float
        real(c_float), parameter :: ONE     = 1.0_c_float
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        h1 = n1/2
        h2 = n2/2
        ! ============================================================
        ! NORM_NOISE (two-pass) on unmasked pixels
        ! ============================================================
        sum_dp    = 0.0_dp
        sum_sq_dp = 0.0_dp
        npix      = 0
        do j = 1, n2
            do i = 1, n1
                if (.not. lmsk(i,j,1)) then
                    xdp = real(self%rmat(i,j,1), dp)
                    sum_dp    = sum_dp    + xdp
                    sum_sq_dp = sum_sq_dp + xdp*xdp
                    npix      = npix + 1
                end if
            end do
        end do
        mean_dp = 0.0_dp
        var_dp  = 0.0_dp
        if (npix > 1) then
            rnpix   = real(npix, dp)
            mean_dp = sum_dp / rnpix
            var_dp  = (sum_sq_dp - sum_dp*sum_dp/rnpix) / real(npix-1, dp)
        end if
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = ONE
        if (var_dp > 0.0_dp) then
            invstd_sp = real(ONE / sqrt(var_dp), c_float)
        end if
        ! ============================================================
        ! FFTSHIFT + (optional) normalization + division-by-instrfun
        !
        ! Division is applied in the *shifted* layout, consistent with
        ! ifft_mask_divwinstrfun_fft (which divides after fftshift).
        ! ============================================================
        if (abs(real(mean_dp, kind=kind(1.0))) > TINY .or. invstd_sp /= ONE) then
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    ! load 4 source pixels
                    v11 = self%rmat(i ,j ,1)
                    v22 = self%rmat(ii,jj,1)
                    v12 = self%rmat(i ,jj,1)
                    v21 = self%rmat(ii,j ,1)
                    ! normalize sources
                    v11 = (v11 - mean_sp) * invstd_sp
                    v22 = (v22 - mean_sp) * invstd_sp
                    v12 = (v12 - mean_sp) * invstd_sp
                    v21 = (v21 - mean_sp) * invstd_sp
                    ! denom values at DESTINATION (shifted) locations
                    u11 = instrfun%rmat(i ,j ,1)   ! dest for v22
                    u22 = instrfun%rmat(ii,jj,1)   ! dest for v11
                    u12 = instrfun%rmat(i ,jj,1)   ! dest for v21
                    u21 = instrfun%rmat(ii,j ,1)   ! dest for v12
                    ! apply safe division at destination
                    if (abs(u22) > EPS_DEN) then
                        invden = ONE / u22
                        v11 = v11 * invden
                    end if
                    if (abs(u11) > EPS_DEN) then
                        invden = ONE / u11
                        v22 = v22 * invden
                    end if
                    if (abs(u21) > EPS_DEN) then
                        invden = ONE / u21
                        v12 = v12 * invden
                    end if
                    if (abs(u12) > EPS_DEN) then
                        invden = ONE / u12
                        v21 = v21 * invden
                    end if
                    ! store swapped (fftshift)
                    self%rmat(i ,j ,1)  = v22
                    self%rmat(ii,jj,1)  = v11
                    self%rmat(i ,jj,1)  = v21
                    self%rmat(ii,j ,1)  = v12
                end do
            end do
        else
            ! mean ~ 0 and invstd ~ 1: only do division + shift
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    v11 = self%rmat(i ,j ,1)
                    v22 = self%rmat(ii,jj,1)
                    v12 = self%rmat(i ,jj,1)
                    v21 = self%rmat(ii,j ,1)
                    u11 = instrfun%rmat(i ,j ,1)
                    u22 = instrfun%rmat(ii,jj,1)
                    u12 = instrfun%rmat(i ,jj,1)
                    u21 = instrfun%rmat(ii,j ,1)
                    ! apply division at destination (same mapping)
                    if (abs(u22) > EPS_DEN) v11 = v11 * (ONE / u22)
                    if (abs(u11) > EPS_DEN) v22 = v22 * (ONE / u11)
                    if (abs(u21) > EPS_DEN) v12 = v12 * (ONE / u21)
                    if (abs(u12) > EPS_DEN) v21 = v21 * (ONE / u12)
                    self%rmat(i ,j ,1) = v22
                    self%rmat(ii,jj,1) = v11
                    self%rmat(i ,jj,1) = v21
                    self%rmat(ii,j ,1) = v12
                end do
            end do
        end if
        ! ============================================================
        ! FFT (FFTW r2c) + scale with reciprocal
        ! ============================================================
        call fftwf_execute_dft_r2c(self%plan_fwd, self%rmat, self%cmat)
        scale_cmat = ONE / real(n1*n2, c_float)
        self%cmat  = self%cmat * scale_cmat
        self%ft    = .true.
    end subroutine norm_noise_divwinstrfun_fft

    !>  \brief  Fused: noise-normalize (on unmasked bg) + fftshift + soft-avg mask + FFT + power spectrum
    !!
    !!  This is a *fuse-only* drop-in building block for later optimization:
    !!    1) norm_noise (two-pass stats on .not. lmsk, then fftshift + normalize)
    !!    2) mask2D_softavg (uses mem_msk_* tables + cosedge_r2_2d)
    !!    3) FFT r2c + scale
    !!    4) power_spectrum (2D path)
    !!
    !!  Assumptions:
    !!    - 2D images: self%ldim(3)=1
    !!    - self has plan_fwd, rmat/cmat allocated
    !!    - mem_msk_cs2/mem_msk_cs2 are already memoized for this image size
    module subroutine norm_noise_mask_fft_powspec( self, lmsk, mskrad, spec )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1), self%ldim(2), self%ldim(3))
        real,         intent(in)    :: mskrad
        real,         intent(inout) :: spec(fdim(self%ldim(1)) - 1)
        ! ---- locals: norm_noise_fft ----
        integer  :: n1, n2, h1, h2, i, j, ii, jj, npix
        real(dp) :: sum_dp, sum_sq_dp, mean_dp, var_dp, rnpix, xdp
        real(c_float) :: mean_sp, invstd_sp, scale_cmat, rswap
        ! ---- locals: mask2D_softavg ----
        real     :: rad_sq, ave, r2, e, cjs2
        integer  :: minlen, np, n1l, n2l
        real(dp) :: sv
        ! ---- locals: power_spectrum (2D only) ----
        integer  :: counts(fdim(self%ldim(1)) - 1)
        real(dp) :: dspec (fdim(self%ldim(1)) - 1)
        integer  :: lims(3,2), filtsz, h, k, sh, hp, kpi
        logical  :: k_neg
        ! ============================================================
        ! 1) NORM_NOISE (two-pass) + fftshift + normalize (fused)
        ! ============================================================
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        h1 = n1/2
        h2 = n2/2
        sum_dp    = 0.0_dp
        sum_sq_dp = 0.0_dp
        npix      = 0
        do j = 1, n2
            do i = 1, n1
                if (.not. lmsk(i,j,1)) then
                    xdp = real(self%rmat(i,j,1), dp)
                    sum_dp    = sum_dp    + xdp
                    sum_sq_dp = sum_sq_dp + xdp*xdp
                    npix      = npix + 1
                end if
            end do
        end do
        mean_dp = 0.0_dp
        var_dp  = 0.0_dp
        if (npix > 1) then
            rnpix   = real(npix, dp)
            mean_dp = sum_dp / rnpix
            var_dp  = (sum_sq_dp - sum_dp*sum_dp/rnpix) / real(npix-1, dp)
        end if
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = 1.0_c_float
        if (var_dp > 0.0_dp) then
            invstd_sp = real(1.0_dp/sqrt(var_dp), c_float)
        end if
        ! ---- fftshift + apply normalization on the fly ----
        if (abs(real(mean_dp, kind=kind(1.0))) > TINY .or. invstd_sp /= 1.0_c_float) then
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    rswap = (self%rmat(i,  j,  1) - mean_sp) * invstd_sp
                    self%rmat(i,  j,  1) = (self%rmat(ii, jj, 1) - mean_sp) * invstd_sp
                    self%rmat(ii, jj, 1) = rswap
                    rswap = (self%rmat(i,  jj, 1) - mean_sp) * invstd_sp
                    self%rmat(i,  jj, 1) = (self%rmat(ii, j,  1) - mean_sp) * invstd_sp
                    self%rmat(ii, j,  1) = rswap
                end do
            end do
        else
            do j = 1, h2
                jj = h2 + j
                do i = 1, h1
                    ii = h1 + i
                    rswap = self%rmat(i,  j,  1)
                    self%rmat(i,  j,  1)   = self%rmat(ii, jj, 1)
                    self%rmat(ii, jj, 1)   = rswap
                    rswap = self%rmat(i,  jj, 1)
                    self%rmat(i,  jj, 1)  = self%rmat(ii, j,  1)
                    self%rmat(ii, j,  1)  = rswap
                end do
            end do
        end if
        ! ============================================================
        ! 2) MASK (softavg) in real-space (same math as mask2D_softavg)
        ! ============================================================
        n1l = n1
        n2l = n2
        minlen = minval(self%ldim(1:2))
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        rad_sq = mskrad * mskrad
        sv = 0.0_dp
        np = 0
        do j = 1, n2l
            cjs2 =  mem_msk_cs2(j)
            do i = 1, n1l
                r2 = mem_msk_cs2(i) + cjs2
                if (r2 > rad_sq) then
                    np = np + 1
                    sv = sv + real(self%rmat(i,j,1), dp)
                endif
            end do
        end do
        if (np <= 0) return
        ave = real(sv / real(np, dp))
        do j = 1, n2l
            cjs2 = mem_msk_cs2(j)
            do i = 1, n1l
                r2 = mem_msk_cs2(i) + cjs2 
                e  = cosedge_r2_2d(r2, minlen, mskrad)
                if (e < 0.0001) then
                    self%rmat(i,j,1) = ave
                else if (e < 0.9999) then
                    self%rmat(i,j,1) = e*self%rmat(i,j,1) + (1.0-e)*ave
                endif
            end do
        end do
        ! ============================================================
        ! 3) FFT (r2c) + scale (same as norm_noise_fft)
        ! ============================================================
        call fftwf_execute_dft_r2c(self%plan_fwd, self%rmat, self%cmat)
        scale_cmat = 1.0_c_float / real(n1*n2, c_float)
        self%cmat  = self%cmat * scale_cmat
        self%ft    = .true.
        ! ============================================================
        ! 4) POWER SPECTRUM (2D path from power_spectrum)
        ! ============================================================
        filtsz = fdim(self%ldim(1)) - 1
        dspec  = 0.0_dp
        counts = 0
        lims   = self%fit%loop_lims(2)
        do k = lims(2,1), lims(2,2)
            k_neg = (k < 0)
            ! k physical index in [1..n2], wrapping negatives to the upper half
            kpi   = merge(k + 1 + n2, k + 1, k_neg)
            do h = lims(1,1), lims(1,2)
                hp = h + 1
                sh = nint(hyp(h,k))
                if (sh == 0 .or. sh > filtsz) cycle
                dspec(sh)  = dspec(sh)  + real(csq_fast(self%cmat(hp, kpi, 1)), dp)
                counts(sh) = counts(sh) + 1
            end do
        end do
        where(counts > 0)
            dspec = dspec / real(counts, dp)
        end where
        spec = real(dspec, kind=sp)
    end subroutine norm_noise_mask_fft_powspec

    module subroutine mask_divwinstrfun_fft( self, mskrad, instrfun )
        class(image), intent(inout) :: self
        real,         intent(in)    :: mskrad
        class(image), intent(in)    :: instrfun
        integer       :: n1, n2, h1, h2, i, j, ii, jj, minlen
        real(c_float) :: temp(self%ldim(1),self%ldim(2)), r2, scale_cmat, cis2, cjs2
        real(c_float) :: v_val, u_val, invden
        real(c_float), parameter :: EPS_DEN   = 1.e-6_c_float
        real(c_float), parameter :: EPS_E     = 1.e-4_c_float
        real(c_float), parameter :: ONE       = 1.0_c_float
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        h1 = n1/2
        h2 = n2/2
        ! ============================================================
        ! FUSED PASS:
        !   - apply soft mask
        !   - divide by instrument function where denom is safe
        !   - fftshift using temp buffer
        ! ============================================================
        ! mask & instrument function
        minlen = minval(self%ldim(1:2))
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        do j = 1, n2
            cjs2 = mem_msk_cs2(j)
            do i = 1, n1
                v_val = self%rmat(i, j, 1)
                u_val = instrfun%rmat(i, j, 1)
                ! soft mask
                r2    = mem_msk_cs2(i) + cjs2
                v_val = v_val * cosedge_r2_2d(r2, minlen, mskrad)
                ! Divide by instrument function if safe
                if (abs(u_val) > EPS_DEN) then
                    invden = ONE / u_val
                    v_val = v_val * invden
                end if
                temp(i, j) = v_val
            end do
        end do
        ! fftshift
        do j = 1, n2
            jj = mod(j + h2 - 1, n2) + 1
            do i = 1, n1
                ii = mod(i + h1 - 1, n1) + 1
                self%rmat(ii, jj, 1) = temp(i, j)
            end do
        end do
        ! ============================================================
        ! FFT (FFTW r2c) + scale with reciprocal
        ! ============================================================
        call fftwf_execute_dft_r2c(self%plan_fwd, self%rmat, self%cmat)
        scale_cmat = 1.0_c_float / real(n1*n2, c_float)
        self%cmat  = self%cmat * scale_cmat
        self%ft    = .true.
    end subroutine mask_divwinstrfun_fft

    module subroutine ifft_mask_divwinstrfun_fft( self, mskrad, instrfun )
        class(image), intent(inout) :: self
        real,         intent(in)    :: mskrad
        class(image), intent(in)    :: instrfun
        integer       :: n1, n2, h1, h2, i, j, ii, jj, minlen, np 
        real(dp)      :: sv
        real(c_float) :: temp(self%ldim(1),self%ldim(2)), rad_sq, ave, r2, e, scale_cmat, cis2, cjs2
        real(c_float) :: v_val, u_val, invden
        real(c_float), parameter :: EPS_DEN = 1.e-6_c_float
        real(c_float), parameter :: EPS_E   = 1.e-4_c_float
        real(c_float), parameter :: ONE     = 1.0_c_float
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        h1 = n1/2
        h2 = n2/2
        ! ============================================================
        ! IFFT
        ! ============================================================
        call fftwf_execute_dft_c2r(self%plan_bwd, self%cmat, self%rmat)
        self%ft = .false.
        ! ============================================================
        ! SHIFT TO PHASE ORIGIN (fftshift)
        ! ============================================================
        do j = 1, h2
            jj = h2 + j
            do i = 1, h1
                ii = h1 + i
                v_val = self%rmat(i ,j ,1);  self%rmat(i ,j ,1)  = self%rmat(ii,jj,1); self%rmat(ii,jj,1) = v_val
                v_val = self%rmat(i ,jj,1);  self%rmat(i ,jj,1)  = self%rmat(ii,j ,1); self%rmat(ii,j ,1) = v_val
            end do
        end do
        ! ============================================================
        ! MASK (softavg): compute average outside radius
        ! ============================================================
        minlen = minval(self%ldim(1:2))
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        rad_sq = real(mskrad*mskrad, c_float)
        sv = 0.0_dp
        np = 0
        do i = 1, n1
            cis2 = mem_msk_cs2(i)
            do j = 1, n2
                r2 = cis2 + mem_msk_cs2(j)
                if (r2 > rad_sq) then
                    np = np + 1
                    sv = sv + real(self%rmat(i,j,1), dp)
                end if
            end do
        end do
        if (np > 0) then
            ave = real(sv / real(np, dp), c_float)
        else
            ave = 0.0_c_float
        end if
        ! ============================================================
        ! FUSED PASS:
        !   - apply softavg mask
        !   - divide by instrument function where denom is safe
        !   - do second fftshift using temp buffer
        ! ============================================================
        ! Process with mask and division
        do j = 1, n2
            cjs2 = mem_msk_cs2(j)
            do i = 1, n1
                v_val = self%rmat(i, j, 1)
                u_val = instrfun%rmat(i, j, 1)
                ! Apply softavg mask
                r2 = mem_msk_cs2(i) + cjs2
                e  = cosedge_r2_2d(r2, minlen, mskrad)
                if (e < EPS_E) then
                    v_val = ave
                else if (e < (ONE - EPS_E)) then
                    v_val = e*v_val + (ONE - e)*ave
                end if
                ! Divide by instrument function if safe
                if (abs(u_val) > EPS_DEN) then
                    invden = ONE / u_val
                    v_val = v_val * invden
                end if
                temp(i, j) = v_val
            end do
        end do
        ! Now apply second fftshift from temp back to self%rmat
        do j = 1, n2
            jj = mod(j + h2 - 1, n2) + 1
            do i = 1, n1
                ii = mod(i + h1 - 1, n1) + 1
                self%rmat(ii, jj, 1) = temp(i, j)
            end do
        end do
        ! ============================================================
        ! FFT (FFTW r2c) + scale with reciprocal
        ! ============================================================
        call fftwf_execute_dft_r2c(self%plan_fwd, self%rmat, self%cmat)
        scale_cmat = 1.0_c_float / real(n1*n2, c_float)
        self%cmat  = self%cmat * scale_cmat
        self%ft    = .true.
    end subroutine ifft_mask_divwinstrfun_fft

    !>  \brief  expand_ft is for getting a Fourier plane using the old SIMPLE logics
    module function expand_ft( self ) result( fplane )
        class(image), intent(in) :: self
        complex, allocatable :: fplane(:,:)
        integer :: xdim, ydim, h, k, phys(3)
        if(is_even(self%ldim(1)))then
            xdim = self%ldim(1)/2
            ydim = self%ldim(2)/2
        else
            xdim = (self%ldim(1)-1)/2
            ydim = (self%ldim(2)-1)/2
        endif
        allocate(fplane(-xdim:xdim,-ydim:ydim))
        fplane = cmplx(0.,0.)
        do h=-xdim,xdim
            do k=-ydim,ydim
                phys = self%comp_addr_phys([h,k,0])
                fplane(h,k) = self%get_fcomp([h,k,0],phys)
            end do
        end do
    end function expand_ft

end submodule simple_image_fft
