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

    module subroutine pad_fft( self, self_out )
        class(image), intent(inout) :: self
        class(image), intent(inout) :: self_out
        integer :: starts(3), stops(3)
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine pad_fft

    subroutine norm_noise_pad_fft( self, lmsk, self_out )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        integer :: npix, starts(3), stops(3)
        real    :: ave, var, ep, sdev_noise
        npix = product(self%ldim) - count(lmsk) ! # background pixels
        ave  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.not. lmsk) / real(npix) ! background average
        if( abs(ave) > TINY ) self%rmat = self%rmat - ave
        ep         = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),      mask=.not. lmsk)
        var        = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))**2.0, mask=.not. lmsk)
        var        = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        sdev_noise = 0.
        if( is_a_number(var) )then
            sdev_noise = sqrt(var)
            if( var > 0. ) self%rmat = self%rmat / sdev_noise
        endif
        starts        = (self_out%ldim - self%ldim) / 2 + 1
        stops         = self_out%ldim - starts + 1
        self_out%ft   = .false.
        self_out%rmat = 0.
        self_out%rmat(starts(1):stops(1),starts(2):stops(2),1)=&
            &self%rmat(:self%ldim(1),:self%ldim(2),1)
        call self_out%fft
    end subroutine norm_noise_pad_fft

    module subroutine norm_noise_fft_clip_shift( self, lmsk, self_out, shvec )
        use, intrinsic :: iso_c_binding, only: c_float, c_float_complex
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        class(image), intent(inout) :: self_out
        real,         intent(in)    :: shvec(2)
        real,          parameter :: SHTHRESH = 0.001
        complex(c_float_complex) :: w1, w2, ph0, ph_h, ph_k, phase
        integer       :: n1, n2, h1, h2,  i, j, ii, jj, npix, lims(3,2), h, hp, k, kpi, kpo 
        real(dp)      :: mean_dp, M2_dp, x_dp, delta, var_dp,  sh(2)
        real(c_float) :: mean_sp, invstd_sp,  rswap, scale_cmat, ratio
        logical       :: k_neg
        ! n3 is always 1 here
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        h1 = n1/2
        h2 = n2/2
        ! ============================================================
        ! NORM_NOISE: one-pass Welford over background pixels (no temps)
        ! ============================================================
        mean_dp = 0.0_dp
        M2_dp   = 0.0_dp
        npix    = 0
        do j = 1, n2
            do i = 1, n1
                if (.not. lmsk(i,j,1)) then
                    npix = npix + 1
                    x_dp = real(self%rmat(i,j,1), dp)
                    delta   = x_dp - mean_dp
                    mean_dp = mean_dp + delta / real(npix, dp)
                    M2_dp   = M2_dp + delta * (x_dp - mean_dp)
                end if
            end do
        end do
        var_dp = 0.0_dp
        if (npix > 1) var_dp = M2_dp / real(npix-1, dp)
        mean_sp   = real(mean_dp, c_float)
        invstd_sp = 1.0_c_float
        if (is_a_number(real(var_dp, kind=kind(1.0))) .and. var_dp > 0.0_dp) then
            invstd_sp = real(1.0_dp / sqrt(var_dp), c_float)
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
                    self%rmat(i, j, 1) = self%rmat(ii, jj, 1)
                    self%rmat(ii, jj, 1) = rswap
                    rswap = self%rmat(i, jj, 1)
                    self%rmat(i, jj, 1) = self%rmat(ii, j, 1)
                    self%rmat(ii, j, 1) = rswap
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
            call self_out%set_smpd(self%smpd * ratio)
            self_out%ft = .true.
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
            call self_out%set_smpd(self%smpd * ratio)
            self_out%ft = .true.
        endif
    end subroutine norm_noise_fft_clip_shift

    module subroutine ifft_mask_divwinstrfun_fft( self, mskrad, instrfun )
        class(image), intent(inout) :: self
        real,         intent(in)    :: mskrad
        class(image), intent(in)    :: instrfun
        integer       :: n1, n2, n3, h1, h2,  i, j, ii, jj, lims(3,2), minlen, npix, i, j, np, n1l, n2l
        real(c_float), parameter :: WWIDTH = 10.
        real(c_float) :: rswap
        real(dp)      :: sumv, sv
        real(c_float) :: rad_sq, ave, r2, e, scale_cmat
        ! n3 is always 1 here
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        n3 = 1
        h1 = n1/2
        h2 = n2/2
        ! ============================================================
        ! IFFT
        ! ============================================================
        call fftwf_execute_dft_c2r(self%plan_bwd,self%cmat,self%rmat)
        self%ft = .false.
        ! ============================================================
        ! SHIFT TO PHASE ORIGIN (fftshift)
        ! ============================================================
        do j = 1, h2
            jj = h2 + j
            do i = 1, h1
                ii = h1 + i
                rswap = self%rmat(i, j, 1)
                self%rmat(i, j, 1) = self%rmat(ii, jj, 1)
                self%rmat(ii, jj, 1) = rswap
                rswap = self%rmat(i, jj, 1)
                self%rmat(i, jj, 1) = self%rmat(ii, j, 1)
                self%rmat(ii, j, 1) = rswap
            end do
        end do
        ! ============================================================
        ! MASK (softavg)
        ! ============================================================
        ! minlen
        minlen = minval(self%ldim(1:2))
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        ! mode
        rad_sq     = mskrad * mskrad
        sumv       = 0.0_dp
        npix       = 0            
        n1l = n1; n2l = n2
        sv = 0.0_dp; np = 0
        ! avg outside radius
        do j = 1, n2l
            do i = 1, n1l
                r2 = mem_msk_cis2(i) + mem_msk_cjs2(j)
                if (r2 > rad_sq) then
                    np = np + 1
                    sv = sv + real(self%rmat(i,j,1), dp)
                endif
            end do
        end do
        if (np > 0) then
            ave = real(sv / real(np, dp))
            ! apply (j,i with i contiguous)
            do j = 1, n2l
                do i = 1, n1l
                    r2 = mem_msk_cis2(i) + mem_msk_cjs2(j)
                    e  = cosedge_r2_2d(r2, minlen, mskrad)
                    if (e < 0.0001) then
                        self%rmat(i,j,1) = ave
                    else if (e < 0.9999) then
                        self%rmat(i,j,1) = e*self%rmat(i,j,1) + (1.0-e)*ave
                    endif
                end do
            end do
        endif
        ! ============================================================
        ! DIVIDE WITH INSTRUMENT FUNCTION (mul w inv of instr)
        ! ============================================================
        where(abs(self%rmat) > 1.e-6) self%rmat = self%rmat/instrfun%rmat
        ! ============================================================
        ! SHIFT TO PHASE ORIGIN (fftshift)
        ! ============================================================
        do j = 1, h2
            jj = h2 + j
            do i = 1, h1
                ii = h1 + i
                rswap = self%rmat(i, j, 1)
                self%rmat(i, j, 1) = self%rmat(ii, jj, 1)
                self%rmat(ii, jj, 1) = rswap
                rswap = self%rmat(i, jj, 1)
                self%rmat(i, jj, 1) = self%rmat(ii, j, 1)
                self%rmat(ii, j, 1) = rswap
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
