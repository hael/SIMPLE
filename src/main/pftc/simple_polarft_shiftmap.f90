!@descr: polarft class submodule: the polar shift maps of a reference/particle pair over every rotation from the memoised spectra
! For every order n and magnitude rho, the all-rotation form of the Fourier-Bessel expansion (simple_polar_shiftmap):
!   a_n(irot, rho) = (1/nrots) sum_m exp(-2 pi i m (irot-1)/nrots) G_n[m, rho],
!   G_n[m, rho] = sum_k w_k J_n(z_k(rho)) R^_k[m-n] conj(Y^_k[m]),
! R^ and Y^ the full-ring spectra of the reference and of the CTF-modulated particle, from the memoised half spectra
! by F[-m] = (-1)^m conj(F[m]). On the grid's directions C(irot, alpha_j, rho) = Re(a_0)/2 + sum_{n>=1} Re(exp(-i n alpha_j) a_n)
! is the real part of one 2D transform of G over (n, m) per magnitude, the orders folded modulo nalpha (exact on the
! sampled directions); the real part of a transform is the transform of the Hermitian-symmetrised input, so it is one
! complex-to-real transform of half the size. Two maps: the full band's, and the screen's over a shorter band.
submodule (simple_polarft_calc) simple_polarft_shiftmap
#include "simple_local_flags.inc"
implicit none

contains

    !> Map imap: the tables of its shift grid (rho in pixels of the box the shift matrices act on, nalpha
    !! directions) over the shells kfromto(1)..kstop, the batched 2D plan and the per-thread work, one per thread of
    !! the memo workspace
    module subroutine new_shiftmap( self, imap, rho, nalpha, kstop )
        use simple_rnd, only: has_global_seed
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: imap
        real(dp),            intent(in)    :: rho(:)
        integer,             intent(in)    :: nalpha, kstop
        integer(c_int) :: n2(2), nin(2), plan_flags
        integer        :: ithr, nmax, nrho, nhalf, k1, k2
        if( .not. self%existence ) THROW_HARD('the calculator must exist; new_shiftmap')
        if( imap < 1 .or. imap > NSHIFTMAPS ) THROW_HARD('no such shift map; new_shiftmap')
        if( kstop <= self%kfromto(1) .or. kstop > self%kfromto(2) ) THROW_HARD('the map''s last shell must lie in the band; new_shiftmap')
        if( .not. allocated(self%shiftmap) ) allocate(self%shiftmap(NSHIFTMAPS), self%smwork(nthr_glob, NSHIFTMAPS))
        call self%kill_one_shiftmap(imap)
        k1 = self%kfromto(1)
        k2 = kstop
        call self%shiftmap(imap)%new(self%pftsz, [k1, k2], self%ldim(1), rho, nalpha)
        nmax  = self%shiftmap(imap)%get_nharm_max()
        nrho  = self%shiftmap(imap)%get_nrho()
        nhalf = nalpha / 2 + 1
        do ithr = 1, nthr_glob
            allocate(self%smwork(ithr,imap)%rhat(self%nrots, k1:k2),&
                    &self%smwork(ithr,imap)%yhat(self%nrots, k1:k2),&
                    &self%smwork(ithr,imap)%prod(self%nrots, k1:k2),&
                    &self%smwork(ithr,imap)%gsum(self%nrots, 0:nmax, nrho),&
                    &self%smwork(ithr,imap)%wj(k1:k2, nrho, 0:nmax),&
                    &self%smwork(ithr,imap)%refpow(self%nrots))
            self%smwork(ithr,imap)%gsum = cmplx(0., 0., kind=sp)
            allocate(self%smwork(ithr,imap)%gin(nalpha, self%nrots, nrho), source=cmplx(0., 0., kind=sp))
            allocate(self%smwork(ithr,imap)%lw(nalpha, self%nrots, nrho))
            self%smwork(ithr,imap)%ghalf_p = fftwf_alloc_complex(int(nhalf * self%nrots * nrho, c_size_t))
            self%smwork(ithr,imap)%gout_p  = fftwf_alloc_real(int(nalpha * self%nrots * nrho, c_size_t))
            call c_f_pointer(self%smwork(ithr,imap)%ghalf_p, self%smwork(ithr,imap)%ghalf, [nhalf, self%nrots, nrho])
            call c_f_pointer(self%smwork(ithr,imap)%gout_p,  self%smwork(ithr,imap)%gout,  [nalpha, self%nrots, nrho])
            self%smwork(ithr,imap)%iref = 0
        enddo
        ! one complex-to-real 2D transform per magnitude over (direction order, ring position), nrho of them in a
        ! batch, out of place; measured plans defeat exact replay, so a fixed seed takes the estimated one (as the
        ! memo workspace)
        plan_flags = FFTW_MEASURE
        if( has_global_seed() ) plan_flags = FFTW_ESTIMATE
        n2(1)   = int(self%nrots, c_int)
        n2(2)   = int(nalpha, c_int)
        nin(1)  = int(self%nrots, c_int)
        nin(2)  = int(nhalf, c_int)
        self%plan_shiftmap(imap) = fftwf_plan_many_dft_c2r(2_c_int, n2, int(nrho, c_int),&
            &self%smwork(1,imap)%ghalf, nin, 1_c_int, int(nhalf * self%nrots, c_int),&
            &self%smwork(1,imap)%gout,  n2,  1_c_int, int(nalpha * self%nrots, c_int), plan_flags)
        if( .not. c_associated(self%plan_shiftmap(imap)) ) THROW_HARD('FFTW failed to create the shift map plan; new_shiftmap')
    end subroutine new_shiftmap

    !> The two maps of the probability tables for a prior of width sigma_px (pixels of the shift matrices): the
    !! full band's and the screen's over the lower half of the band (kscreen, at least two shells), each with
    !! magnitudes a quarter of its band's wavelength apart out to the disk of three sigma and an FFT-friendly
    !! direction count resolving its ring at the outer magnitude; prior_w is the prior's weight 1 / (2 sigma^2)
    module subroutine setup_shiftmaps( self, sigma_px, radius_px, prior_w, kscreen )
        use simple_polar_shiftmap, only: shiftmap_ndirections
        class(polarft_calc), intent(inout) :: self
        real(dp),            intent(in)    :: sigma_px, radius_px
        real(dp),            intent(out)   :: prior_w
        integer,             intent(out)   :: kscreen
        real(dp) :: radius
        if( .not. (sigma_px  > 0.d0) ) THROW_HARD('the prior needs a positive width; setup_shiftmaps')
        if( .not. (radius_px > 0.d0) ) THROW_HARD('the disk needs a positive radius; setup_shiftmaps')
        prior_w = 1.d0 / (2.d0 * sigma_px**2)
        radius  = radius_px
        kscreen = max(self%kfromto(1) + 1, self%kfromto(2) / 2)
        call build_map(SHIFTMAP_FULL,   self%kfromto(2))
        call build_map(SHIFTMAP_SCREEN, kscreen)

    contains

        subroutine build_map( imap, kstop )
            integer, intent(in) :: imap, kstop
            real(dp), allocatable :: rhos(:)
            real(dp) :: step
            integer  :: nrho, i, nalpha
            step = real(self%ldim(1), dp) / real(kstop, dp) / 4.d0
            nrho = floor(radius / step) + 1
            if( nrho < 2 )then
                ! a disk below one step: the origin and its rim
                nrho = 2
                step = radius
            endif
            allocate(rhos(nrho))
            do i = 1, nrho
                rhos(i) = real(i-1, dp) * step
            enddo
            nalpha = shiftmap_ndirections(2 * ceiling(DTWOPI * real(kstop, dp) * rhos(nrho) / real(self%ldim(1), dp)) + 3)
            call self%new_shiftmap(imap, rhos, nalpha, kstop)
            deallocate(rhos)
        end subroutine build_map

    end subroutine setup_shiftmaps

    !> one map's tables, work and plan released (the arrays of maps stay)
    module subroutine kill_one_shiftmap( self, imap )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: imap
        integer :: ithr
        if( .not. allocated(self%shiftmap) ) return
        call self%shiftmap(imap)%kill
        do ithr = 1, size(self%smwork, 1)
            if( allocated(self%smwork(ithr,imap)%rhat) )then
                deallocate(self%smwork(ithr,imap)%rhat, self%smwork(ithr,imap)%yhat, self%smwork(ithr,imap)%prod,&
                    &self%smwork(ithr,imap)%gsum, self%smwork(ithr,imap)%wj, self%smwork(ithr,imap)%refpow,&
                    &self%smwork(ithr,imap)%gin, self%smwork(ithr,imap)%lw)
            endif
            if( c_associated(self%smwork(ithr,imap)%ghalf_p) )then
                call fftwf_free(self%smwork(ithr,imap)%ghalf_p)
                call fftwf_free(self%smwork(ithr,imap)%gout_p)
                self%smwork(ithr,imap)%ghalf_p = c_null_ptr
                self%smwork(ithr,imap)%gout_p  = c_null_ptr
                nullify(self%smwork(ithr,imap)%ghalf, self%smwork(ithr,imap)%gout)
            endif
            self%smwork(ithr,imap)%iref = 0
        enddo
        if( c_associated(self%plan_shiftmap(imap)) )then
            call fftwf_destroy_plan(self%plan_shiftmap(imap))
            self%plan_shiftmap(imap) = c_null_ptr
        endif
    end subroutine kill_one_shiftmap

    module subroutine kill_shiftmap( self )
        class(polarft_calc), intent(inout) :: self
        integer :: imap
        if( .not. allocated(self%shiftmap) ) return
        do imap = 1, NSHIFTMAPS
            call self%kill_one_shiftmap(imap)
        enddo
        deallocate(self%shiftmap, self%smwork)
    end subroutine kill_shiftmap

    module logical function has_shiftmap( self, imap )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: imap
        has_shiftmap = .false.
        if( allocated(self%shiftmap) ) has_shiftmap = self%shiftmap(imap)%exists()
    end function has_shiftmap

    module subroutine shiftmap_dims( self, imap, nrho, nalpha )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: imap
        integer,             intent(out) :: nrho, nalpha
        if( .not. self%has_shiftmap(imap) ) THROW_HARD('no shift map; new_shiftmap first; shiftmap_dims')
        nrho   = self%shiftmap(imap)%get_nrho()
        nalpha = self%shiftmap(imap)%get_nalpha()
    end subroutine shiftmap_dims

    !> the shift of a grid point in the optimizer's frame: the reference's own, the one the shift matrices act in
    module function shiftmap_shift( self, imap, irho, ialpha ) result( s )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: imap, irho, ialpha
        real(dp) :: s(2)
        s = self%shiftmap(imap)%get_shift(irho, ialpha)
    end function shiftmap_shift

    !> the translation area of one grid cell at magnitude irho (pixels^2), the sampler's measure on the polar grid
    module real(dp) function shiftmap_cell_area( self, imap, irho )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: imap, irho
        shiftmap_cell_area = self%shiftmap(imap)%get_cell_area(irho)
    end function shiftmap_cell_area

    !> the particle's weighted power over the shells kfromto(1)..kstop under the current sigma2 (the screen's band)
    module function wsqsum_ptcl_band( self, iptcl, kstop ) result( wsq )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: iptcl, kstop
        real(dp) :: wsq
        integer  :: i, ih, k
        i   = self%pinds(iptcl)
        ih  = merge(1, 2, self%iseven(i))
        wsq = 0.d0
        if( .not. associated(self%sigma2_noise) ) THROW_HARD('no sigma2; wsqsum_ptcl_band')
        do k = self%kfromto(1), min(kstop, self%kfromto(2))
            wsq = wsq + self%kshell_sqsums_ptcls(k,i) / real(self%sigma2_noise(k,ih), dp)
        enddo
    end function wsqsum_ptcl_band

    !> The pair's cross term on map imap's grid, into the calling thread's work: the shell sums G_n[m, rho] over
    !! the map's band for every order and magnitude, then one batched 2D transform; the weighted reference power
    !! under the rotated CTF per rotation alongside, and the particle's power over the band. shiftmap_objfun_grid,
    !! shiftmap_cross and shiftmap_draw read the work on the same thread. The shell weights are the objective's: k
    !! for cc, k/sigma2 for euclid
    module subroutine gen_shiftmap( self, imap, iref, iptcl )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: imap, iref, iptcl
        type(shiftmap_work), pointer :: wk
        real(dp) :: w(self%kfromto(1):self%kfromto(2)), wj(self%kfromto(1):self%kfromto(2))
        real(sp) :: fac
        integer  :: ithr, i, ieo, ih, k, m, n, irho, nrho, nmax, nalpha, nh, jn, sgn, nhalf, mm, nn, kb(2), k1, k2
        if( .not. self%has_shiftmap(imap) ) THROW_HARD('no shift map; new_shiftmap first; gen_shiftmap')
        ithr   =  omp_get_thread_num() + 1
        i      =  self%pinds(iptcl)
        ieo    =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        ih     =  merge(1, 2, self%iseven(i))
        wk     => self%smwork(ithr,imap)
        nrho   =  self%shiftmap(imap)%get_nrho()
        nmax   =  self%shiftmap(imap)%get_nharm_max()
        nalpha =  self%shiftmap(imap)%get_nalpha()
        kb     =  self%shiftmap(imap)%get_kfromto()
        k1     =  kb(1)
        k2     =  kb(2)
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                do k = k1, k2
                    w(k) = real(k, dp)
                enddo
                wk%ptclpow = 0.d0
                do k = k1, k2
                    wk%ptclpow = wk%ptclpow + self%kshell_sqsums_ptcls(k,i)
                enddo
            case(OBJFUN_EUCLID)
                if( .not. allocated(self%sigma2_wk8) ) THROW_HARD('euclid needs the sigma2 weights; gen_shiftmap')
                w(k1:k2)   = self%sigma2_wk8(k1:k2,ih)
                wk%ptclpow = self%wsqsum_ptcl_band(iptcl, k2)
            case DEFAULT
                THROW_HARD('unsupported objective; gen_shiftmap')
        end select
        wk%iref  = iref
        wk%iptcl = iptcl
        ! the Bessel-weighted shell weights of every order and magnitude
        do n = 0, nmax
            do irho = 1, nrho
                if( n > self%shiftmap(imap)%get_nharm(irho) ) cycle
                call self%shiftmap(imap)%shell_weights(w(k1:k2), n, irho, wj(k1:k2))
                wk%wj(:,irho,n) = real(wj(k1:k2), sp)
            enddo
        enddo
        ! the weighted reference power under the rotated CTF: sum_k w_k FT(CTF2) FT(REF2) back-transformed is
        ! 2 nrots times the half-ring sum at every rotation (the reference term of the Euclidean kernels)
        self%crvec1(ithr)%c = cmplx(0., 0., kind=c_float_complex)
        do k = k1, k2
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w(k), sp) * self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        wk%refpow = real(self%crvec1(ithr)%r(1:self%nrots), dp) / real(2 * self%nrots, dp)
        ! the full-ring spectra from the memoised half spectra: F[nrots - m] = (-1)^m conj(F[m])
        do k = k1, k2
            wk%rhat(1:self%pftsz+1,k) =       self%ft_ref(:,k,iref,ieo)
            wk%yhat(1:self%pftsz+1,k) = conjg(self%ft_ptcl_ctf(:,k,i))
            do m = 1, self%pftsz - 1
                sgn = 1 - 2 * mod(m, 2)
                wk%rhat(self%nrots-m+1,k) = real(sgn, sp) * conjg(self%ft_ref(m+1,k,iref,ieo))
                wk%yhat(self%nrots-m+1,k) = real(sgn, sp) *       self%ft_ptcl_ctf(m+1,k,i)
            enddo
        enddo
        ! the shell sums of every order: R^[m-n] conj(Y^[m]) by two contiguous slices, then one weighted sum per
        ! magnitude (single precision, as the spectra)
        do n = 0, nmax
            if( n == 0 )then
                wk%prod = wk%rhat * wk%yhat
            else
                wk%prod(1:n,:)            = wk%rhat(self%nrots-n+1:self%nrots,:) * wk%yhat(1:n,:)
                wk%prod(n+1:self%nrots,:) = wk%rhat(1:self%nrots-n,:)            * wk%yhat(n+1:self%nrots,:)
            endif
            do irho = 1, nrho
                if( n > self%shiftmap(imap)%get_nharm(irho) ) cycle
                wk%gsum(:,n,irho) = wk%wj(k1,irho,n) * wk%prod(:,k1)
                do k = k1 + 1, k2
                    wk%gsum(:,n,irho) = wk%gsum(:,n,irho) + wk%wj(k,irho,n) * wk%prod(:,k)
                enddo
            enddo
        enddo
        ! the orders onto the direction axis, folded modulo nalpha (a_0 enters with a half); the input's other rows
        ! are zero from new_shiftmap and stay so
        do irho = 1, nrho
            nh = self%shiftmap(imap)%get_nharm(irho)
            do n = 0, min(nh, nalpha - 1)
                fac = merge(0.5_sp, 1.0_sp, n == 0)
                wk%gin(n+1,:,irho) = fac * wk%gsum(:,n,irho)
            enddo
            do n = nalpha, nh
                jn = mod(n, nalpha) + 1
                wk%gin(jn,:,irho) = wk%gin(jn,:,irho) + wk%gsum(:,n,irho)
            enddo
        enddo
        ! the real part of the forward transform is the forward transform of the Hermitian half
        ! X_h[n,m] = (X[n,m] + conj(X[-n,-m]))/2; the complex-to-real transform is the backward one, so it is fed
        ! conj(X_h): one batched transform of half the size
        nhalf = nalpha / 2 + 1
        do irho = 1, nrho
            do m = 0, self%nrots - 1
                mm = modulo(-m, self%nrots) + 1
                do n = 0, nhalf - 1
                    nn = modulo(-n, nalpha) + 1
                    wk%ghalf(n+1,m+1,irho) = 0.5_sp * (conjg(wk%gin(n+1,m+1,irho)) + wk%gin(nn,mm,irho))
                enddo
            enddo
        enddo
        call fftwf_execute_dft_c2r(self%plan_shiftmap(imap), wk%ghalf, wk%gout)
    end subroutine gen_shiftmap

    !> the cross term of the thread's pair on map imap at grid rotation irot and magnitude rho(irho), at any direction
    !! alpha (the sub-grid readout): the rotation's coefficients from the shell sums by a direct transform
    module real(dp) function shiftmap_cross( self, imap, irot, irho, alpha )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: imap, irot, irho
        real(dp),            intent(in) :: alpha
        complex(sp) :: a(0:self%shiftmap(imap)%get_nharm_max())
        complex(dp) :: e, em, acc
        real(dp)    :: arg
        integer     :: ithr, n, m, nh
        ithr = omp_get_thread_num() + 1
        if( self%smwork(ithr,imap)%iref < 1 ) THROW_HARD('gen_shiftmap first, on this thread; shiftmap_cross')
        nh  = self%shiftmap(imap)%get_nharm(irho)
        arg = -DTWOPI * real(irot-1, dp) / real(self%nrots, dp)
        e   = cmplx(cos(arg), sin(arg), kind=dp)
        a   = cmplx(0., 0., kind=sp)
        do n = 0, nh
            acc = cmplx(0.d0, 0.d0, kind=dp)
            em  = cmplx(1.d0, 0.d0, kind=dp)
            do m = 1, self%nrots
                acc = acc + cmplx(self%smwork(ithr,imap)%gsum(m,n,irho), kind=dp) * em
                em  = em * e
            enddo
            a(n) = cmplx(acc / real(self%nrots, dp), kind=sp)
        enddo
        shiftmap_cross = self%shiftmap(imap)%eval_rotcoeffs(a, irho, alpha)
    end function shiftmap_cross

    !> The objective of the thread's pair on map imap's grid, vals(ialpha, irot, irho): the correlation for cc, the
    !! normalised Euclidean loss L = 1 + (|CTF R|^2_w - 2 C) / wsqsum (lower is better) for euclid, over the map's
    !! band. The data term only: the sampler adds the untempered prior and the cell measure itself
    module subroutine shiftmap_objfun_grid( self, imap, vals )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: imap
        real(sp),                    intent(out)   :: vals(:,:,:)
        type(shiftmap_work), pointer :: wk
        real(sp) :: scal(self%nrots), offs(self%nrots)
        real(dp) :: invn
        integer  :: ithr, nrho, nalpha, irho, irot
        ithr =  omp_get_thread_num() + 1
        wk   => self%smwork(ithr,imap)
        if( wk%iref < 1 ) THROW_HARD('gen_shiftmap first, on this thread; shiftmap_objfun_grid')
        nrho   = self%shiftmap(imap)%get_nrho()
        nalpha = self%shiftmap(imap)%get_nalpha()
        if( size(vals,1) /= nalpha .or. size(vals,2) /= self%nrots .or. size(vals,3) /= nrho )then
            THROW_HARD('vals must be (nalpha, nrots, nrho); shiftmap_objfun_grid')
        endif
        ! vals = offs(irot) + scal(irot) * C, C = gout / nrots
        invn = 1.d0 / real(self%nrots, dp)
        if( self%p_ptr%cc_objfun == OBJFUN_CC )then
            offs = 0.
            scal = real(invn / sqrt(wk%refpow * wk%ptclpow), sp)
        else
            offs = real(1.d0 + wk%refpow / wk%ptclpow, sp)
            scal = real(-2.d0 * invn / wk%ptclpow, sp)
        endif
        do irho = 1, nrho
            do irot = 1, self%nrots
                vals(:,irot,irho) = offs(irot) + scal(irot) * wk%gout(:,irot,irho)
            enddo
        enddo
    end subroutine shiftmap_objfun_grid

    !> the grid optimum of a shiftmap_objfun_grid (the largest correlation, the smallest loss): its rotation, its shift
    !! in the optimizer's frame, its value and its grid indices
    module subroutine shiftmap_peak( self, imap, vals, irot, xy, val, irho, ialpha )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: imap
        real(sp),            intent(in)  :: vals(:,:,:)
        integer,             intent(out) :: irot
        real(dp),            intent(out) :: xy(2)
        real(sp),            intent(out) :: val
        integer,             intent(out) :: irho, ialpha
        integer :: loc(3)
        if( self%p_ptr%cc_objfun == OBJFUN_CC )then
            loc = maxloc(vals)
        else
            loc = minloc(vals)
        endif
        ialpha = loc(1)
        irot   = loc(2)
        irho   = loc(3)
        val    = vals(ialpha, irot, irho)
        xy     = self%shiftmap(imap)%get_shift(irho, ialpha)
    end subroutine shiftmap_peak

    !> The prior over map imap's grid for a particle whose stored shift is u_prev (particle frame, the pixels of the
    !! shift matrices): a Gaussian of weight prior_w = 1 / (2 sigma^2) on the absolute displacement
    !! |u_prev + R(irot)^T x| of every cell, times the cell's area, over every (direction, rotation, magnitude) cell.
    !! lse_prior is its log normaliser; logp0 the log mass, to that normaliser, of the origin at one rotation (the
    !! central disk: the pose of a class scored at the stored shift)
    module subroutine shiftmap_prior( self, imap, prior_w, u_prev, lse_prior, logp0 )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: imap
        real(dp),            intent(in)  :: prior_w, u_prev(2)
        real(dp),            intent(out) :: lse_prior, logp0
        real(dp) :: c1(self%nrots), c2(self%nrots)
        real(dp) :: ca(self%shiftmap(imap)%get_nalpha()), sa(self%shiftmap(imap)%get_nalpha())
        real(dp) :: usq, rho, loga, l, lmax, s
        integer  :: nrho, nalpha, irho, irot, j, irho0
        if( .not. self%has_shiftmap(imap) ) THROW_HARD('no shift map; new_shiftmap first; shiftmap_prior')
        nrho   = self%shiftmap(imap)%get_nrho()
        nalpha = self%shiftmap(imap)%get_nalpha()
        call prior_tables(self, imap, u_prev, c1, c2, ca, sa)
        usq  = sum(u_prev**2)
        lmax = -huge(lmax)
        do irho = 1, nrho
            rho  = self%shiftmap(imap)%get_rho(irho)
            loga = log(self%shiftmap(imap)%get_cell_area(irho))
            do irot = 1, self%nrots
                do j = 1, nalpha
                    l    = -prior_w * (usq + rho * rho + 2.d0 * rho * (ca(j) * c1(irot) + sa(j) * c2(irot))) + loga
                    lmax = max(lmax, l)
                enddo
            enddo
        enddo
        s = 0.d0
        do irho = 1, nrho
            rho  = self%shiftmap(imap)%get_rho(irho)
            loga = log(self%shiftmap(imap)%get_cell_area(irho))
            do irot = 1, self%nrots
                do j = 1, nalpha
                    l = -prior_w * (usq + rho * rho + 2.d0 * rho * (ca(j) * c1(irot) + sa(j) * c2(irot))) + loga - lmax
                    if( l > -40.d0 ) s = s + exp(l)
                enddo
            enddo
        enddo
        lse_prior = lmax + log(s)
        ! the origin: the smallest magnitude's ring, whose cells all name one point at a rotation
        irho0 = 1
        do irho = 2, nrho
            if( self%shiftmap(imap)%get_rho(irho) < self%shiftmap(imap)%get_rho(irho0) ) irho0 = irho
        enddo
        logp0 = log(real(nalpha,dp) * self%shiftmap(imap)%get_cell_area(irho0)) - prior_w * usq - lse_prior
    end subroutine shiftmap_prior

    !> the prior's centre at every rotation, R(irot) u_prev with R = [[c, s], [-s, c]] (the shift search's frame), and
    !! map imap's directions' cosines and sines
    subroutine prior_tables( self, imap, u_prev, c1, c2, ca, sa )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: imap
        real(dp),            intent(in)  :: u_prev(2)
        real(dp),            intent(out) :: c1(self%nrots), c2(self%nrots), ca(:), sa(:)
        real(dp) :: ang, c, s
        integer  :: irot, j
        do irot = 1, self%nrots
            ang      = real(self%angtab(irot), dp) * DPI / 180.d0
            c        = cos(ang)
            s        = sin(ang)
            c1(irot) =  c * u_prev(1) + s * u_prev(2)
            c2(irot) = -s * u_prev(1) + c * u_prev(2)
        enddo
        do j = 1, size(ca)
            ang   = self%shiftmap(imap)%get_alpha(j)
            ca(j) = cos(ang)
            sa(j) = sin(ang)
        enddo
    end subroutine prior_tables

    !> One pose of the thread's pair drawn from the tempered posterior over map imap's grid, with mass per cell
    !! exp(-temp L) exp(-prior_w |u_prev + R(irot)^T x|^2) area (euclid: L the normalised loss over the map's band,
    !! temp the particle's temperature over that band times the tempering factor), and the pair's evidence
    !! logz = log of the summed mass over the prior's normaliser (shiftmap_prior). l_greedy takes the cell of the
    !! largest mass instead. loss is the data term at the returned cell; the cell's shift in the optimizer's frame
    !! is shiftmap_shift(imap, irho, ialpha)
    !> rho_out with mass_out: the posterior mass of the cells beyond magnitude rho_out (the window's radius: what
    !! the disk still has to find), at temperature temp_mass when given (the physical posterior, beta 1: the draw's
    !! tempering spreads the mass over the disk and would say the disk is always needed), else at temp
    module subroutine shiftmap_draw( self, imap, temp, prior_w, u_prev, lse_prior, l_greedy, logz, irot, irho, ialpha, loss,&
        &rho_out, mass_out, temp_mass )
        use simple_rnd, only: ran3
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: imap
        real(dp),                    intent(in)    :: temp, prior_w, u_prev(2), lse_prior
        logical,                     intent(in)    :: l_greedy
        real(dp),                    intent(out)   :: logz
        integer,                     intent(out)   :: irot, irho, ialpha
        real(sp),                    intent(out)   :: loss
        real(dp), optional,          intent(in)    :: rho_out
        real(dp), optional,          intent(out)   :: mass_out
        real(dp), optional,          intent(in)    :: temp_mass
        type(shiftmap_work), pointer :: wk
        real(dp) :: c1(self%nrots), c2(self%nrots)
        real(dp) :: ca(self%shiftmap(imap)%get_nalpha()), sa(self%shiftmap(imap)%get_nalpha())
        real(dp) :: offs(self%nrots), scal(self%nrots), usq, rho, loga, lmax, s, tgt, l, pen, s_out, s_ring, tm, lm
        integer  :: ithr, nrho, nalpha, ir, ro, j, imax(3)
        ithr =  omp_get_thread_num() + 1
        wk   => self%smwork(ithr,imap)
        if( wk%iref < 1 ) THROW_HARD('gen_shiftmap first, on this thread; shiftmap_draw')
        if( self%p_ptr%cc_objfun /= OBJFUN_EUCLID ) THROW_HARD('the tempered posterior is the Gaussian noise model''s: objfun=euclid; shiftmap_draw')
        nrho   = self%shiftmap(imap)%get_nrho()
        nalpha = self%shiftmap(imap)%get_nalpha()
        call prior_tables(self, imap, u_prev, c1, c2, ca, sa)
        usq  = sum(u_prev**2)
        offs = 1.d0 + wk%refpow / wk%ptclpow
        scal = -2.d0 / (real(self%nrots, dp) * wk%ptclpow)
        ! the log masses, their maximum and its cell
        lmax = -huge(lmax)
        imax = 1
        do ir = 1, nrho
            rho  = self%shiftmap(imap)%get_rho(ir)
            loga = log(self%shiftmap(imap)%get_cell_area(ir))
            do ro = 1, self%nrots
                do j = 1, nalpha
                    pen = prior_w * (usq + rho * rho + 2.d0 * rho * (ca(j) * c1(ro) + sa(j) * c2(ro)))
                    l   = -temp * (offs(ro) + scal(ro) * real(wk%gout(j,ro,ir), dp)) - pen + loga
                    wk%lw(j,ro,ir) = real(l, sp)
                    if( l > lmax )then
                        lmax = l
                        imax = [j, ro, ir]
                    endif
                enddo
            enddo
        enddo
        s     = 0.d0
        s_out = 0.d0
        do ir = 1, nrho
            s_ring = 0.d0
            do ro = 1, self%nrots
                do j = 1, nalpha
                    l = real(wk%lw(j,ro,ir), dp) - lmax
                    if( l > -40.d0 ) s_ring = s_ring + exp(l)
                enddo
            enddo
            s = s + s_ring
            if( present(rho_out) )then
                if( self%shiftmap(imap)%get_rho(ir) > rho_out ) s_out = s_out + s_ring
            endif
        enddo
        logz = lmax + log(s) - lse_prior
        if( present(mass_out) )then
            mass_out = s_out / s
            if( present(temp_mass) )then
                ! the mass beyond at the physical temperature: the data term rescaled by temp_mass / temp
                tm = temp_mass / max(temp, 1.d-30)
                lm = -huge(lm)
                do ir = 1, nrho
                    rho  = self%shiftmap(imap)%get_rho(ir)
                    loga = log(self%shiftmap(imap)%get_cell_area(ir))
                    do ro = 1, self%nrots
                        do j = 1, nalpha
                            pen = prior_w * (usq + rho * rho + 2.d0 * rho * (ca(j) * c1(ro) + sa(j) * c2(ro)))
                            l   = tm * (real(wk%lw(j,ro,ir), dp) + pen - loga) - pen + loga
                            if( l > lm ) lm = l
                        enddo
                    enddo
                enddo
                s     = 0.d0
                s_out = 0.d0
                do ir = 1, nrho
                    rho  = self%shiftmap(imap)%get_rho(ir)
                    loga = log(self%shiftmap(imap)%get_cell_area(ir))
                    s_ring = 0.d0
                    do ro = 1, self%nrots
                        do j = 1, nalpha
                            pen = prior_w * (usq + rho * rho + 2.d0 * rho * (ca(j) * c1(ro) + sa(j) * c2(ro)))
                            l   = tm * (real(wk%lw(j,ro,ir), dp) + pen - loga) - pen + loga - lm
                            if( l > -40.d0 ) s_ring = s_ring + exp(l)
                        enddo
                    enddo
                    s = s + s_ring
                    if( rho > rho_out ) s_out = s_out + s_ring
                enddo
                mass_out = s_out / s
            endif
        endif
        if( l_greedy )then
            ialpha = imax(1)
            irot   = imax(2)
            irho   = imax(3)
        else
            tgt    = real(ran3(), dp) * s
            s      = 0.d0
            ialpha = imax(1)
            irot   = imax(2)
            irho   = imax(3)
            outer: do ir = 1, nrho
                do ro = 1, self%nrots
                    do j = 1, nalpha
                        l = real(wk%lw(j,ro,ir), dp) - lmax
                        if( l <= -40.d0 ) cycle
                        s = s + exp(l)
                        if( s >= tgt )then
                            ialpha = j
                            irot   = ro
                            irho   = ir
                            exit outer
                        endif
                    enddo
                enddo
            enddo outer
        endif
        loss = real(offs(irot) + scal(irot) * real(wk%gout(ialpha,irot,irho), dp), sp)
    end subroutine shiftmap_draw

end submodule simple_polarft_shiftmap
