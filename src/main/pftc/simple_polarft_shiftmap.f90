!@descr: polarft class submodule: the polar shift map of a reference/particle pair over every rotation from the memoised spectra
! For every order n and magnitude rho, the all-rotation form of the Fourier-Bessel expansion (simple_polar_shiftmap):
!   a_n(irot, rho) = (1/nrots) sum_m exp(-2 pi i m (irot-1)/nrots) G_n[m, rho],
!   G_n[m, rho] = sum_k w_k J_n(z_k(rho)) R^_k[m-n] conj(Y^_k[m]),
! R^ and Y^ the full-ring spectra of the reference and of the CTF-modulated particle, from the memoised half spectra
! by F[-m] = (-1)^m conj(F[m]). On the grid's directions C(irot, alpha_j, rho) = Re(a_0)/2 + sum_{n>=1} Re(exp(-i n alpha_j) a_n)
! is the real part of one 2D transform of G over (n, m) per magnitude, the orders folded modulo nalpha (exact on the
! sampled directions); the real part of a transform is the transform of the Hermitian-symmetrised input, so it is one
! complex-to-real transform of half the size.
submodule (simple_polarft_calc) simple_polarft_shiftmap
#include "simple_local_flags.inc"
implicit none

contains

    !> the tables of the shift grid (rho in pixels of the box the shift matrices act on, nalpha directions), the
    !! batched 2D plan and the per-thread work, one per thread of the memo workspace
    module subroutine new_shiftmap( self, rho, nalpha )
        use simple_rnd, only: has_global_seed
        class(polarft_calc), intent(inout) :: self
        real(dp),            intent(in)    :: rho(:)
        integer,             intent(in)    :: nalpha
        integer(c_int) :: n2(2), nin(2), plan_flags
        integer        :: ithr, nmax, nrho, nhalf
        call self%kill_shiftmap
        if( .not. self%existence ) THROW_HARD('the calculator must exist; new_shiftmap')
        allocate(self%shiftmap)
        call self%shiftmap%new(self%pftsz, self%kfromto, self%ldim(1), rho, nalpha)
        nmax  = self%shiftmap%get_nharm_max()
        nrho  = self%shiftmap%get_nrho()
        nhalf = nalpha / 2 + 1
        allocate(self%smwork(nthr_glob))
        do ithr = 1, nthr_glob
            allocate(self%smwork(ithr)%rhat(self%nrots, self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%yhat(self%nrots, self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%prod(self%nrots, self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%gsum(self%nrots, 0:nmax, nrho),&
                    &self%smwork(ithr)%wj(self%kfromto(1):self%kfromto(2), nrho, 0:nmax),&
                    &self%smwork(ithr)%refpow(self%nrots))
            self%smwork(ithr)%gsum = cmplx(0., 0., kind=sp)
            allocate(self%smwork(ithr)%gin(nalpha, self%nrots, nrho), source=cmplx(0., 0., kind=sp))
            self%smwork(ithr)%ghalf_p = fftwf_alloc_complex(int(nhalf * self%nrots * nrho, c_size_t))
            self%smwork(ithr)%gout_p  = fftwf_alloc_real(int(nalpha * self%nrots * nrho, c_size_t))
            call c_f_pointer(self%smwork(ithr)%ghalf_p, self%smwork(ithr)%ghalf, [nhalf, self%nrots, nrho])
            call c_f_pointer(self%smwork(ithr)%gout_p,  self%smwork(ithr)%gout,  [nalpha, self%nrots, nrho])
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
        self%plan_shiftmap = fftwf_plan_many_dft_c2r(2_c_int, n2, int(nrho, c_int),&
            &self%smwork(1)%ghalf, nin, 1_c_int, int(nhalf * self%nrots, c_int),&
            &self%smwork(1)%gout,  n2,  1_c_int, int(nalpha * self%nrots, c_int), plan_flags)
        if( .not. c_associated(self%plan_shiftmap) ) THROW_HARD('FFTW failed to create the shift map plan; new_shiftmap')
    end subroutine new_shiftmap

    module subroutine kill_shiftmap( self )
        class(polarft_calc), intent(inout) :: self
        integer :: ithr
        if( allocated(self%shiftmap) )then
            call self%shiftmap%kill
            deallocate(self%shiftmap)
        endif
        if( allocated(self%smwork) )then
            do ithr = 1, size(self%smwork)
                if( allocated(self%smwork(ithr)%rhat) )then
                    deallocate(self%smwork(ithr)%rhat, self%smwork(ithr)%yhat, self%smwork(ithr)%prod,&
                        &self%smwork(ithr)%gsum, self%smwork(ithr)%wj, self%smwork(ithr)%refpow, self%smwork(ithr)%gin)
                endif
                if( c_associated(self%smwork(ithr)%ghalf_p) )then
                    call fftwf_free(self%smwork(ithr)%ghalf_p)
                    call fftwf_free(self%smwork(ithr)%gout_p)
                    self%smwork(ithr)%ghalf_p = c_null_ptr
                    self%smwork(ithr)%gout_p  = c_null_ptr
                    nullify(self%smwork(ithr)%ghalf, self%smwork(ithr)%gout)
                endif
            enddo
            deallocate(self%smwork)
            call fftwf_destroy_plan(self%plan_shiftmap)
            self%plan_shiftmap = c_null_ptr
        endif
    end subroutine kill_shiftmap

    module logical function has_shiftmap( self )
        class(polarft_calc), intent(in) :: self
        has_shiftmap = allocated(self%shiftmap)
    end function has_shiftmap

    module subroutine shiftmap_dims( self, nrho, nalpha )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(out) :: nrho, nalpha
        if( .not. allocated(self%shiftmap) ) THROW_HARD('no shift map; new_shiftmap first; shiftmap_dims')
        nrho   = self%shiftmap%get_nrho()
        nalpha = self%shiftmap%get_nalpha()
    end subroutine shiftmap_dims

    !> the shift of a grid point in the optimizer's frame: the reference's own, the one the shift matrices act in
    module function shiftmap_shift( self, irho, ialpha ) result( s )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: irho, ialpha
        real(dp) :: s(2)
        s = self%shiftmap%get_shift(irho, ialpha)
    end function shiftmap_shift

    !> the translation area of one grid cell at magnitude irho (pixels^2), the sampler's measure on the polar grid
    module real(dp) function shiftmap_cell_area( self, irho )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: irho
        shiftmap_cell_area = self%shiftmap%get_cell_area(irho)
    end function shiftmap_cell_area

    !> The pair's cross term on the whole grid, into the calling thread's work: the shell sums G_n[m, rho] for every
    !! order and magnitude, then one batched 2D transform; the weighted reference power under the rotated CTF per
    !! rotation alongside. shiftmap_objfun_grid and shiftmap_cross read the work on the same thread. The shell
    !! weights are the objective's: k for cc, k/sigma2 for euclid
    module subroutine gen_shiftmap( self, iref, iptcl )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        type(shiftmap_work), pointer :: wk
        real(dp) :: w(self%kfromto(1):self%kfromto(2)), wj(self%kfromto(1):self%kfromto(2))
        real(sp) :: fac
        integer  :: ithr, i, ieo, ih, k, m, n, irho, nrho, nmax, nalpha, nh, jn, sgn, nhalf, mm, nn
        if( .not. allocated(self%shiftmap) ) THROW_HARD('no shift map; new_shiftmap first; gen_shiftmap')
        ithr   =  omp_get_thread_num() + 1
        i      =  self%pinds(iptcl)
        ieo    =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        ih     =  merge(1, 2, self%iseven(i))
        wk     => self%smwork(ithr)
        nrho   =  self%shiftmap%get_nrho()
        nmax   =  self%shiftmap%get_nharm_max()
        nalpha =  self%shiftmap%get_nalpha()
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                do k = self%kfromto(1), self%kfromto(2)
                    w(k) = real(k, dp)
                enddo
                wk%ptclpow = self%ksqsums_ptcls(i)
            case(OBJFUN_EUCLID)
                if( .not. allocated(self%sigma2_wk8) ) THROW_HARD('euclid needs the sigma2 weights; gen_shiftmap')
                w          = self%sigma2_wk8(:,ih)
                wk%ptclpow = self%wsqsums_ptcls(i)
            case DEFAULT
                THROW_HARD('unsupported objective; gen_shiftmap')
        end select
        wk%iref  = iref
        wk%iptcl = iptcl
        ! the Bessel-weighted shell weights of every order and magnitude
        do n = 0, nmax
            do irho = 1, nrho
                if( n > self%shiftmap%get_nharm(irho) ) cycle
                call self%shiftmap%shell_weights(w, n, irho, wj)
                wk%wj(:,irho,n) = real(wj, sp)
            enddo
        enddo
        ! the weighted reference power under the rotated CTF: sum_k w_k FT(CTF2) FT(REF2) back-transformed is
        ! 2 nrots times the half-ring sum at every rotation (the reference term of the Euclidean kernels)
        self%crvec1(ithr)%c = cmplx(0., 0., kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w(k), sp) * self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        wk%refpow = real(self%crvec1(ithr)%r(1:self%nrots), dp) / real(2 * self%nrots, dp)
        ! the full-ring spectra from the memoised half spectra: F[nrots - m] = (-1)^m conj(F[m])
        do k = self%kfromto(1), self%kfromto(2)
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
                if( n > self%shiftmap%get_nharm(irho) ) cycle
                wk%gsum(:,n,irho) = wk%wj(self%kfromto(1),irho,n) * wk%prod(:,self%kfromto(1))
                do k = self%kfromto(1) + 1, self%kfromto(2)
                    wk%gsum(:,n,irho) = wk%gsum(:,n,irho) + wk%wj(k,irho,n) * wk%prod(:,k)
                enddo
            enddo
        enddo
        ! the orders onto the direction axis, folded modulo nalpha (a_0 enters with a half); the input's other rows
        ! are zero from new_shiftmap and stay so
        do irho = 1, nrho
            nh = self%shiftmap%get_nharm(irho)
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
        call fftwf_execute_dft_c2r(self%plan_shiftmap, wk%ghalf, wk%gout)
    end subroutine gen_shiftmap

    !> the cross term of the thread's pair at grid rotation irot and magnitude rho(irho), at any direction alpha
    !! (the sub-grid readout): the rotation's coefficients from the shell sums by a direct transform
    module real(dp) function shiftmap_cross( self, irot, irho, alpha )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: irot, irho
        real(dp),            intent(in) :: alpha
        complex(sp) :: a(0:self%shiftmap%get_nharm_max())
        complex(dp) :: e, em, acc
        real(dp)    :: arg
        integer     :: ithr, n, m, nh
        ithr = omp_get_thread_num() + 1
        if( self%smwork(ithr)%iref < 1 ) THROW_HARD('gen_shiftmap first, on this thread; shiftmap_cross')
        nh  = self%shiftmap%get_nharm(irho)
        arg = -DTWOPI * real(irot-1, dp) / real(self%nrots, dp)
        e   = cmplx(cos(arg), sin(arg), kind=dp)
        a   = cmplx(0., 0., kind=sp)
        do n = 0, nh
            acc = cmplx(0.d0, 0.d0, kind=dp)
            em  = cmplx(1.d0, 0.d0, kind=dp)
            do m = 1, self%nrots
                acc = acc + cmplx(self%smwork(ithr)%gsum(m,n,irho), kind=dp) * em
                em  = em * e
            enddo
            a(n) = cmplx(acc / real(self%nrots, dp), kind=sp)
        enddo
        shiftmap_cross = self%shiftmap%eval_rotcoeffs(a, irho, alpha)
    end function shiftmap_cross

    !> The objective of the thread's pair on the whole grid, vals(ialpha, irot, irho): the correlation for cc, the
    !! normalised Euclidean loss L = 1 + (|CTF R|^2_w - 2 C) / wsqsum (lower is better) for euclid. The data term
    !! only: the sampler adds the untempered prior and the cell measure (shiftmap_cell_area) itself
    module subroutine shiftmap_objfun_grid( self, vals )
        class(polarft_calc), target, intent(inout) :: self
        real(sp),                    intent(out)   :: vals(:,:,:)
        type(shiftmap_work), pointer :: wk
        real(sp) :: scal(self%nrots), offs(self%nrots)
        real(dp) :: invn
        integer  :: ithr, nrho, nalpha, irho, irot
        ithr =  omp_get_thread_num() + 1
        wk   => self%smwork(ithr)
        if( wk%iref < 1 ) THROW_HARD('gen_shiftmap first, on this thread; shiftmap_objfun_grid')
        nrho   = self%shiftmap%get_nrho()
        nalpha = self%shiftmap%get_nalpha()
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
    module subroutine shiftmap_peak( self, vals, irot, xy, val, irho, ialpha )
        class(polarft_calc), intent(in)  :: self
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
        xy     = self%shiftmap%get_shift(irho, ialpha)
    end subroutine shiftmap_peak

end submodule simple_polarft_shiftmap
