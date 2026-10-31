!@descr: polarft class submodule: the polar shift map of a reference/particle pair over every rotation from the memoised spectra
! For every order n and magnitude rho, the all-rotation form of the Fourier-Bessel expansion (simple_polar_shiftmap):
!   a_n(irot, rho) = (1/nrots) sum_m exp(-2 pi i m (irot-1)/nrots) sum_k w_k J_n(z_k(rho)) R^_k[m-n] conj(Y^_k[m]),
! R^ and Y^ the full-ring spectra of the reference and of the CTF-modulated particle, from the memoised half spectra
! by F[-m] = (-1)^m conj(F[m]); one forward transform per (n, rho) and none of a shifted reference. At a grid point
! C = Re(a_0)/2 + sum_{n>=1} Re(exp(-i n alpha) a_n), the cross term the shift matrix enters every objective through.
submodule (simple_polarft_calc) simple_polarft_shiftmap
#include "simple_local_flags.inc"
implicit none

contains

    !> the tables of the shift grid (rho in pixels of the box the shift matrices act on, nalpha directions) and the
    !! per-thread work, one per thread of the memo workspace
    module subroutine new_shiftmap( self, rho, nalpha )
        class(polarft_calc), intent(inout) :: self
        real(dp),            intent(in)    :: rho(:)
        integer,             intent(in)    :: nalpha
        integer :: ithr, nmax, nrho
        call self%kill_shiftmap
        if( .not. self%existence ) THROW_HARD('the calculator must exist; new_shiftmap')
        allocate(self%shiftmap)
        call self%shiftmap%new(self%pftsz, self%kfromto, self%ldim(1), rho, nalpha)
        nmax = self%shiftmap%get_nharm_max()
        nrho = self%shiftmap%get_nrho()
        allocate(self%smwork(nthr_glob))
        do ithr = 1, nthr_glob
            allocate(self%smwork(ithr)%rhat(self%nrots, self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%yhat(self%nrots, self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%prod(self%nrots, self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%gvec(self%nrots),&
                    &self%smwork(ithr)%a(self%nrots, 0:nmax, nrho),&
                    &self%smwork(ithr)%w(self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%wj(self%kfromto(1):self%kfromto(2)),&
                    &self%smwork(ithr)%refpow(self%nrots),&
                    &self%smwork(ithr)%cmap(self%nrots, nalpha),&
                    &self%smwork(ithr)%pend(2, self%nk))
            self%smwork(ithr)%a = cmplx(0., 0., kind=sp)
        enddo
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
                        &self%smwork(ithr)%gvec, self%smwork(ithr)%a, self%smwork(ithr)%w, self%smwork(ithr)%wj,&
                        &self%smwork(ithr)%refpow, self%smwork(ithr)%cmap, self%smwork(ithr)%pend)
                endif
            enddo
            deallocate(self%smwork)
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

    !> The pair's coefficients a_n(irot, rho) for every order, magnitude and rotation, and the weighted reference power
    !! under the rotated CTF per rotation, into the calling thread's work; shiftmap_cross and shiftmap_objfun_grid
    !! read them on the same thread. The shell weights are the objective's: k for cc, k/sigma2 for euclid
    module subroutine gen_shiftmap( self, iref, iptcl )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        type(shiftmap_work), pointer :: wk
        integer :: ithr, i, ieo, ih, k, m, j, n, irho, nrho, nmax, ncol, sgn
        if( .not. allocated(self%shiftmap) ) THROW_HARD('no shift map; new_shiftmap first; gen_shiftmap')
        ithr =  omp_get_thread_num() + 1
        i    =  self%pinds(iptcl)
        ieo  =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        ih   =  merge(1, 2, self%iseven(i))
        wk   => self%smwork(ithr)
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                do k = self%kfromto(1), self%kfromto(2)
                    wk%w(k) = real(k, dp)
                enddo
                wk%ptclpow = self%ksqsums_ptcls(i)
            case(OBJFUN_EUCLID)
                if( .not. allocated(self%sigma2_wk8) ) THROW_HARD('euclid needs the sigma2 weights; gen_shiftmap')
                wk%w       = self%sigma2_wk8(:,ih)
                wk%ptclpow = self%wsqsums_ptcls(i)
            case DEFAULT
                THROW_HARD('unsupported objective; gen_shiftmap')
        end select
        wk%iref  = iref
        wk%iptcl = iptcl
        ! the weighted reference power under the rotated CTF: sum_k w_k FT(CTF2) FT(REF2) back-transformed is
        ! 2 nrots times the half-ring sum at every rotation (the reference term of the Euclidean kernels)
        self%crvec1(ithr)%c = cmplx(0., 0., kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(wk%w(k), sp) * self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
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
        ! the orders: the shell products at the order, their Bessel-weighted shell sums per magnitude, transformed in
        ! batches of nk (the width of the memo workspace's batched plan)
        nrho = self%shiftmap%get_nrho()
        nmax = self%shiftmap%get_nharm_max()
        ncol = 0
        do n = 0, nmax
            do k = self%kfromto(1), self%kfromto(2)
                do m = 0, self%nrots - 1
                    j = modulo(m - n, self%nrots) + 1
                    wk%prod(m+1,k) = cmplx(wk%rhat(j,k), kind=dp) * cmplx(wk%yhat(m+1,k), kind=dp)
                enddo
            enddo
            do irho = 1, nrho
                if( n > self%shiftmap%get_nharm(irho) ) cycle
                call self%shiftmap%shell_weights(wk%w, n, irho, wk%wj)
                wk%gvec = cmplx(0.d0, 0.d0, kind=dp)
                do k = self%kfromto(1), self%kfromto(2)
                    wk%gvec = wk%gvec + wk%wj(k) * wk%prod(:,k)
                enddo
                ncol = ncol + 1
                self%cmat2_many(ithr)%c(:,ncol) = cmplx(wk%gvec, kind=c_float_complex)
                wk%pend(1,ncol) = n
                wk%pend(2,ncol) = irho
                if( ncol == self%nk )then
                    call self%shiftmap_fft_batch(ithr, ncol)
                    ncol = 0
                endif
            enddo
        enddo
        if( ncol > 0 ) call self%shiftmap_fft_batch(ithr, ncol)
    end subroutine gen_shiftmap

    !> the forward transforms of the ncol pending shell sums of thread ithr, into the coefficients they belong to
    module subroutine shiftmap_fft_batch( self, ithr, ncol )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: ithr, ncol
        type(shiftmap_work), pointer :: wk
        integer :: icol
        wk => self%smwork(ithr)
        call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
        do icol = 1, ncol
            wk%a(:,wk%pend(1,icol),wk%pend(2,icol)) = self%cmat2_many(ithr)%c(:,icol) / real(self%nrots, sp)
        enddo
    end subroutine shiftmap_fft_batch

    !> the cross term of the thread's pair at grid rotation irot and magnitude rho(irho), at any direction alpha
    module real(dp) function shiftmap_cross( self, irot, irho, alpha )
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: irot, irho
        real(dp),            intent(in) :: alpha
        integer :: ithr
        ithr = omp_get_thread_num() + 1
        if( self%smwork(ithr)%iref < 1 ) THROW_HARD('gen_shiftmap first, on this thread; shiftmap_cross')
        shiftmap_cross = self%shiftmap%eval_rotcoeffs(self%smwork(ithr)%a(irot,:,irho), irho, alpha)
    end function shiftmap_cross

    !> The objective of the thread's pair on the whole grid, vals(irot, ialpha, irho): the correlation for cc, the
    !! normalised Euclidean loss L = 1 + (|CTF R|^2_w - 2 C) / wsqsum (lower is better) for euclid. With a prior
    !! (prior_w > 0, u0 the centre as a stored particle-frame shift, the shift search's convention) the loss gains
    !! prior_w |x - R(irot) u0|^2 and the correlation the factor exp(-prior_w |x - R(irot) u0|^2)
    module subroutine shiftmap_objfun_grid( self, vals, prior_w, u0 )
        class(polarft_calc), target, intent(inout) :: self
        real(sp),                    intent(out)   :: vals(:,:,:)
        real(dp),          optional, intent(in)    :: prior_w, u0(2)
        type(shiftmap_work), pointer :: wk
        real(dp), allocatable :: xy0(:,:)
        real(dp) :: ang, c, s, x(2), pen, inv
        integer  :: ithr, nrho, nalpha, irho, j, irot
        logical  :: l_prior, l_cc
        ithr =  omp_get_thread_num() + 1
        wk   => self%smwork(ithr)
        if( wk%iref < 1 ) THROW_HARD('gen_shiftmap first, on this thread; shiftmap_objfun_grid')
        nrho   = self%shiftmap%get_nrho()
        nalpha = self%shiftmap%get_nalpha()
        if( size(vals,1) /= self%nrots .or. size(vals,2) /= nalpha .or. size(vals,3) /= nrho )then
            THROW_HARD('vals must be (nrots, nalpha, nrho); shiftmap_objfun_grid')
        endif
        l_cc    = self%p_ptr%cc_objfun == OBJFUN_CC
        l_prior = .false.
        if( present(prior_w) )then
            if( .not. present(u0) ) THROW_HARD('a prior needs its centre; shiftmap_objfun_grid')
            l_prior = prior_w > 0.d0
        endif
        if( l_prior )then
            ! the centre in the optimizer's frame at every rotation: R(irot) u0, R = [[c, s], [-s, c]]
            allocate(xy0(2, self%nrots))
            do irot = 1, self%nrots
                ang = real(self%angtab(irot), dp) * DPI / 180.d0
                c   = cos(ang)
                s   = sin(ang)
                xy0(1,irot) =  c * u0(1) + s * u0(2)
                xy0(2,irot) = -s * u0(1) + c * u0(2)
            enddo
        endif
        inv = 1.d0 / wk%ptclpow
        do irho = 1, nrho
            call self%shiftmap%map_rotcoeffs(wk%a(:,:,irho), irho, wk%cmap)
            do j = 1, nalpha
                x = self%shiftmap%get_shift(irho, j)
                do irot = 1, self%nrots
                    if( l_cc )then
                        vals(irot,j,irho) = real(wk%cmap(irot,j) / sqrt(wk%refpow(irot) * wk%ptclpow), sp)
                    else
                        vals(irot,j,irho) = real(1.d0 + (wk%refpow(irot) - 2.d0 * wk%cmap(irot,j)) * inv, sp)
                    endif
                    if( l_prior )then
                        pen = prior_w * sum((x - xy0(:,irot))**2)
                        if( l_cc )then
                            vals(irot,j,irho) = vals(irot,j,irho) * real(exp(-pen), sp)
                        else
                            vals(irot,j,irho) = vals(irot,j,irho) + real(pen, sp)
                        endif
                    endif
                enddo
            enddo
        enddo
        if( l_prior ) deallocate(xy0)
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
        irot   = loc(1)
        ialpha = loc(2)
        irho   = loc(3)
        val    = vals(irot, ialpha, irho)
        xy     = self%shiftmap%get_shift(irho, ialpha)
    end subroutine shiftmap_peak

end submodule simple_polarft_shiftmap
