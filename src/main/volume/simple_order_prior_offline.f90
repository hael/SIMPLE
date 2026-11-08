!@descr: offline validation of the PCG order prior on a half-map pair: its closed-form choice against full surrogate solves
! Per half: (F^-1 D F + lambda h diag(1 - w_order)) x = F^-1 D F x0, h the mean of D(k) = 1/sigma_n^2(k) of the pair's
! difference spectrum, the Fourier-diagonal stand-in for the PCG operator (constant D: the closed form). The
! closed-form choice is the loop's (select_order_prior); the same grid solved for real is scored by the same rule.
! Design: doc/implementation_notes/planned/pcg_order_prior.md.
module simple_order_prior_offline
use, intrinsic :: iso_c_binding, only: c_float, c_float_complex
use simple_defs,            only: dp, TINY, logfhandle
use simple_error,           only: simple_exception
use simple_math_ft,         only: fdim, get_resarr
use simple_estimate_ssnr,   only: get_resolution
use simple_image,           only: image
use simple_pcg_solver,      only: pcg_operator, pcg_solver_options, pcg_solver_outcome, pcg_solve
use simple_pcg_order_prior, only: build_density_weight, order_regions, build_order_weights, select_order_prior, &
    &order_prior_stats, choose_order_setting, ORDER_FRAC_GRID, ORDER_LAMBDA_GRID, CORE_SIGNAL_MIN, CORE_MIN_VOX
implicit none

public :: order_prior_offline, run_all_order_prior_offline_tests
private
#include "simple_local_flags.inc"

integer, parameter :: DATA_SCALE_SHELLS = 6      !< low band of the data scale, as reconstructor_pcg
real,    parameter :: SUPPORT_MIN       = 0.5    !< soft-sphere value taken as the evaluation domain
real,    parameter :: EMPTY_SHELL_REL   = 1.e-8  !< difference power below this fraction of the peak carries no data
real,    parameter :: EMPTY_SHELL_FAC   = 100.   !< precision of shells without data, relative to the largest measured
real,    parameter :: CG_RTOL           = 1.e-3
real,    parameter :: CG_XTOL           = 1.e-6  !< below any useful step: the residual decides
integer, parameter :: CG_MAXITS         = 80

type :: order_prior_offline
    private
    type(image)           :: even, odd       !< the prior-free pair
    type(image)           :: support         !< soft mskdiam sphere
    type(image)           :: wimg            !< transform work image
    real,     allocatable :: dgrid(:,:,:)    !< D on the stored half of the Fourier lattice
    real,     allocatable :: dshell(:)       !< D per shell
    logical,  allocatable :: lsup(:,:,:)     !< evaluation domain
    integer :: ldim(3) = 0
    real    :: smpd = 0., res0143 = 0., fsc05 = 0., data_scale = 0., hdiag = 0.
    logical :: l_exists = .false.
  contains
    procedure          :: new
    procedure          :: get_res0143
    procedure          :: density_weights
    procedure          :: run
    procedure          :: kill
    procedure, private :: noise_model
    procedure, private :: fourier_mul
    procedure, private :: closed_form
    procedure, private :: init_op
    procedure, private :: solve_half
end type order_prior_offline

!> one half's surrogate operator for the shared CG engine
type, extends(pcg_operator) :: surrogate_op
    class(order_prior_offline), pointer :: client => null()
    real, allocatable :: pen(:,:,:)       !< lambda_s (1 - w)
    real, allocatable :: sqw(:,:,:)       !< sqrt(w): ordered part of the preconditioner, D^-1
    real, allocatable :: sqs(:,:,:)       !< sqrt(1 - w): ridge part, (D + lambda_s)^-1
    real, allocatable :: t(:,:,:), u(:,:,:)
    real              :: lam_s = 0.
  contains
    procedure :: size    => op_size
    procedure :: apply   => op_apply
    procedure :: precond => op_precond
    procedure :: dot     => op_dot
end type surrogate_op

interface
    !> white-box tests (submodule simple_order_prior_offline_tester, built only with BUILD_TESTS=ON)
    module subroutine run_all_order_prior_offline_tests()
    end subroutine run_all_order_prior_offline_tests
end interface

contains

    !> the prior-free pair, its FSC, the evaluation domain (mskdiam sphere) and the noise model
    subroutine new( self, even, odd, mskdiam )
        class(order_prior_offline), target, intent(inout) :: self
        class(image),                       intent(in)    :: even, odd
        real,                               intent(in)    :: mskdiam
        type(image)            :: e_ft, o_ft
        real, allocatable      :: corrs(:), res(:), ones(:,:,:)
        real(c_float), pointer :: rsup(:,:,:)
        real    :: msk
        integer :: n
        call self%kill
        self%ldim = even%get_ldim()
        self%smpd = even%get_smpd()
        n = self%ldim(1)
        if( any(odd%get_ldim() /= self%ldim) ) THROW_HARD('the halves differ in dimensions; order_prior_offline%new')
        if( any(self%ldim /= n) )              THROW_HARD('cubic volumes only; order_prior_offline%new')
        if( mod(n,2) /= 0 )                    THROW_HARD('even box only; order_prior_offline%new')
        if( even%is_ft() .or. odd%is_ft() )    THROW_HARD('the halves must be in real space; order_prior_offline%new')
        msk = 0.5 * mskdiam / self%smpd
        if( msk < 4. .or. msk > real(n/2) )    THROW_HARD('mskdiam must fit inside the box; order_prior_offline%new')
        call self%even%copy(even)
        call self%odd%copy(odd)
        call e_ft%copy(even)
        call o_ft%copy(odd)
        call e_ft%fft
        call o_ft%fft
        allocate(corrs(fdim(n)-1), source=0.)
        call e_ft%fsc(o_ft, corrs)
        call e_ft%kill
        call o_ft%kill
        res = get_resarr(n, self%smpd)
        call get_resolution(corrs, res, self%fsc05, self%res0143)
        if( self%res0143 <= 2.*self%smpd + TINY ) THROW_HARD('the half-map FSC never falls below 0.143: are the halves independent?')
        allocate(ones(n,n,n), source=1.)
        call self%support%new(self%ldim, self%smpd)
        call self%support%set_rmat(ones, .false.)
        call self%support%mask3D_soft(msk, backgr=0.)
        call self%support%get_rmat_ptr(rsup)
        allocate(self%lsup(n,n,n))
        self%lsup = rsup(:n,:n,:n) >= SUPPORT_MIN
        nullify(rsup)
        call self%noise_model
        self%l_exists = .true.
        write(logfhandle,'(A,I0,A,F7.3,A,F8.3,A,F8.3,A)') '>>> ORDER PRIOR: box ', n, ', smpd ', self%smpd, &
            &' A, half-map FSC=0.5 ', self%fsc05, ' A, FSC=0.143 ', self%res0143, ' A'
        deallocate(corrs, res, ones)
    end subroutine new

    pure real function get_res0143( self )
        class(order_prior_offline), intent(in) :: self
        get_res0143 = self%res0143
    end function get_res0143

    !> the loop's per-half density weights
    subroutine density_weights( self, w )
        class(order_prior_offline), intent(inout) :: self
        type(image),                intent(inout) :: w(2)
        real :: solv(2)
        call build_density_weight(self%even, self%res0143, self%support, w(1), solv(1))
        call build_density_weight(self%odd,  self%res0143, self%support, w(2), solv(2))
        write(logfhandle,'(A,2F7.2,A,F6.3)') '>>> ORDER PRIOR: density weights, solvent % of the support even/odd ', &
            &100.*solv, ', even/odd correlation ', w(1)%real_corr(w(2))
    end subroutine density_weights

    !> the loop's closed-form choice, the same grid solved through the operator and chosen by the same rule, and
    !! the operator pair at the closed-form choice (what the loop ships; the input pair when it is none)
    subroutine run( self, evmap, locres, w_dens, lambda_fixed, w_order, ce, co )
        class(order_prior_offline), target, intent(inout) :: self
        class(image),                       intent(in)    :: evmap, locres
        type(image),                        intent(in)    :: w_dens(2)
        real,                               intent(in)    :: lambda_fixed
        type(image),                        intent(inout) :: w_order(2), ce, co
        integer, parameter :: NF = size(ORDER_FRAC_GRID)
        type(order_prior_stats)   :: st
        type(image)               :: wf(2)
        type(pcg_solver_outcome)  :: oe, oo
        real(c_float), pointer    :: rp(:,:,:)
        real, allocatable, target :: xe(:,:,:), xo(:,:,:), x0e(:,:,:), x0o(:,:,:)
        real,     allocatable     :: lambdas(:), rgrid(:,:), kgrid(:,:)
        logical,  allocatable     :: lcore(:,:,:), lperi(:,:,:)
        real(dp) :: core0, peri0
        logical  :: l_valid(NF)
        integer  :: n, nl, ifr, il, fsel, lsel, icf, lcf
        n = self%ldim(1)
        if( lambda_fixed > 0. )then
            lambdas = [lambda_fixed]
        else
            lambdas = ORDER_LAMBDA_GRID
        endif
        nl = size(lambdas)
        ! the loop's choice, in closed form with the stand-in's diagonal
        call select_order_prior([self%even, self%odd], w_dens, evmap, locres, self%res0143, self%lsup, &
            &spread(self%dshell, 2, 2), [self%hdiag, self%hdiag], lambda_fixed, w_order, st, .true.)
        call print_choice('CLOSED FORM', st%frac, st%lambda_rel, st%removed, st%core_kept, st%l_off)
        write(logfhandle,'(A,F6.2,A,F6.2,A)') '>>> ORDER PRIOR: protein core ', st%pct_core, ' % of the support, dense around it ', &
            &st%pct_peri, ' %'
        ! the same grid through the operator, scored on the solutions' cross-half signal
        call order_regions(w_dens, locres, self%res0143, self%lsup, lcore, lperi)
        call self%even%get_rmat_ptr(rp)
        allocate(x0e(n,n,n), source=rp(:n,:n,:n))
        call self%odd%get_rmat_ptr(rp)
        allocate(x0o(n,n,n), source=rp(:n,:n,:n))
        nullify(rp)
        core0 = cross_signal(x0e, x0o, lcore)
        peri0 = cross_signal(x0e, x0o, lperi)
        allocate(rgrid(NF,nl), kgrid(NF,nl), source=0.)
        allocate(xe(n,n,n), xo(n,n,n))
        do ifr = 1, NF
            call build_order_weights(evmap, self%lsup, ORDER_FRAC_GRID(ifr), wf, l_valid(ifr))
            if( .not. l_valid(ifr) ) cycle
            xe = x0e
            xo = x0o
            do il = 1, nl
                call self%solve_half(1, wf(1), lambdas(il), xe, oe)
                call self%solve_half(2, wf(2), lambdas(il), xo, oo)
                call oe%kill
                call oo%kill
                if( core0 > 0._dp ) kgrid(ifr,il) = real(cross_signal(xe, xo, lcore) / core0)
                if( peri0 > 0._dp ) rgrid(ifr,il) = real(1._dp - cross_signal(xe, xo, lperi) / peri0)
            end do
        end do
        ! the same rule on the operator's numbers
        call choose_order_setting(kgrid, rgrid, l_valid, count(lcore) >= CORE_MIN_VOX .and. core0 > 0._dp, &
            &count(lperi) >= CORE_MIN_VOX .and. peri0 > 0._dp, lambda_fixed > 0., fsel, lsel)
        call print_tables
        if( fsel == 0 )then
            call print_choice('OPERATOR', 0., 0., 0., 1., .true.)
        else
            call print_choice('OPERATOR', ORDER_FRAC_GRID(fsel), lambdas(lsel), rgrid(fsel,lsel), kgrid(fsel,lsel), .false.)
        endif
        icf = 0
        lcf = 0
        if( .not. st%l_off )then
            icf = minloc(abs(ORDER_FRAC_GRID - st%frac), dim=1)
            lcf = minloc(abs(lambdas - st%lambda_rel), dim=1)
            write(logfhandle,'(A,F6.1,A,F6.1,A)') '>>> ORDER PRIOR: the operator at the closed-form choice removes ', &
                &100.*rgrid(icf,lcf), ' % and keeps ', 100.*kgrid(icf,lcf), ' % of the core'
        endif
        if( icf == fsel .and. lcf == lsel )then
            write(logfhandle,'(A)') '>>> ORDER PRIOR: the closed form and the operator agree'
        else
            write(logfhandle,'(A)') '>>> ORDER PRIOR: the closed form and the operator disagree'
        endif
        ! the operator pair at the closed-form choice, as the loop's re-solve
        call self%even%get_rmat_ptr(rp)
        xe = rp(:n,:n,:n)
        call self%odd%get_rmat_ptr(rp)
        xo = rp(:n,:n,:n)
        nullify(rp)
        if( .not. st%l_off )then
            call self%solve_half(1, w_order(1), st%lambda_rel, xe, oe)
            call self%solve_half(2, w_order(2), st%lambda_rel, xo, oo)
            call oe%kill
            call oo%kill
        endif
        call ce%new(self%ldim, self%smpd)
        call co%new(self%ldim, self%smpd)
        call ce%set_rmat(xe, .false.)
        call co%set_rmat(xo, .false.)
        call wf(1)%kill
        call wf(2)%kill
        deallocate(lambdas, rgrid, kgrid, lcore, lperi, xe, xo, x0e, x0o)

    contains

        subroutine print_tables
            integer :: jf, jl
            write(logfhandle,'(A)') '>>> ORDER PRIOR: operator, % of the dense cross-half signal around the core removed,'//&
                &' rows ordered fraction, columns lambda_rel; * keeps 95 % of the core, < chosen'
            write(logfhandle,'(A10,*(F9.2))') 'frac', lambdas
            do jf = 1, NF
                if( .not. l_valid(jf) ) cycle
                write(logfhandle,'(F10.3,*(F8.1,A1))') ORDER_FRAC_GRID(jf), (100.*rgrid(jf,jl), mark(jf,jl), jl=1,nl)
            end do
            write(logfhandle,'(A)') '>>> ORDER PRIOR: operator, % of the core cross-half signal kept'
            write(logfhandle,'(A10,*(F9.2))') 'frac', lambdas
            do jf = 1, NF
                if( .not. l_valid(jf) ) cycle
                write(logfhandle,'(F10.3,*(F8.1,A1))') ORDER_FRAC_GRID(jf), (100.*kgrid(jf,jl), mark(jf,jl), jl=1,nl)
            end do
        end subroutine print_tables

        character(len=1) function mark( jf, jl )
            integer, intent(in) :: jf, jl
            mark = ' '
            if( kgrid(jf,jl) >= CORE_SIGNAL_MIN ) mark = '*'
            if( jf == fsel .and. jl == lsel ) mark = '<'
        end function mark

        !> sum of even x odd over a mask: the independent noise of the halves cancels
        real(dp) function cross_signal( a, b, l ) result( cs )
            real,    intent(in) :: a(:,:,:), b(:,:,:)
            logical, intent(in) :: l(:,:,:)
            integer :: i, j, k
            cs = 0._dp
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        if( l(i,j,k) ) cs = cs + real(a(i,j,k),dp) * real(b(i,j,k),dp)
                    end do
                end do
            end do
        end function cross_signal

    end subroutine run

    subroutine print_choice( tag, frac, lambda_rel, removed, core_kept, l_off )
        character(len=*), intent(in) :: tag
        real,             intent(in) :: frac, lambda_rel, removed, core_kept
        logical,          intent(in) :: l_off
        if( l_off )then
            write(logfhandle,'(A,A,A)') '>>> ORDER PRIOR CHOSEN, ', tag, ': none, no setting keeps 95 % of the core'
        else
            write(logfhandle,'(A,A,A,F6.2,A,F7.2,A,F6.1,A,F6.1,A)') '>>> ORDER PRIOR CHOSEN, ', tag, ': ordered fraction ', &
                &100.*frac, ' %, lambda_rel ', lambda_rel, ', removes ', 100.*removed, ' %, keeps ', 100.*core_kept, &
                &' % of the core'
        endif
    end subroutine print_choice

    !> the per-voxel shrink of one half (1 even, 2 odd): x0 h / (h + lambda_s (1-w))
    subroutine closed_form( self, ihalf, w, lam_rel, out )
        class(order_prior_offline), target, intent(inout) :: self
        integer,                            intent(in)    :: ihalf
        class(image),                       intent(in)    :: w
        real,                               intent(in)    :: lam_rel
        class(image),                       intent(inout) :: out
        real(c_float), pointer :: rx(:,:,:), rw(:,:,:), ro(:,:,:)
        real    :: lam_s
        integer :: n
        n     = self%ldim(1)
        lam_s = lam_rel * self%hdiag
        if( ihalf == 1 )then
            call self%even%get_rmat_ptr(rx)
        else
            call self%odd%get_rmat_ptr(rx)
        endif
        call w%get_rmat_ptr(rw)
        call out%get_rmat_ptr(ro)
        ro(:n,:n,:n) = rx(:n,:n,:n) * self%hdiag / (self%hdiag + lam_s * (1. - min(1., max(0., rw(:n,:n,:n)))))
        nullify(rx, rw, ro)
    end subroutine closed_form

    !> the ridge and the two-region preconditioner of one half's operator
    subroutine init_op( self, op, w, lam_rel )
        class(order_prior_offline), target, intent(inout) :: self
        type(surrogate_op),                 intent(inout) :: op
        class(image),                       intent(in)    :: w
        real,                               intent(in)    :: lam_rel
        real(c_float), pointer :: rw(:,:,:)
        real    :: wv
        integer :: n, i, j, k
        n = self%ldim(1)
        op%client => self
        op%lam_s  = lam_rel * self%hdiag
        allocate(op%pen(n,n,n), op%sqw(n,n,n), op%sqs(n,n,n), op%t(n,n,n), op%u(n,n,n))
        call w%get_rmat_ptr(rw)
        !$omp parallel do collapse(3) default(shared) private(i,j,k,wv) schedule(static) proc_bind(close)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    wv = min(1., max(0., rw(i,j,k)))
                    op%pen(i,j,k) = op%lam_s * (1. - wv)
                    op%sqw(i,j,k) = sqrt(wv)
                    op%sqs(i,j,k) = sqrt(1. - wv)
                end do
            end do
        end do
        !$omp end parallel do
        nullify(rw)
    end subroutine init_op

    !> one half (1 even, 2 odd) through the surrogate operator; x holds the start and returns the solution
    subroutine solve_half( self, ihalf, w, lam_rel, x, outcome )
        class(order_prior_offline), target, intent(inout) :: self
        integer,                            intent(in)    :: ihalf
        class(image),                       intent(in)    :: w
        real,                               intent(in)    :: lam_rel
        real, contiguous, target,           intent(inout) :: x(:,:,:)
        type(pcg_solver_outcome),           intent(out)   :: outcome
        type(surrogate_op)        :: op
        type(pcg_solver_options)  :: opts
        real(c_float), pointer    :: rx0(:,:,:)
        real, pointer, contiguous :: bf(:), xf(:)
        real, allocatable, target :: b(:,:,:)
        integer :: n
        n = self%ldim(1)
        if( ihalf == 1 )then
            call self%even%get_rmat_ptr(rx0)
        else
            call self%odd%get_rmat_ptr(rx0)
        endif
        call self%init_op(op, w, lam_rel)
        if( maxval(op%pen) <= 0. )then
            ! no ridge anywhere: the prior-free half is the solution
            x = rx0(:n,:n,:n)
            outcome%stop_reason = 'no_prior'
            outcome%converged   = .true.
            nullify(rx0)
            return
        endif
        allocate(b(n,n,n))
        call self%fourier_mul(rx0(:n,:n,:n), b, .false., 0.)
        nullify(rx0)
        bf(1:n**3) => b
        xf(1:n**3) => x
        opts%maxits        = CG_MAXITS
        opts%rtol          = CG_RTOL
        opts%xtol          = CG_XTOL
        opts%rescale_start = .true.
        opts%cold_restart  = .true.
        opts%track_preconditioned_residual = .false.
        opts%record_history = .false.
        opts%tag           = 'ORDER PRIOR CG'
        call pcg_solve(op, bf, xf, opts, outcome)
        nullify(bf, xf)
        deallocate(b)
    end subroutine solve_half

    !> out = F^-1 g F v with g = D, or g = 1/(D + shift) when l_inv
    subroutine fourier_mul( self, v, out, l_inv, shift )
        class(order_prior_offline), intent(inout) :: self
        real,                       intent(in)    :: v(:,:,:)
        real,                       intent(out)   :: out(:,:,:)
        logical,                    intent(in)    :: l_inv
        real,                       intent(in)    :: shift
        real(c_float),            pointer :: rp(:,:,:)
        complex(c_float_complex), pointer :: cp(:,:,:)
        integer :: n, cdim(3), i, j, k
        n = self%ldim(1)
        call self%wimg%get_rmat_ptr(rp)
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    rp(i,j,k) = v(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
        nullify(rp)
        call self%wimg%fft
        call self%wimg%get_cmat_ptr(cp)
        cdim = self%wimg%get_array_shape()
        if( l_inv )then
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
            do k = 1, cdim(3)
                do j = 1, cdim(2)
                    do i = 1, cdim(1)
                        cp(i,j,k) = cp(i,j,k) / (self%dgrid(i,j,k) + shift)
                    end do
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
            do k = 1, cdim(3)
                do j = 1, cdim(2)
                    do i = 1, cdim(1)
                        cp(i,j,k) = cp(i,j,k) * self%dgrid(i,j,k)
                    end do
                end do
            end do
            !$omp end parallel do
        endif
        nullify(cp)
        call self%wimg%ifft
        call self%wimg%get_rmat_ptr(rp)
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    out(i,j,k) = rp(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
        nullify(rp)
    end subroutine fourier_mul

    !> D(k) = 1/sigma_n^2(k), the per-half noise variance per voxel from the pair's difference spectrum (fft divides
    !! by N: sigma_n^2 = N <|F(E-O)|^2>/2); shells without data (a Fourier-padded pair) are held at zero
    subroutine noise_model( self )
        class(order_prior_offline), intent(inout) :: self
        type(image)       :: dimg
        real, allocatable :: spec(:)
        real(dp) :: ds_sum, h_sum
        real     :: smax, dmax, nvox
        integer  :: n, nyq, nmax, nlast, sh, a, b, c, cdim(3), logi(3), ds_n, h_n
        n   = self%ldim(1)
        nyq = fdim(n) - 1
        call dimg%copy(self%even)
        call dimg%subtr(self%odd)
        call dimg%fft
        allocate(spec(nyq), source=0.)
        call dimg%power_spectrum(spec)
        call dimg%kill
        smax = maxval(spec)
        if( smax <= 0. ) THROW_HARD('the halves are identical: no noise to model; order_prior_offline%noise_model')
        nlast = 0
        do sh = 1, nyq
            if( spec(sh) > EMPTY_SHELL_REL * smax ) nlast = sh
        end do
        nmax = ceiling(sqrt(3.) * real(nyq)) + 1
        allocate(self%dshell(0:nmax), source=0.)
        nvox = real(n)**3
        do sh = 1, nlast
            self%dshell(sh) = 2. / (nvox * max(spec(sh), EMPTY_SHELL_REL * smax))
        end do
        self%dshell(0) = self%dshell(1)
        dmax = maxval(self%dshell(1:nlast))
        if( nlast == nyq )then
            self%dshell(nyq+1:) = self%dshell(nyq)
        else
            self%dshell(nlast+1:) = EMPTY_SHELL_FAC * dmax
        endif
        call self%wimg%new(self%ldim, self%smpd)
        call self%wimg%zero_and_unflag_ft
        cdim = self%wimg%get_array_shape()
        allocate(self%dgrid(cdim(1),cdim(2),cdim(3)))
        ds_sum = 0._dp
        h_sum  = 0._dp
        ds_n   = 0
        h_n    = 0
        do c = 1, cdim(3)
            do b = 1, cdim(2)
                do a = 1, cdim(1)
                    logi = self%wimg%comp_addr_logi(a,b,c)
                    sh   = min(nmax, nint(sqrt(real(sum(logi**2)))))
                    self%dgrid(a,b,c) = self%dshell(sh)
                    if( sh <= min(DATA_SCALE_SHELLS, nlast) )then
                        ds_sum = ds_sum + real(self%dshell(sh),dp)
                        ds_n   = ds_n + 1
                    endif
                    if( sh <= nlast )then
                        h_sum = h_sum + real(self%dshell(sh),dp)
                        h_n   = h_n + 1
                    endif
                end do
            end do
        end do
        self%data_scale = real(ds_sum / real(ds_n,dp))
        self%hdiag      = real(h_sum  / real(h_n,dp))
        write(logfhandle,'(A)') '>>> ORDER PRIOR: NOISE MODEL D(k) = 1/sigma_n^2(k) FROM THE HALF-PAIR DIFFERENCE SPECTRUM'
        write(logfhandle,'(A,ES10.3,A,ES10.3,A,F8.2)') '    data scale (shells 0-6) ', self%data_scale, &
            &', real-space diagonal ', self%hdiag, ', ratio ', self%data_scale / self%hdiag
        write(logfhandle,'(A)') '        shell     res(A)  D/data_scale'
        sh = 1
        do while( sh <= nlast )
            write(logfhandle,'(I13,F11.2,ES14.3)') sh, real(n) * self%smpd / real(sh), self%dshell(sh) / self%data_scale
            sh = 2 * sh
        end do
        if( nlast < nyq ) write(logfhandle,'(A,I0,A,F7.2,A)') '>>> ORDER PRIOR: no data beyond shell ', nlast, ' (', &
            &real(n) * self%smpd / real(nlast), ' A), a Fourier-padded pair; those shells are held at zero'
        deallocate(spec)
    end subroutine noise_model

    subroutine kill( self )
        class(order_prior_offline), intent(inout) :: self
        call self%even%kill
        call self%odd%kill
        call self%support%kill
        call self%wimg%kill
        if( allocated(self%dgrid)  ) deallocate(self%dgrid)
        if( allocated(self%dshell) ) deallocate(self%dshell)
        if( allocated(self%lsup)   ) deallocate(self%lsup)
        self%ldim       = 0
        self%smpd       = 0.
        self%res0143    = 0.
        self%fsc05      = 0.
        self%data_scale = 0.
        self%hdiag      = 0.
        self%l_exists   = .false.
    end subroutine kill

    ! ---------------- surrogate operator ----------------

    integer function op_size( self ) result( nsz )
        class(surrogate_op), intent(in) :: self
        nsz = product(self%client%ldim)
    end function op_size

    !> y = F^-1 D F x + lambda_s (1-w) x
    subroutine op_apply( self, x, y )
        class(surrogate_op),      intent(inout) :: self
        real, contiguous, target, intent(in)    :: x(:)
        real, contiguous, target, intent(out)   :: y(:)
        real, pointer :: x3(:,:,:), y3(:,:,:)
        integer :: n, i, j, k
        n = self%client%ldim(1)
        x3(1:n,1:n,1:n) => x
        y3(1:n,1:n,1:n) => y
        call self%client%fourier_mul(x3, y3, .false., 0.)
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    y3(i,j,k) = y3(i,j,k) + self%pen(i,j,k) * x3(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
        nullify(x3, y3)
    end subroutine op_apply

    !> sqrt(w) D^-1 sqrt(w) + sqrt(1-w) (D + lambda_s)^-1 sqrt(1-w): SPD as w + (1-w) = 1, exact where w is 0 or 1
    subroutine op_precond( self, r, z )
        class(surrogate_op),      intent(inout) :: self
        real, contiguous, target, intent(in)    :: r(:)
        real, contiguous, target, intent(out)   :: z(:)
        real, pointer :: r3(:,:,:), z3(:,:,:)
        integer :: n
        n = self%client%ldim(1)
        r3(1:n,1:n,1:n) => r
        z3(1:n,1:n,1:n) => z
        z3 = 0.
        call add_part(self%sqw, 0.)
        call add_part(self%sqs, self%lam_s)
        nullify(r3, z3)

    contains

        subroutine add_part( sq, c )
            real, intent(in) :: sq(:,:,:)
            real, intent(in) :: c
            integer :: i, j, k
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        self%t(i,j,k) = sq(i,j,k) * r3(i,j,k)
                    end do
                end do
            end do
            !$omp end parallel do
            call self%client%fourier_mul(self%t, self%u, .true., c)
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        z3(i,j,k) = z3(i,j,k) + sq(i,j,k) * self%u(i,j,k)
                    end do
                end do
            end do
            !$omp end parallel do
        end subroutine add_part

    end subroutine op_precond

    real(dp) function op_dot( self, a, b ) result( val )
        class(surrogate_op),      intent(in) :: self
        real, contiguous, target, intent(in) :: a(:), b(:)
        integer :: i
        val = 0._dp
        !$omp parallel do default(shared) private(i) reduction(+:val) schedule(static) proc_bind(close)
        do i = 1, size(a)
            val = val + real(a(i),dp) * real(b(i),dp)
        end do
        !$omp end parallel do
    end function op_dot

end module simple_order_prior_offline
