!@descr: atomic density model of SINGLE: Gaussian atoms through the measured transfer of the reconstruction, their fit, profiles in the coordination number, noise and rendering
! The design and its notation are doc/implementation_notes/planned/atomic_density_model.md (sections 2, 4 and 12).
! Positions are voxel coordinates (voxel (i,j,k) sits at (i,j,k)), widths B in A**2, intensities are integrals.
! The transfer T(s) is a radial curve per Fourier shell, measured on the map as the projection onto the atoms
! without it and smoothed; an atom through it is K_T(B) = IFT[exp(-B s**2/4) T(s)] = K(B) - h(B), with the halo
! h(B) = IFT[exp(-B s**2/4) (1 - T(s))] tabulated by the radial transform (transfer_kernels). The model is kept as
! two fields: the atoms without the transfer, S = sum_i I_i K(B_i), rendered atom by atom within a few widths of
! each, and their halo H = IFT[F_S (1 - T)], one Fourier transform of S each way (the halo's tails reach across
! most of the box); model = S - H. The per-atom fit reads the table within its sphere and its cell.
! A sweep is Gauss-Seidel over the atoms, cell by cell of a cubic grid coloured by the parity of the cell indices:
! the cells of one colour are fitted in parallel against the fields as they stand (nothing is written meanwhile), the
! atoms of a cell one after the other, each seeing its own atom and halo added back exactly and the changes (atom and
! halo) of the atoms of its cell fitted before it; before the next colour the atoms are rendered again and the halo
! transformed from them. Cells of one colour are a cell apart, farther than a fit sphere reaches, so the result does
! not depend on the thread count beyond the rounding of the halo's threaded transforms; the halos of the other cells
! of a colour enter at the colour's refresh, which leaves the fixed point as it is.
module simple_nano_density
use simple_defs,  only: dp, DPI, PI, nthr_glob, logfhandle
use simple_error, only: simple_exception
implicit none

public :: nano_density, transfer_kernels
public :: gauss_kernel, halo_kernel, render_atoms, fit_atom, regress_profile
public :: gauss_filter3D, robust_spread
public :: CN_FULL
private
#include "simple_local_flags.inc"

integer, parameter :: CN_FULL          = 12     ! close-packed contact number
real,    parameter :: CONTACT_NN_FRAC  = 1.267  ! bond cutoff, in d_NN
real,    parameter :: FIT_NN_FRAC      = 0.5    ! radius of an atom's fit sphere, in d_NN
real,    parameter :: MOVE_NN_FRAC     = 0.35   ! largest move of an atom from its start, in d_NN
real,    parameter :: BSTART_NN_FRAC   = 0.2    ! start width: a Gaussian of this standard deviation, in d_NN
real,    parameter :: BFAC_MIN         = 1.     ! range of the width search, A**2
real,    parameter :: BFAC_MAX         = 400.
integer, parameter :: NGRID_B          = 40     ! log-spaced grids before golden-section refinement
integer, parameter :: NGOLD_B          = 12     ! 0.618**12 of a grid step: B to 0.3 percent
integer, parameter :: NSWEEPS          = 3
! rounds of the measured transfer (atomic_density_model.md, section 12, the sixth and seventh rulings): measure, hold,
! fit; they stop when the curve changes by CURVE_TOL or less over the reliable shells from one round to the next;
! NROUNDS_MAX is a safety limit only. The last round is kept unless its residual over the noise is above the best
! round's by more than ROUND_TOL of it
integer, parameter :: NROUNDS_MAX      = 6
real,    parameter :: CURVE_TOL        = 0.01
real,    parameter :: ROUND_TOL        = 1.e-3
integer, parameter :: NROUND_FIELDS    = 6      ! per round: residual over noise, ln S, ln B_p, beta, gamma, curve change
real,    parameter :: TRANSFER_DS      = 0.02   ! width of the transfer's Fourier shells, 1/A
real,    parameter :: POWER_FRAC       = 0.01   ! a shell is reliable while the atoms' power in it is above this of the peak
real,    parameter :: TAB_DR_VOX       = 0.1    ! radial step of the kernel table, voxels
real,    parameter :: TAB_R_MAX        = 16.    ! radial extent of the kernel table, A
integer, parameter :: NTAB_B           = 161    ! widths of the kernel table, log-spaced over [BFAC_MIN, BFAC_MAX]
real,    parameter :: TAB_DS_INT       = 0.002  ! step of the radial transform's integral, 1/A
real,    parameter :: REACH_FRAC       = 1.e-4  ! a halo is rendered from the table while |K_T| is above this of its peak
integer, parameter :: FFT_NTHR_MAX     = 16     ! threads of the halo's transforms: the system FFTW of the run machine
                                                 ! is many times slower above 16 threads (atomic_density_model.md, phase 5)
! the gauge of the measured transfer: the map fixes the model I exp(-B s**2/4) T(s) only up to T -> T / (c exp(g s**2)),
! I -> c I, B -> B - 4 g; T is held at one, in scale and in Gaussian tilt, over this band, where the pooled transfer of
! the 66 real Pt maps is unity (atomic_density_model.md, section 12 and Progress of phase 5)
real,    parameter :: GAUGE_RES_LO     = 3.1    ! A
real,    parameter :: GAUGE_RES_HI     = 2.0    ! A
real,    parameter :: PRIOR_SD_START   = 1.     ! spread of ln B around the start width before any profile
real,    parameter :: TAU_B_MIN        = 0.1    ! floor of the spread of the width prior
real,    parameter :: STRONG_SNR       = 10.    ! profile atoms: intensity above this many phantom intensity spreads
integer, parameter :: NBIN_MIN         = 20     ! atoms per coordination bin of the intensity spread
integer, parameter :: NPHANTOM         = 1000   ! phantom sites of the intensity noise
real,    parameter :: RENDER_NSD       = 5.     ! render cutoff of a Gaussian, in standard deviations (5: the atoms
                                                 ! rendered without the transfer match the analytic kernel to 1e-5)
real(dp), parameter :: GOLD = (sqrt(5._dp) - 1._dp) / 2._dp
logical,  parameter :: L_CHECK_RESIDUAL = .false. ! diagnostic: every fitted atom's residual against brute force (1 thread)

! the halo h(r; B) = IFT[exp(-B s**2/4) (1 - T(s))] of an atom through the transfer T, and its derivative in ln B,
! tabulated on r (A) and ln B; T is given per Fourier shell of width ds (shell k covers [(k-1) ds, k ds)) and
! interpolated linearly between the shell centres, held beyond them. T = 1 everywhere is the plain Gaussian (no halo)
type :: transfer_kernels
    private
    logical :: l_plain = .true.
    integer :: ntr = 0, nb = 0
    real    :: dr = 0., ds = 0., lnb0 = 0., dlnb = 0.
    real, allocatable :: tsh(:)                    ! T at the shells
    real, allocatable :: htab(:,:), dhtab(:,:)     ! (0:ntr, nb): h and dh/d ln B
    real, allocatable :: reach(:)                  ! (nb): radius (A) beyond which |K_T| < REACH_FRAC of its peak
contains
    procedure :: new => new_kernels
    procedure :: is_plain
    procedure :: halo_at
    procedure :: kernel_at
    procedure :: profile
    procedure :: t_of
    procedure :: kill => kill_kernels
end type transfer_kernels

type :: nano_density
    private
    integer :: ldim(3) = 0, n = 0
    real    :: smpd = 0., d_nn = 0., bstart = 0.
    type(transfer_kernels) :: tk
    real    :: t_ds = TRANSFER_DS                                   ! shell width of t_smooth (1/A)
    real    :: t_change = -1.                                       ! largest difference of the report's curve from the one used
    real    :: gauge_scale = 1., gauge_bfac = 0.                    ! the last measurement's gauge: c and 4 g (A**2)
    integer :: n_gauge = 0                                          ! shells of the gauge band in the last measurement
    integer :: n_held = 0                                           ! shells held by the measurement
    integer :: n_rounds = 0, round_kept = 0                         ! rounds of the transfer taken, the one kept
    logical :: l_rounds_capped = .false.                            ! the rounds reached NROUNDS_MAX still improving
    real,    allocatable :: round_fields(:,:)                       ! (NROUND_FIELDS, n_rounds)
    real,    allocatable :: round_t(:,:)                            ! (shells, n_rounds): the curve each round held
    real,    allocatable :: t_meas(:), t_smooth(:), pow_s(:)        ! the transfer per shell: measured, smoothed; |F_S|**2
    real,    allocatable :: t_report(:)                             ! the report's measurement on the final model, smoothed
    real,    allocatable :: tfilt(:,:,:)                            ! 1 - T at every coefficient of the halo's transform
    real    :: sig2n = 0., s_int = 0., f_corr = 1.
    real    :: ln_s = 0., beta = 0., ln_bp = 0., gamma = 0., tau_b = 0.
    logical :: l_background = .true., l_profiles = .false., l_noise_resid = .false.
    logical :: l_exists = .false.
    real,    allocatable :: map(:,:,:), model(:,:,:), backg(:,:,:)
    real,    allocatable :: smod(:,:,:), halo(:,:,:)              ! the atoms without the transfer, their halo
    real,    allocatable :: xyz(:,:), xyz0(:,:), ints(:), bfacs(:), var_lni(:), var_lnb(:)
    real,    allocatable :: xyz_h(:,:), ints_h(:), bfacs_h(:)     ! the atoms as the halo was last refreshed
    real,    allocatable :: xyz_o(:,:), ints_o(:), bfacs_o(:)     ! the atoms as the atom field was last rendered
    real,    allocatable :: tau_i(:)
    integer, allocatable :: cn(:), bin_lo(:)
contains
    procedure :: new
    procedure :: kill
    procedure :: set_transfer
    procedure :: noise_from_region
    procedure :: noise_from_residual
    procedure :: sweep
    procedure :: measure_transfer
    procedure :: fit_profiles
    procedure :: fit_model
    procedure :: fit_transfer_rounds
    procedure :: get_rounds
    procedure :: resid_over_noise
    procedure :: add_atoms
    procedure :: keep_atoms
    procedure :: get_residual
    procedure :: aperture_intensities
    procedure :: render_sim
    procedure :: render_fit
    procedure :: get_natoms
    procedure :: get_atoms
    procedure :: get_coordination
    procedure :: get_transfer
    procedure :: get_gauge
    procedure :: has_profiles
    procedure :: get_intensity_profile
    procedure :: get_width_profile
    procedure :: get_noise
    procedure, private :: add_model
    procedure, private :: rerender
    procedure, private :: render_narrow
    procedure, private :: refresh_halo
    procedure, private :: set_kernels
    procedure, private :: fit_one
    procedure, private :: residual_correction
    procedure, private :: check_residual
    procedure, private :: gather
    procedure, private :: centroid
    procedure, private :: update_background
    procedure, private :: update_coordination
    procedure, private :: phantom_noise
    procedure, private :: width_prior
    procedure, private :: set_variances
end type nano_density

contains

    ! KERNELS

    ! the unit-integral Gaussian K(r; B) = (4 pi / B)**1.5 exp(-4 pi**2 r**2 / B), r2 in A**2
    elemental real function gauss_kernel( r2, bfac )
        real, intent(in) :: r2, bfac
        gauss_kernel = (4. * PI / bfac)**1.5 * exp(-4. * PI**2 * r2 / bfac)
    end function gauss_kernel

    ! one atom as the reconstruction shows it: K(B) - a K(B + B_h), of integral 1 - a
    elemental real function halo_kernel( r2, bfac, halo_amp, halo_bfac )
        real, intent(in) :: r2, bfac, halo_amp, halo_bfac
        halo_kernel = gauss_kernel(r2, bfac) - halo_amp * gauss_kernel(r2, bfac + halo_bfac)
    end function halo_kernel

    ! add sgn I [K(B) - a K(B + B_h)] of an atom at xyz (voxels) to arr, within RENDER_NSD of the widest term;
    ! the Gaussians factor over the axes, so each is a product of three one-dimensional tables
    subroutine add_atom( arr, smpd, xyz, intens, bfac, halo_amp, halo_bfac, sgn )
        real, intent(inout) :: arr(:,:,:)
        real, intent(in)    :: smpd, xyz(3), intens, bfac, halo_amp, halo_bfac, sgn
        real, allocatable :: d2(:,:), g1(:,:), g2(:,:)
        real    :: c1, e1, c2, e2, rc, rc2, s2
        integer :: n(3), lo(3), hi(3), i, j, k, l, w
        logical :: l_halo
        n      = shape(arr)
        l_halo = halo_amp > 0.
        c1 = sgn * intens * (4. * PI / bfac)**1.5
        e1 = 4. * PI**2 / bfac
        if( l_halo )then
            c2 = sgn * intens * halo_amp * (4. * PI / (bfac + halo_bfac))**1.5
            e2 = 4. * PI**2 / (bfac + halo_bfac)
            rc = RENDER_NSD * sqrt((bfac + halo_bfac) / (8. * PI**2))
        else
            c2 = 0.
            e2 = 0.
            rc = RENDER_NSD * sqrt(bfac / (8. * PI**2))
        endif
        rc2 = rc * rc
        s2  = smpd * smpd
        lo  = max(1, floor(xyz - rc / smpd))
        hi  = min(n, ceiling(xyz + rc / smpd))
        if( any(hi < lo) ) return
        w = maxval(hi - lo) + 1
        allocate(d2(w,3), g1(w,3), g2(w,3), source=0.)
        do l = 1,3
            do i = lo(l),hi(l)
                d2(i-lo(l)+1,l) = (real(i) - xyz(l))**2 * s2
            enddo
            g1(:,l) = exp(-e1 * d2(:,l))
            if( l_halo ) g2(:,l) = exp(-e2 * d2(:,l))
        enddo
        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                    if( d2(i-lo(1)+1,1) + d2(j-lo(2)+1,2) + d2(k-lo(3)+1,3) > rc2 ) cycle
                    arr(i,j,k) = arr(i,j,k) + c1 * g1(i-lo(1)+1,1) * g1(j-lo(2)+1,2) * g1(k-lo(3)+1,3)&
                        &- c2 * g2(i-lo(1)+1,1) * g2(j-lo(2)+1,2) * g2(k-lo(3)+1,3)
                enddo
            enddo
        enddo
    end subroutine add_atom

    ! arr = sum_i I_i [K(B_i) - a K(B_i + B_h)], atoms at xyz(:,i) in voxels
    subroutine render_atoms( arr, smpd, xyz, ints, bfacs, halo_amp, halo_bfac )
        real, intent(out) :: arr(:,:,:)
        real, intent(in)  :: smpd, xyz(:,:), ints(:), bfacs(:), halo_amp, halo_bfac
        integer :: i
        arr = 0.
        do i = 1,size(ints)
            call add_atom(arr, smpd, xyz(:,i), ints(i), bfacs(i), halo_amp, halo_bfac, 1.)
        enddo
    end subroutine render_atoms

    ! THE TRANSFER'S KERNEL TABLE

    ! the table of T (tsh per shell of width ds) for maps of sampling smpd: h and dh/d ln B by the radial transform
    ! h(r) = 4 pi int_0^smax s**2 exp(-B s**2/4) (1 - T(s)) sinc(2 pi s r) ds, smax the outer edge of the last shell
    ! (trapezoids of TAB_DS_INT); then the reach of every width. T = 1 everywhere is the plain Gaussian
    subroutine new_kernels( self, smpd, ds, tsh )
        class(transfer_kernels), intent(inout) :: self
        real,                    intent(in)    :: smpd, ds, tsh(:)
        real(dp), allocatable :: sv(:), w(:,:), dw(:,:), sincv(:)
        integer,  allocatable :: lcut(:)
        real(dp) :: smax, b, peak, kt, hsum, dhsum
        real     :: r
        integer  :: ns, l, j, i
        call self%kill
        if( smpd <= 0. .or. ds <= 0. .or. size(tsh) < 1 ) THROW_HARD('sampling, shell width and shells must be positive; new_kernels')
        self%ds  = ds
        self%tsh = tsh
        self%l_plain = all(tsh == 1.)
        if( self%l_plain ) return
        self%dr   = TAB_DR_VOX * smpd
        self%ntr  = ceiling(TAB_R_MAX / self%dr)
        self%nb   = NTAB_B
        self%lnb0 = log(BFAC_MIN)
        self%dlnb = log(BFAC_MAX / BFAC_MIN) / real(NTAB_B - 1)
        smax = real(size(tsh), dp) * real(ds, dp)
        ns   = ceiling(smax / real(TAB_DS_INT, dp))
        allocate(sv(0:ns), w(0:ns,self%nb), dw(0:ns,self%nb), lcut(self%nb))
        do l = 0,ns
            sv(l) = smax * real(l, dp) / real(ns, dp)
        enddo
        ! the weights of the integral per width: 4 pi s**2 G_B(s) (1 - T(s)) times the trapezoid weight, to where
        ! G_B(s) = exp(-B s**2/4) falls below 1e-8 (lcut)
        do j = 1,self%nb
            b = exp(real(self%lnb0, dp) + real(j - 1, dp) * real(self%dlnb, dp))
            lcut(j) = min(ns, ceiling(sqrt(4._dp * log(1.e8_dp) / b) / smax * real(ns, dp)))
            do l = 0,ns
                w(l,j) = 4._dp * DPI * sv(l)**2 * exp(-b * sv(l)**2 / 4._dp) * (1._dp - t_at(sv(l))) * (smax / real(ns, dp))
                if( l == 0 .or. l == ns ) w(l,j) = 0.5_dp * w(l,j)
                dw(l,j) = -w(l,j) * b * sv(l)**2 / 4._dp
            enddo
        enddo
        allocate(self%htab(0:self%ntr,self%nb), self%dhtab(0:self%ntr,self%nb))
        !$omp parallel do default(shared) private(i,r,l,sincv,j,hsum,dhsum) schedule(static) proc_bind(close)
        do i = 0,self%ntr
            r = real(i) * self%dr
            allocate(sincv(0:ns))
            do l = 0,ns
                if( i == 0 .or. l == 0 )then
                    sincv(l) = 1._dp
                else
                    sincv(l) = sin(2._dp * DPI * sv(l) * real(r, dp)) / (2._dp * DPI * sv(l) * real(r, dp))
                endif
            enddo
            do j = 1,self%nb
                hsum  = sum(w(:lcut(j),j)  * sincv(:lcut(j)))
                dhsum = sum(dw(:lcut(j),j) * sincv(:lcut(j)))
                self%htab(i,j)  = real(hsum)
                self%dhtab(i,j) = real(dhsum)
            enddo
            deallocate(sincv)
        enddo
        !$omp end parallel do
        ! the reach of every width: the largest radius where |K_T| = |K - h| is at least REACH_FRAC of its peak, at most
        ! the table's extent
        allocate(self%reach(self%nb))
        do j = 1,self%nb
            b    = exp(real(self%lnb0, dp) + real(j - 1, dp) * real(self%dlnb, dp))
            peak = 0._dp
            do i = 0,self%ntr
                kt   = real(gauss_kernel((real(i) * self%dr)**2, real(b)), dp) - real(self%htab(i,j), dp)
                peak = max(peak, abs(kt))
            enddo
            self%reach(j) = 0.
            do i = self%ntr,0,-1
                kt = real(gauss_kernel((real(i) * self%dr)**2, real(b)), dp) - real(self%htab(i,j), dp)
                if( abs(kt) >= real(REACH_FRAC, dp) * peak )then
                    self%reach(j) = real(i) * self%dr
                    exit
                endif
            enddo
        enddo

    contains

        ! T at s: linear between the shell centres, held beyond the first and the last
        real(dp) function t_at( sx )
            real(dp), intent(in) :: sx
            real(dp) :: x
            integer  :: k
            x = sx / real(ds, dp) + 0.5_dp
            if( x <= 1._dp )then
                t_at = real(tsh(1), dp)
            else if( x >= real(size(tsh), dp) )then
                t_at = real(tsh(size(tsh)), dp)
            else
                k    = int(x)
                t_at = real(tsh(k), dp) + (x - real(k, dp)) * real(tsh(k+1) - tsh(k), dp)
            endif
        end function t_at

    end subroutine new_kernels

    subroutine kill_kernels( self )
        class(transfer_kernels), intent(inout) :: self
        if( allocated(self%tsh)   ) deallocate(self%tsh)
        if( allocated(self%htab)  ) deallocate(self%htab)
        if( allocated(self%dhtab) ) deallocate(self%dhtab)
        if( allocated(self%reach) ) deallocate(self%reach)
        self%l_plain = .true.
        self%ntr     = 0
        self%nb      = 0
    end subroutine kill_kernels

    logical function is_plain( self )
        class(transfer_kernels), intent(in) :: self
        is_plain = self%l_plain
    end function is_plain

    ! the table's rows at width bfac, interpolated linearly in ln B: h (prof), dh/d ln B (dprof), and the reach (A)
    subroutine profile( self, bfac, prof, dprof, rch )
        class(transfer_kernels), intent(in)  :: self
        real,                    intent(in)  :: bfac
        real,                    intent(out) :: prof(0:), dprof(0:), rch
        real    :: x, wb
        integer :: j
        if( self%l_plain )then
            prof  = 0.
            dprof = 0.
            rch   = 0.
            return
        endif
        x  = min(max((log(bfac) - self%lnb0) / self%dlnb, 0.), real(self%nb - 1) - 1.e-6)
        j  = int(x) + 1
        wb = x - real(j - 1)
        prof  = (1. - wb) * self%htab(:,j)  + wb * self%htab(:,j+1)
        dprof = (1. - wb) * self%dhtab(:,j) + wb * self%dhtab(:,j+1)
        rch   = max(self%reach(j), self%reach(j+1))
    end subroutine profile

    ! T at frequency sx (1/A): linear between the shell centres, held beyond the first and the last
    real function t_of( self, sx )
        class(transfer_kernels), intent(in) :: self
        real,                    intent(in) :: sx
        real    :: x
        integer :: k, nsh
        t_of = 1.
        if( self%l_plain ) return
        nsh = size(self%tsh)
        x = sx / self%ds + 0.5
        if( x <= 1. )then
            t_of = self%tsh(1)
        else if( x >= real(nsh) )then
            t_of = self%tsh(nsh)
        else
            k    = int(x)
            t_of = self%tsh(k) + (x - real(k)) * (self%tsh(k+1) - self%tsh(k))
        endif
    end function t_of

    ! h(r; B) at squared distance r2 (A**2), zero beyond the reach
    real function halo_at( self, r2, bfac )
        class(transfer_kernels), intent(in) :: self
        real,                    intent(in) :: r2, bfac
        real    :: x, wb, r, wr
        integer :: j, i
        halo_at = 0.
        if( self%l_plain ) return
        x  = min(max((log(bfac) - self%lnb0) / self%dlnb, 0.), real(self%nb - 1) - 1.e-6)
        j  = int(x) + 1
        wb = x - real(j - 1)
        r  = sqrt(r2)
        if( r > max(self%reach(j), self%reach(j+1)) ) return
        x  = r / self%dr
        i  = int(x)
        if( i >= self%ntr ) return
        wr = x - real(i)
        halo_at = (1. - wb) * ((1. - wr) * self%htab(i,j)   + wr * self%htab(i+1,j)) +&
                &       wb  * ((1. - wr) * self%htab(i,j+1) + wr * self%htab(i+1,j+1))
    end function halo_at

    ! the atom through the transfer, K_T(r; B) = K(r; B) - h(r; B), at squared distance r2 (A**2)
    real function kernel_at( self, r2, bfac )
        class(transfer_kernels), intent(in) :: self
        real,                    intent(in) :: r2, bfac
        kernel_at = gauss_kernel(r2, bfac) - self%halo_at(r2, bfac)
    end function kernel_at

    ! K_T(r; B) at the squared distances r2 (A**2), the table's rows interpolated once
    function kernel_values( tk, r2, bfac ) result( g )
        class(transfer_kernels), intent(in) :: tk
        real,                    intent(in) :: r2(:), bfac
        real, allocatable :: g(:), prof(:), dprof(:)
        real    :: rch
        integer :: m
        g = gauss_kernel(r2, bfac)
        if( tk%l_plain ) return
        allocate(prof(0:tk%ntr), dprof(0:tk%ntr))
        call tk%profile(bfac, prof, dprof, rch)
        do m = 1,size(r2)
            if( r2(m) > rch * rch ) cycle
            g(m) = g(m) - prof_at(prof, tk%dr, sqrt(r2(m)))
        enddo
    end function kernel_values

    ! value at r (A) of a profile of the table (radial step dr), zero beyond the table
    pure real function prof_at( prof, dr, r )
        real, intent(in) :: prof(0:), dr, r
        real    :: x, wr
        integer :: i
        prof_at = 0.
        x = r / dr
        i = int(x)
        if( i >= ubound(prof, 1) ) return
        wr = x - real(i)
        prof_at = (1. - wr) * prof(i) + wr * prof(i+1)
    end function prof_at

    ! add sgn I h(B) of an atom at xyz (voxels) to arr, within the halo's reach
    subroutine add_halo( arr, smpd, tk, xyz, intens, bfac, sgn )
        real,                   intent(inout) :: arr(:,:,:)
        real,                   intent(in)    :: smpd, xyz(3), intens, bfac, sgn
        class(transfer_kernels), intent(in)   :: tk
        real, allocatable :: prof(:), dprof(:)
        real    :: rch, rv, d2, s2
        integer :: n(3), lo(3), hi(3), i, j, k
        if( tk%l_plain ) return
        allocate(prof(0:tk%ntr), dprof(0:tk%ntr))
        call tk%profile(bfac, prof, dprof, rch)
        prof = sgn * intens * prof
        n  = shape(arr)
        rv = rch / smpd
        lo = max(1, floor(xyz - rv))
        hi = min(n, ceiling(xyz + rv))
        if( any(hi < lo) ) return
        s2 = smpd * smpd
        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                    d2 = sum((real([i,j,k]) - xyz)**2) * s2
                    if( d2 > rch * rch ) cycle
                    arr(i,j,k) = arr(i,j,k) + prof_at(prof, tk%dr, sqrt(d2))
                enddo
            enddo
        enddo
    end subroutine add_halo

    ! PER-ATOM FIT

    ! I and B of I K_T(B) (the atom through the transfer of tk) fitted to the samples y at squared distances r2 (A**2),
    ! noise variance sig2n, with the prior ((ln B - lnb_mean) / lnb_sd)**2 and I >= 0; var_int and var_lnb are the
    ! measurement variances of I and ln B at the solution under white noise (inverse Fisher information of the data alone)
    subroutine fit_atom( y, r2, sig2n, lnb_mean, lnb_sd, tk, intens, bfac, var_int, var_lnb )
        real,                    intent(in)  :: y(:), r2(:), sig2n, lnb_mean, lnb_sd
        class(transfer_kernels), intent(in)  :: tk
        real,                    intent(out) :: intens, bfac, var_int, var_lnb
        real(dp) :: t(NGRID_B), c(NGRID_B), lo, hi, x1, x2, f1, f2, syy, tb, bb
        real(dp) :: g(size(y)), gt(size(y)), k1(size(y)), hh(size(y)), dh(size(y)), rr(size(y)), j11, j12, j22, det
        real, allocatable :: prof(:), dprof(:)
        integer  :: i, ibest
        if( size(y) < 1 .or. size(y) /= size(r2) ) THROW_HARD('samples and distances differ in size; fit_atom')
        if( sig2n <= 0. .or. lnb_sd <= 0. )         THROW_HARD('noise and prior spread must be positive; fit_atom')
        rr  = real(r2, dp)
        syy = sum(real(y, dp)**2)
        allocate(prof(0:max(tk%ntr,0)), dprof(0:max(tk%ntr,0)))
        do i = 1,NGRID_B
            t(i) = log(real(BFAC_MIN, dp)) + real(i-1, dp) * log(real(BFAC_MAX/BFAC_MIN, dp)) / real(NGRID_B-1, dp)
            c(i) = cost(t(i))
        enddo
        ibest = minloc(c, dim=1)
        lo = t(max(ibest-1, 1))
        hi = t(min(ibest+1, NGRID_B))
        x1 = hi - GOLD * (hi - lo)
        x2 = lo + GOLD * (hi - lo)
        f1 = cost(x1)
        f2 = cost(x2)
        do i = 1,NGOLD_B
            if( f1 < f2 )then
                hi = x2
                x2 = x1
                f2 = f1
                x1 = hi - GOLD * (hi - lo)
                f1 = cost(x1)
            else
                lo = x1
                x1 = x2
                f1 = f2
                x2 = lo + GOLD * (hi - lo)
                f2 = cost(x2)
            endif
        enddo
        if( min(f1, f2) > c(ibest) )then
            tb = t(ibest)
        else
            tb = 0.5_dp * (lo + hi)
        endif
        call kernel_on(tb)
        intens = real(int_at())
        bfac = real(exp(tb))
        ! Fisher information of (I, ln B); dK(B)/d ln B = K(B) (4 pi**2 r**2 / B - 1.5), dh/d ln B from the table
        bb = exp(tb)
        gt = k1 * (4._dp * DPI**2 * rr / bb - 1.5_dp) - dh
        j11 = sum(g * g) / sig2n
        j12 = real(intens, dp) * sum(g * gt) / sig2n
        j22 = real(intens, dp)**2 * sum(gt * gt) / sig2n
        det = j11 * j22 - j12 * j12
        if( det > 1.e-12_dp * j11 * j22 .and. det > 0._dp )then
            var_int = real(j22 / det)
            var_lnb = real(j11 / det)
        else
            var_int = huge(1.)
            var_lnb = huge(1.)
        endif

    contains

        ! the kernel K_T at width exp(tt) on the samples (g), with K (k1), h (hh) and dh/d ln B (dh)
        subroutine kernel_on( tt )
            real(dp), intent(in) :: tt
            real(dp) :: b1
            real     :: rch
            integer  :: m
            b1 = exp(tt)
            k1 = (4._dp * DPI / b1)**1.5_dp * exp(-4._dp * DPI**2 * rr / b1)
            hh = 0._dp
            dh = 0._dp
            if( .not. tk%l_plain )then
                call tk%profile(real(b1), prof, dprof, rch)
                do m = 1,size(rr)
                    if( rr(m) > real(rch, dp)**2 ) cycle
                    hh(m) = real(prof_at(prof,  tk%dr, real(sqrt(rr(m)))), dp)
                    dh(m) = real(prof_at(dprof, tk%dr, real(sqrt(rr(m)))), dp)
                enddo
            endif
            g = k1 - hh
        end subroutine kernel_on

        real(dp) function int_at()
            real(dp) :: sgg
            sgg    = sum(g * g)
            int_at = max(sum(real(y, dp) * g) / max(sgg, tiny(1._dp)), 0._dp)
        end function int_at

        real(dp) function cost( tt )
            real(dp), intent(in) :: tt
            real(dp) :: a
            call kernel_on(tt)
            a    = int_at()
            cost = (syy - 2._dp * a * sum(real(y, dp) * g) + a * a * sum(g * g)) / sig2n +&
                &((tt - real(lnb_mean, dp)) / real(lnb_sd, dp))**2
        end function cost

    end subroutine fit_atom

    ! PROFILES

    ! y_i = icept + slope (cn_i - 12) + e_i with e_i ~ N(0, tau(cn_i)**2 + v_i): weighted least squares, the spread
    ! tau in bins of cn (at least NBIN_MIN atoms each, from the interior down; the interior bin is the floor),
    ! tau_all over all atoms. bin_lo(k) is the lowest cn of bin k; bin 1 holds the highest.
    subroutine regress_profile( cn, y, v, slope, icept, bin_lo, tau_bin, tau_all )
        integer,              intent(in)  :: cn(:)
        real,                 intent(in)  :: y(:), v(:)
        real,                 intent(out) :: slope, icept, tau_all
        integer, allocatable, intent(out) :: bin_lo(:)
        real,    allocatable, intent(out) :: tau_bin(:)
        real(dp), allocatable :: w(:), x(:), res(:), tau2(:)
        integer,  allocatable :: ibin(:)
        logical,  allocatable :: finite(:)
        real(dp) :: sw, sx, sy, sxx, sxy, den
        integer  :: m, c, cnt, nb, i, k, ipass
        m = size(y)
        if( m < 2 .or. size(cn) /= m .or. size(v) /= m ) THROW_HARD('profile needs two or more atoms of matching sizes; regress_profile')
        ! bins from the highest coordination down
        allocate(bin_lo(m))
        nb  = 0
        cnt = 0
        do c = maxval(cn),minval(cn),-1
            cnt = cnt + count(cn == c)
            if( cnt >= NBIN_MIN )then
                nb = nb + 1
                bin_lo(nb) = c
                cnt = 0
            endif
        enddo
        if( nb == 0 )then
            nb = 1
            bin_lo(1) = minval(cn)
        else if( cnt > 0 )then
            bin_lo(nb) = minval(cn)
        endif
        bin_lo = bin_lo(:nb)
        allocate(ibin(m))
        do i = 1,m
            do k = 1,nb
                if( cn(i) >= bin_lo(k) )then
                    ibin(i) = k
                    exit
                endif
            enddo
        enddo
        allocate(x(m), w(m), res(m), tau2(nb))
        finite = v < huge(1.)
        x    = real(cn - CN_FULL, dp)
        tau2 = 0._dp
        do ipass = 1,3
            if( ipass == 1 )then
                w = 1._dp
            else
                w = 1._dp / max(real(v, dp) + tau2(ibin), tiny(1._dp))
            endif
            sw  = sum(w)
            sx  = sum(w * x)
            sy  = sum(w * real(y, dp))
            sxx = sum(w * x * x)
            sxy = sum(w * x * real(y, dp))
            den = sw * sxx - sx * sx
            if( den > 1.e-12_dp * sw * sw )then
                slope = real((sw * sxy - sx * sy) / den)
            else
                slope = 0.
            endif
            icept = real((sy - real(slope, dp) * sx) / sw)
            res   = real(y, dp) - real(icept, dp) - real(slope, dp) * x
            do k = 1,nb
                tau2(k) = spread2(ibin == k .and. finite)
            enddo
            tau2 = max(tau2, tau2(1))
        enddo
        tau_bin = real(sqrt(tau2))
        tau_all = real(sqrt(spread2(finite)))

    contains

        ! mean squared residual less mean measurement variance over sel, at least 0
        real(dp) function spread2( sel )
            logical, intent(in) :: sel(:)
            spread2 = 0._dp
            if( count(sel) < 1 ) return
            spread2 = max((sum(res**2, mask=sel) - sum(real(v, dp), mask=sel)) / real(count(sel), dp), 0._dp)
        end function spread2

    end subroutine regress_profile

    ! LIFECYCLE

    ! the model of map (box ldim, sampling smpd) with atoms at xyz (voxels), d_nn the nearest-neighbour distance;
    ! widths start at a Gaussian of BSTART_NN_FRAC d_nn, the transfer at T = 1 (the plain Gaussian), the intensities
    ! at the linear fit of that kernel on the map; l_background (default on) fits the background
    subroutine new( self, map, smpd, d_nn, xyz, l_background )
        class(nano_density), intent(inout) :: self
        real,                intent(in)    :: map(:,:,:), smpd, d_nn, xyz(:,:)
        logical, optional,   intent(in)    :: l_background
        real, allocatable :: y(:), r2(:), g(:)
        integer :: i
        call self%kill
        if( size(xyz,1) /= 3 .or. size(xyz,2) < 1 ) THROW_HARD('positions must be 3 x n, n > 0; new')
        if( smpd <= 0. .or. d_nn <= 0. )             THROW_HARD('sampling and d_NN must be positive; new')
        self%ldim  = shape(map)
        self%smpd  = smpd
        self%d_nn  = d_nn
        self%n     = size(xyz,2)
        self%map   = map
        allocate(self%model(self%ldim(1),self%ldim(2),self%ldim(3)), self%backg(self%ldim(1),self%ldim(2),self%ldim(3)),&
            &self%smod(self%ldim(1),self%ldim(2),self%ldim(3)), self%halo(self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        self%xyz   = xyz
        self%xyz0  = xyz
        self%bstart    = 8. * PI**2 * (BSTART_NN_FRAC * d_nn)**2
        call self%set_kernels(TRANSFER_DS, [1.])
        self%l_background = .true.
        if( present(l_background) ) self%l_background = l_background
        allocate(self%bfacs(self%n), source=self%bstart)
        allocate(self%ints(self%n), self%var_lni(self%n), self%var_lnb(self%n), source=0.)
        allocate(self%cn(self%n), source=0)
        self%xyz_h   = xyz
        self%ints_h  = self%ints
        self%bfacs_h = self%bfacs
        self%xyz_o   = xyz
        self%ints_o  = self%ints
        self%bfacs_o = self%bfacs
        self%l_exists = .true.
        do i = 1,self%n
            call self%gather(self%xyz(:,i), y, r2)
            if( size(y) < 1 ) THROW_HARD('an atom lies outside the map; new')
            g = gauss_kernel(r2, self%bstart)
            self%ints(i) = max(sum(y * g) / sum(g * g), 0.)
        enddo
        call self%rerender
        call self%update_coordination
    end subroutine new

    subroutine kill( self )
        class(nano_density), intent(inout) :: self
        if( .not. self%l_exists ) return
        deallocate(self%map, self%model, self%backg, self%smod, self%halo, self%xyz, self%xyz0, self%ints, self%bfacs,&
            &self%var_lni, self%var_lnb, self%cn, self%xyz_h, self%ints_h, self%bfacs_h, self%xyz_o, self%ints_o,&
            &self%bfacs_o)
        if( allocated(self%tau_i)  ) deallocate(self%tau_i)
        if( allocated(self%bin_lo) ) deallocate(self%bin_lo)
        if( allocated(self%t_meas)   ) deallocate(self%t_meas)
        if( allocated(self%t_smooth) ) deallocate(self%t_smooth)
        if( allocated(self%pow_s)    ) deallocate(self%pow_s)
        if( allocated(self%t_report) ) deallocate(self%t_report)
        if( allocated(self%round_fields) ) deallocate(self%round_fields)
        if( allocated(self%round_t) )      deallocate(self%round_t)
        self%n_rounds        = 0
        self%round_kept      = 0
        self%l_rounds_capped = .false.
        if( allocated(self%tfilt)    ) deallocate(self%tfilt)
        call self%tk%kill
        self%t_change      = -1.
        self%gauge_scale   = 1.
        self%gauge_bfac    = 0.
        self%n_gauge       = 0
        self%n_held        = 0
        self%ldim          = 0
        self%n             = 0
        self%sig2n         = 0.
        self%s_int         = 0.
        self%f_corr        = 1.
        self%ln_s          = 0.
        self%beta          = 0.
        self%ln_bp         = 0.
        self%gamma         = 0.
        self%tau_b         = 0.
        self%l_profiles     = .false.
        self%l_noise_resid  = .false.
        self%l_exists       = .false.
    end subroutine kill

    ! set the transfer to the table t (one value per Fourier shell of width ds from the origin; T = 1 everywhere is the
    ! plain Gaussian); the kernels are tabulated and the model is rendered with it
    subroutine set_transfer( self, ds, t )
        class(nano_density), intent(inout) :: self
        real,                intent(in)    :: ds, t(:)
        if( ds <= 0. .or. size(t) < 1 ) THROW_HARD('the transfer needs a positive shell width and one shell or more; set_transfer')
        call self%set_kernels(ds, t)
        self%t_smooth = t
        self%t_ds     = ds
        call self%rerender
    end subroutine set_transfer

    ! NOISE

    ! sigma_n is the standard deviation of the residual over region (voxels far from every atom); the intensity
    ! noise from phantom sites of the region
    subroutine noise_from_region( self, region )
        class(nano_density), intent(inout) :: self
        logical,             intent(in)    :: region(:,:,:)
        real, allocatable :: vals(:)
        real :: m
        if( any(shape(region) /= self%ldim) ) THROW_HARD('region and map differ in size; noise_from_region')
        if( count(region) < 2 )               THROW_HARD('empty noise region; noise_from_region')
        vals = pack(self%map - self%backg - self%model, region)
        m    = sum(vals) / real(size(vals))
        self%sig2n = sum((vals - m)**2) / real(size(vals) - 1)
        if( self%sig2n <= 0. ) THROW_HARD('no noise in the noise region; noise_from_region')
        self%l_noise_resid = .false.
        call self%phantom_noise(region)
    end subroutine noise_from_region

    ! without an empty region: sigma_n is the robust spread of the fit residual inside envelope farther than
    ! d_NN / 2 from every atom, an upper bound since misfit adds to it; the phantom sites are drawn there
    subroutine noise_from_residual( self, envelope )
        class(nano_density), intent(inout) :: self
        logical,             intent(in)    :: envelope(:,:,:)
        logical, allocatable :: sel(:,:,:)
        real    :: rv, c(3)
        integer :: i, lo(3), hi(3), ii, jj, kk
        if( any(shape(envelope) /= self%ldim) ) THROW_HARD('envelope and map differ in size; noise_from_residual')
        sel = envelope
        rv  = FIT_NN_FRAC * self%d_nn / self%smpd
        do i = 1,self%n
            c  = self%xyz(:,i)
            lo = max(1, floor(c - rv))
            hi = min(self%ldim, ceiling(c + rv))
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        if( sum((real([ii,jj,kk]) - c)**2) <= rv * rv ) sel(ii,jj,kk) = .false.
                    enddo
                enddo
            enddo
        enddo
        if( count(sel) < 2 ) THROW_HARD('no envelope voxel between the atoms; noise_from_residual')
        self%sig2n = robust_spread(pack(self%map - self%backg - self%model, sel))**2
        if( self%sig2n <= 0. ) THROW_HARD('no residual between the atoms; noise_from_residual')
        self%l_noise_resid = .true.
        call self%phantom_noise(sel)
    end subroutine noise_from_residual

    ! the spread s_int of the intensity fitted at fixed width on the residual at phantom sites of sel, and the
    ! factor f_corr by which spatially correlated noise inflates the white-noise variance of that estimate
    subroutine phantom_noise( self, sel )
        class(nano_density), intent(inout) :: self
        logical,             intent(in)    :: sel(:,:,:)
        real, allocatable :: y(:), r2(:), g(:), vals(:)
        real    :: bref, sgg, m
        integer :: nsel, stride, cnt, nsite, ii, jj, kk
        bref = self%bstart
        if( self%l_profiles ) bref = exp(self%ln_bp)
        nsel   = count(sel)
        stride = max(1, nsel / NPHANTOM)
        allocate(vals(nsel / stride + 1))
        cnt   = 0
        nsite = 0
        sgg   = 0.
        do kk = 1,self%ldim(3)
            do jj = 1,self%ldim(2)
                do ii = 1,self%ldim(1)
                    if( .not. sel(ii,jj,kk) ) cycle
                    cnt = cnt + 1
                    if( mod(cnt - 1, stride) /= 0 ) cycle
                    call self%gather(real([ii,jj,kk]), y, r2)
                    g = kernel_values(self%tk, r2, bref)
                    nsite = nsite + 1
                    vals(nsite) = sum(y * g) / sum(g * g)
                    sgg = sgg + sum(g * g)
                enddo
            enddo
        enddo
        if( nsite < 2 ) THROW_HARD('too few phantom sites; phantom_noise')
        m          = sum(vals(:nsite)) / real(nsite)
        self%s_int = sqrt(sum((vals(:nsite) - m)**2) / real(nsite - 1))
        self%f_corr = self%s_int**2 * (sgg / real(nsite)) / self%sig2n
    end subroutine phantom_noise

    ! FIT

    ! one sweep: every atom moves to the centroid of its own residual (at most MOVE_NN_FRAC d_NN from its start) and
    ! is refitted on it within FIT_NN_FRAC d_NN; then the halo is refreshed, the background and the coordination.
    ! The atoms are binned into cubic cells whose side exceeds twice what a fit reads (the fit radius plus the move),
    ! coloured by the parity of the cell indices. The cells of one colour are fitted in parallel against the fields
    ! as they stand: within a cell the atoms go one after the other in index order, each residual corrected
    ! analytically for the atom's own atom and halo and for the changes of the cell's atoms fitted before it (exact
    ! Gauss-Seidel within the cell); cells of one colour are a cell apart, beyond a fit sphere's reach. The atoms are
    ! rendered again before the next colour. Nothing is written while a colour is fitted, so the sweep is the same
    ! whatever the thread count.
    subroutine sweep( self )
        class(nano_density), intent(inout) :: self
        integer, allocatable :: cell_of(:), first(:), next(:), cells(:), upd(:)
        real    :: dmax, rv
        integer :: side, ncell(3), nc, i, ic(3), colour, icell, ncol, j, l, nupd
        if( self%sig2n <= 0. ) THROW_HARD('the noise is not set; sweep')
        dmax  = MOVE_NN_FRAC * self%d_nn / self%smpd
        rv    = FIT_NN_FRAC * self%d_nn / self%smpd
        side  = ceiling(2. * (rv + 2. * dmax)) + 1   ! a fit reads within rv of a position within 2 dmax of the atom's
        ncell = max(1, (self%ldim + side - 1) / side)
        nc    = product(ncell)
        ! the atoms of each cell as linked lists, in atom order
        allocate(cell_of(self%n), first(nc), next(self%n), source=0)
        do i = self%n,1,-1
            ic = min(ncell, max(1, floor((self%xyz(:,i) - 1.) / real(side)) + 1))
            cell_of(i)        = ic(1) + ncell(1) * (ic(2) - 1 + ncell(2) * (ic(3) - 1))
            next(i)           = first(cell_of(i))
            first(cell_of(i)) = i
        enddo
        do colour = 0,7
            ncol = 0
            do icell = 1,nc
                if( first(icell) == 0 ) cycle
                if( cell_colour(icell) /= colour ) cycle
                ncol = ncol + 1
            enddo
            if( ncol == 0 ) cycle
            if( allocated(cells) ) deallocate(cells)
            allocate(cells(ncol))
            l = 0
            do icell = 1,nc
                if( first(icell) == 0 ) cycle
                if( cell_colour(icell) /= colour ) cycle
                l = l + 1
                cells(l) = icell
            enddo
            !$omp parallel do default(shared) private(j,i,upd,nupd) schedule(dynamic) proc_bind(close)
            do j = 1,ncol
                allocate(upd(self%n))
                nupd = 0
                i    = first(cells(j))
                do while( i > 0 )
                    call self%fit_one(i, dmax, upd(:nupd))
                    nupd      = nupd + 1
                    upd(nupd) = i
                    i         = next(i)
                enddo
                deallocate(upd)
            enddo
            !$omp end parallel do
            ! the fields take the fitted atoms: the atoms rendered again, the halo by its transform
            call self%rerender
        enddo
        if( self%l_background ) call self%update_background
        call self%update_coordination

    contains

        integer function cell_colour( icell_ )
            integer, intent(in) :: icell_
            integer :: ic_(3)
            ic_(1) = mod(icell_ - 1, ncell(1)) + 1
            ic_(2) = mod((icell_ - 1) / ncell(1), ncell(2)) + 1
            ic_(3) = (icell_ - 1) / (ncell(1) * ncell(2)) + 1
            cell_colour = mod(ic_(1), 2) + 2 * mod(ic_(2), 2) + 4 * mod(ic_(3), 2)
        end function cell_colour

    end subroutine sweep

    ! atom i against the fields as they stand, its own atom and halo added back and the changes of the atoms upd
    ! (fitted since the atom field was rendered) applied: moved to the centroid of its residual (at most dmax from
    ! its start) and refitted on it; the arrays of atom i change, the fields do not
    subroutine fit_one( self, i, dmax, upd )
        class(nano_density), intent(inout) :: self
        integer,             intent(in)    :: i
        real,                intent(in)    :: dmax
        integer,             intent(in)    :: upd(:)
        real, allocatable :: y(:), r2(:)
        real :: pos(3), d0, lnb_mean, lnb_sd, vi, vb
        pos = self%centroid(self%xyz(:,i), iown=i, upd=upd)
        d0  = sqrt(sum((pos - self%xyz0(:,i))**2))
        if( d0 > dmax ) pos = self%xyz0(:,i) + (pos - self%xyz0(:,i)) * dmax / d0
        call self%gather(pos, y, r2, iown=i, upd=upd)
        if( L_CHECK_RESIDUAL ) call self%check_residual(i, pos, y)
        self%xyz(:,i) = pos
        call self%width_prior(i, lnb_mean, lnb_sd)
        call fit_atom(y, r2, self%sig2n, lnb_mean, lnb_sd, self%tk, self%ints(i), self%bfacs(i), vi, vb)
        call self%set_variances(i, vi, vb)
    end subroutine fit_one

    ! the transfer measured on the map: per Fourier shell of TRANSFER_DS, the projection of the map onto the atoms
    ! without the transfer, T = sum Re(F_map conj(F_S)) / sum |F_S|**2 (the origin left out, the Hermitian half
    ! counted twice off the h = 0 and Nyquist planes), S the atoms as fitted (positions, intensities and widths);
    ! smoothed by a three-shell running median and held at the last reliable value where the atoms' power in a shell,
    ! sum |F_S|**2, is below POWER_FRAC of its peak; then in the gauge: ln T = ln c + g s**2 fitted over the reliable
    ! shells of GAUGE_RES_LO-GAUGE_RES_HI, weighted by the atoms' power, and T divided by c exp(g s**2) (the map fixes
    ! the model only up to that pair). With l_apply the curve is the model's from now on: the atoms take c and -4 g, the
    ! kernels are tabulated and the model rendered. Without, it is the report's measurement of the model as it stands:
    ! kept as t_report beside the curve in use, with their largest difference over the reliable shells (t_change).
    ! The map is the one given, background included; the transforms are single-threaded, the same at any thread count
    subroutine measure_transfer( self, l_apply )
        use simple_image, only: image
        class(nano_density), intent(inout) :: self
        logical,             intent(in)    :: l_apply
        type(image) :: img_m, img_s
        complex,  allocatable :: cm(:,:,:), cs(:,:,:)
        real(dp), allocatable :: num(:), den(:)
        real,     allocatable :: tm(:), tv(:), tsm(:)
        logical,  allocatable :: ok(:)
        integer  :: cshape(3), i, j, k, h, kk, ll, ish, nsh, first, nheld, ngauge
        real     :: sfreq, w, tlast, gscale, gbfac
        call img_m%new(self%ldim, self%smpd, wthreads=.false.)
        call img_s%new(self%ldim, self%smpd, wthreads=.false.)
        call img_m%set_rmat(self%map,  .false.)
        call img_s%set_rmat(self%smod, .false.)
        call img_m%fft()
        call img_s%fft()
        cm = img_m%get_cmat()
        cs = img_s%get_cmat()
        call img_m%kill
        call img_s%kill
        cshape = shape(cm)
        nsh = int(sqrt(sum((0.5 / self%smpd)**2 * [1., 1., 1.])) / TRANSFER_DS) + 1
        allocate(num(nsh), den(nsh), source=0._dp)
        do ll = 1,cshape(3)
            kk = ll - 1
            if( kk > self%ldim(3) / 2 ) kk = kk - self%ldim(3)
            do j = 1,cshape(2)
                k = j - 1
                if( k > self%ldim(2) / 2 ) k = k - self%ldim(2)
                do i = 1,cshape(1)
                    h = i - 1
                    if( h == 0 .and. k == 0 .and. kk == 0 ) cycle
                    sfreq = sqrt((real(h) / real(self%ldim(1)))**2 + (real(k) / real(self%ldim(2)))**2 +&
                        &(real(kk) / real(self%ldim(3)))**2) / self%smpd
                    ish = min(int(sfreq / TRANSFER_DS) + 1, nsh)
                    w = 2.
                    if( h == 0 .or. 2 * h == self%ldim(1) ) w = 1.
                    num(ish) = num(ish) + real(w, dp) * real(cm(i,j,ll) * conjg(cs(i,j,ll)), dp)
                    den(ish) = den(ish) + real(w, dp) * real(abs(cs(i,j,ll))**2, dp)
                enddo
            enddo
        enddo
        allocate(tm(nsh), source=0.)
        ok = den > 0._dp
        where( ok ) tm = real(num / max(den, tiny(1._dp)))
        if( count(ok) < 1 ) THROW_HARD('the atoms carry no power; measure_transfer')
        ! the running median over the shells with power, then the hold where the power is below POWER_FRAC of its peak
        allocate(tv(nsh), tsm(nsh))
        do ish = 1,nsh
            tv(ish) = med3(ish)
        enddo
        ok    = den >= real(POWER_FRAC, dp) * maxval(den)
        first = findloc(ok, .true., dim=1)
        tlast = tv(first)
        nheld = 0
        do ish = 1,nsh
            if( ok(ish) )then
                tlast = tv(ish)
            else
                nheld = nheld + 1
            endif
            tsm(ish) = tlast
        enddo
        call gauge
        if( l_apply )then
            self%t_meas      = tm
            self%t_smooth    = tsm
            self%pow_s       = real(den)
            self%n_held      = nheld
            self%gauge_scale = gscale
            self%gauge_bfac  = gbfac
            self%n_gauge     = ngauge
            self%ints  = self%ints * gscale
            ! the phantom intensity spread is an intensity of the same kernels: it takes the gauge's scale with the atoms
            ! (f_corr, its ratio to the kernel's noise, does not change)
            self%s_int = self%s_int * gscale
            self%bfacs = min(max(self%bfacs - gbfac, BFAC_MIN), BFAC_MAX)
            self%t_ds  = TRANSFER_DS
            call self%set_kernels(TRANSFER_DS, self%t_smooth)
            call self%rerender
            if( self%l_background ) call self%update_background
        else
            self%t_report = tsm
            self%t_change = -1.
            if( allocated(self%t_smooth) )then
                if( size(self%t_smooth) == nsh ) self%t_change = maxval(abs(tsm - self%t_smooth), mask=ok)
            endif
        endif

    contains

        ! ln T = ln c + g s**2 by weighted least squares over the reliable shells whose centres lie in the gauge band (c
        ! alone with one or two, nothing without); the curves (measured and smoothed) divided by c exp(g s**2), the factor
        ! held beyond the last reliable shell and the hold done again
        subroutine gauge
            real(dp) :: sx, sy, sxx, sxy, xk, yk, den_, gg, lnc, wk, sw
            real     :: sc, slast, tl
            sx  = 0._dp
            sy  = 0._dp
            sxx = 0._dp
            sxy = 0._dp
            sw  = 0._dp
            ngauge = 0
            do ish = 1,nsh
                sc = (real(ish) - 0.5) * TRANSFER_DS
                if( .not. ok(ish) ) cycle
                if( sc < 1. / GAUGE_RES_LO .or. sc > 1. / GAUGE_RES_HI ) cycle
                if( .not. tsm(ish) > 0. ) cycle
                xk  = real(sc, dp)**2
                yk  = log(real(tsm(ish), dp))
                wk  = den(ish)
                sw  = sw + wk
                sx  = sx + wk * xk
                sy  = sy + wk * yk
                sxx = sxx + wk * xk * xk
                sxy = sxy + wk * xk * yk
                ngauge = ngauge + 1
            enddo
            gg  = 0._dp
            lnc = 0._dp
            if( ngauge >= 3 )then
                den_ = sw * sxx - sx * sx
                if( den_ > 0._dp ) gg = (sw * sxy - sx * sy) / den_
                lnc = (sy - gg * sx) / sw
            else if( ngauge >= 1 )then
                lnc = sy / sw
            endif
            gscale = real(exp(lnc))
            gbfac  = real(4._dp * gg)
            if( ngauge < 1 ) return
            slast = (real(findloc(ok, .true., dim=1, back=.true.)) - 0.5) * TRANSFER_DS
            do ish = 1,nsh
                sc = min((real(ish) - 0.5) * TRANSFER_DS, slast)
                tsm(ish) = real(real(tsm(ish), dp) / exp(lnc + gg * real(sc, dp)**2))
                tm(ish)  = real(real(tm(ish),  dp) / exp(lnc + gg * real(sc, dp)**2))
            enddo
            tl = tsm(first)
            do ish = 1,nsh
                if( ok(ish) )then
                    tl = tsm(ish)
                else
                    tsm(ish) = tl
                endif
            enddo
        end subroutine gauge

        ! median of the measured values of shells ish-1, ish, ish+1; where the window is incomplete (a neighbour off the
        ! curve or without power) the shell keeps its own value, as a running median keeps its ends
        real function med3( ish_ )
            integer, intent(in) :: ish_
            real    :: v(3), tmp
            integer :: m, l, a, b
            m = 0
            do l = max(1, ish_ - 1),min(nsh, ish_ + 1)
                if( .not. den(l) > 0._dp ) cycle
                m = m + 1
                v(m) = tm(l)
            enddo
            if( m < 3 )then
                med3 = 1.
                if( den(ish_) > 0._dp ) med3 = tm(ish_)
                return
            endif
            do a = 1,m-1
                do b = a+1,m
                    if( v(b) < v(a) )then
                        tmp  = v(a)
                        v(a) = v(b)
                        v(b) = tmp
                    endif
                enddo
            enddo
            med3 = v(2)
        end function med3

    end subroutine measure_transfer

    ! the profiles over the atoms with intensity above STRONG_SNR times the phantom intensity spread: ln I on the
    ! coordination (S, beta and tau_I by bin) and ln B on it (B_p, gamma, tau_B)
    subroutine fit_profiles( self )
        class(nano_density), intent(inout) :: self
        logical, allocatable :: strong(:)
        integer, allocatable :: bins_b(:), cns(:)
        real,    allocatable :: taus_b(:)
        real :: slope_b, tau_all_i
        strong = self%ints > STRONG_SNR * self%s_int
        if( count(strong) < 3 )then
            THROW_WARN('too few strong atoms for the profiles; they are kept')
            return
        endif
        cns = pack(self%cn, strong)
        call regress_profile(cns, log(pack(self%ints, strong)), pack(self%var_lni, strong), self%beta, self%ln_s,&
            &self%bin_lo, self%tau_i, tau_all_i)
        call regress_profile(cns, log(pack(self%bfacs, strong)), pack(self%var_lnb, strong), slope_b, self%ln_bp,&
            &bins_b, taus_b, self%tau_b)
        self%gamma      = -slope_b
        self%l_profiles = .true.
    end subroutine fit_profiles

    ! NSWEEPS sweeps with the profiles after each, the transfer as it is set (held fixed)
    subroutine fit_model( self )
        class(nano_density), intent(inout) :: self
        integer :: isw
        do isw = 1,NSWEEPS
            call self%sweep
            call self%fit_profiles
        enddo
    end subroutine fit_model

    ! the rounds of the measured transfer: each measures the curve from the complete model as it stands, holds it and
    ! fits to convergence (fit_model); they end when the curve changes by CURVE_TOL or less over the reliable shells
    ! from the previous round's, NROUNDS_MAX at most (capped: the limit reached without that). The last round is kept
    ! unless its residual over the noise in envelope is above the best round's by more than ROUND_TOL of it; then the
    ! best is (a round that runs away cannot spoil the map). Then the report's measurement on the round kept
    ! (t_report, t_change). The noise is the one set before the rounds, the same for every round
    subroutine fit_transfer_rounds( self, envelope )
        class(nano_density), intent(inout) :: self
        logical,             intent(in)    :: envelope(:,:,:)
        type(nano_density), allocatable :: best
        real,    allocatable :: fields(:,:), tr(:,:)
        logical, allocatable :: rel(:)
        real    :: rr, rr_best, dt
        integer :: iround, ibest, kept, nr
        logical :: capped
        if( any(shape(envelope) /= self%ldim) ) THROW_HARD('envelope and map differ in size; fit_transfer_rounds')
        allocate(fields(NROUND_FIELDS,NROUNDS_MAX), source=0.)
        rr_best = huge(1.)
        ibest   = 0
        nr      = 0
        capped  = .false.
        do iround = 1,NROUNDS_MAX
            call self%measure_transfer(.true.)
            call self%fit_model
            rr = self%resid_over_noise(envelope)
            nr = iround
            if( iround == 1 ) allocate(tr(size(self%t_smooth),NROUNDS_MAX), source=0.)
            tr(:,iround) = self%t_smooth
            ! the change of the curve from the previous round's, over this round's reliable shells
            dt = -1.
            if( iround > 1 )then
                rel = self%pow_s >= POWER_FRAC * maxval(self%pow_s)
                dt  = maxval(abs(tr(:,iround) - tr(:,iround-1)), mask=rel)
            endif
            fields(:,iround) = [rr, self%ln_s, self%ln_bp, self%beta, self%gamma, dt]
            write(logfhandle,'(A,I0,A,F10.6,A,F8.4,A,F8.4,A,F8.4)') 'transfer round ', iround, ': residual over noise ', rr,&
                &', ln S ', self%ln_s, ', beta ', self%beta, ', curve change ', dt
            if( rr < rr_best )then
                rr_best = rr
                ibest   = iround
                if( allocated(best) ) deallocate(best)
                allocate(best, source=self)
            endif
            if( iround > 1 .and. dt <= CURVE_TOL ) exit
            capped = iround == NROUNDS_MAX
        enddo
        kept = nr
        if( rr > rr_best * (1. + ROUND_TOL) )then
            kept = ibest
            select type( self )
                type is( nano_density )
                    self = best
                class default
                    THROW_HARD('a nano_density is expected; fit_transfer_rounds')
            end select
        endif
        deallocate(best)
        self%n_rounds        = nr
        self%round_kept      = kept
        self%l_rounds_capped = capped
        self%round_fields    = fields(:,:nr)
        self%round_t         = tr(:,:nr)
        call self%measure_transfer(.false.)
    end subroutine fit_transfer_rounds

    ! the rounds taken: per round the residual over the noise, ln S, ln B_p, beta, gamma and the curve's change from the
    ! previous round (-1 for the first) (fields), and the curve it held (curves, per shell of the transfer); the round
    ! kept and whether the rounds reached their limit with the curve still changing
    subroutine get_rounds( self, fields, curves, kept, capped )
        class(nano_density),  intent(in)  :: self
        real,    allocatable, intent(out) :: fields(:,:), curves(:,:)
        integer,              intent(out) :: kept
        logical,              intent(out) :: capped
        if( self%n_rounds > 0 )then
            fields = self%round_fields
            curves = self%round_t
        else
            allocate(fields(NROUND_FIELDS,0), curves(0,0))
        endif
        kept   = self%round_kept
        capped = self%l_rounds_capped
    end subroutine get_rounds

    ! the rms of map - background - model over envelope, in noise standard deviations
    real function resid_over_noise( self, envelope )
        class(nano_density), intent(in) :: self
        logical,             intent(in) :: envelope(:,:,:)
        resid_over_noise = 0.
        if( count(envelope) < 1 .or. .not. self%sig2n > 0. ) return
        resid_over_noise = sqrt(sum((self%map - self%backg - self%model)**2, mask=envelope) / real(count(envelope)) /&
            &self%sig2n)
    end function resid_over_noise

    ! ATOM SET

    ! append atoms at xyz (voxels), each fitted once on its own residual within FIT_NN_FRAC d_NN with the width prior
    ! of the profile (the start width before one); later sweeps move them
    subroutine add_atoms( self, xyz )
        class(nano_density), intent(inout) :: self
        real,                intent(in)    :: xyz(:,:)
        real, allocatable :: y(:), r2(:)
        real    :: lnb_mean, lnb_sd, vi, vb
        integer :: m, n0, i
        m = size(xyz,2)
        if( m < 1 ) return
        if( size(xyz,1) /= 3 ) THROW_HARD('positions must be 3 x m; add_atoms')
        if( self%sig2n <= 0. ) THROW_HARD('the noise is not set; add_atoms')
        n0 = self%n
        self%n       = n0 + m
        self%xyz     = reshape([self%xyz, xyz], [3,self%n])
        self%xyz0    = reshape([self%xyz0, xyz], [3,self%n])
        self%ints    = [self%ints, spread(0., 1, m)]
        self%bfacs   = [self%bfacs, spread(self%bstart, 1, m)]
        self%var_lni = [self%var_lni, spread(huge(1.), 1, m)]
        self%var_lnb = [self%var_lnb, spread(huge(1.), 1, m)]
        self%cn      = [self%cn, spread(0, 1, m)]
        self%xyz_h   = reshape([self%xyz_h, xyz], [3,self%n])
        self%ints_h  = [self%ints_h, spread(0., 1, m)]
        self%bfacs_h = [self%bfacs_h, spread(self%bstart, 1, m)]
        self%xyz_o   = reshape([self%xyz_o, xyz], [3,self%n])
        self%ints_o  = [self%ints_o, spread(0., 1, m)]
        self%bfacs_o = [self%bfacs_o, spread(self%bstart, 1, m)]
        call self%update_coordination
        do i = n0+1,self%n
            call self%gather(self%xyz(:,i), y, r2)
            if( size(y) < 1 ) cycle
            call self%width_prior(i, lnb_mean, lnb_sd)
            call fit_atom(y, r2, self%sig2n, lnb_mean, lnb_sd, self%tk, self%ints(i), self%bfacs(i), vi, vb)
            call self%set_variances(i, vi, vb)
            call self%add_model(i, 1.)
            ! its halo too, rendered from the table (a few atoms), so that the next new atom's residual holds it
            self%ints_h(i)  = self%ints(i)
            self%bfacs_h(i) = self%bfacs(i)
            self%ints_o(i)  = self%ints(i)
            self%bfacs_o(i) = self%bfacs(i)
            if( self%ints(i) > 0. )then
                call add_halo(self%halo,  self%smpd, self%tk, self%xyz(:,i), self%ints(i), self%bfacs(i),  1.)
                call add_halo(self%model, self%smpd, self%tk, self%xyz(:,i), self%ints(i), self%bfacs(i), -1.)
            endif
        enddo
        ! the halos rendered from the table stop at its reach; the transform holds them whole
        call self%rerender
    end subroutine add_atoms

    ! keep the atoms where keep is true, in their order
    subroutine keep_atoms( self, keep )
        class(nano_density), intent(inout) :: self
        logical,             intent(in)    :: keep(:)
        integer :: m
        if( size(keep) /= self%n ) THROW_HARD('mask and atoms differ in size; keep_atoms')
        m = count(keep)
        if( m < 1 ) THROW_HARD('no atom kept; keep_atoms')
        if( m == self%n ) return
        self%xyz     = reshape(pack(self%xyz,  spread(keep, 1, 3)), [3,m])
        self%xyz0    = reshape(pack(self%xyz0, spread(keep, 1, 3)), [3,m])
        self%ints    = pack(self%ints,    keep)
        self%bfacs   = pack(self%bfacs,   keep)
        self%var_lni = pack(self%var_lni, keep)
        self%var_lnb = pack(self%var_lnb, keep)
        self%cn      = pack(self%cn,      keep)
        self%xyz_h   = reshape(pack(self%xyz_h, spread(keep, 1, 3)), [3,m])
        self%ints_h  = pack(self%ints_h,  keep)
        self%bfacs_h = pack(self%bfacs_h, keep)
        self%xyz_o   = reshape(pack(self%xyz_o, spread(keep, 1, 3)), [3,m])
        self%ints_o  = pack(self%ints_o,  keep)
        self%bfacs_o = pack(self%bfacs_o, keep)
        self%n       = m
        call self%rerender
        call self%update_coordination
    end subroutine keep_atoms

    ! arr = map - background - model
    subroutine get_residual( self, arr )
        class(nano_density), intent(in)  :: self
        real,                intent(out) :: arr(:,:,:)
        if( any(shape(arr) /= self%ldim) ) THROW_HARD('array and map differ in size; get_residual')
        arr = self%map - self%backg - self%model
    end subroutine get_residual

    ! the aperture intensity of every atom in src (the map, or a half map on its grid): s**3 times the sum of
    ! src - background - model with the atom's own density added back, within FIT_NN_FRAC d_NN, over the fraction of
    ! its own kernel the aperture holds
    subroutine aperture_intensities( self, src, aper )
        class(nano_density), intent(in)  :: self
        real,                intent(in)  :: src(:,:,:)
        real,                intent(out) :: aper(:)
        real    :: rad, rv, dd, c(3), y, own, frac
        integer :: i, lo(3), hi(3), ii, jj, kk
        if( any(shape(src) /= self%ldim) ) THROW_HARD('map sizes differ; aperture_intensities')
        if( size(aper) /= self%n )         THROW_HARD('one aperture intensity per atom; aperture_intensities')
        rad = FIT_NN_FRAC * self%d_nn
        rv  = rad / self%smpd
        !$omp parallel do default(shared) private(i,c,lo,hi,y,frac,kk,jj,ii,dd,own) schedule(static) proc_bind(close)
        do i = 1,self%n
            c    = self%xyz(:,i)
            lo   = max(1, floor(c - rv))
            hi   = min(self%ldim, ceiling(c + rv))
            y    = 0.
            frac = 0.
            do kk = lo(3),hi(3)
                do jj = lo(2),hi(2)
                    do ii = lo(1),hi(1)
                        dd = sum((real([ii,jj,kk]) - c)**2) * self%smpd**2
                        if( dd > rad * rad ) cycle
                        own  = self%tk%kernel_at(dd, self%bfacs(i))
                        y    = y + src(ii,jj,kk) - self%backg(ii,jj,kk) - self%model(ii,jj,kk) + self%ints(i) * own
                        frac = frac + own
                    enddo
                enddo
            enddo
            aper(i) = 0.
            if( frac > 0. ) aper(i) = y / frac
        enddo
        !$omp end parallel do
    end subroutine aperture_intensities

    ! PRODUCTS AND GETTERS

    ! sum_i I_i K(B_i): the atoms without the transfer and the background
    subroutine render_sim( self, arr )
        class(nano_density), intent(in)  :: self
        real,                intent(out) :: arr(:,:,:)
        if( any(shape(arr) /= self%ldim) ) THROW_HARD('array and map differ in size; render_sim')
        arr = self%smod
    end subroutine render_sim

    ! the full model: the background plus the atoms through the transfer
    subroutine render_fit( self, arr )
        class(nano_density), intent(in)  :: self
        real,                intent(out) :: arr(:,:,:)
        if( any(shape(arr) /= self%ldim) ) THROW_HARD('array and map differ in size; render_fit')
        arr = self%backg + self%model
    end subroutine render_fit

    integer function get_natoms( self )
        class(nano_density), intent(in) :: self
        get_natoms = self%n
    end function get_natoms

    ! positions (voxels), intensities, widths and the measurement variances of ln I and ln B
    subroutine get_atoms( self, xyz, ints, bfacs, var_lni, var_lnb )
        class(nano_density), intent(in)  :: self
        real, allocatable,   intent(out) :: xyz(:,:), ints(:), bfacs(:), var_lni(:), var_lnb(:)
        xyz     = self%xyz
        ints    = self%ints
        bfacs   = self%bfacs
        var_lni = self%var_lni
        var_lnb = self%var_lnb
    end subroutine get_atoms

    function get_coordination( self ) result( cn )
        class(nano_density), intent(in) :: self
        integer, allocatable :: cn(:)
        cn = self%cn
    end function get_coordination

    ! the transfer per Fourier shell of width ds from the origin: as measured (t_meas, with the atoms' power pow_s and
    ! whether it makes the shell reliable; empty when the transfer was set, not measured), as used (t_smooth), and as
    ! the report measured it on the final model (t_report, empty without); the shells the measurement held, and the
    ! largest difference of the report's curve from the one used over the reliable shells (-1 without)
    subroutine get_transfer( self, ds, t_meas, t_smooth, pow_s, reliable, t_report, n_held, t_change )
        class(nano_density),  intent(in)  :: self
        real,                 intent(out) :: ds
        real,    allocatable, intent(out) :: t_meas(:), t_smooth(:), pow_s(:), t_report(:)
        logical, allocatable, intent(out) :: reliable(:)
        integer,              intent(out) :: n_held
        real,                 intent(out) :: t_change
        ds = self%t_ds
        if( allocated(self%t_meas) )then
            t_meas   = self%t_meas
            pow_s    = self%pow_s
            reliable = self%pow_s >= POWER_FRAC * maxval(self%pow_s) .and. self%pow_s > 0.
        else
            allocate(t_meas(0), pow_s(0), reliable(0))
        endif
        if( allocated(self%t_smooth) )then
            t_smooth = self%t_smooth
        else
            allocate(t_smooth(0))
        endif
        if( allocated(self%t_report) )then
            t_report = self%t_report
        else
            allocate(t_report(0))
        endif
        n_held   = self%n_held
        t_change = self%t_change
    end subroutine get_transfer

    ! the gauge of the last measurement: the scale c and the width 4 g (A**2) taken out of the measured transfer and given
    ! to the atoms, and the shells of the gauge band it was fitted on (0: no gauge)
    subroutine get_gauge( self, scale, bfac, nband )
        class(nano_density), intent(in)  :: self
        real,                intent(out) :: scale, bfac
        integer,             intent(out) :: nband
        scale = self%gauge_scale
        bfac  = self%gauge_bfac
        nband = self%n_gauge
    end subroutine get_gauge

    ! whether the profiles were fitted (they need three strong atoms)
    logical function has_profiles( self )
        class(nano_density), intent(in) :: self
        has_profiles = self%l_profiles
    end function has_profiles

    ! ln S, beta, and tau_I by coordination bin (bin_lo the lowest coordination of each bin)
    subroutine get_intensity_profile( self, ln_s, beta, bin_lo, tau_i )
        class(nano_density),  intent(in)  :: self
        real,                 intent(out) :: ln_s, beta
        integer, allocatable, intent(out) :: bin_lo(:)
        real,    allocatable, intent(out) :: tau_i(:)
        if( .not. self%l_profiles ) THROW_HARD('no profile fitted yet; get_intensity_profile')
        ln_s   = self%ln_s
        beta   = self%beta
        bin_lo = self%bin_lo
        tau_i  = self%tau_i
    end subroutine get_intensity_profile

    subroutine get_width_profile( self, ln_bp, gamma, tau_b )
        class(nano_density), intent(in)  :: self
        real,                intent(out) :: ln_bp, gamma, tau_b
        if( .not. self%l_profiles ) THROW_HARD('no profile fitted yet; get_width_profile')
        ln_bp = self%ln_bp
        gamma = self%gamma
        tau_b = self%tau_b
    end subroutine get_width_profile

    ! sigma_n, the phantom intensity spread, the correlation factor, and whether sigma_n is the residual's
    subroutine get_noise( self, sig_n, s_int, f_corr, l_resid )
        class(nano_density), intent(in)  :: self
        real,                intent(out) :: sig_n, s_int, f_corr
        logical,             intent(out) :: l_resid
        sig_n   = sqrt(self%sig2n)
        s_int   = self%s_int
        f_corr  = self%f_corr
        l_resid = self%l_noise_resid
    end subroutine get_noise

    ! PRIVATE

    ! the width prior of atom i: the profile ln B_p + gamma (12 - c_i) with spread tau_B, the start width before one
    subroutine width_prior( self, i, lnb_mean, lnb_sd )
        class(nano_density), intent(in)  :: self
        integer,             intent(in)  :: i
        real,                intent(out) :: lnb_mean, lnb_sd
        if( self%l_profiles )then
            lnb_mean = self%ln_bp + self%gamma * real(CN_FULL - self%cn(i))
            lnb_sd   = max(self%tau_b, TAU_B_MIN)
        else
            lnb_mean = log(self%bstart)
            lnb_sd   = PRIOR_SD_START
        endif
    end subroutine width_prior

    ! the measurement variances of ln I and ln B of atom i from the white-noise variances of the fit, inflated by
    ! f_corr; huge where undefined
    subroutine set_variances( self, i, vi, vb )
        class(nano_density), intent(inout) :: self
        integer,             intent(in)    :: i
        real,                intent(in)    :: vi, vb
        self%var_lni(i) = huge(1.)
        self%var_lnb(i) = huge(1.)
        if( self%ints(i) > sqrt(vi) * 1.e-6 .and. vi < huge(1.) ) self%var_lni(i) = self%f_corr * (vi / self%ints(i)) / self%ints(i)
        if( vb < huge(1.) ) self%var_lnb(i) = self%f_corr * vb
    end subroutine set_variances

    ! atom i's narrow part, sgn I_i K(B_i), into the atoms and the model; its halo lives in the halo field
    subroutine add_model( self, i, sgn )
        class(nano_density), intent(inout) :: self
        integer,             intent(in)    :: i
        real,                intent(in)    :: sgn
        if( self%ints(i) <= 0. ) return
        call add_atom(self%smod,  self%smpd, self%xyz(:,i), self%ints(i), self%bfacs(i), 0., 0., sgn)
        call add_atom(self%model, self%smpd, self%xyz(:,i), self%ints(i), self%bfacs(i), 0., 0., sgn)
    end subroutine add_model

    ! the atoms from scratch, then the halo
    subroutine rerender( self )
        class(nano_density), intent(inout) :: self
        call self%render_narrow
        call self%refresh_halo
    end subroutine rerender

    ! S from every atom, in parallel over planes (every atom adds to every plane it reaches, in atom order), and
    ! model = S - H with the halo as it is
    subroutine render_narrow( self )
        class(nano_density), intent(inout) :: self
        real    :: rc
        integer :: i, k, klo, khi
        self%smod = 0.
        !$omp parallel do default(shared) private(k,i,rc,klo,khi) schedule(static) proc_bind(close)
        do k = 1,self%ldim(3)
            do i = 1,self%n
                if( self%ints(i) <= 0. ) cycle
                rc  = RENDER_NSD * sqrt(self%bfacs(i) / (8. * PI**2)) / self%smpd
                klo = floor(self%xyz(3,i) - rc)
                khi = ceiling(self%xyz(3,i) + rc)
                if( k < klo .or. k > khi ) cycle
                call add_atom(self%smod(:,:,k:k), self%smpd, self%xyz(:,i) - [0., 0., real(k - 1)], self%ints(i),&
                    &self%bfacs(i), 0., 0., 1.)
            enddo
        enddo
        !$omp end parallel do
        self%model   = self%smod - self%halo
        self%xyz_o   = self%xyz
        self%ints_o  = self%ints
        self%bfacs_o = self%bfacs
    end subroutine render_narrow

    ! H = IFT[F_S (1 - T)], model = S - H; the atoms as they are now are what the halo holds until the next refresh
    subroutine refresh_halo( self )
        use simple_image, only: image
        class(nano_density), intent(inout) :: self
        type(image)   :: img
        real, pointer :: p(:,:,:)
        if( self%tk%is_plain() )then
            self%halo = 0.
        else
            call img%new(self%ldim, self%smpd, fft_nthreads=min(nthr_glob, FFT_NTHR_MAX))
            call img%set_rmat(self%smod, .false.)
            call img%fft()
            call img%mul_cmat(self%tfilt)
            call img%ifft()
            call img%get_rmat_ptr(p)
            self%halo = p(:self%ldim(1),:self%ldim(2),:self%ldim(3))
            call img%kill
        endif
        self%model   = self%smod - self%halo
        self%xyz_h   = self%xyz
        self%ints_h  = self%ints
        self%bfacs_h = self%bfacs
    end subroutine refresh_halo

    ! the transfer's kernels and the filter of the halo's transform, 1 - T at every Fourier coefficient of the box
    subroutine set_kernels( self, ds, t )
        use simple_image, only: image
        class(nano_density), intent(inout) :: self
        real,                intent(in)    :: ds, t(:)
        type(image) :: img
        integer     :: cs(3), i, j, l, h, k, m
        call self%tk%new(self%smpd, ds, t)
        if( allocated(self%tfilt) ) deallocate(self%tfilt)
        if( self%tk%is_plain() ) return
        call img%new(self%ldim, self%smpd, wthreads=.false.)
        cs = img%get_array_shape()
        call img%kill
        allocate(self%tfilt(cs(1),cs(2),cs(3)))
        do l = 1,cs(3)
            m = l - 1
            if( m > self%ldim(3) / 2 ) m = m - self%ldim(3)
            do j = 1,cs(2)
                k = j - 1
                if( k > self%ldim(2) / 2 ) k = k - self%ldim(2)
                do i = 1,cs(1)
                    h = i - 1
                    self%tfilt(i,j,l) = 1. - self%tk%t_of(sqrt((real(h) / real(self%ldim(1)))**2 +&
                        &(real(k) / real(self%ldim(2)))**2 + (real(m) / real(self%ldim(3)))**2) / self%smpd)
                enddo
            enddo
        enddo
    end subroutine set_kernels

    ! the correction of the field residual at voxel x for a fit of atom iown: its own atom as the atom field holds
    ! it (the arrays at the last render) added back, its own halo as the halo field holds it (the last refresh)
    ! added back, and the changes of the atoms upd since the atom field was rendered taken out
    real function residual_correction( self, x, iown, upd )
        class(nano_density), intent(in) :: self
        real,                intent(in) :: x(3)
        integer,             intent(in) :: iown, upd(:)
        real    :: dh
        integer :: l, m
        logical :: l_halo
        l_halo = .not. self%tk%is_plain()
        residual_correction = 0.
        if( self%ints_o(iown) > 0. )then
            dh = sum((x - self%xyz_o(:,iown))**2) * self%smpd**2
            residual_correction = residual_correction + self%ints_o(iown) * gauss_kernel(dh, self%bfacs_o(iown))
        endif
        if( l_halo .and. self%ints_h(iown) > 0. )then
            dh = sum((x - self%xyz_h(:,iown))**2) * self%smpd**2
            residual_correction = residual_correction - self%ints_h(iown) * self%tk%halo_at(dh, self%bfacs_h(iown))
        endif
        do l = 1,size(upd)
            m = upd(l)
            if( self%ints(m) > 0. )then
                dh = sum((x - self%xyz(:,m))**2) * self%smpd**2
                residual_correction = residual_correction - self%ints(m) * gauss_kernel(dh, self%bfacs(m))
                if( l_halo ) residual_correction = residual_correction + self%ints(m) * self%tk%halo_at(dh, self%bfacs(m))
            endif
            if( self%ints_o(m) > 0. )then
                dh = sum((x - self%xyz_o(:,m))**2) * self%smpd**2
                residual_correction = residual_correction + self%ints_o(m) * gauss_kernel(dh, self%bfacs_o(m))
            endif
            if( l_halo .and. self%ints_h(m) > 0. )then
                dh = sum((x - self%xyz_h(:,m))**2) * self%smpd**2
                residual_correction = residual_correction - self%ints_h(m) * self%tk%halo_at(dh, self%bfacs_h(m))
            endif
        enddo
    end function residual_correction

    ! diagnostic: the residual gather returned for atom i at c against the map less the background less every other
    ! atom's kernel evaluated directly from the arrays; a mismatch above 1e-3 of the largest sample is reported
    subroutine check_residual( self, i, c, y )
        class(nano_density), intent(in) :: self
        integer,             intent(in) :: i
        real,                intent(in) :: c(3), y(:)
        real    :: rad, rv, dd, yref, ymax, err
        integer :: lo(3), hi(3), ii, jj, kk, m, j
        rad  = FIT_NN_FRAC * self%d_nn
        rv   = rad / self%smpd
        lo   = max(1, floor(c - rv))
        hi   = min(self%ldim, ceiling(c + rv))
        m    = 0
        err  = 0.
        ymax = maxval(abs(y))
        do kk = lo(3),hi(3)
            do jj = lo(2),hi(2)
                do ii = lo(1),hi(1)
                    dd = sum((real([ii,jj,kk]) - c)**2) * self%smpd**2
                    if( dd > rad * rad ) cycle
                    m    = m + 1
                    yref = self%map(ii,jj,kk) - self%backg(ii,jj,kk)
                    do j = 1,self%n
                        if( j == i .or. self%ints(j) <= 0. ) cycle
                        dd   = sum((real([ii,jj,kk]) - self%xyz(:,j))**2) * self%smpd**2
                        yref = yref - self%ints(j) * self%tk%kernel_at(dd, self%bfacs(j))
                    enddo
                    err = max(err, abs(y(m) - yref))
                enddo
            enddo
        enddo
        if( err > 1.e-3 * ymax ) write(*,'(A,I0,A,ES10.3,A,ES10.3)') 'check_residual: atom ', i, ' mismatch ', err, ' of ', ymax
    end subroutine check_residual

    ! the residual map - background - model within FIT_NN_FRAC d_NN of c (voxels), with the squared distances (A**2);
    ! for atom iown corrected by residual_correction (its own atom and halo added back, the changes of upd applied)
    subroutine gather( self, c, y, r2, iown, upd )
        class(nano_density), intent(in)  :: self
        real,                intent(in)  :: c(3)
        real, allocatable,   intent(out) :: y(:), r2(:)
        integer, optional,   intent(in)  :: iown, upd(:)
        real    :: rad, rv, dd
        real    :: ybuf((2 * (ceiling(FIT_NN_FRAC * self%d_nn / self%smpd) + 1) + 1)**3), rbuf(size(ybuf))
        integer :: lo(3), hi(3), ii, jj, kk, m
        if( present(iown) .neqv. present(upd) ) THROW_HARD('iown and upd come together; gather')
        rad = FIT_NN_FRAC * self%d_nn
        rv  = rad / self%smpd
        lo  = max(1, floor(c - rv))
        hi  = min(self%ldim, ceiling(c + rv))
        m   = 0
        do kk = lo(3),hi(3)
            do jj = lo(2),hi(2)
                do ii = lo(1),hi(1)
                    dd = sum((real([ii,jj,kk]) - c)**2) * self%smpd**2
                    if( dd > rad * rad ) cycle
                    m = m + 1
                    rbuf(m) = dd
                    ybuf(m) = self%map(ii,jj,kk) - self%backg(ii,jj,kk) - self%model(ii,jj,kk)
                    if( present(iown) ) ybuf(m) = ybuf(m) + self%residual_correction(real([ii,jj,kk]), iown, upd)
                enddo
            enddo
        enddo
        y  = ybuf(:m)
        r2 = rbuf(:m)
    end subroutine gather

    ! centroid of the positive residual within FIT_NN_FRAC d_NN of c (voxels); for atom iown the residual is
    ! corrected as in gather
    function centroid( self, c, iown, upd ) result( pos )
        class(nano_density), intent(in) :: self
        real,                intent(in) :: c(3)
        integer, optional,   intent(in) :: iown, upd(:)
        real    :: pos(3), rv, w, wsum
        integer :: lo(3), hi(3), ii, jj, kk
        if( present(iown) .neqv. present(upd) ) THROW_HARD('iown and upd come together; centroid')
        rv   = FIT_NN_FRAC * self%d_nn / self%smpd
        lo   = max(1, floor(c - rv))
        hi   = min(self%ldim, ceiling(c + rv))
        pos  = 0.
        wsum = 0.
        do kk = lo(3),hi(3)
            do jj = lo(2),hi(2)
                do ii = lo(1),hi(1)
                    if( sum((real([ii,jj,kk]) - c)**2) > rv * rv ) cycle
                    w = self%map(ii,jj,kk) - self%backg(ii,jj,kk) - self%model(ii,jj,kk)
                    if( present(iown) ) w = w + self%residual_correction(real([ii,jj,kk]), iown, upd)
                    w    = max(w, 0.)
                    pos  = pos + w * real([ii,jj,kk])
                    wsum = wsum + w
                enddo
            enddo
        enddo
        if( wsum > 0. )then
            pos = pos / wsum
        else
            pos = c
        endif
    end function centroid

    ! the background: the residual smoothed with a Gaussian of standard deviation d_NN
    subroutine update_background( self )
        class(nano_density), intent(inout) :: self
        real, allocatable :: resid(:,:,:)
        resid = self%map - self%model
        call gauss_filter3D(resid, self%d_nn / self%smpd, self%backg)
    end subroutine update_background

    ! contact numbers on the current positions: the atoms within CONTACT_NN_FRAC d_NN
    subroutine update_coordination( self )
        class(nano_density), intent(inout) :: self
        real    :: rc2
        integer :: i, j
        rc2 = (CONTACT_NN_FRAC * self%d_nn / self%smpd)**2
        !$omp parallel do default(shared) private(i,j) schedule(static) proc_bind(close)
        do i = 1,self%n
            self%cn(i) = 0
            do j = 1,self%n
                if( j == i ) cycle
                if( sum((self%xyz(:,i) - self%xyz(:,j))**2) <= rc2 ) self%cn(i) = self%cn(i) + 1
            enddo
        enddo
        !$omp end parallel do
    end subroutine update_coordination

    ! FILTERS AND STATISTICS

    ! b = a convolved with the normalised 3D Gaussian of standard deviation sigma_vox voxels, zero outside the box
    subroutine gauss_filter3D( a, sigma_vox, b )
        real, intent(in)  :: a(:,:,:)
        real, intent(in)  :: sigma_vox
        real, intent(out) :: b(:,:,:)
        real, allocatable :: tmp(:,:,:), w(:)
        integer :: n(3), hw, i, j, k, l
        n  = shape(a)
        hw = max(1, ceiling(4. * sigma_vox))
        allocate(w(-hw:hw))
        do l = -hw,hw
            w(l) = exp(-0.5 * (real(l) / sigma_vox)**2)
        enddo
        w = w / sum(w)
        allocate(tmp(n(1),n(2),n(3)))
        !$omp parallel do default(shared) private(i,j,k,l) schedule(static) proc_bind(close)
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    b(i,j,k) = 0.
                    do l = max(-hw, 1-i), min(hw, n(1)-i)
                        b(i,j,k) = b(i,j,k) + w(l) * a(i+l,j,k)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(i,j,k,l) schedule(static) proc_bind(close)
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    tmp(i,j,k) = 0.
                    do l = max(-hw, 1-j), min(hw, n(2)-j)
                        tmp(i,j,k) = tmp(i,j,k) + w(l) * b(i,j+l,k)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(i,j,k,l) schedule(static) proc_bind(close)
        do k = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    b(i,j,k) = 0.
                    do l = max(-hw, 1-k), min(hw, n(3)-k)
                        b(i,j,k) = b(i,j,k) + w(l) * tmp(i,j,k+l)
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine gauss_filter3D

    ! 1.4826 times the median absolute deviation from the median
    real function robust_spread( x )
        use simple_stat, only: median
        real, intent(in) :: x(:)
        real :: med
        med           = median(x)
        robust_spread = 1.4826 * median(abs(x - med))
    end function robust_spread

end module simple_nano_density
