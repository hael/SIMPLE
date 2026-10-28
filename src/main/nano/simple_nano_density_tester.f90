!@descr: unit tests for simple_nano_density: kernels and the transfer table, rendering, per-atom fit, profile regression, noise, the fit of rendered particles with the transfer known and measured, the moved filters
! Expected values are closed forms (kernel integrals, the transfer at zero frequency, the Gaussian convolution identity,
! the notch transfer 1 - a exp(-B_h s**2/4) whose kernel is K(B) - a K(B + B_h)) or the parameters a particle was
! rendered with: a 147-atom cuboctahedral fcc cluster, d_NN 2.75 A at 0.358 A, whose coordination numbers follow from
! its geometry, rendered through the notch, with noise at the real maps' level drawn at fixed seeds.
module simple_nano_density_tester
use simple_test_utils,   only: assert_true, assert_int, assert_real, set_fixed_seed
use simple_defs,         only: dp, PI
use simple_error,        only: simple_exception
use simple_rnd,          only: gasdev
use simple_nano_density, only: nano_density, transfer_kernels, gauss_kernel, halo_kernel, render_atoms, fit_atom,&
    &regress_profile, gauss_filter3D, robust_spread, CN_FULL
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
implicit none
private
public :: run_all_nano_density_tests, rendered_transfer
#include "simple_local_flags.inc"

integer, parameter :: BOX      = 88      ! the particle box, 31.5 A
real,    parameter :: SMPD     = 0.358
real,    parameter :: DNN      = 2.75
integer, parameter :: NSHELL   = 3       ! cuboctahedron of 147 atoms
real,    parameter :: HAMP     = 0.6     ! the transfer the particles are rendered with: 1 - HAMP exp(-HBFAC s**2/4)
real,    parameter :: HBFAC    = 200.
real,    parameter :: DS_FINE  = 0.005   ! shells of the notch table that stands for the known transfer, 1/A
real,    parameter :: DS_MEAS  = 0.02    ! shells of the measured transfer, 1/A (the module's)
! with the transfer measured (ruling 5 of atomic_density_model.md, 2026-10-26): the curve within T_TOL of the rendered
! one where the atoms have power, ln S and ln B_p within MEAS_SCALE_TOL of the truth in the gauge (the scale is a soft
! direction once the curve is measured: the measurement's own noise moves it by about 4% at 147 atoms)
real,    parameter :: T_TOL    = 0.06
real,    parameter :: MEAS_SCALE_TOL = 0.05
! the width of a fully coordinated atom is the real Pt maps' (atomic_density_model.md, Progress of phase 1), whose
! peak over the noise is 15 and whose correlated noise has the smoothing CORR_SD. The fixture is quieter than that
! (PEAK_SNR 40, and half the real maps' intrinsic intensity scatter, INT_SCATTER_IN/OUT) so that the fade beta it
! recovers is determined to about a third of its tolerance across noise draws: the draw depends on the compiler
! (gfortran 15 and 16 seed `random_number` differently), and at the real level the two machines' draws put beta 0.0048
! apart with the kernel known (2026-10-27), one sigma of the 0.005 tolerance. The real level is the 66 maps' business
real,    parameter :: BP       = 28.
real,    parameter :: PEAK_SNR = 40.
real,    parameter :: CORR_SD  = 0.23
real,    parameter :: INT_SCATTER_IN  = 0.015  ! intrinsic scatter of ln I, fully coordinated atoms (real maps 0.03)
real,    parameter :: INT_SCATTER_OUT = 0.03   ! the same at the surface (real maps 0.06)

contains

    subroutine run_all_nano_density_tests()
        write(*,'(A)') '**** running all nano density tests ****'
        call test_kernels()
        call test_transfer_kernels()
        call test_convolution_identity()
        call test_render()
        call test_fit_atom()
        call test_regress_profile()
        call test_noise()
        call test_gauss_filter()
        call test_robust_spread()
        call test_atom_set()
        call test_particle_uniform()
        call test_particle_gradient()
    end subroutine run_all_nano_density_tests

    ! the kernels integrate to 1 and to 1 - a, the transfer at zero frequency (Riemann sums on fine grids)
    subroutine test_kernels()
        real, parameter :: B1 = 20., B2 = 25.
        real(dp) :: s1, s2
        real     :: h, r2
        integer  :: i, j, k, m
        write(*,'(A)') 'test_kernels'
        h  = 0.1
        m  = 30
        s1 = 0._dp
        do k = -m,m
            do j = -m,m
                do i = -m,m
                    s1 = s1 + real(gauss_kernel(real(i*i + j*j + k*k) * h * h, B1), dp)
                enddo
            enddo
        enddo
        call assert_real(1., real(s1) * h**3, 1.e-4, 'gauss_kernel integrates to one')
        h  = 0.25
        m  = 40
        s2 = 0._dp
        do k = -m,m
            do j = -m,m
                do i = -m,m
                    r2 = real(i*i + j*j + k*k) * h * h
                    s2 = s2 + real(halo_kernel(r2, B2, HAMP, HBFAC), dp)
                enddo
            enddo
        enddo
        call assert_real(1. - HAMP, real(s2) * h**3, 1.e-3, 'halo_kernel integrates to 1 - a')
    end subroutine test_kernels

    ! the table of the notch 1 - a exp(-B_h s**2/4) on fine shells reproduces its closed-form kernel K(B) - a K(B + B_h)
    ! and the derivative in ln B of its halo, a B dK(B + B_h)/dB, within 1e-3 of their peaks (the linear interpolation
    ! of the notch between shell centres errs by 2e-4 at most, max |T''| ds**2 / 8); T = 1 is the plain Gaussian exactly
    ! and T = 0 leaves nothing of the atom
    subroutine test_transfer_kernels()
        real, parameter :: BS(3) = [20., 28., 40.]
        type(transfer_kernels) :: tk
        real, allocatable :: t(:), prof(:), dprof(:)
        real    :: r2, err, errd, peak, peakd, rch, bh, dref
        integer :: ib, i
        write(*,'(A)') 'test_transfer_kernels'
        t = notch_table(DS_FINE)
        call tk%new(SMPD, DS_FINE, t)
        call assert_true(.not. tk%is_plain(), 'the notch is not the plain transfer')
        err   = 0.
        errd  = 0.
        allocate(prof(0:1000), dprof(0:1000))
        do ib = 1,3
            peak  = halo_kernel(0., BS(ib), HAMP, HBFAC)
            bh    = BS(ib) + HBFAC
            peakd = 0.
            do i = 0,170
                r2    = (0.05 * real(i))**2
                dref  = HAMP * BS(ib) * gauss_kernel(r2, bh) * (4. * PI**2 * r2 / bh**2 - 1.5 / bh)
                peakd = max(peakd, abs(dref))
            enddo
            do i = 0,170
                r2  = (0.05 * real(i))**2
                err = max(err, abs(tk%kernel_at(r2, BS(ib)) - halo_kernel(r2, BS(ib), HAMP, HBFAC)) / peak)
                dref = HAMP * BS(ib) * gauss_kernel(r2, bh) * (4. * PI**2 * r2 / bh**2 - 1.5 / bh)
                call tk%profile(BS(ib), prof, dprof, rch)
                errd = max(errd, abs(table_value(dprof, sqrt(r2)) - dref) / peakd)
            enddo
        enddo
        write(*,'(A,2ES11.3)') '  notch table: kernel and derivative errors over their peaks ', err, errd
        call assert_true(err  <= 1.e-3, 'the notch table reproduces the closed-form kernel K(B) - a K(B + B_h)')
        call assert_true(errd <= 1.e-3, 'the notch table reproduces the derivative of the halo in ln B')
        call tk%new(SMPD, DS_FINE, spread(1., 1, size(t)))
        call assert_true(tk%is_plain(), 'T = 1 is the plain transfer')
        call assert_real(gauss_kernel(1.3, 28.), tk%kernel_at(1.3, 28.), 0., 'with T = 1 the atom is the Gaussian')
        call tk%new(SMPD, DS_FINE, spread(0., 1, size(t)))
        err = 0.
        do i = 0,170
            r2  = (0.05 * real(i))**2
            err = max(err, abs(tk%kernel_at(r2, 28.)) / gauss_kernel(0., 28.))
        enddo
        call assert_true(err <= 1.e-3, 'with T = 0 nothing of the atom is left')
        call tk%kill

    contains

        real function table_value( pr, r )
            real, intent(in) :: pr(0:), r
            real    :: x, w
            integer :: k
            x = r / (0.1 * SMPD)
            k = int(x)
            w = x - real(k)
            table_value = (1. - w) * pr(k) + w * pr(k+1)
        end function table_value

    end subroutine test_transfer_kernels

    ! K(B1) * K(B2) = K(B1 + B2), with K(B2) applied as fit_transfer does (image%apply_bfac on a zero-padded
    ! cube): the halo atom is K(B1) - a K(B1) * K(B2); render_atoms and atoms%convolve agree on it
    subroutine test_convolution_identity()
        use simple_atoms, only: atoms
        use simple_image, only: image
        integer, parameter :: N = 64
        real,    parameter :: B1 = 25., B2 = 200., Q = 1.3
        real, allocatable :: a1(:,:,:), a2(:,:,:), a3(:,:,:)
        real, pointer     :: p(:,:,:)
        real    :: xyz(3,1)
        type(atoms) :: atm
        type(image) :: vol
        write(*,'(A)') 'test_convolution_identity'
        allocate(a1(N,N,N), a2(N,N,N), a3(N,N,N))
        xyz(:,1) = [32.3, 31.6, 33.2]
        call render_atoms(a1, SMPD, xyz, [Q], [B1], 0., 0.)
        call vol%new([N,N,N], SMPD)
        call vol%set_rmat(a1, .false.)
        call vol%fft()
        call vol%apply_bfac(B2)
        call vol%ifft()
        call vol%get_rmat_ptr(p)
        a2 = p
        call render_atoms(a3, SMPD, xyz, [Q], [B1 + B2], 0., 0.)
        call assert_true(maxval(abs(a2 - a3)) <= 1.e-3 * maxval(a3), 'K(B1) convolved with K(B2) is K(B1 + B2)')
        call render_atoms(a3, SMPD, xyz, [Q], [B1], HAMP, B2)
        call assert_true(maxval(abs(a3 - (a1 - HAMP * a2))) <= 1.e-3 * maxval(a1), 'the halo atom is K(B) - a K(B) * K(B_h)')
        call vol%kill
        call vol%new([N,N,N], SMPD, wthreads=.false.)
        call atm%new(1, dummy=.true.)
        call atm%set_element(1, 'X1')
        call atm%set_coord(1, (xyz(:,1) - 1.) * SMPD)
        call atm%set_occupancy(1, Q)
        call atm%set_beta(1, B1)
        call atm%convolve(vol, 4. * sqrt((B1 + B2) / (8. * PI**2)), halo_amp=HAMP, halo_bfac=B2)
        call assert_true(maxval(abs(vol%get_rmat() - a3)) <= 1.e-4 * maxval(a1), 'render_atoms and atoms%convolve agree on the halo atom')
        call vol%kill
        call atm%kill
    end subroutine test_convolution_identity

    ! rendering adds the atoms; the rendered integral is sum I (1 - a)
    subroutine test_render()
        integer, parameter :: N = 64
        real, allocatable :: arr(:,:,:)
        real    :: xyz(3,2), r2
        write(*,'(A)') 'test_render'
        allocate(arr(N,N,N))
        xyz(:,1) = [30.2, 32.7, 31.1]
        xyz(:,2) = [36.9, 31.4, 33.8]
        call render_atoms(arr, SMPD, xyz, [1.0, 0.7], [25., 32.], HAMP, HBFAC)
        call assert_real(1.7 * (1. - HAMP), sum(arr) * SMPD**3, 0.01 * 1.7 * (1. - HAMP), 'rendered integral is sum I (1 - a)')
        r2 = (sum(([33., 32., 32.] - xyz(:,1))**2) * SMPD**2)
        call assert_real(1.0 * halo_kernel(r2, 25., HAMP, HBFAC) + 0.7 * halo_kernel(sum(([33., 32., 32.] - xyz(:,2))**2) * SMPD**2,&
            &32., HAMP, HBFAC), arr(33,32,32), 1.e-6, 'rendered voxel is the sum of the atoms')
    end subroutine test_render

    ! the per-atom fit returns I and B of a noise-free atom, and in noise is unbiased with the spread it reports
    subroutine test_fit_atom()
        integer, parameter :: NREP = 200
        real,    parameter :: I0 = 1., B0 = 28., SIG = 0.02
        type(transfer_kernels) :: tk
        real, allocatable  :: y(:), r2(:), g(:)
        real    :: fint, fb, vi, vb, ints(NREP), mi, sdi, sdpred
        integer :: i, j, k, m, irep
        write(*,'(A)') 'test_fit_atom'
        call tk%new(SMPD, DS_FINE, notch_table(DS_FINE))
        allocate(r2(729))
        m = 0
        do k = -4,4
            do j = -4,4
                do i = -4,4
                    if( real(i*i + j*j + k*k) * SMPD**2 > (0.5 * DNN)**2 ) cycle
                    m = m + 1
                    r2(m) = real(i*i + j*j + k*k) * SMPD**2
                enddo
            enddo
        enddo
        r2 = r2(:m)
        g  = halo_kernel(r2, B0, HAMP, HBFAC)
        y  = I0 * g
        call fit_atom(y, r2, SIG**2, log(B0), 10., tk, fint, fb, vi, vb)
        call assert_real(I0, fint, 1.e-3 * I0, 'fit_atom: intensity of a noise-free atom')
        call assert_real(B0, fb, 1.e-2 * B0, 'fit_atom: width of a noise-free atom')
        call set_fixed_seed(4711)
        sdpred = 0.
        do irep = 1,NREP
            do i = 1,m
                y(i) = I0 * g(i) + gasdev(0., SIG)
            enddo
            call fit_atom(y, r2, SIG**2, log(B0), 10., tk, ints(irep), fb, vi, vb)
            sdpred = sdpred + vi
        enddo
        sdpred = sqrt(sdpred / real(NREP))
        mi  = sum(ints) / real(NREP)
        sdi = sqrt(sum((ints - mi)**2) / real(NREP - 1))
        ! four standard errors of the mean; the spread within 20%, about four standard errors of an sd from 200
        call assert_real(I0, mi, 4. * sdi / sqrt(real(NREP)), 'fit_atom: intensity unbiased in noise')
        call assert_real(sdi, sdpred, 0.2 * sdi, 'fit_atom: reported spread of the intensity matches its scatter')
    end subroutine test_fit_atom

    ! the regression returns slope and intercept within four standard errors, the interior spread within 25%,
    ! bins of at least 20 atoms from the interior down
    subroutine test_regress_profile()
        integer, parameter :: M = 632
        real,    parameter :: SLOPE0 = 0.05, ICEPT0 = 0.2, TAU_IN = 0.04, TAU_OUT = 0.1
        integer :: cn(M), i, cvals(5) = [5, 7, 8, 9, 12], nper(5) = [100, 200, 100, 12, 220]
        integer, allocatable :: bin_lo(:)
        real,    allocatable :: tau_bin(:)
        real    :: y(M), v(M), w(M), x(M), slope, icept, tau_all, tau, xw, se_slope, se_icept
        integer :: k, j
        write(*,'(A)') 'test_regress_profile'
        call set_fixed_seed(815)
        i = 0
        do k = 1,5
            do j = 1,nper(k)
                i = i + 1
                cn(i) = cvals(k)
            enddo
        enddo
        do i = 1,M
            tau  = TAU_OUT
            if( cn(i) == CN_FULL ) tau = TAU_IN
            v(i) = (0.02 + 0.03 * real(mod(i, 7)) / 6.)**2
            y(i) = ICEPT0 + SLOPE0 * real(cn(i) - CN_FULL) + gasdev(0., sqrt(tau**2 + v(i)))
            w(i) = 1. / (tau**2 + v(i))
        enddo
        call regress_profile(cn, y, v, slope, icept, bin_lo, tau_bin, tau_all)
        x  = real(cn - CN_FULL)
        xw = sum(w * x) / sum(w)
        se_slope = 1. / sqrt(sum(w * (x - xw)**2))
        se_icept = sqrt(1. / sum(w) + xw**2 * se_slope**2)
        call assert_real(SLOPE0, slope, 4. * se_slope, 'regress_profile: slope')
        call assert_real(ICEPT0, icept, 4. * se_icept, 'regress_profile: intercept at full coordination')
        call assert_int(4, size(bin_lo), 'regress_profile: four bins of at least 20 atoms')
        if( size(bin_lo) == 4 ) call assert_true(all(bin_lo == [12, 8, 7, 5]), 'regress_profile: bins 12, 9 with 8, 7 and 5 from the interior down')
        call assert_real(TAU_IN, tau_bin(1), 0.25 * TAU_IN, 'regress_profile: interior spread')
        call assert_real(TAU_OUT, tau_bin(size(tau_bin)), 0.25 * TAU_OUT, 'regress_profile: surface spread')
        cn = CN_FULL
        call regress_profile(cn, y, v, slope, icept, bin_lo, tau_bin, tau_all)
        call assert_real(0., slope, 0., 'regress_profile: no slope without a spread of coordination')
    end subroutine test_regress_profile

    ! on white noise of known sigma: sigma_n within 2%, the phantom intensity spread within 10% of the closed form
    ! sigma / sqrt(sum g**2) and the correlation factor near 1; smoothed noise raises the factor
    subroutine test_noise()
        integer, parameter :: N = 48
        real,    parameter :: SIG = 0.05
        class(nano_density), allocatable :: dens
        real,    allocatable :: map(:,:,:), sm(:,:,:)
        logical, allocatable :: region(:,:,:)
        real    :: xyz(3,1), sig_n, s_int, f_corr, sgg
        logical :: l_resid
        integer :: i, j, k
        write(*,'(A)') 'test_noise'
        allocate(dens)
        call set_fixed_seed(2026)
        allocate(map(N,N,N), sm(N,N,N), region(N,N,N))
        do k = 1,N
            do j = 1,N
                do i = 1,N
                    map(i,j,k) = gasdev(0., SIG)
                enddo
            enddo
        enddo
        region = .false.
        region(10:N-9,10:N-9,10:N-9) = .true.
        xyz(:,1) = [2., 2., 2.]
        call dens%new(map, SMPD, DNN, xyz, l_background=.false.)
        call dens%set_transfer(DS_MEAS, [1.])
        call dens%noise_from_region(region)
        call dens%get_noise(sig_n, s_int, f_corr, l_resid)
        sgg = sum_g2(8. * PI**2 * (0.2 * DNN)**2)
        call assert_real(SIG, sig_n, 0.02 * SIG, 'noise_from_region: sigma of white noise')
        call assert_real(SIG / sqrt(sgg), s_int, 0.1 * SIG / sqrt(sgg), 'phantom intensity spread of white noise')
        call assert_real(1., f_corr, 0.2, 'white noise: correlation factor near 1')
        call assert_true(.not. l_resid, 'noise_from_region: not the residual estimate')
        call gauss_filter3D(map, CORR_SD / SMPD, sm)
        call dens%new(sm, SMPD, DNN, xyz, l_background=.false.)
        call dens%set_transfer(DS_MEAS, [1.])
        call dens%noise_from_region(region)
        call dens%get_noise(sig_n, s_int, f_corr, l_resid)
        call assert_true(f_corr > 2., 'smoothed noise: the correlation factor rises')
        call dens%kill
    end subroutine test_noise

    ! the atom set on two isolated atoms in faint noise: an added atom is fitted to its intensity and width, its aperture
    ! intensity is its intensity, keep_atoms takes an atom out of the model, the residual is the map less the model
    subroutine test_atom_set()
        integer, parameter :: N = 64
        real,    parameter :: I1 = 1.0, I2 = 0.6, B1 = 25., B2 = 32.
        class(nano_density), allocatable :: dens
        real,    allocatable :: map(:,:,:), arr(:,:,:), one(:,:,:), xyz(:,:), ints(:), bfacs(:), vi(:), vb(:), aper(:)
        logical, allocatable :: region(:,:,:)
        real    :: two(3,2), sig
        integer :: i, j, k
        write(*,'(A)') 'test_atom_set'
        two(:,1) = [24.3, 31.6, 32.2]
        two(:,2) = [42.7, 32.4, 31.1]
        allocate(map(N,N,N), arr(N,N,N), one(N,N,N), region(N,N,N))
        call render_atoms(map, SMPD, two, [I1, I2], [B1, B2], HAMP, HBFAC)
        ! noise of 1e-4 of the peak gives the fit its noise variance; the region holds no atom
        call set_fixed_seed(99)
        sig = 1.e-4 * maxval(map)
        do k = 1,N
            do j = 1,N
                do i = 1,N
                    map(i,j,k) = map(i,j,k) + gasdev(0., sig)
                enddo
            enddo
        enddo
        region = .false.
        region(2:10,2:10,2:10) = .true.
        allocate(dens)
        call dens%new(map, SMPD, DNN, two(:,1:1), l_background=.false.)
        call dens%set_transfer(DS_FINE, notch_table(DS_FINE))
        call dens%noise_from_region(region)
        call dens%add_atoms(two(:,2:2))
        call assert_int(2, dens%get_natoms(), 'add_atoms: two atoms')
        call dens%get_atoms(xyz, ints, bfacs, vi, vb)
        call assert_real(I2, ints(2), 2.e-3 * I2, 'add_atoms: the added atom fitted to its intensity')
        call assert_real(B2, bfacs(2), 0.01 * B2, 'add_atoms: the added atom fitted to its width')
        allocate(aper(2))
        call dens%aperture_intensities(map, aper)
        call assert_real(I2, aper(2), 2.e-3 * I2, 'aperture intensity of a fitted isolated atom is its intensity')
        call dens%keep_atoms([.true., .false.])
        call assert_int(1, dens%get_natoms(), 'keep_atoms: one atom kept')
        call dens%get_residual(arr)
        call dens%get_atoms(xyz, ints, bfacs, vi, vb)
        call assert_true(sum((xyz(:,1) - two(:,1))**2) < 1.e-6, 'keep_atoms: the first atom kept')
        ! the model's halo is rendered from the table out to where the atom falls below 1e-4 of its peak, the render is
        ! the closed form: the same to that truncation (an atom's peak is the map's), the interpolation adding less
        call render_atoms(one, SMPD, xyz, ints, bfacs, HAMP, HBFAC)
        call assert_true(maxval(abs(arr - (map - one))) <= 2.e-4 * maxval(map), 'get_residual is the map less the model of the kept atoms')
        call dens%kill
    end subroutine test_atom_set

    ! a unit impulse becomes the normalised separable Gaussian
    subroutine test_gauss_filter()
        integer, parameter :: N = 21, C = 11
        real, parameter    :: SIG = 1.5
        real    :: a(N,N,N), b(N,N,N), w(-6:6)
        integer :: l
        write(*,'(A)') 'test_gauss_filter'
        a = 0.
        a(C,C,C) = 1.
        call gauss_filter3D(a, SIG, b)
        do l = -6,6
            w(l) = exp(-0.5 * (real(l) / SIG)**2)
        enddo
        w = w / sum(w)
        call assert_real(1., sum(b), 1.e-5, 'gauss_filter3D preserves the sum')
        call assert_real(w(2) * w(-1) * w(0), b(C+2,C-1,C), 1.e-7, 'gauss_filter3D is the separable normalised kernel')
    end subroutine test_gauss_filter

    subroutine test_robust_spread()
        real :: x(9)
        write(*,'(A)') 'test_robust_spread'
        x = [1., 2., 3., 4., 5., 6., 7., 8., 100.]
        ! median 5, absolute deviations 4 3 2 1 0 1 2 3 95, their median 2
        call assert_real(1.4826 * 2., robust_spread(x), 1.e-5, 'robust_spread is 1.4826 MAD and ignores the outlier')
    end subroutine test_robust_spread

    ! uniform width, beta = 0.05, white noise: beta within 0.005 and the width rise (12 to 5 neighbours) within 1 A**2
    ! with the transfer known and measured, S within 3% with it known, S and B_p within 5% with it measured, and the
    ! measured transfer within 0.06 of the notch the particle was rendered with; a fit with neighbours subtracted and no
    ! halo or background (the species work's) reads the halos as a width rise of 4.4 A**2 within a factor of two
    subroutine test_particle_uniform()
        real, parameter :: BETA0 = 0.05, GAMMA0 = 0.
        class(nano_density), allocatable :: dens
        real,    allocatable :: map(:,:,:), xyz(:,:), xyz0(:,:), ints(:), bfacs(:), fxyz(:,:), fints(:), fbfacs(:), vi(:), vb(:)
        real,    allocatable :: arr(:,:,:), txyz(:,:), tints(:), tbfacs(:), tvi(:), tvb(:)
        logical, allocatable :: region(:,:,:), envelope(:,:,:)
        integer, allocatable :: cn(:)
        real    :: sig, rms0, rms, sig_n, s_int, f_corr
        integer :: nthr_saved
        logical :: l_resid
        write(*,'(A)') 'test_particle_uniform'
        allocate(dens)
        call make_particle(BETA0, GAMMA0, .false., 31, map, xyz, xyz0, cn, ints, bfacs, region, envelope, sig)
        ! the profiles with the transfer known: what the per-atom fit guarantees
        call dens%new(map, SMPD, DNN, xyz0)
        call dens%set_transfer(DS_FINE, notch_table(DS_FINE))
        call dens%noise_from_region(region)
        call dens%fit_model
        call dens%fit_model
        call check_profiles(dens, 'uniform, transfer known', BETA0, GAMMA0, 0.005)
        ! the transfer measured as well
        call fit_measured(dens, map, xyz0, region, envelope)
        call report_rounds(dens, 'uniform', xyz, ints, bfacs)
        call check_transfer(dens, 'uniform', xyz, ints, bfacs)
        call check_profiles(dens, 'uniform, transfer measured', BETA0, GAMMA0, 0.005, MEAS_SCALE_TOL)
        call dens%get_atoms(fxyz, fints, fbfacs, vi, vb)
        ! the sweep is parallel over cells of a colouring and the transfer's transforms are single-threaded, so a
        ! second run on the same threads is the same fit bit for bit (a race would show), and one thread gives the
        ! same fit up to the order of the threaded sums
        call fit_measured(dens, map, xyz0, region, envelope)
        call dens%get_atoms(txyz, tints, tbfacs, tvi, tvb)
        call assert_true(all(txyz == fxyz) .and. all(tints == fints) .and. all(tbfacs == fbfacs),&
            &'uniform: a second fit on the same threads is the same, bit for bit')
        nthr_saved = 1
        !$ nthr_saved = omp_get_max_threads()
        if( nthr_saved > 1 )then
            !$ call omp_set_num_threads(1)
            call fit_measured(dens, map, xyz0, region, envelope)
            call dens%get_atoms(txyz, tints, tbfacs, tvi, tvb)
            !$ call omp_set_num_threads(nthr_saved)
            call assert_true(maxval(abs(txyz - fxyz)) < 1.e-3 .and. maxval(abs(tints - fints) / fints) < 1.e-3 .and.&
                &maxval(abs(tbfacs - fbfacs) / fbfacs) < 1.e-3, 'uniform: the fit on one thread is the fit on all threads')
        endif
        rms0 = sqrt(sum((xyz0 - xyz)**2) / real(size(cn))) * SMPD
        rms  = sqrt(sum((fxyz - xyz)**2) / real(size(cn))) * SMPD
        call assert_true(rms < 0.05 .and. rms < rms0 / 3., 'uniform: positions within 0.05 A and a third of the start error')
        call assert_true(all(dens%get_coordination() == cn), 'uniform: the coordination of the geometry')
        allocate(arr(BOX,BOX,BOX))
        call dens%render_sim(arr)
        call assert_real(sum(fints), sum(arr) * SMPD**3, 0.01 * sum(fints), 'uniform: the rendered atoms integrate to the sum of the intensities')
        call dens%render_fit(arr)
        call assert_real(sig, sqrt(sum((map - arr)**2, mask=envelope) / real(count(envelope))), 0.15 * sig,&
            &'uniform: the residual of the full model inside the particle is the noise')
        call dens%noise_from_residual(envelope)
        call dens%get_noise(sig_n, s_int, f_corr, l_resid)
        call assert_true(l_resid, 'noise_from_residual: the residual estimate')
        call assert_real(sig, sig_n, 0.15 * sig, 'noise_from_residual: sigma_n from the residual between the atoms')
        ! the same map without the halo term and the background
        call dens%new(map, SMPD, DNN, xyz0, l_background=.false.)
        call dens%set_transfer(DS_MEAS, [1.])
        call dens%noise_from_region(region)
        call dens%fit_model
        call dens%fit_model
        call dens%get_atoms(fxyz, fints, fbfacs, vi, vb)
        call check_artefact(fxyz, fbfacs)
        call dens%kill
    end subroutine test_particle_uniform

    ! width rising towards the surface (gamma = 0.03), intensity falling with coordination (beta = -0.03),
    ! correlated noise: the same checks, beta within the fixture's 0.01
    subroutine test_particle_gradient()
        real, parameter :: BETA0 = -0.03, GAMMA0 = 0.03
        class(nano_density), allocatable :: dens
        real,    allocatable :: map(:,:,:), xyz(:,:), xyz0(:,:), ints(:), bfacs(:)
        logical, allocatable :: region(:,:,:), envelope(:,:,:)
        integer, allocatable :: cn(:)
        real    :: sig
        write(*,'(A)') 'test_particle_gradient'
        allocate(dens)
        call make_particle(BETA0, GAMMA0, .true., 47, map, xyz, xyz0, cn, ints, bfacs, region, envelope, sig)
        call dens%new(map, SMPD, DNN, xyz0)
        call dens%set_transfer(DS_FINE, notch_table(DS_FINE))
        call dens%noise_from_region(region)
        call dens%fit_model
        call dens%fit_model
        call check_profiles(dens, 'gradient, transfer known', BETA0, GAMMA0, 0.01)
        call fit_measured(dens, map, xyz0, region, envelope)
        call report_rounds(dens, 'gradient', xyz, ints, bfacs)
        call check_transfer(dens, 'gradient', xyz, ints, bfacs)
        call check_profiles(dens, 'gradient, transfer measured', BETA0, GAMMA0, 0.01, MEAS_SCALE_TOL)
        call dens%kill
    end subroutine test_particle_gradient

    ! the floors of the fit through the halo against the generating parameters
    ! the profiles against the generating ones; with the transfer known these are what the per-atom fit guarantees.
    ! beta_tol is the fixture's: 0.005 on the white-noise particle and 0.01 on the gradient particle, whose correlated
    ! noise and twelve five-coordinated atoms carry the slope. At the real maps' noise and intrinsic scatter the two
    ! compilers' draws (2026-10-27, kernel known) gave beta 0.0512 and 0.0464 for 0.05, and -0.0288 and -0.0363 for
    ! -0.03: one sigma of the tolerances. The fixture is now 2.7 times quieter and its intrinsic scatter halved, which
    ! by the scaling of both sources puts a draw's scatter near a third of the tolerance; the two machines' builds are
    ! the check of that
    ! scale_tol, with the transfer measured: ln S and ln B_p within it of the truth in the transfer's gauge (ruling 5)
    subroutine check_profiles( dens, tag, beta0, gamma0, beta_tol, scale_tol )
        class(nano_density), intent(in) :: dens
        character(len=*),   intent(in) :: tag
        real,               intent(in) :: beta0, gamma0, beta_tol
        real, optional,     intent(in) :: scale_tol
        integer, allocatable :: bin_lo(:)
        real,    allocatable :: tau_i(:)
        real :: ln_s, beta, ln_bp, gamma, tau_b
        call dens%get_intensity_profile(ln_s, beta, bin_lo, tau_i)
        call dens%get_width_profile(ln_bp, gamma, tau_b)
        write(*,'(A,A,3F9.4,2F9.2)') tag, ': beta, gamma, ln S, B_p, rise ', beta, gamma, ln_s, exp(ln_bp),&
            &width_rise(ln_bp, gamma)
        call assert_real(beta0, beta, beta_tol, tag//': beta within the fixture''s tolerance per neighbour')
        call assert_real(width_rise(log(BP), gamma0), width_rise(ln_bp, gamma), 1., tag//': width rise from 12 to 5 neighbours within 1 A**2')
        if( present(scale_tol) )then
            call assert_real(0., ln_s, scale_tol, tag//': ln S within the tolerance of the measured transfer''s gauge')
            call assert_real(log(BP), ln_bp, scale_tol, tag//': ln B_p within the tolerance of the measured transfer''s gauge')
        else
            call assert_real(0., ln_s, 0.03, tag//': intensity scale S within 3%')
        endif
    end subroutine check_profiles

    ! the production path with the transfer measured: a fit through the plain Gaussian, then the rounds of the transfer
    ! (measured on the complete model, held through a fit, until the curve changes by 0.01 or less), the last round kept
    ! unless clearly worse by the residual over the noise in envelope, and measured again for the report
    subroutine fit_measured( dens, map, xyz0, region, envelope )
        class(nano_density), intent(inout) :: dens
        real,                intent(in)    :: map(:,:,:), xyz0(:,:)
        logical,             intent(in)    :: region(:,:,:), envelope(:,:,:)
        call dens%new(map, SMPD, DNN, xyz0)
        call dens%noise_from_region(region)
        call dens%fit_model
        call dens%fit_transfer_rounds(envelope)
    end subroutine fit_measured

    ! the trajectory of the rounds against the truth: per round ln S, ln B_p, beta, gamma, the largest error of the curve
    ! it held where the atoms have power, the curve's change from the previous round and the residual over the noise;
    ! the limit of rounds is a safety limit, reaching it is reported, not an error (seventh ruling); the round kept is
    ! the last one unless its residual is above the best round's by more than 1e-3 of it
    subroutine report_rounds( dens, tag, xyz, ints, bfacs )
        class(nano_density), intent(in) :: dens
        character(len=*),   intent(in) :: tag
        real,               intent(in) :: xyz(:,:), ints(:), bfacs(:)
        real,    allocatable :: fields(:,:), curves(:,:), t_ref(:)
        logical, allocatable :: have(:)
        real    :: err
        integer :: kept, k, ir, nr
        logical :: capped
        call dens%get_rounds(fields, curves, kept, capped)
        nr = size(fields, 2)
        call assert_true(nr >= 1 .and. kept >= 1 .and. kept <= nr, tag//': the rounds ran and one was kept')
        if( nr < 1 ) return
        call rendered_transfer(xyz, ints, bfacs, BOX, SMPD, HAMP, HBFAC, DS_MEAS, size(curves, 1), t_ref, have)
        do ir = 1,nr
            err = 0.
            do k = 1,size(curves, 1)
                if( have(k) ) err = max(err, abs(curves(k,ir) - t_ref(k)))
            enddo
            write(*,'(A,A,I0,A,F10.6,A,F8.4,A,F8.2,A,F8.4,A,F8.4,A,F8.4,A,F8.4)') tag, ': round ', ir,&
                &': residual over noise ', fields(1,ir), ', ln S ', fields(2,ir), ', B_p ', exp(fields(3,ir)), ', beta ',&
                &fields(4,ir), ', gamma ', fields(5,ir), ', curve error ', err, ', curve change ', fields(6,ir)
        enddo
        write(*,'(A,A,I0,A,I0,A,L1)') tag, ': rounds ', nr, ', kept ', kept, ', capped ', capped
        if( kept /= nr )then
            call assert_true(fields(1,nr) > (1. + 1.e-3) * minval(fields(1,:)) .and. fields(1,kept) == minval(fields(1,:)),&
                &tag//': an earlier round is kept only when the last is clearly worse, and then the best')
        else
            call assert_true(fields(1,nr) <= (1. + 1.e-3) * minval(fields(1,:)),&
                &tag//': the last round is kept when it is not clearly worse than the best')
        endif
    end subroutine report_rounds

    ! the transfer measured on the particle against the one it was rendered with, where its atoms have power (their
    ! power in a shell at least 1% of the peak, the module's criterion, here on the true atoms): the expected value
    ! per shell is the notch averaged over the shell with the true atoms' power as weights, which is what a perfect
    ! measurement returns; within T_TOL
    subroutine check_transfer( dens, tag, xyz, ints, bfacs )
        class(nano_density), intent(in) :: dens
        character(len=*),   intent(in) :: tag
        real,               intent(in) :: xyz(:,:), ints(:), bfacs(:)
        real,    allocatable :: t_meas(:), t_sm(:), pow(:), t_ref(:), t_rep(:)
        logical, allocatable :: rel(:), have(:)
        real    :: ds, t_change, err
        integer :: n_held, k, nchk
        call dens%get_transfer(ds, t_meas, t_sm, pow, rel, t_rep, n_held, t_change)
        call rendered_transfer(xyz, ints, bfacs, BOX, SMPD, HAMP, HBFAC, ds, size(t_sm), t_ref, have)
        err  = 0.
        nchk = 0
        do k = 1,size(t_sm)
            if( .not. have(k) ) cycle
            nchk = nchk + 1
            err  = max(err, abs(t_sm(k) - t_ref(k)))
            write(*,'(A,A,F7.3,A,F8.4,A,F8.4)') tag, ': shell ', (real(k) - 0.5) * ds, ' 1/A, T ', t_sm(k), ', rendered ', t_ref(k)
        enddo
        write(*,'(A,A,I0,A,F8.4,A,F8.4)') tag, ': shells with power ', nchk,&
            &', largest error of the measured transfer ', err, ', change on the final model ', t_change
        call assert_true(nchk >= 10, tag//': the atoms have power in ten shells or more')
        call assert_true(err <= T_TOL, tag//': the measured transfer is the rendered one within 0.06 where the atoms have power')
    end subroutine check_transfer

    ! the transfer a particle was rendered with, as a perfect measurement returns it: the atoms (xyz in voxels, ints,
    ! bfacs) rendered without a transfer in a box of nbox at smpd, then per Fourier shell of width ds (nsh shells from
    ! the origin, the origin left out) the notch 1 - amp exp(-bh s**2/4) averaged over the shell's coefficients with
    ! their power |F|**2 as weights (t_ref); have marks the shells holding at least 1% of the peak power
    subroutine rendered_transfer( xyz, ints, bfacs, nbox, smpd_r, amp, bh, ds, nsh, t_ref, have )
        use simple_image, only: image
        real,                 intent(in)  :: xyz(:,:), ints(:), bfacs(:), smpd_r, amp, bh, ds
        integer,              intent(in)  :: nbox, nsh
        real,    allocatable, intent(out) :: t_ref(:)
        logical, allocatable, intent(out) :: have(:)
        type(image) :: img
        real,     allocatable :: arr(:,:,:)
        complex,  allocatable :: c(:,:,:)
        real(dp), allocatable :: pw(:), pwt(:)
        integer :: cs(3), i, j, l, h, k, m, ish
        real    :: sf, w
        allocate(arr(nbox,nbox,nbox))
        call render_atoms(arr, smpd_r, xyz, ints, bfacs, 0., 0.)
        call img%new([nbox,nbox,nbox], smpd_r, wthreads=.false.)
        call img%set_rmat(arr, .false.)
        call img%fft()
        c  = img%get_cmat()
        cs = shape(c)
        call img%kill
        allocate(pw(nsh), pwt(nsh), source=0._dp)
        do l = 1,cs(3)
            m = l - 1
            if( m > nbox / 2 ) m = m - nbox
            do j = 1,cs(2)
                k = j - 1
                if( k > nbox / 2 ) k = k - nbox
                do i = 1,cs(1)
                    h = i - 1
                    if( h == 0 .and. k == 0 .and. m == 0 ) cycle
                    sf  = sqrt(real(h**2 + k**2 + m**2)) / (real(nbox) * smpd_r)
                    ish = min(int(sf / ds) + 1, nsh)
                    w   = 2.
                    if( h == 0 .or. 2 * h == nbox ) w = 1.
                    pw(ish)  = pw(ish)  + real(w * abs(c(i,j,l))**2, dp)
                    pwt(ish) = pwt(ish) + real(w * abs(c(i,j,l))**2 * (1. - amp * exp(-bh * sf**2 / 4.)), dp)
                enddo
            enddo
        enddo
        have  = pw >= 0.01_dp * maxval(pw) .and. pw > 0._dp
        t_ref = real(pwt / max(pw, tiny(1._dp)))
    end subroutine rendered_transfer

    ! the notch 1 - HAMP exp(-HBFAC s**2/4) at the centres of shells of width ds out to the corner of the box
    function notch_table( ds ) result( t )
        real, intent(in)  :: ds
        real, allocatable :: t(:)
        integer :: k, nsh
        nsh = int(sqrt(3.) * 0.5 / SMPD / ds) + 1
        allocate(t(nsh))
        do k = 1,nsh
            t(k) = 1. - HAMP * exp(-HBFAC * ((real(k) - 0.5) * ds)**2 / 4.)
        enddo
    end function notch_table

    ! B rise of the fit without the halo from the inner to the outer radial fifth: 4.4 A**2 within a factor of two
    subroutine check_artefact( fxyz, fbfacs )
        real, intent(in) :: fxyz(:,:), fbfacs(:)
        real    :: rad(size(fbfacs)), cen(3), rise
        integer :: order(size(fbfacs)), n, nf, i, j, tmp
        n   = size(fbfacs)
        cen = sum(fxyz, dim=2) / real(n)
        do i = 1,n
            rad(i)   = sqrt(sum((fxyz(:,i) - cen)**2))
            order(i) = i
        enddo
        do i = 2,n
            j = i
            do while( j > 1 )
                if( rad(order(j-1)) <= rad(order(j)) ) exit
                tmp = order(j)
                order(j) = order(j-1)
                order(j-1) = tmp
                j = j - 1
            enddo
        enddo
        nf   = n / 5
        rise = sum(fbfacs(order(n-nf+1:n))) / real(nf) - sum(fbfacs(order(1:nf))) / real(nf)
        write(*,'(A,F8.3)') 'no halo term: width rise inner to outer fifth (A**2) ', rise
        call assert_true(rise >= 2.2 .and. rise <= 8.8, 'without the halo term the fit reads a width rise of 4.4 A**2 within a factor of two')
    end subroutine check_artefact

    ! B_p (exp(7 gamma) - 1): the width rise from 12 to 5 neighbours
    real function width_rise( ln_bp, gamma )
        real, intent(in) :: ln_bp, gamma
        width_rise = exp(ln_bp) * (exp(real(CN_FULL - 5) * gamma) - 1.)
    end function width_rise

    ! sum of g**2 over the fit sphere (d_NN / 2) for the plain Gaussian of width b, on the grid
    real function sum_g2( b )
        real, intent(in) :: b
        integer :: i, j, k
        real    :: r2
        sum_g2 = 0.
        do k = -5,5
            do j = -5,5
                do i = -5,5
                    r2 = real(i*i + j*j + k*k) * SMPD**2
                    if( r2 > (0.5 * DNN)**2 ) cycle
                    sum_g2 = sum_g2 + gauss_kernel(r2, b)**2
                enddo
            enddo
        enddo
    end function sum_g2

    ! the cuboctahedral cluster rendered through the transfer with ln I = beta (c - 12) + N(0, tau_I(c)) and
    ! ln B = ln B_p + gamma (12 - c) + N(0, 0.02), plus white or correlated noise; xyz the truth, xyz0 the start
    ! (truth plus 0.1 A per axis), cn from the geometry, region the noise region, envelope the particle's sphere
    subroutine make_particle( beta0, gamma0, l_corr, seed, map, xyz, xyz0, cn, ints, bfacs, region, envelope, sig )
        real,                 intent(in)  :: beta0, gamma0
        logical,              intent(in)  :: l_corr
        integer,              intent(in)  :: seed
        real,    allocatable, intent(out) :: map(:,:,:), xyz(:,:), xyz0(:,:), ints(:), bfacs(:)
        integer, allocatable, intent(out) :: cn(:)
        logical, allocatable, intent(out) :: region(:,:,:), envelope(:,:,:)
        real,                 intent(out) :: sig
        real,    allocatable :: noise(:,:,:), sm(:,:,:)
        real    :: tmp(3,1000), cen(3), step, r, peak, tau, m, rv
        integer :: i, j, k, l, n, cnt(CN_FULL), lo(3), hi(3)
        step = DNN / sqrt(2.)
        cen  = 0.5 * real(BOX) + [0.31, -0.27, 0.13]
        n    = 0
        do k = -NSHELL,NSHELL
            do j = -NSHELL,NSHELL
                do i = -NSHELL,NSHELL
                    if( modulo(i + j + k, 2) /= 0 ) cycle
                    if( abs(i) + abs(j) + abs(k) > 2 * NSHELL ) cycle
                    n = n + 1
                    tmp(:,n) = cen + real([i,j,k]) * step / SMPD
                enddo
            enddo
        enddo
        if( n /= 147 ) THROW_HARD('broken cuboctahedron fixture; make_particle')
        xyz = tmp(:,:n)
        allocate(cn(n), source=0)
        do i = 1,n
            do j = 1,n
                if( j == i ) cycle
                if( sqrt(sum((xyz(:,i) - xyz(:,j))**2)) * SMPD < 1.1 * DNN ) cn(i) = cn(i) + 1
            enddo
        enddo
        cnt = 0
        do i = 1,n
            cnt(cn(i)) = cnt(cn(i)) + 1
        enddo
        if( cnt(5) /= 12 .or. cnt(7) /= 48 .or. cnt(8) /= 24 .or. cnt(9) /= 8 .or. cnt(12) /= 55 )then
            THROW_HARD('cuboctahedron coordination is not 12 x 5, 48 x 7, 24 x 8, 8 x 9, 55 x 12; make_particle')
        endif
        call set_fixed_seed(seed)
        allocate(ints(n), bfacs(n), xyz0(3,n))
        do i = 1,n
            tau = INT_SCATTER_OUT
            if( cn(i) == CN_FULL ) tau = INT_SCATTER_IN
            ints(i)  = exp(beta0 * real(cn(i) - CN_FULL) + gasdev(0., tau))
            bfacs(i) = exp(log(BP) + gamma0 * real(CN_FULL - cn(i)) + gasdev(0., 0.02))
            do j = 1,3
                xyz0(j,i) = xyz(j,i) + gasdev(0., 0.1) / SMPD
            enddo
        enddo
        allocate(map(BOX,BOX,BOX), noise(BOX,BOX,BOX))
        call render_atoms(map, SMPD, xyz, ints, bfacs, HAMP, HBFAC)
        peak = halo_kernel(0., BP, HAMP, HBFAC)
        sig  = peak / PEAK_SNR
        do k = 1,BOX
            do j = 1,BOX
                do i = 1,BOX
                    noise(i,j,k) = gasdev(0., 1.)
                enddo
            enddo
        enddo
        if( l_corr )then
            allocate(sm(BOX,BOX,BOX))
            call gauss_filter3D(noise, CORR_SD / SMPD, sm)
            noise = sm
        endif
        m     = sum(noise) / real(size(noise))
        noise = (noise - m) * sig / sqrt(sum((noise - m)**2) / real(size(noise) - 1))
        map   = map + noise
        ! the noise region: inside the box's inscribed sphere less 2 voxels, farther than 1.5 d_NN from every atom;
        ! the envelope: within the outermost atom's radius plus d_NN / 2
        allocate(region(BOX,BOX,BOX), envelope(BOX,BOX,BOX))
        do k = 1,BOX
            do j = 1,BOX
                do i = 1,BOX
                    r = sqrt(sum((real([i,j,k]) - cen)**2))
                    envelope(i,j,k) = r * SMPD <= real(NSHELL) * DNN + 0.5 * DNN
                    region(i,j,k)   = r <= 0.5 * real(BOX) - 2.
                enddo
            enddo
        enddo
        rv = 1.5 * DNN / SMPD
        do l = 1,n
            lo = max(1, floor(xyz(:,l) - rv))
            hi = min(BOX, ceiling(xyz(:,l) + rv))
            do k = lo(3),hi(3)
                do j = lo(2),hi(2)
                    do i = lo(1),hi(1)
                        if( sum((real([i,j,k]) - xyz(:,l))**2) <= rv * rv ) region(i,j,k) = .false.
                    enddo
                enddo
            enddo
        enddo
    end subroutine make_particle

end module simple_nano_density_tester
