!@descr: library tests of independent half-set PCG reconstruction against gridding (simple_reconstructor_pcg)
! Two disjoint half sets of a deterministic asymmetric phantom, observed through the
! simulate_particles projection path with seeded white noise, are reconstructed with the
! matrix-free PCG operator (fixed soft support) and with production gridding. Pinned:
! ownership of the halves, the realised SNR and noise independence, solve determinism
! across operator instances, the half-map FSC contract (bounded, low-frequency agreement,
! decay towards Nyquist), strong noiseless recovery, and on the 48-view matrix the
! lambda sweep bracketing a noisy-error optimum that beats the gridding control.
! Nightly (lib_reconstruction): about a hundred PCG solves at box 24.
module simple_pcg_halfset_tester
use simple_defs,              only: dp, OSMPL_PAD_FAC
use simple_image,             only: image
use simple_ori,               only: ori
use simple_oris,              only: oris
use simple_projector,         only: projector
use simple_reconstructor_pcg, only: reconstructor_pcg
use simple_reconstructor,     only: reconstructor
use simple_parameters,        only: parameters
use simple_sp_project,        only: sp_project
use simple_sym,               only: sym
use simple_memoize_ft_maps,   only: forget_ft_maps, memoize_ft_maps
use simple_type_defs,         only: CTFFLAG_NO, ctfparams, fplane_type
use simple_test_utils
use ieee_arithmetic,          only: ieee_is_finite
implicit none
private
public :: run_all_pcg_halfset_tests, build_truth_volume, TRUTH_VOLUME_BOX

integer, parameter :: TRUTH_VOLUME_BOX = 24
integer, parameter :: BOX = TRUTH_VOLUME_BOX
integer, parameter :: NBLOBS = 4
real,    parameter :: CTRS(3,NBLOBS) = reshape([&
    &-5.0, -3.0,  2.0, &
    & 4.0,  5.0, -3.0, &
    & 0.0, -6.0, -5.0, &
    & 3.0, -2.0,  6.0], [3,NBLOBS])
real,    parameter :: SIGMAS(NBLOBS) = [2.0, 2.5, 1.8, 2.2]
real,    parameter :: AMPS(NBLOBS)   = [1.0, 0.8, 0.6, 0.5]
real,    parameter :: SMPD           = 1.5
real,    parameter :: LAMBDA         = 1.e-3
real,    parameter :: MSKRAD         = 10.
real,    parameter :: OBS_SNR        = 0.5
integer, parameter :: EVEN_SEED      = 20260815
integer, parameter :: ODD_SEED       = 20260816
integer, parameter :: MAXITS         = 40
integer, parameter :: N_LOW_SHELLS   = 3
integer, parameter :: N_HIGH_SHELLS  = 3
real(dp), parameter :: SNR_RELTOL         = 0.12_dp
real(dp), parameter :: NOISE_CORR_TOL     = 0.08_dp
real(dp), parameter :: CLEAN_CORR_MIN     = 0.85_dp
real(dp), parameter :: LOW_FSC_MIN        = 0.40_dp
real(dp), parameter :: FSC_LOW_HIGH_GAP   = 0.10_dp
real(dp), parameter :: FSC_BOUND_SLACK    = 1.e-5_dp

contains

    subroutine run_all_pcg_halfset_tests()
        write(*,'(A)') '**** running all pcg_halfset tests ****'
        call test_truth_volume_fixture()
        call test_disjoint_half_ownership()
        call test_independent_observations()
        call test_halfset_fsc_24()
        call test_halfset_matrix_48()
    end subroutine run_all_pcg_halfset_tests

    ! ---- fixture ---------------------------------------------------------------

    !> the asymmetric deterministic Gaussian-blob truth shared with test=pcg_recon
    subroutine build_truth_volume( volume )
        real, allocatable, intent(out) :: volume(:,:,:)
        real    :: ctr, dx, dy, dz
        integer :: b, i, j, k
        allocate(volume(BOX,BOX,BOX), source=0.)
        ctr = real(BOX)/2. + 0.5
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, NBLOBS
                        dx = real(i) - ctr - CTRS(1,b)
                        dy = real(j) - ctr - CTRS(2,b)
                        dz = real(k) - ctr - CTRS(3,b)
                        volume(i,j,k) = volume(i,j,k) + AMPS(b) * exp(-(dx*dx + dy*dy + dz*dz) / (2. * SIGMAS(b)**2))
                    enddo
                enddo
            enddo
        enddo
    end subroutine build_truth_volume

    !> deterministic, non-overlapping even and odd orientation sets from one spiral
    subroutine build_disjoint_halves( nhalf, even_oris, odd_oris, even_ids, odd_ids )
        integer,              intent(in)    :: nhalf
        type(oris),           intent(inout) :: even_oris, odd_oris
        integer, allocatable, intent(out)   :: even_ids(:), odd_ids(:)
        type(ori)  :: o
        type(oris) :: all_oris
        integer    :: i
        call all_oris%new(2*nhalf, .false.)
        call all_oris%spiral()
        call even_oris%new(nhalf, .false.)
        call odd_oris%new(nhalf, .false.)
        call o%new(.false.)
        allocate(even_ids(nhalf), odd_ids(nhalf))
        do i = 1, nhalf
            odd_ids(i)  = 2*i - 1
            even_ids(i) = 2*i
            call all_oris%get_ori(odd_ids(i), o)
            call odd_oris%set_ori(i, o)
            call all_oris%get_ori(even_ids(i), o)
            call even_oris%set_ori(i, o)
        enddo
        call o%kill()
        call all_oris%kill()
    end subroutine build_disjoint_halves

    !> clean projections through the simulate_particles path (padded projection, real-space
    !! clip) and seeded white noise at the requested SNR; planes are native PCG observations
    subroutine build_observations( sampler, orientations, seed, snr, planes, clean_planes, &
                                  &clean_images, noisy_images, noise, realized_snr )
        class(reconstructor_pcg), intent(in)    :: sampler
        type(oris),               intent(inout) :: orientations
        integer,                  intent(in)    :: seed
        real,                     intent(in)    :: snr
        complex, allocatable,     intent(out)   :: planes(:,:,:), clean_planes(:,:,:)
        real,    allocatable,     intent(out)   :: clean_images(:,:,:), noisy_images(:,:,:), noise(:,:,:)
        real(dp),                 intent(out)   :: realized_snr
        type(image)     :: source_volume, padded_projection, projection
        type(ori)       :: o
        type(projector) :: truth_projector
        real, allocatable :: clean(:,:,:), noisy(:,:,:), phantom(:,:,:)
        integer  :: i, lims2(2,2), nprojs
        real(dp) :: clean_mean, noise_mean, noise_power, signal_power
        nprojs = orientations%get_noris()
        lims2  = sampler%get_lims2()
        allocate(planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprojs))
        allocate(clean_planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprojs))
        allocate(clean_images(BOX,BOX,nprojs), noisy_images(BOX,BOX,nprojs))
        allocate(noise(BOX,BOX,nprojs), source=0.)
        allocate(clean(BOX,BOX,1), noisy(BOX,BOX,1))
        call build_truth_volume(phantom)
        call source_volume%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call source_volume%set_rmat(phantom, .false.)
        call truth_projector%new(OSMPL_PAD_FAC*[BOX,BOX,BOX], SMPD, wthreads=.false.)
        call source_volume%pad(truth_projector)
        call truth_projector%fft()
        call truth_projector%expand_cmat()
        call padded_projection%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call projection%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call o%new(.false.)
        call set_fixed_seed(seed)
        signal_power = 0._dp
        noise_power  = 0._dp
        do i = 1, nprojs
            call orientations%get_ori(i, o)
            call truth_projector%fproject_serial(o, padded_projection)
            ! simimg performs this no-CTF Fourier round trip before detector noise
            call padded_projection%ifft()
            call padded_projection%fft()
            call padded_projection%ifft()
            call padded_projection%clip(projection)
            call projection%get_rmat_sub(clean)
            clean_images(:,:,i) = clean(:,:,1)
            call projection%fft()
            clean_planes(:,:,i) = sampler%extract_native_plane(projection)
            call projection%ifft()
            call projection%add_gauran(snr)
            call projection%get_rmat_sub(noisy)
            noisy_images(:,:,i) = noisy(:,:,1)
            noise(:,:,i)        = noisy(:,:,1) - clean(:,:,1)
            clean_mean   = sum(real(clean(:,:,1),dp)) / real(BOX*BOX,dp)
            noise_mean   = sum(real(noise(:,:,i),dp)) / real(BOX*BOX,dp)
            signal_power = signal_power + sum((real(clean(:,:,1),dp)-clean_mean)**2)
            noise_power  = noise_power  + sum((real(noise(:,:,i),dp)-noise_mean)**2)
            call projection%fft()
            planes(:,:,i) = sampler%extract_native_plane(projection)
        enddo
        realized_snr = signal_power / noise_power
        call o%kill()
        call projection%kill()
        call padded_projection%kill()
        call truth_projector%kill_expanded()
        call truth_projector%kill()
        call source_volume%kill()
    end subroutine build_observations

    ! ---- solvers ---------------------------------------------------------------

    !> one half with the matrix-free PCG operator: fixed lambda, optional support,
    !! exactly maxits iterations (rtol=0 restarts at zero and never stops early)
    subroutine reconstruct_half( orientations, planes, lambda, maxits, mskrad, reconstruction, niters )
        type(oris),        intent(inout) :: orientations
        complex,           intent(in)    :: planes(-BOX/2:,-BOX/2:,:)
        real,              intent(in)    :: lambda, mskrad
        integer,           intent(in)    :: maxits
        real, allocatable, intent(out)   :: reconstruction(:,:,:)
        integer,           intent(out)   :: niters
        type(reconstructor_pcg) :: pcgop
        call pcgop%new(BOX, SMPD, lambda)
        if( mskrad > 0. ) call pcgop%set_mask(mskrad)
        call pcgop%prep_particles(orientations, use_ctf=.false.)
        call pcgop%build_operators(.false.)
        allocate(reconstruction(BOX,BOX,BOX), source=0.)
        call pcgop%solve(planes, reconstruction, maxits=maxits, rtol=0., niters=niters)
        call pcgop%kill()
    end subroutine reconstruct_half

    !> matched clean and noisy solves of one half with one operator and residuals ||Ax-y||/||y||
    subroutine reconstruct_half_pair( orientations, clean_planes, noisy_planes, lambda, maxits, &
                                     &clean_rec, noisy_rec, niters, residuals )
        type(oris),        intent(inout) :: orientations
        complex,           intent(in)    :: clean_planes(-BOX/2:,-BOX/2:,:), noisy_planes(-BOX/2:,-BOX/2:,:)
        real,              intent(in)    :: lambda
        integer,           intent(in)    :: maxits
        real, allocatable, intent(out)   :: clean_rec(:,:,:), noisy_rec(:,:,:)
        integer,           intent(out)   :: niters(2)
        real(dp),          intent(out)   :: residuals(2)
        type(reconstructor_pcg) :: pcgop
        call pcgop%new(BOX, SMPD, lambda)
        call pcgop%set_mask(MSKRAD)
        call pcgop%prep_particles(orientations, use_ctf=.false.)
        call pcgop%build_operators(.false.)
        allocate(clean_rec(BOX,BOX,BOX), source=0.)
        allocate(noisy_rec(BOX,BOX,BOX), source=0.)
        call pcgop%solve(clean_planes, clean_rec, maxits=maxits, rtol=0., niters=niters(1))
        call pcgop%solve(noisy_planes, noisy_rec, maxits=maxits, rtol=0., niters=niters(2))
        call forward_residual(pcgop, orientations, clean_planes, clean_rec, residuals(1))
        call forward_residual(pcgop, orientations, noisy_planes, noisy_rec, residuals(2))
        call pcgop%kill()
    end subroutine reconstruct_half_pair

    !> ||A x - y|| / ||y|| over the lims2 disk with the same forward operator
    subroutine forward_residual( pcgop, orientations, planes, volume, rel_resid )
        class(reconstructor_pcg), intent(inout) :: pcgop
        type(oris),               intent(inout) :: orientations
        complex,                  intent(in)    :: planes(-BOX/2:,-BOX/2:,:)
        real,                     intent(in)    :: volume(BOX,BOX,BOX)
        real(dp),                 intent(out)   :: rel_resid
        type(ori) :: o
        complex, allocatable :: predicted(:,:)
        real(dp) :: num, den
        integer  :: h, i, k
        allocate(predicted(-BOX/2:BOX/2,-BOX/2:BOX/2))
        call o%new(.false.)
        call pcgop%set_volume(volume)
        num = 0._dp
        den = 0._dp
        do i = 1, orientations%get_noris()
            call orientations%get_ori(i, o)
            call pcgop%forward_plane(o, predicted)
            do k = -BOX/2, BOX/2
                do h = -BOX/2, BOX/2
                    if( h*h + k*k > (BOX/2)**2 ) cycle
                    num = num + abs(cmplx(predicted(h,k) - planes(h,k,i), kind=dp))**2
                    den = den + abs(cmplx(planes(h,k,i), kind=dp))**2
                enddo
            enddo
        enddo
        rel_resid = sqrt(num / max(den, tiny(den)))
        call o%kill()
    end subroutine forward_residual

    !> solutions of one operator at several fixed iteration counts
    subroutine reconstruct_half_trajectory( orientations, planes, its, mskrad, recs, niters )
        type(oris),           intent(inout) :: orientations
        complex,              intent(in)    :: planes(-BOX/2:,-BOX/2:,:)
        integer,              intent(in)    :: its(:)
        real,                 intent(in)    :: mskrad
        real,    allocatable, intent(out)   :: recs(:,:,:,:)
        integer, allocatable, intent(out)   :: niters(:)
        type(reconstructor_pcg) :: pcgop
        integer :: i
        call pcgop%new(BOX, SMPD, LAMBDA)
        if( mskrad > 0. ) call pcgop%set_mask(mskrad)
        call pcgop%prep_particles(orientations, use_ctf=.false.)
        call pcgop%build_operators(.false.)
        allocate(recs(BOX,BOX,BOX,size(its)), source=0.)
        allocate(niters(size(its)), source=0)
        do i = 1, size(its)
            call pcgop%solve(planes, recs(:,:,:,i), maxits=its(i), rtol=0., niters=niters(i))
        enddo
        call pcgop%kill()
    end subroutine reconstruct_half_trajectory

    !> one half through production gridding (KB insertion, no CTF), the A/B control
    subroutine reconstruct_half_gridding( orientations, observations, reconstruction )
        type(oris),        intent(inout) :: orientations
        real,              intent(in)    :: observations(:,:,:)
        real, allocatable, intent(out)   :: reconstruction(:,:,:)
        type(ctfparams)          :: ctfparms
        type(fplane_type)        :: fplane
        type(image)              :: observation, observation_padded, restored
        type(ori)                :: o
        type(parameters), target :: params
        type(reconstructor)      :: gridder
        type(sp_project)         :: project
        type(sym)                :: c1sym
        integer :: i
        params%box        = BOX
        params%box_crop   = BOX
        params%box_croppd = OSMPL_PAD_FAC * BOX
        params%smpd_crop  = SMPD
        params%nstates    = 1
        params%numlen     = 1
        params%oritype    = 'cls3D'   ! class-volume field: no project dependency
        call gridder%new_accumulator(params, project, expand=.true., wthreads=.false.)
        call gridder%set_ft(.true.)
        call restored%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call observation%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call observation_padded%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD, wthreads=.false.)
        call o%new(.false.)
        call c1sym%new('c1')
        call memoize_ft_maps([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], SMPD)
        ctfparms%smpd    = SMPD
        ctfparms%ctfflag = CTFFLAG_NO
        do i = 1, orientations%get_noris()
            call orientations%get_ori(i, o)
            call observation%set_rmat(observations(:,:,i:i), .false.)
            call observation%pad(observation_padded, backgr=0., antialiasing=.false.)
            call observation_padded%fft()
            call observation_padded%gen_fplane4rec([0,BOX/2], SMPD, ctfparms, [0.,0.], fplane)
            call gridder%insert_plane_oversamp(c1sym, o, fplane)
        enddo
        call gridder%compress_exp()
        call gridder%restore_final(restored)
        reconstruction = restored%get_rmat()
        if( allocated(fplane%cmplx_plane)    ) deallocate(fplane%cmplx_plane)
        if( allocated(fplane%ctfsq_plane)    ) deallocate(fplane%ctfsq_plane)
        if( allocated(fplane%transfer_plane) ) deallocate(fplane%transfer_plane)
        call forget_ft_maps()
        call c1sym%kill()
        call o%kill()
        call restored%kill()
        call observation_padded%kill()
        call observation%kill()
        call gridder%kill()
        call project%kill()
    end subroutine reconstruct_half_gridding

    ! ---- measures --------------------------------------------------------------

    !> P x with the soft spherical support the PCG operator applies (x = P u)
    subroutine apply_support( volume )
        real, intent(inout) :: volume(:,:,:)
        type(image) :: img
        call img%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call img%set_rmat(volume, .false.)
        call img%mask3D_soft(MSKRAD, backgr=0.)
        volume = img%get_rmat()
        call img%kill()
    end subroutine apply_support

    subroutine calc_fsc( even_volume, odd_volume, fsc )
        real,              intent(in)  :: even_volume(:,:,:), odd_volume(:,:,:)
        real, allocatable, intent(out) :: fsc(:)
        type(image) :: even_img, odd_img
        call even_img%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call odd_img%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call even_img%set_rmat(even_volume, .false.)
        call odd_img%set_rmat(odd_volume, .false.)
        call even_img%fft()
        call odd_img%fft()
        allocate(fsc(even_img%get_filtsz()))
        call even_img%fsc(odd_img, fsc)
        call even_img%kill()
        call odd_img%kill()
    end subroutine calc_fsc

    !> means of the lowest and highest shells of one FSC curve
    subroutine fsc_low_high( fsc, low, high )
        real,     intent(in)  :: fsc(:)
        real(dp), intent(out) :: low, high
        low  = sum(real(fsc(1:N_LOW_SHELLS),dp)) / real(N_LOW_SHELLS,dp)
        high = sum(real(fsc(size(fsc)-N_HIGH_SHELLS+1:size(fsc)),dp)) / real(N_HIGH_SHELLS,dp)
    end subroutine fsc_low_high

    pure real(dp) function corr3( a, b )
        real, intent(in) :: a(:,:,:), b(:,:,:)
        real(dp) :: amean, bmean, den
        amean = sum(real(a,dp)) / real(size(a),dp)
        bmean = sum(real(b,dp)) / real(size(b),dp)
        den   = sqrt(sum((real(a,dp)-amean)**2) * sum((real(b,dp)-bmean)**2))
        corr3 = sum((real(a,dp)-amean) * (real(b,dp)-bmean)) / max(den, tiny(den))
    end function corr3

    pure real(dp) function rel_l2_error( volume, truth )
        real, intent(in) :: volume(:,:,:), truth(:,:,:)
        real(dp) :: den
        den = sum(real(truth,dp)**2)
        rel_l2_error = sqrt(sum((real(volume,dp)-real(truth,dp))**2) / max(den, tiny(den)))
    end function rel_l2_error

    pure real(dp) function rel_error( actual, expected )
        real(dp), intent(in) :: actual, expected
        rel_error = abs(actual - expected) / max(abs(expected), tiny(expected))
    end function rel_error

    pure logical function interior_minimum( values )
        real(dp), intent(in) :: values(:)
        integer :: loc(1)
        loc = minloc(values)
        interior_minimum = loc(1) > 1 .and. loc(1) < size(values)
    end function interior_minimum

    ! ---- tests -----------------------------------------------------------------

    !> the fixture is reproducible and asymmetric (a symmetric phantom hides orientation bugs)
    subroutine test_truth_volume_fixture()
        real, allocatable :: phantom(:,:,:), again(:,:,:)
        real(dp) :: total, sumsq, den
        write(*,'(A)') 'test_truth_volume_fixture'
        call build_truth_volume(phantom)
        call build_truth_volume(again)
        call assert_true(all(phantom == again), 'the phantom is bit-identical when regenerated')
        call assert_true(minval(phantom) >= 0., 'the Gaussian phantom is non-negative')
        total = sum(real(phantom,dp))
        sumsq = sum(real(phantom,dp)**2)
        ! fingerprints of the established PCG fixture (also test=pcg_recon)
        call assert_true(abs(total - 460.9012_dp) < 1.e-2_dp, 'phantom sum 460.90')
        call assert_true(abs(sumsq - 127.6201_dp) < 1.e-2_dp, 'phantom squared norm 127.62')
        call assert_true(abs(real(maxval(phantom),dp) - 0.910925_dp) < 1.e-4_dp, 'phantom maximum 0.9109')
        den = sumsq
        call assert_true(sqrt(sum((real(phantom,dp) - real(phantom(BOX:1:-1,:,:),dp))**2) / den) > 0.5_dp, &
            &'phantom is asymmetric under x reflection')
        call assert_true(sqrt(sum((real(phantom,dp) - real(phantom(BOX:1:-1,BOX:1:-1,:),dp))**2) / den) > 0.5_dp, &
            &'phantom is asymmetric under a 180-degree z rotation')
    end subroutine test_truth_volume_fixture

    subroutine test_disjoint_half_ownership()
        type(oris) :: even_oris, odd_oris
        integer, allocatable :: even_ids(:), odd_ids(:)
        integer :: i
        logical :: disjoint
        write(*,'(A)') 'test_disjoint_half_ownership'
        call build_disjoint_halves(24, even_oris, odd_oris, even_ids, odd_ids)
        call assert_int(24, even_oris%get_noris(), '24 even orientations')
        call assert_int(24, odd_oris%get_noris(),  '24 odd orientations')
        call assert_true(all(mod(even_ids,2) == 0) .and. all(mod(odd_ids,2) == 1), 'ownership follows the parity split')
        disjoint = .true.
        do i = 1, 24
            if( any(even_ids(i) == odd_ids) ) disjoint = .false.
        enddo
        call assert_true(disjoint, 'even and odd ownership do not overlap')
        call even_oris%kill()
        call odd_oris%kill()
    end subroutine test_disjoint_half_ownership

    !> the observation builder realises the requested SNR, draws independent noise per half
    !! and leaves the other half's data untouched
    subroutine test_independent_observations()
        type(oris) :: even_oris, odd_oris
        type(reconstructor_pcg) :: sampler
        complex, allocatable :: even_planes(:,:,:), even_clean(:,:,:), even_snapshot(:,:,:)
        complex, allocatable :: odd_planes(:,:,:), odd_clean(:,:,:)
        real,    allocatable :: ec_img(:,:,:), en_img(:,:,:), oc_img(:,:,:), on_img(:,:,:)
        real,    allocatable :: even_noise(:,:,:), odd_noise(:,:,:)
        integer, allocatable :: even_ids(:), odd_ids(:)
        real(dp) :: even_snr, odd_snr
        write(*,'(A)') 'test_independent_observations'
        call build_disjoint_halves(24, even_oris, odd_oris, even_ids, odd_ids)
        call sampler%new(BOX, SMPD, LAMBDA)
        call build_observations(sampler, even_oris, EVEN_SEED, OBS_SNR, even_planes, even_clean, ec_img, en_img, even_noise, even_snr)
        even_snapshot = even_planes
        call build_observations(sampler, odd_oris, ODD_SEED, OBS_SNR, odd_planes, odd_clean, oc_img, on_img, odd_noise, odd_snr)
        call sampler%kill()
        call assert_true(all(even_planes == even_snapshot), 'building the odd observations leaves the even data untouched')
        call assert_true(rel_error(even_snr, real(OBS_SNR,dp)) < SNR_RELTOL, 'even half realises the requested SNR (12%)')
        call assert_true(rel_error(odd_snr,  real(OBS_SNR,dp)) < SNR_RELTOL, 'odd half realises the requested SNR (12%)')
        call assert_true(abs(corr3(even_noise, odd_noise)) < NOISE_CORR_TOL, 'even and odd noise are uncorrelated')
        call assert_true(any(even_noise /= 0.), 'noise was added to the even half')
        call assert_true(any(even_planes /= even_clean), 'noisy planes differ from clean planes')
        call even_oris%kill()
        call odd_oris%kill()
    end subroutine test_independent_observations

    !> 24 views per half: PCG halves versus gridding halves, the FSC contract
    subroutine test_halfset_fsc_24()
        integer, parameter :: NHALF = 24
        integer, parameter :: TRAJ_ITS(6) = [1, 2, 4, 8, 16, 40]
        type(oris) :: even_oris, odd_oris
        type(reconstructor_pcg) :: sampler
        complex, allocatable :: even_planes(:,:,:), even_clean(:,:,:), odd_planes(:,:,:), odd_clean(:,:,:)
        real,    allocatable :: ec_img(:,:,:), en_img(:,:,:), oc_img(:,:,:), on_img(:,:,:)
        real,    allocatable :: even_noise(:,:,:), odd_noise(:,:,:), truth(:,:,:), fsc(:), grid_fsc(:)
        real,    allocatable :: even_clean_vol(:,:,:), odd_clean_vol(:,:,:), even_vol(:,:,:), odd_vol(:,:,:)
        real,    allocatable :: traj_even(:,:,:,:), traj_odd(:,:,:,:)
        real,    allocatable :: g_even_clean(:,:,:), g_odd_clean(:,:,:), g_even(:,:,:), g_odd(:,:,:)
        integer, allocatable :: even_ids(:), odd_ids(:), traj_even_niters(:), traj_odd_niters(:)
        real(dp) :: even_snr, odd_snr, low, high, glow, ghigh
        real(dp) :: traj_corr(2,6)
        integer  :: ec_n, oc_n, e_n, o_n, i
        write(*,'(A)') 'test_halfset_fsc_24'
        call build_disjoint_halves(NHALF, even_oris, odd_oris, even_ids, odd_ids)
        call sampler%new(BOX, SMPD, LAMBDA)
        call build_observations(sampler, even_oris, EVEN_SEED, OBS_SNR, even_planes, even_clean, ec_img, en_img, even_noise, even_snr)
        call build_observations(sampler, odd_oris,  ODD_SEED,  OBS_SNR, odd_planes,  odd_clean,  oc_img, on_img, odd_noise,  odd_snr)
        call sampler%kill()
        ! each half owns its operator, observations and solution
        call reconstruct_half(even_oris, even_clean, LAMBDA, MAXITS, MSKRAD, even_clean_vol, ec_n)
        call reconstruct_half(odd_oris,  odd_clean,  LAMBDA, MAXITS, MSKRAD, odd_clean_vol,  oc_n)
        call reconstruct_half(even_oris, even_planes, LAMBDA, MAXITS, MSKRAD, even_vol, e_n)
        call reconstruct_half(odd_oris,  odd_planes,  LAMBDA, MAXITS, MSKRAD, odd_vol,  o_n)
        call assert_true(all([ec_n,oc_n,e_n,o_n] == MAXITS), 'rtol=0 solves run exactly maxits iterations')
        call assert_true(all(ieee_is_finite(even_clean_vol)) .and. all(ieee_is_finite(odd_clean_vol)), 'clean half maps are finite')
        call assert_true(all(ieee_is_finite(even_vol)) .and. all(ieee_is_finite(odd_vol)), 'noisy half maps are finite')
        call assert_true(any(even_vol /= odd_vol), 'independently reconstructed half maps differ')
        ! a second operator instance at the same counts reproduces the solves exactly
        call reconstruct_half_trajectory(even_oris, even_clean, TRAJ_ITS, MSKRAD, traj_even, traj_even_niters)
        call reconstruct_half_trajectory(odd_oris,  odd_clean,  TRAJ_ITS, MSKRAD, traj_odd,  traj_odd_niters)
        call assert_true(all(traj_even_niters == TRAJ_ITS) .and. all(traj_odd_niters == TRAJ_ITS), &
            &'trajectory solves run their requested counts')
        call assert_true(rel_l2_error(traj_even(:,:,:,6), even_clean_vol) < 1.e-5_dp .and. &
            &rel_l2_error(traj_odd(:,:,:,6), odd_clean_vol) < 1.e-5_dp, &
            &'a fresh operator at 40 iterations reproduces the baseline solve (rel L2 < 1e-5)')
        ! gridding control on the same real-space observations
        call reconstruct_half_gridding(even_oris, ec_img, g_even_clean)
        call reconstruct_half_gridding(odd_oris,  oc_img, g_odd_clean)
        call reconstruct_half_gridding(even_oris, en_img, g_even)
        call reconstruct_half_gridding(odd_oris,  on_img, g_odd)
        call assert_true(all(ieee_is_finite(g_even)) .and. all(ieee_is_finite(g_odd)) .and. &
            &all(ieee_is_finite(g_even_clean)) .and. all(ieee_is_finite(g_odd_clean)), 'gridding half maps are finite')
        call assert_true(any(g_even /= g_odd), 'gridding half maps differ')
        ! the constrained solve x = P u is assessed against P x_truth
        call build_truth_volume(truth)
        call apply_support(truth)
        do i = 1, 6
            traj_corr(1,i) = corr3(traj_even(:,:,:,i), truth)
            traj_corr(2,i) = corr3(traj_odd(:,:,:,i),  truth)
        enddo
        call assert_true(all(ieee_is_finite(traj_corr)), 'trajectory correlations are finite')
        call assert_true(maxval(traj_corr(1,:)) > CLEAN_CORR_MIN .and. maxval(traj_corr(2,:)) > CLEAN_CORR_MIN, &
            &'noiseless PCG halves recover the supported truth (corr > 0.85)')
        call calc_fsc(even_vol, odd_vol, fsc)
        call apply_support(g_even)
        call apply_support(g_odd)
        call calc_fsc(g_even, g_odd, grid_fsc)
        call assert_true(size(fsc) >= N_LOW_SHELLS + N_HIGH_SHELLS, 'FSC has enough shells')
        call assert_true(all(ieee_is_finite(fsc)) .and. all(ieee_is_finite(grid_fsc)), 'FSC curves are finite')
        call assert_true(all(real(fsc,dp) >= -1._dp - FSC_BOUND_SLACK) .and. all(real(fsc,dp) <= 1._dp + FSC_BOUND_SLACK), &
            &'FSC lies in [-1,1]')
        call fsc_low_high(fsc, low, high)
        call fsc_low_high(grid_fsc, glow, ghigh)
        call assert_true(low > LOW_FSC_MIN, 'noisy half maps agree at low frequency (FSC > 0.4)')
        call assert_true(low - high > FSC_LOW_HIGH_GAP, 'FSC decays from low to high frequency')
        write(*,'(A,2(ES12.4,1X))') '  realised even/odd SNR:        ', even_snr, odd_snr
        write(*,'(A,6(F7.4,1X))')   '  clean even corr vs iterations:', traj_corr(1,:)
        write(*,'(A,2(F7.4,1X))')   '  PCG low/high FSC:             ', low, high
        write(*,'(A,2(F7.4,1X))')   '  gridding low/high FSC:        ', glow, ghigh
        call even_oris%kill()
        call odd_oris%kill()
    end subroutine test_halfset_fsc_24

    !> 48 views per half: iteration trajectories with and without support, the lambda sweep,
    !! and the gridding control; the noisy raw-L2 optimum must be interior to the sweep and
    !! beat gridding
    subroutine test_halfset_matrix_48()
        integer, parameter :: NHALF = 48
        integer, parameter :: TRAJ_ITS(6) = [1, 2, 4, 8, 16, 40]
        integer, parameter :: OPEN_ITS(3) = [4, 8, 40]
        integer, parameter :: NLAMBDA = 13
        real,    parameter :: LAMBDAS(NLAMBDA) = [1.e-3, 1.e-2, 1.e-1, 1., 10., 100., &
            &1.e3, 2.e3, 3.e3, 5.e3, 1.e4, 1.e5, 1.e6]
        type(oris) :: even_oris, odd_oris
        type(reconstructor_pcg) :: sampler
        complex, allocatable :: even_planes(:,:,:), even_clean(:,:,:), odd_planes(:,:,:), odd_clean(:,:,:)
        real,    allocatable :: ec_img(:,:,:), en_img(:,:,:), oc_img(:,:,:), on_img(:,:,:)
        real,    allocatable :: even_noise(:,:,:), odd_noise(:,:,:), truth(:,:,:), fsc(:)
        real,    allocatable :: ce_traj(:,:,:,:), co_traj(:,:,:,:), ne_traj(:,:,:,:), no_traj(:,:,:,:)
        real,    allocatable :: oce(:,:,:,:), oco(:,:,:,:), one(:,:,:,:), ono(:,:,:,:)
        real,    allocatable :: ce(:,:,:), co(:,:,:), ne(:,:,:), no(:,:,:)
        real,    allocatable :: g_ce(:,:,:), g_co(:,:,:), g_ne(:,:,:), g_no(:,:,:)
        integer, allocatable :: even_ids(:), odd_ids(:), n1(:), n2(:), n3(:), n4(:), m1(:), m2(:), m3(:), m4(:)
        real(dp) :: even_snr, odd_snr, noise_corr, traj_corr(4,6), open_corr(4,3), traj_fsc(2,6), open_fsc(2,3)
        real(dp) :: lambda_corr(4,NLAMBDA), lambda_resid(4,NLAMBDA), lambda_err(4,NLAMBDA), lambda_fsc(2,NLAMBDA)
        real(dp) :: grid_corr(4), grid_err(4), grid_fsc(2), e_res(2), o_res(2)
        integer  :: i, e_n(2), o_n(2), best(2)
        logical  :: its_ok
        write(*,'(A)') 'test_halfset_matrix_48'
        call build_disjoint_halves(NHALF, even_oris, odd_oris, even_ids, odd_ids)
        call sampler%new(BOX, SMPD, LAMBDA)
        call build_observations(sampler, even_oris, EVEN_SEED, OBS_SNR, even_planes, even_clean, ec_img, en_img, even_noise, even_snr)
        call build_observations(sampler, odd_oris,  ODD_SEED,  OBS_SNR, odd_planes,  odd_clean,  oc_img, on_img, odd_noise,  odd_snr)
        call sampler%kill()
        call build_truth_volume(truth)
        call apply_support(truth)
        noise_corr = corr3(even_noise, odd_noise)
        call assert_true(rel_error(even_snr, real(OBS_SNR,dp)) < SNR_RELTOL .and. &
            &rel_error(odd_snr, real(OBS_SNR,dp)) < SNR_RELTOL, '48-view halves realise the requested SNR')
        call assert_true(abs(noise_corr) < NOISE_CORR_TOL, '48-view half noises are uncorrelated')
        ! supported trajectories
        call reconstruct_half_trajectory(even_oris, even_clean,  TRAJ_ITS, MSKRAD, ce_traj, n1)
        call reconstruct_half_trajectory(odd_oris,  odd_clean,   TRAJ_ITS, MSKRAD, co_traj, n2)
        call reconstruct_half_trajectory(even_oris, even_planes, TRAJ_ITS, MSKRAD, ne_traj, n3)
        call reconstruct_half_trajectory(odd_oris,  odd_planes,  TRAJ_ITS, MSKRAD, no_traj, n4)
        call assert_true(all(n1 == TRAJ_ITS) .and. all(n2 == TRAJ_ITS) .and. all(n3 == TRAJ_ITS) .and. all(n4 == TRAJ_ITS), &
            &'supported trajectories run their requested counts')
        do i = 1, 6
            traj_corr(:,i) = [corr3(ce_traj(:,:,:,i), truth), corr3(co_traj(:,:,:,i), truth), &
                             &corr3(ne_traj(:,:,:,i), truth), corr3(no_traj(:,:,:,i), truth)]
            call calc_fsc(ne_traj(:,:,:,i), no_traj(:,:,:,i), fsc)
            call fsc_low_high(fsc, traj_fsc(1,i), traj_fsc(2,i))
        enddo
        ! unsupported trajectories
        call reconstruct_half_trajectory(even_oris, even_clean,  OPEN_ITS, 0., oce, m1)
        call reconstruct_half_trajectory(odd_oris,  odd_clean,   OPEN_ITS, 0., oco, m2)
        call reconstruct_half_trajectory(even_oris, even_planes, OPEN_ITS, 0., one, m3)
        call reconstruct_half_trajectory(odd_oris,  odd_planes,  OPEN_ITS, 0., ono, m4)
        call assert_true(all(m1 == OPEN_ITS) .and. all(m2 == OPEN_ITS) .and. all(m3 == OPEN_ITS) .and. all(m4 == OPEN_ITS), &
            &'unsupported trajectories run their requested counts')
        do i = 1, 3
            open_corr(:,i) = [corr3(oce(:,:,:,i), truth), corr3(oco(:,:,:,i), truth), &
                             &corr3(one(:,:,:,i), truth), corr3(ono(:,:,:,i), truth)]
            call calc_fsc(one(:,:,:,i), ono(:,:,:,i), fsc)
            call fsc_low_high(fsc, open_fsc(1,i), open_fsc(2,i))
        enddo
        ! lambda sweep at a fixed 40-iteration budget
        its_ok = .true.
        do i = 1, NLAMBDA
            call reconstruct_half_pair(even_oris, even_clean, even_planes, LAMBDAS(i), MAXITS, ce, ne, e_n, e_res)
            call reconstruct_half_pair(odd_oris,  odd_clean,  odd_planes,  LAMBDAS(i), MAXITS, co, no, o_n, o_res)
            if( any(e_n /= MAXITS) .or. any(o_n /= MAXITS) ) its_ok = .false.
            lambda_corr(:,i)  = [corr3(ce,truth), corr3(co,truth), corr3(ne,truth), corr3(no,truth)]
            lambda_resid(:,i) = [e_res(1), o_res(1), e_res(2), o_res(2)]
            lambda_err(:,i)   = [rel_l2_error(ce,truth), rel_l2_error(co,truth), rel_l2_error(ne,truth), rel_l2_error(no,truth)]
            call calc_fsc(ne, no, fsc)
            call fsc_low_high(fsc, lambda_fsc(1,i), lambda_fsc(2,i))
        enddo
        call assert_true(its_ok, 'every lambda-sweep solve ran exactly 40 iterations')
        ! gridding control
        call reconstruct_half_gridding(even_oris, ec_img, g_ce)
        call reconstruct_half_gridding(odd_oris,  oc_img, g_co)
        call reconstruct_half_gridding(even_oris, en_img, g_ne)
        call reconstruct_half_gridding(odd_oris,  on_img, g_no)
        call apply_support(g_ce)
        call apply_support(g_co)
        call apply_support(g_ne)
        call apply_support(g_no)
        grid_corr = [corr3(g_ce,truth), corr3(g_co,truth), corr3(g_ne,truth), corr3(g_no,truth)]
        grid_err  = [rel_l2_error(g_ce,truth), rel_l2_error(g_co,truth), rel_l2_error(g_ne,truth), rel_l2_error(g_no,truth)]
        call calc_fsc(g_ne, g_no, fsc)
        call fsc_low_high(fsc, grid_fsc(1), grid_fsc(2))
        best(1) = minloc(lambda_err(3,:), dim=1)
        best(2) = minloc(lambda_err(4,:), dim=1)
        ! the contract
        call assert_true(all(ieee_is_finite(traj_corr)) .and. all(ieee_is_finite(traj_fsc)) .and. &
            &all(ieee_is_finite(open_corr)) .and. all(ieee_is_finite(open_fsc)) .and. &
            &all(ieee_is_finite(lambda_corr)) .and. all(ieee_is_finite(lambda_resid)) .and. &
            &all(ieee_is_finite(lambda_err)) .and. all(ieee_is_finite(lambda_fsc)) .and. &
            &all(ieee_is_finite(grid_corr)) .and. all(ieee_is_finite(grid_err)) .and. all(ieee_is_finite(grid_fsc)), &
            &'every matrix statistic is finite')
        call assert_true(maxval(traj_corr(1,:)) > CLEAN_CORR_MIN .and. maxval(traj_corr(2,:)) > CLEAN_CORR_MIN, &
            &'48-view noiseless PCG halves recover the supported truth (corr > 0.85)')
        call assert_true(interior_minimum(lambda_err(3,:)) .and. interior_minimum(lambda_err(4,:)), &
            &'the noisy raw-L2 lambda optimum is bracketed by the sweep')
        call assert_true(lambda_err(3,best(1)) < grid_err(3) .and. lambda_err(4,best(2)) < grid_err(4), &
            &'the best PCG lambda beats the gridding control on noisy raw L2 error')
        write(*,'(A,2(ES12.4,1X),A,ES12.4)') '  realised SNR: ', even_snr, odd_snr, ' noise corr: ', noise_corr
        do i = 1, 6
            write(*,'(A,I3,A,4(F7.4,1X),A,2(F7.4,1X))') '  its ', TRAJ_ITS(i), ' corr ce/co/ne/no: ', traj_corr(:,i), &
                &' FSC low/high: ', traj_fsc(:,i)
        enddo
        do i = 1, 3
            write(*,'(A,I3,A,4(F7.4,1X),A,2(F7.4,1X))') '  open its ', OPEN_ITS(i), ' corr: ', open_corr(:,i), &
                &' FSC low/high: ', open_fsc(:,i)
        enddo
        do i = 1, NLAMBDA
            write(*,'(A,ES9.2,A,4(F7.4,1X),A,4(F7.4,1X),A,2(F7.4,1X))') '  lambda ', LAMBDAS(i), ' corr: ', lambda_corr(:,i), &
                &' rel L2: ', lambda_err(:,i), ' FSC: ', lambda_fsc(:,i)
        enddo
        write(*,'(A,4(F7.4,1X),A,4(F7.4,1X),A,2(F7.4,1X))') '  gridding corr: ', grid_corr, ' rel L2: ', grid_err, &
            &' FSC: ', grid_fsc
        write(*,'(A,2(ES9.2,1X))') '  best noisy lambda even/odd: ', LAMBDAS(best)
        call even_oris%kill()
        call odd_oris%kill()
    end subroutine test_halfset_matrix_48

end module simple_pcg_halfset_tester
