!@descr: frozen bootstrap smoothing, persistence and matching/reconstruction separation tests
module simple_sigma2_frozen_tester
use, intrinsic :: iso_fortran_env, only: int64, real32, real64
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use simple_sigma2_frozen, only: prepare_frozen_sigma2, read_frozen_sigma2, frozen_sigma2_band_stats
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, &
    &sigma2_state_create_candidate, sigma2_state_write_particles, sigma2_state_write_global, &
    &sigma2_state_publish, sigma2_state_read_global, SIGMA2_PROV_PSPEC, SIGMA2_PROV_RESIDUAL
use simple_string, only: string
use simple_test_utils, only: assert_true, assert_int, assert_real, tests_failed, enter_fixture, leave_fixture
implicit none
private
public :: run_all_sigma2_frozen_tests
character(len=*), parameter :: STORE = 'sigma2_state.bin'
integer, parameter :: NSH = 16

contains

    subroutine run_all_sigma2_frozen_tests()
        type(string) :: saved, root
        integer :: nfail0
        nfail0 = tests_failed
        call enter_fixture('sigma2_frozen', saved, root)
        call test_smoothing_and_persistence()
        call test_invalid_companions()
        call test_band_readout()
        call test_matching_separation()
        call leave_fixture(saved, root, nfail0)
    end subroutine run_all_sigma2_frozen_tests

    subroutine write_store(spectra, generation, provenance, digest)
        real(real32), intent(in) :: spectra(:,:)
        integer(int64), intent(in) :: generation
        integer, intent(in) :: provenance
        integer(int64), optional, intent(in) :: digest
        type(sigma2_state_header) :: header
        integer(int64) :: layout
        integer :: status
        character(len=512) :: message
        layout = 123_int64
        if( present(digest) ) layout = digest
        call sigma2_state_init_header(header, 1, size(spectra,1), 2, 2*size(spectra,1), 1., &
            &generation, layout, provenance)
        call sigma2_state_create_candidate('candidate.bin', header, status, message)
        call assert_int(0, status, 'create canonical fixture: '//trim(message))
        call sigma2_state_write_particles('candidate.bin', 1, spectra, status, message)
        call assert_int(0, status, 'write fixture records: '//trim(message))
        call sigma2_state_write_global('candidate.bin', spectra, status, message)
        call assert_int(0, status, 'write fixture spectra: '//trim(message))
        call sigma2_state_publish('candidate.bin', STORE, status, message)
        call assert_int(0, status, 'publish fixture: '//trim(message))
    end subroutine write_store

    subroutine test_smoothing_and_persistence()
        real(real32) :: input(NSH,2), expected(NSH)
        real(real32), allocatable :: frozen(:,:), measured(:,:), again(:,:)
        integer :: k, status
        character(len=512) :: message
        input(:,1) = 1.
        input(8,1) = exp(1.)
        input(:,2) = [(exp(-2.+0.1*k+0.002*k*k), k=1,NSH)]
        call write_store(input, 1_int64, SIGMA2_PROV_PSPEC)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_int(0, status, 'prepare smoothed bootstrap: '//trim(message))
        call read_frozen_sigma2(STORE, frozen, status, message)
        call assert_int(0, status, 'read frozen bootstrap: '//trim(message))
        if( status /= 0 ) return
        ! Known seven-point quadratic SG impulse response, independent of the fitting implementation.
        expected = 1.
        expected(5:11) = exp(real([-2,3,6,7,6,3,-2])/21.)
        call assert_true(all(abs(frozen(:,1)-expected)<2.e-6), 'SG matches its closed-form impulse response')
        call assert_true(all(abs(frozen(:,2)-input(:,2))<2.e-6), 'quadratic log spectrum is preserved, including endpoints')
        call sigma2_state_read_global(STORE, measured, status, message)
        call assert_int(0, status, 'read measured spectra')
        call assert_true(all(measured==input), 'smoothing does not change the canonical spectra')
        input = 10.*input
        call write_store(input, 2_int64, SIGMA2_PROV_RESIDUAL)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_int(0, status, 'continuation keeps the bootstrap after a residual commit')
        call read_frozen_sigma2(STORE, again, status, message)
        call assert_int(0, status, 'read after residual commit')
        call assert_true(all(again==frozen), 'both frozen halves survive a measured-sigma update exactly')
        call write_store(input, 3_int64, SIGMA2_PROV_PSPEC)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_int(0, status, 'a rebuilt bootstrap replaces the companion')
        call read_frozen_sigma2(STORE, again, status, message)
        call assert_true(all(abs(again-10.*frozen)<2.e-5), 'rebuilt bootstrap supplies the new frozen spectra')
    end subroutine test_smoothing_and_persistence

    subroutine test_invalid_companions()
        real(real32) :: input(NSH,2)
        real(real32), allocatable :: frozen(:,:)
        integer :: status, unit
        character(len=512) :: message
        input = 1.
        call write_store(input, 4_int64, SIGMA2_PROV_RESIDUAL, 456_int64)
        call read_frozen_sigma2(STORE, frozen, status, message)
        call assert_true(status/=0, 'another particle layout is refused')
        open(newunit=unit, file=STORE//'.frozen', status='old')
        close(unit, status='delete')
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_true(status/=0, 'a residual spectrum cannot substitute for a missing bootstrap')
        call write_store(input, 1_int64, SIGMA2_PROV_PSPEC)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_int(0, status, 'restore companion for corruption checks')
        open(newunit=unit, file=STORE//'.frozen', access='stream', form='unformatted', status='old')
        write(unit, pos=57) ieee_value(1.0_real32, ieee_quiet_nan)
        close(unit)
        call read_frozen_sigma2(STORE, frozen, status, message)
        call assert_true(status/=0, 'nonfinite companion spectrum is refused')
        open(newunit=unit, file=STORE//'.frozen', access='stream', form='unformatted', status='replace')
        write(unit) 'short'
        close(unit)
        call read_frozen_sigma2(STORE, frozen, status, message)
        call assert_true(status/=0, 'truncated companion is refused')
        call write_store(input(:6,:), 1_int64, SIGMA2_PROV_PSPEC)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_true(status/=0, 'fewer than seven shells are refused')
        input(3,1) = 0.
        call write_store(input, 1_int64, SIGMA2_PROV_PSPEC)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_true(status/=0, 'nonpositive bootstrap is refused')
    end subroutine test_invalid_companions

    subroutine test_band_readout()
        real(real32) :: input(NSH,2)
        real(real64) :: mean_log(2), edge_log(2)
        integer :: k, status, nabove(2)
        character(len=512) :: message
        input(:,1) = 2.; input(:,2) = 5.
        call write_store(input, 1_int64, SIGMA2_PROV_PSPEC)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_int(0, status, 'prepare readout anchor')
        input(:,1) = [(2.*exp(real(k)/10.),k=1,NSH)]
        input(:,2) = [(5.*exp(-real(k)/20.),k=1,NSH)]
        call write_store(input, 2_int64, SIGMA2_PROV_RESIDUAL)
        call frozen_sigma2_band_stats(STORE,[3,8],mean_log,edge_log,nabove,status,message)
        call assert_int(0,status,'readout uses the requested band')
        call assert_real(0.55,real(mean_log(1)),1.e-6,'even log ratios average uniformly over shells')
        call assert_real(-0.275,real(mean_log(2)),1.e-6,'odd log ratios stay separate')
        call assert_real(0.8,real(edge_log(1)),1.e-6,'even low-pass edge ratio')
        call assert_real(-0.4,real(edge_log(2)),1.e-6,'odd low-pass edge ratio')
        call assert_int(6,nabove(1),'shells above the anchor are counted, not clipped')
        call assert_int(0,nabove(2),'shells below the anchor')
        call frozen_sigma2_band_stats(STORE,[0,8],mean_log,edge_log,nabove,status,message)
        call assert_true(status/=0,'an invalid diagnostic band is refused')
    end subroutine test_band_readout

    !> Analytic radial fields: loss = sum(k/sigma2 * residual**2) / sum(k/sigma2).
    !! Recreate a worker with a wider search band, as solve2D does between stages.
    subroutine test_matching_separation()
        use simple_defs, only: DPI
        use simple_type_defs, only: OBJFUN_EUCLID
        use simple_parameters, only: parameters
        use simple_polarft_calc, only: polarft_calc
        use simple_euclid_sigma2, only: euclid_sigma2
        use simple_ori, only: ori
        class(parameters), allocatable, target :: params
        class(polarft_calc), allocatable :: calc
        class(euclid_sigma2), allocatable, target :: esig
        type(ori) :: pose
        real :: bootstrap(NSH,2), measured(NSH,2), dists(24), weight, denom, loss
        real, allocatable :: plane_before(:,:), residual_before(:), residual_after(:)
        complex :: particle(12,2:NSH), reference(12,2:NSH)
        integer :: stage, half, k, status, last
        character(len=512) :: message
        allocate(params, calc, esig)
        params%box = 2*NSH
        params%box_crop = 2*NSH
        params%pftsz = 12
        params%ctf = 'no'
        params%oritype = 'ptcl2D'
        params%cc_objfun = OBJFUN_EUCLID
        params%fromp = 1
        params%top = 2
        params%prob_beta = 0.3
        bootstrap(:,1) = 2.
        bootstrap(:,2) = 5.
        call write_store(bootstrap, 1_int64, SIGMA2_PROV_PSPEC)
        call prepare_frozen_sigma2(STORE, status, message)
        call assert_int(0, status, 'prepare matching fixture')
        particle = (1.,0.)
        do k=2,NSH
            reference(:,k) = cmplx(0.1*k,0.)
        enddo
        call pose%new(is_ptcl=.true.)
        call pose%set('state',1)
        call pose%set('class',1)
        do stage=1,2
            last = 4+4*stage
            params%kfromto = [2,last]
            do half=1,2
                measured(:,half) = [(real(stage*k+half),k=1,NSH)]
            enddo
            call write_store(measured, int(stage+1,int64), SIGMA2_PROV_RESIDUAL)
            call calc%new(params,1,[1,2],params%kfromto)
            call esig%new(params,calc,string(STORE),2*NSH)
            call esig%read_global()
            call esig%read_part()
            do half=1,2
                call calc%set_eo(half,half==1)
                call calc%set_ptcl_pft(half,particle)
                call calc%set_ref_pft(1,reference,half==1)
            enddo
            call calc%memoize_refs()
            call calc%memoize_ptcls()
            call calc%memoize_sigma2_weights()
            call esig%prep_fplane_grid(4*NSH,2)
            plane_before = esig%fplane%inv_sig2
            call esig%calc_sigma2(calc,1,pose,'class')
            residual_before = esig%get_sigma2_part(1)
            do half=1,2
                denom=0.; loss=0.
                do k=2,last
                    weight=real(k)/measured(k,half)
                    denom=denom+weight
                    loss=loss+weight*(1.-0.1*k)**2
                enddo
                call calc%gen_euclid_dists(1,half,[0.,0.],dists)
                call assert_true(all(abs(dists-loss/denom)<2.e-5), 'default scoring uses measured spectra')
            enddo
            call esig%freeze_for_scoring(calc)
            do half=1,2
                denom=0.; loss=0.
                do k=2,last
                    weight=real(k)/bootstrap(k,half)
                    denom=denom+weight
                    loss=loss+weight*(1.-0.1*k)**2
                enddo
                call calc%gen_euclid_dists(1,half,[0.,0.],dists)
                call assert_true(all(abs(dists-loss/denom)<2.e-5), 'frozen spectrum determines scoring in each stage')
                call assert_real(12.*denom,real(calc%wsqsum_ptcl_now(half)),1.e-4,'frozen half determines self-term')
                call assert_real(real(DPI)*denom*0.3/2.,calc%ptcl_posterior_temp(half),1.e-4, &
                    &'insertion posterior temperature uses the frozen half')
            enddo
            call esig%calc_sigma2(calc,1,pose,'class')
            residual_after = esig%get_sigma2_part(1)
            call assert_true(all(residual_after==residual_before), 'freezing does not reweight residual records')
            call assert_true(all(abs(residual_after-[(0.5*(1.-0.1*k)**2,k=2,last)])<1.e-6), &
                &'residual records retain their unweighted closed-form values')
            call esig%prep_fplane_grid(4*NSH,2)
            call assert_true(all(esig%sigma2_noise==measured), 'measured spectra remain available for reconstruction')
            call assert_true(all(esig%fplane%inv_sig2==plane_before), 'reconstruction inverse weights are unchanged')
            call calc%kill()
            call esig%kill()
        enddo
        call pose%kill()
    end subroutine test_matching_separation

end module simple_sigma2_frozen_tester
