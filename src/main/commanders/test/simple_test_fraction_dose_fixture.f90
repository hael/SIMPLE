!@descr: fraction-dose fixtures shared by testers: a saved motion model with a fraction dose, and the check of a fraction export to its dose-series manifest
module simple_test_fraction_dose_fixture
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,         only: dp, CWD_GLOB
use simple_type_defs,    only: OBJFUN_CC
use simple_string,       only: string
use simple_string_utils, only: int2str, int2str_pad
use simple_syslib,       only: simple_chdir, simple_mkdir
use simple_fileio,       only: file_exists, simple_copy_file, file2rarr
use simple_parameters,   only: parameters
use simple_cmdline,      only: cmdline
use simple_sp_project,   only: sp_project
use simple_motion_model, only: motion_model
use simple_image,        only: image
use simple_dose_series_manifest,       only: dose_series_manifest
use simple_fractionate_model_strategy, only: reconstruct_motion_fractions
use simple_test_utils,   only: assert_true, assert_false, assert_int, assert_char
implicit none
private
public :: check_fraction_dose_export, write_saved_model
public :: MOVIE_NFRAMES, RECON_NPTCLS, SAVED_FRACTION_DOSE

integer, parameter :: MOVIE_NFRAMES       = 12, RECON_NPTCLS = 120
real,    parameter :: SAVED_FRACTION_DOSE = 1.25

contains

    !> Persist realistic model metadata with an EER grouping and target distinct from the saved fraction dose.
    subroutine write_saved_model(params, fname, nframes, dose, kv)
        class(parameters), target, intent(inout) :: params
        type(string),      intent(in)            :: fname
        integer,           intent(in)            :: nframes
        real,              intent(in)            :: dose, kv
        class(motion_model), allocatable :: model
        class(image),        pointer     :: frames(:)
        type(string)                     :: movie
        integer                          :: i, funit
        allocate(model, frames(nframes))
        movie = fname//'.eer'
        open(newunit=funit, file=movie%to_char(), status='replace', access='stream', form='unformatted')
        write(funit) 0
        close(funit)
        do i = 1, nframes
            call frames(i)%new([4,4,1], 1.5, wthreads=.false.)
        enddo
        params%fraction_dose_target = 9.
        params%eer_upsampling = 1
        params%total_dose = real(nframes)*merge(dose,0.,ieee_is_finite(dose))
        call model%new(params, movie, [4,4], 1.5, frames, 8*nframes, 1, kv, merge(dose,0.,ieee_is_finite(dose)), 8)
        model%dose_per_frame = dose
        call model%write(fname//'.star', fname, .false.)
        call model%kill
        do i = 1, nframes
            call frames(i)%kill
        enddo
        deallocate(frames)
    end subroutine write_saved_model

    !> Verify success and metadata-only rejection through the production reconstruction loop.
    subroutine check_fraction_dose_export(root, name, stack, oritab, valid_dose, changed_state)
        use simple_test_flex_phantom_fixture, only: build_flex_pca_phantom_project, &
            &FLEX_PHANTOM_BOX, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
        type(string),     intent(in) :: root, stack, oritab
        character(len=*), intent(in) :: name
        logical,          intent(in) :: valid_dose, changed_state
        class(parameters),           allocatable :: params
        class(cmdline),              allocatable :: cline
        class(sp_project),           allocatable :: project
        class(dose_series_manifest), allocatable :: manifest
        type(string),                allocatable :: projects(:), stacks(:)
        integer,                     allocatable :: ranges(:,:)
        type(string)                             :: arm, directory, model, message, even, odd
        class(image), allocatable :: pproc, lowpass
        real,             pointer :: pproc_pixels(:,:,:), lowpass_pixels(:,:,:)
        real,         allocatable :: fsc(:)
        real(dp)                  :: dlo, dhi
        integer                   :: i, j, nwin, nstates, a, b
        logical                   :: ok
        allocate(params, cline, project, manifest)
        arm = root//'/'//name
        call simple_mkdir(arm)
        call simple_chdir(arm)
        CWD_GLOB = arm%to_char()
        params%box = FLEX_PHANTOM_BOX
        params%smpd = FLEX_PHANTOM_SMPD
        params%mskdiam = FLEX_PHANTOM_MSKDIAM
        params%nparts = 1
        params%nthr = 1
        params%pgrp = 'c1'
        params%cc_objfun = OBJFUN_CC
        params%l_ml_reg = .false.
        params%cc_emit_sigma = 'no'
        model = arm//'/motion.bin'
        call write_saved_model(params, model, MOVIE_NFRAMES, merge(SAVED_FRACTION_DOSE,0.,valid_dose), 300.)
        nwin = 2
        if( valid_dose .and. .not. changed_state ) nwin = 3
        nstates = 1
        if( changed_state ) nstates = 2
        allocate(projects(nwin), stacks(nwin), ranges(2,nwin))
        ranges(:,1) = [1,3]
        ranges(:,2) = [4,6]
        if( nwin == 3 ) ranges(:,3) = [6,8]
        do i = 1, nwin
            directory = arm//'/frames_'//int2str(ranges(1,i))//'-'//int2str(ranges(2,i))
            call simple_mkdir(directory)
            stacks(i) = directory//'/particles.mrcs'
            call simple_copy_file(stack, stacks(i))
            projects(i) = directory//'/fraction.simple'
            call build_flex_pca_phantom_project(projects(i), stacks(i), oritab, nptcls=RECON_NPTCLS)
            call project%read(projects(i))
            call project%os_mic%new(1, is_ptcl=.false.)
            call project%os_mic%set(1, 'mcmodel', model)
            call project%os_mic%set(1, 'nptcls', RECON_NPTCLS)
            call project%os_mic%set_state(1, 1)
            if( changed_state )then
                do j = RECON_NPTCLS-3, RECON_NPTCLS
                    call project%os_ptcl3D%set_state(j, 2)
                enddo
                if( i == 2 )then
                    call assert_int(project%os_ptcl3D%get_eo(1), project%os_ptcl3D%get_eo(RECON_NPTCLS-3), &
                        &'fraction export: state-membership swap keeps the half population')
                    call project%os_ptcl3D%set_state(1, 2)
                    call project%os_ptcl3D%set_state(RECON_NPTCLS-3, 1)
                endif
            endif
            call project%write(projects(i))
            call project%kill
        enddo
        call cline%set('bfac', -120.)
        call cline%set('mirr', 'yes')
        call cline%set('postprocess', 'no')
        params%amsklp = 32.
        params%binwidth = 5
        params%edge = 4
        call reconstruct_motion_fractions(params, cline, projects, ranges, nstates)
        call check_fraction_totals(arm, ranges, nstates, params%box, params%smpd)
        call assert_true(params%amsklp == 32., 'fraction totals: caller low-pass is preserved')
        call assert_int(5, params%binwidth, 'fraction totals: caller binwidth is preserved')
        call assert_int(4, params%edge, 'fraction totals: caller edge is preserved')
        allocate(pproc, lowpass)
        call pproc%new([params%box,params%box,params%box], params%smpd)
        call lowpass%new([params%box,params%box,params%box], params%smpd)
        call assert_true(file_exists(arm//'/fsc_fractions_state01.tsv'), 'fraction export '//name//': raw FSC matrix retained')
        do i = 1, nwin
            call assert_true(file_exists(stacks(i)), 'fraction export '//name//': extracted stack retained')
            directory = arm//'/rec_frames_'//int2str_pad(ranges(1,i),2)//'-'//int2str_pad(ranges(2,i),2)
            call check_cross_fraction_export(directory, ranges, i, nstates, params%box, .not. changed_state)
            call assert_true(file_exists(directory//'/recvol_state01_even_unfil.mrc'), 'fraction export '//name//': raw even half')
            call assert_true(file_exists(directory//'/recvol_state01_odd_unfil.mrc'), 'fraction export '//name//': raw odd half')
            call assert_true(file_exists(directory//'/recvol_state01_pproc.mrc'), &
                &'fraction export '//name//': postprocessed merged map')
            call assert_true(file_exists(directory//'/recvol_state01_lp.mrc'), &
                &'fraction export '//name//': filtered no-B-factor companion')
            call assert_false(file_exists(directory//'/recvol_state01_pproc_mirr.mrc'), &
                &'fraction export '//name//': mirrored output disabled')
            call pproc%read(directory//'/recvol_state01_pproc.mrc')
            call lowpass%read(directory//'/recvol_state01_lp.mrc')
            call pproc%get_rmat_ptr(pproc_pixels)
            call lowpass%get_rmat_ptr(lowpass_pixels)
            call assert_true(all(ieee_is_finite(pproc_pixels)), 'fraction postprocess: finite voxels')
            call assert_true(maxval(abs(pproc_pixels-lowpass_pixels)) <= &
                &1.e-6*max(1.,maxval(abs(lowpass_pixels))), 'fraction postprocess: no sharpening applied')
            fsc = file2rarr(directory//'/fsc_state01.bin')
            call assert_int(FLEX_PHANTOM_BOX/2, size(fsc), 'fraction export '//name//': native FSC shell count')
        enddo
        if( .not. valid_dose .or. changed_state )then
            call assert_false(file_exists(arm//'/dose_series.txt'), 'fraction export '//name//': no invalid manifest published')
        else
            ok = manifest%read(arm//'/dose_series.txt', message)
            call assert_true(ok, 'fraction export: generated manifest is readable')
            if( ok )then
                call assert_int(3, manifest%get_nwindows(), 'fraction export: producer retains all actual intervals')
                ok = manifest%select_windows(message)
                call assert_true(ok, 'fraction export: manifest selects disjoint intervals')
                call assert_int(2, manifest%get_nwindows(), 'fraction export: overlapping predecessor removed')
                call manifest%get_window(2, a, b, dlo, dhi, even, odd)
                call assert_int(6, a, 'fraction export: highest-dose interval retained')
                call assert_true(abs(dlo-6.25_dp) < 1.e-12_dp, 'fraction export: final interval physical dose origin')
                call assert_true(abs(dhi-10._dp) < 1.e-12_dp, 'fraction export: final interval physical dose end')

            endif
        endif
        call pproc%kill
        call lowpass%kill
        call manifest%kill
        call simple_chdir(root)
        CWD_GLOB = root%to_char()
    end subroutine check_fraction_dose_export

    !> Identical input intervals must correlate at unity, independently of their half-map FSC.
    subroutine check_cross_fraction_export(directory, ranges, reference_index, nstates, box, identical)
        type(string), intent(in) :: directory
        integer, intent(in) :: ranges(:,:), reference_index, nstates, box
        logical, intent(in) :: identical
        real, parameter :: FSC_ROUNDTRIP_TOL = 16.*epsilon(1.)
        type(string) :: prefix, expected_header, fname
        character(len=2048) :: header
        real(dp), allocatable :: curves(:)
        real(dp) :: resolution, frequency
        integer :: state, i, unit, ios, shell, j
        allocate(curves(size(ranges,2)-1))
        expected_header = 'shell'//achar(9)//'resolution_A'//achar(9)//'reciprocal_frequency_A^-1'
        do i = 1, size(ranges,2)
            if( i == reference_index ) cycle
            expected_header = expected_header//achar(9)//'FSC_fractions_'//int2str(ranges(1,i))//'-'//int2str(ranges(2,i))
        enddo
        do state = 1, nstates
            prefix = directory//'/fsc_cross_fractions_state'//int2str_pad(state,2)
            call assert_true(file_exists(prefix//'.gp'), 'cross-fraction FSC: reproducible plotting commands')
            call assert_true(file_exists(prefix//'_surface.tsv'), 'cross-fraction FSC: surface data retained')
            fname = prefix//'.tsv'
            open(newunit=unit, file=fname%to_char(), status='old', action='read', iostat=ios)
            call assert_int(0, ios, 'cross-fraction FSC: curve table readable')
            if( ios /= 0 ) cycle
            read(unit,'(A)',iostat=ios) header
            call assert_int(0, ios, 'cross-fraction FSC: header readable')
            call assert_char(expected_header%to_char(), trim(header), 'cross-fraction FSC: ordered targets exclude self')
            do j = 1, box/2
                read(unit,*,iostat=ios) shell, resolution, frequency, curves
                call assert_int(0, ios, 'cross-fraction FSC: complete native shell grid')
                if( ios /= 0 ) exit
                call assert_int(j, shell, 'cross-fraction FSC: native shell indices')
                call assert_true(all(ieee_is_finite(curves)), 'cross-fraction FSC: populated map pairs are finite')
                if( identical ) call assert_true(all(abs(curves-1._dp) < FSC_ROUNDTRIP_TOL), &
                    &'cross-fraction FSC: identical standard maps correlate at unity')
            enddo
            close(unit)
        enddo
    end subroutine check_cross_fraction_export

    !> Array sums independently check half identity, retained overlap and the merged sum's absolute scale.
    subroutine check_fraction_totals(root, ranges, nstates, box, smpd)
        type(string), intent(in) :: root
        integer, intent(in) :: ranges(:,:), nstates, box
        real, intent(in) :: smpd
        real, parameter :: SUM_REL_TOL = 16.*epsilon(1.)
        class(image), allocatable :: volume
        real(dp), allocatable :: expected(:,:,:), half_sum(:,:,:), actual(:,:,:)
        real, allocatable :: envelope(:,:,:), filtered(:,:,:), fsc(:)
        type(string) :: directory, fname, tag, suffix
        integer :: state, half, i
        allocate(volume, expected(box,box,box), half_sum(box,box,box))
        call volume%new([box,box,box], smpd)
        directory = root//'/rec_frames_total'
        call assert_true(file_exists(directory//'/recvol_state01_unfil.mrc'), 'fraction totals: unfiltered merged sum exists')
        do state = 1, nstates
            tag = 'recvol_state'//int2str_pad(state,2)
            half_sum = 0._dp
            do half = 1, 2
                suffix = '_even_unfil.mrc'
                if( half == 2 ) suffix = '_odd_unfil.mrc'
                expected = 0._dp
                do i = 1, size(ranges,2)
                    fname = root//'/rec_frames_'//int2str_pad(ranges(1,i),2)//'-'//int2str_pad(ranges(2,i),2)//'/'//tag//suffix
                    call volume%read(fname)
                    expected = expected + real(volume%get_rmat(),dp)
                enddo
                call volume%read(directory//'/'//tag//suffix)
                actual = real(volume%get_rmat(),dp)
                call assert_true(maxval(abs(actual-expected)) <= SUM_REL_TOL*maxval(abs(expected)), &
                    &'fraction totals: raw half equals the sum of every saved interval without rescaling')
                half_sum = half_sum + actual
            enddo
            call volume%read(directory//'/'//tag//'_unfil.mrc')
            actual = real(volume%get_rmat(),dp)
            call assert_true(maxval(abs(actual-half_sum)) <= SUM_REL_TOL*maxval(abs(half_sum)), &
                &'fraction totals: raw merged volume is even plus odd, without division by two')
            fname = directory//'/automask3D_state'//int2str_pad(state,2)//'.mrc'
            call volume%read(fname)
            envelope = volume%get_rmat()
            call assert_true(all(envelope >= 0. .and. envelope <= 1.), 'fraction totals: density mask range')
            call assert_true(all(envelope == 0. .or. envelope == 1.), 'fraction totals: edge one is binary on the voxel grid')
            call assert_true(any(envelope > 0.) .and. any(envelope == 0.), 'fraction totals: mask includes density and excludes solvent')
            call check_surface_distances(directory, state, envelope, smpd)
            call volume%read(directory//'/'//tag//'_pproc.mrc')
            filtered = volume%get_rmat()
            call assert_true(all(ieee_is_finite(filtered)), 'fraction totals: filtered sum is finite')
            call assert_true(all(pack(filtered,envelope == 0.) == 0.), 'fraction totals: filtered sum retains its density mask')
            call assert_false(file_exists(directory//'/'//tag//'_pproc_mirr.mrc'), 'fraction totals: no mirrored output')
            fname = directory//'/fsc_state'//int2str_pad(state,2)
            fsc = file2rarr(fname//'.bin')
            call assert_int(box/2, size(fsc), 'fraction totals: native FSC shell count')
            call assert_true(all(ieee_is_finite(fsc)), 'fraction totals: spherical FSC is finite')
            call assert_true(file_exists(fname//'.tsv'), 'fraction totals: readable FSC curve is retained')
        enddo
        call volume%kill
    end subroutine check_fraction_totals

    !> Brute-force opposite-label voxel centres verify signs, Euclidean geometry and Angstrom units.
    subroutine check_surface_distances(directory, state, envelope, smpd)
        use simple_test_utils, only: assert_real
        type(string), intent(in) :: directory
        integer, intent(in) :: state
        real, intent(in) :: envelope(:,:,:), smpd
        class(image), allocatable :: volume
        real, allocatable :: distance(:,:,:)
        type(string) :: fname
        integer :: probes(3), n, a, b, c, i, j, k, sqmin, sqdist
        real :: expected, tolerance
        logical :: inside
        n = size(envelope,1)
        allocate(volume)
        call volume%new(shape(envelope), smpd)
        if( state == 1 )then
            fname = directory//'/distance_to_surface.mrc'
        else
            fname = directory//'/distance_to_surface_state'//int2str_pad(state,2)//'.mrc'
        endif
        call assert_true(file_exists(fname), 'surface distance: state-specific output exists')
        if( .not. file_exists(fname) )then
            call volume%kill
            return
        endif
        call volume%read(fname)
        distance = volume%get_rmat()
        call assert_true(all(ieee_is_finite(distance)), 'surface distance: every voxel is finite')
        call assert_true(all((distance > 0.) .eqv. (envelope > 0.5)), 'surface distance: positive inside, negative outside')
        call assert_true(all(distance /= 0.), 'surface distance: opposite voxel centres never coincide')
        tolerance = 16.*epsilon(1.)*sqrt(3.*real(n*n))*smpd
        call assert_real(smpd, minval(abs(distance)), tolerance, 'surface distance: adjacent opposite labels are one voxel apart')
        probes = [1,n/2+1,n]
        do c = 1, 3
            do b = 1, 3
                do a = 1, 3
                    inside = envelope(probes(a),probes(b),probes(c)) > 0.5
                    sqmin = huge(sqmin)
                    do k = 1, n
                        do j = 1, n
                            do i = 1, n
                                if( (envelope(i,j,k) > 0.5) .eqv. inside ) cycle
                                sqdist = (i-probes(a))**2+(j-probes(b))**2+(k-probes(c))**2
                                sqmin = min(sqmin, sqdist)
                            enddo
                        enddo
                    enddo
                    expected = sqrt(real(sqmin))*smpd
                    if( .not. inside ) expected = -expected
                    call assert_real(expected, distance(probes(a),probes(b),probes(c)), tolerance, &
                        &'surface distance: exact nearest opposite voxel, without wrapping or half-voxel offset')
                enddo
            enddo
        enddo
        call volume%kill
    end subroutine check_surface_distances

end module simple_test_fraction_dose_fixture
