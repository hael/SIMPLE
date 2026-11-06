!@descr: tests of fraction dose metadata and the unmasked reconstruction-to-manifest workflow
submodule (simple_refine_motion_model_strategy) simple_refine_motion_model_tester
use simple_defs, only: CWD_GLOB
use simple_image, only: image
use simple_dose_series_manifest, only: dose_series_manifest
use simple_test_fraction_dose_fixture, only: check_fraction_dose_export, write_saved_model, MOVIE_NFRAMES,&
    &RECON_NPTCLS, SAVED_FRACTION_DOSE
use simple_test_utils, only: assert_true, assert_false, assert_int, assert_char, &
    &enter_fixture, leave_fixture, tests_failed
implicit none

integer, parameter :: METADATA_BOX = 8
real, parameter :: METADATA_SMPD = 2.

contains

    !> Exercise metadata from the actual selected state-1 particles and saved decoded-fraction models.
    module subroutine run_all_refine_motion_model_tests()
        class(parameters), allocatable :: params
        class(sp_project), allocatable :: project
        class(dose_series_manifest), allocatable :: manifest
        type(string) :: cwd, root, model1, model2, message, even, odd
        integer, allocatable :: rows(:), halves(:), pinds(:)
        integer :: nfail0, i, a, b, nselected
        integer :: ranges(2,2)
        real(dp) :: dlo, dhi, smpd, kv
        integer :: box
        logical :: ok
        nfail0 = tests_failed
        call enter_fixture('fraction_dose_metadata', cwd, root)
        allocate(params, project, manifest)
        params%box = METADATA_BOX
        params%smpd = METADATA_SMPD
        model1 = root//'/motion1.bin'
        model2 = root//'/motion2.bin'
        call write_saved_model(params, model1, MOVIE_NFRAMES, SAVED_FRACTION_DOSE, 300.)
        call write_saved_model(params, model2, MOVIE_NFRAMES, SAVED_FRACTION_DOSE, 300.)
        call project%os_mic%new(2, is_ptcl=.false.)
        call project%os_stk%new(2, is_ptcl=.false.)
        call project%os_ptcl3D%new(4, is_ptcl=.true.)
        do i = 1, 2
            call project%os_mic%set(i, 'nptcls', 2)
            call project%os_mic%set_state(i, 1)
            call project%os_stk%set(i, 'nptcls', 2)
            call project%os_stk%set_state(i, 1)
        enddo
        call project%os_mic%set(1, 'mcmodel', model1)
        call project%os_mic%set(2, 'mcmodel', model2)
        do i = 1, 4
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'eo', mod(i-1,2))
            call project%os_ptcl3D%set(i, 'stkind', (i-1)/2+1)
        enddo
        pinds = [1,2,3,4]
        ranges = reshape([4,6,7,9], [2,2])
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_true(ok, 'fraction metadata: valid saved models')
        call assert_int(4, size(rows), 'fraction metadata: actual state-1 particle count')
        call assert_true(all(rows == [1,2,3,4]), 'fraction metadata: selected row identity')
        call assert_true(all(halves == [0,1,0,1]), 'fraction metadata: saved half labels')
        if( ok )then
            call manifest%get_geometry(box, smpd, kv)
            call assert_int(METADATA_BOX, box, 'fraction metadata: native box')
            call assert_true(abs(kv-300._dp) < 1.e-12_dp, 'fraction metadata: voltage from saved model')
            call manifest%get_window(1, a, b, dlo, dhi, even, odd)
            call assert_int(4, a, 'fraction metadata: original start fraction')
            call assert_true(abs(dlo-3.75_dp) < 1.e-12_dp, 'fraction metadata: preserves dose before fromf')
            call assert_true(abs(dhi-7.5_dp) < 1.e-12_dp, 'fraction metadata: EER grouping is not applied twice')
            call assert_char('rec_frames_04-06/recvol_state01_even_unfil.mrc', even%to_char(), 'fraction metadata: portable half path')
            call manifest%get_window(2, a, b, dlo, dhi, even, odd)
            call assert_true(abs(dlo-7.5_dp) < 1.e-12_dp, 'fraction metadata: adjacent physical dose boundary')
            call assert_true(abs(dhi-11.25_dp) < 1.e-12_dp, 'fraction metadata: saved dose wins over target dose')
            call check_manifest_partition(manifest, root)
        endif
        call write_saved_model(params, model2, MOVIE_NFRAMES, 2., 300.)
        call project%os_ptcl3D%set_state(3, 2)
        call project%os_ptcl3D%set_state(4, 2)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_true(ok, 'fraction metadata: another state with a different dose is irrelevant')
        call assert_int(2, size(rows), 'fraction metadata: counts describe state 1 only')
        call project%os_ptcl3D%set_state(3, 1)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: contributing micrograph dose mismatch is rejected')
        call assert_true(message%has_substr('dose per fraction differs'), 'fraction metadata: dose rejection explains the reason')
        call project%os_ptcl3D%set_state(4, 1)
        call project%os_ptcl3D%set(1, 'updatecnt', 1)
        call project%os_ptcl3D%set(2, 'updatecnt', 1)
        nselected = 0
        call project%os_ptcl3D%sample4rec([1,4], nselected, pinds)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_true(ok, 'fraction metadata: unselected rows cannot invalidate the dose')
        call assert_int(2, size(rows), 'fraction metadata: uses sample4rec update coverage')
        pinds = [1,2,3,4]
        call write_saved_model(params, model2, MOVIE_NFRAMES, SAVED_FRACTION_DOSE, 200.)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: contributing voltage mismatch is rejected')
        call assert_true(message%has_substr('voltage differs'), 'fraction metadata: voltage rejection explains the reason')
        call write_saved_model(params, model2, MOVIE_NFRAMES, 0., 300.)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: zero saved dose has no physical manifest')
        call write_saved_model(params, model2, MOVIE_NFRAMES, ieee_value(0.,ieee_quiet_nan), 300.)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: nonfinite saved dose is rejected')
        call write_saved_model(params, model2, 6, SAVED_FRACTION_DOSE, 300.)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: model must cover the requested last fraction')
        call assert_true(message%has_substr('does not cover'), 'fraction metadata: coverage rejection explains the reason')
        call project%os_mic%set(2, 'mcmodel', root//'/missing.bin')
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: missing model skips metadata without reading it')
        call project%os_ptcl3D%set_state(3, 2)
        call project%os_ptcl3D%set_state(4, 2)
        call project%os_ptcl3D%set(2, 'eo', 0)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: both state-1 halves need contributors')
        call project%os_ptcl3D%set_state(1, 2)
        call project%os_ptcl3D%set_state(2, 2)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: absent state 1 has no manifest')
        call project%os_ptcl3D%set_state(1, 1)
        call project%os_ptcl3D%set_state(2, 1)
        call project%os_ptcl3D%set(2, 'eo', 1)
        call project%os_mic%new(3, is_ptcl=.false.)
        ok = prepare_fraction_dose_manifest(params, project, pinds, ranges, manifest, rows, halves, message)
        call assert_false(ok, 'fraction metadata: missing particle counts cannot define a sparse stack mapping')
        call assert_false(file_exists(root//'/dose_series.txt'), 'fraction metadata: preparation alone never publishes a manifest')
        call manifest%kill
        call project%kill
        call leave_fixture(cwd, root, nfail0)
    end subroutine run_all_refine_motion_model_tests

    !> Verify serialized particle counts and the independently computed fingerprint of this fixture.
    subroutine check_manifest_partition(manifest, root)
        class(dose_series_manifest), intent(inout) :: manifest
        type(string), intent(in) :: root
        class(image), allocatable :: volume
        type(string) :: even, odd, message, fname
        real(dp) :: dlo, dhi
        integer :: a, b, i, funit, ios, count_value
        character(len=512) :: line
        character(len=32) :: prefix, key
        logical :: ok, found_digest, found_count
        allocate(volume)
        call volume%new([METADATA_BOX,METADATA_BOX,METADATA_BOX], METADATA_SMPD)
        volume = 1.
        do i = 1, manifest%get_nwindows()
            call manifest%get_window(i, a, b, dlo, dhi, even, odd)
            call simple_mkdir(string('rec_frames_'//int2str_pad(a,2)//'-'//int2str_pad(b,2)))
            call volume%write(even)
            call volume%write(odd)
        enddo
        fname = root//'/metadata_roundtrip.txt'
        ok = manifest%write(fname, message)
        call assert_true(ok, 'fraction metadata: complete manifest can be published')
        if( ok )then
            found_digest = .false.
            found_count = .false.
            open(newunit=funit, file=fname%to_char(), status='old', action='read')
            do
                read(funit,'(A)',iostat=ios) line
                if( ios /= 0 ) exit
                if( index(line,'eo_digest') > 0 ) found_digest = index(line,'6C4DD423214CC38D') > 0
                if( index(line,'nptcls') > 0 )then
                    read(line,*,iostat=ios) prefix, key, count_value
                    if( ios == 0 ) found_count = count_value == 4
                endif
            enddo
            close(funit)
            call assert_true(found_digest, 'fraction metadata: stable ordered particle/half fingerprint')
            call assert_true(found_count, 'fraction metadata: serialized state-1 particle count')
        endif
        call volume%kill
    end subroutine check_manifest_partition

    !> Verify retained interval reconstructions and their exported raw FSCs and dose manifest.
    module subroutine run_all_refine_motion_model_lib_tests()
        use simple_test_flex_phantom_fixture, only: create_flex_pca_phantom_fixture
        type(string) :: cwd, root, truth_a, truth_b, truth_mean, truth_diff, oritab, stack
        integer :: nfail0
        nfail0 = tests_failed
        call enter_fixture('fraction_dose_export', cwd, root)
        call create_flex_pca_phantom_fixture(root, 1, truth_a, truth_b, truth_mean, truth_diff, oritab, stack, nptcls=RECON_NPTCLS)
        call check_fraction_dose_export(root, 'valid', stack, oritab, .true., .false.)
        call check_fraction_dose_export(root, 'missing_dose', stack, oritab, .false., .false.)
        call check_fraction_dose_export(root, 'changed_state', stack, oritab, .true., .true.)
        call simple_chdir(root)
        CWD_GLOB = root%to_char()
        call leave_fixture(cwd, root, nfail0)
        CWD_GLOB = cwd%to_char()
    end subroutine run_all_refine_motion_model_lib_tests

end submodule simple_refine_motion_model_tester
