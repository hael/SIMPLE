!@descr: tests of fraction dose metadata and the unmasked reconstruction-to-manifest workflow
module simple_refine_motion_model_tester
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
use simple_defs, only: dp, CWD_GLOB
use simple_type_defs, only: OBJFUN_CC
use simple_string, only: string
use simple_string_utils, only: int2str, int2str_pad
use simple_syslib, only: simple_chdir, simple_mkdir
use simple_fileio, only: file_exists, simple_copy_file, file2rarr
use simple_parameters, only: parameters
use simple_cmdline, only: cmdline
use simple_sp_project, only: sp_project
use simple_motion_model, only: motion_model
use simple_image, only: image
use simple_dose_series_manifest, only: dose_series_manifest
use simple_refine_motion_model_strategy, only: prepare_fraction_dose_manifest, reconstruct_motion_fractions
use simple_test_utils, only: assert_true, assert_false, assert_int, assert_char, &
    &enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_refine_motion_model_tests, run_all_refine_motion_model_lib_tests
public :: check_fraction_dose_export

integer, parameter :: METADATA_BOX = 8, MOVIE_NFRAMES = 12, RECON_NPTCLS = 120
real, parameter :: METADATA_SMPD = 2., SAVED_FRACTION_DOSE = 1.25

contains

    !> Exercise metadata from the actual selected state-1 particles and saved decoded-fraction models.
    subroutine run_all_refine_motion_model_tests()
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

    !> Persist realistic model metadata with an EER grouping and target distinct from the saved fraction dose.
    subroutine write_saved_model(params, fname, nframes, dose, kv)
        class(parameters), target, intent(inout) :: params
        type(string), intent(in) :: fname
        integer, intent(in) :: nframes
        real, intent(in) :: dose, kv
        class(motion_model), allocatable :: model
        class(image), pointer :: frames(:)
        type(string) :: movie
        integer :: i, funit
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

    !> Verify retained interval reconstructions and their exported raw FSCs and dose manifest.
    subroutine run_all_refine_motion_model_lib_tests()
        use simple_flex_pca_application_tester, only: create_flex_pca_phantom_fixture
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

    !> Verify success and metadata-only rejection through the production reconstruction loop.
    subroutine check_fraction_dose_export(root, name, stack, oritab, valid_dose, changed_state)
        use simple_flex_pca_application_tester, only: build_flex_pca_phantom_project, &
            &FLEX_PHANTOM_BOX, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
        type(string), intent(in) :: root, stack, oritab
        character(len=*), intent(in) :: name
        logical, intent(in) :: valid_dose, changed_state
        class(parameters), allocatable :: params
        class(cmdline), allocatable :: cline
        class(sp_project), allocatable :: project
        class(dose_series_manifest), allocatable :: manifest
        type(string), allocatable :: projects(:), stacks(:)
        integer, allocatable :: ranges(:,:)
        type(string) :: arm, directory, model, message, even, odd
        real, allocatable :: fsc(:)
        real(dp) :: dlo, dhi
        integer :: i, j, nwin, nstates, a, b
        logical :: ok
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
        call reconstruct_motion_fractions(params, cline, projects, ranges, nstates)
        call assert_true(file_exists(arm//'/fsc_fractions_state01.tsv'), 'fraction export '//name//': raw FSC matrix retained')
        do i = 1, nwin
            call assert_true(file_exists(stacks(i)), 'fraction export '//name//': extracted stack retained')
            directory = arm//'/rec_frames_'//int2str_pad(ranges(1,i),2)//'-'//int2str_pad(ranges(2,i),2)
            call assert_true(file_exists(directory//'/recvol_state01_even_unfil.mrc'), 'fraction export '//name//': raw even half')
            call assert_true(file_exists(directory//'/recvol_state01_odd_unfil.mrc'), 'fraction export '//name//': raw odd half')
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
        call manifest%kill
        call simple_chdir(root)
        CWD_GLOB = root%to_char()
    end subroutine check_fraction_dose_export

end module simple_refine_motion_model_tester
