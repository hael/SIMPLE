!@descr: library tests of fractional reconstruction through reconstruct3D and the reconstruction service
! On FLEX's two-conformer phantom (simple_flex_pca_application_tester): fractional state weights that
! are not a function of the true conformer alone reconstruct each state as the weight mixture of the
! truth maps (same frame, 8 A low-pass), at least as well as the hard reconstruction of the same
! particles matches its truth, minus 0.01; 0/1 weights equal to the labels reproduce the hard maps
! byte for byte; gridding and PCG trailing seeds represent the same applied mass; a trailing chain of
! another weight set generation is re-seeded, one of the same is blended (both backends); under PCG a
! state left without particles is dropped and its previous maps are carried forward.
module simple_rec3D_service_tester
use simple_defs,                        only: logfhandle
use simple_string,                      only: string
use simple_string_utils,                only: int2str_pad
use simple_syslib,                      only: simple_mkdir, simple_chdir, simple_rename
use simple_fileio,                      only: files_identical, simple_copy_file
use simple_image,                       only: image
use simple_cmdline,                     only: cmdline
use simple_sp_project,                  only: sp_project
use simple_commanders_rec,              only: commander_rec3D
use simple_state_weight_set,            only: state_weight_set
use simple_trail_chain_manifest,        only: trail_chain_manifest, TRAIL_MANIFEST_OK
use simple_refine3D_fnames,             only: refine3D_state_vol_fname, refine3D_trail_manifest_fname, &
    &refine3D_pcg_trail_manifest_fname
use simple_test_truth_metrics,          only: compare_to_truth
use simple_flex_pca_application_tester, only: create_flex_pca_phantom_fixture, build_flex_pca_phantom_project, &
    &FLEX_PHANTOM_BOX, FLEX_PHANTOM_NPTCLS, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
use simple_error,                       only: simple_exception
use simple_test_utils,                  only: assert_false, assert_int, assert_true, enter_fixture, leave_fixture, tests_failed
implicit none
#include "simple_local_flags.inc"
private
public :: run_all_rec3D_service_lib_tests

integer, parameter :: NTHR       = 4
integer, parameter :: NSTATES    = 2
real,    parameter :: CORR_LP    = 8.0   !< same-frame low-pass of the map correlations (A)
! the margin of the flex-on-SIMPLE plan, measured weighted minus hard -0.0001 to +0.0004; MASS_RELTOL is not derived
real,    parameter :: CORR_MARGIN = 0.01
real,    parameter :: MASS_RELTOL = 1.e-5

type(string) :: root, truth(NSTATES), oritab, stack

contains

    subroutine run_all_rec3D_service_lib_tests()
        type(string) :: cwd_saved, truth_mean, truth_diff
        integer      :: nfail0
        write(*,'(A)') '**** running fractional reconstruction library tests ****'
        nfail0 = tests_failed
        call enter_fixture('rec3D_service', cwd_saved, root)
        call create_flex_pca_phantom_fixture(root, NTHR, truth(1), truth(2), truth_mean, truth_diff, oritab, stack)
        call test_labels_equal_hard('gridding')
        call test_unmasked_gridding_export()
        call test_labels_equal_hard('pcg')
        call test_fractional_phantom('gridding')
        call test_fractional_phantom('pcg')
        call test_trailing_weight_identity()
        call test_empty_state_carried()
        call test_request_settings()
        call simple_chdir(root)
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine run_all_rec3D_service_lib_tests

    !> conformer A (odd rows, state 1) or B (even rows, state 2) of the phantom
    pure integer function conformer( i )
        integer, intent(in) :: i
        conformer = merge(1, 2, mod(i,2) == 1)
    end function conformer

    !> the fixed fractional rule: mostly the particle's conformer, modulated by its index, so a
    !! state's expected map is a mixture of both truth maps rather than either of them
    subroutine fractional_weights( weights, labels )
        real,    intent(out) :: weights(FLEX_PHANTOM_NPTCLS,NSTATES)
        integer, intent(out) :: labels(FLEX_PHANTOM_NPTCLS)
        integer :: i
        do i = 1, FLEX_PHANTOM_NPTCLS
            weights(i,1) = merge(0.75, 0.30, conformer(i) == 1) + 0.1 * cos(0.37 * real(i))
            weights(i,2) = 1. - weights(i,1)
            labels(i)    = maxloc(weights(i,:), dim=1)
        enddo
    end subroutine fractional_weights

    !> a fresh phantom project in dir (entered), states set to labels
    subroutine make_project( dir, labels )
        character(len=*), intent(in) :: dir
        integer,          intent(in) :: labels(FLEX_PHANTOM_NPTCLS)
        class(sp_project), allocatable:: project
        integer :: i
        allocate(project)
        call simple_chdir(root)
        call simple_mkdir(string(dir))
        call build_flex_pca_phantom_project(root//'/'//dir//'/phantom.simple', stack, oritab)
        call simple_chdir(string(dir))
        call project%read(string('phantom.simple'))
        do i = 1, FLEX_PHANTOM_NPTCLS
            call project%os_ptcl3D%set_state(i, labels(i))
        enddo
        call project%write(string('phantom.simple'))
        call project%kill
    end subroutine make_project

    subroutine publish_weights( weights, labels )
        real,    intent(in) :: weights(FLEX_PHANTOM_NPTCLS,NSTATES)
        integer, intent(in) :: labels(FLEX_PHANTOM_NPTCLS)
        class(sp_project), allocatable:: project
        type(state_weight_set) :: wset
        integer :: i
        allocate(project)
        call project%read(string('phantom.simple'))
        call wset%publish(project, project%os_ptcl3D, string('phantom.simple'), &
            &[(i, i = 1, FLEX_PHANTOM_NPTCLS)], weights, labels, 'tester')
        call wset%kill
        call project%kill
    end subroutine publish_weights

    !> a request's settings (states, one part, no filter, every state) replace the caller's parameters for the
    !! reconstruction only: one weight column reconstructs one state under a caller of two states and state=2,
    !! whose parameters come back unchanged
    subroutine test_request_settings()
        use simple_rec3D_service, only: rec3D_service, rec3D_request, REC3D_WEIGHTS_TABLE, REC3D_OUTPUT_RAW, &
            &REC3D_DISPATCH_INPROC, REC3D_BACKEND_GRIDDING
        use simple_parameters,    only: parameters
        use simple_builder,       only: builder
        use simple_eo_partition,  only: ensure_field_halves
        use simple_fileio,        only: file_exists
        class(parameters),   allocatable :: params
        class(builder),      allocatable :: build
        class(cmdline),      allocatable :: cline, cline_rec
        type(rec3D_service), allocatable :: service
        type(rec3D_request) :: request
        integer :: truth_labels(FLEX_PHANTOM_NPTCLS), i
        write(*,'(A)') 'test_request_settings'
        truth_labels = [(conformer(i), i = 1, FLEX_PHANTOM_NPTCLS)]
        call make_project('request', truth_labels)
        allocate(params, build, cline, cline_rec, service)
        call cline%set('prg',         'reconstruct3D')
        call cline%set('projfile',    'phantom.simple')
        call cline%set('mkdir',       'no')
        call cline%set('nstates',     NSTATES)
        call cline%set('state',       2)
        call cline%set('mskdiam',     FLEX_PHANTOM_MSKDIAM)
        call cline%set('pgrp',        'c1')
        call cline%set('rec_backend', 'gridding')
        call cline%set('postprocess', 'no')
        call cline%set('oritype',     'ptcl3D')
        call cline%set('objfun',      'cc')
        call cline%set('nthr',        NTHR)
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
        call ensure_field_halves(build%spproj, params%oritype, params%projfile)
        ! the command line volume assembly rebuilds its parameters from, as the FLEX caller gives it
        cline_rec = cline
        call cline_rec%set('nstates',   1)
        call cline_rec%set('rho_floor', 'yes')
        call cline_rec%set('filt_mode', 'none')
        call cline_rec%set('automsk',   'no')
        call cline_rec%delete('state')
        request%pinds = [(i, i = 1, FLEX_PHANTOM_NPTCLS)]
        allocate(request%weights_table(FLEX_PHANTOM_NPTCLS,1), source=1.)
        request%weights       = REC3D_WEIGHTS_TABLE
        request%backend       = REC3D_BACKEND_GRIDDING
        request%output        = REC3D_OUTPUT_RAW
        request%prefix        = 'request'
        request%register      = .false.
        request%dispatch      = REC3D_DISPATCH_INPROC
        request%nstates       = 1
        request%one_part      = .true.
        request%no_nonuniform = .true.
        request%no_trail_rec  = .true.
        request%all_states    = .true.
        request%no_automsk    = .true.
        call service%new(params, build, cline_rec, REC3D_DISPATCH_INPROC, NTHR)
        call service%execute(params, build, cline_rec, request)
        call service%kill
        call assert_true(file_exists(string('request_state01_even.mrc')), 'the request''s one state is reconstructed')
        call assert_false(file_exists(string('request_state02_even.mrc')), 'no second state is reconstructed')
        call assert_int(NSTATES, params%nstates, 'the caller''s state count comes back')
        call assert_true(params%l_state_defined, 'the caller''s state selection comes back')
        call assert_int(2, params%state, 'the caller''s selected state comes back')
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call cline_rec%kill
        call cline%kill
    end subroutine test_request_settings

    !> reconstruct3D in this process and directory
    subroutine run_rec3D( backend, l_weighted, extra_key, extra_val )
        character(len=*),           intent(in) :: backend
        logical,                    intent(in) :: l_weighted
        character(len=*), optional, intent(in) :: extra_key(:), extra_val(:)
        type(commander_rec3D) :: xrec3D
        class(cmdline), allocatable:: cline
        integer :: i
        allocate(cline)
        call cline%set('prg',         'reconstruct3D')
        call cline%set('projfile',    'phantom.simple')
        call cline%set('mkdir',       'no')
        call cline%set('nstates',     NSTATES)
        call cline%set('mskdiam',     FLEX_PHANTOM_MSKDIAM)
        call cline%set('pgrp',        'c1')
        call cline%set('rec_backend', backend)
        call cline%set('postprocess', 'no')
        call cline%set('nthr',        NTHR)
        if( l_weighted ) call cline%set('m_estimator', 'flex')
        if( present(extra_key) )then
            do i = 1, size(extra_key)
                call cline%set_from_text(trim(extra_key(i)), trim(extra_val(i)))
            enddo
        endif
        call xrec3D%execute(cline)
        call cline%kill
    end subroutine run_rec3D

    !> Raw export preserves the standard maps/FSC while retaining density beyond their spherical support.
    subroutine test_unmasked_gridding_export()
        use simple_defs, only: COSMSKHALFWIDTH
        use simple_parameters, only: parameters
        use simple_refine3D_fnames, only: refine3D_fsc_fname
        use simple_halfmap_diagnostics, only: read_support_provenance
        use simple_test_utils, only: assert_char, assert_real
        character(len=5), parameter :: SUFFIXES(3) = [character(len=5) :: '', '_even', '_odd']
        character(len=16), parameter :: KEYS(4) = [character(len=16) :: 'objfun', 'ml_reg', 'automsk', 'msk_crop']
        character(len=16), parameter :: VALS(4) = [character(len=16) :: 'cc', 'no', 'no', '8']
        real, parameter :: MASK_RADIUS = 8., REL_TOL = 2.e-6
        class(parameters), allocatable :: defaults
        class(image), allocatable :: raw, standard
        real, allocatable :: rraw(:,:,:), rstandard(:,:,:)
        type(string) :: fname, rawname, baseline
        character(len=32) :: support
        integer :: labels(FLEX_PHANTOM_NPTCLS), i, s, h, x, y, z, c
        real :: radius2, raw_outside, standard_outside, scale, err
        logical :: constrained, found
        allocate(defaults, raw, standard)
        call assert_char('no', defaults%unfil_unmasked, 'raw export is disabled by default')
        labels = [(conformer(i), i=1,FLEX_PHANTOM_NPTCLS)]
        call make_project('unmasked_gridding', labels)
        call run_rec3D('gridding', .false., KEYS, VALS)
        do s = 1, NSTATES
            do h = 1, size(SUFFIXES)
                fname = state_map_fname(s, trim(SUFFIXES(h)))
                baseline = fname//'.baseline'
                call simple_copy_file(fname, baseline)
                if( h > 1 )then
                    rawname = state_map_fname(s, trim(SUFFIXES(h))//'_unfil')
                    call assert_true(files_identical(fname, rawname), 'default raw half retains spherical support')
                endif
            enddo
            fname = refine3D_fsc_fname(s)
            call simple_copy_file(fname, fname//'.baseline')
        enddo
        call run_rec3D('gridding', .false., [character(len=16) :: KEYS, 'unfil_unmasked'], &
            &[character(len=16) :: VALS, 'yes'])
        call raw%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD, wthreads=.false.)
        call standard%new([FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX,FLEX_PHANTOM_BOX], FLEX_PHANTOM_SMPD, wthreads=.false.)
        c = FLEX_PHANTOM_BOX / 2 + 1
        do s = 1, NSTATES
            do h = 1, size(SUFFIXES)
                fname = state_map_fname(s, trim(SUFFIXES(h)))
                call assert_true(files_identical(fname, fname//'.baseline'), 'raw export leaves standard map unchanged')
                if( h == 1 ) cycle
                rawname = state_map_fname(s, trim(SUFFIXES(h))//'_unfil')
                call read_support_provenance(rawname, constrained, found, support_kind=support)
                call assert_true(found, 'unmasked half has support provenance')
                call assert_false(constrained, 'unmasked half has no envelope constraint')
                call assert_char('none', trim(support), 'unmasked half records none support')
                call raw%read(rawname)
                call standard%read(fname)
                rraw = raw%get_rmat()
                rstandard = standard%get_rmat()
                scale = max(maxval(abs(rraw)), tiny(1.))
                raw_outside = 0.
                standard_outside = 0.
                do z = 1, FLEX_PHANTOM_BOX
                    do y = 1, FLEX_PHANTOM_BOX
                        do x = 1, FLEX_PHANTOM_BOX
                            radius2 = real((x-c)**2 + (y-c)**2 + (z-c)**2)
                            if( radius2 <= (MASK_RADIUS + COSMSKHALFWIDTH)**2 ) cycle
                            raw_outside = max(raw_outside, abs(rraw(x,y,z)))
                            standard_outside = max(standard_outside, abs(rstandard(x,y,z)))
                        enddo
                    enddo
                enddo
                call assert_true(raw_outside > REL_TOL * scale, 'unmasked half retains density outside the sphere')
                call assert_real(0., standard_outside, REL_TOL * scale, 'standard half is zero outside the sphere')
                call raw%mask3D_soft(MASK_RADIUS, backgr=0.)
                rraw = raw%get_rmat()
                err = maxval(abs(rraw-rstandard)) / scale
                call assert_real(0., err, REL_TOL, 'masking the raw export reproduces its standard half')
            enddo
            fname = refine3D_fsc_fname(s)
            call assert_true(files_identical(fname, fname//'.baseline'), 'raw export leaves FSC unchanged')
        enddo
        ! Reusing the names in default mode must replace the explicit unmasked provenance.
        call run_rec3D('gridding', .false., KEYS, VALS)
        do s = 1, NSTATES
            do h = 2, size(SUFFIXES)
                fname = state_map_fname(s, trim(SUFFIXES(h)))
                rawname = state_map_fname(s, trim(SUFFIXES(h))//'_unfil')
                call assert_true(files_identical(fname, rawname), 'default rerun restores masked raw exports')
                call read_support_provenance(rawname, constrained, found, support_kind=support)
                call assert_true(found, 'rewritten half retains support provenance')
                call assert_char('sphere', trim(support), 'default rerun replaces stale none support')
            enddo
        enddo
        call raw%kill
        call standard%kill
        call fname%kill
        call rawname%kill
        call baseline%kill
    end subroutine test_unmasked_gridding_export

    !> 0/1 weights equal to the hard labels give the hard run's state maps and half maps byte for byte
    subroutine test_labels_equal_hard( backend )
        character(len=*), intent(in) :: backend
        character(len=5), parameter :: SUFFIXES(3) = [character(len=5) :: '', '_even', '_odd']
        real    :: weights(FLEX_PHANTOM_NPTCLS,NSTATES)
        integer :: truth_labels(FLEX_PHANTOM_NPTCLS), i, s, k
        type(string) :: fname, fname_hard
        logical :: l_equal
        write(*,'(A)') 'test_labels_equal_hard '//backend
        truth_labels = [(conformer(i), i = 1, FLEX_PHANTOM_NPTCLS)]
        call make_project('eq_'//backend, truth_labels)
        call run_rec3D(backend, .false.)
        do s = 1, NSTATES
            do k = 1, size(SUFFIXES)
                fname      = state_map_fname(s, trim(SUFFIXES(k)))
                fname_hard = fname//'.hard'
                call simple_rename(fname, fname_hard)
            enddo
        enddo
        weights = 0.
        do i = 1, FLEX_PHANTOM_NPTCLS
            weights(i,truth_labels(i)) = 1.
        enddo
        call publish_weights(weights, truth_labels)
        call run_rec3D(backend, .true.)
        l_equal = .true.
        do s = 1, NSTATES
            do k = 1, size(SUFFIXES)
                fname      = state_map_fname(s, trim(SUFFIXES(k)))
                fname_hard = fname//'.hard'
                l_equal = l_equal .and. files_identical(fname, fname_hard)
            enddo
        enddo
        call assert_true(l_equal, backend//': 0/1 weights equal to the labels reproduce the hard maps bit for bit')
        call fname%kill
        call fname_hard%kill
    end subroutine test_labels_equal_hard

    !> PCG with a state that lost every particle while the other has them (the gridding rule): the
    !! run completes, the empty state's previous maps stay byte for byte, the other state is rewritten
    subroutine test_empty_state_carried()
        character(len=5), parameter :: SUFFIXES(3) = [character(len=5) :: '', '_even', '_odd']
        class(sp_project), allocatable:: project
        type(string)     :: fname, fname_prev
        integer :: truth_labels(FLEX_PHANTOM_NPTCLS), i, k, status
        logical :: l_kept
        allocate(project)
        write(*,'(A)') 'test_empty_state_carried pcg'
        truth_labels = [(conformer(i), i = 1, FLEX_PHANTOM_NPTCLS)]
        call make_project('empty_pcg', truth_labels)
        call run_rec3D('pcg', .false.)
        do k = 1, size(SUFFIXES)
            fname      = state_map_fname(2, trim(SUFFIXES(k)))
            fname_prev = fname//'.prev'
            call simple_copy_file(fname, fname_prev, status)
            if( status /= 0 ) THROW_HARD('could not keep '//fname%to_char())
        enddo
        fname      = state_map_fname(1, '')
        fname_prev = fname//'.prev'
        call simple_copy_file(fname, fname_prev, status)
        if( status /= 0 ) THROW_HARD('could not keep '//fname%to_char())
        ! every particle to state 1: state 2 has none
        call project%read(string('phantom.simple'))
        do i = 1, FLEX_PHANTOM_NPTCLS
            call project%os_ptcl3D%set_state(i, 1)
        enddo
        call project%write(string('phantom.simple'))
        call project%kill
        call run_rec3D('pcg', .false.)
        l_kept = .true.
        do k = 1, size(SUFFIXES)
            fname      = state_map_fname(2, trim(SUFFIXES(k)))
            fname_prev = fname//'.prev'
            l_kept = l_kept .and. files_identical(fname, fname_prev)
        enddo
        call assert_true(l_kept, 'pcg: the empty state keeps its previous maps byte for byte')
        fname      = state_map_fname(1, '')
        fname_prev = fname//'.prev'
        call assert_false(files_identical(fname, fname_prev), 'pcg: the populated state is reconstructed again')
        call fname%kill
        call fname_prev%kill
    end subroutine test_empty_state_carried

    function state_map_fname( s, suffix ) result( fname )
        integer,          intent(in) :: s
        character(len=*), intent(in) :: suffix
        type(string) :: fname
        fname = string('recvol_state')//int2str_pad(s,2)//suffix//'.mrc'
    end function state_map_fname

    !> per state: the weighted map against its expected mixture, the hard map against its truth
    subroutine test_fractional_phantom( backend )
        character(len=*), intent(in) :: backend
        real    :: weights(FLEX_PHANTOM_NPTCLS,NSTATES), corr_frac(NSTATES), corr_hard(NSTATES)
        real    :: masses(2,NSTATES), f05, f0143
        integer :: labels(FLEX_PHANTOM_NPTCLS), truth_labels(FLEX_PHANTOM_NPTCLS), i, s
        type(string) :: mixture
        write(*,'(A)') 'test_fractional_phantom '//backend
        call fractional_weights(weights, labels)
        truth_labels = [(conformer(i), i = 1, FLEX_PHANTOM_NPTCLS)]
        ! hard: the same particles labelled by their conformer
        call make_project('hard_'//backend, truth_labels)
        call run_rec3D(backend, .false.)
        do s = 1, NSTATES
            call compare_to_truth(truth(s), refine3D_state_vol_fname(s), FLEX_PHANTOM_MSKDIAM, corr_hard(s), &
                &f05, f0143, corr_lp=CORR_LP)
        enddo
        ! fractional
        call make_project('frac_'//backend, labels)
        call publish_weights(weights, labels)
        call run_rec3D(backend, .true.)
        masses = 0.
        do i = 1, FLEX_PHANTOM_NPTCLS
            masses(conformer(i),:) = masses(conformer(i),:) + weights(i,:)
        enddo
        do s = 1, NSTATES
            mixture = string('mixture_state')//int2str_pad(s,2)//'.mrc'
            call write_mixture(mixture, masses(:,s))
            call compare_to_truth(mixture, refine3D_state_vol_fname(s), FLEX_PHANTOM_MSKDIAM, corr_frac(s), &
                &f05, f0143, corr_lp=CORR_LP)
            write(logfhandle,'(A,A,A,I0,A,2F9.1,A,F7.4,A,F7.4)') '>>> FRACTIONAL PHANTOM ', backend, ' STATE ', s, &
                &' MASS A/B', masses(:,s), ' CORR(WEIGHTED, MIXTURE)=', corr_frac(s), ' CORR(HARD, TRUTH)=', corr_hard(s)
            call assert_true(corr_frac(s) >= corr_hard(s) - CORR_MARGIN, backend// &
                &': a weighted state map matches its expected mixture as well as the hard map its truth')
        enddo
        call mixture%kill
    end subroutine test_fractional_phantom

    !> the expected map of a state: the truth maps mixed by the state's mass from each conformer
    subroutine write_mixture( fname, masses )
        class(string), intent(in) :: fname
        real,          intent(in) :: masses(2)
        type(image) :: a, b
        integer :: ldim(3)
        ldim = [FLEX_PHANTOM_BOX, FLEX_PHANTOM_BOX, FLEX_PHANTOM_BOX]
        call a%new(ldim, FLEX_PHANTOM_SMPD, wthreads=.false.)
        call b%new(ldim, FLEX_PHANTOM_SMPD, wthreads=.false.)
        call a%read(truth(1))
        call b%read(truth(2))
        call a%mul(masses(1) / sum(masses))
        call b%mul(masses(2) / sum(masses))
        call a%add(b)
        call a%write(fname)
        call a%kill
        call b%kill
    end subroutine write_mixture

    !> trailing chains under fractional weights: the two backends' seeds represent the same applied mass,
    !! and a chain of another weight set generation is re-seeded (generation 1) while one of the same
    !! generation is blended (generation 2)
    subroutine test_trailing_weight_identity()
        type(trail_chain_manifest) :: man
        real    :: weights(FLEX_PHANTOM_NPTCLS,NSTATES), mass_grid, mass_pcg, mass_expected
        integer :: labels(FLEX_PHANTOM_NPTCLS), status
        integer(kind=8) :: id_first(2), id_second(2)
        write(*,'(A)') 'test_trailing_weight_identity'
        call fractional_weights(weights, labels)
        mass_expected = sum(weights(:,1))
        ! PCG seed (the weights exceed the PCG membership threshold, so both backends see every row)
        call make_project('trail_pcg', labels)
        call publish_weights(weights, labels)
        call run_rec3D('pcg', .true., [character(len=16) :: 'trail_seed'], [character(len=16) :: 'yes'])
        call man%read(refine3D_pcg_trail_manifest_fname(1), status)
        call assert_int(TRAIL_MANIFEST_OK, status, 'the PCG trailing seed publishes a manifest')
        mass_pcg = man%get_mrep()
        id_first = man%get_wset_id()
        call fractional_update_then_new_generation(weights, labels)
        call run_rec3D('pcg', .true., [character(len=16) :: 'update_frac', 'trail_rec'], [character(len=16) :: '0.5', 'yes'])
        call man%read(refine3D_pcg_trail_manifest_fname(1), status)
        id_second = man%get_wset_id()
        call assert_true(any(id_second /= id_first), 'PCG: the new weight set generation is recorded with the chain')
        call assert_int(1, man%get_gen(), 'PCG: a chain of another weight set generation is re-seeded, never blended')
        call run_rec3D('pcg', .true., [character(len=16) :: 'update_frac', 'trail_rec'], [character(len=16) :: '0.5', 'yes'])
        call man%read(refine3D_pcg_trail_manifest_fname(1), status)
        call assert_int(2, man%get_gen(), 'PCG: a chain of the same weight set generation is blended')
        ! gridding seed
        call make_project('trail_grid', labels)
        call publish_weights(weights, labels)
        call run_rec3D('gridding', .true., [character(len=16) :: 'trail_seed'], [character(len=16) :: 'yes'])
        call man%read(refine3D_trail_manifest_fname(1), status)
        call assert_int(TRAIL_MANIFEST_OK, status, 'the gridding trailing seed publishes a manifest')
        mass_grid = man%get_mrep()
        id_first  = man%get_wset_id()
        write(logfhandle,'(A,3F12.3)') '>>> TRAILING SEED MASS EXPECTED/GRIDDING/PCG ', mass_expected, mass_grid, mass_pcg
        call assert_true(abs(mass_grid - mass_expected) <= MASS_RELTOL * mass_expected, &
            &'the gridding chain represents the applied mass of the state')
        call assert_true(abs(mass_pcg - mass_grid) <= MASS_RELTOL * mass_grid, &
            &'PCG and gridding chains report the same applied mass for the same input')
        ! a fractional update of half the rows under a NEW weight set generation: re-seeded, not blended
        call fractional_update_then_new_generation(weights, labels)
        call run_rec3D('gridding', .true., [character(len=16) :: 'update_frac', 'trail_rec'], &
            &[character(len=16) :: '0.5', 'yes'])
        call man%read(refine3D_trail_manifest_fname(1), status)
        id_second = man%get_wset_id()
        call assert_true(any(id_second /= id_first), 'the new weight set generation is recorded with the chain')
        call assert_int(1, man%get_gen(), 'a chain of another weight set generation is re-seeded, never blended')
        ! the same update under the same generation: blended
        call run_rec3D('gridding', .true., [character(len=16) :: 'update_frac', 'trail_rec'], &
            &[character(len=16) :: '0.5', 'yes'])
        call man%read(refine3D_trail_manifest_fname(1), status)
        call assert_int(2, man%get_gen(), 'a chain of the same weight set generation is blended')
        call assert_true(all(man%get_wset_id() == id_second), 'the blended chain keeps its weight set identity')
        call man%kill
    end subroutine test_trailing_weight_identity

    !> every row updated once, half of them in the current sample; then the weights republished, so the
    !! set's generation (its identity) changes
    subroutine fractional_update_then_new_generation( weights, labels )
        real,    intent(in) :: weights(FLEX_PHANTOM_NPTCLS,NSTATES)
        integer, intent(in) :: labels(FLEX_PHANTOM_NPTCLS)
        class(sp_project), allocatable:: project
        integer :: i
        allocate(project)
        call project%read(string('phantom.simple'))
        do i = 1, FLEX_PHANTOM_NPTCLS
            call project%os_ptcl3D%set(i, 'updatecnt', 1.)
            call project%os_ptcl3D%set(i, 'sampled', merge(1., 0., mod(i,4) < 2))
        enddo
        call project%write(string('phantom.simple'))
        call project%kill
        call publish_weights(weights, labels)
    end subroutine fractional_update_then_new_generation

end module simple_rec3D_service_tester
