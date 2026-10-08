!@descr: library tests of fractional reconstruction through reconstruct3D and the reconstruction service
! On FLEX's two-conformer phantom (simple_flex_pca_application_tester): fractional state weights that
! are not a function of the true conformer alone reconstruct each state as the weight mixture of the
! truth maps (same frame, 8 A low-pass), at least as well as the hard reconstruction of the same
! particles matches its truth, minus 0.01; 0/1 weights equal to the labels reproduce the hard maps
! byte for byte; gridding and PCG trailing seeds represent the same applied mass; a trailing chain of
! another weight set generation is re-seeded, one of the same is blended (both backends).
module simple_rec3D_service_tester
use simple_defs,                        only: logfhandle
use simple_string,                      only: string
use simple_string_utils,                only: int2str_pad
use simple_syslib,                      only: simple_mkdir, simple_chdir, simple_rename
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
use simple_test_utils
implicit none
private
public :: run_all_rec3D_service_lib_tests

integer, parameter :: NTHR       = 4
integer, parameter :: NSTATES    = 2
real,    parameter :: CORR_LP    = 8.0   !< same-frame low-pass of the map correlations (A)
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
        call test_labels_equal_hard('pcg')
        call test_fractional_phantom('gridding')
        call test_fractional_phantom('pcg')
        call test_trailing_weight_identity()
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
        type(sp_project) :: project
        integer :: i
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
        type(sp_project)       :: project
        type(state_weight_set) :: wset
        integer :: i
        call project%read(string('phantom.simple'))
        call wset%publish(project, project%os_ptcl3D, string('phantom.simple'), &
            &[(i, i = 1, FLEX_PHANTOM_NPTCLS)], weights, labels, 'tester')
        call wset%kill
        call project%kill
    end subroutine publish_weights

    !> reconstruct3D in this process and directory
    subroutine run_rec3D( backend, l_weighted, extra_key, extra_val )
        character(len=*),           intent(in) :: backend
        logical,                    intent(in) :: l_weighted
        character(len=*), optional, intent(in) :: extra_key(:), extra_val(:)
        type(commander_rec3D) :: xrec3D
        type(cmdline)         :: cline
        integer :: i
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
                call cline%set(trim(extra_key(i)), trim(extra_val(i)))
            enddo
        endif
        call xrec3D%execute(cline)
        call cline%kill
    end subroutine run_rec3D

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
                l_equal = l_equal .and. same_bytes(fname, fname_hard)
            enddo
        enddo
        call assert_true(l_equal, backend//': 0/1 weights equal to the labels reproduce the hard maps bit for bit')
        call fname%kill
        call fname_hard%kill
    end subroutine test_labels_equal_hard

    function state_map_fname( s, suffix ) result( fname )
        integer,          intent(in) :: s
        character(len=*), intent(in) :: suffix
        type(string) :: fname
        fname = string('recvol_state')//int2str_pad(s,2)//suffix//'.mrc'
    end function state_map_fname

    logical function same_bytes( fname1, fname2 ) result( l_same )
        class(string), intent(in) :: fname1, fname2
        character(len=1), allocatable :: b1(:), b2(:)
        integer :: u1, u2, n1, n2
        l_same = .false.
        inquire(file=fname1%to_char(), size=n1)
        inquire(file=fname2%to_char(), size=n2)
        if( n1 <= 0 .or. n1 /= n2 ) return
        allocate(b1(n1), b2(n2))
        open(newunit=u1, file=fname1%to_char(), access='stream', form='unformatted', status='old', action='read')
        open(newunit=u2, file=fname2%to_char(), access='stream', form='unformatted', status='old', action='read')
        read(u1) b1
        read(u2) b2
        close(u1)
        close(u2)
        l_same = all(b1 == b2)
        deallocate(b1, b2)
    end function same_bytes

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
        type(sp_project) :: project
        integer :: i
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
