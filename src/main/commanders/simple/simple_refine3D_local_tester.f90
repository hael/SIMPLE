!@descr: local-refinement startup rejection and artifact-provenance checks
module simple_refine3D_local_tester
use, intrinsic :: iso_fortran_env, only: int8
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use simple_string, only: string
use simple_syslib, only: simple_mkdir, simple_chdir, simple_list_dirs, file_exists
use simple_cmdline, only: cmdline
use simple_image, only: image
use simple_sp_project, only: sp_project
use simple_commanders_refine3D, only: commander_refine3D_local
use simple_test_dispatch, only: test_dispatcher, test_execution, EXEC_FAILED, proc_available
use simple_test_utils, only: assert_true, assert_false, assert_int, assert_char, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_refine3D_local_tests

integer, parameter :: BOX = 16
real, parameter :: SMPD = 1.5
character(len=24), parameter :: CASES(16) = [character(len=24) :: &
    &'missing_even', 'missing_odd', 'bad_dimensions', 'bad_sampling', 'bad_extent', 'bad_projection', &
    &'bad_angle', 'bad_shift', 'fractional_half', 'unassigned_half', 'one_half', 'no_active', &
    &'merged_lp', 'merged_nu', 'reference_halves', 'reference_state']
character(len=20), parameter :: INPUTS(4) = [character(len=20) :: &
    &'source.simple', 'reference.mrc', 'reference_even.mrc', 'reference_odd.mrc']
type :: file_snapshot
    integer(int8), allocatable :: bytes(:)
    logical :: exists = .false.
end type file_snapshot

contains

    ! Fatal checks run in isolated children; the parent checks exit, reason and source preservation.
    subroutine run_all_refine3D_local_tests()
        type(test_dispatcher) :: dispatcher
        type(test_execution) :: execution
        type(file_snapshot) :: before(4), after
        type(string) :: saved, root, case_dir, name
        type(string), allocatable :: dirs(:)
        integer :: icase, i, nfail0
        logical :: interrupted
#if !defined(__linux__)
        write(*,'(A)') 'SKIP: local startup rejection tests require the Linux test dispatcher'
        return
#endif
        if( .not. proc_available() )then
            write(*,'(A)') 'SKIP: local startup rejection tests require /proc'
            return
        endif
        nfail0 = tests_failed
        call enter_fixture('local_startup', saved, root)
        call test_suite_binding()
        do icase = 1, size(CASES)
            case_dir = root%to_char()//'/'//trim(CASES(icase))
            call simple_mkdir(case_dir)
            call simple_chdir(case_dir)
            call write_inputs(trim(CASES(icase)))
            do i = 1, size(INPUTS)
                call snapshot(string(trim(INPUTS(i))), before(i))
            enddo
            call simple_chdir(root)
            call dispatcher%new(1, 1)
            execution = test_execution()
            execution%prg = trim(CASES(icase))
            execution%label = trim(CASES(icase))
            execution%dir = case_dir%to_char()
            execution%deadline = 10
            call dispatcher%add(execution)
            call dispatcher%run(rejection_route, interrupted)
            execution = dispatcher%get_execution(1)
            call assert_false(interrupted, trim(CASES(icase))//': child was not interrupted')
            call assert_int(EXEC_FAILED, execution%state, trim(CASES(icase))//': input was rejected')
            call assert_true(execution%exit_code /= 0, trim(CASES(icase))//': rejection has nonzero exit')
            name = case_dir%to_char()//'/run.log'
            call assert_true(log_has(name, expected_reason(trim(CASES(icase)))), &
                &trim(CASES(icase))//': intended preflight diagnostic')
            dirs = simple_list_dirs(case_dir)
            if( allocated(dirs) )then
                call assert_int(0, size(dirs), trim(CASES(icase))//': no parameter-created run directory')
                deallocate(dirs)
            else
                call assert_true(.true., trim(CASES(icase))//': no parameter-created run directory')
            endif
            do i = 1, size(INPUTS)
                name = case_dir%to_char()//'/'//trim(INPUTS(i))
                call snapshot(name, after)
                call assert_true(after%exists .eqv. before(i)%exists, &
                    &trim(CASES(icase))//': input existence unchanged '//trim(INPUTS(i)))
                if( .not. before(i)%exists .or. .not. after%exists ) cycle
                call assert_int(size(before(i)%bytes), size(after%bytes), 'source byte count unchanged')
                if( size(before(i)%bytes) == size(after%bytes) )then
                    call assert_true(all(before(i)%bytes == after%bytes), &
                        &trim(CASES(icase))//': source bytes unchanged '//trim(INPUTS(i)))
                endif
            enddo
            call dispatcher%kill()
        enddo
        call test_provenance()
        call test_startup_staging()
        call leave_fixture(saved, root, nfail0)
    end subroutine run_all_refine3D_local_tests

    ! The real parameter constructor must retain the selector before the expensive simulation.
    subroutine test_suite_binding()
        use simple_parameters, only: parameters
        class(parameters), allocatable :: params
        class(cmdline), allocatable :: cline
        allocate(params, cline)
        call cline%set('prg', 'cont_refine3D_1jxy')
        call cline%set('mkdir', 'no')
        call cline%set('box', 32)
        call cline%set('smpd', SMPD)
        call cline%set('mskdiam', 12.)
        call cline%set('nthr', 1)
        call params%new(cline, silent=.true.)
        call assert_char('', trim(params%suite), 'suite parser: default selects the full gate')
        call cline%set('suite', 'nu_handoff')
        call params%new(cline, silent=.true.)
        call assert_char('nu_handoff', trim(params%suite), 'suite parser: focused NU handoff survives construction')
        call cline%set('suite', 'all')
        call params%new(cline, silent=.true.)
        call assert_char('all', trim(params%suite), 'suite parser: explicit full gate survives construction')
        call cline%kill()
    end subroutine test_suite_binding

    ! Small files exercise identity and interrupted publication without a refinement run.
    subroutine test_provenance()
        use simple_parameters, only: parameters
        use simple_fileio, only: arr2file
        use simple_syslib, only: del_file
        use simple_halfmap_diagnostics, only: support_provenance_fname
        use simple_refine3D_provenance, only: publish_reuse_provenance, reuse_provenance_matches, &
            &begin_reference_set, publish_reference_set, validate_reference_set, reference_matches_project, provenance_file, &
            &REFSET_OK, REFSET_INVALID, REFSET_LEGACY
        class(sp_project), allocatable :: project
        class(parameters), allocatable :: params
        class(image), allocatable :: volume
        type(string) :: sigma, chain, files(3), halves(2), merged, fsc, checked_fsc, record
        character(len=256) :: message
        integer :: i, status
        logical :: valid, support_verified
        allocate(project, params, volume)
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'provenance')
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'stk', 'particles.mrc')
        call project%os_ptcl3D%new(4, is_ptcl=.true.)
        do i = 1,4
            call project%os_ptcl3D%set(i, 'stkind', 1)
            call project%os_ptcl3D%set(i, 'indstk', i)
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'eo', mod(i,2))
            call project%os_ptcl3D%set(i, 'updatecnt', 1)
            call project%os_ptcl3D%set_euler(i, [real(i), 20., 30.])
            call project%os_ptcl3D%set_shift(i, [0., 0.])
        enddo
        params%oritype = 'ptcl3D'
        params%projfile = 'producer.simple'
        params%box_crop = BOX
        params%smpd_crop = SMPD
        params%mskdiam = 12.
        params%l_ml_reg = .false.
        params%l_nonuniform = .false.
        sigma = 'provenance_sigma.bin'
        call arr2file([1.,2.], sigma)
        call publish_reuse_provenance(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma])
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_true(valid, 'provenance: unchanged sigma producer identity')
        call project%os_ptcl3D%set(1, 'eo', 0)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed half assignment')
        call project%os_ptcl3D%set(1, 'eo', 1)
        call project%os_ptcl3D%set_state(1, 0)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed active population')
        call project%os_ptcl3D%set_state(1, 1)
        call project%os_ptcl3D%set(1, 'indstk', 2)
        call project%os_ptcl3D%set(2, 'indstk', 1)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: reordered particle rows')
        call project%os_ptcl3D%set(1, 'indstk', 1)
        call project%os_ptcl3D%set(2, 'indstk', 2)
        params%mskdiam = 14.
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed noise mask setting')
        params%mskdiam = 12.
        call project%os_ptcl3D%set(1, 'dfx', 1.)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_false(valid, 'provenance: changed observation CTF')
        call project%os_ptcl3D%set(1, 'dfx', 0.)
        call project%os_ptcl3D%set_euler(1, [2.,20.,30.])
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'sigma', sigma, [sigma], message)
        call assert_true(valid, 'provenance: sigma generation may lag poses')
        call project%os_ptcl3D%set_euler(1, [1.,20.,30.])
        chain = 'provenance_chain.txt'
        files(1) = chain
        files(2) = 'provenance_chain_even.bin'
        files(3) = 'provenance_chain_odd.bin'
        do i = 1,3
            call arr2file([real(i),0.], files(i))
        enddo
        call publish_reuse_provenance(params, project, project%os_ptcl3D, 'chain', chain, files)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_true(valid, 'provenance: unchanged chain producer identity')
        call project%os_ptcl3D%set(1, 'updatecnt', 0)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_false(valid, 'provenance: changed represented population')
        call project%os_ptcl3D%set(1, 'updatecnt', 1)
        call project%os_ptcl3D%set_euler(1, [2.,20.,30.])
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_false(valid, 'provenance: changed chain producer poses')
        call project%os_ptcl3D%set_euler(1, [1.,20.,30.])
        call arr2file([0.,1.], chain)
        valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', chain, files, message)
        call assert_false(valid, 'provenance: changed chain generation with equal byte count')
        merged = 'provenance_reference.mrc'
        halves(1) = 'provenance_reference_even.mrc'
        halves(2) = 'provenance_reference_odd.mrc'
        fsc = 'provenance_fsc.bin'
        call volume%new([BOX,BOX,BOX], SMPD)
        call volume%write(merged)
        call volume%write(halves(1))
        call volume%write(halves(2))
        call arr2file([1.,0.5], fsc)
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 1, merged, halves, fsc, 0.)
        status = validate_reference_set(merged, checked_fsc, support_verified, message)
        call assert_int(REFSET_OK, status, 'provenance: complete reference set')
        valid = reference_matches_project(merged, project, 1, .false., message)
        call assert_true(valid, 'reference: recorded continuation identity')
        call project%os_ptcl3D%set(1, 'eo', 0)
        valid = reference_matches_project(merged, project, 1, .true., message)
        call assert_false(valid, 'reference: explicit input cannot bypass same-lineage half changes')
        call project%os_ptcl3D%set(1, 'eo', 1)
        call project%os_ptcl3D%set_state(1, 0)
        valid = reference_matches_project(merged, project, 1, .false., message)
        call assert_false(valid, 'reference: changed active population')
        call project%os_ptcl3D%set_state(1, 1)
        valid = reference_matches_project(merged, project, 2, .false., message)
        call assert_false(valid, 'reference: different selected producer state')
        call project%os_ptcl3D%set_euler(1, [2.,20.,30.])
        valid = reference_matches_project(merged, project, 1, .false., message)
        call assert_true(valid, 'reference: starting poses may differ from lagged reference poses')
        call project%os_ptcl3D%set_euler(1, [1.,20.,30.])
        call project%os_ptcl3D%set(1, 'indstk', 2)
        call project%os_ptcl3D%set(2, 'indstk', 1)
        valid = reference_matches_project(merged, project, 1, .true., message)
        call assert_false(valid, 'reference: reordered rows are not an independent external lineage')
        call project%os_ptcl3D%set(1, 'indstk', 1)
        call project%os_ptcl3D%set(2, 'indstk', 2)
        call project%projinfo%set(1, 'projname', 'independent')
        valid = reference_matches_project(merged, project, 1, .false., message)
        call assert_false(valid, 'reference: registered input requires its particle lineage')
        valid = reference_matches_project(merged, project, 1, .true., message)
        call assert_true(valid, 'reference: independent explicit input keeps its trusted external contract')
        call project%projinfo%set(1, 'projname', 'provenance')
        record = support_provenance_fname(merged)
        call arr2file([1.], record)
        status = validate_reference_set(merged, checked_fsc, support_verified, message)
        call assert_int(REFSET_OK, status, 'provenance: extra sidecar does not alter the recorded base set')
        call assert_false(support_verified, 'provenance: unrecorded support sidecar is not trusted')
        call arr2file([0.5,1.], fsc)
        status = validate_reference_set(merged, checked_fsc, support_verified, message)
        call assert_int(REFSET_INVALID, status, 'provenance: stale FSC with equal byte count')
        call begin_reference_set(merged)
        status = validate_reference_set(merged, checked_fsc, support_verified, message)
        call assert_int(REFSET_INVALID, status, 'provenance: interrupted publication is not legacy')
        record = provenance_file(merged)//'.pending'
        call del_file(record)
        status = validate_reference_set(merged, checked_fsc, support_verified, message)
        call assert_int(REFSET_LEGACY, status, 'provenance: unrecorded legacy base pair')
        call test_reference_derivatives(params, project)
        call project%kill()
        call volume%kill()
    end subroutine test_provenance

    subroutine test_reference_derivatives( params, project )
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_parameters, only: parameters
        use simple_fileio, only: arr2file, file2rarr, add2fbody
        use simple_syslib, only: del_file
        use simple_defs_fname, only: MRC_EXT, NUFILT_SUFFIX, NULOCRES_SUFFIX, AUTOMASK_FBODY, NU_ENVMASK_FBODY
        use simple_refine3D_provenance, only: begin_reference_set, publish_reference_set, stage_reference_derivatives
        use simple_sigma2_state_file, only: sigma2_state_digest_file
        class(parameters), intent(inout) :: params
        class(sp_project), intent(in) :: project
        type(string) :: merged, halves(2), fsc, source(6), staged(6), target, name
        integer(int64) :: digests(6)
        real, allocatable :: values(:)
        real :: align_lp
        integer :: i
        logical :: reused
        character(len=256) :: message
        merged = 'nu_reference.mrc'
        halves = [string('nu_reference_even.mrc'), string('nu_reference_odd.mrc')]
        fsc = 'nu_fsc.bin'
        target = 'staged_reference.mrc'
        call arr2file([1.,0.], merged)
        do i = 1,2
            call arr2file([real(i),0.], halves(i))
        enddo
        call arr2file([1.,0.5], fsc)
        source(1) = add2fbody(merged, MRC_EXT, NUFILT_SUFFIX)
        source(2) = add2fbody(halves(1), MRC_EXT, NUFILT_SUFFIX)
        source(3) = add2fbody(halves(2), MRC_EXT, NUFILT_SUFFIX)
        source(4) = add2fbody(merged, MRC_EXT, NULOCRES_SUFFIX)
        source(5) = AUTOMASK_FBODY//'02'//MRC_EXT
        source(6) = NU_ENVMASK_FBODY//'02'//MRC_EXT
        staged(1) = add2fbody(target, MRC_EXT, NUFILT_SUFFIX)
        name = add2fbody(target, MRC_EXT, '_even')
        staged(2) = add2fbody(name, MRC_EXT, NUFILT_SUFFIX)
        name = add2fbody(target, MRC_EXT, '_odd')
        staged(3) = add2fbody(name, MRC_EXT, NUFILT_SUFFIX)
        staged(4) = add2fbody(target, MRC_EXT, NULOCRES_SUFFIX)
        staged(5) = AUTOMASK_FBODY//'01'//MRC_EXT
        staged(6) = NU_ENVMASK_FBODY//'01'//MRC_EXT
        do i = 1,6
            call arr2file([real(i),0.], source(i))
            digests(i) = sigma2_state_digest_file(source(i))
        enddo
        params%l_nonuniform = .true.
        params%filt_mode = 'nonuniform'
        params%automsk = 'nu'
        params%lp = 9.
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 8.)
        ! The live LP is an output of matching, not a fixed preparation setting.
        params%lp = 20.
        reused = stage_reference_derivatives(params, merged, target, align_lp, message)
        call assert_true(reused, 'NU bundle: verified derivatives are staged')
        call assert_true(abs(align_lp - 8.) < 1.e-6, 'NU bundle: recorded matching LP is restored')
        do i = 1,6
            call assert_true(file_exists(staged(i)), 'NU bundle: recorded derivative exists at destination')
            if( .not. file_exists(staged(i)) ) cycle
            values = file2rarr(staged(i))
            call assert_true(size(values) == 2, 'NU bundle: destination shape')
            if( size(values) == 2 ) call assert_true(all(values == [real(i),0.]), 'NU bundle: copied independent values')
            call assert_true(sigma2_state_digest_file(source(i)) == digests(i), 'NU bundle: source bytes preserved')
        enddo
        reused = stage_reference_derivatives(params, merged, merged, align_lp, message)
        call assert_true(reused, 'NU bundle: in-place map copies are accepted')
        do i = 1,6
            call assert_true(sigma2_state_digest_file(source(i)) == digests(i), 'NU bundle: in-place copies preserve bytes')
        enddo
        call assert_true(file_exists(staged(6)), 'NU bundle: old destination envelope exists before staging')
        call del_file(source(6))
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 8.)
        reused = stage_reference_derivatives(params, merged, target, align_lp, message)
        call assert_true(reused, 'NU bundle: missing recorded envelope still stages verified constituents')
        call assert_false(file_exists(staged(6)), 'NU bundle: unrecorded destination envelope is removed')
        call assert_true(file_exists(staged(5)), 'NU bundle: recorded destination mask is retained')
        call assert_true(sigma2_state_digest_file(source(5)) == digests(5), 'NU bundle: recorded source mask is preserved')
        params%amsklp = params%amsklp + 1.
        reused = stage_reference_derivatives(params, merged, string('incompatible.mrc'), align_lp, message)
        call assert_false(reused, 'NU bundle: changed preparation settings use the fallback')
        name = add2fbody(string('incompatible.mrc'), MRC_EXT, NUFILT_SUFFIX)
        call assert_false(file_exists(name), 'NU bundle: rejected derivatives are not copied')
        params%amsklp = params%amsklp - 1.
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 0.)
        reused = stage_reference_derivatives(params, merged, target, align_lp, message)
        call assert_false(reused, 'NU bundle: absent matching LP uses the fallback')
        params%l_nonuniform = .false.
        call begin_reference_set(merged)
        call publish_reference_set(params, project, project%os_ptcl3D, 2, merged, halves, fsc, 8.)
        params%l_nonuniform = .true.
        reused = stage_reference_derivatives(params, merged, target, align_lp, message)
        call assert_false(reused, 'NU bundle: unrecorded existing derivatives are not trusted')
        reused = stage_reference_derivatives(params, string('legacy.mrc'), target, align_lp, message)
        call assert_false(reused, 'NU bundle: absent provenance uses the fallback')
    end subroutine test_reference_derivatives

    ! Impulses have nonzero power in every shell: equal/opposite pairs have FSC +1/-1.
    subroutine test_startup_staging()
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_parameters, only: parameters
        use simple_fileio, only: file2rarr
        use simple_halfmap_diagnostics, only: write_base_halfmap_fsc
        use simple_syslib, only: simple_abspath
        use simple_refine3D_strategy, only: stage_refine3D_local_references
        use simple_sigma2_state_file, only: sigma2_state_digest_file
        class(parameters), allocatable :: params
        class(sp_project), allocatable :: project
        class(image), allocatable :: volume
        type(string) :: source(3), target, target_path, fsc_file, missing_fsc, registered
        integer(int64) :: digests(3)
        real, allocatable :: fsc(:)
        real :: registered_smpd
        integer :: i, registered_box
        logical :: reused
        allocate(params, project, volume)
        params%projfile = 'ownership_not_written.simple'
        params%l_nonuniform = .false.
        source = [string('ownership.mrc'), string('ownership_even.mrc'), string('ownership_odd.mrc')]
        fsc_file = 'ownership_fsc.bin'
        missing_fsc = 'ownership_missing_fsc.bin'
        call volume%new([BOX,BOX,BOX], SMPD)
        call volume%set_rmat_at(1, 1, 1, 1.)
        do i = 1,3
            call volume%write(source(i))
        enddo
        call write_base_halfmap_fsc(source(2), source(3), [BOX,BOX,BOX], SMPD, fsc_file)
        fsc = file2rarr(fsc_file)
        call assert_int(BOX/2, size(fsc), 'startup FSC: one value per nonzero shell')
        call assert_true(all(abs(fsc-1.) < 1.e-6), 'startup FSC: identical impulses give one in every shell')
        call volume%set_rmat_at(1, 1, 1, -1.)
        call volume%write(source(3))
        do i = 1,3
            digests(i) = sigma2_state_digest_file(source(i))
        enddo
        ! No project file exists: staging must use and update the supplied project in memory.
        call project%os_ptcl3D%new(1, is_ptcl=.true.)
        target = source(1)
        call stage_refine3D_local_references(params, project, target, missing_fsc, .false., reused)
        call assert_false(reused, 'startup staging: legacy pair takes the base-half fallback')
        call assert_true(project%isthere_in_osout('vol', 1), 'startup staging: caller receives the volume registration')
        call project%get_vol('vol', 1, registered, registered_smpd, registered_box)
        target_path = simple_abspath(target)
        call assert_char(target_path%to_char(), registered%to_char(), 'startup staging: caller volume is the staged volume')
        call assert_int(BOX, registered_box, 'startup staging: registered native box')
        call assert_true(abs(registered_smpd-SMPD) < 1.e-6, 'startup staging: registered native sampling')
        call project%get_fsc(1, fsc_file, registered_box)
        fsc = file2rarr(fsc_file)
        call assert_true(all(abs(fsc+1.) < 1.e-6), 'startup staging: missing FSC is rebuilt from opposite impulses')
        do i = 1,3
            call assert_true(sigma2_state_digest_file(source(i)) == digests(i), 'startup staging: source map bytes unchanged')
        enddo
        call assert_false(file_exists(params%projfile), 'startup staging: project publication remains with the caller')
        call test_chain_staging(params, project)
        call project%kill()
        call volume%kill()
    end subroutine test_startup_staging

    subroutine test_chain_staging( params, project )
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_parameters, only: parameters
        use simple_fileio, only: arr2file
        use simple_refine3D_strategy, only: stage_refine3D_local_chain
        use simple_refine3D_fnames, only: refine3D_trail_rec_fname, refine3D_trail_rho_fname, &
            &refine3D_trail_manifest_fname, refine3D_pcg_trail_accum_fname, refine3D_pcg_trail_manifest_fname
        use simple_trail_chain_manifest, only: trail_chain_manifest
        use simple_refine3D_provenance, only: publish_reuse_provenance, reuse_provenance_matches
        use simple_sigma2_state_file, only: sigma2_state_digest_file
        class(parameters), intent(inout) :: params
        class(sp_project), intent(inout) :: project
        type(trail_chain_manifest) :: manifest
        type(string) :: files(5), anchor
        integer(int64) :: sizes(4), digests(5)
        integer :: backend, ncomponents, i, status
        logical :: valid
        character(len=256) :: message
        params%l_update_frac = .true.
        params%l_trail_rec = .true.
        params%l_m_estimator_flex = .false.
        params%oritype = 'ptcl3D'
        params%box_crop = BOX
        params%smpd_crop = SMPD
        call project%projinfo%new(1, is_ptcl=.false.)
        call project%projinfo%set(1, 'projname', 'ownership')
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'stk', 'ownership_particles.mrc')
        call project%os_ptcl3D%set(1, 'stkind', 1)
        call project%os_ptcl3D%set(1, 'indstk', 1)
        call project%os_ptcl3D%set_state(1, 1)
        call project%os_ptcl3D%set(1, 'eo', 0)
        call project%os_ptcl3D%set(1, 'updatecnt', 1)
        call project%os_ptcl3D%set_euler(1, [0.,0.,0.])
        call project%os_ptcl3D%set_shift(1, [0.,0.])
        do backend = 1,2
            if( backend == 1 )then
                params%rec_backend = 'gridding'
                files(1) = refine3D_trail_manifest_fname(1)
                files(2) = refine3D_trail_rec_fname(1, 'even')
                files(3) = refine3D_trail_rho_fname(1, 'even')
                files(4) = refine3D_trail_rec_fname(1, 'odd')
                files(5) = refine3D_trail_rho_fname(1, 'odd')
                ncomponents = 4
            else
                params%rec_backend = 'pcg'
                files(1) = refine3D_pcg_trail_manifest_fname(1)
                files(2) = refine3D_pcg_trail_accum_fname(1, 'even')
                files(3) = refine3D_pcg_trail_accum_fname(1, 'odd')
                ncomponents = 2
            endif
            sizes = 0_int64
            do i = 1,ncomponents
                call arr2file([real(i),0.], files(i+1))
                inquire(file=files(i+1)%to_char(), size=sizes(i))
            enddo
            call manifest%new(BOX, SMPD, 1, 1, 1, 1, sizes, 1., 1, [0_int64,0_int64])
            call manifest%write(files(1), status)
            call assert_int(0, status, 'startup chain: fixture manifest published')
            call publish_reuse_provenance(params, project, project%os_ptcl3D, 'chain', files(1), files(:ncomponents+1))
            valid = reuse_provenance_matches(params, project, project%os_ptcl3D, 'chain', &
                &files(1), files(:ncomponents+1), message)
            call assert_true(valid, 'startup chain: fixture provenance is consumable: '//trim(message))
            if( .not. valid ) cycle
            do i = 1,ncomponents+1
                digests(i) = sigma2_state_digest_file(files(i))
            enddo
            anchor = 'ownership.mrc'
            valid = stage_refine3D_local_chain(params, project, anchor, BOX, SMPD)
            call assert_true(valid, 'startup chain: compatible in-place chain accepted')
            do i = 1,ncomponents+1
                call assert_true(sigma2_state_digest_file(files(i)) == digests(i), 'startup chain: in-place file bytes unchanged')
            enddo
            call project%os_ptcl3D%set_euler(1, [1.,0.,0.])
            valid = stage_refine3D_local_chain(params, project, anchor, BOX, SMPD)
            call assert_false(valid, 'startup chain: changed producer poses reject reuse')
            call project%os_ptcl3D%set_euler(1, [0.,0.,0.])
            call assert_false(file_exists(params%projfile), 'startup chain: no project file read or write needed')
        enddo
        call manifest%kill()
    end subroutine test_chain_staging

    ! The real public commander must reject before params%new, not a test-only copy of its checks.
    subroutine rejection_route( which, cline )
        character(len=*), intent(in) :: which
        class(cmdline), intent(inout) :: cline
        class(commander_refine3D_local), allocatable :: command
        allocate(command)
        call cline%delete('test')
        call cline%set('prg', 'refine3D_local')
        call cline%set('projfile', 'source.simple')
        call cline%set('vol1', 'reference.mrc')
        call cline%set('mskdiam', 12.)
        call cline%set('nthr', 1)
        call cline%set('maxits', 1)
        call cline%set('mkdir', 'yes')
        select case(which)
            case('merged_lp')
                call cline%set('filt_mode', 'none')
                call cline%set('lp', 8.)
            case('merged_nu')
                call cline%set('filt_mode', 'nonuniform_lpset')
        end select
        call command%execute(cline)
    end subroutine rejection_route

    subroutine write_inputs( which )
        use simple_parameters, only: parameters
        use simple_refine3D_provenance, only: begin_reference_set, publish_reference_set
        use simple_fileio, only: arr2file
        character(len=*), intent(in) :: which
        class(sp_project), allocatable :: project
        class(image), allocatable :: volume
        class(parameters), allocatable :: params
        type(string) :: halves(2)
        integer :: i, dims(3), producer_state
        real :: sampling, nan
        allocate(project, volume)
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'box', BOX)
        call project%os_stk%set(1, 'smpd', SMPD)
        call project%os_stk%set(1, 'stk', 'particles.mrc')
        call project%os_ptcl3D%new(2, is_ptcl=.true.)
        do i = 1,2
            call project%os_ptcl3D%set_euler(i, [0., 0., 0.])
            call project%os_ptcl3D%set_shift(i, [0., 0.])
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'proj', 1.)
            call project%os_ptcl3D%set(i, 'eo', real(i-1))
            call project%os_ptcl3D%set(i, 'stkind', 1)
            call project%os_ptcl3D%set(i, 'indstk', i)
        enddo
        nan = ieee_value(0., ieee_quiet_nan)
        select case(which)
            case('bad_projection')
                call project%os_ptcl3D%set(1, 'proj', 0.)
            case('bad_angle')
                call project%os_ptcl3D%set_euler(1, [nan, 0., 0.])
            case('bad_shift')
                call project%os_ptcl3D%set_shift(1, [nan, 0.])
            case('fractional_half')
                call project%os_ptcl3D%set(1, 'eo', 0.5)
            case('unassigned_half')
                call project%os_ptcl3D%set(1, 'eo', -1.)
            case('one_half')
                call project%os_ptcl3D%set(1, 'eo', 1.)
            case('no_active')
                call project%os_ptcl3D%set_all2single('state', 0.)
        end select
        call project%write(string('source.simple'))
        call project%kill()
        do i = 2, size(INPUTS)
            if( which == 'missing_even' .and. i == 3 ) cycle
            if( which == 'missing_odd' .and. i == 4 ) cycle
            dims = BOX
            sampling = SMPD
            if( which == 'bad_dimensions' .and. i == 3 ) dims = BOX/2
            if( which == 'bad_sampling' .and. i == 4 ) sampling = 2.
            if( which == 'bad_extent' ) sampling = 2.
            call volume%new(dims, sampling)
            call volume%write(string(trim(INPUTS(i))))
            call volume%kill()
        enddo
        if( which == 'reference_halves' .or. which == 'reference_state' )then
            allocate(params)
            call project%read(string('source.simple'))
            call project%update_projinfo(string('source.simple'))
            params%oritype = 'ptcl3D'
            params%projfile = 'source.simple'
            params%box_crop = BOX
            params%smpd_crop = SMPD
            params%l_ml_reg = .false.
            params%l_nonuniform = .false.
            halves = [string('reference_even.mrc'), string('reference_odd.mrc')]
            call arr2file([1.,0.5], string('reference_fsc.bin'))
            producer_state = 1
            if( which == 'reference_state' ) producer_state = 2
            call begin_reference_set(string('reference.mrc'))
            call publish_reference_set(params, project, project%os_ptcl3D, producer_state, &
                &string('reference.mrc'), halves, string('reference_fsc.bin'), 0.)
            if( which == 'reference_halves' )then
                call project%os_ptcl3D%set(1, 'eo', 1)
                call project%os_ptcl3D%set(2, 'eo', 0)
                call project%write(string('source.simple'))
            endif
            call project%kill()
        endif
    end subroutine write_inputs

    function expected_reason( which ) result(reason)
        character(len=*), intent(in) :: which
        character(len=96) :: reason
        select case(which)
            case('missing_even')
                reason = 'Continuous startup requires an existing reference: reference_even.mrc'
            case('missing_odd')
                reason = 'Continuous startup requires an existing reference: reference_odd.mrc'
            case('bad_dimensions', 'bad_sampling', 'bad_extent')
                reason = 'Continuous startup reference dimensions, sampling, or physical extent are incompatible'
            case('merged_lp')
                reason = 'refine3D_local requires half-specific references; use lpstop instead of lp'
            case('merged_nu')
                reason = 'refine3D_local does not support filt_mode=nonuniform_lpset'
            case('reference_halves')
                reason = 'Continuous startup reference provenance failed: reference half assignments'
            case('reference_state')
                reason = 'Continuous startup reference provenance failed: reference producer state'
            case('one_half')
                reason = 'Continuous startup requires particles in both half-sets'
            case('no_active')
                reason = 'refine=cont found no active particle'
            case default
                reason = 'refine=cont continues from the 3D poses'
        end select
    end function expected_reason

    logical function log_has( name, expected ) result(found)
        class(string), intent(in) :: name
        character(len=*), intent(in) :: expected
        character(len=2048) :: line
        integer :: unit, ios
        found = .false.
        open(newunit=unit, file=name%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(unit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, trim(expected)) > 0 ) found = .true.
        enddo
        close(unit)
    end function log_has

    subroutine snapshot( name, saved )
        class(string), intent(in) :: name
        type(file_snapshot), intent(inout) :: saved
        integer :: unit, nbytes
        if( allocated(saved%bytes) ) deallocate(saved%bytes)
        saved%exists = file_exists(name)
        if( .not. saved%exists ) return
        inquire(file=name%to_char(), size=nbytes)
        allocate(saved%bytes(nbytes))
        open(newunit=unit, file=name%to_char(), access='stream', form='unformatted', status='old', action='read')
        read(unit) saved%bytes
        close(unit)
    end subroutine snapshot

end module simple_refine3D_local_tester
