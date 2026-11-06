!@descr: unit tests for the continuous in-plane state of the refine3D search (simple_strategy3D_srch, _alloc, _utils)
! The search-state contract behind inpl_cont: the grid seed of a continuous candidate, candidate
! storage with and without the continuous arrays, the joint and discrete-seed routes, the invalid
! joint-evaluation predicate, resolve_inplane_e3, the inpl_cont policy (internal, default yes, no
! under refine=cont, never a command-line key) and the stored shift as the search's previous
! shift (N10; a small phantom written to the run directory and removed).
module simple_strategy3D_inplane_tester
use simple_core_module_api,   only: sp, dp
use simple_rnd,               only: gasdev
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_string,            only: string
use simple_linked_list,       only: list_iterator
use simple_ui,                only: make_ui, get_prg_ptr
use simple_ui_program,        only: ui_program, ui_program_input
use simple_refine3D_strategy, only: strip_refine3D_search_only_args
use simple_strategy3D_alloc,  only: clean_strategy3D, prep_strategy3D, s3D, seed_continuous_inplane_candidate
use simple_strategy3D_srch,   only: strategy3D_srch, strategy3D_spec
use simple_polarft_calc,      only: SHIFTMAP_WINDOW
use simple_builder,           only: builder
use simple_image,             only: image
use simple_syslib,            only: del_file
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_strategy3D_utils,  only: resolve_inplane_e3, candidate_pose
use simple_pftc_shsrch_grad,  only: pftc_shsrch_grad, shift_ref2ptcl
use simple_matcher_2Dprep,    only: prepimg4align
use simple_matcher_ptcl_io,   only: prep_rec_observation, prep_rec_observation_ptcl, pass_team, killimgbatch, &
                                    &rec_state_half
use simple_memoize_ft_maps,   only: memoize_ft_maps
use simple_matcher_ptcl_batch, only: check_pass_team, alloc_ptcl_imgs
use simple_forked_process,    only: forked_process, FORK_STATUS_FAILED
use simple_defs,              only: nthr_glob
use simple_test_utils
!$ use omp_lib,               only: omp_get_max_threads, omp_set_num_threads
implicit none
private
public :: run_all_strategy3D_inplane_tests

!> a pass on a team larger than nthr_glob, which must stop
type, extends(forked_process) :: oversized_team_pass
contains
    procedure :: execute => run_oversized_team_pass
end type oversized_team_pass

!> selects a state half among particles of which one is in no half (eo=-1), which must stop
type, extends(forked_process) :: halfless_selection
contains
    procedure :: execute => select_halfless
end type halfless_selection

contains

    subroutine run_all_strategy3D_inplane_tests()
        write(*,'(A)') '**** running all strategy3D in-plane tests ****'
        call test_grid_seed()
        call test_storage_without_continuous_arrays()
        call test_storage_with_continuous_arrays()
        call test_joint_storage()
        call test_discrete_seed_storage_and_invalid_predicate()
        call test_resolve_inplane_e3()
        call test_inpl_cont_policy()
        call test_inpl_cont_derivation()
        call test_key_derivation()
        call test_previous_shift_is_stored_shift()
        call test_table_increment_axes()
        call test_map_initialised_polish()
        call test_continuous_seed_beyond_window()
        call test_rec_observation_mask()
        call test_pass_team_check()
        call test_pass_state_team_sized()
        call test_rec_state_half_stops_on_eo()
    end subroutine run_all_strategy3D_inplane_tests

    !> the reconstruction's state halves (rec_state_half, also the PCG backend's) take the used particles' half
    !! through get_half: a used particle with eo=-1 stops; a rejected one (state 0) is not looked at
    subroutine test_rec_state_half_stops_on_eo()
        class(builder), allocatable, target :: build
        type(halfless_selection) :: child
        integer, allocatable :: selected(:)
        write(*,'(A)') 'test_rec_state_half_stops_on_eo'
        allocate(build)
        call make_half_field(build, 0)
        call rec_state_half(build, 1, 0, [1,2,3,4], selected)
        call assert_int(2, size(selected), 'the even half of state 1 without the rejected particle')
        call child%start(name=string('REC_STATE_HALF_EO'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'a used particle in no half stops the selection')
        call child%destroy
        call build%spproj%kill
    end subroutine test_rec_state_half_stops_on_eo

    !> four particles, eo [0,1,-1,0]; particle 3 has state @p state3, the others state 1
    subroutine make_half_field( build, state3 )
        class(builder), target, intent(inout) :: build
        integer,                intent(in)    :: state3
        integer :: i
        integer, parameter :: EOS(4) = [0, 1, -1, 0]
        call build%spproj%os_ptcl3D%new(4, is_ptcl=.true.)
        do i = 1, 4
            call build%spproj%os_ptcl3D%set_state(i, 1)
            call build%spproj%os_ptcl3D%set(i, 'eo', EOS(i))
        enddo
        call build%spproj%os_ptcl3D%set_state(3, state3)
        build%spproj_field => build%spproj%os_ptcl3D
    end subroutine make_half_field

    subroutine select_halfless( self, cline )
        class(halfless_selection), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(builder), allocatable, target :: build
        integer, allocatable :: selected(:)
        allocate(build)
        call make_half_field(build, 1)
        call rec_state_half(build, 1, 0, [1,2,3,4], selected)
    end subroutine select_halfless

    !> a pass's per-thread particle images are sized by the team its regions run, not by nthr_glob
    !! (OpenMP builds only)
    subroutine test_pass_state_team_sized()
        class(parameters), allocatable :: params
        class(builder),    allocatable :: build
        type(image),       allocatable :: ptcl_imgs(:), ptcl_imgs_pad(:)
        integer :: nthr_saved, nteam_saved, i
        logical :: l_omp
        write(*,'(A)') 'test_pass_state_team_sized'
        l_omp = .false.
        !$ l_omp = .true.
        if( .not. l_omp ) return
        nthr_saved  = nthr_glob
        nteam_saved = 1
        !$ nteam_saved = omp_get_max_threads()
        allocate(params, build)
        params%box        = 16
        params%smpd       = 2.
        params%box_crop   = 16
        params%box_croppd = 32
        params%smpd_crop  = 2.
        nthr_glob = 1
        !$ call omp_set_num_threads(2)
        call assert_int(2, pass_team(), 'the pass team is the OpenMP team')
        call alloc_ptcl_imgs(params, build, ptcl_imgs, ptcl_imgs_pad, 2)
        call assert_int(2, size(ptcl_imgs),     'the particle images follow the team, not nthr_glob')
        call assert_int(2, size(ptcl_imgs_pad), 'the padded particle images follow the team, not nthr_glob')
        do i = 1, size(ptcl_imgs)
            call ptcl_imgs(i)%kill
            call ptcl_imgs_pad(i)%kill
        enddo
        deallocate(ptcl_imgs, ptcl_imgs_pad)
        call killimgbatch(build)
        nthr_glob = nthr_saved
        !$ call omp_set_num_threads(nteam_saved)
    end subroutine test_pass_state_team_sized

    !> a team of at most nthr_glob passes; a larger one stops (in a child; OpenMP builds only)
    subroutine test_pass_team_check()
        type(oversized_team_pass) :: child
        integer :: nthr_saved, nteam_saved
        logical :: l_omp
        write(*,'(A)') 'test_pass_team_check'
        l_omp = .false.
        !$ l_omp = .true.
        if( .not. l_omp ) return
        nthr_saved  = nthr_glob
        nteam_saved = 1
        !$ nteam_saved = omp_get_max_threads()
        nthr_glob = 2
        !$ call omp_set_num_threads(2)
        call check_pass_team('test_pass_team_check')
        !$ call omp_set_num_threads(1)
        call check_pass_team('test_pass_team_check')
        call child%start(name=string('PASS_TEAM_OVERSIZED'))
        call child%await_final_status()
        call assert_int(FORK_STATUS_FAILED, child%status(), 'a pass on a team larger than nthr_glob stops')
        call child%destroy
        nthr_glob = nthr_saved
        !$ call omp_set_num_threads(nteam_saved)
    end subroutine test_pass_team_check

    subroutine run_oversized_team_pass( self, cline )
        class(oversized_team_pass), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        nthr_glob = 1
        !$ call omp_set_num_threads(2)
        call check_pass_team('run_oversized_team_pass')
    end subroutine run_oversized_team_pass

    !> nrefs candidates, one thread; the continuous arrays only on request
    subroutine alloc_state( nrefs, continuous )
        integer, intent(in) :: nrefs
        logical, intent(in) :: continuous
        allocate(s3D%proj_space_shift(2,nrefs,1), s3D%proj_space_corrs(nrefs,1), s3D%proj_space_inplinds(nrefs,1))
        s3D%proj_space_shift    = 0.
        s3D%proj_space_corrs    = -huge(1.)
        s3D%proj_space_inplinds = 0
        if( continuous )then
            allocate(s3D%proj_space_inplcoords(nrefs,1), s3D%proj_space_inplvalid(nrefs,1))
            s3D%proj_space_inplcoords = 0.
            s3D%proj_space_inplvalid  = .false.
        endif
    end subroutine alloc_state

    subroutine test_grid_seed()
        real    :: coordinate
        logical :: valid
        write(*,'(A)') 'test_grid_seed'
        call seed_continuous_inplane_candidate(1, coordinate, valid)
        call assert_real(1., coordinate, 1.e-6, 'the first grid index is preserved as the continuous seed')
        call assert_false(valid, 'a grid seed is not a valid continuous result')
        call seed_continuous_inplane_candidate(288, coordinate, valid)
        call assert_real(288., coordinate, 1.e-6, 'the upper grid index is preserved as the continuous seed')
        call assert_false(valid, 'the upper grid seed is not a valid continuous result')
    end subroutine test_grid_seed

    !> default-off storage holds only score, shift and integer angle; storing must not
    !! require or create continuous-only state
    subroutine test_storage_without_continuous_arrays()
        type(strategy3D_srch) :: srch
        write(*,'(A)') 'test_storage_without_continuous_arrays'
        call alloc_state(2, .false.)
        srch%ithr = 1
        call srch%store_solution(1, 9, 0.25, sh=[0.5, -0.25])
        call assert_int(9, s3D%proj_space_inplinds(1,1), 'integer in-plane index stored')
        call assert_real(0.25, s3D%proj_space_corrs(1,1), 1.e-6, 'score stored')
        call assert_true(all(abs(s3D%proj_space_shift(:,1,1) - [0.5, -0.25]) <= 1.e-6), 'shift stored')
        call assert_false(allocated(s3D%proj_space_inplcoords) .or. allocated(s3D%proj_space_inplvalid), &
            &'default-off storage creates no continuous-only state')
        call clean_strategy3D
    end subroutine test_storage_without_continuous_arrays

    subroutine test_storage_with_continuous_arrays()
        type(strategy3D_srch) :: srch
        write(*,'(A)') 'test_storage_with_continuous_arrays'
        call alloc_state(2, .true.)
        srch%ithr   = 1
        srch%nsolns = 0
        call srch%store_solution(2, 17, 0.5, sh=[1.25, -0.75])
        call assert_int(1, srch%nsolns, 'the first stored candidate increments the solution count')
        call assert_int(17, s3D%proj_space_inplinds(2,1), 'stored integer in-plane index')
        call assert_real(17., s3D%proj_space_inplcoords(2,1), 1.e-6, 'the continuous coordinate is the grid seed')
        call assert_false(s3D%proj_space_inplvalid(2,1), 'a legacy stored candidate is not continuous-valid')
        call assert_true(all(s3D%proj_space_shift(:,2,1) == [1.25, -0.75]), 'stored shift')
        ! a rejected (worse) candidate leaves an accepted continuous result alone
        s3D%proj_space_inplcoords(2,1) = 17.25
        s3D%proj_space_inplvalid(2,1)  = .true.
        call srch%store_solution(2, 18, 0.4, sh=[2., 2.])
        call assert_real(17.25, s3D%proj_space_inplcoords(2,1), 1.e-6, 'a rejected candidate keeps the continuous coordinate')
        call assert_true(s3D%proj_space_inplvalid(2,1), 'a rejected candidate keeps the continuous validity')
        ! an improving grid candidate replaces and invalidates stale continuous state
        call srch%store_solution(2, 19, 0.6, sh=[0.5, -0.5])
        call assert_int(19, s3D%proj_space_inplinds(2,1), 'an improving candidate replaces the index')
        call assert_real(19., s3D%proj_space_inplcoords(2,1), 1.e-6, 'an improving candidate resets the coordinate to its grid')
        call assert_false(s3D%proj_space_inplvalid(2,1), 'an improving candidate invalidates the stale continuous result')
        call clean_strategy3D
    end subroutine test_storage_with_continuous_arrays

    subroutine test_joint_storage()
        type(strategy3D_srch) :: srch
        write(*,'(A)') 'test_joint_storage'
        call alloc_state(2, .true.)
        srch%ithr = 1
        call srch%store_solution(2, 17, 0.5, sh=[1.25, -0.75])
        call srch%store_continuous_solution(2, 18, 17.625_dp, 0.6, [1.1, -0.7])
        call assert_int(18, s3D%proj_space_inplinds(2,1), 'the joint result keeps its nearest grid index')
        call assert_real(17.625, s3D%proj_space_inplcoords(2,1), 1.e-6, 'the joint result keeps its continuous coordinate')
        call assert_true(s3D%proj_space_inplvalid(2,1), 'the joint result is marked valid')
        call assert_real(0.6, s3D%proj_space_corrs(2,1), 1.e-6, 'the joint score is committed')
        call assert_true(all(abs(s3D%proj_space_shift(:,2,1) - [1.1, -0.7]) <= 1.e-6), 'the joint shift is committed')
        call assert_true(s3D%proj_space_corrs(1,1) <= -huge(1.)/2., 'another reference candidate is untouched')
        call clean_strategy3D
    end subroutine test_joint_storage

    !> polish-only policy: the joint solve only polishes the committed pose; a finite
    !! no-improvement result stores the selected seed and is not continuous-valid
    subroutine test_discrete_seed_storage_and_invalid_predicate()
        type(strategy3D_srch) :: srch
        real, parameter :: SEED_SHIFT(2) = [1.25, -0.75]
        write(*,'(A)') 'test_discrete_seed_storage_and_invalid_predicate'
        srch%continuous_active = .true.
        call assert_false(srch%joint_evaluation_invalid(.true.), 'a valid joint no-improvement result is not invalid')
        call assert_true(srch%joint_evaluation_invalid(.false.), 'an invalid joint result is identified for retention')
        call alloc_state(2, .true.)
        s3D%proj_space_corrs(2,1)      = 0.625
        s3D%proj_space_inplcoords(2,1) = 17.
        s3D%proj_space_inplinds(2,1)   = 17
        srch%ithr = 1
        call srch%store_discrete_seed_solution(2, 23, 0.75, SEED_SHIFT)
        call assert_true(all(abs(s3D%proj_space_shift(:,2,1) - SEED_SHIFT) <= 1.e-6), 'the selected seed shift is stored')
        call assert_real(0.75, s3D%proj_space_corrs(2,1), 1.e-6, 'the selected seed score is stored')
        call assert_int(23, s3D%proj_space_inplinds(2,1), 'the incoming grid angle is replaced')
        call assert_real(23., s3D%proj_space_inplcoords(2,1), 1.e-6, 'the coordinate follows the new grid angle')
        call assert_false(s3D%proj_space_inplvalid(2,1), 'a finite no-improvement is not continuous-valid')
        call assert_true(all(s3D%proj_space_shift(:,1,1) == 0.), 'another reference is untouched')
        call clean_strategy3D
    end subroutine test_discrete_seed_storage_and_invalid_predicate

    subroutine test_resolve_inplane_e3()
        real, parameter :: DANG = 1.25
        write(*,'(A)') 'test_resolve_inplane_e3'
        call assert_real(340., resolve_inplane_e3(17, 17.4, .true.,  .false., DANG), 1.e-5, &
            &'inactive continuous mode keeps the integer grid angle')
        call assert_real(340., resolve_inplane_e3(17, 17.4, .false., .true.,  DANG), 1.e-5, &
            &'an invalid continuous coordinate keeps the integer grid angle')
        call assert_real(339.5, resolve_inplane_e3(17, 17.4, .true.,  .true.,  DANG), 1.e-5, &
            &'an accepted continuous coordinate is converted to e3')
        call assert_real(0., resolve_inplane_e3(1, 1., .true., .true., DANG), 1.e-5, &
            &'the first grid angle maps to e3 = 0 in continuous mode')
    end subroutine test_resolve_inplane_e3

    !> inpl_cont: internal, default yes (refine=cont turns it off), not a command-line key
    !! since release 4, so no refine program exposes it; pose_cont is stripped from a
    !! child command line as a matcher-only option
    subroutine test_inpl_cont_policy()
        type(cmdline)    :: cline, child_cline
        type(parameters) :: defaults
        write(*,'(A)') 'test_inpl_cont_policy'
        call assert_char('yes', trim(defaults%inpl_cont), 'inpl_cont defaults to yes')
        call make_ui
        call assert_false(has_search_input('refine3D',        'inpl_cont'), 'refine3D does not expose inpl_cont')
        call assert_false(has_search_input('refine3D_auto',   'inpl_cont'), 'refine3D_auto does not expose inpl_cont')
        call assert_false(has_search_input('refine3D_states', 'inpl_cont'), 'refine3D_states does not expose inpl_cont')
        call assert_false(has_search_input('solve3D',         'inpl_cont'), 'solve3D does not expose inpl_cont')
        call assert_false(has_search_input('solve2D',         'inpl_cont'), 'solve2D does not expose inpl_cont')
        call cline%set('prg', 'refine3D')
        call cline%set('pose_cont', 'yes')
        child_cline = cline
        call strip_refine3D_search_only_args(child_cline)
        call assert_false(child_cline%defined('pose_cont'), 'pose_cont is stripped from the child command line')
        call assert_true(cline%defined('pose_cont'), 'the parent command line keeps pose_cont')
        call child_cline%kill
        call cline%kill
    end subroutine test_inpl_cont_policy

    !> the phase derivation: refine=cont turns inpl_cont off and is a Cartesian pass, the polish
    !! with pose_cont=yes; a discrete mode keeps the internal default
    subroutine test_inpl_cont_derivation()
        class(parameters), allocatable :: params
        type(cmdline) :: cline
        write(*,'(A)') 'test_inpl_cont_derivation'
        call cline%set('smpd',    1.5)
        call cline%set('box',     32)
        call cline%set('mskdiam', 36.)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        call cline%set('nthr',    1)
        call cline%set('refine',  'cont')
        allocate(params)
        call params%new(cline)
        call assert_char('no', trim(params%inpl_cont), 'refine=cont turns inpl_cont off')
        call assert_true(params%l_cart_refine .and. .not. params%l_cont_polish, 'refine=cont is a Cartesian pass')
        deallocate(params)
        call cline%set('pose_cont', 'yes')
        allocate(params)
        call params%new(cline)
        call assert_true(params%l_cart_refine .and. params%l_cont_polish, 'refine=cont pose_cont=yes is the polish')
        deallocate(params)
        call cline%delete('pose_cont')
        call cline%set('refine', 'shc')
        allocate(params)
        call params%new(cline)
        call assert_char('yes', trim(params%inpl_cont), 'a discrete mode keeps inpl_cont=yes')
        call assert_false(params%l_cart_refine .or. params%l_cont_polish, 'a discrete mode is not a Cartesian pass')
        deallocate(params)
        call cline%kill
    end subroutine test_inpl_cont_derivation

    !> the polar shift map's key: on by default under euclid in a polar mode with shifts, the classic modes
    !! included (ruling 20); the legacy increment prior off under it (the review's T3); zero switches the key off and
    !! gives the legacy prior back, as the stage controllers pass it
    subroutine test_key_derivation()
        class(parameters), allocatable :: params
        type(cmdline) :: cline
        character(len=12) :: modes(4) = [character(len=12) :: 'shc', 'neigh', 'prob', 'greedy']
        integer :: i
        write(*,'(A)') 'test_key_derivation'
        call cline%set('smpd',    1.5)
        call cline%set('box',     32)
        call cline%set('mskdiam', 36.)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call cline%set('nthr',    1)
        call cline%set('trs',     5.)
        do i = 1, size(modes)
            call cline%set('refine', trim(modes(i)))
            allocate(params)
            call params%new(cline)
            call assert_true(params%l_centring, 'the key centres under refine='//trim(modes(i)))
            call assert_false(params%l_shift_prior, 'no increment prior under the key, refine='//trim(modes(i)))
            deallocate(params)
        enddo
        call cline%set('shift_sigma', 0.)
        allocate(params)
        call params%new(cline)
        call assert_false(params%l_shift_map .or. params%l_centring, 'shift_sigma=0 switches the key off')
        call assert_true(params%l_shift_prior, 'shift_sigma=0 gives the increment prior back')
        deallocate(params)
        call cline%delete('shift_sigma')
        call cline%set('refine', 'cont')
        allocate(params)
        call params%new(cline)
        call assert_false(params%l_centring, 'the Cartesian pass does not centre')
        deallocate(params)
        call cline%kill
    end subroutine test_key_derivation

    !> a search control of the named UI program; with expected_default, one with that default
    logical function has_search_input( program_name, key, expected_default ) result( found )
        character(len=*),           intent(in) :: program_name, key
        character(len=*), optional, intent(in) :: expected_default
        type(ui_program), pointer :: program => null()
        type(list_iterator) :: iterator
        type(string) :: name
        class(*), allocatable :: value
        found = .false.
        name  = program_name
        call get_prg_ptr(name, program)
        call name%kill
        if( .not. associated(program) ) return
        iterator = program%srch_ctrls%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type(input => value)
                type is(ui_program_input)
                    if( input%param%key%to_char() == key )then
                        found = .true.
                        if( present(expected_default) ) found = input%param%has_default .and. &
                            &input%param%cval_default%to_char() == expected_default
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
    end function has_search_input

    !> N10: the matcher takes the shift increment of a search from the project orientation
    !! before and after srch, where it formerly read the search object's previous shift
    !! (s%prev_shvec). That shift is the one stored in the project before the search: both
    !! preparations of the search object (prep4srch, prep4prob) set it from the stored record.
    !! Fixture: a one-particle strategy3D toolbox on a small phantom. Expected value: the
    !! stored shift.
    subroutine test_previous_shift_is_stored_shift()
        character(len=*), parameter :: PHANTOM = 'tmp_strategy3D_inplane_tester_phantom.mrc'
        real,             parameter :: STORED_SHIFT(2) = [1.25, -0.75]
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)         :: cline
        type(strategy3D_srch) :: srch
        type(strategy3D_spec) :: spec
        write(*,'(A)') 'test_previous_shift_is_stored_shift'
        call one_particle_fixture(PHANTOM, STORED_SHIFT, b, p, cline, srch, spec)
        call assert_true(all(srch%prev_shvec == STORED_SHIFT), &
            &'the previous shift of the search object is not the shift stored in the project')
        call assert_true(all(b%spproj_field%get_2Dshift(1) == STORED_SHIFT), &
            &'preparing the search changed the stored shift')
        call srch%kill
        call clean_strategy3D
        call b%pftc%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        call cline%kill
        call del_file(PHANTOM)
    end subroutine test_previous_shift_is_stored_shift

    !> a one-particle strategy3D toolbox on a small phantom (two directions), the particle at a stored shift, the
    !! search object prepared for the probabilistic route; trs a shift search's window (none by default)
    subroutine one_particle_fixture( phantom, stored_shift, b, p, cline, srch, spec, trs, objfun, box_crop, mskdiam )
        character(len=*),              intent(in)    :: phantom
        real,                          intent(in)    :: stored_shift(2)
        type(builder),         target, intent(inout) :: b
        type(parameters),      target, intent(inout) :: p
        type(cmdline),                 intent(inout) :: cline
        type(strategy3D_srch),         intent(inout) :: srch
        type(strategy3D_spec),         intent(inout) :: spec
        real,                optional, intent(in)    :: trs
        character(len=*),    optional, intent(in)    :: objfun
        integer,             optional, intent(in)    :: box_crop
        real,                optional, intent(in)    :: mskdiam
        integer, parameter :: BOX = 32
        type(image)        :: vol
        real, allocatable  :: rmat(:,:,:)
        integer :: i, j, k
        allocate(rmat(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    rmat(i,j,k) = exp(-real((i-14)**2 + (j-18)**2 + (k-16)**2)/18.) + &
                        &0.5*exp(-real((i-20)**2 + (j-13)**2 + (k-19)**2)/12.)
                end do
            end do
        end do
        call vol%new([BOX,BOX,BOX], 1.5)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(phantom), del_if_exists=.true.)
        call vol%kill
        call cline%set('vol1',    phantom)
        call cline%set('mskdiam', 36.)
        call cline%set('smpd',    1.5)
        call cline%set('lp',      6.)
        call cline%set('nptcls',  1.)
        call cline%set('nspace',  2)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'cc')
        if( present(objfun) )   call cline%set('objfun',   objfun)
        if( present(trs) )      call cline%set('trs',      trs)
        if( present(box_crop) ) call cline%set('box_crop', box_crop)
        if( present(mskdiam) )  call cline%set('mskdiam',  mskdiam)
        call cline%check
        call b%init_params_and_build_strategy3D_tbox(cline, p)
        call set_bp_range3D(p, b, cline)
        call b%pftc%new(p, p%nspace, [1,1], p%kfromto)
        call b%spproj_field%set_euler(1, [31., 57., 113.])
        call b%spproj_field%set_shift(1, stored_shift)
        call b%spproj_field%set_state(1, 1)
        call prep_strategy3D(p, b)
        spec%iptcl     = 1
        spec%iptcl_map = 1
        call srch%new(p, spec, b)
        srch%ithr = 1
        call srch%prep4prob
    end subroutine one_particle_fixture

    !> One convention for the table's increments (review R1): every increment is in the particle's axes. The window
    !! path writes minimize(sh_rot=.true.), the optimizer's shift rotated from the reference's axes at the rotation it
    !! was found at; the map path writes the drawn cell rotated the same way; the consumer adds the increment to the
    !! stored shift unconverted
    subroutine test_table_increment_axes()
        character(len=*), parameter :: PHANTOM = 'tmp_strategy3D_inplane_tester_phantom.mrc'
        real,             parameter :: STORED_SHIFT(2) = [1.25, -0.75], S_TRUE(2) = [1.5, -2.0], TRS = 5.
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)            :: cline
        type(strategy3D_srch)    :: srch
        type(strategy3D_spec)    :: spec
        type(pftc_shsrch_grad)   :: grad_native, grad_ref
        complex(sp), allocatable :: ref(:,:), ptcl(:,:)
        real    :: lims(2,2), cxy_native(3), cxy_ref(3), euler(3), shvec(2), u(2)
        integer :: irot_n, irot_r, nrots, inpl, nchecked, pdim(3), kfromto(2), pftsz, k, pp
        write(*,'(A)') 'test_table_increment_axes'
        call set_fixed_seed(20261031)
        call one_particle_fixture(PHANTOM, STORED_SHIFT, b, p, cline, srch, spec)
        nrots   = b%pftc%get_nrots()
        pftsz   = b%pftc%get_pftsz()
        kfromto = b%pftc%get_kfromto()
        pdim    = b%pftc%get_pdim_interp()
        ! a white reference on the band; the particle is that reference rotated, then shifted by S_TRUE about its own axes
        allocate(ref(pdim(1),pdim(2):pdim(3)), ptcl(pdim(1),pdim(2):pdim(3)), source=cmplx(0.,0.))
        do k = kfromto(1), kfromto(2)
            do pp = 1, pftsz
                ref(pp,k) = cmplx(gasdev(), gasdev())
            enddo
        enddo
        call b%pftc%set_ref_pft(1, ref, iseven=.true.)
        call b%pftc%set_eo(1, .true.)
        call b%pftc%memoize_refs
        lims(:,1) = -TRS
        lims(:,2) =  TRS
        call grad_native%new_alternating(b, lims, lims_init=lims, maxits=100, coarse_init=.true.)
        call grad_ref%new_alternating(b, lims, lims_init=lims, maxits=100, coarse_init=.true.)
        nchecked = 0
        do inpl = 1, nrots, max(1, nrots / 5)
            call rotate_half_ring(pftsz, kfromto(2)-kfromto(1)+1, ref(:,kfromto(1):kfromto(2)), inpl, ptcl(:,kfromto(1):kfromto(2)))
            call b%pftc%set_ptcl_pft(1, ptcl)
            call b%pftc%shift_ptcl(1, S_TRUE)
            call b%pftc%memoize_ptcls
            irot_n = inpl
            irot_r = inpl
            call grad_native%set_indices(1, 1)
            cxy_native = grad_native%minimize(irot=irot_n, sh_rot=.true.)
            call grad_ref%set_indices(1, 1)
            cxy_ref    = grad_ref%minimize(irot=irot_r, sh_rot=.false.)
            call assert_true(irot_n > 0 .and. irot_r > 0, 'the shift search found a better pose than the start')
            if( irot_n < 1 .or. irot_r < 1 ) cycle
            if( irot_n /= irot_r ) cycle
            nchecked = nchecked + 1
            call assert_true(abs(norm2(cxy_ref(2:3)) - norm2(S_TRUE)) <= 0.1, 'the shift search recovers the shift''s magnitude')
            u = shift_ref2ptcl(cxy_ref(2:3), b%pftc%get_rot(irot_r))
            call assert_true(all(abs(cxy_native(2:3) - u) <= 1.e-4), 'the window path''s table increment is the optimizer''s shift rotated into the particle''s axes')
            call candidate_pose(srch, 1, irot_n, cxy_native(2:3), euler, shvec)
            call assert_true(all(abs(shvec - (STORED_SHIFT + cxy_native(2:3))) <= 1.e-5), 'the consumer adds the native increment to the stored shift unconverted')
        enddo
        call assert_true(nchecked >= 3, 'the two searches agreed on the rotation at three starts or more')
        call grad_native%kill
        call grad_ref%kill
        call srch%kill
        call clean_strategy3D
        call b%pftc%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        call cline%kill
        call del_file(PHANTOM)
    end subroutine test_table_increment_axes

    !> Phase 4's pin: the map-initialised polish (the coarse window map's cell, then the L-BFGS from it in a box of
    !! one map step about it, minimize_in_cell) reaches the lattice-initialised search's optimum on the
    !! white-reference fixture, at the same rotation
    subroutine test_map_initialised_polish()
        character(len=*), parameter :: PHANTOM = 'tmp_strategy3D_inplane_tester_phantom.mrc'
        real,             parameter :: STORED_SHIFT(2) = [0., 0.], S_TRUE(2) = [1.5, -2.0], TRS = 5.
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)            :: cline
        type(strategy3D_srch)    :: srch
        type(strategy3D_spec)    :: spec
        type(pftc_shsrch_grad)   :: grad_lattice, grad_map
        complex(sp), allocatable :: ref(:,:), ptcl(:,:)
        real(sp), allocatable :: vals(:,:,:)
        real     :: lims(2,2), cxy_l(3), cxy_m(3), xy_cell(2), step, val
        real(dp) :: prior_w, xy_dp(2)
        integer  :: irot_l, irot_m, nrots, inpl, nchecked, pdim(3), kfromto(2), pftsz, k, pp, kscreen, irho, ialpha
        integer  :: nrho, nalpha
        logical  :: l_polished
        write(*,'(A)') 'test_map_initialised_polish'
        call set_fixed_seed(20261101)
        call one_particle_fixture(PHANTOM, STORED_SHIFT, b, p, cline, srch, spec)
        nrots   = b%pftc%get_nrots()
        pftsz   = b%pftc%get_pftsz()
        kfromto = b%pftc%get_kfromto()
        pdim    = b%pftc%get_pdim_interp()
        allocate(ref(pdim(1),pdim(2):pdim(3)), ptcl(pdim(1),pdim(2):pdim(3)), source=cmplx(0.,0.))
        do k = kfromto(1), kfromto(2)
            do pp = 1, pftsz
                ref(pp,k) = cmplx(gasdev(), gasdev())
            enddo
        enddo
        call b%pftc%set_ref_pft(1, ref, iseven=.true.)
        call b%pftc%set_eo(1, .true.)
        call b%pftc%memoize_refs
        lims(:,1) = -TRS
        lims(:,2) =  TRS
        call grad_lattice%new_alternating(b, lims, lims_init=lims, maxits=100, coarse_init=.true.)
        call grad_map%new_alternating(b, lims, lims_init=lims, maxits=100, coarse_init=.true.)
        call b%pftc%setup_shiftmaps(3.d0, real(TRS, dp), real(TRS, dp), prior_w, kscreen)
        step = real(norm2(b%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 2, 1) - b%pftc%shiftmap_shift(SHIFTMAP_WINDOW, 1, 1)))
        call b%pftc%shiftmap_dims(SHIFTMAP_WINDOW, nrho, nalpha)
        allocate(vals(nalpha, nrots, nrho))
        nchecked = 0
        do inpl = 1, nrots, max(1, nrots / 5)
            call rotate_half_ring(pftsz, kfromto(2)-kfromto(1)+1, ref(:,kfromto(1):kfromto(2)), inpl, ptcl(:,kfromto(1):kfromto(2)))
            call b%pftc%set_ptcl_pft(1, ptcl)
            call b%pftc%shift_ptcl(1, S_TRUE)
            call b%pftc%memoize_ptcls
            irot_l = inpl
            call grad_lattice%set_indices(1, 1)
            cxy_l = grad_lattice%minimize(irot=irot_l, sh_rot=.true.)
            ! the map's peak is the initialiser's cell (objective-agnostic, as the fills' greedy draw under euclid)
            call b%pftc%gen_shiftmap(SHIFTMAP_WINDOW, 1, 1)
            call b%pftc%shiftmap_objfun_grid(SHIFTMAP_WINDOW, vals)
            call b%pftc%shiftmap_peak(SHIFTMAP_WINDOW, vals, irot_m, xy_dp, val, irho, ialpha)
            xy_cell = real(xy_dp)
            call grad_map%set_indices(1, 1)
            call grad_map%minimize_in_cell(irot_m, xy_cell, step, cxy_m, l_polished)
            call assert_true(irot_l > 0, 'the lattice-initialised search found a better pose than the start')
            call assert_true(l_polished, 'the map-initialised polish beat the origin')
            if( irot_l < 1 ) cycle
            nchecked = nchecked + 1
            call assert_true(abs(norm2(cxy_m(2:3)) - norm2(S_TRUE)) <= 0.15, 'the map-initialised polish recovers the shift''s magnitude')
            call assert_true(all(abs(cxy_m(2:3) - cxy_l(2:3)) <= 0.15), 'the map-initialised polish reaches the lattice-initialised optimum')
            call assert_true(irot_m == irot_l, 'the two starts agree on the rotation')
        enddo
        call assert_true(nchecked >= 3, 'the identity held at three starts or more')
        call grad_lattice%kill
        call grad_map%kill
        call srch%kill
        call clean_strategy3D
        call b%pftc%kill_shiftmap
        call b%pftc%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        call cline%kill
        call del_file(PHANTOM)
    end subroutine test_map_initialised_polish

    !> The review's T1: the continuous route's joint domain contains the incoming pose. A pose found on a disk beyond
    !! the window (the white-reference fixture shifted by S_TRUE, polished within +-8 px) goes through the committed
    !! pose's continuous route under trs = 2: it comes back unclipped and scored no worse than itself
    subroutine test_continuous_seed_beyond_window()
        character(len=*), parameter :: PHANTOM = 'tmp_strategy3D_inplane_tester_phantom.mrc'
        real,             parameter :: STORED_SHIFT(2) = [0., 0.], S_TRUE(2) = [4.0, -3.5], RDISK = 8.
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)            :: cline
        type(strategy3D_srch)    :: srch
        type(strategy3D_spec)    :: spec
        type(pftc_shsrch_grad)   :: grad_disk
        complex(sp), allocatable :: ref(:,:), ptcl(:,:)
        real    :: lims(2,2), cxy(3), sh(2), corr, inpl_coord
        integer :: irot, inpl, nrots, pdim(3), kfromto(2), pftsz, k, pp, nchecked
        logical :: inpl_valid
        write(*,'(A)') 'test_continuous_seed_beyond_window'
        call set_fixed_seed(20261103)
        call one_particle_fixture(PHANTOM, STORED_SHIFT, b, p, cline, srch, spec, trs=2.)
        call assert_true(srch%uses_continuous_refinement(), 'the fixture runs the continuous route')
        nrots   = b%pftc%get_nrots()
        pftsz   = b%pftc%get_pftsz()
        kfromto = b%pftc%get_kfromto()
        pdim    = b%pftc%get_pdim_interp()
        allocate(ref(pdim(1),pdim(2):pdim(3)), ptcl(pdim(1),pdim(2):pdim(3)), source=cmplx(0.,0.))
        do k = kfromto(1), kfromto(2)
            do pp = 1, pftsz
                ref(pp,k) = cmplx(gasdev(), gasdev())
            enddo
        enddo
        call b%pftc%set_ref_pft(1, ref, iseven=.true.)
        call b%pftc%set_eo(1, .true.)
        call b%pftc%memoize_refs
        lims(:,1) = -RDISK
        lims(:,2) =  RDISK
        call grad_disk%new_alternating(b, lims, lims_init=lims, maxits=100, coarse_init=.true.)
        nchecked = 0
        do inpl = 1, nrots, max(1, nrots / 5)
            call rotate_half_ring(pftsz, kfromto(2)-kfromto(1)+1, ref(:,kfromto(1):kfromto(2)), inpl, ptcl(:,kfromto(1):kfromto(2)))
            call b%pftc%set_ptcl_pft(1, ptcl)
            call b%pftc%shift_ptcl(1, S_TRUE)
            call b%pftc%memoize_ptcls
            irot = inpl
            call grad_disk%set_indices(1, 1)
            cxy = grad_disk%minimize(irot=irot, sh_rot=.true.)
            if( irot < 1 ) cycle
            if( maxval(abs(cxy(2:3))) <= 2. * sqrt(2.) ) cycle
            nchecked = nchecked + 1
            sh   = cxy(2:3)
            corr = 0.
            call srch%refine_assignment_continuously(1, irot, corr, sh, inpl_coord, inpl_valid)
            call assert_true(all(abs(sh - cxy(2:3)) <= 0.25), 'the continuous route keeps a seed beyond the window unclipped')
            call assert_true(corr >= cxy(1) - 1.e-4, 'the continuous route scores the seed itself, not its clipped pose')
        enddo
        call assert_true(nchecked >= 3, 'the seed lay beyond the window at three starts or more')
        call grad_disk%kill
        call srch%kill
        call clean_strategy3D
        call b%pftc%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        call cline%kill
        call del_file(PHANTOM)
    end subroutine test_continuous_seed_beyond_window

    !> The review's M2: a reconstruction observation is the same whether its record is kept or the particle is read
    !! again. The search preparation of a particle on the disk normalises it under the box's mask and notes the choice
    !! (note_obs_mask); its record (the prepared transform's unshifted crop) equals the reread through
    !! prep_rec_observation_ptcl, which the ordinary mask would not give (an off-centred Gaussian on a structured
    !! background, as the review's probe); a particle in its window notes the ordinary mask
    subroutine test_rec_observation_mask()
        character(len=*), parameter :: PHANTOM = 'tmp_strategy3D_inplane_tester_phantom.mrc'
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)            :: cline
        type(strategy3D_srch)    :: srch
        type(strategy3D_spec)    :: spec
        type(image)              :: raw, work, out, out_pd, rec, obs, obs_l
        real, allocatable        :: rmat(:,:,:), a(:,:,:), c(:,:,:)
        real    :: d_box, d_lmsk, scale
        integer :: i, j
        write(*,'(A)') 'test_rec_observation_mask'
        call one_particle_fixture(PHANTOM, [0., 0.], b, p, cline, srch, spec, trs=3., objfun='euclid', box_crop=28, mskdiam=18.)
        call assert_true(allocated(b%lmsk_box), 'the key builds the box''s mask')
        call assert_true(count(b%lmsk_box) > count(b%lmsk), 'the box''s mask is wider than the molecular one')
        call assert_false(b%obs_on_box(1), 'without a preparation the ordinary mask')
        ! the stack record the preparation reads the CTF parameters from
        call b%spproj_field%set(1, 'stkind', 1)
        call b%spproj_field%set(1, 'indstk', 1)
        call b%spproj%os_stk%new(1, is_ptcl=.false.)
        call b%spproj%os_stk%set(1, 'nptcls_stk', 1)
        call b%spproj%os_stk%set(1, 'smpd',       p%smpd)
        call b%spproj%os_stk%set(1, 'ctf',        'no')
        call b%spproj%os_stk%set(1, 'kv',         300.)
        call b%spproj%os_stk%set(1, 'cs',         2.7)
        call b%spproj%os_stk%set(1, 'fraca',      0.1)
        allocate(rmat(p%box,p%box,1))
        do j = 1, p%box
            do i = 1, p%box
                rmat(i,j,1) = sin(0.37 * real(i)) + cos(0.23 * real(j)) + 6. * exp(-real((i-24)**2 + (j-22)**2) / 10.)
            enddo
        enddo
        call raw%new([p%box,p%box,1], p%smpd)
        call raw%set_rmat(rmat, .false.)
        call work%new([p%box,p%box,1], p%smpd)
        call out%new([p%box_crop,p%box_crop,1], p%smpd_crop)
        call out%memoize_mask_coords
        call memoize_ft_maps(out%get_ldim(), out%get_smpd())
        call out_pd%new([p%box_croppd,p%box_croppd,1], p%smpd_crop)
        call rec%new(  [p%box_crop,p%box_crop,1], p%smpd_crop)
        call obs%new(  [p%box_crop,p%box_crop,1], p%smpd_crop)
        call obs_l%new([p%box_crop,p%box_crop,1], p%smpd_crop)
        ! on the disk (not yet centred)
        call b%spproj_field%set(1, 'shsig', -1.)
        call record_and_reread
        call assert_true(b%obs_on_box(1), 'the disk particle''s preparation notes the box''s mask')
        call work%copy(raw)
        call prep_rec_observation(work, b%lmsk, obs_l, .false.)
        c      = obs_l%get_rmat()
        d_lmsk = maxval(abs(a - c)) / scale
        write(*,'(A,ES10.2,A,ES10.2)') '  disk: |record - reread| ', d_box, ', |record - ordinary-mask reread| ', d_lmsk
        call assert_true(d_box <= 1.e-5, 'the reread equals the record under the noted mask')
        call assert_true(d_lmsk > 1.e-2, 'the ordinary mask would give another observation')
        ! in its window
        call b%spproj_field%set(1, 'shsig', 0.1)
        call record_and_reread
        call assert_false(b%obs_on_box(1), 'the window particle''s preparation notes the ordinary mask')
        call assert_true(d_box <= 1.e-5, 'its reread equals its record too')
        call raw%kill
        call work%kill
        call out%kill
        call out_pd%kill
        call rec%kill
        call obs%kill
        call obs_l%kill
        call srch%kill
        call clean_strategy3D
        call b%pftc%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        call cline%kill
        call del_file(PHANTOM)

    contains

        !> the record writer's transform (the search preparation, then the unshifted crop) against the reread
        subroutine record_and_reread
            call work%copy(raw)
            call prepimg4align(p, b, 1, work, out, out_pd)
            call work%clip(rec)
            call rec%ifft
            call work%copy(raw)
            call prep_rec_observation_ptcl(b, 1, work, obs, .false.)
            a     = rec%get_rmat()
            c     = obs%get_rmat()
            scale = maxval(abs(a))
            d_box = maxval(abs(a - c)) / scale
        end subroutine record_and_reread

    end subroutine test_rec_observation_mask

    !> a polar section rotated by irot in the calculator's convention (the polar shift map tester's helper, checked there
    !! against the calculator's own rotation)
    pure subroutine rotate_half_ring( pftsz, nk, pft, irot, pft_rot )
        integer,     intent(in)  :: pftsz, nk, irot
        complex(sp), intent(in)  :: pft(pftsz,nk)
        complex(sp), intent(out) :: pft_rot(pftsz,nk)
        integer :: p, q
        do p = 1, pftsz
            q = modulo(p - irot, 2*pftsz) + 1
            if( q <= pftsz )then
                pft_rot(p,:) = pft(q,:)
            else
                pft_rot(p,:) = conjg(pft(q-pftsz,:))
            endif
        enddo
    end subroutine rotate_half_ring

end module simple_strategy3D_inplane_tester
