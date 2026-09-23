!@descr: unit tests for the continuous in-plane state of the refine3D search (simple_strategy3D_srch, _alloc, _utils)
! The search-state contract behind inpl_cont in refine3D, with no fixture: the grid seed
! of a continuous candidate, candidate storage with and without the continuous arrays
! (a rejected candidate leaves an accepted continuous result alone, an improving grid
! candidate resets it), the joint and the discrete-seed storage routes, the invalid
! joint-evaluation predicate, resolve_inplane_e3, and the inpl_cont policy: default
! yes, exposed with that default on the three refine3D programs, stripped from the
! child command line.
module simple_strategy3D_inplane_tester
use simple_core_module_api,   only: dp
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_string,            only: string
use simple_linked_list,       only: list_iterator
use simple_ui,                only: make_ui, get_prg_ptr
use simple_ui_program,        only: ui_program, ui_program_input
use simple_refine3D_strategy, only: strip_refine3D_search_only_args
use simple_strategy3D_alloc,  only: clean_strategy3D, s3D, seed_continuous_inplane_candidate
use simple_strategy3D_srch,   only: strategy3D_srch
use simple_strategy3D_utils,  only: resolve_inplane_e3
use simple_test_utils
implicit none
private
public :: run_all_strategy3D_inplane_tests

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
    end subroutine run_all_strategy3D_inplane_tests

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

    !> inpl_cont: default yes, exposed with that default on the refine3D programs,
    !! stripped from a child command line as a matcher-only option
    subroutine test_inpl_cont_policy()
        type(cmdline)    :: cline, child_cline
        type(parameters) :: defaults
        write(*,'(A)') 'test_inpl_cont_policy'
        call assert_char('yes', trim(defaults%inpl_cont), 'inpl_cont defaults to yes')
        call make_ui
        call assert_true(has_search_input('refine3D',        'inpl_cont', 'yes'), 'refine3D exposes inpl_cont=yes')
        call assert_true(has_search_input('refine3D_auto',   'inpl_cont', 'yes'), 'refine3D_auto exposes inpl_cont=yes')
        call assert_true(has_search_input('refine3D_states', 'inpl_cont', 'yes'), 'refine3D_states exposes inpl_cont=yes')
        call cline%set('prg', 'refine3D')
        call cline%set('inpl_cont', 'yes')
        child_cline = cline
        call strip_refine3D_search_only_args(child_cline)
        call assert_false(child_cline%defined('inpl_cont'), 'inpl_cont is stripped from the child command line')
        call assert_true(cline%defined('inpl_cont'), 'the parent command line keeps inpl_cont')
        call child_cline%kill
        call cline%kill
    end subroutine test_inpl_cont_policy

    !> a search control of the named UI program with the given default
    logical function has_search_input( program_name, key, expected_default ) result( found )
        character(len=*), intent(in) :: program_name, key, expected_default
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
                        found = input%param%has_default .and. input%param%cval_default%to_char() == expected_default
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
    end function has_search_input

end module simple_strategy3D_inplane_tester
