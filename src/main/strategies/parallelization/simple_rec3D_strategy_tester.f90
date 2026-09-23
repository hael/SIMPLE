!@descr: unit test routines for the rec3D backend selector (simple_rec3D_strategy)
! Pins the parameter defaults the selector relies on (rec_backend=gridding, maxits_pcg=2),
! backend-name resolution and wiring, and the dynamic type the factory returns for every
! branch: shared-memory gridding, distributed gridding (nparts without part), shared-memory
! PCG, distributed PCG (nparts > 1) and the worker part of either backend.
module simple_rec3D_strategy_tester
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_rec3D_strategy, only: rec3D_strategy, rec3D_inmem_strategy, rec3D_distr_strategy, &
    &rec3D_pcg_inmem_strategy, create_rec3D_strategy, rec3D_backend_id, rec3D_backend_is_wired, &
    &REC3D_BACKEND_INVALID, REC3D_BACKEND_GRIDDING, REC3D_BACKEND_PCG
use simple_test_utils
implicit none
private
public :: run_all_rec3D_strategy_tests

integer, parameter :: KIND_INMEM = 1, KIND_DISTR = 2, KIND_PCG_INMEM = 3, KIND_OTHER = 0

contains

    subroutine run_all_rec3D_strategy_tests()
        write(*,'(A)') '**** running all rec3D_strategy tests ****'
        call test_defaults()
        call test_backend_resolution()
        call test_factory_branches()
    end subroutine run_all_rec3D_strategy_tests

    !> the dynamic type of a strategy as a small integer, exact type (pcg_inmem extends inmem)
    integer function strategy_kind( strategy )
        class(rec3D_strategy), intent(in) :: strategy
        select type(strategy)
            type is(rec3D_inmem_strategy)
                strategy_kind = KIND_INMEM
            type is(rec3D_distr_strategy)
                strategy_kind = KIND_DISTR
            type is(rec3D_pcg_inmem_strategy)
                strategy_kind = KIND_PCG_INMEM
            class default
                strategy_kind = KIND_OTHER
        end select
    end function strategy_kind

    subroutine test_defaults()
        type(parameters) :: params
        write(*,'(A)') 'test_defaults'
        call assert_char('gridding', trim(params%rec_backend), 'rec_backend defaults to gridding')
        call assert_int(2, params%maxits_pcg, 'maxits_pcg defaults to 2')
        call assert_true(params%rtol <= 0., 'rtol defaults to <= 0 (exactly maxits_pcg iterations)')
    end subroutine test_defaults

    subroutine test_backend_resolution()
        write(*,'(A)') 'test_backend_resolution'
        call assert_int(REC3D_BACKEND_GRIDDING, rec3D_backend_id('gridding'), 'gridding resolves')
        call assert_int(REC3D_BACKEND_PCG,      rec3D_backend_id('pcg'),      'pcg resolves')
        call assert_int(REC3D_BACKEND_GRIDDING, rec3D_backend_id('gridding   '), 'trailing blanks are ignored')
        call assert_int(REC3D_BACKEND_INVALID,  rec3D_backend_id('invalid'),  'unknown name is invalid')
        call assert_int(REC3D_BACKEND_INVALID,  rec3D_backend_id('PCG'),      'names are case-sensitive')
        call assert_int(REC3D_BACKEND_INVALID,  rec3D_backend_id(''),         'empty name is invalid')
        call assert_true(rec3D_backend_is_wired(REC3D_BACKEND_GRIDDING), 'gridding is wired')
        call assert_true(rec3D_backend_is_wired(REC3D_BACKEND_PCG),      'pcg is wired')
        call assert_false(rec3D_backend_is_wired(REC3D_BACKEND_INVALID), 'invalid is not wired')
    end subroutine test_backend_resolution

    subroutine test_factory_branches()
        type(cmdline) :: cline
        class(rec3D_strategy), allocatable :: strategy
        write(*,'(A)') 'test_factory_branches'
        ! bare command line: shared-memory gridding
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_INMEM, strategy_kind(strategy), 'default is the gridding in-memory strategy')
        deallocate(strategy)
        ! nparts without part: the distributed master
        call cline%set('nparts', 4)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_DISTR, strategy_kind(strategy), 'gridding with nparts is the distributed strategy')
        deallocate(strategy)
        ! nparts with part: a worker runs in memory
        call cline%set('part', 1)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_INMEM, strategy_kind(strategy), 'a gridding worker part runs the in-memory strategy')
        deallocate(strategy)
        call cline%delete('part')
        call cline%delete('nparts')
        ! pcg backend, shared memory
        call cline%set('rec_backend', 'pcg')
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_PCG_INMEM, strategy_kind(strategy), 'pcg is the PCG in-memory strategy')
        deallocate(strategy)
        ! pcg with nparts > 1: distributed
        call cline%set('nparts', 2)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_DISTR, strategy_kind(strategy), 'pcg with nparts > 1 is the distributed strategy')
        deallocate(strategy)
        ! pcg with nparts = 1: stays in memory
        call cline%set('nparts', 1)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_PCG_INMEM, strategy_kind(strategy), 'pcg with nparts = 1 stays the PCG in-memory strategy')
        deallocate(strategy)
        ! pcg worker part
        call cline%set('nparts', 2)
        call cline%set('part', 2)
        strategy = create_rec3D_strategy(cline)
        call assert_int(KIND_PCG_INMEM, strategy_kind(strategy), 'a pcg worker part runs the PCG in-memory strategy')
        deallocate(strategy)
        call cline%kill
    end subroutine test_factory_branches

end module simple_rec3D_strategy_tester
