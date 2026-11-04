!@descr: unit test routines for the rec3D backend selector (simple_rec3D_strategy, simple_rec3D_service)
! Pins the parameter defaults the selector relies on (rec_backend=gridding, maxits_pcg=2),
! backend-name resolution and wiring, and the dynamic type the factory returns for every
! branch: shared-memory gridding, distributed gridding (nparts without part), shared-memory
! PCG, distributed PCG (nparts > 1) and the worker part of either backend; the PCG worker's entries of the
! weighted insertion (topk > 1).
module simple_rec3D_strategy_tester
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_rec3D_strategy, only: rec3D_strategy, rec3D_inmem_strategy, rec3D_distr_strategy, &
    &rec3D_pcg_inmem_strategy, create_rec3D_strategy
use simple_rec3D_service,  only: rec3D_backend_id, rec3D_backend_is_wired, REC3D_BACKEND_INVALID, &
    &REC3D_BACKEND_GRIDDING, REC3D_BACKEND_PCG
use simple_reconstructor_pcg, only: reconstructor_pcg, PCG_LAMBDA
use simple_rec3D_pcg_strategy, only: pcg_extra_entries, accumulate_pcg_extras
use simple_eul_prob_tab_utils, only: prob_insert_poses
use simple_ori,                only: ori
use simple_oris,               only: oris
use simple_type_defs,          only: ctfparams, CTFFLAG_NO
use simple_rnd,                only: ran3, seed_rnd_fixed
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
        call test_pcg_extra_entries()
    end subroutine run_all_rec3D_strategy_tests

    !> The PCG worker's weighted insertion: the extra poses as entries after the committed ones, each from its
    !! particle's observation. Extras at the committed pose with weights summing to one reproduce the plain
    !! accumulation, in one batch and in two; an extra moved off the committed pose changes it
    subroutine test_pcg_extra_entries()
        integer, parameter :: NP = 4, BOX = 32
        real,    parameter :: SMPD = 2.
        type(reconstructor_pcg) :: op
        type(oris)              :: sel, entries
        type(ori)               :: o
        type(ctfparams)         :: ctfp
        type(prob_insert_poses) :: extra
        complex, allocatable    :: planes(:,:,:), y_x(:,:,:), b(:,:,:), b0(:,:,:), b2(:,:,:)
        real,    allocatable    :: d(:,:,:), d0(:,:,:), d2(:,:,:), wts(:)
        integer, allocatable    :: xptcl(:), xfirst(:)
        integer :: pinds(NP), lims2(2,2), i, h, k
        write(*,'(A)') 'test_pcg_extra_entries'
        call seed_rnd_fixed(20261103)
        pinds        = [11, 12, 13, 14]
        ctfp%ctfflag = CTFFLAG_NO
        ctfp%smpd    = SMPD
        call sel%new(NP, .true.)
        call o%new(.true.)
        do i = 1, NP
            call o%set_euler([30.*real(i), 40. + 10.*real(i), 15.*real(i)])
            call o%set_ctfvars(ctfp)
            call o%set_shift([0.5*real(i), -0.3*real(i)])
            call sel%set_ori(i, o)
        enddo
        call op%new(BOX, SMPD, PCG_LAMBDA)
        lims2 = op%get_lims2()
        allocate(planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NP))
        do i = 1, NP
            do k = lims2(2,1), lims2(2,2)
                do h = lims2(1,1), lims2(1,2)
                    planes(h,k,i) = cmplx(ran3() - 0.5, ran3() - 0.5)
                enddo
            enddo
        enddo
        call op%kill
        ! the plain accumulation
        call accumulate_plain(b0, d0)
        ! particle 12: one extra at its committed pose (0.6 / 0.4); 14: two (0.5 / 0.3 / 0.2); 11 and 13 none
        extra%nptcls = NP
        extra%npose  = 3
        allocate(extra%loc(maxval(pinds)), extra%nextra(NP), extra%first(NP), source=0)
        allocate(extra%wc(NP), source=1.)
        allocate(extra%euler(3,3), extra%shift(2,3), extra%w(3), source=0.)
        extra%loc(pinds) = [1, 2, 3, 4]
        extra%first      = [1, 1, 2, 2]
        extra%nextra     = [0, 1, 0, 2]
        extra%wc         = [1., 0.6, 1., 0.5]
        extra%w          = [0.4, 0.3, 0.2]
        call sel%get_ori(2, o)
        extra%euler(:,1) = o%get_euler()
        extra%shift(:,1) = o%get_2Dshift()
        call sel%get_ori(4, o)
        extra%euler(:,2) = o%get_euler()
        extra%euler(:,3) = o%get_euler()
        extra%shift(:,2) = o%get_2Dshift()
        extra%shift(:,3) = o%get_2Dshift()
        call pcg_extra_entries(sel, pinds, [(1., i = 1, NP)], extra, 1., entries, wts, xptcl, xfirst)
        call assert_int(NP + 3, entries%get_noris(), 'one entry per committed and per extra pose')
        call assert_true(all(xptcl == [2, 4, 4]) .and. all(xfirst == [1, 1, 2, 2, 4]), 'the extras follow in particle order')
        call assert_true(all(abs(wts - [1., 0.6, 1., 0.5, 0.4, 0.3, 0.2]) < 1.e-6), 'the entries carry the poses'' weights')
        allocate(y_x(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), 3))
        call accumulate_entries([1, NP], b, d)
        call assert_true(maxval(abs(b - b0)) <= 1.e-4 * maxval(abs(b0)), 'extras at the committed pose reproduce the plain right-hand side')
        call assert_true(maxval(abs(d - d0)) <= 1.e-4 * maxval(abs(d0)), 'extras at the committed pose reproduce the plain density')
        call accumulate_entries([2, NP], b2, d2)
        call assert_true(maxval(abs(b2 - b)) <= 1.e-5 * maxval(abs(b)) .and. maxval(abs(d2 - d)) <= 1.e-5 * maxval(abs(d)),&
            &'two batches accumulate as one')
        ! an extra moved off the committed pose enters at its own pose
        extra%euler(:,3) = extra%euler(:,3) + [20., 15., 0.]
        call entries%kill
        call pcg_extra_entries(sel, pinds, [(1., i = 1, NP)], extra, 1., entries, wts, xptcl, xfirst)
        call accumulate_entries([1, NP], b2, d2)
        call assert_true(maxval(abs(b2 - b)) > 1.e-3 * maxval(abs(b)), 'an extra pose of its own changes the right-hand side')
        call o%kill
        call sel%kill
        call entries%kill
        call extra%kill

    contains

        subroutine accumulate_plain( bp, dp_ )
            complex, allocatable, intent(out) :: bp(:,:,:)
            real,    allocatable, intent(out) :: dp_(:,:,:)
            call op%new(BOX, SMPD, PCG_LAMBDA)
            call op%prep_particles(sel, use_ctf=.true.)
            call op%begin_accum
            call op%accumulate_batch(planes, NP, 1)
            call op%get_raw_accum(bp, dp_)
            call op%kill
        end subroutine accumulate_plain

        !> the committed and the extra entries in batches of the given size (bsz(1)) over the NP particles
        subroutine accumulate_entries( bsz, be, de )
            integer,              intent(in)  :: bsz(2)
            complex, allocatable, intent(out) :: be(:,:,:)
            real,    allocatable, intent(out) :: de(:,:,:)
            integer :: ib, ie
            call op%new(BOX, SMPD, PCG_LAMBDA)
            call op%prep_particles(entries, use_ctf=.true., weights=wts)
            call op%begin_accum
            do ib = 1, bsz(2), bsz(1)
                ie = min(bsz(2), ib + bsz(1) - 1)
                call op%accumulate_batch(planes(:,:,ib:ie), ie - ib + 1, ib)
                call accumulate_pcg_extras(op, planes(:,:,ib:ie), [ib, ie], NP, xptcl, xfirst, y_x)
            enddo
            call op%get_raw_accum(be, de)
            call op%kill
        end subroutine accumulate_entries

    end subroutine test_pcg_extra_entries

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
