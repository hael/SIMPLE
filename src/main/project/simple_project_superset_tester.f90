!@descr: unit tests of the abinitio3D_addon superset relation (simple_project_superset)
! Current project 20 rows / 2 stacks, frozen project 14 rows (its first stack): identity refusals
! naming the first offender, membership and floors, mask and restore.
module simple_project_superset_tester
use simple_defs,             only: STDLEN
use simple_string_utils,     only: int2str
use simple_sp_project,       only: sp_project
use simple_project_superset, only: project_superset
use simple_test_utils
implicit none
private
public :: run_all_project_superset_tests

integer, parameter :: NPTCLS     = 20 !< rows of the current project
integer, parameter :: NPTCLS_FRZ = 14 !< rows of the frozen project: the current project's first stack

contains

    subroutine run_all_project_superset_tests()
        write(*,'(A)') '**** running all project superset tests ****'
        call test_identity_accepts_superset()
        call test_identity_refusals()
        call test_membership()
        call test_mask_and_restore()
    end subroutine run_all_project_superset_tests

    ! ---- fixtures -------------------------------------------------------------

    !> NPTCLS rows in two stacks, every row active, with CTF parameters: stack 1
    !! holds rows 1-NPTCLS_FRZ, stack 2 the rows appended after the base run
    subroutine make_current( spproj )
        type(sp_project), intent(inout) :: spproj
        integer :: i, istk
        call spproj%kill
        call spproj%projinfo%new(1, is_ptcl=.false.)
        call spproj%projinfo%set(1, 'projname', 'current')
        call spproj%os_stk%new(2, is_ptcl=.false.)
        do istk = 1, 2
            call spproj%os_stk%set(istk, 'stk',   '/data/stacks/stack_'//char(48+istk)//'.mrcs')
            call spproj%os_stk%set(istk, 'fromp', merge(1, NPTCLS_FRZ + 1, istk == 1))
            call spproj%os_stk%set(istk, 'top',   merge(NPTCLS_FRZ, NPTCLS, istk == 1))
            call spproj%os_stk%set(istk, 'box',   64)
            call spproj%os_stk%set(istk, 'smpd',  1.3)
            call spproj%os_stk%set(istk, 'ctf',   'yes')
            call spproj%os_stk%set(istk, 'kv',    300.)
            call spproj%os_stk%set(istk, 'cs',    2.7)
            call spproj%os_stk%set(istk, 'fraca', 0.1)
        enddo
        call spproj%os_ptcl3D%new(NPTCLS, is_ptcl=.true.)
        do i = 1, NPTCLS
            istk = merge(1, 2, i <= NPTCLS_FRZ)
            call spproj%os_ptcl3D%set(i, 'stkind', istk)
            call spproj%os_ptcl3D%set(i, 'indstk', i - (istk-1)*NPTCLS_FRZ)
            call spproj%os_ptcl3D%set(i, 'dfx',    1.0 + 0.05*real(i))
            call spproj%os_ptcl3D%set(i, 'dfy',    1.1 + 0.05*real(i))
            call spproj%os_ptcl3D%set(i, 'angast', 15.)
            call spproj%os_ptcl3D%set_state(i, 1)
        enddo
        spproj%os_ptcl2D = spproj%os_ptcl3D
    end subroutine make_current

    !> the frozen project: the current project's first stack, rows
    !! 1-NPTCLS_FRZ; rows 1-12 selected (2D and 3D), 3D-aligned with state
    !! labels 1/2; rows 1-10 were updated, rows 11-12 never were (sampled run)
    subroutine make_frozen( cur, frozen, nstates )
        type(sp_project), intent(in)    :: cur
        type(sp_project), intent(inout) :: frozen
        integer,          intent(in)    :: nstates
        integer :: i
        frozen = cur
        call frozen%projinfo%set(1, 'projname', 'frozen')
        frozen%os_stk    = cur%os_stk%extract_subset(1, 1)
        frozen%os_ptcl2D = cur%os_ptcl2D%extract_subset(1, NPTCLS_FRZ)
        frozen%os_ptcl3D = cur%os_ptcl3D%extract_subset(1, NPTCLS_FRZ)
        do i = 1, NPTCLS_FRZ
            if( i <= 12 )then
                call frozen%os_ptcl2D%set_state(i, 1)
                call frozen%os_ptcl3D%set_state(i, 1 + mod(i, nstates))
                call frozen%os_ptcl3D%set_euler(i, [real(10*i), real(5*i), real(3*i)])
                call frozen%os_ptcl3D%set_shift(i, [0.5*real(i), -0.25*real(i)])
                call frozen%os_ptcl3D%set(i, 'corr', 0.01*real(i))
                call frozen%os_ptcl3D%set(i, 'eo',   mod(i,2))
                if( i <= 10 )then
                    call frozen%os_ptcl3D%set(i, 'updatecnt', 3)
                else
                    call frozen%os_ptcl3D%set(i, 'updatecnt', 0)
                endif
            else
                call frozen%os_ptcl2D%set_state(i, 0)
                call frozen%os_ptcl3D%set_state(i, 0)
            endif
        enddo
    end subroutine make_frozen

    !> an identity refusal of the two-state fixture that names the particle
    subroutine expect_refusal( cur, frozen, what, expected_first )
        type(sp_project), intent(inout) :: cur, frozen
        character(len=*), intent(in)    :: what
        integer,          intent(in)    :: expected_first
        type(project_superset) :: superset
        character(len=STDLEN)  :: msg
        integer :: status
        call superset%new(cur, frozen, 2, status, msg)
        call assert_true(status /= 0, what//' is refused')
        call assert_true(names_particle(msg, expected_first), what//': the refusal names the first offending particle')
        call assert_int(0, superset%get_nfrozen(), what//': a refused relation is left empty')
    end subroutine expect_refusal

    !> the refusal ends by naming particle iptcl
    logical function names_particle( msg, iptcl ) result( l_names )
        character(len=*), intent(in) :: msg
        integer,          intent(in) :: iptcl
        character(len=:), allocatable :: tail
        integer :: n
        tail    = '; first offending particle: '//int2str(iptcl)
        n       = len_trim(msg)
        l_names = .false.
        if( n >= len(tail) ) l_names = msg(n-len(tail)+1:n) == tail
    end function names_particle

    ! ---- tests ----------------------------------------------------------------

    subroutine test_identity_accepts_superset()
        type(sp_project)       :: cur, frozen
        type(project_superset) :: superset
        character(len=STDLEN)  :: msg
        integer :: status
        write(*,'(A)') 'test_identity_accepts_superset'
        call make_current(cur)
        call make_frozen(cur, frozen, 2)
        call superset%new(cur, frozen, 2, status, msg)
        call assert_int(0, status, 'a shorter frozen project on the shared indices is valid: '//trim(msg))
        call superset%kill
        call cur%kill
        call frozen%kill
    end subroutine test_identity_accepts_superset

    subroutine test_identity_refusals()
        type(sp_project)       :: cur, frozen, ref
        type(project_superset) :: superset
        character(len=STDLEN)  :: msg
        integer :: status
        write(*,'(A)') 'test_identity_refusals'
        call make_current(ref)
        ! permuted rows: rows 3 and 4 swap their images in the current project
        cur = ref
        call cur%os_ptcl3D%set(3, 'indstk', 4)
        call cur%os_ptcl3D%set(4, 'indstk', 3)
        cur%os_ptcl2D = cur%os_ptcl3D
        call make_frozen(ref, frozen, 2)
        call expect_refusal(cur, frozen, 'a row permutation', 3)
        ! a changed stack source
        cur = ref
        call cur%os_stk%set(1, 'stk', '/data/stacks/other.mrcs')
        call expect_refusal(cur, frozen, 'a changed stack source', 1)
        ! appended rows from a stack the frozen project holds would enter the union twice
        cur = ref
        call cur%os_stk%set(2, 'stk', '/data/stacks/stack_1.mrcs')
        call expect_refusal(cur, frozen, 'appended rows from a stack of the frozen project', NPTCLS_FRZ + 1)
        ! ptcl2D and ptcl3D naming different images in the current project
        cur = ref
        call cur%os_ptcl2D%set(7, 'indstk', 8)
        call expect_refusal(cur, frozen, 'a ptcl2D/ptcl3D mismatch', 7)
        ! a frozen particle deselected in the frozen project's ptcl2D
        cur = ref
        call make_frozen(ref, frozen, 2)
        call frozen%os_ptcl2D%set_state(5, 0)
        call expect_refusal(cur, frozen, 'a ptcl2D/ptcl3D selection mismatch of a frozen particle', 5)
        ! a frozen member inactive in the current project
        call make_frozen(ref, frozen, 2)
        cur = ref
        call cur%os_ptcl2D%set_state(9, 0)
        call expect_refusal(cur, frozen, 'a frozen member inactive in the current project', 9)
        ! a CTF mismatch on a frozen particle; on a cohort row it does not matter
        cur = ref
        call cur%os_ptcl3D%set(15, 'dfx', 3.3)
        call superset%new(cur, frozen, 2, status, msg)
        call assert_int(0, status, 'a changed CTF of a cohort particle is accepted: '//trim(msg))
        call cur%os_ptcl3D%set(6, 'dfx', 3.3)
        call expect_refusal(cur, frozen, 'a changed CTF of a frozen particle', 6)
        ! another optics group of a frozen particle
        cur = ref
        call cur%os_ptcl3D%set(4, 'ogid', 2)
        call expect_refusal(cur, frozen, 'another optics group of a frozen particle', 4)
        ! another sampling of a stack
        cur = ref
        call cur%os_stk%set(1, 'smpd', 1.31)
        call expect_refusal(cur, frozen, 'another stack sampling', 1)
        ! another box
        cur = ref
        call cur%os_stk%set(1, 'box', 72)
        call expect_refusal(cur, frozen, 'another stack box', 1)
        ! a current project that ends before frozen particle 10
        cur = ref
        cur%os_stk    = ref%os_stk%extract_subset(1, 1)
        cur%os_ptcl2D = ref%os_ptcl2D%extract_subset(1, 9)
        cur%os_ptcl3D = ref%os_ptcl3D%extract_subset(1, 9)
        call expect_refusal(cur, frozen, 'a frozen particle missing from the current project', 10)
        call superset%kill
        call cur%kill
        call frozen%kill
        call ref%kill
    end subroutine test_identity_refusals

    subroutine test_membership()
        type(sp_project)       :: cur, frozen, work
        type(project_superset) :: superset
        character(len=STDLEN)  :: msg
        integer :: status, i
        logical :: l_ok
        write(*,'(A)') 'test_membership'
        call make_current(cur)
        call make_frozen(cur, frozen, 2)
        call superset%new(cur, frozen, 2, status, msg)
        call assert_int(0, status, 'the membership of a valid pair is computed: '//trim(msg))
        ! rows 1-10 frozen; rows 11-12 were never updated and join rows 13-14
        ! (deselected in the frozen project) and 15-20 (appended)
        call assert_int(10, superset%get_nfrozen(),        'frozen = frozen state > 0 and updatecnt > 0')
        call assert_int(10, superset%get_ncohort(),        'cohort = current active and not frozen')
        call assert_int(2,  superset%get_nnever_updated(), 'never-updated frozen-project rows are counted')
        call assert_int(5,  superset%get_nfrozen_state(1), 'frozen particles of state 1 (even rows)')
        call assert_int(5,  superset%get_nfrozen_state(2), 'frozen particles of state 2 (odd rows)')
        call assert_false(superset%is_small_cohort(), 'a cohort as large as the frozen population raises no warning')
        ! the rows themselves, as the masked working copy shows them
        work = cur
        call superset%mask(work)
        l_ok = .true.
        do i = 1, NPTCLS
            if( i <= 10 )then
                l_ok = l_ok .and. work%os_ptcl2D%get_state(i) == 0
            else
                l_ok = l_ok .and. work%os_ptcl2D%get_state(i) > 0
            endif
        enddo
        call assert_true(l_ok, 'rows 1-10 are the frozen rows and rows 11-20 the cohort')
        ! an inactive current row belongs to neither set (one state: floor 5)
        call cur%os_ptcl2D%set_state(20, 0)
        call make_frozen(cur, frozen, 1)
        call superset%new(cur, frozen, 1, status, msg)
        call assert_int(0, status, 'a cohort of 9 is above the floor of one state: '//trim(msg))
        call assert_int(9, superset%get_ncohort(), 'an inactive current particle is not in the cohort')
        ! the per-state floor: 2 states need at least 10 cohort particles
        call make_frozen(cur, frozen, 2)
        call superset%new(cur, frozen, 2, status, msg)
        call assert_true(status /= 0, 'a cohort below 5 particles per inherited state is refused')
        call assert_int(0, superset%get_ncohort(), 'a refused relation is left empty')
        call superset%new(cur, frozen, 1, status, msg)
        call assert_true(status /= 0, 'a frozen state label above nstates is refused')
        ! an empty cohort
        do i = 11, NPTCLS
            call cur%os_ptcl2D%set_state(i, 0)
        enddo
        call superset%new(cur, frozen, 2, status, msg)
        call assert_true(status /= 0, 'an empty cohort is refused')
        ! an empty inherited state (every frozen particle in state 1 of 2)
        call make_current(cur)
        call make_frozen(cur, frozen, 1)
        call superset%new(cur, frozen, 2, status, msg)
        call assert_true(status /= 0, 'an empty inherited state is refused')
        call superset%new(cur, frozen, 1, status, msg)
        call assert_int(0, status, 'the same frozen project is valid for its own single state')
        call superset%kill
        call cur%kill
        call frozen%kill
        call work%kill
    end subroutine test_membership

    !> one inherited state, so that a deselected row leaves the cohort above its floor
    subroutine test_mask_and_restore()
        type(sp_project)       :: cur, frozen, work
        type(project_superset) :: superset
        character(len=STDLEN)  :: msg
        real    :: e_work(3), e_frz(3), sh_work(2), sh_frz(2)
        integer :: status, i
        logical :: l_ok
        write(*,'(A)') 'test_mask_and_restore'
        call make_current(cur)
        call cur%os_ptcl2D%set_state(20, 0)   ! a row the user deselected stays deselected
        call make_frozen(cur, frozen, 1)
        call superset%new(cur, frozen, 1, status, msg)
        call assert_int(0, status, 'the relation to mask is valid: '//trim(msg))
        work = cur
        call superset%mask(work)
        l_ok = .true.
        do i = 1, NPTCLS
            if( i <= 10 )then
                l_ok = l_ok .and. work%os_ptcl2D%get_state(i) == 0 .and. work%os_ptcl3D%get_state(i) == 0
            else
                l_ok = l_ok .and. work%os_ptcl2D%get_state(i) == cur%os_ptcl2D%get_state(i)
            endif
        enddo
        call assert_true(l_ok, 'masking zeroes the frozen rows in ptcl2D and ptcl3D and nothing else')
        call assert_int(9, work%count_state_gt_zero(), 'the masked working copy counts the cohort alone')
        call superset%validate_cohort_states(work, status, msg)
        call assert_int(0, status, 'the labelled cohort holds the floor of its state: '//trim(msg))
        ! ptcl3D rows 11-20 hold the state (row 20 is deselected in ptcl2D only)
        do i = 11, 16
            call work%os_ptcl3D%set_state(i, 0)
        enddo
        call superset%validate_cohort_states(work, status, msg)
        call assert_true(status /= 0, 'a state below the cohort floor after labelling is refused')
        do i = 11, 16
            call work%os_ptcl3D%set_state(i, 1)
        enddo
        ! the refinement moves the cohort; the frozen rows come back from the frozen project
        do i = 11, NPTCLS
            call work%os_ptcl3D%set_euler(i, [1., 2., 3.])
        enddo
        call superset%restore(work, frozen)
        l_ok = .true.
        do i = 1, NPTCLS
            l_ok = l_ok .and. work%os_ptcl2D%get_state(i) == cur%os_ptcl2D%get_state(i)
            if( i > 10 ) cycle
            e_work  = work%os_ptcl3D%get_euler(i)
            e_frz   = frozen%os_ptcl3D%get_euler(i)
            sh_work = work%os_ptcl3D%get_2Dshift(i)
            sh_frz  = frozen%os_ptcl3D%get_2Dshift(i)
            l_ok = l_ok .and. all(abs(e_work - e_frz) < 1.e-4) .and. all(abs(sh_work - sh_frz) < 1.e-6) .and. &
                &work%os_ptcl3D%get_state(i) == frozen%os_ptcl3D%get_state(i) .and. &
                &work%os_ptcl3D%get_updatecnt(i) == frozen%os_ptcl3D%get_updatecnt(i) .and. &
                &work%os_ptcl3D%get_eo(i) == frozen%os_ptcl3D%get_eo(i) .and. &
                &abs(work%os_ptcl3D%get(i, 'corr') - frozen%os_ptcl3D%get(i, 'corr')) < 1.e-6
        enddo
        call assert_true(l_ok, 'restoration brings back every frozen 3D record and every saved ptcl2D state')
        e_work = work%os_ptcl3D%get_euler(15)
        call assert_true(all(abs(e_work - [1., 2., 3.]) < 1.e-4), 'restoration leaves cohort poses as refined')
        call superset%kill
        call cur%kill
        call frozen%kill
        call work%kill
    end subroutine test_mask_and_restore

end module simple_project_superset_tester
