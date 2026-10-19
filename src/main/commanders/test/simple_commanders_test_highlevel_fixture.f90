!@descr: the multi-state fixture checks of simple_commanders_test_highlevel: the cache, the modes and the declared point groups of the truths
! The blob family in its three modes, each built twice (the second build must come from the cache),
! with the exact state counts and what each mode promises its project holds; the compositional and
! conformational beta-galactosidase states, each truth holding its declared point group and not the
! next larger one. Its own entry (ruling of 2026-10-17): these checks are not simulate_particles'.
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_fixture
use simple_test_gate,          only: test_gate
use simple_test_state_fixture, only: state_fixture, FIX_SEARCH, FIX_KNOWN, FIX_VIEWS, point_group_check
use simple_test_truth_metrics, only: group_operators, axis_operator
implicit none
#include "simple_local_flags.inc"

!> the declared group holds when the truth correlates with itself above this under its operators
real, parameter :: MIN_CC_GROUP = 0.99
!> the next larger group is broken when its extra operator lowers the self-correlation by at least
!! this below the declared group's: relative, so that a state broken by a small rotation (whose
!! self-correlation under the extra operator stays high) is still told from one that holds it
real, parameter :: MIN_CC_DROP  = 0.01

contains

!> simple_test_exec test=state_fixture [nthr=<n>]
module subroutine exec_test_state_fixture( self, cline )
    use simple_ui, only: make_ui
    class(commander_test_state_fixture), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    class(parameters), allocatable :: params
    type(test_gate) :: gate
    type(string)    :: cwd_saved, root
    integer :: nthr, status
    logical :: all_ok
    call make_ui
    call simple_getcwd(cwd_saved)
    allocate(params)
    call params%new(cline)
    nthr = params%nthr
    deallocate(params)
    call gate%new(filepath(cwd_saved, 'metrics.tsv'))
    root = filepath(cwd_saved, 'test_state_fixture_'//int2str(get_process_id()))
    if( dir_exists(root) ) call simple_rmdir(root)
    call simple_mkdir(root)
    call simple_chdir(root, status)
    if( status /= 0 ) THROW_HARD('could not enter '//root%to_char())
    CWD_GLOB = root%to_char()
    call check_blob_modes(nthr, gate)
    call check_betagal('compositional', 3, 60, nthr, gate)
    call check_betagal('conformational', 2, 40, nthr, gate)
    call simple_chdir(cwd_saved, status)
    if( status /= 0 ) THROW_HARD('could not return to '//cwd_saved%to_char())
    CWD_GLOB = cwd_saved%to_char()
    all_ok = gate%passed()
    call gate%kill
    if( all_ok )then
        call simple_rmdir(root)
        write(logfhandle,'(a)') 'PASS: the multi-state fixture holds its modes, its cache and its point groups'
        call simple_end('**** SIMPLE_TEST_STATE_FIXTURE NORMAL STOP ****')
    else
        THROW_HARD('TEST_STATE_FIXTURE FAILED')
    endif
end subroutine exec_test_state_fixture

!> the blob family in its three modes: each built twice (the second from the cache), the exact
!! counts of populations 50/30/20 % of 300 particles, the project each mode promises; and the blob
!! truths broken under C2
subroutine check_blob_modes( nthr, gate )
    use simple_sp_project, only: sp_project
    use simple_oris,       only: oris
    integer,         intent(in)    :: nthr
    type(test_gate), intent(inout) :: gate
    character(len=8), parameter :: MODES(3) = ['search', 'known ', 'views ']
    type(state_fixture) :: fix
    class(sp_project), allocatable :: proj
    type(oris) :: truth
    integer, allocatable :: states(:)
    integer :: imode, k, nrep, i
    real    :: ccg, ccl
    logical :: ok
    allocate(proj)
    do imode = FIX_SEARCH, FIX_VIEWS
        fix           = state_fixture()
        fix%family    = 'blob'
        fix%mode      = imode
        fix%nptcls    = 300
        fix%nviews    = 10
        fix%pops(1:3) = [0.5, 0.3, 0.2]
        fix%nthr      = nthr
        call fix%build
        call fix%build
        call gate%check('fixture blob '//trim(MODES(imode))//': the second build reads the cache', fix%reused)
        call truth%new(fix%nptcls, is_ptcl=.true.)
        call truth%read(fix%truth_doc(), [1, fix%nptcls])
        states = nint(truth%get_all('state'))
        call gate%check('fixture blob '//trim(MODES(imode))//': state counts 150, 90 and 60', &
            &count(states == 1) == 150 .and. count(states == 2) == 90 .and. count(states == 3) == 60)
        call proj%read(fix%project_fname())
        ! the records are compared only once their numbers agree
        ok = proj%os_ptcl3D%get_noris() == fix%nptcls .and. proj%os_ptcl2D%get_noris() == fix%nptcls
        if( ok )then
            ok = all(nint(proj%os_ptcl3D%get_all('state')) == 1)
            select case(imode)
                case(FIX_SEARCH)
                    ok = ok .and. all(abs(proj%os_ptcl3D%get_all('e2')) < 1.e-6)
                    ok = ok .and. maxval(abs(truth%get_all('x'))) > 0.5 .and. &
                        &maxval(abs(truth%get_all('x'))) <= fix%sherr + 1.e-4
                case(FIX_KNOWN)
                    ok = ok .and. all(abs(proj%os_ptcl3D%get_all('e2') - truth%get_all('e2')) < 1.e-3)
                    ok = ok .and. all(abs(proj%os_ptcl3D%get_all('x') - truth%get_all('x')) < 1.e-3)
                case(FIX_VIEWS)
                    nrep = fix%nptcls / fix%nviews
                    do i = 1, fix%nptcls
                        ok = ok .and. proj%os_ptcl2D%get_class(i) == (i - 1) / nrep + 1
                        ok = ok .and. abs(truth%get(i, 'e2') - truth%get(((i - 1) / nrep) * nrep + 1, 'e2')) < 1.e-3
                    enddo
            end select
        endif
        call gate%check('fixture blob '//trim(MODES(imode))//': the project holds what its mode promises', ok)
        call proj%kill
        call truth%kill
        if( imode == FIX_SEARCH )then
            do k = 1, 3
                call point_group_check(fix%truth_fname(k), group_operators('c1'), axis_operator([0.,0.,1.], 180.), &
                    &ccg, ccl)
                call gate%metric('fixture blob truth '//int2str(k)//': C2 lowers its self-correlation by', ccg - ccl, &
                    &MIN_CC_DROP, ccg - ccl >= MIN_CC_DROP)
            enddo
        endif
    enddo
    deallocate(proj)
end subroutine check_blob_modes

!> the beta-galactosidase states of kind: each truth correlates with itself above MIN_CC_GROUP
!! under its declared group, and the extra operator of the next larger one (c1 -> c2, c2 -> c4,
!! d2 -> d4: about z) lowers that by at least MIN_CC_DROP
subroutine check_betagal( kind, nstates, nptcls, nthr, gate )
    character(len=*), intent(in)    :: kind
    integer,          intent(in)    :: nstates, nptcls, nthr
    type(test_gate),  intent(inout) :: gate
    type(state_fixture) :: fix
    real    :: ccg, ccl, extra(3,3)
    integer :: k
    fix         = state_fixture()
    fix%family  = 'betagal'
    fix%kind    = kind
    fix%nstates = nstates
    fix%nptcls  = nptcls
    fix%nthr    = nthr
    call fix%build
    do k = 1, fix%nstates
        select case(trim(fix%pgrps(k)))
            case('c1');   extra = axis_operator([0.,0.,1.], 180.)
            case default; extra = axis_operator([0.,0.,1.], 90.)
        end select
        call point_group_check(fix%truth_fname(k), group_operators(trim(fix%pgrps(k))), extra, ccg, ccl)
        call gate%metric('fixture betagal '//kind//' truth '//int2str(k)//' under '//trim(fix%pgrps(k)), ccg, &
            &MIN_CC_GROUP, ccg > MIN_CC_GROUP)
        call gate%metric('fixture betagal '//kind//' truth '//int2str(k)//': the next larger group lowers it by', &
            &ccg - ccl, MIN_CC_DROP, ccg - ccl >= MIN_CC_DROP)
    enddo
end subroutine check_betagal

end submodule simple_commanders_test_highlevel_fixture
