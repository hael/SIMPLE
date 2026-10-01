!@descr: abinitio3D_addon superset relation of a current and a frozen project: identity, frozen/cohort membership, masking and restoration
! Current and frozen projects share row indices; shared rows must name the same image (and CTF/optics
! for frozen rows), appended rows come from new stacks. frozen = frozen ptcl3D state>0 & updatecnt>0;
! current-inactive frozen rows are retired. mask zeroes frozen/retired rows; restore keeps retired rows inactive.
! Contract: doc/policies/3D/abinitio3D_addon_policy.md sec. 4.
module simple_project_superset
use simple_core_module_api
use simple_sp_project, only: sp_project
implicit none

public :: project_superset, COHORT_WARN_FRAC
private
#include "simple_local_flags.inc"

integer, parameter :: MIN_COHORT_STATE_POP = 5    !< hard floor: cohort particles per inherited state
real,    parameter :: COHORT_WARN_FRAC     = 0.05 !< warn below this fraction of the frozen population

!> the validated superset relation of a current and a frozen project: the
!! active frozen rows, retired rows, the cohort, and saved current ptcl2D states
type :: project_superset
    private
    integer :: nrows = 0, nrows_frozen = 0, nstates = 0 !< row counts of the current and the frozen project
    integer :: nfrozen = 0, ncohort = 0, nnever_updated = 0, nretired = 0
    logical, allocatable :: l_frozen(:), l_retired(:)
    integer, allocatable :: nfrozen_state(:)
    integer, allocatable :: saved_state2D(:)   !< current ptcl2D states before masking
  contains
    procedure :: new
    procedure :: validate_cohort_states
    procedure :: mask
    procedure :: retire_from_frozen
    procedure :: restore
    procedure :: get_nfrozen
    procedure :: get_ncohort
    procedure :: get_nnever_updated
    procedure :: get_nretired
    procedure :: get_nfrozen_state
    procedure :: is_small_cohort
    procedure :: kill
end type project_superset

contains

    !> Validate the relation and define the membership, with the refusals
    !! that precede any write: a row-wise identity defect (the message names
    !! the first offending particle), an empty inherited state, an empty
    !! cohort, and a cohort that gives an inherited state fewer than
    !! MIN_COHORT_STATE_POP particles. A cohort below COHORT_WARN_FRAC of the
    !! frozen population is allowed (is_small_cohort, warned by the caller).
    !! On refusal status /= 0 and the object is left empty.
    subroutine new( self, cur, frozen, nstates, status, msg )
        class(project_superset), intent(inout) :: self
        class(sp_project),       intent(inout) :: cur, frozen
        integer,                 intent(in)    :: nstates
        integer,                 intent(out)   :: status
        character(len=*),        intent(out)   :: msg
        logical, allocatable :: l_cohort(:)
        integer :: n, i, s, nper_state
        call self%kill
        call validate_identity(cur, frozen, status, msg)
        if( status /= 0 ) return
        status = 1
        n = cur%os_ptcl3D%get_noris()
        self%nrows        = n
        self%nrows_frozen = frozen%os_ptcl3D%get_noris()
        self%nstates      = nstates
        allocate(self%l_frozen(n), self%l_retired(n), l_cohort(n), source=.false.)
        allocate(self%nfrozen_state(nstates), source=0)
        do i = 1, n
            if( i <= self%nrows_frozen ) self%l_frozen(i) = is_frozen_row(frozen, i)
            if( self%l_frozen(i) .and. cur%os_ptcl2D%get_state(i) <= 0 )then
                self%l_frozen(i)  = .false.
                self%l_retired(i) = .true.
                self%nretired     = self%nretired + 1
            endif
            if( self%l_frozen(i) )then
                s = frozen%os_ptcl3D%get_state(i)
                if( s > nstates )then
                    write(msg,'(A,I0,A,I0)') 'frozen particle ', i, ' carries a state label above nstates: ', s
                    call self%kill
                    return
                endif
                self%nfrozen_state(s) = self%nfrozen_state(s) + 1
            else
                l_cohort(i) = cur%os_ptcl2D%get_state(i) > 0
                if( l_cohort(i) .and. i <= self%nrows_frozen )then
                    if( frozen%os_ptcl3D%get_state(i) > 0 ) self%nnever_updated = self%nnever_updated + 1
                endif
            endif
        enddo
        self%nfrozen = count(self%l_frozen)
        self%ncohort = count(l_cohort)
        if( self%nfrozen < 1 )then
            msg = 'the frozen project has no frozen particles (state > 0 and updatecnt > 0)'
        else if( any(self%nfrozen_state < 1) )then
            msg = 'an inherited state has no frozen particles; the state layout must be contiguous 1..nstates'
        else if( self%ncohort < 1 )then
            msg = 'the current project adds no particles to the frozen solution (empty cohort)'
        else
            ! balanced labelling gives every inherited state ncohort/nstates
            ! cohort particles or one more
            nper_state = self%ncohort / nstates
            if( nper_state < MIN_COHORT_STATE_POP )then
                write(msg,'(A,I0,A,I0,A,I0,A)') 'the cohort (', self%ncohort, ' particles) gives an inherited state ', &
                    &nper_state, ', below the floor of ', MIN_COHORT_STATE_POP, ' per state'
            else
                status = 0
            endif
        endif
        if( status /= 0 ) call self%kill
    end subroutine new

    !> After the cohort is labelled in the masked working copy (frozen rows at
    !! state 0), every inherited state holds at least MIN_COHORT_STATE_POP
    !! cohort particles
    subroutine validate_cohort_states( self, spproj, status, msg )
        class(project_superset), intent(in)    :: self
        class(sp_project),       intent(inout) :: spproj
        integer,                 intent(out)   :: status
        character(len=*),        intent(out)   :: msg
        integer :: s, pop
        status = 1
        msg    = ''
        if( .not. allocated(self%saved_state2D) ) THROW_HARD('the cohort is validated only after masking')
        do s = 1, self%nstates
            pop = spproj%os_ptcl3D%get_pop(s, 'state')
            if( pop < MIN_COHORT_STATE_POP )then
                write(msg,'(A,I0,A,I0,A,I0)') 'inherited state ', s, ' holds ', pop, &
                    &' cohort particles, below the floor of ', MIN_COHORT_STATE_POP
                return
            endif
        enddo
        status = 0
    end subroutine validate_cohort_states

    !> Save the working copy's ptcl2D states and set state 0 in ptcl2D and
    !! ptcl3D for every frozen row, so every counting, sampling and labelling
    !! routine of the established workflow sees the cohort alone
    subroutine mask( self, spproj )
        class(project_superset), intent(inout) :: self
        class(sp_project),       intent(inout) :: spproj
        integer :: i
        if( .not. allocated(self%l_frozen) ) THROW_HARD('the superset relation was never established')
        if( spproj%os_ptcl2D%get_noris() /= self%nrows .or. spproj%os_ptcl3D%get_noris() /= self%nrows ) &
            &THROW_HARD('frozen-row mask does not match the working project')
        if( allocated(self%saved_state2D) ) deallocate(self%saved_state2D)
        allocate(self%saved_state2D(self%nrows))
        do i = 1, self%nrows
            self%saved_state2D(i) = spproj%os_ptcl2D%get_state(i)
            if( self%l_frozen(i) .or. self%l_retired(i) )then
                call spproj%os_ptcl2D%set_state(i, 0)
                call spproj%os_ptcl3D%set_state(i, 0)
            endif
        enddo
    end subroutine mask

    !> Remove current-project rejections from the private frozen copy before
    !! its accumulators are generated; row identity and numbering stay intact.
    subroutine retire_from_frozen( self, frozen )
        class(project_superset), intent(in)    :: self
        class(sp_project),       intent(inout) :: frozen
        integer :: i
        if( .not. allocated(self%l_retired) ) THROW_HARD('the superset relation was never established')
        if( frozen%os_ptcl2D%get_noris() /= self%nrows_frozen .or. &
            &frozen%os_ptcl3D%get_noris() /= self%nrows_frozen ) THROW_HARD('retired-row mask does not match the frozen project')
        do i = 1, self%nrows_frozen
            if( .not. self%l_retired(i) ) cycle
            call frozen%os_ptcl2D%set_state(i, 0)
            call frozen%os_ptcl3D%set_state(i, 0)
        enddo
    end subroutine retire_from_frozen

    !> Restore active frozen rows and every saved ptcl2D state; retired rows
    !! stay inactive and cohort 3D records stay as refined.
    subroutine restore( self, spproj, frozen )
        class(project_superset), intent(in)    :: self
        class(sp_project),       intent(inout) :: spproj
        class(sp_project),       intent(in)    :: frozen
        integer :: i
        if( .not. allocated(self%saved_state2D) ) THROW_HARD('frozen rows were never masked')
        if( spproj%os_ptcl3D%get_noris() /= self%nrows .or. frozen%os_ptcl3D%get_noris() /= self%nrows_frozen ) &
            &THROW_HARD('frozen-row restore does not match the working project')
        do i = 1, self%nrows
            call spproj%os_ptcl2D%set_state(i, self%saved_state2D(i))
            if( self%l_retired(i) )then
                call spproj%os_ptcl3D%set_state(i, 0)
                cycle
            endif
            if( .not. self%l_frozen(i) ) cycle
            call spproj%os_ptcl3D%transfer_3Dparams(i, frozen%os_ptcl3D, i)
            call spproj%os_ptcl3D%set_state(i, frozen%os_ptcl3D%get_state(i))
        enddo
    end subroutine restore

    integer function get_nfrozen( self ) result( n )
        class(project_superset), intent(in) :: self
        n = self%nfrozen
    end function get_nfrozen

    integer function get_ncohort( self ) result( n )
        class(project_superset), intent(in) :: self
        n = self%ncohort
    end function get_ncohort

    !> cohort rows the frozen project had selected but never updated
    integer function get_nnever_updated( self ) result( n )
        class(project_superset), intent(in) :: self
        n = self%nnever_updated
    end function get_nnever_updated

    integer function get_nretired( self ) result( n )
        class(project_superset), intent(in) :: self
        n = self%nretired
    end function get_nretired

    !> frozen particles of one inherited state
    integer function get_nfrozen_state( self, state ) result( n )
        class(project_superset), intent(in) :: self
        integer,                 intent(in) :: state
        if( state < 1 .or. state > self%nstates ) THROW_HARD('state is outside the inherited state layout')
        n = self%nfrozen_state(state)
    end function get_nfrozen_state

    !> a cohort below COHORT_WARN_FRAC of the frozen population
    logical function is_small_cohort( self ) result( l_small )
        class(project_superset), intent(in) :: self
        l_small = real(self%ncohort) < COHORT_WARN_FRAC * real(self%nfrozen)
    end function is_small_cohort

    subroutine kill( self )
        class(project_superset), intent(inout) :: self
        self%nrows = 0; self%nrows_frozen = 0; self%nstates = 0
        self%nfrozen = 0; self%ncohort = 0; self%nnever_updated = 0; self%nretired = 0
        if( allocated(self%l_frozen)      ) deallocate(self%l_frozen)
        if( allocated(self%l_retired)     ) deallocate(self%l_retired)
        if( allocated(self%nfrozen_state) ) deallocate(self%nfrozen_state)
        if( allocated(self%saved_state2D) ) deallocate(self%saved_state2D)
    end subroutine kill

    ! PRIVATE HELPERS

    logical function is_frozen_row( frozen, i ) result( l_frozen )
        class(sp_project), intent(in) :: frozen
        integer,           intent(in) :: i
        l_frozen = frozen%os_ptcl3D%get_state(i) > 0
        if( l_frozen ) l_frozen = frozen%os_ptcl3D%get_updatecnt(i) > 0
    end function is_frozen_row

    !> Row-wise physical identity of the current and the frozen project on the
    !! rows both hold, and the superset relation: every frozen particle lies
    !! within the current project's rows, and appended rows (beyond the frozen
    !! project's last row) come from stacks the frozen project does not hold.
    !! status /= 0 names the defect and, for a row defect, the first offending
    !! particle index.
    subroutine validate_identity( cur, frozen, status, msg )
        class(sp_project), intent(inout) :: cur, frozen
        integer,           intent(out)   :: status
        character(len=*),  intent(out)   :: msg
        type(ctfparams) :: ctf_cur, ctf_frz
        integer         :: n, nf, i
        status = 1
        msg    = ''
        n = cur%os_ptcl3D%get_noris()
        if( n < 1 )then
            msg = 'the current project has no particles'
            return
        endif
        if( cur%os_ptcl2D%get_noris() /= n )then
            msg = 'the current project ptcl2D and ptcl3D segments differ in length'
            return
        endif
        nf = frozen%os_ptcl3D%get_noris()
        if( frozen%os_ptcl2D%get_noris() /= nf )then
            msg = 'the frozen project ptcl2D and ptcl3D segments differ in length'
            return
        endif
        do i = 1, n
            if( i > nf )then
                ! an appended row: consistent in the current project
                if( image_id(cur, 'ptcl2D', i) /= image_id(cur, 'ptcl3D', i) )then
                    msg = 'the current project ptcl2D and ptcl3D rows name different images'
                    call name_particle(i)
                    return
                endif
                cycle
            endif
            ! the same physical image in every segment of both projects
            if( .not. same_image(cur, 'ptcl3D', frozen, 'ptcl3D', i) )then
                msg = 'ptcl3D rows name different images (permuted rows or a changed stack source)'
                call name_particle(i)
                return
            endif
            if( .not. same_image(cur, 'ptcl2D', frozen, 'ptcl2D', i) )then
                msg = 'ptcl2D rows name different images (permuted rows or a changed stack source)'
                call name_particle(i)
                return
            endif
            if( image_id(cur, 'ptcl2D', i) /= image_id(cur, 'ptcl3D', i) )then
                msg = 'the current project ptcl2D and ptcl3D rows name different images'
                call name_particle(i)
                return
            endif
            if( image_id(frozen, 'ptcl2D', i) /= image_id(frozen, 'ptcl3D', i) )then
                msg = 'the frozen project ptcl2D and ptcl3D rows name different images'
                call name_particle(i)
                return
            endif
            if( .not. is_frozen_row(frozen, i) ) cycle
            ! a frozen member must be active where the fresh-start selection is made
            if( frozen%os_ptcl2D%get_state(i) <= 0 )then
                msg = 'a frozen particle is deselected in the frozen project ptcl2D (ptcl2D/ptcl3D selection mismatch)'
                call name_particle(i)
                return
            endif
            ! the optics and CTF identity that reproduces the frozen contribution
            ctf_cur = cur%get_ctfparams('ptcl3D', i)
            ctf_frz = frozen%get_ctfparams('ptcl3D', i)
            if( .not. same_ctf(ctf_cur, ctf_frz) )then
                msg = 'a frozen particle has other CTF or optics parameters in the current project'
                call name_particle(i)
                return
            endif
            if( cur%os_ptcl3D%isthere(i, 'ogid') .or. frozen%os_ptcl3D%isthere(i, 'ogid') )then
                if( cur%os_ptcl3D%get_int(i, 'ogid') /= frozen%os_ptcl3D%get_int(i, 'ogid') )then
                    msg = 'a frozen particle belongs to another optics group in the current project'
                    call name_particle(i)
                    return
                endif
            endif
        enddo
        ! frozen-project rows past the current project's last row hold no frozen particle
        do i = n + 1, nf
            if( is_frozen_row(frozen, i) )then
                msg = 'a frozen particle is missing from the current project (its row lies past the last)'
                call name_particle(i)
                return
            endif
        enddo
        if( n > nf )then
            call check_appended_stacks
            if( len_trim(msg) > 0 ) return
        endif
        status = 0

    contains

        subroutine name_particle( iptcl )
            integer, intent(in) :: iptcl
            msg = trim(msg)//'; first offending particle: '//int2str(iptcl)
        end subroutine name_particle

        !> appended rows come from stacks the frozen project does not hold: an
        !! appended copy of a frozen-project image would enter the union twice
        subroutine check_appended_stacks
            type(string)         :: stk_cur, stk_frz
            logical, allocatable :: l_app(:)
            integer              :: j, k, istk, stkind, ind
            allocate(l_app(cur%os_stk%get_noris()), source=.false.)
            do j = nf + 1, n
                call cur%map_ptcl_ind2stk_ind('ptcl3D', j, stkind, ind)
                l_app(stkind) = .true.
            enddo
            do k = 1, size(l_app)
                if( .not. l_app(k) ) cycle
                stk_cur = cur%os_stk%get_str(k, 'stk')
                do istk = 1, frozen%os_stk%get_noris()
                    stk_frz = frozen%os_stk%get_str(istk, 'stk')
                    if( stk_cur%to_char() /= stk_frz%to_char() ) cycle
                    msg = 'appended particles come from a stack the frozen project holds: '//stk_cur%to_char()
                    do j = nf + 1, n
                        call cur%map_ptcl_ind2stk_ind('ptcl3D', j, stkind, ind)
                        if( stkind == k ) exit
                    enddo
                    call name_particle(j)
                    call stk_cur%kill
                    call stk_frz%kill
                    return
                enddo
            enddo
            call stk_cur%kill
            call stk_frz%kill
        end subroutine check_appended_stacks

    end subroutine validate_identity

    !> the same stack file, physical image and stack geometry for row i
    logical function same_image( a, seg_a, b, seg_b, i ) result( l_same )
        class(sp_project), intent(inout) :: a, b
        character(len=*),  intent(in)    :: seg_a, seg_b
        integer,           intent(in)    :: i
        l_same = image_id(a, seg_a, i) == image_id(b, seg_b, i)
    end function same_image

    !> row i's physical identity: stack file, image index, stack box and sampling
    function image_id( p, seg, i ) result( id )
        class(sp_project), intent(inout) :: p
        character(len=*),  intent(in)    :: seg
        integer,           intent(in)    :: i
        character(len=:), allocatable :: id
        type(string) :: stk
        integer :: stkind, ind
        call p%map_ptcl_ind2stk_ind(seg, i, stkind, ind)
        stk = p%os_stk%get_str(stkind, 'stk')
        id  = trim(stk%to_char())//'|'//int2str(ind)//'|'//int2str(p%os_stk%get_int(stkind, 'box'))// &
            &'|'//trim(real2str(p%os_stk%get(stkind, 'smpd')))
        call stk%kill
    end function image_id

    logical function same_ctf( a, b ) result( l_same )
        type(ctfparams), intent(in) :: a, b
        l_same = a%ctfflag == b%ctfflag .and. a%smpd == b%smpd .and. a%kv == b%kv .and. a%cs == b%cs .and. &
            &a%fraca == b%fraca .and. a%dfx == b%dfx .and. a%dfy == b%dfy .and. a%angast == b%angast .and. &
            &a%phshift == b%phshift
    end function same_ctf

end module simple_project_superset
