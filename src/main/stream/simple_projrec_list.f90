module simple_projrec_list
include 'simple_lib.f08'
implicit none

public :: projrec_list
private
#include "simple_local_flags.inc"

interface projrec_list
    module procedure constructor
end interface projrec_list

! Convenience type to hold information about individual project files
type projrec
    type(string) :: projname             ! project file name
    integer      :: micind     = 0       ! index of micrograph in project
    integer      :: nptcls     = 0       ! # of particles
    integer      :: nptcls_sel = 0       ! # of particles (state=1)
    logical      :: included   = .false. ! whether record has been imported
end type projrec

type projrec_list
    private
    type(projrec), allocatable :: list(:)
contains
    ! constructors/lifecycle
    procedure          :: new
    procedure          :: insert
    procedure          :: replace
    procedure          :: kill
    ! overloaded assignment and append
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: append
    generic            :: operator(//)  => append
    ! checkers
    ! procedure          :: exists
    ! procedure          :: is_project
    ! ! getters
    ! procedure          :: get_projname,   get_projname_arr
    ! procedure          :: get_micind,     get_micind_arr
    ! procedure          :: get_nptcls,     get_nptcl_arr
    ! procedure          :: get_nptcls_sel, get_nptcls_sel_arr
    ! ! calculators
    ! procedure          :: get_nptcls_tot
    ! procedure          :: get_nptcls_sel_tot
end type projrec_list

! getter/setter
! call move_alloc(projrecords, old_records)

! size()

! reallocate global project
                    ! if( nprev_imports == 0 )then
                    !     call spproj_glob%os_mic%new(nmics, is_ptcl=.false.) ! first import
                    !     allocate(projrecords(nmics))
                    ! else
                    !     call spproj_glob%os_mic%reallocate(n_completed)
                    !     call move_alloc(projrecords, old_records)
                    !     allocate(projrecords(n_completed))
                    !     if( n_old > 0 ) projrecords(1:n_old) = old_records(:)
                    !     call kill_projrecords(old_records)
                    ! endif

! CHECK
! update_projects_list


contains

    ! constructors/lifecycle

    function constructor( n ) result( self )
        integer, optional, intent(in) :: n
        type(projrec_list)  :: self
        if( present(n) )then
            allocate(self%list(n))
        else
            allocate(self%list(0))
        endif
    end function constructor

    subroutine new( self, n )
        class(projrec_list), intent(inout) :: self
        integer, optional,   intent(in)    :: n
        call self%kill
        if( present(n) )then
            allocate(self%list(n))
        else
            allocate(self%list(0))
        endif
    end subroutine new

    subroutine insert( self, self_old )
        class(projrec_list), intent(inout) :: self
        class(projrec_list), intent(in)    :: self_old
        type(projrec), allocatable :: tmp_list(:)
        integer :: n, n_old
        if( .not. allocated(self_old%list) ) THROW_HARD('cannot insert unalloacted projrec_list')
        if( .not. allocated(self%list) ) allocate(self%list(0)) 
        n     = size(self%list)
        n_old = size(self_old%list)
        if( n_old <= n )then
            self%list(1:n_old) = self_old%list(:)
        else
            allocate(tmp_list(n_old))
            tmp_list = self_old%list
            call move_alloc(tmp_list, self%list)
        endif
    end subroutine insert

    subroutine replace( self, self_old )
        class(projrec_list), intent(inout) :: self
        class(projrec_list), intent(inout) :: self_old
        if( .not. allocated(self_old%list) ) THROW_HARD('cannot replace with unalloacted projrec_list')
        if( .not. allocated(self%list) ) allocate(self%list(0)) 
        call move_alloc(self_old%list, self%list)
    end subroutine replace

    subroutine kill( self )
        class(projrec_list), intent(inout) :: self
        if( allocated(self%list) )then
            call self%list(:)%projname%kill
            deallocate(self%list)
        endif
    end subroutine kill

    ! overloaded assignment and append

    subroutine assign(lhs, rhs)
        class(projrec_list), intent(inout) :: lhs
        class(projrec_list), intent(in)    :: rhs
        call lhs%kill
        if( allocated(rhs%list) ) allocate(lhs%list(size(rhs%list)), source=rhs%list)
    end subroutine assign

    function append(lhs, rhs) result(res)
        class(projrec_list), intent(in) :: lhs
        class(projrec_list), intent(in) :: rhs
        type(projrec_list) :: res
        type(projrec), allocatable :: tmp_list(:)

        ! 2 be implemented

    end function append

    !> convert a list of projects into one project
    !  previous mic/stk/ptcl2D,ptcl3D are wipped, other fields untouched
!     subroutine projrecords2proj( records, spproj )
!         class(projrecord), intent(in)    :: records(:)
!         class(sp_project), intent(inout) :: spproj
!         type(sp_project)              :: tmpproj
!         type(string) :: stack_name, projname, prev_projname
!         integer :: iptcl, fromp, ifromp, itop, jptcl, nptcls_tot
!         integer :: nrecs, nmics, nptcls, imic, micind
!         logical :: has_ptcl
!         call spproj%os_mic%kill
!         call spproj%os_stk%kill
!         call spproj%os_ptcl2D%kill
!         call spproj%os_ptcl3D%kill
!         nrecs      = size(records)
!         if( nrecs == 0 ) return 
!         nmics      = nrecs
!         nptcls_tot = sum(records(:)%nptcls)
!         has_ptcl   = nptcls_tot > 0
!         call spproj%os_mic%new(nmics,is_ptcl=.false.)
!         call spproj%os_stk%new(nmics,is_ptcl=.false.)
!         if( has_ptcl ) call spproj%os_ptcl2D%new(nptcls_tot,is_ptcl=.true.)
!         prev_projname = ''
!         jptcl = 0
!         fromp = 1
!         do imic = 1,nmics
!             ! read individual project (up to STREAM_NMOVS_SET entries)
!             projname = records(imic)%projname
!             if( projname /= prev_projname )then
!                 call tmpproj%kill
!                 call tmpproj%read_mic_stk_ptcl2D_segments(projname)
!                 prev_projname = projname
!             endif
!             ! mic
!             micind = records(imic)%micind
!             call spproj%os_mic%transfer_ori(imic, tmpproj%os_mic, micind)
!             ! stack
!             nptcls = records(imic)%nptcls
!             if( nptcls == 0 )cycle
!             call spproj%os_stk%transfer_ori(imic, tmpproj%os_stk, micind)
!             ! update stack path to absolute
!             stack_name = spproj%get_stkname(imic)
!             if( stack_name%to_char([1,1]) == '/' )then
!                 ! already absolute path, should always be the case
!             else if( stack_name%to_char([1,3]) == '../' )then
!                 stack_name = simple_abspath(stack_name)
!                 call spproj%os_stk%set(imic, 'stk', stack_name)
!             else
!                 THROW_HARD('Unexpected file path format for: '//stack_name%to_char())
!             endif
!             ! particles
!             ifromp = spproj%os_stk%get_fromp(imic)
!             itop   = spproj%os_stk%get_top(imic)
!             do iptcl = ifromp,itop
!                 jptcl = jptcl+1 ! global index
!                 call spproj%os_ptcl2D%transfer_ori(jptcl, tmpproj%os_ptcl2D, iptcl)
!                 call spproj%os_ptcl2D%set_stkind(jptcl, imic)
!             enddo
!             call spproj%os_stk%set(imic, 'fromp', fromp)
!             call spproj%os_stk%set(imic, 'top',   fromp+nptcls-1)
!             fromp = fromp + nptcls
!         enddo
!         call tmpproj%kill
!         if( has_ptcl ) spproj%os_ptcl3D = spproj%os_ptcl2D
!     end subroutine projrecords2proj

end module simple_projrec_list
