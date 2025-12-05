module simple_projrec_list
include 'simple_lib.f08'
implicit none

public :: projrec, projrec_list
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
    procedure          :: copy_from
    procedure          :: replace
    procedure          :: kill
    ! overloaded assignment and append
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: append
    generic            :: operator(//)  => append
    ! checkers
    procedure          :: is_allocated
    procedure, private :: is_project_1, is_project_2
    generic            :: is_project => is_project_1, is_project_2
    procedure          :: len
    ! setter
    procedure          :: push
    procedure          :: set
    ! getters
    procedure          :: get_projname,   get_projname_arr
    procedure          :: get_micind,     get_micind_arr
    procedure          :: get_nptcls,     get_nptcls_arr
    procedure          :: get_nptcls_sel, get_nptcls_sel_arr
    procedure          :: get_nptcls_tot
    procedure          :: get_nptcls_sel_tot
end type projrec_list

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

    subroutine copy_from( self, self_old )
        class(projrec_list), intent(inout) :: self
        class(projrec_list), intent(in)    :: self_old
        type(projrec), allocatable :: tmp_list(:)
        integer :: n, n_old
        if( .not. allocated(self_old%list) ) return
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
    end subroutine copy_from

    subroutine replace( self, self_old )
        class(projrec_list), intent(inout) :: self
        class(projrec_list), intent(inout) :: self_old
        if( .not. allocated(self_old%list) ) return
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
        type(projrec), allocatable :: tmp_list(:)
        type(projrec_list) :: res
        integer :: n_lhs, n_rhs
        n_lhs = 0
        if( allocated(lhs%list) ) n_lhs = size(lhs%list) 
        n_rhs = 0
        if( allocated(rhs%list) ) n_rhs = size(rhs%list)
        if( n_lhs > 0 .and. n_rhs > 0  )then
            allocate(tmp_list(n_lhs + n_rhs))
            tmp_list(1:n_lhs)    = lhs%list(:)
            tmp_list(n_lhs + 1:) = rhs%list(:)
            call move_alloc(tmp_list, res%list)
        else if( n_lhs > 0 .and. n_rhs == 0  )then
            allocate(res%list(n_lhs), source=lhs%list)
        else if( n_lhs == 0 .and. n_rhs > 0 )then
            allocate(res%list(n_rhs), source=rhs%list)
        else
            allocate(res%list(0))
        endif
    end function append

    ! checkers

    pure logical function is_allocated( self )
        class(projrec_list), intent(in) :: self
        is_allocated = allocated(self%list)
    end function is_allocated

    pure logical function is_project_1(self, ind, projname)
        class(projrec_list), intent(in) :: self
        integer,             intent(in) :: ind
        class(string),       intent(in) :: projname
        if (.not. allocated(self%list)) then
            is_project_1 = .false.
        else if (ind < 1 .or. ind > size(self%list)) then
            is_project_1 = .false.
        else
            is_project_1 = self%list(ind)%projname .eq. projname
        end if
    end function is_project_1

    pure logical function is_project_2(self, ind, projname)
        class(projrec_list), intent(in) :: self
        integer,             intent(in) :: ind
        character(len=*),    intent(in) :: projname
        if (.not. allocated(self%list)) then
            is_project_2 = .false.
        else if (ind < 1 .or. ind > size(self%list)) then
            is_project_2 = .false.
        else
            is_project_2 = self%list(ind)%projname .eq. projname
        end if
    end function is_project_2

    pure function len( self ) result ( sz )
        class(projrec_list), intent(in) :: self
        integer :: sz
        if (allocated(self%list)) then
            sz = size(self%list)
        else
            sz = 0
        end if
    end function len

    ! setters

    subroutine push( self, record )
        class(projrec_list), intent(inout) :: self
        type(projrec),       intent(in) :: record
        type(projrec_list) :: one_elem
        one_elem = projrec_list(1)
        one_elem%list(1) = record
        self = self // one_elem
    end subroutine push

    subroutine set( self, ind, record )
        class(projrec_list), intent(inout) :: self
        integer,             intent(in)    :: ind
        type(projrec),       intent(in)    :: record
        if (.not. allocated(self%list) .and. ind == 1 ) then
            call self%push(record)
        else if ( ind == size(self%list) + 1 )then
            call self%push(record)            
        else if( ind > 0 .and. ind <= size(self%list) )then
            self%list(ind) = record
        else
            THROW_HARD('index ind out of bounds')
        end if
    end subroutine set

    ! getters

    type(string) function get_projname( self, ind )
        class(projrec_list), intent(in) :: self
        integer,             intent(in) :: ind
        if (.not. allocated(self%list)) then
            get_projname = ''
        else if (ind < 1 .or. ind > size(self%list)) then
            get_projname = ''
        else
            get_projname = self%list(ind)%projname
        end if
    end function get_projname

    pure integer function get_micind( self, ind )
        class(projrec_list), intent(in) :: self
        integer,             intent(in) :: ind
        if (.not. allocated(self%list)) then
            get_micind = 0
        else if (ind < 1 .or. ind > size(self%list)) then
            get_micind = 0
        else
            get_micind = self%list(ind)%micind
        end if
    end function get_micind

    pure integer function get_nptcls( self, ind )
        class(projrec_list), intent(in) :: self
        integer,             intent(in) :: ind
        if (.not. allocated(self%list)) then
            get_nptcls = 0
        else if (ind < 1 .or. ind > size(self%list)) then
            get_nptcls = 0
        else
            get_nptcls = self%list(ind)%nptcls
        end if
    end function get_nptcls

    pure integer function get_nptcls_sel( self, ind )
        class(projrec_list), intent(in) :: self
        integer,             intent(in) :: ind
        if (.not. allocated(self%list)) then
            get_nptcls_sel = 0
        else if (ind < 1 .or. ind > size(self%list)) then
            get_nptcls_sel = 0
        else
            get_nptcls_sel = self%list(ind)%nptcls_sel
        end if
    end function get_nptcls_sel

    pure function get_projname_arr( self ) result( pnames )
        class(projrec_list), intent(in) :: self
        type(string), allocatable :: pnames(:)
        integer :: n
        n = 0
        if( allocated(self%list) ) n = size(self%list)
        if( n > 0 )then
            pnames = self%list(:)%projname
        else    
            allocate(pnames(0))
        endif
    end function get_projname_arr

    pure function get_micind_arr( self ) result( micinds )
        class(projrec_list), intent(in) :: self
        integer, allocatable :: micinds(:)
        integer :: n
        n = 0
        if( allocated(self%list) ) n = size(self%list)
        if( n > 0 )then
            micinds = self%list(:)%micind
        else    
            allocate(micinds(0))
        endif
    end function get_micind_arr

    pure function get_nptcls_arr( self ) result( nptcls )
        class(projrec_list), intent(in) :: self
        integer, allocatable :: nptcls(:)
        integer :: n
        n = 0
        if( allocated(self%list) ) n = size(self%list)
        if( n > 0 )then
            nptcls = self%list(:)%nptcls
        else    
            allocate(nptcls(0))
        endif
    end function get_nptcls_arr

    pure function get_nptcls_sel_arr( self ) result( nptcls_sel )
        class(projrec_list), intent(in) :: self
        integer, allocatable :: nptcls_sel(:)
        integer :: n
        n = 0
        if( allocated(self%list) ) n = size(self%list)
        if( n > 0 )then
            nptcls_sel = self%list(:)%nptcls_sel
        else    
            allocate(nptcls_sel(0))
        endif
    end function get_nptcls_sel_arr

    integer function get_nptcls_tot(self)
        class(projrec_list), intent(in) :: self
        get_nptcls_tot = 0
        if( allocated(self%list) )then
            if( size(self%list) > 0 ) get_nptcls_tot = sum(self%list(:)%nptcls)
        end if
    end function get_nptcls_tot

    integer function get_nptcls_sel_tot(self)
        class(projrec_list), intent(in) :: self
        get_nptcls_sel_tot = 0
         if( allocated(self%list) )then
            if( size(self%list) > 0 ) get_nptcls_sel_tot = sum(self%list(:)%nptcls_sel)            
        end if
    end function get_nptcls_sel_tot

    ! THIS ONE BELONGS IN THE SP_PROJECT AREA
    ! SHOULD BE A METHOD THAT IS CONSTRUCTING A SP_PROJECT INSTANCE FROM THE RECORDS PROVIDED BY THIS CLASS
    !> convert a list of projects into one project
    !  previous mic/stk/ptcl2D,ptcl3D are wiped, other fields untouched
    ! subroutine projrecords2proj( records, spproj )
    !     class(projrecord), intent(in)    :: records(:)
    !     class(sp_project), intent(inout) :: spproj
    !     type(sp_project)              :: tmpproj
    !     type(string) :: stack_name, projname, prev_projname
    !     integer :: iptcl, fromp, ifromp, itop, jptcl, nptcls_tot
    !     integer :: nrecs, nmics, nptcls, imic, micind
    !     logical :: has_ptcl
    !     call spproj%os_mic%kill
    !     call spproj%os_stk%kill
    !     call spproj%os_ptcl2D%kill
    !     call spproj%os_ptcl3D%kill
    !     nrecs      = size(records)
    !     if( nrecs == 0 ) return 
    !     nmics      = nrecs
    !     nptcls_tot = sum(records(:)%nptcls)
    !     has_ptcl   = nptcls_tot > 0
    !     call spproj%os_mic%new(nmics,is_ptcl=.false.)
    !     call spproj%os_stk%new(nmics,is_ptcl=.false.)
    !     if( has_ptcl ) call spproj%os_ptcl2D%new(nptcls_tot,is_ptcl=.true.)
    !     prev_projname = ''
    !     jptcl = 0
    !     fromp = 1
    !     do imic = 1,nmics
    !         ! read individual project (up to STREAM_NMOVS_SET entries)
    !         projname = records(imic)%projname
    !         if( projname /= prev_projname )then
    !             call tmpproj%kill
    !             call tmpproj%read_mic_stk_ptcl2D_segments(projname)
    !             prev_projname = projname
    !         endif
    !         ! mic
    !         micind = records(imic)%micind
    !         call spproj%os_mic%transfer_ori(imic, tmpproj%os_mic, micind)
    !         ! stack
    !         nptcls = records(imic)%nptcls
    !         if( nptcls == 0 )cycle
    !         call spproj%os_stk%transfer_ori(imic, tmpproj%os_stk, micind)
    !         ! update stack path to absolute
    !         stack_name = spproj%get_stkname(imic)
    !         if( stack_name%to_char([1,1]) == '/' )then
    !             ! already absolute path, should always be the case
    !         else if( stack_name%to_char([1,3]) == '../' )then
    !             stack_name = simple_abspath(stack_name)
    !             call spproj%os_stk%set(imic, 'stk', stack_name)
    !         else
    !             THROW_HARD('Unexpected file path format for: '//stack_name%to_char())
    !         endif
    !         ! particles
    !         ifromp = spproj%os_stk%get_fromp(imic)
    !         itop   = spproj%os_stk%get_top(imic)
    !         do iptcl = ifromp,itop
    !             jptcl = jptcl+1 ! global index
    !             call spproj%os_ptcl2D%transfer_ori(jptcl, tmpproj%os_ptcl2D, iptcl)
    !             call spproj%os_ptcl2D%set_stkind(jptcl, imic)
    !         enddo
    !         call spproj%os_stk%set(imic, 'fromp', fromp)
    !         call spproj%os_stk%set(imic, 'top',   fromp+nptcls-1)
    !         fromp = fromp + nptcls
    !     enddo
    !     call tmpproj%kill
    !     if( has_ptcl ) spproj%os_ptcl3D = spproj%os_ptcl2D
    ! end subroutine projrecords2proj

end module simple_projrec_list
