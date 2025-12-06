module simple_poly_rec_list
include 'simple_lib.f08'
use simple_linked_list
implicit none

public :: project_rec
public :: process_rec
public :: chunk_rec
public :: poly_rec_list
private
#include "simple_local_flags.inc"

! Convenience type to hold information about individual project files
type project_rec
    type(string) :: projname             ! project file name
    integer      :: micind     = 0       ! index of micrograph in project
    integer      :: nptcls     = 0       ! # of particles
    integer      :: nptcls_sel = 0       ! # of particles (state=1)
    logical      :: included   = .false. ! whether the record has been imported
end type project_rec

! Convenience type to hold information about processes
type process_rec
    type(string) :: str_id               ! unique string ID
    type(string) :: folder               ! location
    type(string) :: projfile             ! project filename
    type(string) :: volume               ! volume filename
    type(string) :: alnvolume            ! aligned volume filename
    logical      :: submitted = .false.  ! process has been submitted (running)
    logical      :: completed = .false.  ! has completed
    logical      :: included  = .false.  ! whether the record has been post-processed/analyzed
end type process_rec

! Convenience type to keep track of converged chunks
type chunk_rec
    integer      :: id        = 0        ! unique ID
    type(string) :: projfile             ! project filename
    logical      :: busy      = .false.  ! true after submission and until completion is detected
    logical      :: processed = .false.  ! chunk: has converged; set: has been clustered/selected/matched
    logical      :: included  = .false.  ! whether the set has been imported into the pool
end type chunk_rec

type poly_rec_list
    private
    type(linked_list) :: rec_list
contains
    ! constructors/lifecycle
    procedure          :: push
    procedure          :: copy_from
    procedure          :: replace
    procedure          :: kill
    ! overloaded assignment and append
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure, private :: append
    generic            :: operator(//)  => append
    ! checkers
    procedure          :: size
    procedure          :: is_empty
    procedure, private :: is_project
    ! ! setter
    ! procedure          :: push
    ! procedure          :: set
    ! ! getters
    ! procedure          :: get
end type poly_rec_list

contains

    subroutine push( self, rec )
        class(poly_rec_list), intent(inout) :: self
        class(*),             intent(in)    :: rec
        call self%rec_list%push_back(rec)
    end subroutine push

    subroutine copy_from( self, self_old )
        class(poly_rec_list), intent(inout) :: self
        class(poly_rec_list), intent(in)    :: self_old
        type(list_iterator)   :: self_iter, self_old_iter
        class(*), allocatable :: rec
        integer :: n, n_old
        if( self_old%rec_list%size() == 0 )then
            THROW_HARD('Nothing to copy_from, self_old%rec_list is empty')
        endif
        n = self%rec_list%size()
        n_old = self_old%rec_list%size()
        if( n_old <= n )then
            self_iter = self%rec_list%begin()
            self_old_iter = self_old%rec_list%begin()
            do while (self_old_iter%has_value())
                call self_old_iter%get(rec)
                call self_iter%set(rec)
                call self_iter%next()
                call self_old_iter%next()
            end do
        else
            call self_old_iter%advance(n)
            do while (self_old_iter%has_value())
                call self_old_iter%get(rec)
                call self%rec_list%push_back(rec)
                call self_old_iter%next()
            end do
        endif
    end subroutine copy_from

    subroutine replace( self, self_old )
        class(poly_rec_list), intent(inout) :: self
        class(poly_rec_list), intent(inout) :: self_old
        if( self_old%rec_list%size() == 0 )then
            THROW_HARD('Nothing to relace with, self_old%rec_list is empty')
        endif
        call self%rec_list%replace_with(self_old%rec_list)
    end subroutine replace

    subroutine kill( self )
        class(poly_rec_list), intent(inout) :: self
        call self%rec_list%kill
    end subroutine kill

    ! overloaded assignment and append

    subroutine assign(lhs, rhs)
        class(poly_rec_list), intent(inout) :: lhs
        class(poly_rec_list), intent(in)    :: rhs
        lhs%rec_list = rhs%rec_list
    end subroutine assign

    function append(lhs, rhs) result(res)
        class(poly_rec_list), intent(in) :: lhs
        class(poly_rec_list), intent(in) :: rhs
        type(poly_rec_list) :: res
        integer             :: n_lhs, n_rhs
        n_lhs = lhs%rec_list%size() 
        n_rhs = rhs%rec_list%size()
        if( n_lhs > 0 .and. n_rhs > 0  )then
            res%rec_list = lhs%rec_list // rhs%rec_list
        else if( n_lhs > 0 .and. n_rhs == 0  )then
            res%rec_list = lhs%rec_list
        else if( n_lhs == 0 .and. n_rhs > 0 )then
            res%rec_list = rhs%rec_list
        endif
    end function append

    ! checkers

    pure elemental integer function size(self) result(k)
        class(poly_rec_list), intent(in) :: self
        k = self%rec_list%size()
    end function size

    pure logical function is_empty(self) result(tf)
        class(poly_rec_list), intent(in) :: self
        tf = (self%rec_list%size() == 0)
    end function is_empty

    logical function is_project(self, ind, projname)
        class(poly_rec_list), intent(in) :: self
        integer,             intent(in) :: ind
        class(string),       intent(in) :: projname
        class(*), allocatable :: any
        integer :: n
        n = self%rec_list%size()
        if( n == 0 )then
            is_project = .false.
        else if( ind < 1 .or. ind > n )then
            is_project = .false.
        else
            call self%rec_list%at(ind, any)
            select type (any)
                type is (project_rec)
                    is_project = any%projname .eq. projname
                type is (process_rec)
                    is_project = .false.
                type is (chunk_rec)
                    is_project = .false.
                class default
                    THROW_HARD('unsupported type')
            end select
        end if
    end function is_project

    ! setters

    subroutine set( self, ind, rec )
        class(poly_rec_list), intent(inout) :: self
        integer,              intent(in)    :: ind
        class(*),             intent(in)    :: rec
        type(list_iterator) :: self_iter
        integer :: i, n
        n = self%rec_list%size() 
        if( n == 0 .and. ind == 1 )then
            call self%push(rec)
        else if( ind == n + 1 )then
            call self%push(rec)            
        else if( ind > 0 .and. ind <= n )then          
            i = 0
            self_iter = self%rec_list%begin()
            do while( self_iter%has_value() )
                i = i + 1
                if( i == ind )then
                    call self_iter%set(rec)
                    exit
                endif
                call self_iter%next()
            enddo
        else
            THROW_HARD('index ind out of bounds')
        end if
    end subroutine set

    ! ! getters

    ! type(string) function get_projname( self, ind )
    !     class(poly_rec_list), intent(in) :: self
    !     integer,             intent(in) :: ind
    !     if (.not. allocated(self%list)) then
    !         get_projname = ''
    !     else if (ind < 1 .or. ind > size(self%list)) then
    !         get_projname = ''
    !     else
    !         get_projname = self%list(ind)%projname
    !     end if
    ! end function get_projname

    ! pure integer function get_micind( self, ind )
    !     class(poly_rec_list), intent(in) :: self
    !     integer,             intent(in) :: ind
    !     if (.not. allocated(self%list)) then
    !         get_micind = 0
    !     else if (ind < 1 .or. ind > size(self%list)) then
    !         get_micind = 0
    !     else
    !         get_micind = self%list(ind)%micind
    !     end if
    ! end function get_micind

    ! pure integer function get_nptcls( self, ind )
    !     class(poly_rec_list), intent(in) :: self
    !     integer,             intent(in) :: ind
    !     if (.not. allocated(self%list)) then
    !         get_nptcls = 0
    !     else if (ind < 1 .or. ind > size(self%list)) then
    !         get_nptcls = 0
    !     else
    !         get_nptcls = self%list(ind)%nptcls
    !     end if
    ! end function get_nptcls

    ! pure integer function get_nptcls_sel( self, ind )
    !     class(poly_rec_list), intent(in) :: self
    !     integer,             intent(in) :: ind
    !     if (.not. allocated(self%list)) then
    !         get_nptcls_sel = 0
    !     else if (ind < 1 .or. ind > size(self%list)) then
    !         get_nptcls_sel = 0
    !     else
    !         get_nptcls_sel = self%list(ind)%nptcls_sel
    !     end if
    ! end function get_nptcls_sel

    ! pure function get_projname_arr( self ) result( pnames )
    !     class(poly_rec_list), intent(in) :: self
    !     type(string), allocatable :: pnames(:)
    !     integer :: n
    !     n = 0
    !     if( allocated(self%list) ) n = size(self%list)
    !     if( n > 0 )then
    !         pnames = self%list(:)%projname
    !     else    
    !         allocate(pnames(0))
    !     endif
    ! end function get_projname_arr

    ! pure function get_micind_arr( self ) result( micinds )
    !     class(poly_rec_list), intent(in) :: self
    !     integer, allocatable :: micinds(:)
    !     integer :: n
    !     n = 0
    !     if( allocated(self%list) ) n = size(self%list)
    !     if( n > 0 )then
    !         micinds = self%list(:)%micind
    !     else    
    !         allocate(micinds(0))
    !     endif
    ! end function get_micind_arr

    ! pure function get_nptcls_arr( self ) result( nptcls )
    !     class(poly_rec_list), intent(in) :: self
    !     integer, allocatable :: nptcls(:)
    !     integer :: n
    !     n = 0
    !     if( allocated(self%list) ) n = size(self%list)
    !     if( n > 0 )then
    !         nptcls = self%list(:)%nptcls
    !     else    
    !         allocate(nptcls(0))
    !     endif
    ! end function get_nptcls_arr

    ! pure function get_nptcls_sel_arr( self ) result( nptcls_sel )
    !     class(poly_rec_list), intent(in) :: self
    !     integer, allocatable :: nptcls_sel(:)
    !     integer :: n
    !     n = 0
    !     if( allocated(self%list) ) n = size(self%list)
    !     if( n > 0 )then
    !         nptcls_sel = self%list(:)%nptcls_sel
    !     else    
    !         allocate(nptcls_sel(0))
    !     endif
    ! end function get_nptcls_sel_arr

    ! integer function get_nptcls_tot(self)
    !     class(poly_rec_list), intent(in) :: self
    !     get_nptcls_tot = 0
    !     if( allocated(self%list) )then
    !         if( size(self%list) > 0 ) get_nptcls_tot = sum(self%list(:)%nptcls)
    !     end if
    ! end function get_nptcls_tot

    ! integer function get_nptcls_sel_tot(self)
    !     class(poly_rec_list), intent(in) :: self
    !     get_nptcls_sel_tot = 0
    !      if( allocated(self%list) )then
    !         if( size(self%list) > 0 ) get_nptcls_sel_tot = sum(self%list(:)%nptcls_sel)            
    !     end if
    ! end function get_nptcls_sel_tot

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

end module simple_poly_rec_list
