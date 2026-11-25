!@descr: a class-average selection session on a project's cls2D segment: visibility, display order, selection and its write
! Reads cls2D pop, res and state; a class is shown unless its pop is present and 0, and starts selected when
! its state is above 0 (or absent). save writes the states as prg=selection oritype=cls2D does
! (report_state2cls2D), the four changed segments to the project file and a 0/1 selection.txt.
module simple_cls2D_selection
use simple_defs,       only: dp, logfhandle
use simple_error,      only: simple_exception
use simple_string,     only: string
use simple_sp_project, only: sp_project
implicit none

public :: cls2D_selection, SEL_SORT_RES, SEL_SORT_POP, SEL_SORT_CLASS
private
#include "simple_local_flags.inc"

integer,          parameter :: SEL_SORT_RES   = 1   ! resolution ascending (the start)
integer,          parameter :: SEL_SORT_POP   = 2   ! population descending
integer,          parameter :: SEL_SORT_CLASS = 3   ! class index ascending
character(len=*), parameter :: SELECTION_FNAME = 'selection.txt'

type :: cls2D_selection
    private
    class(sp_project), pointer :: spproj => null()
    type(string)               :: projfile
    real(dp), allocatable      :: pop(:), res(:)
    logical,  allocatable      :: l_pop(:), l_res(:), visible(:), selected(:)
    integer,  allocatable      :: order(:)        ! the visible classes in display order
    integer :: ncls     = 0
    integer :: nvis     = 0
    integer :: anchor   = 0                       ! the last clicked class (0: none)
    integer :: sort_key = SEL_SORT_RES
    logical :: append   = .false.
    logical :: dirty    = .false.
    logical :: l_exists = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: get_ncls
    procedure :: get_nvisible
    procedure :: get_class
    procedure :: get_sort
    procedure :: get_visible_mask
    procedure :: has_field
    procedure :: has_pop
    procedure :: has_res
    procedure :: get_pop
    procedure :: get_res
    procedure :: is_selected
    procedure :: is_dirty
    procedure :: counts
    procedure :: sort
    procedure :: click
    procedure :: select_all
    procedure :: deselect_all
    procedure :: save
end type cls2D_selection

contains

    ! the session on spproj (kept by pointer), written to projfile; append as prg=selection's
    subroutine new( self, spproj, projfile, append )
        class(cls2D_selection),    intent(inout) :: self
        class(sp_project), target, intent(inout) :: spproj
        class(string),             intent(in)    :: projfile
        logical,                   intent(in)    :: append
        integer :: icls
        call self%kill
        self%spproj   => spproj
        self%projfile = projfile
        self%append   = append
        self%ncls     = spproj%os_cls2D%get_noris()
        if( self%ncls == 0 ) THROW_HARD('empty cls2D segment; cls2D_selection%new')
        allocate(self%pop(self%ncls), self%res(self%ncls), self%l_pop(self%ncls), self%l_res(self%ncls),&
            &self%visible(self%ncls), self%selected(self%ncls))
        do icls = 1,self%ncls
            self%l_pop(icls) = spproj%os_cls2D%isthere(icls, 'pop')
            self%l_res(icls) = spproj%os_cls2D%isthere(icls, 'res')
            self%pop(icls)   = real(spproj%os_cls2D%get(icls, 'pop'), dp)
            self%res(icls)   = real(spproj%os_cls2D%get(icls, 'res'), dp)
            self%visible(icls) = .true.
            if( self%l_pop(icls) ) self%visible(icls) = self%pop(icls) /= 0._dp
            self%selected(icls) = self%visible(icls)
            if( self%visible(icls) )then
                if( spproj%os_cls2D%isthere(icls, 'state') ) self%selected(icls) = spproj%os_cls2D%get_state(icls) > 0
            endif
        end do
        self%nvis = count(self%visible)
        allocate(self%order(self%nvis))
        self%order = pack([(icls, icls = 1,self%ncls)], self%visible)
        self%anchor = 0
        self%dirty  = .false.
        call self%sort(SEL_SORT_RES)
        self%l_exists = .true.
    end subroutine new

    subroutine kill( self )
        class(cls2D_selection), intent(inout) :: self
        if( allocated(self%pop)      ) deallocate(self%pop)
        if( allocated(self%res)      ) deallocate(self%res)
        if( allocated(self%l_pop)    ) deallocate(self%l_pop)
        if( allocated(self%l_res)    ) deallocate(self%l_res)
        if( allocated(self%visible)  ) deallocate(self%visible)
        if( allocated(self%selected) ) deallocate(self%selected)
        if( allocated(self%order)    ) deallocate(self%order)
        call self%projfile%kill
        self%spproj   => null()
        self%ncls     = 0
        self%nvis     = 0
        self%anchor   = 0
        self%l_exists = .false.
    end subroutine kill

    integer function get_ncls( self )
        class(cls2D_selection), intent(in) :: self
        get_ncls = self%ncls
    end function get_ncls

    integer function get_nvisible( self )
        class(cls2D_selection), intent(in) :: self
        get_nvisible = self%nvis
    end function get_nvisible

    ! the class shown at display position pos (1 .. get_nvisible)
    integer function get_class( self, pos )
        class(cls2D_selection), intent(in) :: self
        integer,                intent(in) :: pos
        get_class = self%order(pos)
    end function get_class

    integer function get_sort( self )
        class(cls2D_selection), intent(in) :: self
        get_sort = self%sort_key
    end function get_sort

    function get_visible_mask( self ) result( mask )
        class(cls2D_selection), intent(in) :: self
        logical :: mask(self%ncls)
        mask = self%visible
    end function get_visible_mask

    ! whether any class carries the field ('pop' or 'res')
    logical function has_field( self, key )
        class(cls2D_selection), intent(in) :: self
        character(len=*),       intent(in) :: key
        select case(key)
            case('pop')
                has_field = any(self%l_pop)
            case('res')
                has_field = any(self%l_res)
            case DEFAULT
                has_field = .false.
        end select
    end function has_field

    logical function has_pop( self, icls )
        class(cls2D_selection), intent(in) :: self
        integer,                intent(in) :: icls
        has_pop = self%l_pop(icls)
    end function has_pop

    logical function has_res( self, icls )
        class(cls2D_selection), intent(in) :: self
        integer,                intent(in) :: icls
        has_res = self%l_res(icls)
    end function has_res

    real(dp) function get_pop( self, icls )
        class(cls2D_selection), intent(in) :: self
        integer,                intent(in) :: icls
        get_pop = self%pop(icls)
    end function get_pop

    real(dp) function get_res( self, icls )
        class(cls2D_selection), intent(in) :: self
        integer,                intent(in) :: icls
        get_res = self%res(icls)
    end function get_res

    logical function is_selected( self, icls )
        class(cls2D_selection), intent(in) :: self
        integer,                intent(in) :: icls
        is_selected = self%selected(icls)
    end function is_selected

    logical function is_dirty( self )
        class(cls2D_selection), intent(in) :: self
        is_dirty = self%dirty
    end function is_dirty

    ! selected and visible classes; particles (pop) of the selected classes and of all classes
    subroutine counts( self, nsel, nvis, psel, ptot )
        class(cls2D_selection), intent(in)  :: self
        integer,                intent(out) :: nsel, nvis
        integer(kind=8),        intent(out) :: psel, ptot
        integer :: icls
        nsel = count(self%selected)
        nvis = self%nvis
        psel = 0_8
        ptot = 0_8
        do icls = 1,self%ncls
            if( .not. self%l_pop(icls) ) cycle
            ptot = ptot + nint(self%pop(icls), kind=8)
            if( self%selected(icls) ) psel = psel + nint(self%pop(icls), kind=8)
        end do
    end subroutine counts

    ! re-sorts the display order, stable on the current order; a class without the key goes last
    subroutine sort( self, key )
        use simple_srch_sort_loc, only: scores2order_stable
        class(cls2D_selection), intent(inout) :: self
        integer,                intent(in)    :: key
        real(dp) :: scores(self%nvis)
        integer  :: perm(self%nvis), pos, icls
        select case(key)
            case(SEL_SORT_CLASS)
                self%order = pack([(icls, icls = 1,self%ncls)], self%visible)
            case(SEL_SORT_RES, SEL_SORT_POP)
                do pos = 1,self%nvis
                    icls = self%order(pos)
                    scores(pos) = -huge(1._dp)
                    if( key == SEL_SORT_RES )then
                        if( self%l_res(icls) ) scores(pos) = -self%res(icls)
                    else
                        if( self%l_pop(icls) ) scores(pos) = self%pop(icls)
                    endif
                end do
                perm = scores2order_stable(scores)
                self%order = self%order(perm)
            case DEFAULT
                THROW_HARD('unknown sort key; cls2D_selection%sort')
        end select
        self%sort_key = key
    end subroutine sort

    ! a click toggles icls; with shift and an earlier click, the range between them in display order
    ! is selected; every click becomes the next anchor
    subroutine click( self, icls, shift )
        class(cls2D_selection), intent(inout) :: self
        integer,                intent(in)    :: icls
        logical,                intent(in)    :: shift
        integer :: p1, p2
        if( .not. self%visible(icls) ) return
        p1 = 0
        if( shift .and. self%anchor > 0 ) p1 = findloc(self%order, self%anchor, dim=1)
        if( p1 > 0 )then
            p2 = findloc(self%order, icls, dim=1)
            self%selected(self%order(min(p1,p2):max(p1,p2))) = .true.
        else
            self%selected(icls) = .not. self%selected(icls)
        endif
        self%anchor = icls
        self%dirty  = .true.
    end subroutine click

    subroutine select_all( self )
        class(cls2D_selection), intent(inout) :: self
        self%selected = self%visible
        self%anchor   = 0
        if( self%nvis > 0 ) self%anchor = self%order(self%nvis)
        self%dirty = .true.
    end subroutine select_all

    subroutine deselect_all( self )
        class(cls2D_selection), intent(inout) :: self
        self%selected = .false.
        self%anchor   = 0
        self%dirty    = .true.
    end subroutine deselect_all

    ! the selection into the project in memory and on disk (cls2D, cls3D, ptcl2D, ptcl3D), and selection.txt
    subroutine save( self )
        use simple_fileio, only: fopen, fclose
        class(cls2D_selection), intent(inout) :: self
        integer, allocatable :: states(:)
        integer :: funit, icls
        if( .not. self%l_exists ) THROW_HARD('no selection session; cls2D_selection%save')
        allocate(states(self%ncls), source=0)
        where( self%selected ) states = 1
        call self%spproj%report_state2cls2D(states, self%append)
        call self%spproj%write_segment_inside('cls2D', self%projfile)
        if( self%spproj%os_cls3D%get_noris()  > 0 ) call self%spproj%write_segment_inside('cls3D',  self%projfile)
        if( self%spproj%os_ptcl2D%get_noris() > 0 ) call self%spproj%write_segment_inside('ptcl2D', self%projfile)
        if( self%spproj%os_ptcl3D%get_noris() > 0 ) call self%spproj%write_segment_inside('ptcl3D', self%projfile)
        call fopen(funit, file=string(SELECTION_FNAME), status='replace', action='write')
        do icls = 1,self%ncls
            write(funit,'(i0)') states(icls)
        end do
        call fclose(funit)
        self%dirty = .false.
        write(logfhandle,'(a,i0,a,i0,a)') '>>> SELECTION SAVED: ', count(states == 1), ' classes kept, ', &
            &count(states == 0), ' rejected, to '//self%projfile%to_char()//' and '//SELECTION_FNAME
    end subroutine save

end module simple_cls2D_selection
