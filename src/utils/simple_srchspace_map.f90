module simple_srchspace_map
use simple_error,         only: simple_exception
use simple_srch_sort_loc, only: mask2inds
implicit none

public :: srchspace_map
private
#include "simple_local_flags.inc"

type :: srchspace_map
    private
    integer              :: nspace = 0
    integer              :: nsub   = 0
    integer, allocatable :: full2sub_map(:)  ! size(nspace)
    integer, allocatable :: sub2full_map(:)  ! size(nsub)
contains
    procedure :: new
    procedure :: get_inds_in_full
    procedure :: get_sub2full_map
    procedure :: get_full2sub_map
    procedure :: full2sub
    procedure :: sub2full
    procedure :: get_nspace
    procedure :: get_nsub
    procedure :: kill
end type srchspace_map

contains

    subroutine new(self, nspace, nspace_sub, distmat)
        class(srchspace_map), intent(inout) :: self
        integer,              intent(in)    :: nspace, nspace_sub
        real,                 intent(in)    :: distmat(nspace_sub, nspace)
        call self%kill()
        if (nspace_sub > nspace) THROW_HARD('subspace size cannot be larger than full space size')
        self%nspace = nspace
        self%nsub   = nspace_sub
        allocate(self%full2sub_map(nspace), self%sub2full_map(nspace_sub))
        ! For each full-space point (column), pick closest sub-space point (row)
        self%full2sub_map = minloc(distmat, dim=1)
        ! For each sub-space point (row), pick closest full-space point (column)
        self%sub2full_map = minloc(distmat, dim=2)
    end subroutine new

    pure function get_inds_in_full(self, sub_idx) result(inds)
        class(srchspace_map), intent(in) :: self
        integer,              intent(in) :: sub_idx
        integer, allocatable :: inds(:)
        logical :: mask(self%nspace)
        ! default empty
        allocate(inds(0))
        if (.not. allocated(self%full2sub_map)) return
        if (sub_idx < 1 .or. sub_idx > self%nsub) return
        mask = (self%full2sub_map == sub_idx)
        inds = mask2inds(mask)
    end function get_inds_in_full

    pure function get_sub2full_map(self) result(arr)
        class(srchspace_map), intent(in) :: self
        integer, allocatable :: arr(:)
        if (.not. allocated(self%sub2full_map)) then
            allocate(arr(0))
        else
            allocate(arr(size(self%sub2full_map)), source=self%sub2full_map)
        end if
    end function get_sub2full_map

    pure function get_full2sub_map(self) result(arr)
        class(srchspace_map), intent(in) :: self
        integer, allocatable :: arr(:)
        if (.not. allocated(self%full2sub_map)) then
            allocate(arr(0))
        else
            allocate(arr(size(self%full2sub_map)), source=self%full2sub_map)
        end if
    end function get_full2sub_map

    pure integer function full2sub(self, full_idx)
        class(srchspace_map), intent(in) :: self
        integer,              intent(in) :: full_idx
        full2sub = 0
        if (.not. allocated(self%full2sub_map)) return
        if (full_idx < 1 .or. full_idx > size(self%full2sub_map)) return
        full2sub = self%full2sub_map(full_idx)
    end function full2sub

    pure integer function sub2full(self, sub_idx)
        class(srchspace_map), intent(in) :: self
        integer,              intent(in) :: sub_idx
        sub2full = 0
        if (.not. allocated(self%sub2full_map)) return
        if (sub_idx < 1 .or. sub_idx > size(self%sub2full_map)) return
        sub2full = self%sub2full_map(sub_idx)
    end function sub2full

    pure integer function get_nspace(self)
        class(srchspace_map), intent(in) :: self
        get_nspace = self%nspace
    end function get_nspace

    pure integer function get_nsub(self)
        class(srchspace_map), intent(in) :: self
        get_nsub = self%nsub
    end function get_nsub

    subroutine kill(self)
        class(srchspace_map), intent(inout) :: self
        if (allocated(self%full2sub_map)) deallocate(self%full2sub_map)
        if (allocated(self%sub2full_map)) deallocate(self%sub2full_map)
        self%nspace = 0
        self%nsub   = 0
    end subroutine kill

end module simple_srchspace_map
