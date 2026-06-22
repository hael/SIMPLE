module simple_eulspace_neigh_map
use simple_core_module_api
use simple_srchspace_map, only: srchspace_map
implicit none

public :: eulspace_neigh_map
public :: eulspace_neigh_groups
private
#include "simple_local_flags.inc"

type :: eulspace_neigh_groups
    private
    integer :: nsub   = 0
    integer :: nitems = 0
    integer, allocatable :: rowptr(:)
    integer, allocatable :: item_inds(:)
    integer, allocatable :: pops(:)
contains
    procedure :: get_group  => groups_get_group
    procedure :: get_pops   => groups_get_pops
    procedure :: get_pop    => groups_get_pop
    procedure :: get_nsub   => groups_get_nsub
    procedure :: get_nitems => groups_get_nitems
    procedure :: has_data   => groups_has_data
    procedure :: kill       => groups_kill
end type eulspace_neigh_groups

type :: eulspace_neigh_map
    private
    integer :: nspace = 0
    integer :: nsub   = 0
    integer, allocatable :: full2sub_map(:) ! size(nspace), values in [1..nsub]
contains
    procedure, private :: new_from_labels
    procedure, private :: new_from_spaces
    generic            :: new => new_from_labels, new_from_spaces
    procedure          :: get_neighbors_mask
    procedure          :: get_neighbors_list
    procedure          :: get_neighbors_mask_pooled
    procedure          :: group_by_full_inds
    procedure          :: get_full2sub_map
    procedure          :: get_nspace
    procedure          :: get_nsub
    procedure          :: kill
end type eulspace_neigh_map

contains

    subroutine new_from_labels(self, full2sub_map, nspace_sub)
        class(eulspace_neigh_map), intent(inout) :: self
        integer,                   intent(in)    :: full2sub_map(:)
        integer,                   intent(in)    :: nspace_sub
        call self%kill
        if (nspace_sub < 1) THROW_HARD('new_from_labels: nspace_sub must be >= 1')
        if (size(full2sub_map) < 1) THROW_HARD('new_from_labels: full2sub_map cannot be empty')
        if (minval(full2sub_map) < 1 .or. maxval(full2sub_map) > nspace_sub) then
            THROW_HARD('new_from_labels: full2sub_map has values outside [1,nspace_sub]')
        endif
        self%nspace = size(full2sub_map)
        self%nsub   = nspace_sub
        allocate(self%full2sub_map(self%nspace), source=full2sub_map)
    end subroutine new_from_labels

    subroutine new_from_spaces(self, eulspace, eulspace_sub, pgrpsym, normalize_dist)
        class(eulspace_neigh_map), intent(inout) :: self
        class(oris),               intent(in)    :: eulspace, eulspace_sub
        class(sym),                intent(inout) :: pgrpsym
        logical, optional,         intent(in)    :: normalize_dist
        type(srchspace_map) :: mapper
        type(ori) :: o, o_sub, osym
        integer, allocatable :: labels(:)
        integer :: i, j, nspace, nspace_sub
        real, allocatable :: distmat(:,:)
        real :: dtmp, inplrotdist
        logical :: l_norm
        l_norm = .true.
        if (present(normalize_dist)) l_norm = normalize_dist
        nspace     = eulspace%get_noris()
        nspace_sub = eulspace_sub%get_noris()
        if (nspace_sub < 1 .or. nspace < 1) THROW_HARD('new_from_spaces: empty eulspace/eulspace_sub')
        allocate(distmat(nspace_sub, nspace))
        !$omp parallel do default(shared) private(i,j,o,o_sub,osym,dtmp,inplrotdist) proc_bind(close) schedule(static)
        do i = 1, nspace_sub
            call eulspace_sub%get_ori(i, o_sub)
            do j = 1, nspace
                call eulspace%get_ori(j, o)
                call pgrpsym%sym_dists(o_sub, o, osym, dtmp, inplrotdist)
                distmat(i,j) = dtmp
            enddo
        enddo
        !$omp end parallel do
        if (l_norm) call normalize_minmax(distmat)
        call mapper%new(nspace, nspace_sub, distmat)
        labels = mapper%get_full2sub_map()
        call self%new(labels, nspace_sub)
        if (allocated(labels)) deallocate(labels)
        if (allocated(distmat)) deallocate(distmat)
        call mapper%kill
    end subroutine new_from_spaces

    subroutine get_neighbors_mask(self, sub_idx, mask)
        class(eulspace_neigh_map), intent(in)  :: self
        integer,                   intent(in)  :: sub_idx
        logical,                   intent(out) :: mask(self%nspace)
        mask = (self%full2sub_map == sub_idx)
    end subroutine get_neighbors_mask

    subroutine get_neighbors_list(self, sub_idx, proj_list)
        class(eulspace_neigh_map), intent(in)    :: self
        integer,                   intent(in)    :: sub_idx
        integer, allocatable,      intent(inout) :: proj_list(:)
        integer :: nsel, i
        if (allocated(proj_list)) deallocate(proj_list)
        nsel = count(self%full2sub_map == sub_idx)
        allocate(proj_list(nsel))
        if (nsel > 0) proj_list = pack([(i, i=1,self%nspace)], self%full2sub_map == sub_idx)
    end subroutine get_neighbors_list

    subroutine get_neighbors_mask_pooled(self, sub_idxs, mask)
        class(eulspace_neigh_map), intent(in)  :: self
        integer,                   intent(in)  :: sub_idxs(:)
        logical,                   intent(out) :: mask(self%nspace)
        integer :: i
        mask = .false.
        do i = 1, size(sub_idxs)
            mask = mask .or. (self%full2sub_map == sub_idxs(i))
        enddo
    end subroutine get_neighbors_mask_pooled

    subroutine group_by_full_inds(self, full_inds, groups, item_inds)
        class(eulspace_neigh_map),    intent(in)    :: self
        integer,                      intent(in)    :: full_inds(:)
        type(eulspace_neigh_groups),  intent(inout) :: groups
        integer, optional,            intent(in)    :: item_inds(:)
        integer, allocatable :: cursor(:)
        integer :: i, isub, nitems, pos
        call groups%kill
        if (.not. allocated(self%full2sub_map)) THROW_HARD('group_by_full_inds: empty eulspace_neigh_map')
        nitems = size(full_inds)
        if (present(item_inds)) then
            if (size(item_inds) /= nitems) THROW_HARD('group_by_full_inds: item_inds/full_inds size mismatch')
        endif
        groups%nsub   = self%nsub
        groups%nitems = nitems
        allocate(groups%pops(self%nsub), source=0)
        do i = 1, nitems
            if (full_inds(i) < 1 .or. full_inds(i) > self%nspace) then
                THROW_HARD('group_by_full_inds: full index outside eulspace_neigh_map range')
            endif
            isub = self%full2sub_map(full_inds(i))
            groups%pops(isub) = groups%pops(isub) + 1
        enddo
        allocate(groups%rowptr(self%nsub + 1), source=1)
        do i = 1, self%nsub
            groups%rowptr(i + 1) = groups%rowptr(i) + groups%pops(i)
        enddo
        allocate(groups%item_inds(nitems))
        allocate(cursor(self%nsub), source=groups%rowptr(1:self%nsub))
        do i = 1, nitems
            isub = self%full2sub_map(full_inds(i))
            pos  = cursor(isub)
            if (present(item_inds)) then
                groups%item_inds(pos) = item_inds(i)
            else
                groups%item_inds(pos) = i
            endif
            cursor(isub) = cursor(isub) + 1
        enddo
        if (allocated(cursor)) deallocate(cursor)
    end subroutine group_by_full_inds

    function get_full2sub_map(self) result(labels)
        class(eulspace_neigh_map), intent(in) :: self
        integer, allocatable :: labels(:)
        allocate(labels(self%nspace), source=self%full2sub_map)
    end function get_full2sub_map

    pure integer function get_nspace(self)
        class(eulspace_neigh_map), intent(in) :: self
        get_nspace = self%nspace
    end function get_nspace

    pure integer function get_nsub(self)
        class(eulspace_neigh_map), intent(in) :: self
        get_nsub = self%nsub
    end function get_nsub

    subroutine kill(self)
        class(eulspace_neigh_map), intent(inout) :: self
        if (allocated(self%full2sub_map)) deallocate(self%full2sub_map)
        self%nspace = 0
        self%nsub   = 0
    end subroutine kill

    function groups_get_group(self, sub_idx) result(items)
        class(eulspace_neigh_groups), intent(in) :: self
        integer,                       intent(in) :: sub_idx
        integer, allocatable :: items(:)
        integer :: nsel
        allocate(items(0))
        if (.not. allocated(self%rowptr)) return
        if (.not. allocated(self%item_inds)) return
        if (sub_idx < 1 .or. sub_idx > self%nsub) return
        nsel = self%rowptr(sub_idx + 1) - self%rowptr(sub_idx)
        deallocate(items)
        allocate(items(nsel))
        if (nsel > 0) items = self%item_inds(self%rowptr(sub_idx):self%rowptr(sub_idx + 1) - 1)
    end function groups_get_group

    function groups_get_pops(self) result(arr)
        class(eulspace_neigh_groups), intent(in) :: self
        integer, allocatable :: arr(:)
        if (.not. allocated(self%pops)) then
            allocate(arr(0))
        else
            allocate(arr(size(self%pops)), source=self%pops)
        endif
    end function groups_get_pops

    pure integer function groups_get_pop(self, sub_idx)
        class(eulspace_neigh_groups), intent(in) :: self
        integer,                       intent(in) :: sub_idx
        groups_get_pop = 0
        if (.not. allocated(self%pops)) return
        if (sub_idx < 1 .or. sub_idx > self%nsub) return
        groups_get_pop = self%pops(sub_idx)
    end function groups_get_pop

    pure integer function groups_get_nsub(self)
        class(eulspace_neigh_groups), intent(in) :: self
        groups_get_nsub = self%nsub
    end function groups_get_nsub

    pure integer function groups_get_nitems(self)
        class(eulspace_neigh_groups), intent(in) :: self
        groups_get_nitems = self%nitems
    end function groups_get_nitems

    pure logical function groups_has_data(self)
        class(eulspace_neigh_groups), intent(in) :: self
        groups_has_data = allocated(self%rowptr) .and. allocated(self%item_inds) .and. allocated(self%pops)
    end function groups_has_data

    subroutine groups_kill(self)
        class(eulspace_neigh_groups), intent(inout) :: self
        if (allocated(self%rowptr))    deallocate(self%rowptr)
        if (allocated(self%item_inds)) deallocate(self%item_inds)
        if (allocated(self%pops))      deallocate(self%pops)
        self%nsub   = 0
        self%nitems = 0
    end subroutine groups_kill

end module simple_eulspace_neigh_map
