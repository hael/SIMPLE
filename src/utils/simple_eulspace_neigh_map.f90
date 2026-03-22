module simple_eulspace_neigh_map
use simple_core_module_api
use simple_srchspace_map, only: srchspace_map
implicit none

public :: eulspace_neigh_map
private
#include "simple_local_flags.inc"

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
        if (.not. allocated(self%full2sub_map)) THROW_HARD('get_neighbors_mask: map not initialized')
        if (sub_idx < 1 .or. sub_idx > self%nsub) THROW_HARD('get_neighbors_mask: sub_idx out of range')
        mask = (self%full2sub_map == sub_idx)
    end subroutine get_neighbors_mask

    subroutine get_neighbors_list(self, sub_idx, proj_list)
        class(eulspace_neigh_map), intent(in)  :: self
        integer,                   intent(in)  :: sub_idx
        integer, allocatable,      intent(out) :: proj_list(:)
        integer :: nsel, i

        if (.not. allocated(self%full2sub_map)) THROW_HARD('get_neighbors_list: map not initialized')
        if (sub_idx < 1 .or. sub_idx > self%nsub) THROW_HARD('get_neighbors_list: sub_idx out of range')

        nsel = count(self%full2sub_map == sub_idx)
        allocate(proj_list(nsel))
        if (nsel > 0) proj_list = pack([(i, i=1,self%nspace)], self%full2sub_map == sub_idx)
    end subroutine get_neighbors_list

    subroutine get_neighbors_mask_pooled(self, sub_idxs, mask)
        class(eulspace_neigh_map), intent(in)  :: self
        integer,                   intent(in)  :: sub_idxs(:)
        logical,                   intent(out) :: mask(self%nspace)
        integer :: i
        if (.not. allocated(self%full2sub_map)) THROW_HARD('get_neighbors_mask_pooled: map not initialized')
        if (size(sub_idxs) < 1) THROW_HARD('get_neighbors_mask_pooled: sub_idxs cannot be empty')
        mask = .false.
        do i = 1, size(sub_idxs)
            if (sub_idxs(i) < 1 .or. sub_idxs(i) > self%nsub) THROW_HARD('get_neighbors_mask_pooled: sub_idx out of range')
            mask = mask .or. (self%full2sub_map == sub_idxs(i))
        enddo
    end subroutine get_neighbors_mask_pooled

    function get_full2sub_map(self) result(labels)
        class(eulspace_neigh_map), intent(in) :: self
        integer, allocatable :: labels(:)
        if (.not. allocated(self%full2sub_map)) THROW_HARD('get_full2sub_map: map not initialized')
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

end module simple_eulspace_neigh_map
