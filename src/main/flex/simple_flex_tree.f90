!@descr: indexed volume hierarchy with project, fitted-model and decision records
module simple_flex_tree
use, intrinsic :: iso_fortran_env, only: int64
use simple_indexed_tree, only: indexed_tree
use simple_string,       only: string
use simple_error,        only: simple_exception
use simple_defs,         only: dp
implicit none
private
public :: flex_tree, flex_tree_node, flex_tree_model, flex_tree_decision
public :: TREE_UNASSESSED, TREE_SPLIT, TREE_MERGE, TREE_LEAF, TREE_JUNK
#include "simple_local_flags.inc"

integer, parameter :: TREE_UNASSESSED = 0, TREE_SPLIT = 1, TREE_MERGE = 2, TREE_LEAF = 3, TREE_JUNK = 4

! One named fit, e.g. initial/focused/paired. The embedding owns latents, eigenvalues and precisions.
! Its particle indices are rows of this fit's project; that project's pind maps to the source.
type :: flex_tree_model
    type(string) :: name, project_file, embedding_file
    type(string) :: mean_file, mask_file, noise_file
    type(string) :: mean_scale_file ! optional scale applied to mean_file for this fit
    type(string), allocatable :: basis_files(:)
    type(string) :: weights_file ! state-weight manifest, distinct from latent coefficients
    type(string) :: targets_file ! state coordinates, bandwidths and effective populations
    type(string) :: provenance_file
    type(string), allocatable :: diagnostics(:)
    integer :: ncomp = 0, box = 0
    real    :: smpd = 0. ! working lattice; zero box/smpd means not recorded
    logical :: unit_noise = .false. ! true records the unit-spectrum fallback
    integer(int64) :: layout_digest = 0_int64 ! zero means not recorded
end type flex_tree_model

! Reports carry split/site evidence; small state-tree merge statistics stay directly accessible.
type :: flex_tree_decision
    integer :: outcome = TREE_UNASSESSED
    type(string) :: reason
    type(string), allocatable :: reports(:)
    logical :: has_merge_stats = .false.
    real(dp) :: merge_d2 = 0._dp, merge_z = 0._dp
    real(dp), allocatable :: child_populations(:) ! get_children() order
end type flex_tree_decision

! Copyable payload: absent fields describe results that have not been produced.
! Explicit maps describe this node's result (including an aggregate); otherwise use the project registry.
type :: flex_tree_node
    type(string) :: name, project_file
    integer, allocatable :: project_states(:) ! absent/empty selects the whole node project
    type(string) :: map_file, even_file, odd_file, mask_file, fsc_file
    type(string), allocatable :: map_metadata_files(:) ! support/provenance sidecars
    type(string), allocatable :: site_files(:)
    type(flex_tree_model),    allocatable :: models(:)
    type(flex_tree_decision), allocatable :: decision
end type flex_tree_node

! Topology and payloads share stable node IDs. Projects own membership and canonical pind values.
! Callers create and retain referenced files; this container performs no project I/O or remapping.
type, extends(indexed_tree) :: flex_tree
    private
    type(flex_tree_node), allocatable :: nodes(:)
    type(string) :: source_project
  contains
    procedure :: new
    procedure :: kill
    procedure :: add_root
    procedure :: add_child
    procedure :: set_source_project
    procedure :: get_source_project
    procedure :: set_node
    procedure :: get_node
    procedure, private :: check_payload_node
    procedure, private :: validate_node
    procedure, private :: reserve_nodes
end type flex_tree

contains

    ! Keep the base constructor's interface; every payload and the source reference start empty.
    subroutine new( self, nroots )
        class(flex_tree), intent(inout) :: self
        integer,            intent(in)    :: nroots
        if( nroots < 0 ) THROW_HARD('flex_tree: negative root count; new')
        call self%kill
        call self%indexed_tree%new(nroots)
        call self%reserve_nodes(nroots)
    end subroutine new

    subroutine kill( self )
        class(flex_tree), intent(inout) :: self
        if( allocated(self%nodes) ) deallocate(self%nodes)
        call self%source_project%kill
        call self%indexed_tree%kill
    end subroutine kill

    ! Inherited join_roots dispatches here; a joined parent has no inferred project, model or decision.
    integer function add_root( self ) result( inode )
        class(flex_tree), intent(inout) :: self
        inode = self%indexed_tree%add_root()
        call self%reserve_nodes(inode)
    end function add_root

    ! A changed child list invalidates the parent's earlier split/merge/terminal decision.
    integer function add_child( self, iparent ) result( inode )
        class(flex_tree), intent(inout) :: self
        integer,            intent(in)    :: iparent
        inode = self%indexed_tree%add_child(iparent)
        call self%reserve_nodes(inode)
        if( allocated(self%nodes(iparent)%decision) ) deallocate(self%nodes(iparent)%decision)
    end function add_child

    ! Changing the source would reinterpret every stored pind; new/kill starts a new namespace.
    subroutine set_source_project( self, project_file )
        class(flex_tree), intent(inout) :: self
        class(string),      intent(in)    :: project_file
        if( project_file%is_blank() ) THROW_HARD('flex_tree: empty source project')
        if( .not. self%source_project%is_blank() )then
            if( self%source_project%to_char() /= project_file%to_char() ) THROW_HARD('flex_tree: source project already set')
        endif
        self%source_project = project_file%to_char()
    end subroutine set_source_project

    function get_source_project( self ) result( project_file )
        class(flex_tree), intent(in) :: self
        type(string) :: project_file
        project_file = self%source_project%to_char()
    end function get_source_project

    ! Replace the complete payload after validation. Intrinsic assignment copies all nested allocations.
    subroutine set_node( self, inode, data )
        class(flex_tree),     intent(inout) :: self
        integer,                intent(in)    :: inode
        type(flex_tree_node), intent(in)    :: data
        call self%check_payload_node(inode)
        call self%validate_node(inode, data)
        self%nodes(inode) = data
    end subroutine set_node

    ! Return a copy: editing it changes the tree only after a successful set_node.
    function get_node( self, inode ) result( data )
        class(flex_tree), intent(in) :: self
        integer,            intent(in) :: inode
        type(flex_tree_node) :: data
        call self%check_payload_node(inode)
        data = self%nodes(inode)
    end function get_node

    subroutine validate_node( self, inode, data )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        class(flex_tree),     intent(in) :: self
        integer,                intent(in) :: inode
        type(flex_tree_node), intent(in) :: data
        integer, allocatable :: children(:)
        integer :: i, j
        if( allocated(data%project_states) )then
            if( size(data%project_states) > 0 )then
                if( data%project_file%is_blank() ) THROW_HARD('flex_tree: state selection has no project')
            endif
            do i = 1, size(data%project_states)
                if( data%project_states(i) < 1 ) THROW_HARD('flex_tree: project states must be positive')
                if( any(data%project_states(:i-1) == data%project_states(i)) ) THROW_HARD('flex_tree: duplicate project state')
            end do
        endif
        if( data%even_file%is_blank() .neqv. data%odd_file%is_blank() ) THROW_HARD('flex_tree: supply both half maps')
        if( allocated(data%map_metadata_files) ) call check_files(data%map_metadata_files)
        if( allocated(data%site_files) ) call check_files(data%site_files)
        if( allocated(data%models) )then
            do i = 1, size(data%models)
                call check_model(data%models(i))
                do j = 1, i-1
                    if( data%models(i)%name%to_char() == data%models(j)%name%to_char() ) &
                        &THROW_HARD('flex_tree: duplicate model name')
                end do
            end do
        endif
        if( .not. allocated(data%decision) ) return
        associate( decision => data%decision )
            children = self%get_children(inode)
            select case(decision%outcome)
                case(TREE_UNASSESSED)
                case(TREE_SPLIT, TREE_MERGE)
                    if( size(children) < 2 ) THROW_HARD('flex_tree: split/merge decision needs at least two children')
                case(TREE_LEAF, TREE_JUNK)
                    if( size(children) /= 0 ) THROW_HARD('flex_tree: terminal decision on a node with children')
                case DEFAULT
                    THROW_HARD('flex_tree: invalid decision outcome')
            end select
            if( allocated(decision%reports) ) call check_files(decision%reports)
            if( decision%has_merge_stats )then
                if( decision%outcome /= TREE_MERGE ) THROW_HARD('flex_tree: merge statistics on a non-merge decision')
                if( .not. allocated(decision%child_populations) ) THROW_HARD('flex_tree: missing merge populations')
                if( size(decision%child_populations) /= size(children) ) THROW_HARD('flex_tree: merge population count mismatch')
                if( .not. ieee_is_finite(decision%merge_d2) ) THROW_HARD('flex_tree: nonfinite merge distance')
                if( .not. ieee_is_finite(decision%merge_z) ) THROW_HARD('flex_tree: nonfinite merge significance')
                if( any(.not. ieee_is_finite(decision%child_populations)) ) THROW_HARD('flex_tree: nonfinite merge populations')
                if( any(decision%child_populations < 0._dp) ) THROW_HARD('flex_tree: negative merge populations')
            else
                if( allocated(decision%child_populations) ) THROW_HARD('flex_tree: populations supplied without merge statistics')
            endif
        end associate
    end subroutine validate_node

    subroutine check_model( model )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        type(flex_tree_model), intent(in) :: model
        if( model%name%is_blank() ) THROW_HARD('flex_tree: unnamed model')
        if( model%project_file%is_blank() ) THROW_HARD('flex_tree: model has no fitting project')
        if( model%ncomp < 0 .or. model%box < 0 ) THROW_HARD('flex_tree: negative model dimensions')
        if( .not. ieee_is_finite(model%smpd) ) THROW_HARD('flex_tree: nonfinite model sampling')
        if( model%smpd < 0. ) THROW_HARD('flex_tree: negative model sampling')
        if( (model%box == 0) .neqv. (model%smpd == 0.) ) THROW_HARD('flex_tree: supply model box and sampling together')
        if( model%unit_noise )then
            if( .not. model%noise_file%is_blank() ) THROW_HARD('flex_tree: both unit noise and a noise file supplied')
        endif
        if( allocated(model%basis_files) )then
            call check_files(model%basis_files)
            if( size(model%basis_files) > 0 )then
                if( size(model%basis_files) /= model%ncomp ) THROW_HARD('flex_tree: basis count differs from model rank')
                if( model%mean_file%is_blank() ) THROW_HARD('flex_tree: basis has no matching mean')
                if( model%box == 0 ) THROW_HARD('flex_tree: basis has no working lattice')
            endif
        endif
        if( allocated(model%diagnostics) ) call check_files(model%diagnostics)
    end subroutine check_model

    subroutine check_files( files )
        type(string), intent(in) :: files(:)
        integer :: i
        do i = 1, size(files)
            if( files(i)%is_blank() ) THROW_HARD('flex_tree: blank entry in file list')
        end do
    end subroutine check_files

    subroutine check_payload_node( self, inode )
        class(flex_tree), intent(in) :: self
        integer,            intent(in) :: inode
        if( inode < 1 .or. inode > self%get_nnodes() ) THROW_HARD('flex_tree: node index out of range')
    end subroutine check_payload_node

    subroutine reserve_nodes( self, needed )
        class(flex_tree), intent(inout) :: self
        integer,            intent(in)    :: needed
        type(flex_tree_node), allocatable :: tmp(:)
        integer :: capacity, newcap
        capacity = 0
        if( allocated(self%nodes) ) capacity = size(self%nodes)
        if( needed <= capacity ) return
        newcap = max(16, needed)
        if( capacity <= shiftr(huge(capacity), 1) ) newcap = max(newcap, 2 * capacity)
        allocate(tmp(newcap))
        if( capacity > 0 ) tmp(:capacity) = self%nodes
        call move_alloc(tmp, self%nodes)
    end subroutine reserve_nodes

end module simple_flex_tree
