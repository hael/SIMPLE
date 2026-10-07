!@descr: flex_pca state-reconstruction contract, gridding/PCG backends, delivery and service
module simple_flex_pca_states_backend
use simple_core_module_api, only: simple_exception
use simple_builder,         only: builder
use simple_image,           only: image
use simple_parameters,      only: parameters
use simple_flex_pca_rounds, only: flex_pca_rounds
implicit none

public :: flex_states_backend, flex_state_maps, flex_state_delivery_policy
public :: flex_rec_box, flex_rec_smpd
private

!> the views of one reconstructed state, as the backend finalizes them
type :: flex_state_maps
    type(image) :: combined, even, odd
    logical     :: l_halves = .false.   !< even/odd present (a halfset-split round)
  contains
    procedure :: kill => flex_state_maps_kill
end type flex_state_maps

!> what the common delivery applies to a backend's maps. These are the backends' explicit
!! differences (doc/refactoring_notes/completed/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18.md
!! 6.7): each backend declares its own values, nothing is normalised between them.
type :: flex_state_delivery_policy
    logical :: l_mask         = .true.    !< background removal + soft spherical mask, on the FSC copies and the delivered views
    logical :: l_project_fsc_fallback = .true.  !< the project FSC low-pass when the state eo-FSC is unmeasurable, and on single-set rounds
    character(len=8) :: tag = ''          !< log tag after 'FLEX STATE'
end type flex_state_delivery_policy

!> One state per weight column, per halfset when l_fuse. The service drives every backend the
!! same way: begin, then on a worker accumulate_local_or_write_part; on the distributed master
!! fold_parts, in shared memory accumulate_local_or_write_part; then finalize_maps per state,
!! then kill. The selection (particle rows and weight table) is held here for the backends.
type, abstract :: flex_states_backend
    integer, allocatable :: pinds(:)
    real,    allocatable :: state_weights(:,:)
    integer :: nstates = 0, box_rec = 0
    real    :: smpd_rec = 0.
    logical :: l_fuse = .false., l_floor_rho = .false.
  contains
    procedure(begin_iface),    deferred :: begin
    procedure(pass_iface),     deferred :: accumulate_local_or_write_part
    procedure(pass_iface),     deferred :: fold_parts
    procedure(finalize_iface), deferred :: finalize_maps
    procedure(policy_iface),   deferred :: delivery_policy
    procedure(kill_iface),     deferred :: kill
    procedure :: set_selection  => backend_set_selection
    procedure :: kill_selection => backend_kill_selection
end type flex_states_backend

abstract interface
    subroutine begin_iface( self, params, build, rounds, pinds, state_weights, nstates, l_fuse, l_floor_rho, box_rec, smpd_rec )
        import :: flex_states_backend, parameters, builder, flex_pca_rounds
        class(flex_states_backend), intent(inout) :: self
        class(parameters),          intent(inout) :: params
        class(builder),             intent(inout) :: build
        class(flex_pca_rounds),     intent(inout) :: rounds
        integer,                    intent(in)    :: pinds(:), nstates, box_rec
        real,                       intent(in)    :: state_weights(:,:), smpd_rec
        logical,                    intent(in)    :: l_fuse, l_floor_rho
    end subroutine begin_iface

    subroutine pass_iface( self, params, build, rounds )
        import :: flex_states_backend, parameters, builder, flex_pca_rounds
        class(flex_states_backend), intent(inout) :: self
        class(parameters),          intent(inout) :: params
        class(builder),             intent(inout) :: build
        class(flex_pca_rounds),     intent(inout) :: rounds
    end subroutine pass_iface

    subroutine finalize_iface( self, params, build, rounds, state, maps )
        import :: flex_states_backend, flex_state_maps, parameters, builder, flex_pca_rounds
        class(flex_states_backend), intent(inout) :: self
        class(parameters),          intent(inout) :: params
        class(builder),             intent(inout) :: build
        class(flex_pca_rounds),     intent(inout) :: rounds
        integer,                    intent(in)    :: state
        type(flex_state_maps),      intent(inout) :: maps
    end subroutine finalize_iface

    function policy_iface( self ) result( policy )
        import :: flex_states_backend, flex_state_delivery_policy
        class(flex_states_backend), intent(in) :: self
        type(flex_state_delivery_policy) :: policy
    end function policy_iface

    subroutine kill_iface( self )
        import :: flex_states_backend
        class(flex_states_backend), intent(inout) :: self
    end subroutine kill_iface
end interface

contains

    !> Box/sampling of the delivered state maps. Decoupled from box_crop because the embedding is a
    !! low-frequency object while the state maps are plain backprojections that carry signal beyond
    !! the covariance band. With box_rec==box_crop this is a no-op.
    pure integer function flex_rec_box( params ) result( box_rec )
        class(parameters), intent(in) :: params
        box_rec = params%box_crop
        if( params%box_rec >= 1 ) box_rec = params%box_rec
    end function flex_rec_box

    pure real function flex_rec_smpd( params ) result( smpd_rec )
        class(parameters), intent(in) :: params
        smpd_rec = params%smpd_crop
        if( params%box_rec >= 1 .and. params%smpd_rec > 0. ) smpd_rec = params%smpd_rec
    end function flex_rec_smpd

    !> The shared part of begin: the selection and the lattice.
    subroutine backend_set_selection( self, pinds, state_weights, nstates, l_fuse, l_floor_rho, box_rec, smpd_rec )
        class(flex_states_backend), intent(inout) :: self
        integer,                    intent(in)    :: pinds(:), nstates, box_rec
        real,                       intent(in)    :: state_weights(:,:), smpd_rec
        logical,                    intent(in)    :: l_fuse, l_floor_rho
        call self%kill_selection
        allocate(self%pinds(size(pinds)), source=pinds)
        allocate(self%state_weights(size(pinds),nstates), source=state_weights(:,1:nstates))
        self%nstates     = nstates
        self%l_fuse      = l_fuse
        self%l_floor_rho = l_floor_rho
        self%box_rec     = box_rec
        self%smpd_rec    = smpd_rec
    end subroutine backend_set_selection

    subroutine backend_kill_selection( self )
        class(flex_states_backend), intent(inout) :: self
        if( allocated(self%pinds) )         deallocate(self%pinds)
        if( allocated(self%state_weights) ) deallocate(self%state_weights)
        self%nstates = 0; self%box_rec = 0; self%smpd_rec = 0.
        self%l_fuse = .false.; self%l_floor_rho = .false.
    end subroutine backend_kill_selection

    subroutine flex_state_maps_kill( self )
        class(flex_state_maps), intent(inout) :: self
        call self%combined%kill
        call self%even%kill
        call self%odd%kill
        self%l_halves = .false.
    end subroutine flex_state_maps_kill

end module simple_flex_pca_states_backend
