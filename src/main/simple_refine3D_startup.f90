!@descr: continuation startup decision and effective first-iteration sampling
module simple_refine3D_startup
implicit none
private
public :: refine3D_startup, STARTUP_ORDINARY, STARTUP_REUSE_CHAIN, STARTUP_SEED_CHAIN

integer, parameter :: STARTUP_ORDINARY    = 0
integer, parameter :: STARTUP_REUSE_CHAIN = 1
integer, parameter :: STARTUP_SEED_CHAIN  = 2

type :: refine3D_startup
    private
    integer :: mode = STARTUP_ORDINARY
    real    :: requested_fraction = 1.
    logical :: requested_fractional = .false.
    logical :: requested_trailing = .false.
    logical :: full_update_active = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: read_cline
    procedure :: write_cline
    procedure :: preserves_history
    procedure :: begin_iteration
    procedure :: finish_iteration
end type refine3D_startup

contains

    subroutine new( self, mode )
        class(refine3D_startup), intent(inout) :: self
        integer,                intent(in)     :: mode
        call self%kill
        self%mode = mode
    end subroutine new

    subroutine kill( self )
        class(refine3D_startup), intent(inout) :: self
        self%mode = STARTUP_ORDINARY
        self%requested_fraction = 1.
        self%requested_fractional = .false.
        self%requested_trailing = .false.
        self%full_update_active = .false.
    end subroutine kill

    ! Decode the child-command handoff once; subsequent decisions use the typed state.
    subroutine read_cline( self, cline )
        use simple_cmdline, only: cmdline
        use simple_string,  only: string
        class(refine3D_startup), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(string) :: value
        call self%new(STARTUP_ORDINARY)
        if( .not. cline%defined('trail_seed') ) return
        value = cline%get_carg('trail_seed')
        select case(value%to_char())
            case('no')
                self%mode = STARTUP_REUSE_CHAIN
            case('yes')
                self%mode = STARTUP_SEED_CHAIN
        end select
        call value%kill
        call cline%delete('trail_seed')
    end subroutine read_cline

    subroutine write_cline( self, cline )
        use simple_cmdline, only: cmdline
        class(refine3D_startup), intent(in)    :: self
        class(cmdline),          intent(inout) :: cline
        select case(self%mode)
            case(STARTUP_REUSE_CHAIN)
                call cline%set('trail_seed', 'no')
            case(STARTUP_SEED_CHAIN)
                call cline%set('trail_seed', 'yes')
            case default
                call cline%delete('trail_seed')
        end select
    end subroutine write_cline

    logical function preserves_history( self ) result(preserve)
        class(refine3D_startup), intent(in) :: self
        preserve = self%mode == STARTUP_REUSE_CHAIN
    end function preserves_history

    ! Keep the request separate from the full update that produces a new chain.
    logical function begin_iteration( self, params, cline ) result(full_update)
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        class(refine3D_startup), intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(cmdline),          intent(inout) :: cline
        full_update = self%mode == STARTUP_SEED_CHAIN .and. &
            &params%which_iter == params%startit .and. params%l_update_frac
        if( .not. full_update ) return
        self%requested_fraction = params%update_frac
        self%requested_fractional = params%l_update_frac
        self%requested_trailing = params%l_trail_rec
        self%full_update_active = .true.
        params%update_frac = 1.
        params%l_update_frac = .false.
        params%l_trail_rec   = .false.
        call cline%set('update_frac', 1.)
        ! Reconstruction alone consumes this handshake to publish the full-mass chain.
        call cline%set('trail_seed', 'yes')
    end function begin_iteration

    subroutine finish_iteration( self, params, cline )
        use simple_parameters, only: parameters
        use simple_cmdline,    only: cmdline
        class(refine3D_startup), intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(cmdline),          intent(inout) :: cline
        if( .not. self%full_update_active ) return
        params%update_frac = self%requested_fraction
        params%l_update_frac = self%requested_fractional
        params%l_trail_rec = self%requested_trailing
        call cline%set('update_frac', self%requested_fraction)
        call cline%delete('trail_seed')
        call self%kill
    end subroutine finish_iteration

end module simple_refine3D_startup
