module simple_pftcc_bfacsrch
use simple_opt_spec,          only: opt_spec
use simple_pftcc_opt,         only: pftcc_opt
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_simplex_pftcc_opt, only: simplex_pftcc_opt
use simple_defs               ! use all in there
implicit none

public :: pftcc_bfacsrch
private

type, extends(pftcc_opt) :: pftcc_bfacsrch
    private
    type(opt_spec)                   :: ospec                 !< optimizer specification object
    type(simplex_pftcc_opt)          :: nlopt                 !< optimizer object
    class(polarft_corrcalc), pointer :: pftcc_ptr   =>null()  !< pointer to pftcc object
    integer                          :: reference   = 0       !< reference pft
    integer                          :: particle    = 0       !< particle pft
    integer                          :: rot         = 1       !< in-plane rotation
    integer                          :: nrestarts   = 5       !< simplex restarts (randomized bounds)
  contains
    procedure :: new         => bfacsrch_new
    procedure :: set_indices => bfacsrch_set_indices
    procedure :: set_inipop  => bfacsrch_set_inipop
    procedure :: costfun     => bfacsrch_costfun
    procedure :: minimize    => bfacsrch_minimize
    procedure :: get_nevals  => bfacsrch_get_nevals
    procedure :: get_peaks   => bfacsrch_get_peaks
    procedure :: kill
end type pftcc_bfacsrch

contains

    subroutine bfacsrch_new( self, pftcc, lims, shbarrier, nrestarts, npeaks, maxits, vols )
        use simple_projector, only: projector
        class(pftcc_bfacsrch),              intent(inout) :: self
        class(polarft_corrcalc),    target, intent(in)    :: pftcc
        real,                               intent(in)    :: lims(:,:)
        character(len=*), optional,         intent(in)    :: shbarrier
        integer,          optional,         intent(in)    :: nrestarts, npeaks, maxits
        class(projector), optional, target, intent(in)    :: vols(:)
        self%nrestarts = 10
        if( present(nrestarts) ) self%nrestarts = nrestarts
        ! make optimizer spec
        call self%ospec%specify('simplex', 1, ftol=1e-4,&
        &gtol=1e-4, limits=lims, nrestarts=self%nrestarts)
         ! generate the simplex optimizer object 
        call self%nlopt%new(self%ospec)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
    end subroutine bfacsrch_new

    subroutine bfacsrch_set_indices( self, ref, ptcl, rot, state )
        class(pftcc_bfacsrch), intent(inout) :: self
        integer,               intent(in)    :: ref, ptcl
        integer, optional,     intent(in)    :: rot, state
        self%reference = ref 
        self%particle  = ptcl
        self%rot = 1
        if( present(rot) ) self%rot = rot
    end subroutine bfacsrch_set_indices

    subroutine bfacsrch_set_inipop( self, inipop )
        class(pftcc_bfacsrch), intent(inout) :: self
        real,                  intent(in)    :: inipop(:,:)
        stop 'Not for simplex use; simple_pftcc_bfacsrch%srch_set_inipop'
    end subroutine bfacsrch_set_inipop

    function bfacsrch_costfun( self, vec, D ) result( cost )
        class(pftcc_bfacsrch), intent(inout) :: self
        integer,               intent(in)    :: D
        real,                  intent(in)    :: vec(D)
        real :: cost
        cost = self%pftcc_ptr%euclid(self%reference, self%particle, self%rot, vec(1))
    end function bfacsrch_costfun

    function bfacsrch_minimize( self, irot, shvec, rxy, fromto ) result( cb )
        use simple_math, only: rotmat2d
        class(pftcc_bfacsrch), intent(inout) :: self
        integer, optional,     intent(in)    :: irot
        real,    optional,     intent(in)    :: shvec(:)
        real,    optional,     intent(in)    :: rxy(:)
        integer, optional,     intent(in)    :: fromto(2)
        real              :: cost, cost_init
        real, allocatable :: cb(:)
        allocate(cb(2))
        ! minimisation
        self%ospec%x = 0.
        self%ospec%nevals = 0
        cost_init = self%costfun(self%ospec%x, self%ospec%ndim)
        call self%nlopt%minimize(self%ospec, self, cost)
        if( cost < cost_init )then
            cb(1) = -cost      ! correlation
            cb(2) = self%ospec%x(1)
        else
            cb(1) = -cost_init ! correlation
            cb(2) = 0.
        endif
    end function bfacsrch_minimize

    function bfacsrch_get_nevals( self ) result( nevals )
        class(pftcc_bfacsrch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function bfacsrch_get_nevals

    subroutine bfacsrch_get_peaks( self, peaks )
        class(pftcc_bfacsrch), intent(inout) :: self
        real, allocatable,     intent(out)   :: peaks(:,:)
        allocate(peaks(1,2))
        peaks(1,:) = self%ospec%x
    end subroutine bfacsrch_get_peaks

    subroutine kill( self )
        class(pftcc_bfacsrch), intent(inout) :: self
        self%pftcc_ptr => null()
    end subroutine kill

end module simple_pftcc_bfacsrch
