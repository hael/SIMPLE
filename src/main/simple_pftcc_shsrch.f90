! rotational origin shift alignment of band-pass limited polar projections in the Fourier domain
module simple_pftcc_shsrch
#include "simple_lib.f08"
use simple_opt_spec,          only: opt_spec
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_opt_factory,       only: opt_factory
use simple_optimizer,         only: optimizer
implicit none

public :: pftcc_shsrch
private

type :: pftcc_shsrch
    private
    type(opt_spec)                   :: ospec                  !< optimizer specification object
    class(optimizer),        pointer :: nlopt                  !< optimizer object
    class(polarft_corrcalc), pointer :: pftcc_ptr   => null()  !< pointer to pftcc object
    integer                          :: reference   =  0       !< reference pft
    integer                          :: particle    =  0       !< particle pft
    integer                          :: nrots       =  0       !< # rotations
    integer                          :: maxits      =  100     !< max # iterations
    integer                          :: nrestarts   =   5      !< simplex restarts (randomized bounds)
    logical                          :: shbarr      =  .true.  !< shift barrier constraint or not
  contains
    procedure :: new
    procedure :: set_indices
    procedure :: minimize
    procedure :: costfun
end type pftcc_shsrch

contains

    !> Shift search constructor
    subroutine new( self, pftcc, lims, lims_init, shbarrier, nrestarts, maxits )
        use simple_projector, only: projector
        class(pftcc_shsrch),             intent(inout) :: self           !< instance
        class(polarft_corrcalc), target, intent(in)    :: pftcc          !< correlator
        real,                            intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,             optional,      intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        character(len=*), optional,      intent(in)    :: shbarrier      !< shift barrier constraint or not
        integer,          optional,      intent(in)    :: nrestarts      !< simplex restarts (randomized bounds)
        integer,          optional,      intent(in)    :: maxits         !< maximum iterations
        type(opt_factory) :: opt_fact
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%nrestarts = 5
        if( present(nrestarts) ) self%nrestarts = nrestarts
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        ! make optimizer spec
        if( present(lims_init) )then
            call self%ospec%specify('simplex', 2, ftol=1e-4, gtol=1e-4, limits=lims,&
                &limits_init=lims_init, nrestarts=self%nrestarts, maxits=self%maxits)
        else
            call self%ospec%specify('simplex', 2, ftol=1e-4, gtol=1e-4,&
                &limits=lims, nrestarts=self%nrestarts, maxits=self%maxits)
        endif
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%nlopt)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! get # rotations
        self%nrots = pftcc%get_nrots()
        ! associate costfun
        self%ospec%costfun => costfun_wrapper
    end subroutine new

    !> shsrch_set_indices Set indicies for shift search
    !! \param ref reference
    !! \param ptcl particle index
    !! \param rot rotational index
    !! \param state current state
    !!
    subroutine set_indices( self, ref, ptcl )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: ref, ptcl
        self%reference = ref
        self%particle  = ptcl
    end subroutine set_indices

    !> Wrapper for cost function (gcc7+)
    function costfun_wrapper(self, vec, D) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real :: cost
        select type (self)
        class is (pftcc_shsrch)
            cost = self%costfun(vec, D)
        class default
            write (*,*) 'error in simple_pftcc_shsrch, costfun_wrapper: unknown type'
            stop
        end select
    end function costfun_wrapper

    !> Cost function
    function costfun( self, vec, D ) result( cost )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: D
        real,                intent(in)    :: vec(D)
        real :: cost
        real :: corrs(self%nrots)
        if( self%shbarr )then
            if( any(vec(:) < self%ospec%limits(:,1)) .or.&
                &any(vec(:) > self%ospec%limits(:,2)) )then
                cost = 1.
                return
            endif
        endif
        call self%pftcc_ptr%gencorrs(self%reference, self%particle, vec, corrs)
        cost = -maxval(corrs)
    end function costfun

    !> minimisation routine
    function minimize( self, irot ) result( cxy )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(out)   :: irot
        real, allocatable :: cxy(:)
        real              :: cost, cost_init, corrs(self%nrots)
        integer           :: loc(1)
        allocate(cxy(3))
        ! minimisation
        self%ospec%x      = 0.
        self%ospec%nevals = 0
        cost_init = self%costfun(self%ospec%x, self%ospec%ndim)
        call self%nlopt%minimize(self%ospec, self, cost)
        if( cost <= cost_init )then
            ! get rotation index
            call self%pftcc_ptr%gencorrs(self%reference, self%particle, self%ospec%x, corrs)
            loc   = maxloc(corrs)
            irot  = loc(1)
            ! set output corr & shift
            cxy(1)  = -cost        ! correlation
            cxy(2:) = self%ospec%x ! shift
            ! rotate the shift vector to the frame of reference
            cxy(2:) = matmul(cxy(2:), rotmat2d(self%pftcc_ptr%get_rot(irot)))
        else
            irot    = 0
            cxy(1)  = -cost_init ! correlation
            cxy(2:) = 0.
        endif
    end function minimize

    function get_nevals( self ) result( nevals )
        class(pftcc_shsrch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function get_nevals

    subroutine get_peaks( self, peaks )
        class(pftcc_shsrch), intent(inout) :: self
        real, allocatable,   intent(out)   :: peaks(:,:) !< output peak matrix
        allocate(peaks(1,2))
        peaks(1,:) = self%ospec%x
    end subroutine get_peaks

end module simple_pftcc_shsrch
