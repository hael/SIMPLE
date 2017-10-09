! rotational origin shift alignment of band-pass limited polar projections in the Fourier domain
module simple_pftcc_shsrch
#include "simple_lib.f08"
use simple_opt_spec,          only: opt_spec
use simple_pftcc_opt,         only: pftcc_opt
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_simplex_pftcc_opt, only: simplex_pftcc_opt
implicit none

public :: pftcc_shsrch
private

type, extends(pftcc_opt) :: pftcc_shsrch
    private
    type(opt_spec)                   :: ospec                    !< optimizer specification object
    type(simplex_pftcc_opt)          :: nlopt                    !< optimizer object
    class(polarft_corrcalc), pointer :: pftcc_ptr      =>null()  !< pointer to pftcc object
    integer                          :: reference      = 0       !< reference pft
    integer                          :: particle       = 0       !< particle pft
    integer                          :: ldim(3)        = [0,0,0] !< logical dimension of Cartesian image
    integer                          :: nrots          = 0       !< # rotations
    integer                          :: maxits         = 100     !< max # iterations
    real                             :: maxshift       = 0.      !< maximal shift
    logical                          :: shbarr         = .true.  !< shift barrier constraint or not
    integer                          :: nrestarts      =  5      !< simplex restarts (randomized bounds)
  contains
    procedure :: new         => shsrch_new
    procedure :: set_indices => shsrch_set_indices
    procedure :: costfun     => shsrch_costfun
    procedure :: minimize    => shsrch_minimize
end type pftcc_shsrch

contains

    !> Shift search constructor
    subroutine shsrch_new( self, pftcc, lims, lims_init, shbarrier, nrestarts, maxits )
        use simple_projector, only: projector
        class(pftcc_shsrch),                intent(inout) :: self           !< instance
        class(polarft_corrcalc),    target, intent(in)    :: pftcc          !< correlator
        real,                               intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,             optional,         intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        character(len=*), optional,         intent(in)    :: shbarrier      !< shift barrier constraint or not
        integer,          optional,         intent(in)    :: nrestarts      !< simplex restarts (randomized bounds)
        integer,          optional,         intent(in)    :: maxits         !< maximum iterations
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
        ! generate the simplex optimizer object
        call self%nlopt%new(self%ospec)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! get logical dimension
        self%ldim = self%pftcc_ptr%get_ldim()
        ! get # rotations
        self%nrots = pftcc%get_nrots() 
        ! set maxshift
        self%maxshift = real(maxval(self%ldim))/2.
    end subroutine shsrch_new

    !> shsrch_set_indices Set indicies for shift search
    !! \param ref reference
    !! \param ptcl particle index
    !! \param rot rotational index
    !! \param state current state
    !!
    subroutine shsrch_set_indices( self, ref, ptcl )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: ref, ptcl
        self%reference = ref
        self%particle  = ptcl
    end subroutine shsrch_set_indices

    !> Cost function
    function shsrch_costfun( self, vec, D ) result( cost )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: D          !< size of vec
        real,                intent(in)    :: vec(D)     !< input search values
        real :: vec_here(2) !< current set of values
        real :: cost, corrs(self%nrots)
        vec_here = vec
        where( abs(vec) < 1.e-6 ) vec_here = 0.
        if( self%shbarr )then
            if( any(vec_here(:) < self%ospec%limits(:,1)) .or.&
               &any(vec_here(:) > self%ospec%limits(:,2)) )then
                cost = 1.
                return
            endif
        endif
        corrs = self%pftcc_ptr%gencorrs_fft(self%reference, self%particle, vec_here)
        cost = -maxval(corrs)
    end function shsrch_costfun

    !> minimisation routine
    function shsrch_minimize( self, irot ) result( cxy )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(out)   :: irot
        real, allocatable :: cxy(:)
        real              :: cost, cost_init, corrs(self%nrots)
        integer           :: loc(1)
        allocate(cxy(3))
        ! minimisation
        self%ospec%x      = 0.
        self%ospec%nevals = 0
        cost_init         = self%costfun(self%ospec%x, self%ospec%ndim)
        call self%nlopt%minimize(self%ospec, self, cost)
        if( cost <= cost_init )then
            ! get rotation index
            corrs = self%pftcc_ptr%gencorrs_fft(self%reference, self%particle, self%ospec%x)
            loc   = maxloc(corrs)
            irot  = loc(1)
            ! set output corr & shift
            cxy(1)  = -cost        ! correlation
            cxy(2:) = self%ospec%x ! shift
            ! rotate the shift vector to the frame of reference
            cxy(2:) = matmul(cxy(2:), rotmat2d(self%pftcc_ptr%get_rot(irot)))
            if( any(cxy(2:) > self%maxshift) .or. any(cxy(2:) < -self%maxshift) )then
                cxy(1)  = -1.
                cxy(2:) = 0.
            endif
        else
            irot    = 0
            cxy(1)  = -cost_init ! correlation
            cxy(2:) = 0.
        endif
    end function shsrch_minimize

    function shsrch_get_nevals( self ) result( nevals )
        class(pftcc_shsrch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function shsrch_get_nevals

    subroutine shsrch_get_peaks( self, peaks )
        class(pftcc_shsrch), intent(inout) :: self
        real, allocatable,   intent(out)   :: peaks(:,:) !< output peak matrix
        allocate(peaks(1,2))
        peaks(1,:) = self%ospec%x
    end subroutine shsrch_get_peaks

end module simple_pftcc_shsrch
