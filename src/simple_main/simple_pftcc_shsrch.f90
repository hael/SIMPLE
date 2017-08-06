!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple optimisation method: Shift search of pftcc objects
module simple_pftcc_shsrch
use simple_opt_spec,          only: opt_spec
use simple_pftcc_opt,         only: pftcc_opt
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_simplex_pftcc_opt, only: simplex_pftcc_opt
use simple_defs               ! use all in there
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
    integer                          :: rot            = 1       !< in-plane rotation
    integer                          :: ldim(3)        = [0,0,0] !< logical dimension of Cartesian image
    integer                          :: maxits         =  100    !< max nr of iterations
    real                             :: maxshift       = 0.      !< maximal shift
    logical                          :: shbarr         = .true.  !< shift barrier constraint or not
    integer                          :: nrestarts      =  5      !< simplex restarts (randomized bounds)
  contains
    procedure :: new         => shsrch_new
    procedure :: set_indices => shsrch_set_indices
    procedure :: set_inipop  => shsrch_set_inipop
    procedure :: costfun     => shsrch_costfun
    procedure :: minimize    => shsrch_minimize
    procedure :: get_nevals  => shsrch_get_nevals
    procedure :: get_peaks   => shsrch_get_peaks
    procedure :: kill
end type pftcc_shsrch

contains

    !> Shift search constructor
    subroutine shsrch_new( self, pftcc, lims, shbarrier, nrestarts, npeaks, maxits, vols )
        use simple_projector, only: projector
        class(pftcc_shsrch),                intent(inout) :: self
        class(polarft_corrcalc),    target, intent(in)    :: pftcc     !< pointer to pftcc object
        real,                               intent(in)    :: lims(:,:) !< logical dimension of Cartesian image
        character(len=*), optional,         intent(in)    :: shbarrier !< shift barrier constraint or not
        integer,          optional,         intent(in)    :: nrestarts !< simplex restarts (randomized bounds)
        integer,          optional,         intent(in)    :: npeaks    !< # peaks
        integer,          optional,         intent(in)    :: maxits    !< maximum iterations
        class(projector), optional, target, intent(in)    :: vols(:)   !< projector volume objects
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
        call self%ospec%specify('simplex', 2, ftol=1e-4, gtol=1e-4,&
            limits=lims, nrestarts=self%nrestarts, maxits=self%maxits)
        ! generate the simplex optimizer object
        call self%nlopt%new(self%ospec)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! get logical dimension
        self%ldim = self%pftcc_ptr%get_ldim()
        ! set maxshift
        self%maxshift = real(maxval(self%ldim))/2.
    end subroutine shsrch_new

    !> shsrch_set_indices Set indicies for shift search
    !! \param ref reference
    !! \param ptcl particle index
    !! \param rot rotational index
    !! \param state current state
    !!
    subroutine shsrch_set_indices( self, ref, ptcl, rot, state )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: ref, ptcl
        integer, optional,   intent(in)    :: rot, state
        self%reference = ref
        self%particle  = ptcl
        self%rot = 1
        if( present(rot) ) self%rot = rot
    end subroutine shsrch_set_indices

    !> shsrch_set_inipop Set init population for shift search
    subroutine shsrch_set_inipop( self, inipop )
        class(pftcc_shsrch), intent(inout) :: self
        real,                intent(in)    :: inipop(:,:) !< initial population
        stop 'Not for simplex use; simple_pftcc_shsrch%srch_set_inipop'
    end subroutine shsrch_set_inipop

    !> Cost function
    function shsrch_costfun( self, vec, D ) result( cost )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: D          !< size of vec
        real,                intent(in)    :: vec(D)     !< input search values
        real    :: vec_here(2)    !< current set of values
        real    :: cost
        vec_here = vec
        where( abs(vec) < 1.e-6 ) vec_here = 0.
        if( self%shbarr )then
            if( vec_here(1) < self%ospec%limits(1,1) .or.&
               &vec_here(1) > self%ospec%limits(1,2) )then
                cost = 1.
                return
            else if( vec_here(2) < self%ospec%limits(2,1) .or.&
                    &vec_here(2) > self%ospec%limits(2,2) )then
                cost = 1.
                return
            endif
        endif
        cost = -self%pftcc_ptr%corr(self%reference, self%particle, self%rot, vec_here)
    end function shsrch_costfun

    !> minimisation routine
    function shsrch_minimize( self, irot, shvec, rxy, fromto ) result( cxy )
        use simple_math, only: rotmat2d
        class(pftcc_shsrch), intent(inout) :: self
        integer, optional,   intent(in)    :: irot        !< index of rotation (obsolete)
        real,    optional,   intent(in)    :: shvec(:)    !< search values vector (obsolete)
        real,    optional,   intent(in)    :: rxy(:)      !< (obsolete)
        integer, optional,   intent(in)    :: fromto(2)   !< (obsolete)
        real, allocatable :: cxy(:)
        real              :: cost, cost_init
        allocate(cxy(3))
        ! minimisation
        self%ospec%x      = 0.
        self%ospec%nevals = 0
        cost_init         = self%costfun(self%ospec%x, self%ospec%ndim)
        call self%nlopt%minimize(self%ospec, self, cost)
        if( cost < cost_init )then
            cxy(1)  = -cost        ! correlation
            cxy(2:) = self%ospec%x ! shift
            ! rotate the shift vector to the frame of reference
            cxy(2:) = matmul(cxy(2:), rotmat2d(self%pftcc_ptr%get_rot(self%rot)))
            if( any(cxy(2:) > self%maxshift) .or. any(cxy(2:) < -self%maxshift) )then
                cxy(1)  = -1.
                cxy(2:) = 0.
            endif
        else
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

    !> Destructor
    subroutine kill( self )
        class(pftcc_shsrch), intent(inout) :: self
        self%pftcc_ptr => null()
    end subroutine kill

end module simple_pftcc_shsrch
