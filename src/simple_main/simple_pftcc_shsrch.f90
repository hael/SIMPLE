!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
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
    type(opt_spec)                   :: ospec                 !< optimizer specification object
    type(simplex_pftcc_opt)          :: nlopt                 !< optimizer object
    class(polarft_corrcalc), pointer :: pftcc_ptr   =>null()  !< pointer to pftcc object
    integer                          :: reference   = 0       !< reference pft
    integer                          :: particle    = 0       !< particle pft
    integer                          :: rot         = 1       !< in-plane rotation
    integer                          :: ldim(3)     = [0,0,0] !< logical dimension of Cartesian image
    real                             :: rotmat(2,2) = 0.      !< rotation matrix for checking limits
    real                             :: maxshift    = 0.      !< maximal shift
    logical                          :: shbarr      = .true.  !< shift barrier constraint or not
    integer                          :: nrestarts   =  5      !< simplex restarts (randomized bounds)
  contains
    procedure :: new         => shsrch_new
    procedure :: set_indices => shsrch_set_indices
    procedure :: costfun     => shsrch_costfun
    procedure :: minimize    => shsrch_minimize
    procedure :: get_nevals  => shsrch_get_nevals
    procedure :: kill
end type pftcc_shsrch

contains
    !> shift search method constructor
    subroutine shsrch_new( self, pftcc, lims, shbarrier, nrestarts, vols )
        use simple_projector,        only: projector
        class(pftcc_shsrch),                intent(inout) :: self
        class(polarft_corrcalc),    target, intent(in)    :: pftcc     !< pointer to pftcc object              
        real,                               intent(in)    :: lims(:,:) !< logical dimension of Cartesian image
        character(len=*), optional,         intent(in)    :: shbarrier !< shift barrier constraint or not       
        integer,          optional,         intent(in)    :: nrestarts !< simplex restarts (randomized bounds)  
        class(projector), optional, target, intent(in)    :: vols(:)   !< projector volume objects
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%nrestarts = 5
        if( present(nrestarts) ) self%nrestarts = nrestarts
        ! make optimizer spec
        call self%ospec%specify('simplex', 2, ftol=1e-4,&
        &gtol=1e-4, limits=lims, nrestarts=self%nrestarts)
        ! generate the simplex optimizer object
        call self%nlopt%new(self%ospec)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! get logical dimension
        self%ldim = self%pftcc_ptr%get_ldim()
        ! set maxshift
        self%maxshift = real(maxval(self%ldim))/2.
        ! rotmat init
        self%rotmat      = 0.
        self%rotmat(1,1) = 1.
        self%rotmat(2,2) = 1.
    end subroutine shsrch_new

    subroutine shsrch_set_indices( self, ref, ptcl, rot, state )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: ref, ptcl
        integer, optional,   intent(in)    :: rot, state
        self%reference = ref
        self%particle  = ptcl
        if( present(rot) ) self%rot = rot
    end subroutine shsrch_set_indices

    function shsrch_costfun( self, vec, D ) result( cost )
        class(pftcc_shsrch), intent(inout) :: self
        integer,             intent(in)    :: D          !< size of vec
        real,                intent(in)    :: vec(D)     !< input search values
        real    :: vec_here(2)    !< current set of values
        real    :: rotvec_here(2) !< current set of values rotated to frame of reference
        real    :: cost
        vec_here = vec
        if( abs(vec(1)) < 1e-6 ) vec_here(1) = 0.
        if( abs(vec(2)) < 1e-6 ) vec_here(2) = 0.
        rotvec_here = matmul(vec_here,self%rotmat)
        if( self%shbarr )then
            if( rotvec_here(1) < self%ospec%limits(1,1) .or.&
               &rotvec_here(1) > self%ospec%limits(1,2) )then
                cost = 1.
                return
            else if( rotvec_here(2) < self%ospec%limits(2,1) .or.&
                    &rotvec_here(2) > self%ospec%limits(2,2) )then
                cost = 1.
                return
            endif
        endif
        cost = -self%pftcc_ptr%corr(self%reference, self%particle, self%rot, vec_here)
    end function shsrch_costfun

    !> \todo This function needs fixing - only self is modified
    function shsrch_minimize( self, irot, shvec, rxy, fromto ) result( cxy )
        use simple_math, only: rotmat2d
        class(pftcc_shsrch), intent(inout) :: self
        integer, optional,   intent(in)    :: irot        !< index of rotation (obsolete)
        real,    optional,   intent(in)    :: shvec(:)    !< search values vector (obsolete)
        real,    optional,   intent(in)    :: rxy(:)      !< (obsolete)
        integer, optional,   intent(in)    :: fromto(2)   !< (obsolete)
        real              :: cost, cost_init
        real, allocatable :: cxy(:)
        allocate(cxy(3))
        ! set rotmat for boudary checking and final rotation
        self%rotmat = rotmat2d( self%pftcc_ptr%get_rot(self%rot) )
        ! minimisation
        self%ospec%x = 0.
        self%ospec%nevals = 0
        cost_init = self%costfun(self%ospec%x, self%ospec%ndim)
        call self%nlopt%minimize(self%ospec, self, cost)
        if( cost < cost_init )then
            cxy(1)  = -cost ! correlation
            ! rotate the shift vector to the frame of reference
            cxy(2:) = self%ospec%x ! shift
            cxy(2:) = matmul(cxy(2:),self%rotmat)
            if( any(cxy(2:) > self%maxshift) .or. any(cxy(2:) < -self%maxshift) )then
                cxy(1)  = -1.
                cxy(2:) = 0.
            endif
        else
             cxy(1)  = -cost_init ! correlation
             cxy(2:) = 0.
        endif
        ! clean exit
        self%rotmat      = 0.
        self%rotmat(1,1) = 1.
        self%rotmat(2,2) = 1.
    end function shsrch_minimize

    function shsrch_get_nevals( self ) result( nevals )
        class(pftcc_shsrch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function shsrch_get_nevals

    subroutine kill( self )
        class(pftcc_shsrch), intent(inout) :: self
        self%pftcc_ptr => null()
    end subroutine kill

end module simple_pftcc_shsrch
