module simple_pftcc_srch
use simple_opt_spec,          only: opt_spec
use simple_pftcc_opt,         only: pftcc_opt
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_projector,         only: projector
use simple_simplex_pftcc_opt, only: simplex_pftcc_opt
use simple_defs               ! use all in there
implicit none

public :: pftcc_srch
private

type, extends(pftcc_opt) :: pftcc_srch
    private
    type(opt_spec)                   :: ospec                  !< optimizer specification object
    type(simplex_pftcc_opt)          :: nlopt                  !< optimizer object
    class(polarft_corrcalc), pointer :: pftcc_ptr   => null()  !< pointer to pftcc object
    class(projector),        pointer :: vols_ptr(:) => null()  !< pointer to pftcc object
    integer                          :: reference   =  0       !< reference pft
    integer                          :: particle    =  0       !< particle pft
    integer                          :: state       =  0       !< state & index to vols array
    integer                          :: ldim(3)     =  [0,0,0] !< logical dimension of Cartesian image
    real                             :: maxshift    =  0.      !< maximal shift
    logical                          :: shbarr      =  .true.  !< shift barrier constraint or not
    integer                          :: nrestarts   =  5       !< simplex restarts (randomized bounds)
  contains
    procedure :: new         => srch_new
    procedure :: set_indices => srch_set_indices
    procedure :: costfun     => srch_costfun
    procedure :: minimize    => srch_minimize
    procedure :: get_nevals  => srch_get_nevals
    procedure :: kill
end type pftcc_srch

contains

    subroutine srch_new( self, pftcc, lims, shbarrier, nrestarts, vols )
        class(pftcc_srch),                  intent(inout) :: self
        class(polarft_corrcalc),    target, intent(in)    :: pftcc
        real,                               intent(in)    :: lims(:,:)
        character(len=*), optional,         intent(in)    :: shbarrier
        integer,          optional,         intent(in)    :: nrestarts
        class(projector), optional, target, intent(in)    :: vols(:)
        real :: srchlims(5,2)
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%nrestarts = 5
        if( present(nrestarts) ) self%nrestarts = nrestarts 
        ! make optimizer spec
        srchlims = lims
        call self%ospec%specify('simplex', 5, ftol=1e-4,&
        &gtol=1e-4, limits=srchlims, nrestarts=self%nrestarts)
        ! generate the simplex optimizer object 
        call self%nlopt%new(self%ospec)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! set pointer to volume
        self%vols_ptr => vols
        ! get logical dimension
        self%ldim = self%pftcc_ptr%get_ldim()
        ! set maxshift
        self%maxshift = real(maxval(self%ldim))/2.
    end subroutine srch_new
    
    subroutine srch_set_indices( self, ref, ptcl, rot, state )
        class(pftcc_srch), intent(inout) :: self
        integer,               intent(in)    :: ref, ptcl
        integer, optional,     intent(in)    :: rot, state
        self%reference = ref 
        self%particle  = ptcl
        self%state     = state
    end subroutine srch_set_indices

    function srch_costfun( self, vec, D ) result( cost )
        use simple_math, only: enforce_cyclic_limit
        use simple_ori,  only: ori
        class(pftcc_srch), intent(inout) :: self
        integer,           intent(in)    :: D
        real,              intent(in)    :: vec(D)
        type(ori) :: o
        real      :: vec_here(5)    ! current set of values
        real      :: cost
        integer   :: rot
        vec_here = vec
        ! check so that the in-plane rotation is within the limit
        call enforce_cyclic_limit(vec_here(1), 360.)
        call enforce_cyclic_limit(vec_here(2), 360.)
        call enforce_cyclic_limit(vec_here(3), 360.)
        ! zero small shifts
        if( abs(vec(4)) < 1e-6 ) vec_here(4) = 0.
        if( abs(vec(5)) < 1e-6 ) vec_here(5) = 0.
        ! check shift boundaries
        if( self%shbarr )then
            if(vec_here(4) < self%ospec%limits(4,1) .or. vec_here(4) > self%ospec%limits(4,2))then
                cost = 1.
                return
            else if(vec_here(5) < self%ospec%limits(5,1) .or. vec_here(5) > self%ospec%limits(5,2))then
                cost = 1.
                return
            endif
        endif
        ! projection
        call o%new
        call o%set_euler(vec_here(1:3))
        call self%vols_ptr(self%state)%fproject_polar(self%reference, o, self%pftcc_ptr,&
        &expanded=.true., serial=.true.)
        ! correlation. CTF matrix MUST be created prior to using the minimizer
        cost = -self%pftcc_ptr%corr(self%reference, self%particle, 1, vec_here(4:5))
    end function srch_costfun
    
    function srch_minimize( self, irot, shvec, rxy ) result( crxy )
        use simple_math, only: rotmat2d, enforce_cyclic_limit
        use simple_rnd,  only: ran3
        class(pftcc_srch),     intent(inout) :: self
        integer, optional,     intent(in)    :: irot
        real,    optional,     intent(in)    :: shvec(:)
        real,    optional,     intent(in)    :: rxy(:)
        real, allocatable :: crxy(:)
        real    :: cost_init, cost
        logical :: irot_here, shvec_here, rxy_here
        allocate(crxy(6))
        irot_here  = present(irot)
        shvec_here = present(shvec)
        rxy_here   = present(rxy)
        if( shvec_here ) stop 'shvec is not taken into account; srch_minimize'
        if( irot_here ) stop 'irot is not taken into account; srch_minimize'
        if( .not.rxy_here ) stop 'rxy_here argument strictly required; srch_minimize'
        if( size(rxy).ne.3 )stop 'input rxy is euler angles only; srch_minimize'
        if( self%state == 0 .or. self%state > size(self%vols_ptr) )stop 'Incompatible state; srch_minimize'
        ! initialisation
        self%ospec%x(1:3) = rxy
        self%ospec%x(4:5) = 0.  ! particle is pre-shifted
        ! previous correlation
        cost_init = self%costfun(self%ospec%x, self%ospec%ndim)
        ! minimisation
        call self%nlopt%minimize(self%ospec, self, cost)
        if(cost < cost_init)then
            ! correlation improvement
            crxy(1) = -cost               ! correlation
            crxy(2:4) = self%ospec%x(1:3) ! euler angles
            ! check so that all rotations are within the limit
            call enforce_cyclic_limit(crxy(2), 360.)
            call enforce_cyclic_limit(crxy(3), 360.)
            call enforce_cyclic_limit(crxy(4), 360.)
            ! shifts by vector addition must be done in the driver
            crxy(5:) = self%ospec%x(4:)
            if( any(crxy(5:) > self%maxshift) .or. any(crxy(5:) < -self%maxshift) )then
                crxy(1)  = -1.
                crxy(5:) = 0.
            endif
        else
            ! no improvement
            crxy(1)   = -cost_init  ! previous correlation
            crxy(2:4) = rxy         ! previous euler
            crxy(5:)  = 0.          ! previous shift
        endif
    end function srch_minimize

    function  srch_get_nevals( self ) result( nevals )
        class(pftcc_srch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function  srch_get_nevals

    subroutine kill( self )
        class(pftcc_srch), intent(inout) :: self
        self%pftcc_ptr => null()
        self%vols_ptr  => null()         
    end subroutine kill

end module simple_pftcc_srch
