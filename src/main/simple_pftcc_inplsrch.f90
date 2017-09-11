! in-plane alignment of band-pass limited polar projections in the Fourier domain 
module simple_pftcc_inplsrch
use simple_opt_spec,          only: opt_spec
use simple_pftcc_opt,         only: pftcc_opt
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_simplex_pftcc_opt, only: simplex_pftcc_opt
use simple_defs               ! use all in there
implicit none

public :: pftcc_inplsrch
private

type, extends(pftcc_opt) :: pftcc_inplsrch
    private
    type(opt_spec)                   :: ospec                  !< optimizer specification object
    type(simplex_pftcc_opt)          :: nlopt                  !< optimizer object
    class(polarft_corrcalc), pointer :: pftcc_ptr   => null()  !< pointer to pftcc object
    integer                          :: reference   =  0       !< reference pft
    integer                          :: particle    =  0       !< particle pft
    integer                          :: ldim(3)     =  [0,0,0] !< logical dimension of Cartesian image
    integer                          :: maxits      =  100     !< max nr of iterations
    real                             :: shift_scale =  1.      !< shift scale factor
    real                             :: maxshift    =  0.      !< maximal shift
    logical                          :: shbarr      =  .true.  !< shift barrier constraint or not
    integer                          :: nrestarts   =  5       !< simplex restarts (randomized bounds)
  contains
    procedure :: new         => inplsrch_new
    procedure :: set_indices => inplsrch_set_indices
    procedure :: set_inipop  => inplsrch_set_inipop
    procedure :: costfun     => inplsrch_costfun
    procedure :: minimize    => inplsrch_minimize
    procedure :: get_nevals  => inplsrch_get_nevals
    procedure :: get_peaks   => inplsrch_get_peaks
end type pftcc_inplsrch

contains

    subroutine inplsrch_new( self, pftcc, lims, shbarrier, nrestarts, npeaks, maxits, vols )
        use simple_projector, only: projector
        class(pftcc_inplsrch),              intent(inout) :: self
        class(polarft_corrcalc),    target, intent(in)    :: pftcc
        real,                               intent(in)    :: lims(:,:)
        character(len=*), optional,         intent(in)    :: shbarrier
        integer,          optional,         intent(in)    :: nrestarts, npeaks, maxits
        class(projector), optional, target, intent(in)    :: vols(:)
        real :: inpllims(3,2)
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
        inpllims(1,1)  = 0.
        inpllims(1,2)  = 360.
        inpllims(2:,:) = lims
        call self%ospec%specify('simplex', 3, ftol=1e-4, gtol=1e-4,&
            limits=inpllims, nrestarts=self%nrestarts, maxits=self%maxits)
        ! generate the simplex optimizer object 
        call self%nlopt%new(self%ospec)
        ! set pointer to corrcalc object
        self%pftcc_ptr => pftcc
        ! get logical dimension
        self%ldim = self%pftcc_ptr%get_ldim()
        ! set maxshift
        self%maxshift = real(maxval(self%ldim))/2.
    end subroutine inplsrch_new
    
    subroutine inplsrch_set_indices( self, ref, ptcl, rot, state )
        class(pftcc_inplsrch), intent(inout) :: self
        integer,               intent(in)    :: ref, ptcl
        integer, optional,     intent(in)    :: rot, state
        self%reference = ref 
        self%particle  = ptcl
    end subroutine inplsrch_set_indices

    !>  \brief  is a setter
    subroutine inplsrch_set_inipop( self, inipop )
        class(pftcc_inplsrch), intent(inout) :: self
        real,                  intent(in)    :: inipop(:,:)
        stop 'Not for simplex use; simple_pftcc_inplsrch%srch_set_inipop'
    end subroutine inplsrch_set_inipop

    function inplsrch_costfun( self, vec, D ) result( cost )
        use simple_math, only: enforce_cyclic_limit
        class(pftcc_inplsrch), intent(inout) :: self
        integer,               intent(in)    :: D
        real,                  intent(in)    :: vec(D)
        real    :: vec_here(3)    ! current set of values
        real    :: shift(2), cost
        integer :: rot
        vec_here = vec
        ! check so that the in-plane rotation is within the limit
        call enforce_cyclic_limit(vec_here(1), 360.)
        ! check boundary
        if( self%shbarr )then
            if( any(vec_here(2:3) < self%ospec%limits(2:3,1)) .or.&
               &any(vec_here(2:3) > self%ospec%limits(2:3,2)) )then
                cost = 1.
                return
            endif
        endif
        ! unscale shift
        shift = vec_here(2:3) / self%shift_scale
        ! zero small shift
        where( abs(shift) < 1.e-6 ) shift = 0.
        ! cost
        rot  =  self%pftcc_ptr%get_roind(vec_here(1))
        cost = -self%pftcc_ptr%corr(self%reference, self%particle, rot, shift)
    end function inplsrch_costfun
    
    function inplsrch_minimize( self, irot, shvec, rxy, fromto ) result( crxy )
        use simple_math, only: rotmat2d, enforce_cyclic_limit
        use simple_rnd,  only: ran3
        class(pftcc_inplsrch), intent(inout) :: self
        integer, optional,     intent(in)    :: irot
        real,    optional,     intent(in)    :: shvec(:)
        real,    optional,     intent(in)    :: rxy(:)
        integer, optional,     intent(in)    :: fromto(2)
        real, allocatable :: crxy(:)
        logical           :: irot_here, shvec_here, rxy_here
        allocate(crxy(4))         !! fixed size, not checking
        irot_here  = present(irot)
        shvec_here = present(shvec)
        rxy_here   = present(rxy)
        if( irot_here  .and. rxy_here ) stop 'irot & rxy cannot be simultaneously present; inplsrch_minimize'
        if( shvec_here .and. rxy_here ) stop 'shvec_here & rxy_here cannot be simultaneously present; inplsrch_minimize'
        ! initialisation
        self%ospec%x      = 0.
        self%ospec%x(1)   = ran3()*360.
        self%ospec%nevals = 0
        if( rxy_here   ) self%ospec%x(1:3) = rxy(1:3)
        if( irot_here  ) self%ospec%x(1)   = self%pftcc_ptr%get_rot(irot)
        if( shvec_here ) self%ospec%x(2:3) = shvec(1:2)
        ! determines & applies shift scaling
        self%shift_scale = (self%ospec%limits(1,2) - self%ospec%limits(1,1)) /&
        &( maxval(self%ospec%limits(2:,2)) - minval(self%ospec%limits(2:,1)) )
        self%ospec%limits(2:,:) = self%ospec%limits(2:,:) * self%shift_scale
        self%ospec%x(2:3)       = self%ospec%x(2:3)       * self%shift_scale
        ! minimisation
        call self%nlopt%minimize(self%ospec, self, crxy(1))
        ! solution and un-scaling
        crxy(1)           = -crxy(1)                             ! correlation
        crxy(2)           = self%ospec%x(1)                      ! in-plane angle
        self%ospec%x(2:3) = self%ospec%x(2:3) / self%shift_scale ! shift unscaling
        crxy(3:4)         = self%ospec%x(2:3)                    ! shift
        self%ospec%limits(2:3,:) = self%ospec%limits(2:3,:) / self%shift_scale
        ! check so that the in-plane rotation is within the limit
        call enforce_cyclic_limit(crxy(2), 360.)
        ! rotate the shift vector to the frame of reference
        crxy(3:) = matmul(crxy(3:), rotmat2d(crxy(2)))
        if( any(crxy(3:) > self%maxshift) .or. any(crxy(3:) < -self%maxshift) )then
            crxy(1)  = -1.
            crxy(3:) = 0.
        endif
    end function inplsrch_minimize

    function  inplsrch_get_nevals( self ) result( nevals )
        class(pftcc_inplsrch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function  inplsrch_get_nevals

    subroutine inplsrch_get_peaks( self, peaks )
        class(pftcc_inplsrch), intent(inout) :: self
        real, allocatable,     intent(out)   :: peaks(:,:)
        allocate(peaks(1,3))   !! fixed size, not checking
        peaks(1,:) = self%ospec%x
    end subroutine inplsrch_get_peaks
    
end module simple_pftcc_inplsrch
