! continuous 3D alignment of band-pass limited polar projections in the Fourier domain
module simple_pftcc_srch
use simple_opt_spec,          only: opt_spec
use simple_pftcc_opt,         only: pftcc_opt
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_projector,         only: projector
use simple_de_pftcc_opt,      only: de_pftcc_opt
use simple_defs               ! use all in there
implicit none

public :: pftcc_srch
private

type, extends(pftcc_opt) :: pftcc_srch
    private
    type(opt_spec)                   :: ospec                  !< optimizer specification object
    type(de_pftcc_opt)               :: nlopt                  !< simplex optimizer object
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
    procedure           :: new         => srch_new
    procedure           :: set_indices => srch_set_indices
    procedure           :: set_inipop  => srch_set_inipop
    procedure           :: costfun     => srch_costfun
    procedure           :: minimize    => srch_minimize
    procedure           :: get_nevals  => srch_get_nevals
    procedure           :: get_peaks   => srch_get_peaks
    procedure, private  :: check_lims
    procedure           :: kill
end type pftcc_srch

contains

    subroutine srch_new( self, pftcc, lims, shbarrier, nrestarts, npeaks, maxits, vols )
        class(pftcc_srch),                  intent(inout) :: self
        class(polarft_corrcalc),    target, intent(in)    :: pftcc
        real,                               intent(in)    :: lims(:,:)
        character(len=*), optional,         intent(in)    :: shbarrier
        integer,          optional,         intent(in)    :: nrestarts, npeaks, maxits
        class(projector), optional, target, intent(in)    :: vols(:)
        real    :: srchlims(3,2)
        integer :: ndim, npeaks_here, maxits_here
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%nrestarts = 1
        ! make optimizer spec
        srchlims = lims
        ndim     = 3
        if( present(nrestarts) )self%nrestarts = nrestarts 
        maxits_here = 1000 * ndim
        if(present(maxits))maxits_here = maxits
        ! de
        if(present(npeaks))then
            npeaks_here = npeaks
            call self%ospec%specify('de', ndim, limits=srchlims,&
            &nrestarts=1, maxits=maxits_here, npeaks=npeaks_here)
        else
            call self%ospec%specify('de', ndim, limits=srchlims,&
            &nrestarts=1, maxits=maxits_here)
        endif
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
    
    !>  \brief  is a setter
    subroutine srch_set_indices( self, ref, ptcl, rot, state )
        class(pftcc_srch), intent(inout) :: self
        integer,           intent(in)    :: ref, ptcl
        integer, optional, intent(in)    :: rot, state
        self%reference = ref 
        self%particle  = ptcl
        self%state     = state
    end subroutine srch_set_indices

    !>  \brief  is a setter
    subroutine srch_set_inipop( self, inipop )
        class(pftcc_srch), intent(inout) :: self
        real,              intent(in)    :: inipop(:,:)
        integer :: i
        call self%ospec%set_inipop(inipop)
        do i = 1, size(self%ospec%inipopulation, dim=1)
            call self%check_lims(self%ospec%inipopulation(i,:))
        enddo
    end subroutine srch_set_inipop

    function srch_costfun( self, vec, D ) result( cost )
        use simple_math, only: enforce_cyclic_limit
        use simple_ori,  only: ori
        class(pftcc_srch), intent(inout) :: self
        integer,           intent(in)    :: D
        real,              intent(in)    :: vec(D)
        type(ori) :: o
        real      :: vec_here(self%ospec%ndim), cost
        integer   :: i
        vec_here = vec
        call self%check_lims(vec_here)
        do i = 1,self%ospec%ndim
            ! euler angles & shift boundaries
            if(vec_here(i) < self%ospec%limits(i,1) .or.&
              &vec_here(i) > self%ospec%limits(i,2))then
                cost = 1.
                return
            endif
        enddo
        ! projection
        call o%new
        call o%set_euler(vec_here)
        call self%vols_ptr(self%state)%fproject_polar(self%reference, o, self%pftcc_ptr,&
        &serial=.true.)
        ! correlation
        cost = -self%pftcc_ptr%corr(self%reference, self%particle, 1)
    end function srch_costfun
    
    function srch_minimize( self, irot, shvec, rxy, fromto ) result( crxy )
        class(pftcc_srch),     intent(inout) :: self
        integer, optional,     intent(in)    :: irot
        real,    optional,     intent(in)    :: shvec(:)
        real,    optional,     intent(in)    :: rxy(:)
        integer, optional,     intent(in)    :: fromto(2)
        real, allocatable :: crxy(:)
        real    :: cost_init, cost
        integer :: i
        logical :: irot_here, shvec_here, rxy_here
        allocate(crxy(4))
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
        ! previous correlation
        cost_init     = self%costfun(self%ospec%x, self%ospec%ndim)
        self%ospec%yb = cost_init
        ! minimisation
        call self%nlopt%minimize(self%ospec, self, cost)
        ! updates best solution
        crxy(1)   = -cost    ! correlation
        crxy(2:4) = self%ospec%x
        call self%check_lims(crxy(2:4))
        ! updates population
        self%ospec%peaks(:,4) = -self%ospec%peaks(:,4)  ! cost to correlation
        ! improvement
        do i = 1, self%ospec%npop
            call self%check_lims(self%ospec%peaks(i,1:3))
        enddo
    end function srch_minimize

    function  srch_get_nevals( self ) result( nevals )
        class(pftcc_srch), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function  srch_get_nevals

    subroutine srch_get_peaks( self, peaks )
        class(pftcc_srch), intent(inout) :: self
        real, allocatable, intent(out)   :: peaks(:,:)
        peaks = self%ospec%peaks
    end subroutine srch_get_peaks

    subroutine check_lims( self, individual )
        use simple_math, only: enforce_cyclic_limit
        class(pftcc_srch), intent(inout) :: self
        real,              intent(inout) :: individual(3)
        call enforce_cyclic_limit(individual(1), 360.)
        call enforce_cyclic_limit(individual(2), 360.)
        call enforce_cyclic_limit(individual(3), 360.)
    end subroutine check_lims

    subroutine kill( self )
        class(pftcc_srch), intent(inout) :: self
        self%pftcc_ptr => null()
        self%vols_ptr  => null()
        call self%ospec%kill
    end subroutine kill

end module simple_pftcc_srch
