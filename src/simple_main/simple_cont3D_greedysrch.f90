module simple_cont3D_greedysrch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_projector,        only: projector
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_pftcc_srch,      only: pftcc_srch
use simple_math              ! use all in there
implicit none

public :: cont3D_greedysrch
private

logical, parameter :: debug = .false.

type cont3D_greedysrch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr   => null()   !< polar fourier correlation calculator
    class(projector),        pointer :: vols_ptr(:) => null()   !< volumes for projection
    type(pftcc_srch)                 :: srch_obj                !< shift search object
    type(ori)                        :: o_in                    !< input orientation
    type(ori)                        :: o_out                   !< best orientation found
    logical,             allocatable :: state_exists(:)         !< indicates whether each state is populated
    real                             :: lims(2,2)     = 0.     !< shift search limits
    real                             :: prev_corr     = -1.    !< previous correlation
    real                             :: prev_shift(2) = 0.     !< previous correlation
    real                             :: angthresh     = 0.     !< angular threshold
    real                             :: specscore     = 0.     !< previous spectral score
    integer                          :: iptcl         = 0      !< orientation general index
    integer                          :: ref           = 0      !< previous ref
    integer                          :: prev_ref      = 0      !< previous ref
    integer                          :: state         = 0      !< previous state
    integer                          :: prev_state    = 0      !< previous state
    integer                          :: nstates       = 0      !< number of states
    character(len=STDLEN)            :: shbarr        = 'yes'  !< shift barrier flag
    logical                          :: exists = .false.

  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PREP ROUTINES
    procedure, private :: prep_srch
    procedure, private :: prep_ori
    ! SEARCH ROUTINES
    procedure          :: exec_srch
    procedure, private :: do_srch
    ! GETTERS/SETTERS
    procedure          :: get_ori
    procedure          :: does_exist
    ! DESTRUCTOR
    procedure          :: kill
end type cont3D_greedysrch

contains

    !>  \brief  is a constructor
    subroutine new( self, p, pftcc, vols)
        class(cont3D_greedysrch),              intent(inout) :: self  !< instance
        class(params),                   intent(in)    :: p     !< parameters
        class(polarft_corrcalc), target, intent(in)    :: pftcc
        class(projector),        target, intent(in)    :: vols(:)    !< references
        call self%kill
        ! set constants
        self%pftcc_ptr => pftcc
        self%vols_ptr  => vols
        self%lims(:,1) = -p%trs
        self%lims(:,2) = p%trs
        self%nstates   = p%nstates
        if( size(vols).ne.self%nstates )stop 'Inconsistent number of volumes; cont3D_greedysrch::new'
        self%angthresh = p%athres
        self%shbarr    = p%shbarrier
        ! done
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> cont3D_greedysrch::CONSTRUCTED NEW SIMPLE_cont3D_greedysrch OBJECT'
    end subroutine new

    ! PREP ROUTINES

    !>  \brief  is the master search routine
    subroutine prep_srch(self, a, iptcl, iref, istate)
        class(cont3D_greedysrch), intent(inout) :: self
        class(oris),        intent(inout) :: a
        integer,            intent(in)    :: iptcl, iref, istate
        real, allocatable :: frc(:)
        real :: lims(5,2)
        self%iptcl = iptcl
        self%ref   = iref
        self%state = istate
        self%o_in  = a%get_ori(self%iptcl)
        self%prev_shift = self%o_in%get_shift()
        ! state
        allocate(self%state_exists(self%nstates))
        self%state_exists = a%get_state_exist(self%nstates)
        self%state = istate
        if( .not.self%state_exists(self%state) )stop 'state is empty; cont3D_greedysrch::prep_srch'
        ! correlation
        call self%vols_ptr(self%state)%fproject_polar(self%ref, self%o_in,&
        &self%pftcc_ptr, expanded=.true., serial=.true.)
        self%prev_corr = self%pftcc_ptr%corr(self%ref, self%iptcl, 1, [0.,0.])
        ! specscore
        frc = self%pftcc_ptr%genfrc(self%ref, self%iptcl, 1)
        self%specscore = max(0., median_nocopy(frc))
        ! search object
        lims(:3,1) = 0.
        lims(:3,2) = 360.
        lims(4:,:) = self%lims
        call self%srch_obj%new(self%pftcc_ptr, lims, shbarrier=self%shbarr,&
        &vols=self%vols_ptr)
        deallocate(frc)
        if( debug ) write(*,'(A)') '>>> cont3D_greedysrch::END OF PREP_SRCH'
    end subroutine prep_srch

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine exec_srch( self, a, iptcl, iref, istate )
        class(cont3D_greedysrch), intent(inout) :: self
        class(oris),         intent(inout) :: a
        integer,             intent(in)    :: iptcl, iref, istate
        if(nint(a%get(iptcl,'state')) > 0)then
            call self%prep_srch(a, iptcl, iref, istate)
            call self%do_srch
            call self%prep_ori
            call a%set_ori(self%iptcl, self%o_out)
        else
            call a%reject(iptcl)
        endif
        if( debug ) write(*,'(A)') '>>> cont3D_greedysrch::END OF SRCH'
    end subroutine exec_srch

    !>  \brief  performs the shift search
    subroutine do_srch( self )
        class(cont3D_greedysrch), intent(inout) :: self
        real, allocatable :: solution(:)
        real :: euls(3)
        euls = self%o_in%get_euler()
        call self%srch_obj%set_indices(self%ref, self%iptcl, state=self%state)
        solution = self%srch_obj%minimize(rxy=euls)
        ! updates correlation, euler angles and shift
        self%o_out = self%o_in
        if(solution(1) >= self%prev_corr)then
            call self%o_out%set('corr', solution(1))
            call self%o_out%set_euler(solution(2:4))
            call self%o_out%set_shift(self%prev_shift + solution(5:6))
        else
            call self%o_out%set('corr', self%prev_corr)
        endif
        deallocate(solution)
    end subroutine do_srch

    !>  \brief  updates solutions orientations
    subroutine prep_ori( self )
        class(cont3D_greedysrch), intent(inout) :: self
        real      :: u(2), mat(2,2), x1(2), x2(2)
        real      :: frac, euldist, mi_inpl, mi_proj, mi_state, mi_joint
        integer   :: roind, prev_roind
        call self%o_out%set('ow', 1.)
        call self%o_out%set('specscore', self%specscore)
        call self%o_out%set('sdev', 0.)
        call self%o_out%set('proj', 0.)
        call self%o_out%set('state', real(self%state))
        ! dist
        euldist = rad2deg(self%o_in.euldist.self%o_out)
        call self%o_out%set('dist', euldist)
        ! overlap between distributions
        roind      = self%pftcc_ptr%get_roind(360.-self%o_out%e3get())
        prev_roind = self%pftcc_ptr%get_roind(360.-self%o_in%e3get())
        mi_inpl  = 0.
        mi_state = 0.
        mi_joint = 0.
        mi_proj  = 0.
        if(prev_roind == roind)then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif
        ! if(self%nstates > 1)then
        !     if(self%prev_state == self%state)then
        !         mi_state = mi_state + 1.
        !         mi_joint = mi_joint + 1.
        !     endif
        !     mi_joint = mi_joint/3.
        ! else
        !     mi_state = 1.
        !     mi_joint = mi_joint/2.
        ! endif
        call self%o_out%set('mi_proj',  mi_proj)
        call self%o_out%set('mi_inpl',  mi_inpl)
        call self%o_out%set('mi_state', 1.) ! ignored for now
        call self%o_out%set('mi_joint', 1.) ! ignored for now
        ! in-plane distance
        ! make in-plane unit vector
        u = [0., 1.]
        ! calculate previous vec
        mat = rotmat2d(self%o_in%e3get())
        x1  = matmul(u,mat)
        ! calculate new vec
        mat = rotmat2d(self%o_out%e3get())
        x2  = matmul(u, mat)
        call self%o_out%set('dist_inpl', rad2deg(myacos(dot_product(x1,x2))))
        ! frac
        call self%o_out%set('frac', 0.)        
        !frac = 100.*(1.-min(euldist/euldist_thresh, 1.))
        !call self%o_out%set('frac', frac)
        ! done
        if(debug)write(*,*)'simple_cont3D_greedysrch::prep_ori done'
    end subroutine prep_ori

    ! GETTERS/SETTERS

    !>  \brief  returns the best solution orientation
    function get_ori( self )result( o )
        class(cont3D_greedysrch), intent(inout) :: self
        type(ori) :: o
        if(.not.self%exists)stop 'search has not been performed; cont3D_greedysrch::get_ori'
        o = self%o_out
    end function get_ori

    !>  \brief  whether object has been initialized
    logical function does_exist( self )
        class(cont3D_greedysrch), intent(inout) :: self
        does_exist = self%exists
    end function does_exist

    ! DESTRUCTOR

    !>  \brief  is the destructor
    subroutine kill( self )
        class(cont3D_greedysrch), intent(inout) :: self
        call self%srch_obj%kill
        self%pftcc_ptr => null()
        self%vols_ptr  => null()
        call self%o_in%kill
        call self%o_out%kill
        if( allocated(self%state_exists) )deallocate(self%state_exists)
        self%exists = .false.
    end subroutine kill

end module simple_cont3D_greedysrch