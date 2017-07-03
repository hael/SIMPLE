module simple_cont3D_de_srch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_projector,        only: projector
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_pftcc_srch,       only: pftcc_srch
use simple_math              ! use all in there
implicit none

public :: cont3D_de_srch
private

integer, parameter :: NDOF           = 5
integer, parameter :: maxits_per_dof = 400
integer, parameter :: MAXITS         = NDOF * maxits_per_dof
integer, parameter :: NINIPOP        = 50
integer, parameter :: NINIPOP_HIGH   = 80 ! < 107
logical, parameter :: debug = .false.

type cont3D_de_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr   => null()  !< polar fourier correlation calculator
    class(projector),        pointer :: vols_ptr(:) => null()  !< volumes for projection
    type(pftcc_srch)                 :: srch_obj               !< shift search object
    type(ori)                        :: o_in                   !< input orientation
    type(ori)                        :: o_out                  !< best orientation found
    type(oris)                       :: o_peaks                !< best orientation found
    logical,             allocatable :: state_exists(:)        !< indicates whether each state is populated
    real                             :: lims(5,2)     = 0.     !< shift search limits
    real                             :: prev_corr     = -1.    !< previous correlation
    real                             :: prev_shift(2) = 0.     !< previous correlation
    real                             :: angthresh     = 0.     !< angular threshold
    real                             :: specscore     = 0.     !< previous spectral score
    real                             :: trs           = 0.     !< shift limit
    integer                          :: npeaks        = 0      !< number of peaks
    integer                          :: iptcl         = 0      !< orientation general index
    integer                          :: ref           = 0      !< previous ref
    integer                          :: prev_ref      = 0      !< previous ref
    integer                          :: state         = 0      !< previous state
    integer                          :: prev_state    = 0      !< previous state
    integer                          :: nstates       = 0      !< number of states
    character(len=STDLEN)            :: shbarr        = 'yes'  !< shift barrier flag
    character(len=STDLEN)            :: refine        = ''     !< 
    logical                          :: exists = .false.

  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PREP ROUTINES
    procedure, private :: prep_srch
    procedure, private :: prep_oris
    procedure, private :: stochastic_weights_L2norm
    ! SEARCH ROUTINES
    procedure          :: exec_srch
    procedure, private :: do_srch
    ! GETTERS/SETTERS
    procedure          :: get_best_ori
    procedure          :: get_softoris
    procedure          :: does_exist
    ! DESTRUCTOR
    procedure          :: kill
end type cont3D_de_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, p, pftcc, vols)
        class(cont3D_de_srch),              intent(inout) :: self  !< instance
        class(params),                   intent(in)    :: p     !< parameters
        class(polarft_corrcalc), target, intent(in)    :: pftcc
        class(projector),        target, intent(in)    :: vols(:)    !< references
        call self%kill
        ! set constants
        self%pftcc_ptr  => pftcc
        self%vols_ptr   => vols
        self%lims(:3,:) = p%eullims
        self%trs        = p%trs
        self%lims(4,:)  = [-self%trs, self%trs]
        self%lims(5,:)  = [-self%trs, self%trs]
        self%nstates    = p%nstates
        if( size(vols).ne.self%nstates )stop 'Inconsistent number of volumes; cont3D_de_srch::new'
        self%angthresh  = p%athres
        self%shbarr     = p%shbarrier
        self%npeaks     = p%npeaks
        self%refine     = p%refine
        self%nstates    = p%nstates
        ! done
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> cont3D_de_srch::CONSTRUCTED NEW SIMPLE_cont3D_de_srch OBJECT'
    end subroutine new

    ! PREP ROUTINES

    !>  \brief  is the master search routine
    subroutine prep_srch(self, a, iptcl, iref, istate)
        class(cont3D_de_srch), intent(inout) :: self
        class(oris),        intent(inout) :: a
        integer,            intent(in)    :: iptcl, iref, istate
        real, allocatable :: frc(:), inipop(:,:)
        self%iptcl      = iptcl
        self%ref        = iref
        self%state      = istate
        self%o_in       = a%get_ori(self%iptcl)
        self%prev_shift = self%o_in%get_shift()
        ! state
        allocate(self%state_exists(self%nstates))
        self%state_exists = a%get_state_exist(self%nstates)
        self%state        = istate
        if( .not.self%state_exists(self%state) )stop 'state is empty; cont3D_de_srch::prep_srch'
        ! target correlation
        call self%vols_ptr(self%state)%fproject_polar(self%ref, self%o_in,&
        &self%pftcc_ptr, serial=.true.)
        self%prev_corr = self%pftcc_ptr%corr(self%ref, self%iptcl, self%state, [0.,0.])
        ! spectral score
        frc = self%pftcc_ptr%genfrc(self%ref, self%iptcl, 1)
        self%specscore = max(0., median_nocopy(frc))
        ! DE search object
        call self%srch_obj%new(self%pftcc_ptr, self%lims, shbarrier=self%shbarr,&
        &npeaks=self%npeaks, maxits=MAXITS, vols=self%vols_ptr)
        call self%srch_obj%set_indices(self%ref, self%iptcl, state=self%state)
        call gen_inipop()
        call self%srch_obj%set_inipop(inipop)
        ! cleanup
        deallocate(frc, inipop)
        if( debug ) write(*,'(A)') '>>> cont3D_de_srch::END OF PREP_SRCH'

        contains

            subroutine gen_inipop()
                type(oris) :: orefs
                type(ori)  :: o
                integer    :: i, npop
                if( self%npeaks == 1 )then
                    npop = NINIPOP
                else
                    npop = NINIPOP
                    ! increases local initial population upon previous failure
                    if( self%o_in%get('ow') > 0.99 )npop = NINIPOP_HIGH
                endif
                allocate( inipop(npop, NDOF) )
                call orefs%rnd_gau_neighbors(npop, self%o_in, self%angthresh, self%lims(:3,:))
                do i = 1, npop
                    if(i == 1)then
                        o = self%o_in
                        call o%set_shift([0.,0.])
                    else
                        o = orefs%get_ori(i)
                        call o%rnd_shift(self%trs)
                    endif
                    inipop(i,:3)  = o%get_euler()   ! angles
                    inipop(i,4:5) = o%get_shift()   ! shift
                enddo
            end subroutine gen_inipop
    end subroutine prep_srch

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine exec_srch( self, a, iptcl, iref, istate )
        class(cont3D_de_srch), intent(inout) :: self
        class(oris),           intent(inout) :: a
        integer,               intent(in)    :: iptcl, iref, istate
        real, allocatable :: peaks(:,:)
        if(nint(a%get(iptcl,'state')) > 0)then
            call self%prep_srch(a, iptcl, iref, istate)
            call self%do_srch(peaks)
            call self%prep_oris(peaks)
            deallocate(peaks)
            call a%set_ori(self%iptcl, self%o_out)
        else
            call a%reject(iptcl)
        endif
        if( debug ) write(*,'(A)') '>>> cont3D_de_srch::END OF SRCH'
    end subroutine exec_srch

    !>  \brief  performs the shift search
    subroutine do_srch( self, peaks )
        class(cont3D_de_srch), intent(inout) :: self
        real, allocatable,     intent(out)   :: peaks(:,:)
        real, allocatable :: solution(:)
        real :: euls(3)
        euls = self%o_in%get_euler()
        solution = self%srch_obj%minimize(rxy=euls)
        call self%srch_obj%get_peaks(peaks)      ! search outcome
        deallocate(solution)
    end subroutine do_srch

    !>  \brief  updates solutions orientations
    subroutine prep_oris( self, peaks )
        use simple_strings, only: int2str
        class(cont3D_de_srch), intent(inout) :: self
        real,                  intent(inout) :: peaks(:,:)
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        type(ori)  :: o
        type(oris) :: os
        real       :: ang_sdev, dist_inpl, wcorr
        real       :: euldist, mi_proj, mi_inpl, mi_state, mi_joint
        integer    :: i, state, roind, prev_state, prev_roind, cnt
        ! UPDATES SCOPE NPEAKS VALUE
        self%npeaks = count(peaks(:,6) > 0.)
        ! init
        call self%o_peaks%new(self%npeaks)
        call os%new(self%npeaks)
        ! o_peaks <- ospec peaks
        cnt = 0
        do i = 1, size(peaks(:,6))
            if( peaks(i,6) <= 0. )cycle
            cnt = cnt + 1
            o   = self%o_in
            ! no in-plane convention as taken care of at online extraction
            call o%set_euler(peaks(i,1:3))
             ! shift addition
            call o%set_shift(peaks(i,4:5) + self%prev_shift)
            call o%set('corr', peaks(i,6))
            call os%set_ori(cnt, o)
        enddo
        ! sorting
        if(self%npeaks > 1)then
            allocate(inds(self%npeaks))
            inds  = (/ (i, i=1,self%npeaks) /)
            corrs = os%get_all('corr')
            call hpsort(self%npeaks, corrs, inds)
            do i = 1, self%npeaks
                o = os%get_ori(inds(i))
                call self%o_peaks%set_ori(i, o)
            enddo
            deallocate(corrs, inds)
        else
            self%o_peaks = os
        endif
        ! spectral score
        call self%o_peaks%set_all2single('specscore', self%specscore)
        ! best orientation & stochastic weights and weighted correlation
        call self%stochastic_weights_L2norm(wcorr)
        self%o_out = self%o_peaks%get_ori(self%npeaks)
        call self%o_out%set('corr', wcorr)  
        ! angular distances & deviation
        euldist   = rad2deg( self%o_in.euldist.self%o_out )
        dist_inpl = rad2deg( self%o_in.inplrotdist.self%o_out )
        ang_sdev  = self%o_peaks%ang_sdev(self%refine, self%nstates, self%npeaks)
        call self%o_out%set('dist', euldist)
        call self%o_out%set('dist_inpl', dist_inpl)
        call self%o_out%set('sdev', ang_sdev)
        ! overlap between distributions
        mi_proj  = 0.
        mi_state = 0.
        if( euldist < 0.25 )then
            mi_proj  = mi_proj + 1.
        endif
        if(self%nstates > 1)then
            state      = nint(self%o_out%get('state'))
            prev_state = nint(self%o_in%get('state'))
            if(prev_state == state)mi_state = mi_state + 1.
        endif
        call self%o_out%set('proj', 0.)
        call self%o_out%set('mi_proj',  mi_proj)
        call self%o_out%set('mi_inpl',  1.)
        call self%o_out%set('mi_state', mi_state)
        call self%o_out%set('mi_joint', 1.)       
        ! frac (unused in convergence)
        call self%o_out%set('frac', 100.)
        if(debug)write(*,*)'simple_cont3D_srch::prep_softoris done'
    end subroutine prep_oris

    !>  \brief  determines and updates stochastic weights
    subroutine stochastic_weights_L2norm( self, wcorr )
        class(cont3D_de_srch), intent(inout) :: self
        real,                    intent(out)   :: wcorr
        real              :: ws(self%npeaks), dists(self%npeaks)
        real, allocatable :: corrs(:)
        type(ori)         :: o
        integer           :: ipeak
        if( self%npeaks == 1 )then
            call self%o_peaks%set(1, 'ow', 1.)
            wcorr = self%o_peaks%get(1, 'corr')
            return
        endif
        ! get correlations
        corrs = self%o_peaks%get_all('corr')
        ! calculate L1 norms
        do ipeak = 1, self%npeaks
            o = self%o_peaks%get_ori(ipeak)
            call self%vols_ptr(self%state)%fproject_polar(self%ref, o,&
            &self%pftcc_ptr, serial=.true.)
            dists(ipeak) = self%pftcc_ptr%euclid(self%ref, self%iptcl, 1)
        end do
        ! calculate normalised weights and weighted corr
        ws    = exp(-dists)
        ws    = ws/sum(ws)
        wcorr = sum(ws*corrs)
        ! update npeaks individual weights
        call self%o_peaks%set_all('ow', ws)
        ! cleanup
        deallocate(corrs)
    end subroutine stochastic_weights_L2norm

    ! GETTERS

    !>  \brief  returns the solution set of orientations
    function get_softoris( self )result( os )
        class(cont3D_de_srch), intent(inout) :: self
        type(oris) :: os
        if(.not.self%exists)stop 'search has not been performed; cont3D_de_srch::get_softoris'
        os = self%o_peaks
    end function get_softoris

    !>  \brief  returns the best solution orientation
    function get_best_ori( self )result( o )
        class(cont3D_de_srch), intent(inout) :: self
        type(ori) :: o
        if(.not.self%exists)stop 'search has not been performed; cont3D_de_srch::get_best_ori'
        o = self%o_out
    end function get_best_ori

    ! GETTERS/SETTERS

    !>  \brief  whether object has been initialized
    logical function does_exist( self )
        class(cont3D_de_srch), intent(inout) :: self
        does_exist = self%exists
    end function does_exist

    ! DESTRUCTOR

    !>  \brief  is the destructor
    subroutine kill( self )
        class(cont3D_de_srch), intent(inout) :: self
        call self%srch_obj%kill
        self%pftcc_ptr => null()
        self%vols_ptr  => null()
        call self%o_in%kill
        call self%o_out%kill
        call self%o_peaks%kill
        if( allocated(self%state_exists) )deallocate(self%state_exists)
        self%exists = .false.
    end subroutine kill

end module simple_cont3D_de_srch