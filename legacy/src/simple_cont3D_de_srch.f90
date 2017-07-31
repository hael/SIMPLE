module simple_cont3D_de_srch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_projector,        only: projector
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_sym,              only: sym
use simple_pftcc_srch,       only: pftcc_srch
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_math              ! use all in there
implicit none

public :: cont3D_de_srch
private

real,    parameter :: FACTWEIGHTS_THRESH = 0.001    !< threshold for factorial weights
integer, parameter :: NDOF           = 3
integer, parameter :: maxits_per_dof = 300
integer, parameter :: MAXITS         = NDOF * maxits_per_dof
integer, parameter :: NINIPOP        = 103 ! <= 103
logical, parameter :: debug = .false.

type cont3D_de_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr   => null()  !< polar fourier correlation calculator
    class(projector),        pointer :: vols_ptr(:) => null()  !< volumes for projection
    type(pftcc_shsrch)               :: shsrch_obj             !< origin shift search object
    type(pftcc_srch)                 :: srch_obj               !< shift search object
    type(sym)                        :: se
    type(ori)                        :: o_in                   !< input orientation
    type(ori)                        :: o_out                  !< best orientation found
    type(oris)                       :: o_peaks                !< best orientation found
    logical,             allocatable :: state_exists(:)        !< indicates whether each state is populated
    real                             :: lims(5,2)     = 0.     !< shift search limits
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
    procedure, private :: prep_inipop
    procedure, private :: prep_oris
    procedure, private :: stochastic_weights
    ! SEARCH ROUTINES
    procedure          :: exec_srch
    procedure, private :: euler_srch
    procedure, private :: shift_srch
    ! GETTERS/SETTERS
    procedure          :: get_best_ori
    procedure          :: get_peaks
    procedure          :: does_exist
    ! DESTRUCTOR
    procedure          :: kill
end type cont3D_de_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, p, pftcc, vols, se )
        class(cont3D_de_srch),           intent(inout) :: self     !< instance
        class(params),                   intent(in)    :: p        !< parameters
        class(polarft_corrcalc), target, intent(in)    :: pftcc    !< corrcalc obj
        class(projector),        target, intent(in)    :: vols(:)  !< references
        class(sym),                      intent(in)    :: se
        call self%kill
        self%pftcc_ptr  => pftcc
        self%vols_ptr   => vols
        self%lims(:3,:) = p%eullims
        self%trs        = p%trs
        self%lims(4,:)  = [-self%trs, self%trs]
        self%lims(5,:)  = self%lims(4,:)
        self%nstates    = p%nstates
        if( size(vols).ne.self%nstates )stop 'Inconsistent number of volumes; cont3D_de_srch::new'
        self%angthresh  = p%athres
        self%shbarr     = p%shbarrier
        self%npeaks     = 2 * p%npeaks
        self%refine     = p%refine
        self%nstates    = p%nstates
        self%se         = se
        self%exists     = .true.
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
        call self%o_peaks%new(self%npeaks)
        ! state
        allocate(self%state_exists(self%nstates))
        self%state_exists = a%get_state_exist(self%nstates)
        self%state        = istate
        if( .not.self%state_exists(self%state) )stop 'state is empty; cont3D_de_srch::prep_srch'
        ! target correlation
        call self%vols_ptr(self%state)%fproject_polar(self%ref, self%o_in,&
        &self%pftcc_ptr, serial=.true.)
        ! spectral score
        frc = self%pftcc_ptr%genfrc(self%ref, self%iptcl, 1)
        self%specscore = max(0., median_nocopy(frc))
        ! DE search object
        call self%srch_obj%new(self%pftcc_ptr, self%lims(:3,:), shbarrier=self%shbarr,&
        &maxits=MAXITS, vols=self%vols_ptr)
        call self%srch_obj%set_indices(self%ref, self%iptcl, state=self%state)
        call self%prep_inipop( inipop )
        call self%srch_obj%set_inipop( inipop )
        ! shift search object
        call self%shsrch_obj%new(self%pftcc_ptr, self%lims(4:5,:), shbarrier=self%shbarr)
        ! cleanup
        deallocate(frc, inipop)
        if( debug ) write(*,'(A)') '>>> cont3D_de_srch::END OF PREP_SRCH'
    end subroutine prep_srch

    subroutine prep_inipop( self, inipop )
        use simple_rnd, only: gasdev, ran3
        class(cont3D_de_srch), intent(inout) :: self
        real, allocatable,     intent(out)   :: inipop(:,:)
        type(ori)  :: o, o_transform
        real       :: euler(3)
        integer    :: i
        allocate( inipop(NINIPOP, NDOF) )
        ! previous best
        inipop(1,:) = self%o_in%get_euler()
        ! others...
        call o_transform%new
        do i = 2, NINIPOP
            euler    = 0.
            euler(2) = gasdev(0., 2. * self%angthresh) ! double angular threshold!!
            call o_transform%set_euler( euler )
            o = self%o_in.compose.o_transform
            euler(3) = self%o_in%e3get() + (ran3()-0.5)*90.
            call enforce_cyclic_limit(euler(3), 360.)
            call o%e3set( euler(3) )
            call self%se%rot_to_asym( o )
            ! set individual
            inipop(i,:) = o%get_euler()
        enddo
    end subroutine prep_inipop

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine exec_srch( self, a, iptcl, iref, istate )
        class(cont3D_de_srch), intent(inout) :: self
        class(oris),           intent(inout) :: a
        integer,               intent(in)    :: iptcl, iref, istate
        if(nint(a%get(iptcl,'state')) > 0)then
            call self%prep_srch(a, iptcl, iref, istate)
            call self%euler_srch
            call self%shift_srch
            call self%prep_oris
            call a%set_ori(self%iptcl, self%o_out)
        else
            call a%reject(iptcl)
        endif
        if( debug ) write(*,'(A)') '>>> cont3D_de_srch::END OF SRCH'
    end subroutine exec_srch

    !>  \brief  performs the shift search
    subroutine euler_srch( self )
        class(cont3D_de_srch), intent(inout) :: self
        type(ori)            :: o
        real,    allocatable :: peaks(:,:), corrs(:), solution(:)
        integer, allocatable :: inds(:)
        real                 :: euls(3)
        integer              :: ipeak, cnt, ind, n
        euls     = self%o_in%get_euler()
        solution = self%srch_obj%minimize(rxy=euls)
        call self%srch_obj%get_peaks(peaks)      ! search outcome
        deallocate(solution)
        ! o_peaks <- ospec peaks
        n     = size(peaks(:,1))
        corrs = peaks(:,NDOF+1)
        inds  = (/(ipeak,ipeak=1,n)/)
        call hpsort(n, corrs, inds)
        cnt = 0
        do ipeak = n, n-self%npeaks+1, -1
            ind = inds(ipeak)
            cnt = cnt + 1
            o   = self%o_in
            call o%set_euler(peaks(ind,1:3))
            call o%set('corr', peaks(ind,NDOF+1))
            call self%o_peaks%set_ori(cnt, o)
        enddo
    end subroutine euler_srch

    !>  \brief  performs the shift search
    subroutine shift_srch( self )
        class(cont3D_de_srch), intent(inout) :: self
        real, allocatable :: cxy(:)
        type(ori) :: o
        real      :: cc
        integer   :: i, inpl_ind, state
        do i = 1, self%npeaks
                o = self%o_peaks%get_ori( i )
                call o%e3set(0.)
                cc = o%get('corr')
                inpl_ind = self%pftcc_ptr%get_roind(360.-o%e3get())
                call self%vols_ptr(self%state)%fproject_polar(1, o, self%pftcc_ptr, serial=.true.)
                call self%shsrch_obj%set_indices(1, self%iptcl, inpl_ind)
                cxy = self%shsrch_obj%minimize()
                if( cxy(1) >= cc )then
                    call self%o_peaks%set(i, 'corr', cxy(1))
                    call self%o_peaks%set(i, 'x', self%prev_shift(1) + cxy(2))
                    call self%o_peaks%set(i, 'y', self%prev_shift(2) + cxy(3))
                endif
            end do
    end subroutine shift_srch

    !>  \brief  updates solutions orientations
    subroutine prep_oris( self )
        class(cont3D_de_srch), intent(inout) :: self
        real,    allocatable :: corrs(:)
        type(ori)  :: o
        real       :: euldist, mi_proj, mi_state, ang_sdev, dist_inpl, wcorr
        integer    :: i, state, prev_state, cnt, best_loc(1)
        ! stochastic weights and weighted correlation
        call self%stochastic_weights(wcorr)
        ! best orientation
        corrs      = self%o_peaks%get_all('corr')
        best_loc   = maxloc(corrs)
        self%o_out = self%o_peaks%get_ori(best_loc(1))
        call self%o_out%set('corr', wcorr)
        ! angular distances & deviation
        o = self%o_out
        call self%se%sym_euldist(self%o_in, o, euldist)
        dist_inpl = rad2deg( self%o_in.inplrotdist.o )
        ang_sdev  = self%o_peaks%ang_sdev(self%refine, self%nstates, self%npeaks)
        call self%o_out%set('dist', euldist)
        call self%o_out%set('dist_inpl', dist_inpl)
        call self%o_out%set('sdev', ang_sdev)
        ! spectral score
        call self%o_out%set('specscore', self%specscore)
        ! overlap between distributions
        mi_proj  = 0.
        mi_state = 0.
        if( euldist < 1.0 )then
            mi_proj  = mi_proj + 1.
        endif
        if(self%nstates > 1)then
            state      = nint(self%o_out%get('state'))
            prev_state = nint(self%o_in%get('state'))
            if(prev_state == state)mi_state = mi_state + 1.
        endif
        if(self%o_in%isthere('proj'))     call self%o_out%set('proj', 0.)
        if(self%o_in%isthere('mi_inpl'))  call self%o_out%set('mi_inpl',  1.)
        if(self%o_in%isthere('mi_joint')) call self%o_out%set('mi_joint', 1.)       
        call self%o_out%set('mi_proj',  mi_proj)
        call self%o_out%set('mi_state', mi_state)
        call self%o_out%set('frac', 100.)
        if(debug)write(*,*)'simple_cont3D_srch::prep_oris done'
    end subroutine prep_oris

    !>  \brief  determines, thresholds and updates stochastic weights
    subroutine stochastic_weights( self, wcorr )
        class(cont3D_de_srch), intent(inout) :: self
        real,                  intent(out)   :: wcorr
        type(oris)        :: os
        real, allocatable :: corrs(:)
        real              :: logws(self%npeaks), ws(self%npeaks)
        integer           :: order(self%npeaks), ipeak, cnt
        logical           :: included(self%npeaks)
        if( self%npeaks == 1 )then
            call self%o_peaks%set(1, 'ow', 1.)
            wcorr = self%o_peaks%get(1, 'corr')
            return
        endif
        corrs = self%o_peaks%get_all('corr')
        ws    = exp(-(1.-corrs))
        logws = log(ws)
        order = (/(ipeak,ipeak=1,self%npeaks)/)
        call hpsort(self%npeaks, logws, order)
        call reverse(order)
        call reverse(logws)
        forall(ipeak = 1:self%npeaks) ws(order(ipeak)) = exp(sum(logws(:ipeak))) 
        ! thresholding of the weights
        included = (ws >= FACTWEIGHTS_THRESH)
        where( .not.included ) ws = 0.
        ! weighted correlation
        wcorr = sum(ws*corrs, mask=included) / sum(ws, mask=included)
        call self%o_peaks%set_all('ow', ws)
        ! removes zero weighted peaks and update global npeaks
        os = self%o_peaks
        self%npeaks = count(ws > TINY)
        call self%o_peaks%new(self%npeaks)
        cnt  = 0
        do ipeak = 1, os%get_noris()
            if( included(ipeak) )then
                cnt = cnt + 1
                call self%o_peaks%set_ori(cnt, os%get_ori(ipeak))
            endif
        enddo
    end subroutine stochastic_weights

    ! GETTERS

    !>  \brief  returns the solution set of orientations
    function get_peaks( self )result( os )
        class(cont3D_de_srch), intent(inout) :: self
        type(oris) :: os
        if(.not.self%exists)stop 'search has not been performed; cont3D_de_srch::get_softoris'
        os = self%o_peaks
    end function get_peaks

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