module simple_cont3D_srch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_math              ! use all in there
implicit none

public :: cont3D_srch
private

real,    parameter :: E3HALFWINSZ = 90. !< in-plane angle half window size
logical, parameter :: debug = .false.

type cont3D_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr => null() !< polar fourier correlation calculator
    type(pftcc_shsrch)               :: shsrch_obj          !< shift search object
    type(oris)                       :: reforis             !< per ptcl search space and result of euler angles search
    type(oris)                       :: softoris            !< references returned
    type(oris)                       :: shiftedoris         !< references whose shifts are searched
    type(ori)                        :: o_in                !< input orientation
    type(ori)                        :: o_out               !< best orientation found
    logical,             allocatable :: state_exists(:)     !< indicates whether each state is populated
    real                    :: lims(2,2)  = 0.     !< shift search limits
    real                    :: prev_corr  = -1.    !< previous correlation
    real                    :: angthresh  = 0.     !< angular threshold
    real                    :: specscore  = 0.     !< previous spectral score
    integer                 :: iptcl      = 0      !< orientation general index
    integer                 :: prev_ref   = 0      !< previous reference
    integer                 :: prev_roind = 0      !< previous in-plane rotational index
    integer                 :: prev_state = 0      !< previous state
    integer                 :: nstates    = 0      !< number of states
    integer                 :: nbetter    = 0      !< number of improving references found
    integer                 :: neval      = 0      !< number of references evaluated
    integer                 :: npeaks     = 0      !< number of references returned
    integer                 :: nrefs      = 0      !< number of references in search space
    integer                 :: nrots      = 0      !< number of in-plane rotations
    character(len=STDLEN)   :: refine = ''
    logical                 :: exists = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PREP ROUTINES
    procedure, private :: prep_srch
    procedure, private :: prep_softoris
    ! SEARCH ROUTINES
    procedure          :: exec_srch
    procedure, private :: do_euler_srch
    procedure, private :: do_shift_srch
    ! GETTERS
    procedure          :: get_best_ori
    procedure          :: get_softoris
    procedure          :: does_exist
    ! DESTRUCTOR
    procedure          :: kill
end type cont3D_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, p, e, pftcc)
        class(cont3D_srch),             intent(inout) :: self  !< instance
        class(params),                   intent(in)    :: p     !< parameters
        class(oris),                     intent(in)    :: e     !< references
        class(polarft_corrcalc), target, intent(in)    :: pftcc
        call self%kill
        ! set constants
        self%pftcc_ptr    => pftcc
        self%reforis   = e
        self%npeaks    = p%npeaks
        self%lims(:,1) = -p%trs
        self%lims(:,2) =  p%trs
        self%angthresh = p%athres
        self%nstates   = p%nstates
        self%nrefs     = self%pftcc_ptr%get_nrefs()
        self%nrots     = self%pftcc_ptr%get_nrots()
        if(self%reforis%get_noris().ne.self%nrefs)&
            &stop 'Inconsistent number of references & orientations'
        ! done
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> cont3D_srch::CONSTRUCTED NEW SIMPLE_cont3D_srch OBJECT'
    end subroutine new

    ! PREP ROUTINES

    !>  \brief  is the master search routine
    subroutine prep_srch(self, a, iptcl)
        class(cont3D_srch), intent(inout) :: self
        class(oris),         intent(inout) :: a
        integer,             intent(in)    :: iptcl
        real, allocatable :: frc(:)
        if(iptcl == 0)stop 'ptcl index mismatch; cont3D_srch::do_srch'
        self%iptcl      = iptcl
        self%o_in       = a%get_ori(self%iptcl)
        self%prev_state = nint(self%o_in%get('state'))
        allocate(self%state_exists(self%nstates))
        self%state_exists = a%get_state_exist(self%nstates)
        if( .not.self%state_exists(self%prev_state) )stop 'state is empty; cont3D_srch::prep_srch'
        self%prev_roind = self%pftcc_ptr%get_roind(360. - self%o_in%e3get())
        self%prev_ref   = self%reforis%find_closest_proj(self%o_in) ! state not taken into account
        self%prev_corr  = self%pftcc_ptr%corr(self%prev_ref, self%iptcl, self%prev_roind)
        ! specscore
        frc = self%pftcc_ptr%genfrc(self%prev_ref, self%iptcl, self%prev_roind)
        self%specscore  = max(0., median_nocopy(frc))
        ! shift search object
        call self%shsrch_obj%new(self%pftcc_ptr, self%lims, nrestarts=3)
        deallocate(frc)
        if( debug ) write(*,'(A)') '>>> cont3D_srch::END OF PREP_SRCH'
    end subroutine prep_srch

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine exec_srch( self, a, iptcl )
        class(cont3D_srch), intent(inout) :: self
        class(oris),        intent(inout) :: a
        integer,            intent(in)    :: iptcl
        if(nint(a%get(iptcl,'state')) > 0)then
            ! INIT
            call self%prep_srch(a, iptcl)
            ! EULER ANGLES SEARCH
            call self%do_euler_srch
            ! SHIFT SEARCH
            call self%do_shift_srch
            ! outcome prep
            call self%prep_softoris
            ! update
            call a%set_ori(self%iptcl, self%o_out)
        else
            call a%reject(iptcl)
        endif
        if( debug ) write(*,'(A)') '>>> cont3D_srch::END OF SRCH'
    end subroutine exec_srch

    !>  \brief  performs euler angles search
    subroutine do_euler_srch( self )
        class(cont3D_srch), intent(inout) :: self
        integer, allocatable :: roind_vec(:)    ! slice of in-plane angles
        real                 :: inpl_corr
        integer              :: iref
        ! init
        self%nbetter = 0
        self%neval   = 0
        roind_vec    = self%pftcc_ptr%get_win_roind(360.-self%o_in%e3get(), E3HALFWINSZ)
        ! search
        ! the input search space in stochastic, so no need for a randomized search order
        do iref = 1,self%nrefs
            if(iref == self%prev_ref)cycle
            call greedy_inpl_srch(iref, inpl_corr)
            self%neval = self%neval+1
            if(inpl_corr >= self%prev_corr)self%nbetter = self%nbetter+1
            if(self%nbetter >= self%npeaks)exit
        enddo
        if( self%nbetter<self%npeaks )then
            ! previous reference considered last
            call greedy_inpl_srch(self%prev_ref, inpl_corr)
            self%nbetter = self%nbetter + 1
            self%neval   = self%nrefs
        endif            
        deallocate(roind_vec)
        if(debug)write(*,*)'simple_cont3D_srch::do_refs_srch done'

        contains

            subroutine greedy_inpl_srch(iref_here, corr_here)
                integer, intent(in)    :: iref_here
                real,    intent(inout) :: corr_here
                real    :: corrs(self%nrots), e3
                integer :: loc(1), inpl_ind
                corrs     = self%pftcc_ptr%gencorrs(iref_here, self%iptcl, roind_vec=roind_vec)
                loc       = maxloc(corrs)
                inpl_ind  = loc(1)
                corr_here = corrs(inpl_ind)
                e3 = 360. - self%pftcc_ptr%get_rot(inpl_ind)
                call self%reforis%e3set(iref_here, e3)
                call self%reforis%set(iref_here, 'corr', corr_here)
            end subroutine greedy_inpl_srch
            
    end subroutine do_euler_srch

    !>  \brief  performs the shift search
    subroutine do_shift_srch( self )
        use simple_pftcc_shsrch      ! use all in there
        class(cont3D_srch), intent(inout) :: self
        real, allocatable :: corrs(:)
        type(ori)         :: o
        real, allocatable :: cxy(:)
        real              :: prev_shift_vec(2), corr
        integer           :: ref_inds(self%nrefs), i, iref, inpl_ind, cnt
        call self%shiftedoris%new(self%npeaks)
        ! SORT ORIENTATIONS
        ref_inds = (/ (i, i=1, self%nrefs) /)
        corrs    = self%reforis%get_all('corr')
        call hpsort(self%nrefs, corrs, ref_inds)
        ! SHIFT SEARCH
        prev_shift_vec = self%o_in%get_shift()
        cnt = 0
        do i = self%nrefs-self%npeaks+1,self%nrefs
            cnt      = cnt+1
            iref     = ref_inds(i)
            o        = self%reforis%get_ori(iref)
            corr     = o%get('corr')
            inpl_ind = self%pftcc_ptr%get_roind(360.-o%e3get())
            call self%shsrch_obj%set_indices(iref, self%iptcl, inpl_ind)
            cxy = self%shsrch_obj%minimize()
            if(cxy(1) >= corr)then
                call o%set('corr', cxy(1))
                call o%set_shift(cxy(2:) + prev_shift_vec)
            else
                call o%set_shift(prev_shift_vec)
            endif
            call self%shiftedoris%set_ori(cnt, o)
        enddo
        ! done        
        deallocate(corrs)
        if(debug)write(*,*)'simple_cont3d_src::do_shift_srch done'
    end subroutine do_shift_srch

    !>  \brief  updates solutions orientations
    subroutine prep_softoris( self )
        class(cont3D_srch), intent(inout) :: self
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        type(ori) :: o
        real      :: euldist_thresh, ang_sdev, dist_inpl
        real      :: frac, wcorr, euldist, mi_proj, mi_inpl, mi_state, mi_joint
        integer   :: i, state, roind, prev_state
        call self%softoris%new(self%npeaks)
        ! sort oris
        if(self%npeaks > 1)then
            allocate(inds(self%npeaks))
            inds  = (/ (i, i=1,self%npeaks) /)
            corrs = self%shiftedoris%get_all('corr')
            call hpsort(self%npeaks, corrs, inds)
            do i = 1,self%npeaks
                o = self%shiftedoris%get_ori(inds(i))
                call self%softoris%set_ori(i, o)
            enddo
            deallocate(corrs, inds)
        else
            self%softoris = self%shiftedoris
        endif
        ! stochastic weights and weighted correlation
        if(self%npeaks > 1)then
            call self%softoris%stochastic_weights( wcorr )
        else
            call self%softoris%set(1, 'ow', 1.)
            wcorr = self%softoris%get(1, 'corr')
        endif
        ! variables inherited from input ori: CTF, w, lp
        if(self%o_in%isthere('dfx'))   call self%softoris%set_all2single('dfx',self%o_in%get('dfx'))
        if(self%o_in%isthere('dfy'))   call self%softoris%set_all2single('dfy',self%o_in%get('dfy'))
        if(self%o_in%isthere('angast'))call self%softoris%set_all2single('angast',self%o_in%get('angast'))
        if(self%o_in%isthere('cs'))    call self%softoris%set_all2single('cs',self%o_in%get('cs'))
        if(self%o_in%isthere('kv'))    call self%softoris%set_all2single('kv',self%o_in%get('kv'))
        if(self%o_in%isthere('fraca')) call self%softoris%set_all2single('fraca',self%o_in%get('fraca'))
        if(self%o_in%isthere('lp'))    call self%softoris%set_all2single('lp',self%o_in%get('lp'))
        call self%softoris%set_all2single('w', self%o_in%get('w'))
        call self%softoris%set_all2single('specscore', self%specscore)
        ! best orientation
        self%o_out = self%softoris%get_ori(self%npeaks)
        call self%o_out%set('corr', wcorr)  
        ! angular distances & deviation
        euldist   = rad2deg( self%o_in.euldist.self%o_out )
        dist_inpl = rad2deg( self%o_in.inplrotdist.self%o_out )
        ang_sdev  = self%softoris%ang_sdev(self%refine, self%nstates, self%npeaks)
        call self%o_out%set('dist', euldist)
        call self%o_out%set('dist_inpl', dist_inpl)
        call self%o_out%set('sdev', ang_sdev)
        ! overlap between distributions
        euldist_thresh = max(0.1, self%angthresh/10.)
        roind    = self%pftcc_ptr%get_roind(360.-self%o_out%e3get())
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_state = 0.
        mi_joint = 0.
        if( euldist < 0.1 )then
            mi_proj  = mi_proj + 1.
            mi_joint = mi_joint + 1.
        endif
        if(self%prev_roind == roind)then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif
        if(self%nstates > 1)then
            state      = nint(self%o_out%get('state'))
            prev_state = nint(self%o_in%get('state'))
            if(prev_state == state)then
                mi_state = mi_state + 1.
                mi_joint = mi_joint + 1.
            endif
            mi_joint = mi_joint/3.
        else
            mi_joint = mi_joint/2.
        endif
        call self%o_out%set('proj', 0.)
        call self%o_out%set('mi_proj',  mi_proj)
        call self%o_out%set('mi_inpl',  mi_inpl)
        call self%o_out%set('mi_state', mi_state)
        call self%o_out%set('mi_joint', mi_joint)       
        ! frac
        frac = 100. * (real(self%neval)/real(self%nrefs))
        call self%o_out%set('frac', frac)
        ! done
        if(debug)write(*,*)'simple_cont3D_srch::prep_softoris done'
    end subroutine prep_softoris

    ! GETTERS

    !>  \brief  returns the solution set of orientations
    function get_softoris( self )result( os )
        class(cont3D_srch), intent(inout) :: self
        type(oris) :: os
        if(.not.self%exists)stop 'search has not been performed; cont3D_srch::get_softoris'
        os = self%softoris
    end function get_softoris

    !>  \brief  returns the best solution orientation
    function get_best_ori( self )result( o )
        class(cont3D_srch), intent(inout) :: self
        type(ori) :: o
        if(.not.self%exists)stop 'search has not been performed; cont3D_srch::get_softoris'
        o = self%o_out
    end function get_best_ori

    !>  \brief  whether object has been initialized
    function does_exist( self )result( l )
        class(cont3D_srch), intent(inout) :: self
        logical :: l
        l = self%exists
    end function does_exist

    ! DESTRUCTOR

    !>  \brief  is the destructor
    subroutine kill( self )
        class(cont3D_srch), intent(inout) :: self
        call self%shsrch_obj%kill
        self%pftcc_ptr => null()
        call self%shiftedoris%kill
        call self%softoris%kill
        call self%reforis%kill
        call self%o_in%kill
        call self%o_out%kill
        self%iptcl  = 0
        self%exists = .false.
    end subroutine kill

end module simple_cont3D_srch
