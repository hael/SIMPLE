!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple search class: 3D continous algorithm
module simple_cont3D_srch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_pftcc_inplsrch,   only: pftcc_inplsrch
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_math              ! use all in there
implicit none

public :: cont3D_srch
private

real,    parameter :: E3HALFWINSZ = 90. !< in-plane angle half window size
#include "simple_local_flags.inc"

type cont3D_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr => null() !< polar fourier correlation calculator
    type(pftcc_shsrch)               :: shsrch_obj          !< shift search object
    type(pftcc_inplsrch)             :: inplsrch_obj        !< in-plane search object
    type(oris)                       :: reforis             !< per ptcl search space and result of euler angles search
    type(oris)                       :: softoris            !< references returned
    type(oris)                       :: shiftedoris         !< references whose shifts are searched
    type(ori)                        :: o_in                !< input orientation
    type(ori)                        :: o_out               !< best orientation found
    logical,             allocatable :: state_exists(:)     !< indicates whether each state is populated
    real                             :: lims(2,2)  = 0.     !< shift search limits
    real                             :: prev_corr  = -1.    !< previous correlation
    real                             :: angthresh  = 0.     !< angular threshold
    real                             :: specscore  = 0.     !< previous spectral score
    integer                          :: iptcl      = 0      !< orientation general index
    integer                          :: prev_ref   = 0      !< previous reference
    integer                          :: prev_roind = 0      !< previous in-plane rotational index
    integer                          :: prev_state = 0      !< previous state
    integer                          :: nstates    = 0      !< number of states
    integer                          :: npeaks     = 0      !< number of references returned
    integer                          :: nrefs      = 0      !< number of references in search space
    integer                          :: nrots      = 0      !< number of in-plane rotations
    character(len=STDLEN)            :: shbarr = ''         !< shift barrier flag
    character(len=STDLEN)            :: refine = ''
    logical                          :: exists = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PREP ROUTINES
    procedure, private :: prep_srch
    procedure, private :: prep_softoris
    procedure, private :: stochastic_weights_L2norm
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
        class(cont3D_srch),              intent(inout) :: self  !< instance
        class(params),                   intent(in)    :: p     !< parameters
        class(oris),                     intent(in)    :: e     !< references
        class(polarft_corrcalc), target, intent(in)    :: pftcc !< correlation calculator obj
        call self%kill
        ! set constants
        self%pftcc_ptr => pftcc
        self%reforis   = e
        self%npeaks    = p%npeaks
        self%lims(:,1) = -p%trs
        self%lims(:,2) =  p%trs
        self%angthresh = p%athres
        self%nstates   = p%nstates
        self%shbarr    = p%shbarrier
        self%nrefs     = self%pftcc_ptr%get_nrefs()
        self%nrots     = self%pftcc_ptr%get_nrots()
        if(self%reforis%get_noris().ne.self%nrefs)&
            &stop 'Inconsistent number of references & orientations'
        ! done
        self%exists = .true.
        DebugPrint  '>>> cont3D_srch::CONSTRUCTED NEW SIMPLE_cont3D_srch OBJECT'
    end subroutine new

    ! PREP ROUTINES

    !>  \brief  is the master search routine
    subroutine prep_srch(self, a, iptcl)
        class(cont3D_srch),  intent(inout) :: self
        class(oris),         intent(inout) :: a     !< oris object        
        integer,             intent(in)    :: iptcl !< index to ori in a  
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
        deallocate(frc)
        ! shift search object
        call self%shsrch_obj%new(self%pftcc_ptr, self%lims,&
        &shbarrier=self%shbarr, nrestarts=5)
        call self%inplsrch_obj%new(self%pftcc_ptr, self%lims,&
        &shbarrier=self%shbarr, nrestarts=5)
        DebugPrint '>>> cont3D_srch::END OF PREP_SRCH'
    end subroutine prep_srch

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine exec_srch( self, a, iptcl )
        class(cont3D_srch), intent(inout) :: self
        class(oris),        intent(inout) :: a       !< oris object - references
        integer,            intent(in)    :: iptcl   !< index to ori in a
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
        DebugPrint  '>>> cont3D_srch::END OF SRCH'
    end subroutine exec_srch

    !>  \brief  performs euler angles search
    subroutine do_euler_srch( self )
        class(cont3D_srch), intent(inout) :: self
        !integer, allocatable :: roind_vec(:)    ! slice of in-plane angles
        real                 :: inpl_corr
        integer              :: iref
        ! roind_vec    = self%pftcc_ptr%get_win_roind(360.-self%o_in%e3get(), E3HALFWINSZ)
        ! search
        ! the input search space in stochastic, so no need for a randomized search order
        do iref = 1,self%nrefs
            call greedy_inpl_srch(iref, inpl_corr)
        enddo         
        ! deallocate(roind_vec)
        if(debug)write(*,*)'simple_cont3D_srch::do_refs_srch done'

        contains
            !> greedy inplane search
            subroutine greedy_inpl_srch(iref_here, corr_here)
                integer, intent(in)    :: iref_here
                real,    intent(inout) :: corr_here
                real    :: corrs(self%nrots), e3
                integer :: loc(1), inpl_ind
                ! corrs     = self%pftcc_ptr%gencorrs(iref_here, self%iptcl, roind_vec=roind_vec)
                corrs     = self%pftcc_ptr%gencorrs(iref_here, self%iptcl)
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
        class(cont3D_srch), intent(inout) :: self
        real, allocatable :: corrs(:)
        type(ori)         :: o
        real, allocatable :: cxy(:), crxy(:)
        real              :: prev_shift_vec(2), corr
        integer           :: ref_inds(self%nrefs), i, iref, inpl_ind, cnt
        logical :: greedy_inpl = .true.
        call self%shiftedoris%new(self%npeaks)
        prev_shift_vec = self%o_in%get_shift()
        ! SORT ORIENTATIONS
        ref_inds = (/ (i, i=1, self%nrefs) /)
        corrs    = self%reforis%get_all('corr')
        call hpsort(self%nrefs, corrs, ref_inds)
        deallocate(corrs)
        ! IN-PLANE/SHIFT SEARCH
        cnt = 0
        do i = self%nrefs-self%npeaks+1,self%nrefs
            cnt      = cnt+1
            iref     = ref_inds(i)
            o        = self%reforis%get_ori(iref)
            corr     = o%get('corr')
            inpl_ind = self%pftcc_ptr%get_roind(360.-o%e3get())
            if( greedy_inpl )then
                ! in-plane search
                call self%inplsrch_obj%set_indices(iref, self%iptcl)
                crxy = self%inplsrch_obj%minimize(irot=inpl_ind)
                if(crxy(1) >= corr)then
                    call o%set('corr', crxy(1))
                    call o%e3set(360. - crxy(2))
                    call o%set_shift(crxy(3:4) + prev_shift_vec)
                else
                    call o%set_shift(prev_shift_vec)
                endif
                deallocate(crxy)
            else
                ! shift search
                call self%shsrch_obj%set_indices(iref, self%iptcl, inpl_ind)
                cxy = self%shsrch_obj%minimize()
                if(cxy(1) >= corr)then
                    call o%set('corr', cxy(1))
                    call o%set_shift(cxy(2:) + prev_shift_vec)
                else
                    call o%set_shift(prev_shift_vec)
                endif
                deallocate(cxy)
            endif
            ! stores solution
            call self%shiftedoris%set_ori(cnt, o)
        enddo
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
        call self%stochastic_weights_L2norm(wcorr)
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
        roind    = self%pftcc_ptr%get_roind(360.-self%o_out%e3get())
        mi_proj  = 0.
        mi_state = 0.
        if( euldist < 0.8 )mi_proj  = mi_proj + 1.
        if(self%nstates > 1)then
            state      = nint(self%o_out%get('state'))
            prev_state = nint(self%o_in%get('state'))
            if(prev_state == state)mi_state = mi_state + 1.
        endif
        call self%o_out%set('mi_proj',  mi_proj)
        call self%o_out%set('mi_state', mi_state)
        if(self%o_out%isthere('proj'))call self%o_out%set('proj', 0.)
        if(self%o_in%isthere('mi_inpl')) call self%o_out%set('mi_inpl', 1.)
        if(self%o_in%isthere('mi_joint'))call self%o_out%set('mi_joint', 1.)
        ! fractional search space
        call self%o_out%set('frac', 100.)
        ! done
        if(debug)write(*,*)'simple_cont3D_srch::prep_softoris done'
    end subroutine prep_softoris

    !>  \brief  determines and updates stochastic weights
    subroutine stochastic_weights_L2norm( self, wcorr )
        class(cont3D_srch),      intent(inout) :: self
        real,                    intent(out)   :: wcorr
        real              :: ws(self%npeaks), dists(self%npeaks)
        real, allocatable :: corrs(:)
        type(ori)         :: o
        integer           :: roind, ref, ipeak
        if( self%npeaks == 1 )then
            call self%softoris%set(1, 'ow', 1.)
            wcorr = self%softoris%get(1, 'corr')
            return
        endif
        ! get correlations
        corrs = self%softoris%get_all('corr')
        ! calculate L1 norms
        do ipeak = 1, self%npeaks
            o     = self%softoris%get_ori(ipeak)
            roind = self%pftcc_ptr%get_roind(360.-o%e3get())
            ref   = self%reforis%find_closest_proj(o, 1)
            dists(ipeak) = self%pftcc_ptr%euclid(ref, self%iptcl, roind)
        end do
        ! calculate normalised weights and weighted corr
        ! ws    = exp(-dists)
        ! ws    = ws/sum(ws)
        ! wcorr = sum(ws*corrs)
        ! calculate weights and weighted corr
        ws    = exp(-dists)
        wcorr = sum(ws*corrs) / sum(ws)
        ! update npeaks individual weights
        call self%softoris%set_all('ow', ws)
        ! cleanup
        deallocate(corrs)
    end subroutine stochastic_weights_L2norm

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
        if(allocated(self%state_exists))deallocate(self%state_exists)
        self%iptcl  = 0
        self%exists = .false.
    end subroutine kill

end module simple_cont3D_srch
