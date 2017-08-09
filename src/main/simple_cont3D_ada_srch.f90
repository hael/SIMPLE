! continuous adaptive 3D orientation search
module simple_cont3D_ada_srch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_projector,        only: projector
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_sym,              only: sym
use simple_pftcc_inplsrch,   only: pftcc_inplsrch
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_math              ! use all in there
implicit none

public :: cont3D_ada_srch
private
! #include "simple_local_flags.inc"
real,    parameter :: FACTWEIGHTS_THRESH = 0.001    !< threshold for factorial weights
real,    parameter :: E3HALFWINSZ = 90.             !< in-plane angle half window size
integer, parameter :: NREFS       = 80
logical, parameter :: debug = .false.

type cont3D_ada_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr   => null()  !< polar fourier correlation calculator
    class(projector),        pointer :: vols_ptr(:) => null()  !< volumes for projection
    type(pftcc_inplsrch)             :: inplsrch_obj           !< in-plane search object
    type(pftcc_shsrch)               :: shiftsrch_obj          !< in-plane search object
    type(sym)                        :: se                     !< symmery object
    type(ori)                        :: o_in                   !< input orientation
    type(ori)                        :: o_out                  !< best orientation found
    type(oris)                       :: o_peaks                !< peaks
    type(oris)                       :: o_srch                 !< all search orientations
    type(ori)                        :: o_best                 !< current best orientation
    logical,             allocatable :: state_exists(:)        !< indicates whether each state is populated
    real                             :: lims(5,2)     = 0.     !< shift search limits
    real                             :: best_corr     = -1.    !< current best correlation
    real                             :: prev_shift(2) = 0.     !< previous shift
    real                             :: angthresh     = 0.     !< angular threshold
    real                             :: specscore     = 0.     !< previous spectral score
    integer                          :: nrefs         = 0      !< total number of references
    integer                          :: nrots         = 0      !< number of in-plane rotations
    integer                          :: nbetter       = 0      !< number of time an improving reference is found
    integer                          :: npeaks        = 0      !< number of peaks
    integer                          :: npeaks_inpl   = 0      !< number of peaks
    integer                          :: iptcl         = 0      !< orientation general index
    integer                          :: ref           = 0      !< previous ref
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
    procedure, private :: prep_peaks
    procedure, private :: stochastic_weights
    ! SEARCH ROUTINES
    procedure          :: exec_srch
    procedure, private :: do_euler_srch
    procedure, private :: do_inpl_srch
    procedure, private :: gen_ori
    ! GETTERS/SETTERS
    procedure          :: get_best_ori
    procedure          :: get_peaks
    procedure          :: does_exist
    ! DESTRUCTOR
    procedure          :: kill
end type cont3D_ada_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, p, pftcc, vols, se )
        class(cont3D_ada_srch),           intent(inout) :: self
        class(params),                   intent(in)    :: p        !< parameters
        class(polarft_corrcalc), target, intent(in)    :: pftcc    !< corrcalc object
        class(projector),        target, intent(in)    :: vols(:)  !< volume for references
        class(sym),                      intent(in)    :: se       !< symmetry object
        real :: trs
        call self%kill
        ! set constants
        self%pftcc_ptr => pftcc
        self%vols_ptr  => vols
        self%se        = se
        self%nstates   = p%nstates
        if( size(vols).ne.self%nstates )stop 'Inconsistent number of volumes; cont3D_ada_srch%new'
        self%angthresh   = p%athres
        self%shbarr      = p%shbarrier
        self%npeaks      = p%npeaks
        self%npeaks_inpl = self%npeaks * 3
        self%refine      = p%refine
        self%nstates     = p%nstates
        self%nrefs       = NREFS !!! 1 state assumed here
        self%nrots       = self%pftcc_ptr%get_nrots()
        self%state       = 1     !!! 1 state assumed here
        self%lims(:3,:)  = self%se%srchrange()
        trs = p%trs
        self%lims(4,:)  = [-trs, trs]
        self%lims(5,:)  = self%lims(4,:)
        if( pftcc%get_nrefs() .ne. self%nrefs )&
        &stop 'Non-congruent number of references in pftcc; cont3D_ada_srch%new'
        ! done
        self%exists = .true.
        if( debug )write(*,'(A)') '>>> cont3D_ada_srch%CONSTRUCTED NEW SIMPLE_cont3D_ada_srch OBJECT'
    end subroutine new

    ! PREP ROUTINES

    !>  \brief  is the master search routine
    subroutine prep_srch(self, a, iptcl)
        class(cont3D_ada_srch), intent(inout) :: self
        class(oris),         intent(inout) :: a     !< search orientations
        integer,             intent(in)    :: iptcl !< input particle
        type(ori)         :: o
        real, allocatable :: frc(:)
        integer           :: inpl_ind
        if(iptcl == 0)stop 'ptcl index mismatch; cont3D_ada_srch::prep_srch'
        ! init
        self%iptcl      = iptcl
        self%o_in       = a%get_ori(self%iptcl)
        self%o_best     = self%o_in
        call self%o_best%set_shift([0.,0.])
        self%prev_state = nint(self%o_in%get('state'))
        self%prev_shift = self%o_in%get_shift()
        call self%o_srch%new( self%nrefs )
        call self%o_peaks%new( self%npeaks )
        ! state init
        allocate(self%state_exists(self%nstates))
        self%state_exists = a%get_state_exist(self%nstates)
        if( .not.self%state_exists(self%prev_state) )stop 'state is empty; cont3D_ada_srch::prep_srch'
        ! current reference put first
        o = self%o_in
        call o%e3set(0.)
        call self%vols_ptr(self%state)%fproject_polar(1, o, self%pftcc_ptr, serial=.true.)
        inpl_ind       = self%pftcc_ptr%get_roind(360. - self%o_best%e3get())
        self%best_corr = self%pftcc_ptr%corr(1, self%iptcl, inpl_ind)
        ! specscore
        frc = self%pftcc_ptr%genfrc(1, self%iptcl, inpl_ind)
        self%specscore = max(0., median_nocopy(frc))
        deallocate(frc)
        ! in-plane search object
        call self%inplsrch_obj%new(self%pftcc_ptr, self%lims(4:5,:), shbarrier=self%shbarr, nrestarts=4)
        if( debug ) write(*,'(A)') '>>> cont3D_ada_srch::END OF PREP_SRCH'
    end subroutine prep_srch

    !>  \brief  updates solutions orientations
    subroutine prep_peaks( self )
        class(cont3D_ada_srch), intent(inout) :: self
        real, allocatable :: corrs(:)
        type(ori) :: o
        real      :: ang_sdev, dist_inpl, frac, wcorr, euldist, mi_proj,mi_state
        integer   :: i, state, prev_state, cnt, ref_inds(self%nrefs), iref, best_loc(1)
        ! sort oris
        ref_inds = (/ (iref, iref=1,self%nrefs) /)
        corrs    = self%o_srch%get_all('corr')
        call hpsort(self%nrefs, corrs, ref_inds)
        cnt = 0
        do i = self%nrefs, self%nrefs-self%npeaks+1, -1
            cnt  = cnt + 1
            iref = ref_inds(i)
            o    = self%o_srch%get_ori(iref)
            call self%o_peaks%set_ori(cnt, o)
        enddo
        deallocate(corrs)
        ! best reference index
        corrs    = self%o_peaks%get_all('corr')
        best_loc = maxloc(corrs)
        deallocate(corrs)
        ! weights and weighted correlation
        call self%stochastic_weights( wcorr )
        ! variables inherited from input ori: CTF, w, lp
        if(self%o_in%isthere('dfx'))   call self%o_peaks%set_all2single('dfx',self%o_in%get('dfx'))
        if(self%o_in%isthere('dfy'))   call self%o_peaks%set_all2single('dfy',self%o_in%get('dfy'))
        if(self%o_in%isthere('angast'))call self%o_peaks%set_all2single('angast',self%o_in%get('angast'))
        if(self%o_in%isthere('cs'))    call self%o_peaks%set_all2single('cs',self%o_in%get('cs'))
        if(self%o_in%isthere('kv'))    call self%o_peaks%set_all2single('kv',self%o_in%get('kv'))
        if(self%o_in%isthere('fraca')) call self%o_peaks%set_all2single('fraca',self%o_in%get('fraca'))
        if(self%o_in%isthere('lp'))    call self%o_peaks%set_all2single('lp',self%o_in%get('lp'))
        call self%o_peaks%set_all2single('w', self%o_in%get('w'))
        call self%o_peaks%set_all2single('specscore', self%specscore)
        ! best orientation
        self%o_out = self%o_peaks%get_ori(best_loc(1))
        call self%o_out%set('corr', wcorr)  
        ! angular distances & deviation
        o = self%o_out
        call self%se%sym_euldist( self%o_in, o, euldist )
        dist_inpl = rad2deg( self%o_in.inplrotdist.o )
        ang_sdev  = self%o_peaks%ang_sdev(self%refine, self%nstates, self%npeaks)
        call self%o_out%set('dist',      euldist)
        call self%o_out%set('dist_inpl', dist_inpl)
        call self%o_out%set('sdev',      ang_sdev)
        ! overlap between distributions
        mi_proj  = 0.
        mi_state = 0.
        if( euldist < 0.5 )mi_proj  = mi_proj + 1.
        if(self%nstates > 1)then
            state      = nint(self%o_out%get('state'))
            prev_state = nint(self%o_in%get('state'))
            if(prev_state == state)mi_state = mi_state + 1.
        else
            mi_state = 1.
        endif
        call self%o_out%set('mi_proj',  mi_proj)
        call self%o_out%set('mi_state', mi_state)
        if(self%o_out%isthere('proj'))call self%o_out%set('proj', 0.)
        if(self%o_in%isthere('mi_inpl')) call self%o_out%set('mi_inpl', 1.)
        if(self%o_in%isthere('mi_joint'))call self%o_out%set('mi_joint', 1.)
        ! fractional search space
        frac = 100. * real(self%nrefs - self%nbetter) / real(self%nrefs)
        call self%o_out%set('frac', frac)
        ! done
        if(debug)write(*,*)'simple_cont3D_ada_srch::prep_o_peaks done'
    end subroutine prep_peaks

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine exec_srch( self, a, iptcl )
        class(cont3D_ada_srch), intent(inout) :: self
        class(oris),        intent(inout) :: a      !< search orientations
        integer,            intent(in)    :: iptcl  !< input particle
        if(nint(a%get(iptcl,'state')) > 0)then
            call self%prep_srch(a, iptcl)
            call self%do_euler_srch
            call self%do_inpl_srch
            call self%prep_peaks
            call a%set_ori(self%iptcl, self%o_out)
        else
            call a%reject(iptcl)
        endif
        if( debug ) write(*,'(A)') '>>> cont3D_ada_srch::END OF SRCH'
    end subroutine exec_srch

    !>  \brief  is for generating a random orientation in the vicinity of
    !>          the current best
    subroutine gen_ori( self, o )
        use simple_rnd, only: gasdev
        class(cont3D_ada_srch), intent(inout) :: self
        class(ori),             intent(out)   :: o !< rand orientation
        type(ori) :: o_transform
        real      :: val
        call o_transform%new
        call o_transform%rnd_euler
        val = gasdev(0., self%angthresh)
        call enforce_cyclic_limit(val, 180.)
        call o_transform%e2set(val)
        call o_transform%e3set(0.)
        o = self%o_best.compose.o_transform
        call self%se%rot_to_asym( o )
        call o%e3set(0.)
    end subroutine gen_ori

    !>  \brief  performs euler angles search
    subroutine do_euler_srch( self )
        class(cont3D_ada_srch), intent(inout) :: self
        type(ori)            :: o
        integer, allocatable :: roind_vec(:)    ! slice of in-plane angles
        real                 :: corrs(self%nrots), inpl_corr
        integer              :: iref, loc(1), inpl_ind
        self%nbetter = 0
        do iref = 1, self%nrefs
            ! reference extraction
            if( iref.eq.1 )then
                ! previous best orientation: projection done in self%prep_srch
                roind_vec = self%pftcc_ptr%get_win_roind(360.-self%o_best%e3get(), E3HALFWINSZ)
                o = self%o_best
            else
                call self%gen_ori( o )
                call self%vols_ptr(self%state)%fproject_polar(iref, o, self%pftcc_ptr, serial=.true.)
            endif
            ! in-plane search
            corrs     = self%pftcc_ptr%gencorrs(iref, self%iptcl, roind_vec=roind_vec)
            loc       = maxloc(corrs)
            inpl_ind  = loc(1)
            inpl_corr = corrs(inpl_ind)
            call o%e3set(360. - self%pftcc_ptr%get_rot(inpl_ind))
            call o%set('corr', inpl_corr)
            if(inpl_corr > self%best_corr)then
                ! updates new centre
                self%best_corr = inpl_corr
                self%nbetter   = self%nbetter + 1
                self%o_best    = o
                deallocate(roind_vec)
                roind_vec = self%pftcc_ptr%get_win_roind(360.-self%o_best%e3get(), E3HALFWINSZ)
            endif
            ! stash
            call self%o_srch%set_ori(iref, o)
        enddo
        ! done
        deallocate(roind_vec)
        if(debug)write(*,*)'simple_cont3D_srch::do_refs_srch done'
    end subroutine do_euler_srch

    !>  \brief  performs the in-plane search
    subroutine do_inpl_srch( self )
        class(cont3D_ada_srch), intent(inout) :: self
        real, allocatable :: corrs(:)
        type(ori)         :: o
        real, allocatable :: crxy(:)
        real              :: corr
        integer           :: i, iref, inpl_ind, ref_inds(self%nrefs)
        ! sort searched orientations
        ref_inds = (/(iref,iref=1,self%nrefs)/)
        corrs = self%o_srch%get_all('corr')
        call hpsort(self%nrefs, corrs, ref_inds)
        deallocate(corrs)
        ! in-plane search
        do i = self%nrefs-self%npeaks+1, self%nrefs
            iref     = ref_inds(i)
            o        = self%o_srch%get_ori(iref)
            corr     = o%get('corr')
            inpl_ind = self%pftcc_ptr%get_roind(360.-o%e3get())
            ! in-plane search
            call self%inplsrch_obj%set_indices(iref, self%iptcl)
            crxy = self%inplsrch_obj%minimize(irot=inpl_ind)
            if(crxy(1) >= corr)then
                call o%set('corr', crxy(1))
                call o%e3set(360. - crxy(2))
                ! shift addition
                call o%set_shift(crxy(3:4) + self%prev_shift)
            else
                call o%set_shift(self%prev_shift)
            endif
            deallocate(crxy)
            ! stores solution
            call self%o_srch%set_ori(iref, o)
        enddo
        if(debug)write(*,*)'simple_cont3d_srch::do_inpl_srch done'
    end subroutine do_inpl_srch

    !>  \brief  determines and updates stochastic weights
    subroutine stochastic_weights( self, wcorr )
        class(cont3D_ada_srch), intent(inout) :: self
        real,                   intent(out)   :: wcorr   !<  stochastic weights
        type(ori)         :: o
        real, allocatable :: frc(:)
        real              :: ws(self%npeaks), corrs(self%npeaks), logws(self%npeaks), frcmed
        integer           :: ipeak, iref, iroind, order(self%npeaks)
        logical           :: included(self%npeaks)
        if( self%npeaks == 1 )then
            call self%o_peaks%set(1, 'ow', 1.)
            wcorr = self%o_peaks%get(1, 'corr')
            return
        endif
        corrs = self%o_peaks%get_all('corr')
        ! calculate L1 norms
        ! do ipeak = 1, self%npeaks
        !     o      = self%o_peaks%get_ori(ipeak)
        !     iroind = self%pftcc_ptr%get_roind(360.-o%e3get())
        !     !!!!!!!!!!!!!!! note state default 1 here
        !     iref   = self%o_srch%find_closest_proj(o, 1)
        !     ! frc    = self%pftcc_ptr%genfrc(iref, self%iptcl, iroind)
        !     ! frcmed = max(0., median_nocopy(frc))
        !     ! ws(ipeak) = exp(-(1.-frcmed))
        !     ! deallocate(frc)
        ! end do
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
    end subroutine stochastic_weights

    ! GETTERS

    !>  \brief  returns the solution set of orientations
    function get_peaks( self )result( os )
        class(cont3D_ada_srch), intent(inout) :: self
        type(oris) :: os
        if(.not.self%exists)stop 'search has not been performed; cont3D_ada_srch::get_o_peaks'
        os = self%o_peaks
    end function get_peaks

    !>  \brief  returns the best solution orientation
    function get_best_ori( self )result( o )
        class(cont3D_ada_srch), intent(inout) :: self
        type(ori) :: o
        if(.not.self%exists)stop 'search has not been performed; cont3D_ada_srch::get_best_ori'
        o = self%o_out
    end function get_best_ori

    ! GETTERS/SETTERS

    !>  \brief  whether object has been initialized
    logical function does_exist( self )
        class(cont3D_ada_srch), intent(inout) :: self
        does_exist = self%exists
    end function does_exist

    ! DESTRUCTOR

    !>  \brief  is the destructor
    subroutine kill( self )
        class(cont3D_ada_srch), intent(inout) :: self
        self%pftcc_ptr => null()
        self%vols_ptr  => null()
        call self%o_in%kill
        call self%o_out%kill
        call self%o_best%kill
        call self%o_peaks%kill
        call self%o_srch%kill
        if( allocated(self%state_exists) )deallocate(self%state_exists)
        self%exists = .false.
    end subroutine kill

end module simple_cont3D_ada_srch
