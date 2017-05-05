module simple_pcont3D_srch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_math              ! use all in there
implicit none

public :: pcont3D_srch
private

real,    parameter :: E3HALFWINSZ = 90. !< in-plane angle half window size
logical, parameter :: debug = .false.

type pcont3D_srch
    private
    class(polarft_corrcalc), pointer :: ppftcc => null()  !< polar fourier correlation calculator
    type(pftcc_shsrch)               :: shsrch_obj        !< shift search object
    type(oris) :: reforis             !< per ptcl search space and result of euler angles search
    type(oris) :: softoris            !< references returned
    type(oris) :: shiftedoris         !< references whose shifts are searched
    type(ori)  :: orientation_in      !< input orientation
    type(ori)  :: orientation_out     !< best orientation found

    real       :: lims(2,2)  = 0.     !< shift search limits
    real       :: prev_corr  = -1.    !< previous correlation
    real       :: specscore  = 0.     !< previous spectral score
    integer    :: iptcl      = 0      !< orientation general index
    integer    :: prev_ref   = 0      !< previous reference
    integer    :: prev_roind = 0      !< previous in-plane rotational index
    integer    :: prev_state = 0      !< previous state
    integer    :: nstates    = 0      !< number of states
    integer    :: nbetter    = 0      !< number of improving references found
    integer    :: neval      = 0      !< number of references evaluated
    integer    :: npeaks     = 0      !< number of references returned
    integer    :: nrefs      = 0      !< number of references in search space
    integer    :: nrots      = 0      !< number of in-plane rotations
    logical    :: exists = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PREP ROUTINES
    procedure, private :: prep_srch
    procedure, private :: prep_softoris
    ! SEARCH ROUTINES
    procedure          :: do_srch
    procedure, private :: do_refs_srch
    procedure, private :: do_shift_srch
    procedure, private :: ang_sdev
    ! GETTERS
    procedure          :: get_best_ori
    procedure          :: get_softoris
    procedure          :: does_exist
    ! DESTRUCTOR
    procedure          :: kill
end type pcont3D_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, p, a, e, pftcc, iptcl )
        class(pcont3D_srch),             intent(inout) :: self  !< instance
        class(params),                   intent(in)    :: p     !< parameters
        class(oris),                     intent(inout) :: a     !< ptcls orientations 
        class(oris),                     intent(in)    :: e     !< references
        class(polarft_corrcalc), target, intent(in)    :: pftcc
        integer,                         intent(in)    :: iptcl !< general particle index
        call self%kill
        ! particle index
        self%iptcl = iptcl
        ! input orientation
        self%orientation_in = a%get_ori(self%iptcl)
        self%prev_state     = nint(self%orientation_in%get('state'))
        if(self%prev_state == 0)return
        ! set constants
        self%ppftcc    => pftcc
        self%reforis   = e
        self%npeaks    = p%npeaks
        self%lims(:,1) = -p%trs
        self%lims(:,2) =  p%trs
        self%nstates   = p%nstates
        self%nrefs     = self%ppftcc%get_nrefs()
        self%nrots     = self%ppftcc%get_nrots()
        if(self%reforis%get_noris().ne.self%nrefs)&
            &stop 'Inconsistent number of references & orientations'
        ! done
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> PCONT3D_SRCH::CONSTRUCTED NEW SIMPLE_PCONT3D_SRCH OBJECT'
    end subroutine new

    ! PREP ROUTINES

    !>  \brief  is the master search routine
    subroutine prep_srch( self )
        class(pcont3D_srch), intent(inout) :: self
        real, allocatable :: frc(:)
        ! INIT
        call self%shsrch_obj%new(self%ppftcc, self%lims)
        self%prev_roind = self%ppftcc%get_roind(360.-self%orientation_in%e3get())
        self%prev_ref   = self%reforis%find_closest_proj(self%orientation_in) ! state not taken into account
        self%prev_corr  = self%ppftcc%corr(self%prev_ref, self%iptcl, self%prev_roind)
        frc             = self%ppftcc%genfrc(self%prev_ref, self%iptcl, self%prev_roind)
        self%specscore  = max(0., median_nocopy(frc))
        deallocate(frc)
        if( debug ) write(*,'(A)') '>>> PCONT3D_SRCH::END OF PREP_SRCH'
    end subroutine prep_srch

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine do_srch( self, a )
        class(pcont3D_srch), intent(inout) :: self
        class(oris),         intent(inout) :: a
        if(self%exists)then
            ! INIT
            call self%prep_srch
            ! EULER ANGLES SEARCH
            call self%do_refs_srch
            ! SHIFT SEARCH
            call self%do_shift_srch
            ! outcome prep
            call self%prep_softoris
        else
            call self%orientation_out%reject
        endif
        call a%set_ori(self%iptcl, self%orientation_out)
        if( debug ) write(*,'(A)') '>>> PCONT3D_SRCH::END OF SRCH'
    end subroutine do_srch

    !>  \brief  performs euler angles search
    subroutine do_refs_srch( self )
        class(pcont3D_srch), intent(inout) :: self
        integer, allocatable :: roind_vec(:)    ! slice of in-plane angles
        real                 :: inpl_corr
        integer              :: iref
        ! init
        self%nbetter = 0
        self%neval   = 0
        roind_vec = self%ppftcc%get_win_roind(360.-self%orientation_in%e3get(), E3HALFWINSZ)
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
        if(debug)write(*,*)'simple_pcont3d_srch::do_refs_srch done'

        contains

            ! Unused at the moment, for testing purpose
            subroutine shc_inpl_srch(iref_here, corr_here)
                use simple_rnd, only: shcloc
                integer, intent(in)    :: iref_here
                real,    intent(inout) :: corr_here
                real    :: corrs(self%nrots), e3
                integer :: inpl_ind
                corrs    = self%ppftcc%gencorrs(iref_here, self%iptcl, roind_vec=roind_vec)
                inpl_ind = shcloc(self%nrots, corrs, self%prev_corr)
                if(inpl_ind > 0)then
                    corr_here = corrs(inpl_ind)
                    e3 = 360. - self%ppftcc%get_rot(inpl_ind)
                    call self%reforis%e3set(iref_here, e3)
                else
                    corr_here = 0.
                endif
                call self%reforis%set(iref_here, 'corr', corr_here)
            end subroutine shc_inpl_srch

            subroutine greedy_inpl_srch(iref_here, corr_here)
                integer, intent(in)    :: iref_here
                real,    intent(inout) :: corr_here
                real    :: corrs(self%nrots), e3
                integer :: loc(1), inpl_ind
                corrs     = self%ppftcc%gencorrs(iref_here, self%iptcl, roind_vec=roind_vec)
                loc       = maxloc(corrs)
                inpl_ind  = loc(1)
                corr_here = corrs(inpl_ind)
                e3 = 360. - self%ppftcc%get_rot(inpl_ind)
                call self%reforis%e3set(iref_here, e3)
                call self%reforis%set(iref_here, 'corr', corr_here)
            end subroutine greedy_inpl_srch
    end subroutine do_refs_srch

    !>  \brief  performs the shift search
    subroutine do_shift_srch( self )
        use simple_pftcc_shsrch      ! use all in there
        class(pcont3D_srch), intent(inout) :: self
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
        prev_shift_vec = self%orientation_in%get_shift()
        cnt = 0
        do i = self%nrefs-self%npeaks+1,self%nrefs
            cnt      = cnt+1
            iref     = ref_inds(i)
            o        = self%reforis%get_ori(iref)
            corr     = o%get('corr')
            inpl_ind = self%ppftcc%get_roind(360.-o%e3get())
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
        if(debug)write(*,*)'simple_pcont3d_src::do_shift_srch done'
    end subroutine do_shift_srch

    !>  \brief  updates solutions orientations
    subroutine prep_softoris( self )
        class(pcont3D_srch), intent(inout) :: self
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        type(ori) :: o
        real      :: ws(self%npeaks), u(2), mat(2,2), x1(2), x2(2)
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
            ws    = 0.
            corrs = self%softoris%get_all('corr')
            where(corrs > TINY)ws = exp(corrs)
            ws    = ws/sum(ws)
            wcorr = sum(ws*corrs) 
            call self%softoris%set_all('ow', ws)
            deallocate( corrs )
        else
            call self%softoris%set(1, 'ow', 1.)
            wcorr = self%softoris%get(1, 'corr')
        endif
        ! variables inherited from input ori: CTF, w, lp
        o = self%orientation_in
        if(o%isthere('dfx'))   call self%softoris%set_all2single('dfx',o%get('dfx'))
        if(o%isthere('dfy'))   call self%softoris%set_all2single('dfy',o%get('dfy'))
        if(o%isthere('angast'))call self%softoris%set_all2single('angast',o%get('angast'))
        if(o%isthere('cs'))    call self%softoris%set_all2single('cs',o%get('cs'))
        if(o%isthere('kv'))    call self%softoris%set_all2single('kv',o%get('kv'))
        if(o%isthere('fraca')) call self%softoris%set_all2single('fraca',o%get('fraca'))
        if(o%isthere('lp'))    call self%softoris%set_all2single('lp',o%get('lp'))
        call self%softoris%set_all2single('w', o%get('w'))
        call self%softoris%set_all2single('specscore', self%specscore)
        call o%kill
        ! best orientation
        self%orientation_out = self%softoris%get_ori(self%npeaks) ! best ori
        call self%orientation_out%set('corr', wcorr)  
        ! angular standard deviation
        call self%orientation_out%set('sdev', self%ang_sdev())
        ! dist
        euldist = rad2deg(self%orientation_in.euldist.self%orientation_out)
        call self%orientation_out%set('dist',euldist)
        ! overlap between distributions
        roind    = self%ppftcc%get_roind(360.-self%orientation_out%e3get())
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
            state      = nint(self%orientation_out%get('state'))
            prev_state = nint(self%orientation_in%get('state'))
            if(prev_state == state)then
                mi_state = mi_state + 1.
                mi_joint = mi_joint + 1.
            endif
            mi_joint = mi_joint/3.
        else
            mi_joint = mi_joint/2.
        endif
        call self%orientation_out%set('mi_proj',  mi_proj)
        call self%orientation_out%set('mi_inpl',  mi_inpl)
        call self%orientation_out%set('mi_state', mi_state)
        call self%orientation_out%set('mi_joint', mi_joint)
        ! in-plane distance
        ! make in-plane unit vector
        u(1) = 0.
        u(2) = 1.
        ! calculate previous vec
        mat  = rotmat2d(self%orientation_in%e3get())
        x1   = matmul(u,mat)
        ! calculate new vec
        mat     = rotmat2d(self%orientation_out%e3get())
        x2      = matmul(u,mat)
        call self%orientation_out%set('dist_inpl', rad2deg(myacos(dot_product(x1,x2))))
        ! frac
        frac = 100. * (real(self%neval)/real(self%nrefs))
        call self%orientation_out%set('frac', frac)
        ! done
        if(debug)write(*,*)'simple_pcont3d_srch::prep_softoris done'
    end subroutine prep_softoris

    !>  \brief  standard deviation
    function ang_sdev( self )result( sdev )
        class(pcont3D_srch), intent(inout) :: self
        real    :: sdev
        integer :: nstates, state, pop
        sdev = 0.
        if(self%npeaks < 3)return
        nstates = 0
        do state=1,self%nstates
            pop = self%softoris%get_statepop(state)
            if( pop > 0 )then
                sdev = sdev + ang_sdev_state(state)
                nstates = nstates + 1
            endif
        enddo
        sdev = sdev / real( nstates )
        if( debug ) write(*,'(A)') '>>> PCONT3D_SRCH::CALCULATED ANG_SDEV'
    
        contains
            
            function ang_sdev_state( istate )result( isdev )
                use simple_stat, only: moment
                integer, intent(in)  :: istate
                type(ori)            :: o_best, o
                type(oris)           :: os
                real, allocatable    :: dists(:), ws(:)
                integer, allocatable :: inds(:)
                real                 :: ave, isdev, var
                integer              :: loc(1), alloc_stat, i, ind, n, cnt
                logical              :: err
                isdev = 0.
                inds = self%softoris%get_state(istate)
                n = size(inds)
                if( n < 3 )return ! because one is excluded in the next step & moment needs at least 2 objs
                call os%new( n )
                allocate(ws(n), dists(n-1), stat=alloc_stat)
                call alloc_err('ang_sdev_state; simple_pcont3D_srch', alloc_stat)
                ws    = 0.
                dists = 0.
                ! get best ori
                do i=1,n
                    ind = inds(i)
                    ws(i) = self%softoris%get( ind,'ow')
                    call os%set_ori(i, self%softoris%get_ori(ind))
                enddo
                loc    = maxloc(ws)
                o_best = os%get_ori(loc(1))
                ! build distance vector
                cnt = 0
                do i=1,n
                    if( i==loc(1) )cycle
                    cnt = cnt+1
                    o = os%get_ori(i)
                    dists(cnt) = rad2deg(o.euldist.o_best)
                enddo
                call moment(dists, ave, isdev, var, err)
                deallocate(ws, dists, inds)
                call os%kill
            end function ang_sdev_state
    end function ang_sdev

    ! GETTERS

    !>  \brief  returns the solution set of orientations
    function get_softoris( self )result( os )
        class(pcont3D_srch), intent(inout) :: self
        type(oris) :: os
        if(.not.self%exists)stop 'search has not been performed; pcont3d_srch::get_softoris'
        os = self%softoris
    end function get_softoris

    !>  \brief  returns the best solution orientation
    function get_best_ori( self )result( o )
        class(pcont3D_srch), intent(inout) :: self
        type(ori) :: o
        if(.not.self%exists)stop 'search has not been performed; pcont3d_srch::get_softoris'
        o = self%orientation_out
    end function get_best_ori

    !>  \brief  whether object has been initialized
    function does_exist( self )result( l )
        class(pcont3D_srch), intent(inout) :: self
        logical :: l
        l = self%exists
    end function does_exist

    ! DESTRUCTOR

    !>  \brief  is the destructor
    subroutine kill( self )
        class(pcont3D_srch), intent(inout) :: self
        call self%shsrch_obj%kill
        self%ppftcc => null()
        call self%shiftedoris%kill
        call self%softoris%kill
        call self%reforis%kill
        call self%orientation_in%kill
        call self%orientation_out%kill
        self%iptcl  = 0
        self%exists = .false.
    end subroutine kill

end module simple_pcont3D_srch
