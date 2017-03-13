module simple_pcont3D_srch
use simple_defs
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_pftcc_shsrch      ! use all in there
use simple_prime_srch,       only: prime_srch
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_ctf,              only: ctf
use simple_math              ! use all in there
implicit none

public :: pcont3D_srch
private

real,    parameter :: E3HALFWINSZ = 60.      !< in-plane angle half window size 
logical, parameter :: debug = .false.

type pcont3D_srch
    private
    class(params),           pointer :: pp     => null()            !< pointer to parameters
    class(oris),             pointer :: pe     => null()            !< references orientations
    class(polarft_corrcalc), pointer :: ppftcc => null()
    type(prime_srch) :: srch_common             !< functionalities common to primesrch2D/3D
    type(oris)       :: reforis                 !< per ptcl search space and result of euler angles search
    type(oris)       :: softoris
    type(oris)       :: shiftedoris
    type(ori)        :: orientation_in
    type(ori)        :: orientation_out
    real             :: lims(2,2)      = 0.
    integer          :: nstates        = 0      !< number states
    integer          :: neff_states    = 0      !< number of existing states (eg non empty)
    integer          :: npeaks         = 0
    integer          :: nrefs_per_ptcl = 0      !< number of references per particle
    logical          :: exists = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SEARCH ROUTINES
    procedure          :: do_srch
    procedure, private :: do_refs_srch
    procedure, private :: do_shift_srch
    procedure, private :: prep_softoris
    procedure, private :: ang_sdev
    ! GETTERS
    procedure          :: get_best_ori
    procedure          :: get_softoris
    ! DESTRUCTOR
    procedure          :: kill
end type pcont3D_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, p, o_in, e, pftcc )
        class(pcont3D_srch),             intent(inout) :: self  !< instance
        class(params),           target, intent(in)    :: p     !< parameters
        class(ori),              target, intent(inout) :: o_in  !< builder
        class(oris),             target, intent(in)    :: e     !< builder
        class(polarft_corrcalc), target, intent(in)    :: pftcc
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%pp     => p 
        self%pe     => e
        self%ppftcc => pftcc
        self%npeaks         = self%pp%npeaks
        self%lims(:,1)      = -self%pp%trs
        self%lims(:,2)      = self%pp%trs
        self%nstates        = self%pp%nstates
        self%nrefs_per_ptcl = self%ppftcc%get_nrefs()
        ! orientations
        self%orientation_in  = o_in
        self%orientation_out = self%orientation_in
        call self%reforis%new( self%nrefs_per_ptcl )
        call self%shiftedoris%new(self%npeaks)
        call self%softoris%new( self%npeaks )
        ! other stuff
        self%srch_common = prime_srch(self%pp, self%nrefs_per_ptcl, self%ppftcc%get_nrots())
        ! done
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> PCONT3D_SRCH::CONSTRUCTED NEW SIMPLE_PCONT3D_SRCH OBJECT'
    end subroutine new

    ! SEARCH ROUTINES

    !>  \brief  is the master search routine
    subroutine do_srch( self )
        class(pcont3D_srch), intent(inout) :: self
        ! EULER ANGLES SEARCH
        call self%do_refs_srch
        ! SHIFT SEARCH
        call self%do_shift_srch
        ! outcome prep
        call self%prep_softoris
        if( debug ) write(*,'(A)') '>>> PCONT3D_SRCH::END OF SRCH'
    end subroutine do_srch

    !>  \brief  performs euler angles search
    subroutine do_refs_srch( self )
        class(pcont3D_srch), intent(inout) :: self
        integer, allocatable :: roind_vec(:)    ! slice of in-plane angles
        type(ori) :: o
        real      :: corrs( self%ppftcc%get_nrots() )
        integer   :: loc(1), iref, inpl_ind
        ! init
        corrs = -1.
        !roind_vec = self%ppftcc%get_win_roind(360.-self%orientation_in%e3get(), E3HALFWINSZ)
        ! search
        do iref = 1,self%nrefs_per_ptcl
            o        = self%pe%get_ori( iref )
            corrs    = self%ppftcc%gencorrs(iref, 1)
            loc      = maxloc(corrs)     
            inpl_ind = loc(1)                                  
            call o%set( 'corr', corrs(inpl_ind) )
            call o%e3set( 360.-self%srch_common%rot(inpl_ind) )
            call self%reforis%set_ori( iref, o )
        enddo
        ! done
        !deallocate( roind_vec )
        if(debug)write(*,*)'simple_pcont3d_srch::do_refs_srch done'
    end subroutine do_refs_srch

    !>  \brief  performs the shift search
    subroutine do_shift_srch( self )
        use simple_pftcc_shsrch      ! use all in there
        class(pcont3D_srch), intent(inout) :: self
        integer, allocatable :: inds(:)
        real,    allocatable :: corrs(:)
        type(ori) :: o
        real      :: cxy(3), shift_vec(2), corr
        integer   :: i, iref, inpl_ind
        shift_vec = self%orientation_in%get_shift()
        call pftcc_shsrch_init(self%ppftcc, self%lims)
        ! SELECT BEST NPEAKS ORIENTATIONS
        allocate(inds(self%nrefs_per_ptcl))
        inds  = (/ (i, i=1,self%nrefs_per_ptcl) /)
        corrs = self%reforis%get_all('corr')
        if( any(corrs<-.5) )then
            print *, 'Invalid correlation in simple_pcont3D_srch::do_shift_srch', corrs
        endif
        call hpsort( self%nrefs_per_ptcl, corrs, inds )
        ! SHIFT SEARCH
        do i = 1,self%npeaks
            iref     = inds( self%nrefs_per_ptcl-i+1 )   ! per ptcl ref index
            o        = self%reforis%get_ori( iref )
            corr     = o%get('corr')
            !inpl_ind = self%srch_common%roind( o%e3get() )
            inpl_ind = self%srch_common%roind( 360.-o%e3get() )
            call pftcc_shsrch_set_indices(iref, 1, inpl_ind)
            cxy = pftcc_shsrch_minimize()
            if( cxy(1) >= corr )then
                call o%set('corr', cxy(1))
                call o%set_shift(cxy(2:) + shift_vec)
            else
                call o%set_shift(shift_vec)
            endif
            call self%shiftedoris%set_ori( i, o )
        enddo
        ! done
        deallocate( inds, corrs )
        if(debug)write(*,*)'simple_pcont3d_src::do_shift_srch done'
    end subroutine do_shift_srch

    !>  \brief  updates solutions orientations
    subroutine prep_softoris( self )
        class(pcont3D_srch), intent(inout) :: self
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        type(ori) :: o
        real      :: ws(self%npeaks)
        real      :: frac, dist, maxdist, wcorr, euldist, sdev
        real      :: mi_class, mi_inpl, mi_state, mi_joint
        integer   :: i, cnt, state, roind, prev_roind, prev_state
        ! sort oris
        allocate( inds(self%npeaks) )
        inds  = (/ (i, i=1,self%npeaks) /)
        corrs = self%shiftedoris%get_all('corr')
        call hpsort(self%npeaks, corrs, inds)
        do i = 1,self%npeaks
            o = self%shiftedoris%get_ori( inds(i) )
            call self%softoris%set_ori(i, o)
        enddo
        deallocate( corrs, inds )
        ! stochastic weights and weighted correlation
        ws = 0.
        corrs = self%softoris%get_all('corr')
        where( corrs > TINY ) ws = exp(corrs)
        ws    = ws/sum(ws)
        wcorr = sum(ws*corrs) 
        call self%softoris%set_all( 'ow', ws )
        ! best orientation
        o = self%softoris%get_ori(self%npeaks)
        call self%orientation_out%set('corr', wcorr )  
        call self%orientation_out%set('ow', o%get('ow') )
        call self%orientation_out%set_euler( o%get_euler() )
        call self%orientation_out%set_shift( o%get_shift() )
        call self%orientation_out%set('state', o%get('state') )
        ! angular standard deviation
        sdev = self%ang_sdev()
        call self%orientation_out%set('sdev', sdev)
        call self%softoris%set_all2single('sdev', sdev)
        !
        ! make unit vector
        ! u(1) = 0.
        ! u(2) = 1.
        ! ! calculate previous vec
        ! mat  = rotmat2d(o2update%e3get())
        ! x1   = matmul(u,mat)
        ! ! get new orientation
        ! o = self%o_npeaks%get_ori( self%npeaks )
        ! ! calculate new vec
        ! mat     = rotmat2d(o%e3get())
        ! x2      = matmul(u,mat)
        ! state   = nint( o%get('state') )
        ! if( .not.self%state_exists(state) )then
        !     print *,'Empty state in simple_prime3d_srch; get_ori_best'
        !     stop
        ! endif
        ! dist
        euldist = rad2deg( self%orientation_in.euldist.self%orientation_out )
        call self%orientation_out%set('dist',euldist)
        call self%softoris%set_all2single('sdev', sdev)
        ! calculate overlap between distributions
        roind      = self%srch_common%roind(360.-self%orientation_out%e3get() )
        prev_roind = self%srch_common%roind(360.-self%orientation_in%e3get())
        state      = nint( self%orientation_out%get('state') )
        prev_state = nint( self%orientation_in%get('state') )
        mi_class = 1.
        mi_inpl  = 0.
        mi_state = 0.
        mi_joint = 1.
        if( prev_roind == roind )then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%nstates > 1 )then
            if( prev_state == state )then
                mi_state = mi_state + 1.
                mi_joint = mi_joint + 1.
            endif
            mi_joint = mi_joint/3.
        else
            mi_joint = mi_joint/2.
        endif
        ! overlaps
        call self%orientation_out%set('mi_class', mi_class)
        call self%orientation_out%set('mi_inpl',  mi_inpl)
        call self%orientation_out%set('mi_state', mi_state)
        call self%orientation_out%set('mi_joint', mi_joint)
        ! ! set the distances before we update the orientation
        ! call o2update%set('dist', 0.5*euldist + 0.5*o2update%get('dist'))
        ! call o2update%set('dist_inpl', rad2deg(myacos(dot_product(x1,x2))))
        ! frac
        maxdist = -huge(maxdist)
        do i=1,self%npeaks-1
            o = self%softoris%get_ori(i)
            dist = o.euldist.self%orientation_out
            if( dist>maxdist )maxdist = dist
        enddo
        cnt = 0
        do i=1,self%nrefs_per_ptcl
            o = self%reforis%get_ori(i)
            dist = o.euldist.self%orientation_out
            if( dist<=maxdist )cnt = cnt + 1
        enddo
        cnt  = cnt - self%npeaks
        frac = 100.*real(self%nrefs_per_ptcl-cnt)/real(self%nrefs_per_ptcl)
        call self%orientation_out%set('frac', frac)
        call self%softoris%set_all2single('frac', frac)
        ! CTF
        o = self%orientation_out
        if(o%isthere('dfx'))call self%softoris%set_all2single('dfx',o%get('dfx'))
        if(o%isthere('dfy'))call self%softoris%set_all2single('dfy',o%get('dfy'))
        if(o%isthere('angast'))call self%softoris%set_all2single('angast',o%get('anagst'))
        if(o%isthere('cs'))call self%softoris%set_all2single('cs',o%get('cs'))
        if(o%isthere('kv'))call self%softoris%set_all2single('kv',o%get('kv'))
        if(o%isthere('fraca'))call self%softoris%set_all2single('fraca',o%get('fraca'))
        ! done
        deallocate( corrs )
        if(debug)write(*,*)'simple_pcont3d_src::prep_softoris done'
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
            pop = self%softoris%get_statepop( state )
            if( pop > 0 )then
                sdev = sdev + ang_sdev_state( state )
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
                inds = self%softoris%get_state( istate )
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
                    call os%set_ori( i, self%softoris%get_ori(ind) )
                enddo
                loc    = maxloc( ws )
                o_best = os%get_ori( loc(1) )
                ! build distance vector
                cnt = 0
                do i=1,n
                    if( i==loc(1) )cycle
                    cnt = cnt+1
                    o = os%get_ori( i )
                    dists( cnt ) = rad2deg( o.euldist.o_best )
                enddo
                call moment(dists, ave, isdev, var, err)
                deallocate( ws, dists, inds )
                call os%kill
            end function ang_sdev_state
    end function ang_sdev

    ! GETTERS

    !>  \brief  returns the solution set of orientations
    function get_softoris( self )result( os )
        class(pcont3D_srch), intent(inout) :: self
        type(oris) :: os
        os = self%softoris
    end function get_softoris

    !>  \brief  returns the best solution orientation
    function get_best_ori( self )result( o )
        class(pcont3D_srch), intent(inout) :: self
        type(ori) :: o
        o = self%orientation_out
    end function get_best_ori

    ! DESTRUCTOR

    !>  \brief  is the destructor
    subroutine kill( self )
        class(pcont3D_srch), intent(inout) :: self
        self%pp     => null()
        self%pe     => null()
        self%ppftcc => null()
        call self%srch_common%kill
        call self%shiftedoris%kill
        call self%softoris%kill
        call self%reforis%kill
        call self%orientation_in%kill
        call self%orientation_out%kill
        self%exists = .false.
    end subroutine kill

end module simple_pcont3D_srch
