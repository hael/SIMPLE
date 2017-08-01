module simple_prime3D_srch
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
use simple_ran_tabu,         only: ran_tabu
use simple_jiffys,           only: alloc_err
use simple_params,           only: params
use simple_math,             only: hpsort, rad2deg, round2even, put_last
use simple_pftcc_shsrch      ! singleton
use simple_defs              ! singleton
use simple_cuda_defs
implicit none

public :: prime3D_srch
private

logical, parameter :: debug=.false.

type reference
    integer :: proj  = 1     !< projection direction index
    integer :: state = 1     !< state of reference (conformational or compositional)
    real    :: inpl_corr     !< best in-plane correlation
    integer :: inpl_ind      !< best in-plane index
    real    :: shvec(2) = 0. !< shift vector
    real    :: w = 0.        !< orientation weight
end type

type prime3D_srch
    private
    integer                      :: nrefs      = 0        !< total number of references (nstates*nprojs)
    integer                      :: nstates    = 0        !< number of states
    integer                      :: nprojs     = 0        !< number of projections (same as number of oris in o_refs)
    integer                      :: nrots      = 0        !< number of in-plane rotations in polar representation
    integer                      :: npeaks     = 0        !< number of peaks (nonzero orientation weights)
    integer                      :: nbetter    = 0        !< nr of better orientations identified
    integer                      :: nrefs_eval = 0        !< nr of references evaluated
    integer                      :: proj       = 0        !< previous projection index
    integer                      :: rot        = 0        !< previous in-plane rotation index
    integer                      :: astep      = 1        !< angular step size
    integer                      :: state      = 0        !< previous state index
    integer                      :: ref        = 0        !< previous reference index
    integer                      :: nnn        = 0        !< number of nearest neighbors
    real                         :: atres      = 0.       !< angular treshold
    real                         :: trs        = 0.       !< shift range parameter [-trs,trs]
    real                         :: shvec(2)   = 0.       !< previous origin shift vector
    real                         :: corr_prev  = 1.       !< previous best correlation
    class(oris),   pointer       :: o_refs => null()      !< pointer to orientations of references in proj_space
    class(params), pointer       :: pp => null()          !< pointer to parameters singleton
    type(reference), allocatable :: proj_space(:)         !< array of references (the space of projection directions)
    integer, allocatable         :: proj_space_inds(:)    !< projection space index array
    integer, allocatable         :: refinds(:,:)          !< table of reference indices
    integer, allocatable         :: srch_order(:)         !< stochastoc search order
    integer, allocatable         :: state_range(:,:)      !< state mapping to reference range
    integer, allocatable         :: inplmat(:,:)          !< in-plane rotation index matrix for GPU-based search
    real,    allocatable         :: proj_space_corrs(:)   !< projection space correlations (in-plane maxima)
    real,    allocatable         :: angtab(:)             !< list of in-plane angles
    real,    allocatable         :: prev_corrs(:)         !< previous particle-reference correlations
    real,    allocatable         :: corrmat3d(:,:,:)      !< 3D correlation matrix for GPU-based search
    real,    allocatable         :: corrmat2d(:,:)        !< 2D correlation matrix for GPU-based search
    type(ran_tabu)               :: rt                    !< object for random number generation
    character(len=STDLEN)        :: refine = ''           !< refinement flag
    character(len=STDLEN)        :: ctf = ''              !< ctf flag
    logical                      :: use_cpu  = .true.     !< indicates if CPU logics is being used
    logical                      :: bench_gpu = .false.   !< indicates if GPU logics should be tested
    logical                      :: use_gpu  = .false.    !< indicates if GPU logics should be used
    logical                      :: doshift  = .true.     !< origin shift search indicator
    logical                      :: exists   = .false.    !< 2 indicate existence
    ! benchmark control logicals data structure
    type(t_bench)                :: s_bench
  contains
    ! CONSTRUCTOR & INITIALIZER
    procedure          :: new
    procedure          :: init
    ! CALCULATORS
    procedure, private :: calc_roind
    procedure, private :: calc_stateind
    procedure, private :: calc_rot
    ! GETTERS & SETTERS
    procedure          :: get_ori_best
    procedure          :: get_state
    procedure          :: get_proj
    procedure          :: get_class
    procedure          :: get_nbetter
    procedure          :: get_ori
    procedure          :: get_euler
    procedure          :: ang_sdev
    ! PREPARATION ROUTINES
    procedure, private :: prep4srch
    procedure          :: calc_corrs_on_gpu
    ! SEARCH ROUTINES
    procedure          :: exec_prime3D_srch
    procedure          :: exec_prime3D_qcont_srch
    procedure          :: exec_prime3D_shc_srch
    procedure, private :: stochastic_srch
    procedure, private :: stochastic_srch_qcont
    procedure, private :: stochastic_srch_shc
    procedure, private :: shift_srch
    procedure, private :: shift_srch_qcont
    procedure, private :: shift_srch_shc
    procedure, private :: greedy_srch
    procedure, private :: stochastic_weights
    ! DESTRUCTOR
    procedure :: kill
end type prime3D_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, e, p )
        class(prime3D_srch), intent(inout) :: self !< instance
        class(params), intent(in), target  :: p    !< parameters
        class(oris), intent(in), target    :: e    !< Euler angles
        integer :: alloc_stat, cnt, i, s
        real    :: dang, ang, x
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%pp         => p 
        self%nrefs      =  p%nspace*p%nstates
        self%nstates    =  p%nstates
        self%nprojs     =  p%nspace
        self%nrots      =  round2even(twopi*real(p%ring2))
        self%npeaks     =  p%npeaks
        self%atres      =  p%thres
        self%nbetter    =  0
        self%nrefs_eval =  0
        self%trs        =  p%trs
        self%o_refs     => e
        self%doshift    =  p%doshift
        self%refine     =  p%refine
        self%ctf        =  p%ctf
        self%nnn        =  p%nnn
        self%bench_gpu  =  .false.
        self%use_gpu    =  .false.
        self%use_cpu    =  .true.
        if( self%refine .eq. 'no' )then 
            self%astep = nint(atan(p%lp/(4.*real(p%ring2)*p%smpd))/(twopi/real(self%nrots)))
            self%astep = max(1,self%astep)
            self%astep = min(5,self%astep)
        else
            self%astep = 1
        endif
        if( p%bench_gpu .eq. 'yes' .or. p%use_gpu .eq. 'yes' )then
            if( p%bench_gpu .eq. 'yes' ) self%bench_gpu = .true.
            if( p%use_gpu  .eq. 'yes' ) self%use_gpu  = .true.
            if( p%top-p%fromp+1 /= self%nrefs )&
            &stop 'the particle chunk is not correctly balanced for GPU execution!'
             self%use_cpu = .false.
        endif
        ! construct composites
        allocate( self%proj_space(self%nrefs), self%refinds(self%nstates,self%nprojs),&
                  self%proj_space_inds(self%nrefs), self%proj_space_corrs(self%nrefs),&
                  self%angtab(self%nrots), self%state_range(self%nstates,2), stat=alloc_stat )
        call alloc_err('In: new; simple_prime3D_srch, 1', alloc_stat)
        self%proj_space_inds  = 0
        self%proj_space_corrs = -1.
        if( self%refine.eq.'neigh' .or. self%refine.eq.'shcneigh' )then ! local refinement
            allocate(self%srch_order(self%nnn), stat=alloc_stat)
            call alloc_err('In: new; simple_prime3D_srch, 2', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu(self%nnn)
        else
            allocate(self%srch_order(self%nrefs), stat=alloc_stat)
            call alloc_err('In: new; simple_prime3D_srch, 3', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu(self%nrefs)
        endif
        ! generate the list of in-plane angles and indices
        dang = twopi/real(self%nrots)
        do i=1,self%nrots
            self%angtab(i) = rad2deg((i-1)*dang)
        end do
        ! initialize structures
        call self%init
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine
    
    !>  \brief  is an initializer
    subroutine init( self )
        class(prime3D_srch), intent(inout) :: self !< instance
        integer :: cnt, i, s
        ! create index transformation & output data (correlation) storage structure
        cnt = 0
        self%state_range(1,1) = 1
        do s=1,self%nstates
            if( s > 1 ) self%state_range(s,1) = cnt+1
            do i=1,self%nprojs
                cnt = cnt+1
                self%proj_space(cnt)%state     = s
                self%proj_space(cnt)%proj      = i
                self%refinds(s,i)              = cnt
                self%proj_space(cnt)%inpl_corr = -1.
                self%proj_space(cnt)%inpl_ind  = 0
                self%proj_space_inds(cnt)      = i
            end do
            self%state_range(s,2) = cnt
        end do
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::INITIALISED SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine
    
    ! CALCULATORS
    
    !>  \brief calculate the in-plane rotational index for the rot in-plane angle
    function calc_roind( self, rot ) result( ind )
        class(prime3D_srch), intent(in) :: self
        real, intent(in)                :: rot
        integer :: ind, i, loc(1)
        real    :: dists(self%nrots)
        do i=1,self%nrots
            dists(i) = sqrt((self%angtab(i)-rot)**2.)
        end do
        loc = minloc(dists)
        ind = loc(1)
    end function
    
    !>  \brief calculate the state index for a given reference index
    function calc_stateind( self, ref ) result( state_ind )
        class(prime3D_srch), intent(in) :: self
        integer, intent(in)             :: ref
        integer :: state_ind, s
        do s=1,self%nstates
            if( ref >= self%state_range(s,1) .and. ref <= self%state_range(s,2) )then
                state_ind = s
                return
            endif
        end do
    end function
        
    !>  \brief calculate the in-plane angle of in-plane rotational index rot
    function calc_rot( self, ind ) result( rot )
        class(prime3D_srch), intent(in) :: self
        integer, intent(in) :: ind
        real :: rot
        if( ind < 1 .or. ind > self%nrots )then
            write(*,*) 'rotational index is: ', ind, ' which is out of range; calc_rot; simple_prime3D_srch'
            stop
        endif
        rot = self%angtab(ind)
    end function
    
    ! GETTERS & SETTERS
    
    !>  \brief  to get the best orientation from the discrete space
    subroutine get_ori_best( self, o )
        use simple_ori,  only: ori
        class(prime3D_srch), intent(in) :: self
        class(ori), intent(inout)       :: o
        type(ori)                       :: o_copy
        real                            :: euls(3), mi, frac, dist
        character(len=STDLEN)           :: dig, key
        integer                         :: proj_new, rot_new, state_new, class, s, ref
        real                            :: x, y, corr, w
        ! copy the input orientation
        o_copy = o
        ! set the reference index
        if( self%refine.eq.'shc' .or. self%refine.eq.'shcneigh' )then
            ref = self%nrefs
        else
            ref = self%proj_space_inds(self%nrefs)
        endif
        ! get new indices
        proj_new  = self%proj_space(ref)%proj
        rot_new   = self%proj_space(ref)%inpl_ind
        state_new = self%proj_space(ref)%state
        class     = ref
        if( self%refine.eq.'shc' .or. self%refine.eq.'shcneigh' .or. (self%refine.eq.'no'.and.self%nstates>1))then
            class = self%refinds(state_new,proj_new)
        endif
        euls(1)   = self%o_refs%e1get(proj_new)
        euls(2)   = self%o_refs%e2get(proj_new)
        if( rot_new < 1 .or. rot_new > self%nrots )then
            write(*,*) 'rot_new  in simple_prime3D_srch::get_ori_best is: ', rot_new, ' which is out of range'
            stop
        endif
        euls(3) = 360.-self%calc_rot(rot_new) ! change sgn to fit convention
        if( euls(3) == 360. ) euls(3) = 0.
        ! calculate overlap between distributions
        mi = 0.
        if( self%proj  == proj_new  ) mi = mi+1
        if( self%rot   == rot_new   ) mi = mi+1
        if( self%state == state_new ) mi = mi+1
        mi = mi/3.
        ! set parameters
        call o%set_euler(euls)
        ! shifts must be obtained by vector addition
        x = self%shvec(1)+self%proj_space(ref)%shvec(1)
        y = self%shvec(2)+self%proj_space(ref)%shvec(2)
        call o%set('x',x)
        call o%set('y',y)
        corr = self%proj_space(ref)%inpl_corr
        w    = self%proj_space(ref)%w
        if( self%refine.eq.'neigh' .or. self%refine.eq.'shcneigh' )then
            frac = 100.*(real(self%nrefs_eval)/real(self%nnn))
        else
            frac = 100.*(real(self%nrefs_eval)/real(self%nrefs))
        endif
        dist = 0.5*rad2deg(o_copy.euldist.o)+0.5*o_copy%get('dist')
        call o%set_list(['state  ','class  ','corr   ','dist   ','mi_hard','frac   ', 'ow     '],&
        [real(state_new),real(class),corr,dist,mi,frac,w])
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::GOT BEST ORI'
    end subroutine
    
    !>  \brief  get the state corresponding to class
    function get_state( self, class ) result( s )
        class(prime3D_srch), intent(inout) :: self
        integer, intent(in)              :: class
        integer                          :: s
        s = self%proj_space(class)%state
    end function
    
    !>  \brief  get the projection direction corresponding to class
    function get_proj( self, class ) result( p )
        class(prime3D_srch), intent(inout) :: self
        integer, intent(in)                :: class
        integer                            :: p
        p = self%proj_space(class)%proj
    end function
    
    !>  \brief  get the projection direction corresponding to class
    function get_class( self, state, proj ) result( class )
        class(prime3D_srch), intent(inout) :: self
        integer, intent(in)                :: state, proj
        integer                            :: class
        class = self%refinds(state,proj)
    end function
    
    !>  \brief  the number of better orientations
    function get_nbetter( self ) result( nbetter )
        class(prime3D_srch), intent(inout) :: self
        integer :: nbetter
        nbetter = max(1,self%nbetter)
    end function
    
    !>  \brief  to get one orientation from the discrete space
    subroutine get_ori( self, i, o )
        use simple_ori, only: ori
        class(prime3D_srch), intent(inout) :: self
        integer, intent(in)              :: i
        type(ori), intent(inout)         :: o
        real                             :: euls(3), x, y
        integer                          :: ref
        if( self%refine.eq.'shc' .or. self%refine.eq.'shcneigh' )then
            stop 'get_ori not for shc-modes; simple_prime3D_srch'
        endif 
        ref = self%proj_space_inds(self%nrefs-i+1)
        ! shifts must be obtained by vector addition
        x = self%shvec(1)+self%proj_space(ref)%shvec(1)
        y = self%shvec(2)+self%proj_space(ref)%shvec(2)
        if( abs(x) < 1e-6 ) x = 0.
        if( abs(y) < 1e-6 ) y = 0.
        call o%set('x',x)
        call o%set('y',y)
        euls(1) = self%o_refs%e1get(self%proj_space(ref)%proj)
        euls(2) = self%o_refs%e2get(self%proj_space(ref)%proj)
        euls(3) = 360.-self%calc_rot(self%proj_space(ref)%inpl_ind) ! change sgn to fit convention
        call o%set_euler(euls)
        call o%set('state',real(self%proj_space(ref)%state))
        call o%set('corr',self%proj_space(ref)%inpl_corr)
        call o%set('ow',self%proj_space(ref)%w)
    end subroutine
    
    !>  \brief  to get one orientation from the discrete space
    subroutine get_euler( self, i, o )
        use simple_ori, only: ori
        class(prime3D_srch), intent(inout) :: self
        integer, intent(in)              :: i
        type(ori), intent(inout)         :: o
        real                             :: euls(3)
        integer                          :: ref
        ref = self%proj_space_inds(self%nrefs-i+1)
        if( self%refine.eq.'shc' .or. self%refine.eq.'shcneigh' )then
            stop 'get_euler not for shc-modes; simple_prime3D_srch'
        endif 
        euls(1) = self%o_refs%e1get(self%proj_space(ref)%proj)
        euls(2) = self%o_refs%e2get(self%proj_space(ref)%proj)
        euls(3) = 360.-self%calc_rot(self%proj_space(ref)%inpl_ind) ! change sgn to fit convention
        call o%set_euler(euls)
        call o%set('state',real(self%proj_space(ref)%state))
        call o%set('corr',self%proj_space(ref)%inpl_corr)
    end subroutine
    
    !>  \brief  standard deviation
    function ang_sdev( self ) result( sdev )
        use simple_ori, only: ori
        class(prime3D_srch), intent(inout) :: self
        real                             :: sdev
        integer                          :: nstates, s, i, cnt, alloc_stat, ref
        type(reference)                  :: best_ones(self%npeaks)
        if( self%refine.eq.'shc' .or. self%refine.eq.'shcneigh' )then
            sdev = 0.
            return
        endif
        cnt = 0
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            cnt = cnt+1
            ref = self%proj_space_inds(i)
            best_ones(cnt) = self%proj_space(ref)   
        end do
        nstates = maxval(best_ones%state)
        sdev = 0.
        do s=1,nstates
            sdev = sdev+ang_sdev_state()
        end do
        sdev = sdev/real(nstates)
        cnt = 0
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::CALCULATED ANG_SDEV'
    
        contains

            function ang_sdev_state( ) result( sdev )
                use simple_stat, only: moment
                real                        :: ave, sdev, var
                type(ori)                   :: o_best, o
                integer                     :: loc(1), alloc_stat, i, cnt
                real, allocatable           :: dists(:)
                logical                     :: err
                real                        :: w
                call o_best%new
                call o%new
                ! count the number of occurances of state
                cnt = 0
                do i=1,size(best_ones)
                    if( best_ones(i)%state == s ) cnt = cnt+1
                end do
                if( cnt <= 3 ) then ! because one is excluded in the next step & moment needs at least 2 objs
                    sdev = 0.
                    return
                endif
                allocate( dists(cnt-1), stat=alloc_stat )
                call alloc_err( 'ang_sdev_state; simple_prime3D_srch', alloc_stat )
                ! copy data
                cnt = 0
                w = 0.
                loc(1) = 0
                do i=1,size(best_ones)
                    if( best_ones(i)%state == s )then
                        cnt = cnt+1
                        if( best_ones(i)%w > w )then
                            loc(1) = i
                            w = best_ones(i)%w
                        endif
                    endif
                end do
                if( loc(1) == 0 )then
                    sdev = 0.
                    return
                else
                    call self%get_euler(loc(1), o_best)
                    cnt = 0
                    do i=1,size(best_ones)
                        if( i == loc(1) ) cycle
                        if( best_ones(i)%state == s )then
                            cnt = cnt+1
                            call self%get_euler(i, o)
                            dists(cnt) = rad2deg(o_best.euldist.o) 
                        endif
                    end do
                    call moment(dists(:cnt), ave, sdev, var, err)
                endif
                deallocate(dists)
                call o_best%kill
                call o%kill
            end function
            
    end function
    
    ! PREPARATION ROUTINES
    
    !>  \brief  prepares for the search on CPU
    subroutine prep4srch( self, pftcc, iptcl, lp, o_prev )
        use simple_ori, only: ori
        class(prime3D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        class(ori), intent(inout), optional    :: o_prev
        logical :: did_set_corr_prev
        integer :: s
        real    :: lims(2,2), cc_t, cc_t_min_1
        ! make 'random' projection direction order
        call self%rt%reset
        call self%rt%ne_ran_iarr(self%srch_order)
        ! initialize origin shift search object
        lims(1,1) = -self%trs
        lims(1,2) = self%trs
        lims(2,1) = -self%trs
        lims(2,2) = self%trs
        ! we are using simplex with barrier constraint
        call pftcc_shsrch_init(pftcc, 'simplex', lims, 3)
        if( present(o_prev) )then
            ! find previous discrete alignment parameters
            self%proj  = self%o_refs%find_closest_proj(o_prev) ! projection index
            self%rot   = self%calc_roind(360.-o_prev%e3get())  ! in-plane angle index
            self%state = nint(o_prev%get('state'))             ! state index
            self%ref   = self%refinds(self%state,self%proj)    ! reference index
            self%shvec = [o_prev%get('x'),o_prev%get('y')]     ! shift vector
            if( self%state > self%nstates )then
                stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
            endif
            ! put prev_best last to avoid cycling
            call put_last(self%ref, self%srch_order)
            ! calculate previous best corr (treshold for better)
            did_set_corr_prev = .false.
            cc_t = pftcc%corr(self%ref, iptcl, self%rot)
            if( self%refine.eq.'no' .and. self%ctf.eq.'no' )then        
                cc_t_min_1 = -1.
                if( o_prev%isthere('corr') ) cc_t_min_1 = o_prev%get('corr')
                if( o_prev%isthere('lp') )then
                    if( abs(o_prev%get('lp')-lp) < 0.5 )then ! previous and present correlation values are comparable
                        if( cc_t_min_1 > 0. )then
                            self%corr_prev = 0.5*cc_t+0.5*cc_t_min_1 ! diversifying limit
                            did_set_corr_prev = .true.
                        endif
                    endif
                endif
                if( .not. did_set_corr_prev ) self%corr_prev = cc_t ! greedy limit
            else
                self%corr_prev = cc_t
            endif
        else
            self%proj      = 1
            self%rot       = 1
            self%state     = 1
            self%ref       = 1
            self%shvec     = 0.
            self%corr_prev = 1.
        endif
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::PREPARED FOR SIMPLE_PRIME3D_SRCH'
    end subroutine
    
    !>  \brief  prepares the matrices for PRIME3D search using Frederic's GPU code
    subroutine calc_corrs_on_gpu( self, pftcc, a, testflag )
        use simple_oris,    only: oris
        use simple_ori,     only: ori
        use simple_jiffys,  only: alloc_err
        use simple_corrmat, only: project_corrmat3D_greedy, project_corrmat3D_shc
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(oris),             intent(inout) :: a
        logical, optional,       intent(in)    :: testflag
        integer :: alloc_stat
        logical :: ttestflag
        ttestflag = .false.
        if( present(testflag) ) ttestflag = testflag
        ! the resulting corrmat2d and inplmat 2D matrices are indexed (iptcl,iref)        
        allocate(self%corrmat3d(self%nrefs,self%nrefs,self%nrots),&
                 self%corrmat2d(self%nrefs,self%nrefs),&
                 self%inplmat(self%nrefs,self%nrefs), stat=alloc_stat)
        call alloc_err("simple_prime3D_srch::prep_matrices", alloc_stat)
        ! we need to expand the dimensions in pftcc before calculating correlations on GPU
        call pftcc%expand_dim
        ! calculate all correlations needed for search
        if( ttestflag )then
            call pftcc%gencorrs_all_tester(self%corrmat2d, self%inplmat)
            deallocate(self%corrmat3d)
        else
            call pftcc%gencorrs_all(self%pp, self%corrmat3d)
            ! process corrmat3d to accomplish the different modes of stochastic search
            if( self%refine .eq. 'no' )then
                call project_corrmat3D_greedy(self%nrefs, self%corrmat3d, self%corrmat2d, self%inplmat)
                deallocate(self%corrmat3d)
            else if( self%refine .eq. 'shc' )then
                call project_corrmat3D_shc(self%nrefs, self%corrmat3d, self%prev_corrs, self%corrmat2d, self%inplmat)
                deallocate(self%corrmat3d)
            else
                write(*,*) 'Refinement mode: ', trim(self%refine)
                stop 'This refinement mode is not currently supported on GPU'
            endif
        endif            
    end subroutine calc_corrs_on_gpu

    ! SEARCH ROUTINES
    
    !>  \brief a master prime search routine 4 CPU
    subroutine exec_prime3D_srch( self, pftcc, iptcl, lp, wcorr, o, nnmat, cnt_glob )
        use simple_ori, only: ori
        class(prime3D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        real, intent(out)                      :: wcorr
        class(ori), intent(inout), optional    :: o
        integer, intent(in), optional          :: nnmat(self%nprojs,self%nnn)
        integer, intent(in), optional          :: cnt_glob
        call self%prep4srch(pftcc, iptcl, lp, o)
        call self%stochastic_srch(pftcc, iptcl, nnmat, cnt_glob)
        call self%stochastic_weights(wcorr)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine
    
    !>  \brief a master prime search routine 4 CPU
    subroutine exec_prime3D_qcont_srch( self, refvols, proj, tfun, pftcc, iptcl, lp, wcorr, o )
        use simple_image,     only: image
        use simple_projector, only: projector
        use simple_ori,       only: ori
        use simple_ctf,       only: ctf
        class(prime3D_srch), intent(inout)     :: self
        class(image), intent(inout)            :: refvols(self%nstates)
        class(projector), intent(inout)        :: proj
        class(ctf), intent(inout)              :: tfun
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        real, intent(out)                      :: wcorr
        class(ori), intent(inout)              :: o
        call self%prep4srch(pftcc, iptcl, lp, o)
        call self%stochastic_srch_qcont(refvols, proj, tfun, pftcc, iptcl, o)
        call self%stochastic_weights(wcorr)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine
    
    !>  \brief a master prime search routine 4 CPU
    subroutine exec_prime3D_shc_srch( self, pftcc, iptcl, lp, corr_best, o, nnmat, cnt_glob )
        use simple_ori, only: ori
        class(prime3D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        real, intent(out)                      :: corr_best
        class(ori), intent(inout), optional    :: o
        integer, intent(in), optional          :: nnmat(self%nprojs,self%nnn)
        integer, intent(in), optional          :: cnt_glob
        call self%prep4srch(pftcc, iptcl, lp, o)
        if( present(o) )then
            call self%stochastic_srch_shc(pftcc, iptcl, corr_best, nnmat, cnt_glob)
        else
            call self%greedy_srch(pftcc, iptcl, corr_best, cnt_glob)
        endif
        if( debug ) write(*,'(A)') '>>> PRIME3D_SHC_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine
    
    !>  \brief  executes the stochastic soft orientation search on CPU
    subroutine stochastic_srch( self, pftcc, iptcl, nnmat, cnt_glob )
        class(prime3D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        integer, intent(in), optional          :: nnmat(self%nprojs,self%nnn)
        integer, intent(in), optional          :: cnt_glob
        integer :: i, loc(1), ref, endit
        real    :: corr, corrs(self%nrots)
        if( present(nnmat) )then
            if( self%bench_gpu .or. self%use_gpu ) stop 'Neigh modes not implemented for GPU; simple_prime3D :: stochastic_srch'
        endif
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime3D :: stochastic_srch'
        endif
        ! initialize
        self%nbetter          = 0
        self%proj_space_corrs = -1.
        self%nrefs_eval       = 0
        self%proj_space_inds  = 0
        ! search
        endit = self%nrefs
        if( present(nnmat) )endit = self%nnn
        do i=1,endit
            ! set the stochastic reference index
            if( present(nnmat) )then
                ref = nnmat(self%proj,self%srch_order(i))
            else
                ref = self%srch_order(i)
            endif
            if( ref < 1 .or. ref > self%nrefs ) stop 'ref index out of bound; simple_prime3D_srch::stochastic_srch'
            ! stash the stochastic reference index
            self%proj_space_inds(ref) = ref
            if( self%use_cpu )then
                corrs = pftcc%gencorrs(ref, iptcl, self%astep)
                loc   = maxloc(corrs)
                self%proj_space(ref)%inpl_ind  = loc(1)
                self%proj_space(ref)%inpl_corr = corrs(loc(1))
            else
                self%proj_space(ref)%inpl_ind  = self%inplmat(cnt_glob,ref)
                self%proj_space(ref)%inpl_corr = self%corrmat2d(cnt_glob,ref)
            endif
            ! update the correlation
            self%proj_space_corrs(ref) = self%proj_space(ref)%inpl_corr
            ! update nbetter to keep track of how many improving solutions we have identified
            if( self%npeaks == 1 )then
                if( self%proj_space_corrs(ref) > self%corr_prev ) self%nbetter = self%nbetter+1
            else
                if( self%proj_space_corrs(ref) >= self%corr_prev ) self%nbetter = self%nbetter+1
            endif
            ! keep track of how many references we are evaluating
            self%nrefs_eval = self%nrefs_eval+1
            ! exit condition
            if( self%nbetter == self%npeaks ) exit
        end do
        ! sort in projection direction space (the max inpl corrs)
        call hpsort(self%nrefs, self%proj_space_corrs, self%proj_space_inds)
        ! search shifts
        call self%shift_srch(iptcl)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'
    end subroutine

    !>  \brief  executes the stochastic quasi-continuous soft orientation search on CPU
    subroutine stochastic_srch_qcont( self, refvols, proj, tfun, pftcc, iptcl, o )
        use simple_image,     only: image
        use simple_projector, only: projector
        use simple_ori,       only: ori
        use simple_ctf,       only: ctf
        use simple_rnd,       only: shcloc
        class(prime3D_srch), intent(inout)     :: self
        class(image), intent(inout)            :: refvols(self%nstates)
        class(projector), intent(inout)        :: proj
        class(ctf), intent(inout)              :: tfun
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        type(ori), intent(inout)               :: o
        type(ori) :: osample
        integer   :: isample, loc(1), mystate
        real      :: corr, corrs(self%nrots), dfx, dfy, angast
        ! initialize
        mystate               = 1
        self%nbetter          = 0
        self%proj_space_corrs = -1.
        self%nrefs_eval       = 0
        self%proj_space_inds  = 0
        if( self%pp%ctf .ne. 'no' )then
            if( self%pp%ctfmode .eq. 'astig' )then ! astigmatic CTF
                dfx = o%get('dfx')
                dfy = o%get('dfy')
                angast = o%get('angast')
            else if( self%pp%ctfmode .eq. 'noastig' )then
                dfx = o%get('dfx')
                dfy = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; stochastic_srch_qcont; simple_prime3D_srch'
            endif
        endif
        ! generate stochastic search space
        call self%o_refs%rnd_proj_space(self%nrefs)
        ! search
        do isample=1,self%nrefs
            ! get random ori sample
            osample = self%o_refs%get_ori(isample)
            ! online generation of polar central section
            call proj%fproject_polar(1,refvols(mystate),osample,self%pp,pftcc)
            if( self%pp%ctf .ne. 'no' ) call pftcc%apply_ctf(tfun, dfx, dfy, angast )
            ! stash the index
            self%proj_space_inds(isample) = isample
            ! now the logics becomes the same as above
            corrs = pftcc%gencorrs(1, iptcl)
            loc   = maxloc(corrs)
            self%proj_space(isample)%inpl_ind  = loc(1)
            self%proj_space(isample)%inpl_corr = corrs(loc(1))
            ! update the correlation
            self%proj_space_corrs(isample) = self%proj_space(isample)%inpl_corr
            ! update nbetter to keep track of how many improving solutions we have identified
            if( self%npeaks == 1 )then
                if( self%proj_space_corrs(isample) > self%corr_prev ) self%nbetter = self%nbetter+1
            else
                if( self%proj_space_corrs(isample) >= self%corr_prev ) self%nbetter = self%nbetter+1
            endif
            ! keep track of how many references we are evaluating
            self%nrefs_eval = self%nrefs_eval+1
            ! exit condition
            if( self%nbetter == self%npeaks ) exit
        end do
        ! sort in projection direction space (the max inpl corrs)
        call hpsort(self%nrefs, self%proj_space_corrs, self%proj_space_inds)
        ! search shifts
        call self%shift_srch_qcont(refvols, proj, pftcc, mystate, iptcl)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED QUASI-CONTINUOUS STOCHASTIC SEARCH'
    end subroutine
    
    !>  \brief  executes the stochastic hard orientation search using pure stochastic hill climbing
    !!          (no probabilistic weighting + stochastic search of in-plane angles) on CPU
    subroutine stochastic_srch_shc( self, pftcc, iptcl, corr_best, nnmat, cnt_glob )
        use simple_rnd, only: irnd_uni, shcloc
        class(prime3D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(out)                      :: corr_best
        integer, intent(in), optional          :: nnmat(self%nprojs,self%nnn)
        integer, intent(in), optional          :: cnt_glob
        integer :: i, ref, endit, proj, inpl
        real    :: corrs(self%nrots), corr
        logical :: found_better
        if( present(nnmat) )then
            if( self%bench_gpu .or. self%use_gpu ) stop 'Neigh modes not implemented for GPU; simple_prime3D :: stochastic_srch'
        endif
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime3D :: stochastic_srch'
        endif
        ! initialize
        self%nrefs_eval = 0
        ! search
        endit = self%nrefs
        if( present(nnmat) ) endit = self%nnn
        found_better = .false.
        if( debug ) print *,'self%nrefs',self%nrefs
        do i=1,endit
            ! set the stochastic reference index
            if( present(nnmat) )then
                if( self%nstates > 1 )then
                    ref = self%refinds(irnd_uni(self%nstates),nnmat(self%proj,self%srch_order(i)))
                else
                    ref = nnmat(self%proj,self%srch_order(i))
                endif
            else
                ref = self%srch_order(i)
            endif
            if( ref < 1 .or. ref > self%nrefs ) stop 'ref index out of bound; simple_prime3D_srch::stochastic_srch_shc'
            ! align using SHC
            if( self%use_cpu )then
                corrs = pftcc%gencorrs(ref, iptcl)
                inpl  = shcloc(corrs, self%corr_prev)
                corr = 0.
                if( inpl > 0 ) corr = corrs(inpl)
            else
                inpl = self%inplmat(cnt_glob,ref)
                corr = self%corrmat2d(cnt_glob,ref)
            endif
            self%nrefs_eval = self%nrefs_eval+1
            if( inpl == 0 ) cycle
            ! update reference index
            self%proj_space(self%nrefs)%proj = self%proj_space_inds(ref)
            ! update correlation
            self%proj_space(self%nrefs)%inpl_corr = corr
            ! update in-plane angle index
            self%proj_space(self%nrefs)%inpl_ind = inpl
            ! update state index
            self%proj_space(self%nrefs)%state = self%calc_stateind(ref)
            ! set best weight to one
            self%proj_space(self%nrefs)%w = 1.
            ! indicate that we found a better solution
            found_better = .true.
            exit ! first-improvement heuristic, characteristic of stochastic hill climbing
        end do
        if( found_better )then
            ! best ref has already been updated
        else
            ! keep the old parameters
            self%proj_space(self%nrefs)%proj      = self%proj
            self%proj_space(self%nrefs)%inpl_corr = self%corr_prev
            self%proj_space(self%nrefs)%inpl_ind  = self%rot
            self%proj_space(self%nrefs)%state     = self%state
            self%proj_space_inds(self%nrefs)      = self%nprojs ! proj_space_inds must have values <= nprojs (=nspace)
        endif
        ! set corr_best
        corr_best = self%proj_space(self%nrefs)%inpl_corr
        ! search shifts
        call self%shift_srch_shc(iptcl)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH (REFINE=SHC)'
    end subroutine stochastic_srch_shc
    
    !>  \brief  greedy hill climbing (4 initialisation)
    subroutine greedy_srch( self, pftcc, iptcl, corr_best, cnt_glob )
        use simple_rnd, only: irnd_uni, shcloc
        class(prime3D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(out)                      :: corr_best
        integer, intent(in), optional          :: cnt_glob
        integer :: i, ref, proj, inpl, loc(1)
        real    :: corrs(self%nrots), corr
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime3D :: stochastic_srch'
        endif
        ! initialize
        self%nrefs_eval = 0
        corr_best       = -1.
        do i=1,self%nrefs
            ! set reference index
            ref = i
            ! align
            if( self%use_cpu )then
                corrs = pftcc%gencorrs(ref, iptcl, self%astep)
                loc   = maxloc(corrs)
                inpl  = loc(1)
                corr  = corrs(inpl)
            else
                inpl = self%inplmat(cnt_glob,ref)
                corr = self%corrmat2d(cnt_glob,ref)
            endif
            self%nrefs_eval = self%nrefs_eval+1
            if( corr > corr_best )then
                corr_best = corr
                ! update reference index
                self%proj_space(self%nrefs)%proj = self%proj_space_inds(ref)
                ! update correlation
                self%proj_space(self%nrefs)%inpl_corr = corr
                ! update in-plane angle index
                self%proj_space(self%nrefs)%inpl_ind = inpl
                ! update state index
                self%proj_space(self%nrefs)%state = self%calc_stateind(ref)
                ! set best weight to one
                self%proj_space(self%nrefs)%w = 1.
            endif
        end do
        ! search shifts
        call self%shift_srch_shc(iptcl)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH 4 INITI (REFINE=SHC)'
    end subroutine greedy_srch
    
    !>  \brief  executes the simplex-based shift search
    subroutine shift_srch( self, iptcl )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        integer :: i, ref
        real    :: cxy(3)
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            if( self%doshift )then
                call pftcc_shsrch_set_indices(ref, iptcl, self%proj_space(ref)%inpl_ind)
                cxy = pftcc_shsrch_minimize()
                self%proj_space(ref)%shvec = cxy(2:3)
                self%proj_space(ref)%inpl_corr = cxy(1)
            else
                self%proj_space(ref)%shvec = [0.,0.]
            endif
        end do
    end subroutine
    
    !>  \brief  executes the simplex-based shift search
    subroutine shift_srch_qcont( self, refvols, proj, pftcc, state, iptcl )
        use simple_image,     only: image
        use simple_projector, only: projector
        use simple_ori,       only: ori
        class(prime3D_srch),     intent(inout) :: self
        class(image),            intent(inout) :: refvols(self%nstates)
        class(projector),        intent(inout) :: proj
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: state, iptcl
        integer   :: i, ref
        real      :: cxy(3)
        type(ori) :: o
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            if( self%doshift )then
                ! get random ori sample
                o = self%o_refs%get_ori(ref)
                ! online generation of polar central section
                call proj%fproject_polar(1,refvols(state),o,self%pp,pftcc)
                call pftcc_shsrch_set_indices(1,iptcl,self%proj_space(ref)%inpl_ind)
                cxy = pftcc_shsrch_minimize()
                self%proj_space(ref)%shvec = cxy(2:3)
                self%proj_space(ref)%inpl_corr = cxy(1)
            else
                self%proj_space(ref)%shvec = [0.,0.]
            endif
        end do
    end subroutine
    
    !>  \brief  executes the simplex-based shift search
    subroutine shift_srch_shc( self, iptcl )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real :: cxy(3)
        ! search shifts
        if( self%doshift )then
            call pftcc_shsrch_set_indices(self%proj_space(self%nrefs)%proj, iptcl, self%proj_space(self%nrefs)%inpl_ind)
            cxy = pftcc_shsrch_minimize()
            self%proj_space(self%nrefs)%shvec = cxy(2:3)
            self%proj_space(self%nrefs)%inpl_corr = cxy(1)
        else
            self%proj_space(self%nrefs)%shvec = [0.,0.]
        endif
    end subroutine
    
    !>  \brief  calculates stochastic weights
    subroutine stochastic_weights( self, wcorr )
        class(prime3D_srch), intent(inout) :: self
        real, intent(out)                :: wcorr
        real                             :: wsum, spsum
        integer                          :: i, j, nbetter_inpl, alloc_stat, ref
        ! set unnormalised weights and calculate weight sum
        self%proj_space(:)%w = 0.
        wsum                 = 0.
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            if( ref > 0 )then
                if( self%proj_space(ref)%inpl_corr > 0. )then
                    self%proj_space(ref)%w = exp(self%proj_space(ref)%inpl_corr)
                    wsum = wsum+self%proj_space(ref)%w
                endif
            endif
        end do
        ! normalise weights & calculate weighted correlation
        wcorr = 0.     
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            if( ref > 0 )then
                if( self%proj_space(ref)%inpl_corr > 0. )then
                    self%proj_space(ref)%w = self%proj_space(ref)%w/wsum
                    wcorr = wcorr+self%proj_space(ref)%inpl_corr*self%proj_space(ref)%w
                endif
            endif
        end do
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH :: CALCULATED STOCHASTIC WEIGHTS PRIME3D_SRCH'        
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime3D_srch), intent(inout) :: self !< instance
        if( self%exists )then
            self%o_refs => null()
            self%pp     => null()
            deallocate(self%proj_space, self%proj_space_inds, self%refinds,&
            self%srch_order, self%state_range, self%proj_space_corrs, self%angtab )
            if( allocated(self%inplmat) )    deallocate(self%inplmat)
            if( allocated(self%prev_corrs) ) deallocate(self%prev_corrs)
            if( allocated(self%corrmat3d) )  deallocate(self%corrmat3d)
            if( allocated(self%corrmat2d) )  deallocate(self%corrmat2d)
            call self%rt%kill
            self%exists = .false.
        endif        
    end subroutine

end module simple_prime3D_srch
