! PRIME2D stochastic search routines

module simple_prime2D_srch
#include "simple_lib.f08"
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch,      only: pftcc_shsrch      ! simplex-based angle and shift search
use simple_pftcc_grad_shsrch, only: pftcc_grad_shsrch ! gradient-based angle and shift search
use simple_oris,              only: oris
use simple_timer              ! use all in there
implicit none

public :: prep4prime2D_srch, prime2D_srch
private

logical, parameter :: L_BENCH = .true.
logical, parameter :: DEBUG   = .false.

! allocatables for prime3D_srch are class variables to improve caching and reduce alloc overheads 
integer, allocatable :: cls_pops(:)
integer, allocatable :: srch_order(:,:)

type prime2D_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr => null()  !< pointer to pftcc (corrcalc) object
    class(oris),             pointer :: a_ptr     => null()  !< pointer to b%a (primary particle orientation table)
    type(pftcc_shsrch)               :: shsrch_obj           !< origin shift search object
    type(pftcc_grad_shsrch)          :: grad_shsrch_obj      !< origin shift search object, L-BFGS with gradient
    integer                          :: nrefs         =  0   !< number of references
    integer                          :: nrots         =  0   !< number of in-plane rotations in polar representation
    integer                          :: nrefs_eval    =  0   !< nr of references evaluated
    integer                          :: prev_class    =  0   !< previous class index
    integer                          :: best_class    =  0   !< best class index found by search
    integer                          :: prev_rot      =  0   !< previous in-plane rotation index
    integer                          :: best_rot      =  0   !< best in-plane rotation found by search
    integer                          :: nthr          =  0   !< number of threads
    integer                          :: fromp         =  1   !< from particle index
    integer                          :: top           =  1   !< to particle index
    integer                          :: nnn           =  0   !< # nearest neighbors
    integer                          :: iptcl         =  0   !< global particle index 
    real                             :: trs           =  0.  !< shift range parameter [-trs,trs]
    real                             :: prev_shvec(2) =  0.  !< previous origin shift vector
    real                             :: best_shvec(2) =  0.  !< best shift vector found by search
    real                             :: prev_corr     = -1.  !< previous best correlation
    real                             :: best_corr     = -1.  !< best corr found by search
    real                             :: specscore     =  0.  !< spectral score
    character(len=STDLEN)            :: refine        = ''   !< refinement flag
    character(len=STDLEN)            :: opt           = ''   !< optimizer flag
    ! timer vars
    real(timer_int_kind)             :: rt_refloop, rt_inpl, rt_tot
    integer(timer_int_kind)          :: t_refloop, t_inpl, t_tot
    ! logical flags
    logical                          :: dyncls  = .true.     !< whether to turn on dynamic class update (use of low population threshold)
    logical                          :: doshift = .true.     !< origin shift search indicator
  contains
    procedure          :: new
    procedure          :: exec_prime2D_srch
    procedure          :: stochastic_srch
    procedure          :: greedy_srch
    procedure          :: nn_srch
    procedure, private :: inpl_srch
    procedure          :: prep4srch
    procedure          :: update_best
    procedure          :: get_times
    procedure          :: kill
end type prime2D_srch

contains

    subroutine prep4prime2D_srch( b, p, which_iter )
        use simple_ran_tabu,  only: ran_tabu
        use simple_params,    only: params
        use simple_build,     only: build
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        integer,        intent(in)    :: which_iter
        type(ran_tabu) :: rt
        integer        :: iptcl, prev_class
        if( allocated(cls_pops)   ) deallocate(cls_pops)
        if( allocated(srch_order) ) deallocate(srch_order)
        ! gather class populations
        if( b%a%isthere('class') )then
            if( p%weights2D .eq. 'yes' .and. which_iter > 3 )then
                call b%a%get_pops(cls_pops, 'class', consider_w=.true., maxn=p%ncls)
            else
                call b%a%get_pops(cls_pops, 'class', consider_w=.false., maxn=p%ncls)
            endif
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(cls_pops(p%ncls), source=MINCLSPOPLIM+1, stat=alloc_stat)
            allocchk("simple_hadamard2D_matcher; prime2D_exec cls_pops")
        endif
        ! make random reference direction order
        allocate(srch_order(p%fromp:p%top,p%ncls), source=0)
        rt = ran_tabu(p%ncls)
        do iptcl = p%fromp, p%top
            call rt%reset
            call rt%ne_ran_iarr( srch_order(iptcl,:) )
            prev_class = nint(b%a%get(iptcl,'class'))
            call put_last(prev_class, srch_order(iptcl,:))
        end do
        if( any(srch_order == 0) ) stop 'Invalid index in srch_order; simple_prime2D_srch :: prep4prime2D_srch'
        call rt%kill
    end subroutine prep4prime2D_srch
    
    !>  \brief  is a constructor
    subroutine new( self, iptcl, pftcc, a, p )
        use simple_params, only: params
        class(prime2D_srch),             intent(inout) :: self   !< instance
        integer,                         intent(in)    :: iptcl  !< global particle index
        class(polarft_corrcalc), target, intent(inout) :: pftcc  !< correlator
        class(oris),             target, intent(in)    :: a      !< primary particle orientation table
        class(params),                   intent(in)    :: p      !< parameters
        integer :: alloc_stat, i
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%pftcc_ptr  => pftcc
        self%a_ptr      => a
        self%iptcl      =  iptcl
        self%nrefs      =  p%ncls
        self%nrots      =  round2even(twopi*real(p%ring2))
        self%nrefs_eval =  0
        self%trs        =  p%trs
        self%doshift    =  p%doshift
        self%refine     =  p%refine
        self%nthr       =  p%nthr
        self%fromp      =  p%fromp
        self%top        =  p%top
        self%nnn        =  p%nnn
        self%dyncls     =  (p%dyncls.eq.'yes')
        self%opt        =  p%opt
        ! construct composites
        lims(:,1)       = -p%trs
        lims(:,2)       =  p%trs
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        call self%shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init, nrestarts=3, maxits=60)
        call self%grad_shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init, maxits=60)
        if( L_BENCH )then
            ! init timers
            self%rt_refloop = 0.
            self%rt_inpl    = 0.
            self%rt_tot     = 0.  
        endif
        if( DEBUG ) print *, '>>> PRIME2D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME2D_SRCH OBJECT'
    end subroutine new

    !>  \brief the master prime search routine
    subroutine exec_prime2D_srch( self, greedy, extr_bound )
        class(prime2D_srch), intent(inout) :: self
        logical, optional,   intent(in)    :: greedy
        real,    optional,   intent(in)    :: extr_bound
        logical :: ggreedy
        if( self%refine .eq. 'yes' )then
            ggreedy = .true.
        else
            ggreedy = .false.
            if( present(greedy) ) ggreedy = greedy
        endif
        if( ggreedy )then
            call self%greedy_srch
        else
            call self%stochastic_srch(extr_bound)
        endif
        if( DEBUG ) print *, '>>> PRIME2D_SRCH::EXECUTED PRIME2D_SRCH'
    end subroutine exec_prime2D_srch

    !>  \brief  executes greedy rotational search (for initialisation)
    subroutine greedy_srch( self )
        class(prime2D_srch), intent(inout) :: self
        integer :: iref,loc(1),inpl_ind
        real    :: corrs(self%nrots),inpl_corr,corr
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            call self%prep4srch
            corr = self%prev_corr
            do iref=1,self%nrefs
                if( cls_pops(iref) == 0 )cycle
                call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs)
                loc       = maxloc(corrs)
                inpl_ind  = loc(1)
                inpl_corr = corrs(inpl_ind)         
                if( inpl_corr >= corr )then
                    corr            = inpl_corr
                    self%best_class = iref
                    self%best_corr  = inpl_corr
                    self%best_rot   = inpl_ind
                endif    
            end do
            self%nrefs_eval = self%nrefs
            call self%inpl_srch
            call self%update_best
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *, '>>> PRIME2D_SRCH::FINISHED STOCHASTIC SEARCH'
    end subroutine greedy_srch

    !>  \brief  executes stochastic rotational search
    subroutine stochastic_srch( self, extr_bound )
        class(prime2D_srch), intent(inout) :: self
        real, optional,      intent(in)    :: extr_bound
        integer :: iref, loc(1), isample, inpl_ind, nptcls, class_glob, inpl_glob
        real    :: corrs(self%nrots), inpl_corr, corr_bound, cc_glob
        logical :: found_better, do_inplsrch, glob_best_set
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            if( L_BENCH ) self%t_tot = tic()
            do_inplsrch   = .true.
            corr_bound    = -1.
            cc_glob       = -1.
            glob_best_set = .false. 
            if( present(extr_bound) ) corr_bound = extr_bound
            call self%prep4srch
            if( corr_bound < 0. .or. self%prev_corr > corr_bound )then
                ! SHC move
                found_better = .false.
                self%nrefs_eval = 0
                if( L_BENCH ) self%t_refloop = tic()
                do isample=1,self%nrefs
                    iref = srch_order(self%iptcl, isample)
                    ! keep track of how many references we are evaluating
                    self%nrefs_eval = self%nrefs_eval + 1
                    ! passes empty classes
                    if( cls_pops(iref) == 0 )cycle
                    ! shc update
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs)
                    inpl_ind  = shcloc(self%nrots, corrs, self%prev_corr)
                    if( inpl_ind == 0 )then
                        ! update inpl_ind & inpl_corr to greedy best
                        loc       = maxloc(corrs)
                        inpl_ind  = loc(1)
                        inpl_corr = corrs(inpl_ind)
                    else
                        ! use the parameters selected by SHC condition
                        inpl_corr       = corrs(inpl_ind)
                        self%best_class = iref
                        self%best_corr  = inpl_corr
                        self%best_rot   = inpl_ind
                        found_better    = .true.
                    endif
                    ! keep track of global best
                    if( inpl_corr > cc_glob )then
                        cc_glob       = inpl_corr
                        class_glob    = iref
                        inpl_glob     = inpl_ind
                        glob_best_set = .true.
                    endif
                    ! first improvement heuristic
                    if( found_better ) exit
                end do
                if( L_BENCH ) self%rt_refloop = self%rt_refloop + toc(self%t_refloop)
            else
                ! random move
                self%nrefs_eval = 1 ! evaluate one random ref
                isample = 1         ! random .ne. prev
                iref    = srch_order(self%iptcl, isample)
                if( self%dyncls )then
                    ! all good
                else
                    ! makes sure the ptcl does not land in an empty class
                    ! such that a search is performed
                    do while( cls_pops(iref) == 0 )
                        isample = isample + 1
                        iref    = srch_order(self%iptcl, isample)
                        if( isample.eq.self%nrefs )exit
                    enddo
                endif
                if( cls_pops(iref) == 0 )then
                    ! empty class
                    do_inplsrch = .false.               ! no in-plane search
                    inpl_ind    = irnd_uni(self%nrots)  ! random in-plane
                    nptcls      = self%a_ptr%get_noris()
                    inpl_corr   = -1.
                    do while( inpl_corr < TINY )
                        inpl_corr = self%a_ptr%get(irnd_uni(nptcls), 'corr') ! random correlation
                    enddo
                else
                    ! populated class
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs) 
                    loc       = maxloc(corrs)
                    inpl_ind  = loc(1)
                    inpl_corr = corrs(inpl_ind)
                endif
                self%best_class = iref
                self%best_corr  = inpl_corr
                self%best_rot   = inpl_ind
                found_better    = .true.
            endif
            if( found_better )then
                ! best ref has already been updated
            else
                if( glob_best_set )then
                    ! use the globally best parameters
                    self%best_class = class_glob
                    self%best_corr  = cc_glob
                    self%best_rot   = inpl_glob
                else
                    ! keep the old parameters
                    self%best_class = self%prev_class 
                    self%best_corr  = self%prev_corr
                    self%best_rot   = self%prev_rot
                endif
            endif
            if( do_inplsrch )then
                if( L_BENCH ) self%t_inpl = tic()
                call self%inpl_srch
                if( L_BENCH ) self%rt_inpl = self%rt_inpl + toc(self%t_inpl)
            endif
            if( .not. is_a_number(self%best_corr) )then
                print *, 'FLOATING POINT EXCEPTION ALARM; simple_prime2D_srch :: stochastic_srch'
                print *, self%iptcl, self%best_class, self%best_corr, self%best_rot
                print *, (corr_bound < 0. .or. self%prev_corr > corr_bound)
            endif
            call self%update_best
            if( L_BENCH ) self%rt_tot = self%rt_tot + toc(self%t_tot)
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *, '>>> PRIME2D_SRCH::FINISHED STOCHASTIC SEARCH'
    end subroutine stochastic_srch

    !>  \brief  executes nearest-neighbor rotational search
    subroutine nn_srch( self, nnmat )
        class(prime2D_srch), intent(inout) :: self
        integer,             intent(in)    :: nnmat(self%nrefs,self%nnn)
        integer :: iref,loc(1),inpl_ind,inn
        real    :: corrs(self%nrots),inpl_corr,corr
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            ! find previous discrete alignment parameters
            self%prev_class = nint(self%a_ptr%get(self%iptcl,'class'))                        ! class index
            self%prev_rot   = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl))     ! in-plane angle index
            self%prev_shvec = [self%a_ptr%get(self%iptcl,'x'),self%a_ptr%get(self%iptcl,'y')] ! shift vector
            ! set best to previous best by default
            self%best_class = self%prev_class         
            self%best_rot   = self%prev_rot
            ! calculate spectral score
            self%specscore  = self%pftcc_ptr%specscore(self%prev_class, self%iptcl, self%prev_rot)
            corr            = -1.
            ! evaluate neighbors (greedy selection)
            do inn=1,self%nnn
                iref      = nnmat(self%prev_class,inn)
                if( cls_pops(iref) == 0 )cycle
                call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs) 
                loc       = maxloc(corrs)
                inpl_ind  = loc(1)
                inpl_corr = corrs(inpl_ind)         
                if( inpl_corr >= corr )then
                    corr            = inpl_corr
                    self%best_class = iref
                    self%best_corr  = inpl_corr
                    self%best_rot   = inpl_ind
                endif    
            end do
            self%nrefs_eval = self%nrefs
            call self%inpl_srch
            call self%update_best
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *, '>>> PRIME2D_SRCH::FINISHED NEAREST-NEIGHBOR SEARCH'
    end subroutine nn_srch

    !>  \brief  executes the shift search over the best matching reference
    subroutine inpl_srch( self )
        class(prime2D_srch), intent(inout) :: self
        real, allocatable :: cxy(:)
        integer           :: irot
        self%best_shvec = [0.,0.]
        if( self%doshift )then
            select case(trim(self%opt))
                case('bfgs')
                    call self%grad_shsrch_obj%set_indices(self%best_class, self%iptcl)
                    cxy = self%grad_shsrch_obj%minimize(irot=irot)
                    if( irot > 0 )then
                        self%best_corr  = cxy(1)
                        self%best_rot   = irot
                        self%best_shvec = cxy(2:3)
                    endif
                case DEFAULT
                    call self%shsrch_obj%set_indices(self%best_class, self%iptcl)
                    cxy = self%shsrch_obj%minimize(irot=irot)
                    if( irot > 0 )then
                        self%best_corr  = cxy(1)
                        self%best_rot   = irot
                        self%best_shvec = cxy(2:3)
                    endif
            end select
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED SHIFT SEARCH'
    end subroutine inpl_srch

    !>  \brief  prepares for the search
    subroutine prep4srch( self )
        class(prime2D_srch), intent(inout) :: self
        real             :: corrs(self%pftcc_ptr%get_nrots())
        ! find previous discrete alignment parameters
        self%prev_class = nint(self%a_ptr%get(self%iptcl,'class')) ! class index
        if( self%dyncls )then
            if( self%prev_class > 0 )then
                ! reassignement to a class with higher population
                do while( cls_pops(self%prev_class) <= MINCLSPOPLIM )
                   self%prev_class = irnd_uni(self%nrefs)
                enddo
            endif
        endif
        self%prev_rot   = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl))     ! in-plane angle index
        self%prev_shvec = [self%a_ptr%get(self%iptcl,'x'),self%a_ptr%get(self%iptcl,'y')] ! shift vector
        ! set best to previous best by default
        self%best_class = self%prev_class
        self%best_rot   = self%prev_rot
        ! calculate previous best corr (treshold for better)
        if( self%prev_class > 0 )then
            call self%pftcc_ptr%gencorrs(self%prev_class, self%iptcl, corrs)
            self%prev_corr  = max(0., corrs(self%prev_rot))
            self%best_corr  = self%prev_corr
        else
            self%prev_class = irnd_uni(self%nrefs)
            self%prev_corr  = 0.
            self%best_corr  = 0.
            self%specscore  = 0.
        endif
        ! calculate spectral score
        self%specscore = self%pftcc_ptr%specscore(self%prev_class, self%iptcl, self%prev_rot)
        if( DEBUG ) print *, '>>> PRIME2D_SRCH::PREPARED FOR SIMPLE_PRIME2D_SRCH'
    end subroutine prep4srch

    !>  \brief  to get the class
    subroutine update_best( self )
        use simple_ori,  only: ori
        class(prime2D_srch), intent(in) :: self
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        real :: e3, mi_class, mi_inpl, mi_joint
        ! get in-plane angle
        e3   = 360. - self%pftcc_ptr%get_rot(self%best_rot) ! change sgn to fit convention
        ! calculate in-plane rot dist (radians)
        u(1) = 0.
        u(2) = 1.
        mat  = rotmat2d(e3)
        x1   = matmul(u,mat)
        mat  = rotmat2d(self%a_ptr%e3get(self%iptcl))
        x2   = matmul(u,mat)
        dist = myacos(dot_product(x1,x2))
        ! calculate overlap between distributions
        mi_class = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( self%prev_class == self%best_class )then
            mi_class = mi_class + 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%prev_rot ==  self%best_rot )then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif 
        mi_joint = mi_joint / 2.
        ! update parameters
        call self%a_ptr%e3set(self%iptcl,e3)
        call self%a_ptr%set_shift(self%iptcl, self%prev_shvec + self%best_shvec) 
        call self%a_ptr%set(self%iptcl, 'inpl',       real(self%best_rot))
        call self%a_ptr%set(self%iptcl, 'class',      real(self%best_class))
        call self%a_ptr%set(self%iptcl, 'corr',       self%best_corr)
        call self%a_ptr%set(self%iptcl, 'specscore',  self%specscore)
        call self%a_ptr%set(self%iptcl, 'dist_inpl',  rad2deg(dist))
        call self%a_ptr%set(self%iptcl, 'mi_class',   mi_class)
        call self%a_ptr%set(self%iptcl, 'mi_inpl',    mi_inpl)
        call self%a_ptr%set(self%iptcl, 'mi_joint',   mi_joint)
        call self%a_ptr%set(self%iptcl, 'frac',       100.*(real(self%nrefs_eval)/real(self%nrefs)))
        if( DEBUG ) print *, '>>> PRIME2D_SRCH::GOT BEST ORI'
    end subroutine update_best

    subroutine get_times( self, rt_refloop, rt_inpl, rt_tot )
        class(prime2D_srch),  intent(in) :: self
        real(timer_int_kind), intent(out):: rt_refloop, rt_inpl, rt_tot
        rt_refloop = self%rt_refloop
        rt_inpl    = self%rt_inpl
        rt_tot     = self%rt_tot
    end subroutine get_times

    subroutine kill( self )
        class(prime2D_srch),  intent(inout) :: self !< instance
        call self%grad_shsrch_obj%kill
    end subroutine kill

end module simple_prime2D_srch
