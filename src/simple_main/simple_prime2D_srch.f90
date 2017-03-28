module simple_prime2D_srch
use simple_defs              ! use all in there
use simple_math              ! use all in there
use simple_pftcc_shsrch      ! use all in there
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime_srch,       only: prime_srch
use simple_strings,          only: str_has_substr
implicit none

public :: prime2D_srch
private

logical, parameter :: DEBUG=.false.

type prime2D_srch
    private
    type(prime_srch)      :: srch_common          !< functionalities common to primesrch2D/3D
    integer               :: nrefs         = 0    !< number of references
    integer               :: nrots         = 0    !< number of in-plane rotations in polar representation
    integer               :: nrefs_eval    = 0    !< nr of references evaluated
    integer               :: prev_class    = 0    !< previous class index
    integer               :: best_class    = 0    !< best class index found by search
    integer               :: prev_rot      = 0    !< previous in-plane rotation index
    integer               :: best_rot      = 0    !< best in-plane rotation found by search
    integer               :: nthr          = 0    !< number of threads
    real                  :: trs           = 0.   !< shift range parameter [-trs,trs]
    real                  :: prev_shvec(2) = 0.   !< previous origin shift vector
    real                  :: best_shvec(2) = 0.   !< best ishift vector found by search
    real                  :: prev_corr     = -1.  !< previous best correlation
    real                  :: best_corr     = -1.  !< best corr found by search
    integer, allocatable  :: srch_order(:)        !< stochastic search order
    integer, allocatable  :: parts(:,:)           !< balanced partitions over references
    integer, allocatable  :: inplmat(:,:)         !< in-plane indices in matrix formulated search
    real,    allocatable  :: corrmat2d(:,:)       !< correlations in matrix formulated search
    character(len=STDLEN) :: refine               !< refinement flag
    logical               :: doshift = .true.     !< origin shift search indicator
    logical               :: exists  = .false.    !< 2 indicate existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS
    procedure :: get_nrots
    procedure :: get_corr
    procedure :: get_inpl
    procedure :: get_cls
    procedure :: get_roind
    ! PREPARATION ROUTINE
    procedure :: prep4srch
    ! SEARCH ROUTINES
    procedure :: exec_prime2D_srch
    procedure :: greedy_srch
    procedure :: stochastic_srch
    procedure :: shift_srch
    ! DESTRUCTOR
    procedure :: kill
end type prime2D_srch

contains

    ! CONSTRUCTOR
    
    !>  \brief  is a constructor
    subroutine new( self, p, ncls )
        use simple_params,     only: params
        use simple_map_reduce, only: split_nobjs_even
        class(prime2D_srch), intent(inout) :: self !< instance
        class(params),       intent(in)    :: p    !< parameters
        integer, optional,   intent(in)    :: ncls !< overriding dummy (ncls)
        integer :: alloc_stat, i
        real    :: dang
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        if( present(ncls) )then
            self%nrefs = ncls
        else
            self%nrefs = p%ncls
        endif
        self%nrots      = round2even(twopi*real(p%ring2))
        self%refine     = p%refine
        self%nrefs_eval = 0
        self%trs        = p%trs
        self%doshift    = p%doshift
        self%nthr       = p%nthr
        if( self%nrefs < self%nthr ) stop 'ncls < nthr not allowed; simple_prime2D_srch :: new'
        ! construct composites
        self%srch_common = prime_srch(p)
        ! find number of threads & create the same number of balanced partitions
        self%parts = split_nobjs_even(self%nrefs,self%nthr)
        ! the instance now exists
        self%exists = .true.
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::CONSTRUCTED NEW SIMPLE_prime2D_srch OBJECT'
    end subroutine new

    ! GETTERS

    !>  \brief nrots getter
    pure integer function get_nrots( self )
        class(prime2D_srch), intent(in) :: self
        get_nrots = self%nrots
    end function get_nrots

    !>  \brief correlation getter
    real function get_corr( self, iptcl, iref )
        class(prime2D_srch), intent(in) :: self
        integer,             intent(in) :: iptcl, iref
        get_corr = self%corrmat2d(iptcl, iref)
    end function get_corr
    
    !>  \brief in-plane index getter
    integer function get_inpl( self, iptcl, iref )
        class(prime2D_srch), intent(in) :: self
        integer,             intent(in) :: iptcl, iref
        get_inpl = self%inplmat(iptcl, iref)
    end function get_inpl
    
    !>  \brief  to get the class
    subroutine get_cls( self, o )
        use simple_math, only: myacos, rotmat2d, rad2deg
        use simple_ori,  only: ori
        class(prime2D_srch), intent(in)    :: self
        class(ori),          intent(inout) :: o
        real    :: euls(3), mi_class, mi_inpl, mi_joint
        integer :: class, rot
        real    :: x, y, mat(2,2), u(2), x1(2), x2(2)
        ! make unit vector
        u(1)     = 0.
        u(2)     = 1.
        ! calculate previous vec
        mat      = rotmat2d(o%e3get())
        x1       = matmul(u,mat)
        ! get new indices
        class    = self%best_class
        rot      = self%best_rot
        ! get in-plane angle
        euls     = 0.
        euls(3)  = 360.-self%srch_common%rot(rot) ! change sgn to fit convention
        if( euls(3) == 360. ) euls(3) = 0.
        call o%set_euler(euls)
        ! calculate new vec & distance (in degrees)
        mat      = rotmat2d(o%e3get())
        x2       = matmul(u,mat)
        ! calculate overlap between distributions
        mi_class = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( self%prev_class == class )then
            mi_class = mi_class + 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%prev_rot == rot )then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif 
        mi_joint = mi_joint / 2.
        ! set parameters
        x = self%prev_shvec(1)
        y = self%prev_shvec(2)
        if( self%doshift )then
            ! shifts must be obtained by vector addition
            x = x + self%best_shvec(1)
            y = y + self%best_shvec(2)
        endif
        call o%set('x',         x)
        call o%set('y',         y)
        call o%set('class',     real(class))
        call o%set('corr',      self%best_corr)
        call o%set('dist_inpl', rad2deg(myacos(dot_product(x1,x2))))
        call o%set('mi_class',  mi_class)
        call o%set('mi_inpl',   mi_inpl)
        call o%set('mi_joint',  mi_joint)
        call o%set('frac', 100.*(real(self%nrefs_eval)/real(self%nrefs)))
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::GOT BEST ORI'
    end subroutine get_cls

    !>  \brief returns the in-plane rotational index for the rot in-plane angle
    integer function get_roind( self, rot )
        class(prime2D_srch), intent(in) :: self
        real,                intent(in) :: rot
        get_roind = self%srch_common%roind(rot)
    end function get_roind

    ! PREPARATION ROUTINES

    !>  \brief  prepares for the search
    subroutine prep4srch( self, pftcc, iptcl, o_prev, corrmat )
        use simple_ori,      only: ori
        use simple_ran_tabu, only: ran_tabu
        use simple_rnd,      only: irnd_uni
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori), optional,    intent(inout) :: o_prev
        real,       optional,    intent(in)    :: corrmat(self%nrefs,self%nrots)
        type(ran_tabu) :: rt
        real           :: lims(2,2)
        ! initialize in-plane search classes
        lims(1,1) = -self%trs
        lims(1,2) =  self%trs
        lims(2,1) = -self%trs
        lims(2,2) =  self%trs
        call pftcc_shsrch_init(pftcc, lims)
        if( present(o_prev) )then
            ! find previous discrete alignment parameters
            self%prev_class = nint(o_prev%get('class'))                    ! class index
            self%prev_rot   = self%srch_common%roind(360.-o_prev%e3get())  ! in-plane angle index
            self%prev_shvec = [o_prev%get('x'),o_prev%get('y')]            ! shift vector
            ! set best to previous best by default
            self%best_class = self%prev_class         
            self%best_rot   = self%prev_rot
            if( present(corrmat) )then
                self%prev_corr = corrmat(self%prev_class,self%prev_rot)
                self%best_corr = corrmat(self%prev_class,self%prev_rot)
            else
                ! calculate previous best corr (treshold for better)
                self%prev_corr  = pftcc%corr(self%prev_class, iptcl, self%prev_rot)
                self%best_corr  = self%prev_corr
            endif
        else
            self%prev_class = irnd_uni(self%nrefs)
            self%prev_rot   = 1
            self%prev_shvec = 0.
            self%prev_corr  = 1.
        endif
        rt = ran_tabu(self%nrefs)
        if( allocated(self%srch_order) ) deallocate(self%srch_order)
        allocate(self%srch_order(self%nrefs))
        ! make random reference direction order
        call rt%ne_ran_iarr(self%srch_order)
        ! put prev_best last to avoid cycling
        call put_last(self%prev_class, self%srch_order)
        if( any(self%srch_order == 0) ) stop 'Invalid index in srch_order; simple_prime2D_srch :: prep4srch'
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::PREPARED FOR SIMPLE_prime2D_srch'
    end subroutine prep4srch

    ! SEARCH ROUTINES

    !>  \brief a master prime search routine
    subroutine exec_prime2D_srch( self, pftcc, a, pfromto, greedy)
        use simple_oris, only: oris
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(oris),             intent(inout) :: a
        integer,                 intent(in)    :: pfromto(2)
        logical, optional,       intent(in)    :: greedy
        real    :: lims(2,2)
        logical :: ggreedy
        ggreedy = .false.
        if( present(greedy) ) ggreedy = greedy
        if( self%refine .eq. 'greedy' .or. ggreedy )then
            call self%greedy_srch(pftcc, a, pfromto)
        else
            call self%stochastic_srch(pftcc, a, pfromto)
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::EXECUTED PRIME2D_SRCH'
    end subroutine exec_prime2D_srch

    !>  \brief  executes the greedy rotational search
    subroutine greedy_srch( self, pftcc, a, pfromto )
        use simple_oris, only: oris
        use simple_ori,  only:  ori
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(oris),             intent(inout) :: a
        integer,                 intent(in)    :: pfromto(2)
        type(ori)         :: orientation
        integer           :: iref, iptcl, loc(2)
        real              :: corrs(self%nrots)
        real, allocatable :: corrmat3d(:,:,:)
        ! calculate all corrs
        allocate(corrmat3d(pfromto(1):pfromto(2),self%nrefs,self%nrots))
        call pftcc%gencorrs_all_cpu(corrmat3d)
        ! search in-plane
        do iptcl=pfromto(1),pfromto(2)
            orientation = a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                ! initialize
                call self%prep4srch(pftcc, iptcl, orientation, corrmat=corrmat3d(iptcl,:,:))
                ! greedy selection
                loc = maxloc(corrmat3d(iptcl,:,:))
                self%best_class = loc(1)
                self%best_rot = loc(2)
                self%best_corr = corrmat3d(iptcl,loc(1),loc(2))
                ! search shifts
                call self%shift_srch(iptcl)
                ! we always evaluate all references using the greedy approach
                self%nrefs_eval = self%nrefs
                ! output info
                call self%get_cls(orientation)
                call a%set_ori(iptcl,orientation)
            else
                call orientation%reject
            endif
        end do
        deallocate(corrmat3d)
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED GREEDY SEARCH'
    end subroutine greedy_srch

    !>  \brief  executes the greedy rotational search
    subroutine stochastic_srch( self, pftcc, a, pfromto )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_oris, only: oris
        use simple_ori,  only: ori
        use simple_rnd,  only: shcloc
        use simple_syscalls
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(oris),             intent(inout) :: a
        integer,                 intent(in)    :: pfromto(2)
        type(ori)          :: orientation
        integer            :: iptcl,iref,loc(1),i,inpl_ind
        real               :: corrs(self%nrots),inpl_corr,corr
        logical            :: found_better,nnsrch
        real, allocatable  :: corrmat3d(:,:,:)
        ! calculate all corrs
        allocate(corrmat3d(pfromto(1):pfromto(2),self%nrefs,self%nrots))
        call pftcc%gencorrs_all_cpu(corrmat3d)
        ! execute search
        do iptcl=pfromto(1),pfromto(2)
            orientation = a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                ! initialize
                call self%prep4srch(pftcc, iptcl, orientation, corrmat=corrmat3d(iptcl,:,:))
                found_better = .false.
                self%nrefs_eval = 0
                do i=1,self%nrefs
                    iref = self%srch_order(i)
                    ! keep track of how many references we are evaluating
                    self%nrefs_eval = self%nrefs_eval + 1
                    corrs = corrmat3d(iptcl,iref,:)
                    inpl_ind = shcloc(self%nrots, corrs, self%prev_corr)
                    inpl_corr = 0.
                    if( inpl_ind > 0 ) inpl_corr = corrs(inpl_ind)             
                    if( inpl_ind > 0 .and. inpl_corr >= self%prev_corr )then
                        ! update the class
                        self%best_class = iref
                        ! update the correlation
                        self%best_corr = inpl_corr
                        ! update the in-plane angle
                        self%best_rot = inpl_ind
                        ! indicate that we found a better solution
                        found_better = .true.
                        exit ! first-improvement heuristic
                    endif    
                end do
                if( found_better )then
                    ! best ref has already been updated
                else
                    ! keep the old parameters
                    self%best_class = self%prev_class 
                    self%best_corr  = self%prev_corr
                    self%best_rot   = self%prev_rot
                endif
                ! search shifts
                call self%shift_srch(iptcl)
                ! output info
                call self%get_cls(orientation)
                call a%set_ori(iptcl,orientation)
            else
                call orientation%reject
            endif
        end do
        deallocate(corrmat3d)
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED STOCHASTIC SEARCH'
    end subroutine stochastic_srch

    !>  \brief  executes the in-plane search over one reference
    subroutine shift_srch( self, iptcl )
        class(prime2D_srch), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real :: cxy(3), corr_before, corr_after
        self%best_shvec = [0.,0.]
        if( self%doshift )then
            corr_before = self%best_corr
            call pftcc_shsrch_set_indices(self%best_class, iptcl, self%best_rot)
            cxy = pftcc_shsrch_minimize()
            corr_after = cxy(1)
            if( corr_after > corr_before )then
                self%best_corr  = cxy(1)
                self%best_shvec = cxy(2:3)
            endif
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED SHIFT SEARCH'
    end subroutine shift_srch

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime2D_srch), intent(inout) :: self !< instance
        if( self%exists )then
            if( allocated(self%corrmat2d)  ) deallocate(self%corrmat2d)
            if( allocated(self%inplmat)    ) deallocate(self%inplmat)
            if( allocated(self%srch_order) ) deallocate(self%srch_order)
            if( allocated(self%parts)      ) deallocate(self%parts)
            call self%srch_common%kill
            self%exists = .false.
        endif
    end subroutine kill

    !>  \brief a master prime search routine
    ! subroutine exec_prime2D_srch_old( self, pftcc, iptcl, o, nnmat )
    !     use simple_ori, only: ori
    !     class(prime2D_srch),     intent(inout) :: self
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     integer,                 intent(in)    :: iptcl
    !     class(ori), optional,    intent(inout) :: o
    !     integer,    optional,    intent(in)    :: nnmat(:,:)
    !     call self%prep4srch(pftcc, iptcl, o, nnmat=nnmat)
    !     if( .not.present(o) .or. self%refine.eq.'greedy' )then
    !         call self%greedy_srch_old(pftcc, iptcl)
    !     else
    !         call self%stochastic_srch_shc_old(pftcc, iptcl)
    !     endif
    !     if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::EXECUTED PRIME2D_SRCH'
    ! end subroutine exec_prime2D_srch_old

    !>  \brief  executes the greedy rotational search
    ! subroutine greedy_srch_old( self, pftcc, iptcl )
    !     class(prime2D_srch),     intent(inout) :: self
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     integer,                 intent(in)    :: iptcl
    !     integer :: ref, loc(1), inpl_ind, i, endit
    !     real    :: cxy(3), corrs(self%nrots), inpl_corr
    !     self%prev_corr = -1.
    !     endit = self%nrefs
    !     if( str_has_substr(self%refine,'neigh') ) endit = self%nnn
    !     do i=1,endit
    !         ref       = self%srch_order(i)
    !         corrs     = pftcc%gencorrs(ref, iptcl)
    !         loc       = maxloc(corrs)
    !         inpl_ind  = loc(1)
    !         inpl_corr = corrs(inpl_ind)
    !         if( inpl_corr >= self%prev_corr )then
    !             ! update the class
    !             self%best_class = ref
    !             ! update the correlations
    !             self%best_corr = inpl_corr
    !             self%prev_corr = self%best_corr
    !             ! update the in-plane angle
    !             self%best_rot = inpl_ind
    !         endif
    !     end do
    !     ! we always evaluate all references using the greedy approach
    !     self%nrefs_eval = self%nrefs
    !     if( str_has_substr(self%refine,'neigh') ) self%nrefs_eval = self%nnn
    !     ! search in-plane
    !     call self%shift_srch(iptcl)
    !     if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED GREEDY SEARCH'
    ! end subroutine greedy_srch_old

    !>  \brief  executes the stochastic rotational search
    ! subroutine stochastic_srch_shc_old( self, pftcc, iptcl )
    !     use simple_rnd, only: shcloc
    !     class(prime2D_srch),     intent(inout) :: self
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     integer,                 intent(in)    :: iptcl
    !     integer :: i, iref, inpl_ind, loc(1), endit
    !     real    :: corr_new, cxy(3), corrs(self%nrots), inpl_corr
    !     logical :: found_better
    !     ! initialize
    !     found_better    = .false.
    !     self%nrefs_eval = 0
    !     endit           = self%nrefs
    !     if( str_has_substr(self%refine,'neigh') ) endit = self%nnn
    !     ! search
    !     do i=1,endit
    !         iref = self%srch_order(i)
    !         ! keep track of how many references we are evaluating
    !         self%nrefs_eval = self%nrefs_eval + 1
    !         corrs           = pftcc%gencorrs(iref, iptcl)
    !         inpl_ind        = shcloc(self%nrots, corrs, self%prev_corr)
    !         inpl_corr       = 0.
    !         if( inpl_ind > 0 ) inpl_corr = corrs(inpl_ind)
    !         if( inpl_ind > 0 )then
    !             ! update the class
    !             self%best_class = iref
    !             ! update the correlation
    !             self%best_corr  = inpl_corr
    !             ! update the in-plane angle
    !             self%best_rot   = inpl_ind
    !             ! indicate that we found a better solution
    !             found_better    = .true.
    !             exit ! first-improvement heuristic
    !         endif    
    !     end do
    !     if( found_better )then
    !         ! best ref has already been updated
    !     else
    !         ! keep the old parameters
    !         self%best_class = self%prev_class 
    !         self%best_corr  = self%prev_corr
    !         self%best_rot   = self%prev_rot
    !     endif
    !     call self%shift_srch(iptcl)
    !     if( DEBUG ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED STOCHASTIC SEARCH'
    ! end subroutine stochastic_srch_shc_old

end module simple_prime2D_srch
