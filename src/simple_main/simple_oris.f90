!==Class simple_oris
!
! simple_oris handles orientation information. The code is distributed with the hope that it will be useful,
! but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2009-05-26
!
!==Changes and comments below
!* new class, HE 2012-03-27
!* could easily implement alignment document comparsions with this new class (convergence criterion)
!* all methods tested
!
module simple_oris
use simple_ori,    only: ori
use simple_math,   only: hpsort
use simple_defs    ! singleton
use simple_jiffys  ! singleton
implicit none

public :: oris, test_oris
private

!>  \brief  aggregates ori obejcts
type :: oris
    private
    type(ori), allocatable :: o(:)
    integer :: n=0
  contains
    ! CONSTRUCTORS
    procedure          :: new
    ! GETTERS
    procedure          :: e1get
    procedure          :: e2get
    procedure          :: e3get
    procedure          :: get_euler
    procedure          :: get_noris
    procedure          :: get_ori
    procedure          :: get
    procedure          :: get_all
    procedure          :: get_mat
    procedure          :: get_normal
    procedure, private :: isthere_1
    procedure, private :: isthere_2
    generic            :: isthere => isthere_1, isthere_2
    procedure          :: get_ncls
    procedure          :: get_pop
    procedure          :: get_clspop
    procedure          :: get_statepop
    procedure          :: get_ptcls_in_state
    procedure          :: get_nstates
    procedure          :: get_nlabels
    procedure          :: split_state
    procedure          :: split_class
    procedure          :: expand_classes
    procedure          :: remap_classes
    procedure          :: get_cls_pinds
    procedure          :: get_cls_oris
    procedure          :: get_state
    procedure          :: get_arr
    procedure, private :: calc_sum
    procedure          :: get_sum
    procedure          :: get_avg
    procedure, private :: calc_nonzero_sum
    procedure          :: get_nonzero_sum
    procedure          :: get_nonzero_avg
    procedure          :: get_ctfparams
    procedure          :: print
    procedure          :: print_mats
    ! procedure          :: distance_transform
    ! SETTERS
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: set_euler
    procedure          :: e1set
    procedure          :: e2set
    procedure          :: e3set
    procedure          :: set
    procedure          :: set_ori
    procedure          :: set_all
    procedure          :: set_all2single
    procedure          :: e3swapsgn
    procedure          :: swape1e3
    procedure          :: zero
    procedure          :: zero_projs
    procedure          :: zero_shifts
    procedure          :: mul_shifts
    procedure          :: rnd_oris
    procedure          :: rnd_proj_space
    procedure          :: rnd_neighbors
    procedure          :: rnd_ori
    procedure          :: rnd_cls
    procedure          :: rnd_oris_discrete
    procedure          :: rnd_inpls
    procedure          :: rnd_trs
    procedure          :: rnd_ctf
    procedure          :: rnd_states
    procedure          :: rnd_classes
    procedure          :: rnd_lps
    procedure          :: rnd_corrs
    procedure          :: revshsgn
    procedure          :: revorisgn
    procedure          :: merge_classes
    procedure          :: symmetrize
    procedure          :: merge
    procedure          :: merge_files
    ! I/O
    procedure          :: read
    procedure, private :: write_1
    procedure, private :: write_2
    generic            :: write => write_1, write_2
    ! CALCULATORS
    procedure, private :: homogeneity_1
    procedure, private :: homogeneity_2
    generic            :: homogeneity => homogeneity_1, homogeneity_2
    procedure          :: cohesion_norm
    procedure          :: cohesion
    procedure          :: cohesion_ideal
    procedure          :: separation_norm
    procedure          :: separation
    procedure          :: separation_ideal
    procedure          :: histogram
    procedure          :: round_shifts
    procedure          :: introd_alig_err
    procedure          :: introd_ctf_err
    procedure, private :: rot_1
    procedure, private :: rot_2
    generic            :: rot => rot_1, rot_2
    procedure, private :: median_1
    generic            :: median => median_1
    procedure          :: stats
    procedure          :: minmax
    procedure, private :: spiral_1
    procedure, private :: spiral_2
    generic            :: spiral => spiral_1, spiral_2
    procedure          :: qspiral
    procedure          :: order
    procedure          :: order_cls
    procedure          :: calc_hard_ptcl_weights
    procedure          :: find_closest_proj
    procedure          :: find_closest_projs
    procedure          :: find_closest_ori
    procedure          :: find_closest_oris
    procedure          :: discretize
    procedure          :: nearest_neighbors
    procedure          :: find_angres
    procedure          :: find_angres_geod
    procedure          :: find_npeaks
    procedure          :: class_calc_frac
    procedure          :: class_dist_stat
    procedure          :: class_corravg
    procedure          :: best_in_class
    procedure, private :: calc_avg_class_shifts_1
    procedure, private :: calc_avg_class_shifts_2
    generic            :: calc_avg_class_shifts => calc_avg_class_shifts_1, calc_avg_class_shifts_2
    procedure, private :: map3dshift22d_1
    procedure, private :: map3dshift22d_2
    generic            :: map3dshift22d => map3dshift22d_1, map3dshift22d_2
    procedure          :: add_shift2class
    procedure, nopass  :: corr_oris
    procedure          :: gen_diverse
    procedure          :: gen_diversity_score
    procedure          :: gen_diversity_scores
    procedure          :: ori_generator
    procedure          :: gen_smat
    procedure          :: geodesic_dist
    procedure          :: gen_subset
    procedure, private :: diststat_1
    procedure, private :: diststat_2
    generic            :: diststat => diststat_1, diststat_2
    procedure          :: cluster_diststat
    procedure          :: mirror3d
    procedure          :: mirror2d
    procedure          :: cls_overlap
    ! DESTRUCTOR
    procedure          :: kill
end type

interface oris
    module procedure constructor
end interface

type(ori), pointer   :: op(:)=>null()
type(oris), pointer  :: ops  =>null()
integer, allocatable :: classpops(:)

contains

    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    function constructor( n ) result( self )
        integer, intent(in) :: n
        type(oris) :: self
        call self%new(n)
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, n )
        class(oris), intent(inout) :: self
        integer, intent(in) :: n
        integer :: alloc_stat, i
        call self%kill
        self%n = n
        allocate( self%o(self%n), stat=alloc_stat )
        call alloc_err('new; simple_oris', alloc_stat)
        do i=1,n
            call self%o(i)%new
        end do
    end subroutine new
    
    ! GETTERS
    
    !>  \brief  returns the i:th first Euler angle 
    pure function e1get( self, i ) result( e1 )
        class(oris), intent(in) :: self
        integer, intent(in)     :: i
        real                    :: e1
        e1 = self%o(i)%e1get()
    end function e1get
    
    !>  \brief  returns the i:th second Euler angle 
    pure function e2get( self, i ) result( e2 )
        class(oris), intent(in) :: self
        integer, intent(in)     :: i
        real                    :: e2
        e2 = self%o(i)%e2get()
    end function e2get
    
    !>  \brief  returns the i:th third Euler angle 
    pure function e3get( self, i ) result( e3 )
        class(oris), intent(in) :: self
        integer, intent(in)     :: i
        real                    :: e3
        e3 = self%o(i)%e3get()
    end function e3get
    
    !>  \brief  returns the i:th Euler triplet
    pure function get_euler( self, i ) result( euls )
        class(oris), intent(in) :: self
        integer, intent(in)     :: i
        real                    :: euls(3)
        euls = self%o(i)%get_euler()
    end function get_euler
    
    !>  \brief  is for getting the number of orientations
    pure function get_noris( self ) result( n )
        class(oris), intent(in) :: self
        integer :: n
        n = self%n
    end function get_noris
    
    !>  \brief  returns an individual orientation
    function get_ori( self, i ) result( o )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: i
        type(ori)                  :: o
        if( self%n == 0 ) stop 'oris object does not exist; get_ori; simple_oris'
        if( i > self%n .or. i < 1 )then
            write(*,*) 'trying to get ori: ', i, ' among: ', self%n, ' oris'
            stop 'i out of range; get_ori; simple_oris'
        endif
        o = self%o(i)
    end function get_ori
    
    !>  \brief  is a multiparameter getter
    function get( self, i, key ) result( val ) 
        class(oris), intent(inout)     :: self
        integer, intent(in)            :: i
        character(len=*), intent(in)   :: key
        real                           :: val
        val = self%o(i)%get(key)
    end function get

    !>  \brief  is for getting an array of 'key' values
    function get_all( self, key ) result( arr ) 
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        real, allocatable :: arr(:)
        integer           :: i, alloc_stat
        allocate( arr(self%n), stat=alloc_stat)
        call alloc_err('get_all; simple_oris', alloc_stat)
        do i=1,self%n
            arr(i) = self%o(i)%get(key)
        enddo
    end function get_all
    
    !>  \brief  is for getting the i:th rotation matrix
    pure function get_mat( self, i ) result( mat )
        class(oris), intent(in) :: self
        integer, intent(in)     :: i
        real :: mat(3,3)
        mat = self%o(i)%get_mat()
    end function get_mat
    
    !>  \brief  is for getting the i:th normal
    pure function get_normal( self, i ) result( normal )
        class(oris), intent(in) :: self
        integer, intent(in)     :: i
        real :: normal(3)
        normal = self%o(i)%get_normal()
    end function get_normal
    
    !>  \brief  is for checking if parameter is present
    function isthere_1( self, key ) result( is )
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: key
        logical                      :: is
        integer                      :: i
        is = .false.
        do i=1,self%n
            is = self%o(i)%isthere(key)
            if( is ) exit
        end do
    end function isthere_1
    
    !>  \brief  is for checking if parameter is present
    function isthere_2( self, i, key ) result( is )
        class(oris), intent(inout)   :: self
        integer, intent(in)          :: i 
        character(len=*), intent(in) :: key
        logical                      :: is
        is = self%o(i)%isthere(key)
    end function isthere_2
    
    !>  \brief  is for checking the number of classes
    function get_ncls( self, minpop ) result( ncls )
        class(oris), intent(inout)    :: self
        integer, intent(in), optional :: minpop
        integer :: i, ncls, mycls, ncls_here
        ncls = 1
        do i=1,self%n
            mycls = nint(self%o(i)%get('class'))
            if( mycls > ncls ) ncls = mycls
        end do
        if( present(minpop) )then
            ncls_here = 0
            do i=1,ncls
                if( self%get_clspop(i) >= minpop ) ncls_here = ncls_here+1
            end do
            ncls = ncls_here
        endif
    end function get_ncls
    
    !>  \brief  is for checking label population
    function get_pop( self, ind, label ) result( pop )
        class(oris), intent(inout)   :: self
        integer, intent(in)          :: ind
        character(len=*), intent(in) :: label
        integer :: mycls, pop, i, mystate
        pop = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate > 0 )then
                mycls = nint(self%o(i)%get(label))
                if( mycls == ind )then
                    pop = pop+1
                endif
            endif
        end do
    end function get_pop
    
    !>  \brief  is for checking class population
    function get_clspop( self, class ) result( pop )
        class(oris), intent(inout) :: self
        integer, intent(in) :: class
        integer :: mycls, pop, i, mystate
        pop = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate > 0 )then
                mycls   = nint(self%o(i)%get('class'))
                if( mycls == class )then
                    pop = pop+1
                endif
            endif
        end do
    end function get_clspop
    
    !>  \brief  is for checking state population
    function get_statepop( self, state ) result( pop )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: state
        integer :: mystate, pop, i
        pop = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate == state ) pop = pop+1
        end do
    end function get_statepop
    
    !>  \brief  4 getting the particle indices corresponding to state
    function get_ptcls_in_state( self, state ) result( ptcls )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: state
        real, allocatable          :: ptcls(:)
        integer                    :: pop, mystate, cnt, alloc_stat, i
        pop = self%get_statepop(state)
        allocate( ptcls(pop), stat=alloc_stat )
        call alloc_err("In: get_ptcls_in_state; simple_oris", alloc_stat)
        cnt = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate == state )then
                cnt = cnt+1
                ptcls(cnt) = i
            endif
        end do
    end function get_ptcls_in_state
    
    !>  \brief  is for getting the number of states
    function get_nstates( self ) result( nstates )
        class(oris), intent(inout) :: self
        integer                    :: nstates, i, mystate
        nstates = 1
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate > nstates ) nstates = mystate
        end do
    end function get_nstates
    
    !>  \brief  is for getting the number of labels
    function get_nlabels( self ) result( nlabels )
        class(oris), intent(inout) :: self
        integer                    :: nlabels, i, mylabel
        nlabels = 1
        do i=1,self%n
            mylabel = nint(self%o(i)%get('label'))
            if( mylabel > nlabels ) nlabels = mylabel
        end do
    end function get_nlabels
    
    !>  \brief  for balanced split of a state group
    subroutine split_state( self, which )
        use simple_ran_tabu, only: ran_tabu
        class(oris), intent(inout) :: self
        integer, intent(in)        :: which
        integer, allocatable :: ptcls_in_which(:)
        integer, allocatable :: states(:)
        type(ran_tabu)       :: rt
        integer              :: alloc_stat, n, nstates, iptcl
        nstates = self%get_nstates()
        if( which < 1 .or. which > nstates )then
            stop 'which (state) is out of range; simple_oris::split_state'
        endif
        ptcls_in_which = self%get_state(which)
        n = size(ptcls_in_which)
        allocate(states(n), stat=alloc_stat)
        call alloc_err("simple_oris::split_state", alloc_stat)
        rt = ran_tabu(n)
        call rt%balanced(2, states)
        do iptcl=1,n
            if( states(iptcl) == 1 )then
                ! do nothing, leave this state as is
            else
                call self%o(ptcls_in_which(iptcl))%set('state', real(nstates+1))
            endif
        end do
        call rt%kill
        deallocate(ptcls_in_which,states)
    end subroutine split_state
    
    !>  \brief  for balanced split of a class
    subroutine split_class( self, which )
        use simple_ran_tabu, only: ran_tabu
        class(oris), intent(inout) :: self
        integer, intent(in)        :: which
        integer, allocatable :: ptcls_in_which(:)
        integer, allocatable :: members(:)
        type(ran_tabu)       :: rt
        integer              :: alloc_stat, n, nmembers, iptcl
        nmembers = self%get_ncls()
        if( which < 1 .or. which > nmembers )then
            stop 'which (member is out of range; simple_oris::split_class'
        endif
        ptcls_in_which = self%get_cls_pinds(which)
        n = size(ptcls_in_which)
        allocate(members(n), stat=alloc_stat)
        call alloc_err("simple_oris::split_class", alloc_stat)
        rt = ran_tabu(n)
        call rt%balanced(2, members)
        do iptcl=1,n
            if( members(iptcl) == 1 )then
                ! do nothing, leave this state as is
            else
                call self%o(ptcls_in_which(iptcl))%set('class', real(nmembers+1))
            endif
        end do
        call rt%kill
        deallocate(ptcls_in_which,members)
    end subroutine split_class
    
    !>  \brief  for expanding the number of classes using balanced splitting 
    subroutine expand_classes( self, ncls_target )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: ncls_target
        integer, allocatable       :: pops(:)
        integer :: ncls, loc(1), myncls, icls
        ncls = self%get_ncls()
        if( ncls_target <= ncls ) stop 'Number of target classes cannot be <= original number'
        ! calculate class populations
        allocate(pops(ncls_target))
        pops = 0
        do icls=1,ncls
            pops(icls) = self%get_clspop(icls)
        end do
        myncls = ncls
        do while( myncls < ncls_target )
            ! split largest
            loc = maxloc(pops)
            call self%split_class(loc(1))
            ! update number of classes
            myncls = myncls+1
            ! update pops
            pops(loc(1)) = self%get_clspop(loc(1))
            pops(myncls) = self%get_clspop(myncls)
        end do
    end subroutine expand_classes
    
    !>  \brief  for remapping classes after exclusion
    subroutine remap_classes( self )
        class(oris), intent(inout) :: self
        integer :: ncls, clsind_remap, pop, icls, iptcl, old_cls
        integer , allocatable :: clspops(:)
        ncls = self%get_ncls()
        allocate(clspops(ncls))
        do icls=1,ncls
            clspops(icls) = self%get_clspop(icls)
        end do
        if( any(clspops == 0) )then
            clsind_remap = ncls
            do icls=1,ncls
                pop = clspops(icls)
                if( pop > 1 )then
                    clsind_remap = clsind_remap+1
                    do iptcl=1,self%n
                        old_cls = nint(self%o(iptcl)%get('class'))
                        if( old_cls == icls ) call self%o(iptcl)%set('class', real(clsind_remap))
                    end do
                else
                    do iptcl=1,self%n
                        old_cls = nint(self%o(iptcl)%get('class'))
                        if( old_cls == icls )then
                            call self%o(iptcl)%set('class', 0.)
                            call self%o(iptcl)%set('state', 0.)
                        endif
                    end do
                endif
            end do
            do iptcl=1,self%n
                old_cls = nint(self%o(iptcl)%get('class'))
                if( old_cls /= 0 ) call self%o(iptcl)%set('class', real(old_cls-ncls))
            end do
        endif
        deallocate(clspops)
    end subroutine remap_classes
    
    !>  \brief  is for getting an allocatable array with ptcl indices of the class 'class'
    function get_cls_pinds( self, class ) result( clsarr )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: class
        integer, allocatable       :: clsarr(:)
        integer                    :: clsnr, i, pop, alloc_stat, cnt
        pop = self%get_clspop( class )
        if( pop > 0 )then
            allocate( clsarr(pop), stat=alloc_stat )
            call alloc_err('get_cls; simple_oris', alloc_stat)
            cnt = 0
            do i=1,self%n
                clsnr = nint(self%get( i, 'class'))
                if( clsnr == class )then
                    cnt = cnt+1
                    clsarr(cnt) = i
                endif
            end do
        endif
    end function get_cls_pinds

    !>  \brief  is for extracting the subset of oris with class label class
    function get_cls_oris( self, class ) result( clsoris )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: class
        type(oris)                 :: clsoris
        integer                    :: clsnr, i, pop, cnt
        pop = self%get_clspop( class )
        if( pop > 0 )then
            call clsoris%new(pop)
            cnt = 0
            do i=1,self%n
                clsnr = nint(self%get(i, 'class'))
                if( clsnr == class )then
                    cnt = cnt+1
                    clsoris%o(cnt) = self%o(i)
                endif
            end do
        endif
    end function get_cls_oris
    
    !>  \brief  is for getting an allocatable array with ptcl indices of the stategroup 'state'
    function get_state( self, state ) result( statearr )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: state
        integer, allocatable       :: statearr(:)
        integer                    :: statenr, i, pop, alloc_stat, cnt
        pop = self%get_statepop( state )
        if( pop > 0 )then
            allocate( statearr(pop), stat=alloc_stat )
            call alloc_err('get_state; simple_oris', alloc_stat)
            cnt = 0
            do i=1,self%n
                statenr = nint(self%get( i, 'state' ))
                if( statenr == state )then
                    cnt = cnt+1
                    statearr(cnt) = i
                endif
            end do
        endif
    end function get_state
    
    !>  \brief  is for getting an array of 'which' variables with 
    !!          filtering based on class/state
    function get_arr( self, which, class, state ) result( vals )
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which 
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        real, allocatable             :: vals(:)
        integer                       :: pop, cnt, clsnr, i, alloc_stat, mystate
        real                          :: val
        if( present(class) )then
            pop = self%get_clspop(class)
        else if( present(state) )then
            pop = self%get_statepop(state)
        else
            pop = self%n
        endif
        if( pop > 0 )then
            allocate( vals(pop), stat=alloc_stat )
            call alloc_err('get_arr; simple_oris', alloc_stat)
            cnt = 0
            do i=1,self%n
                val = self%get(i, which)
                if( present(class) )then
                    mystate = nint(self%get( i, 'state'))
                    if( mystate > 0 )then
                        clsnr   = nint(self%get( i, 'class'))
                        if( clsnr == class )then
                            cnt = cnt+1
                            vals(cnt) = val
                        endif
                    endif
                else if( present(state) )then
                    mystate = nint(self%get( i, 'state'))
                    if( mystate == state )then
                        cnt = cnt+1
                        vals(cnt) = val
                    endif
                else
                    vals(i) = val
                endif
            end do
        endif
    end function get_arr
    
    !>  \brief  is for calculating the sum of 'which' variables with 
    !!          filtering based on class/state/fromto
    subroutine calc_sum( self, which, sum, cnt, class, state, fromto, mask )
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which
        real, intent(out)             :: sum
        integer, intent(out)          :: cnt
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        integer, intent(in), optional :: fromto(2)
        logical, intent(in), optional :: mask(self%n)
        integer                       :: clsnr, i, mystate, istart, istop
        real                          :: val
        logical                       :: proceed
        cnt = 0
        sum = 0.
        if( present(fromto) )then
            istart = fromto(1)
            istop  = fromto(2)
        else
            istart = 1
            istop  = self%n
        endif
        do i=istart,istop
            mystate = nint(self%get( i, 'state'))
            if( mystate == 0 ) cycle
            if( present(mask) )then
                proceed = mask(i)
            else    
                proceed = .true.
            endif
            if( proceed )then
                val = self%get(i, which)
                if( isnan(val) ) val=0.
                if( present(class) )then
                    clsnr = nint(self%get( i, 'class'))
                    if( clsnr == class )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else if( present(state) )then
                    if( mystate == state )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else
                    cnt = cnt+1
                    sum = sum+val
                endif
            endif
        end do
    end subroutine calc_sum
    
    !>  \brief  is for getting the sum of 'which' variables with 
    !!          filtering based on class/state/fromto
    function get_sum( self, which, class, state, fromto, mask) result( sum )
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which 
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        integer, intent(in), optional :: fromto(2)
        logical, intent(in), optional :: mask(self%n)
        integer                       :: cnt
        real                          :: sum
        call self%calc_sum(which, sum, cnt, class, state, fromto, mask)
    end function get_sum
    
    !>  \brief  is for getting the average of 'which' variables with 
    !!          filtering based on class/state/fromto
    function get_avg( self, which, class, state, fromto, mask ) result( avg )
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which 
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        integer, intent(in), optional :: fromto(2)
        logical, intent(in), optional :: mask(self%n)
        integer                       :: cnt
        real                          :: avg, sum
        call self%calc_sum(which, sum, cnt, class, state, fromto, mask)
        avg = sum/real(cnt)
    end function get_avg
    
    !>  \brief  is for calculating the nonzero sum of 'which' variables with 
    !!          filtering based on class/state/fromto
    subroutine calc_nonzero_sum( self, which, sum, cnt, class, state, fromto )
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which
        real, intent(out)             :: sum
        integer, intent(out)          :: cnt
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        integer, intent(in), optional :: fromto(2)
        integer                       :: clsnr, i, mystate, istart, istop
        real                          :: val
        cnt = 0
        sum = 0.
        if( present(fromto) )then
            istart = fromto(1)
            istop  = fromto(2)
        else
            istart = 1
            istop  = self%n
        endif
        do i=istart,istop
            mystate = nint(self%get( i, 'state'))
            if( mystate == 0 ) cycle
            val = self%get(i, which)
            if( isnan(val) ) val=0.
            if( val > 0. )then
                if( present(class) )then
                    clsnr = nint(self%get( i, 'class'))
                    if( clsnr == class )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else if( present(state) )then
                    
                    if( mystate == state )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else
                    cnt = cnt+1
                    sum = sum+val
                endif
            endif
        end do
    end subroutine calc_nonzero_sum
    
    !>  \brief  is for getting the sum of 'which' variables with 
    !!          filtering based on class/state/fromto
    function get_nonzero_sum( self, which, class, state, fromto ) result( sum )
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which 
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        integer, intent(in), optional :: fromto(2)
        integer                       :: cnt
        real                          :: sum
        call self%calc_nonzero_sum(which, sum, cnt, class, state, fromto)
    end function get_nonzero_sum
    
    !>  \brief  is for getting the average of 'which' variables with 
    !!          filtering based on class/state/fromto
    function get_nonzero_avg( self, which, class, state, fromto ) result( avg )
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which 
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        integer, intent(in), optional :: fromto(2)
        integer                       :: cnt
        real                          :: avg, sum
        call self%calc_nonzero_sum(which, sum, cnt, class, state, fromto)
        avg = sum/real(cnt)
    end function get_nonzero_avg

    !>  \brief  is for getting a matrix of CTF parameters
    function get_ctfparams( self, pfromto ) result( ctfparams )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: pfromto(2)
        real, allocatable :: ctfparams(:,:)
        integer :: nparams, cnt, iptcl
        logical :: astig
        nparams = pfromto(2)-pfromto(1)+1
        allocate(ctfparams(nparams,6))
        ctfparams = 0.
        astig = .false.
        if( self%isthere('dfx') )then
            ! ok
            if( self%isthere('dfy') ) astig = .true.
        else
            stop 'CTF parameters not in oris; simple_oris :: get_ctfparams'
        endif
        cnt = 0
        do iptcl=pfromto(1),pfromto(2)
            cnt = cnt + 1
            ctfparams(cnt,1) = self%get(iptcl,'kv')
            ctfparams(cnt,2) = self%get(iptcl,'cs')
            ctfparams(cnt,3) = self%get(iptcl,'fraca')
            ctfparams(cnt,4) = self%get(iptcl,'dfx')
            if( astig )then
                ctfparams(cnt,5) = self%get(iptcl,'dfy')
                ctfparams(cnt,6) = self%get(iptcl,'angast')
            else
                ctfparams(cnt,5) = self%get(iptcl,'dfx')
                ctfparams(cnt,6) = 0.
            endif
        end do
    end function get_ctfparams
    
    !>  \brief  is for printing
    subroutine print( self, i )
        class(oris), intent(inout) :: self
        integer, intent(in) :: i
        call self%o(i)%print
    end subroutine print
    
    !>  \brief  is for printing
    subroutine print_mats( self )
        class(oris), intent(inout) :: self
        integer                    :: i
        do i=1,self%n
            write(*,*) i
            call self%o(i)%print_mat
        end do
    end subroutine print_mats
    
    ! SETTERS
    
    !>  \brief is polymorphic assignment, overloaded as (=)
    subroutine assign( self_out, self_in )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        integer   :: i
        call self_out%new(self_in%n)
        do i=1,self_in%n
            self_out%o(i) = self_in%o(i)
        end do
    end subroutine assign
    
    !>  \brief  is a setter
    subroutine set_euler( self, i, euls )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: i
        real, intent(in)           :: euls(3)
        call self%o(i)%set_euler(euls)
    end subroutine set_euler
    
    !>  \brief  is a setter
    subroutine e1set( self, i, e1 )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: i
        real, intent(in)           :: e1
        call self%o(i)%e1set(e1)
    end subroutine e1set
    
    !>  \brief  is a setter
    subroutine e2set( self, i, e2 )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: i
        real, intent(in)           :: e2
        call self%o(i)%e2set(e2)
    end subroutine e2set
    
    !>  \brief  is a setter 
    subroutine e3set( self, i, e3 )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: i
        real, intent(in)           :: e3
        call self%o(i)%e3set(e3)
    end subroutine e3set
    
    !>  \brief  is a setter
    subroutine set( self, i, key, val )
        class(oris), intent(inout)   :: self
        integer, intent(in)          :: i
        character(len=*), intent(in) :: key
        real, intent(in)             :: val
        call self%o(i)%set(key, val)
    end subroutine set
    
    !>  \brief  is a setter
    subroutine set_ori( self, i, o )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: i
        class(ori), intent(in)     :: o
        self%o(i) = o
    end subroutine set_ori
    
    !>  \brief  is a setter
    subroutine set_all( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all

     !>  \brief  is a setter
    subroutine set_all2single( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, val)
        end do
    end subroutine set_all2single

    !>  \brief  is a setter
    subroutine e3swapsgn( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%e3set(i,360.-self%e3get(i))
        end do
    end subroutine e3swapsgn
    
    !>  \brief  is a setter
    subroutine swape1e3( self )
        class(oris), intent(inout) :: self
        integer :: i
        real :: e, euls(3)
        do i=1,self%n
            euls = self%get_euler(i)
            e = euls(1)
            euls(1) = euls(3)
            euls(3) = e
            call self%set_euler(i,euls)
        end do
    end subroutine swape1e3
    
    !>  \brief  zero the 'which' var
    subroutine zero( self, which )
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, 0.)
        end do
    end subroutine zero
    
    !>  \brief  zero the projection directions
    subroutine zero_projs( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%e1set(0.)
            call self%o(i)%e2set(0.)
        end do
    end subroutine zero_projs
    
    !>  \brief  zero the shifts
    subroutine zero_shifts( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', 0.)
            call self%o(i)%set('y', 0.)
        end do
    end subroutine zero_shifts
    
    !>  \brief  zero the shifts
    subroutine mul_shifts( self, mul )
        class(oris), intent(inout) :: self
        real, intent(in)           :: mul
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', mul*self%o(i)%get('x'))
            call self%o(i)%set('y', mul*self%o(i)%get('y'))
        end do
    end subroutine mul_shifts

    !>  \brief  randomizes eulers in oris
    subroutine rnd_oris( self, trs, eullims)
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
        integer :: i
        do i=1,self%n
            call self%rnd_ori(i, trs, eullims)
        end do 
    end subroutine rnd_oris
    
    !>  \brief  generates nnn stochastic neighbors to o_prev with angular threshold athres
    !!          optional proj indicates projection direction only threshold or not
    subroutine rnd_neighbors( self, nnn, o_prev, athres, proj )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: nnn
        class(ori),        intent(in)    :: o_prev
        real,              intent(in)    :: athres
        logical, optional, intent(in)    :: proj
        integer :: i
        call self%new(nnn)
        do i=1,nnn
            call self%o(i)%rnd_euler(o_prev, athres, proj)
            if( present(proj) )then
                if( proj ) call self%e3set(i,0.)
            endif
        end do
    end subroutine rnd_neighbors

    !>  \brief  generate random projection direction space around a given one
    subroutine rnd_proj_space( self, nsample, o_prev, thres, eullims )
        use simple_math, only: rad2deg
        class(oris), intent(inout)           :: self
        integer,     intent(in)              :: nsample
        class(ori),  intent(inout), optional :: o_prev
        real,        intent(inout), optional :: eullims(3,2)
        real,        intent(in),    optional :: thres
        type(ori) :: o_stoch
        integer   :: i
        logical   :: within_lims, found
        if( present(o_prev) .and. .not.present(thres) ) &
            & stop 'Missing angular threshold in simple_oris::rnd_proj_space'
        if( .not.present(o_prev) .and. present(thres) ) &
            & stop 'Missing orientation in simple_oris::rnd_proj_space'
        within_lims = .false.
        if( present(eullims) )within_lims = .true.
        call self%new(nsample)
        if( present(o_prev).and.present(thres) )then
            do i=1,self%n
                found = .false.
                do while( .not.found )
                    call o_stoch%new
                    if( within_lims )then
                        call o_stoch%rnd_euler( eullims )
                    else
                        call o_stoch%rnd_euler
                    endif
                    if( rad2deg( o_stoch.euldist.o_prev ) > 2.*thres )cycle
                    found = .true.
                    call o_stoch%e3set( 0.)
                    self%o( i ) = o_stoch
                    ! call o_stoch%kill
                end do
            end do
        else
            do i=1,self%n
                if( within_lims)then
                    call self%rnd_ori( i, eullims=eullims )
                else
                    call self%rnd_ori( i )
                endif
                call self%e3set(i,0.)
            end do
        endif
    end subroutine rnd_proj_space

    !>  \brief  randomizes eulers in oris
    subroutine rnd_ori( self, i, trs, eullims )
        use simple_rnd, only: ran3
        class(oris), intent(inout)         :: self
        integer,     intent(in)            :: i
        real,        intent(in), optional  :: trs
        real,        intent(inout), optional  :: eullims(3,2)
        call self%o(i)%rnd_ori( trs, eullims )
    end subroutine rnd_ori
    
    !>  \brief  for generating random clustering
    subroutine rnd_cls( self, ncls, srch_inpl )
        class(oris), intent(inout)    :: self
        integer, intent(in)           :: ncls
        logical, intent(in), optional :: srch_inpl
        logical :: ssrch_inpl
        integer :: i
        call self%rnd_classes(ncls)
        ssrch_inpl = .true.
        if( present(srch_inpl) ) ssrch_inpl = srch_inpl
        do i=1,self%n
            if( ssrch_inpl )then
                call self%o(i)%rnd_inpl
            endif
        end do
    end subroutine rnd_cls
    
    !>  \brief  randomizes eulers in oris
    subroutine rnd_oris_discrete( self, discrete, nsym, eullims )
        use simple_rnd, only: irnd_uni
        class(oris), intent(inout) :: self
        integer, intent(in)        :: discrete, nsym
        real, intent(in)           :: eullims(3,2)
        type(oris)                 :: tmp
        real                       :: euls(3)
        integer                    :: i, irnd
        euls(3) = 0.
        tmp = oris(discrete)
        call tmp%spiral(nsym, eullims)
        do i=1,self%n
            irnd    = irnd_uni(discrete)
            euls    = tmp%get_euler( irnd )
            euls(3) = self%o( i )%e3get()
            call self%o(i)%set_euler( euls )
        end do
        call tmp%kill
    end subroutine rnd_oris_discrete
    
    !>  \brief  randomizes the in-plane degrees of freedom
    subroutine rnd_inpls( self, trs )
        use simple_rnd, only: ran3
        class(oris), intent(inout) :: self
        real, intent(in), optional :: trs
        integer :: i
        real :: x, y
        do i=1,self%n
            if( present(trs) )then
                if( trs == 0. )then
                    x = 0.
                    y = 0. 
                else
                    x = ran3()*2.0*trs-trs
                    y = ran3()*2.0*trs-trs
                endif
            else
                x = 0.
                y = 0.
            endif
            call self%o(i)%e3set(ran3()*359.99)
            call self%o(i)%set('x', x)
            call self%o(i)%set('y', y)
        end do
    end subroutine rnd_inpls

    !>  \brief  randomizes the origin shifts
    subroutine rnd_trs( self, trs )
        use simple_rnd, only: ran3
        class(oris), intent(inout) :: self
        real,        intent(in)    :: trs
        integer :: i
        real :: x, y
        do i=1,self%n
            x = ran3()*2.0*trs-trs
            y = ran3()*2.0*trs-trs
            call self%o(i)%set('x', x)
            call self%o(i)%set('y', y)
        end do
    end subroutine rnd_trs
    
    !>  \brief  randomizes the CTF parameters
    subroutine rnd_ctf( self, kv, cs, fraca, defocus, deferr, astigerr )
        use simple_rnd, only: ran3
        class(oris),    intent(inout) :: self
        real,           intent(in)    :: kv, cs, fraca, defocus, deferr
        real, optional, intent(in)    :: astigerr
        integer :: i
        real    :: dfx, dfy, angast, err
        do i=1,self%n
            call self%o(i)%set('kv',    kv   )
            call self%o(i)%set('cs',    cs   )
            call self%o(i)%set('fraca', fraca)
            do
                err = ran3()*deferr
                if( ran3() < 0.5 )then
                    dfx = defocus-err
                else
                    dfx = defocus+err
                endif
                if( dfx > 0. ) exit
            end do
            call self%o(i)%set('dfx', dfx)
            if( present(astigerr) )then
                do
                    err = ran3()*astigerr
                    if( ran3() < 0.5 )then
                        dfy = dfx-err
                    else
                        dfy = dfx+err
                    endif
                    if( dfy > 0. ) exit
                end do
                angast = ran3()*359.99
                call self%o(i)%set('dfy', dfy)
                call self%o(i)%set('angast', angast)
            endif
        end do
    end subroutine rnd_ctf
    
    !>  \brief  balanced randomisation of states in oris 
    subroutine rnd_states( self, nstates )
        use simple_ran_tabu, only: ran_tabu
        class(oris), intent(inout) :: self
        integer, intent(in)        :: nstates
        integer, allocatable       :: states(:)
        type(ran_tabu)             :: rt
        integer :: alloc_stat, i, state
        if( nstates > 1 )then
            allocate( states(self%n), stat=alloc_stat )
            call alloc_err("simple_oris::rnd_states", alloc_stat)
            rt = ran_tabu(self%n)
            call rt%balanced(nstates, states) 
            do i=1,self%n
                state = nint(self%o(i)%get('state'))
                if( state /= 0 )then
                    call self%o(i)%set('state', real(states(i)))
                endif
            end do
            call rt%kill
            deallocate(states)
        else if( nstates<=0)then
            stop 'Invalid value for nstates; simple_oris :: rnd_states'
        else
            ! nstates = 1; zero-preserving
            do i=1,self%n
                state = nint(self%o(i)%get('state'))
                if(  state /= 0 )call self%o(i)%set('state', 1.)
            end do
        endif            
    end subroutine rnd_states

    !>  \brief  randomizes classes in oris
    subroutine rnd_classes( self, ncls )
        use simple_ran_tabu, only: ran_tabu
        class(oris), intent(inout) :: self
        integer, intent(in)        :: ncls
        integer, allocatable       :: classes(:)
        type(ran_tabu)             :: rt
        integer :: alloc_stat, i
        if( ncls > 1 )then
            allocate(classes(self%n), stat=alloc_stat)
            call alloc_err("simple_oris::rnd_classes", alloc_stat)
            rt = ran_tabu(self%n)
            call rt%balanced(ncls, classes)
            do i=1,self%n
                call self%o(i)%set('class', real(classes(i)))
            end do
            call rt%kill
            deallocate(classes)
        endif
    end subroutine rnd_classes
    
    !>  \brief  randomizes low-pass limits in oris
    subroutine rnd_lps( self )
        use simple_rnd, only: ran3
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('lp', ran3()*100.)
        end do
    end subroutine rnd_lps
    
    !>  \brief  randomizes correlations in oris 
    subroutine rnd_corrs( self )
        use simple_rnd, only: ran3
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('corr', ran3())
        end do
    end subroutine rnd_corrs
    
    !>  \brief  reverses the sign of the shifts in oris 
    subroutine revshsgn( self ) 
        class(oris), intent(inout) :: self
        integer :: i
        real :: x, y
        do i=1,self%n
            x = self%o(i)%get('x')
            y = self%o(i)%get('y')
            call self%o(i)%set('x', -x)
            call self%o(i)%set('y', -y)
        end do
    end subroutine revshsgn
    
    !>  \brief  reverses the sign of the shifts and rotations in oris 
    subroutine revorisgn( self ) 
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', -self%o(i)%get('x'))
            call self%o(i)%set('y', -self%o(i)%get('y'))
            call self%o(i)%set_euler(-self%o(i)%get_euler())
        end do
    end subroutine revorisgn
    
    !>  \brief  is for merging class class into class_merged
    subroutine merge_classes( self, class_merged, class )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: class_merged, class
        integer                    :: i, clsnr
        do i=1,self%n
            clsnr = nint(self%get(i, 'class'))
            if(clsnr == class) call self%set(i, 'class', real(class_merged))
        end do
    end subroutine merge_classes
    
    !>  \brief  for extending algndoc according to nr of symmetry ops
    subroutine symmetrize( self, nsym )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: nsym
        type(oris)                 :: tmp
        integer                    :: cnt, i, j
        tmp = oris(self%get_noris()*nsym)
        cnt = 0
        do i=1,self%get_noris()
            do j=1,nsym
                cnt = cnt+1
                call tmp%set_ori(cnt, self%get_ori(i))
            end do
        end do
        self = tmp
        call tmp%kill
    end subroutine symmetrize
    
    !>  \brief  for merging two oris objects into one
    subroutine merge( self1, self2add )
        class(oris), intent(inout) :: self1, self2add
        type(oris)                 :: self
        integer                    :: ntot, cnt, i
        ntot = self1%n+self2add%n
        self = oris(ntot)
        cnt  = 0
        if( self1%n > 0 )then
            do i=1,self1%n
                cnt = cnt+1
                self%o(cnt) = self1%o(i)
            end do
        endif
        if( self2add%n > 0 )then
            do i=1,self2add%n
                cnt = cnt+1
                self%o(cnt) = self2add%o(i)
            end do
        endif
        self1 = self 
        call self%kill
        call self2add%kill    
    end subroutine merge
    
    !>  \brief  for merging two oris files into one object
    subroutine merge_files( self, fname1, fname2 )
        use simple_jiffys, only: nlines, file_exists
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: fname1, fname2
        integer    :: nl1, nl2
        type(oris) :: o2
        logical    :: here1, here2
        call self%kill
        here1 = file_exists(fname1)
        here2 = file_exists(fname2)
        nl1   = 0
        nl2   = 0
        if(here1) nl1 = nlines(fname1)
        if(here2) nl2 = nlines(fname2)        
        if( here1 )then
            call self%new(nl1)
            call self%read(fname1)
            if( .not. here2 ) return
        else
            call self%new(nl2)
            call self%read(fname2)
            return
        endif
        call o2%new(nl2)
        call o2%read(fname2)
        call self%merge(o2)
    end subroutine merge_files
    
    ! I/O
    
    !>  \brief  reads orientation info from file
    subroutine read( self, orifile, nst )
        use simple_jiffys, only: get_fileunit
        class(oris), intent(inout)     :: self
        character(len=*), intent(in)   :: orifile
        integer, optional, intent(out) :: nst
        integer                        :: file_stat, i, fnr, state
        fnr = get_fileunit( )
        open(unit=fnr, FILE=orifile, STATUS='OLD', action='READ', iostat=file_stat)
        call fopen_err( 'In: read_oris, module: simple_oris.f90', file_stat )
        if( present(nst) ) nst = 0
        do i=1,self%n
            call self%o(i)%read(fnr)
            if( present(nst) )then
                state = int(self%o(i)%get('state'))
                nst = max(1,max(state,nst))
            endif
        end do
        close(fnr)
    end subroutine read
    
    !>  \brief  writes orientation info to file
    subroutine write_1( self, orifile )
        use simple_jiffys, only: get_fileunit
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: orifile
        integer                      :: file_stat, fnr, i
        fnr = get_fileunit( )
        open(unit=fnr, FILE=orifile, STATUS='REPLACE', action='WRITE', iostat=file_stat)
        call fopen_err( 'In: write_1, module: simple_oris.f90', file_stat )
        do i=1,self%n
            call self%o(i)%write(fnr)
        end do
        close(fnr)
    end subroutine write_1
    
    !>  \brief  writes orientation info to file
    subroutine write_2( self, i, orifile  )
        use simple_jiffys, only: get_fileunit
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: orifile
        integer, intent(in)          :: i
        integer                      :: fnr, file_stat
        fnr = get_fileunit( )
        open(unit=fnr, FILE=orifile, STATUS='UNKNOWN', action='WRITE', position='APPEND', iostat=file_stat)
        call fopen_err( 'In: write_2, module: simple_oris.f90', file_stat )
        call self%o(i)%write(fnr)
        close(fnr)
    end subroutine write_2
    
    ! CALCULATORS
    
    !>  \brief  for calculating homogeneity of labels within a state group or class
    subroutine homogeneity_1( self, icls, which, minp, thres, cnt, homo_cnt, homo_avg )
        class(oris), intent(inout)   :: self
        integer, intent(in)          :: icls
        character(len=*), intent(in) :: which
        integer, intent(in)          :: minp
        real, intent(in)             :: thres
        real, intent(inout)          :: cnt, homo_cnt, homo_avg
        integer                      :: pop, i, nlabels, medlab(1)
        real                         :: homo_cls
        real, allocatable            :: vals(:)
        integer, allocatable         :: labelpops(:)
        logical, parameter           :: debug=.false.
        nlabels = self%get_nlabels()
        if( debug ) print *, 'number of labels: ', nlabels
        allocate(labelpops(nlabels))
        if( debug ) print *, 'allocated labelpops'
        if( debug ) print *, '>>> PROCESSING CLASS: ', icls
        ! get pop 
        if( which .eq. 'class' )then
            pop = self%get_clspop(icls)
        else
            pop = self%get_statepop(icls)
        endif
        ! apply population treshold
        if( pop < minp ) return
        ! increment counter
        cnt = cnt+1.
        ! get label values
        if( which .eq. 'class' )then
            vals = self%get_arr('label', class=icls)
        else
            vals = self%get_arr('label', state=icls)
        endif
        if( debug ) print *, nint(vals)
        ! count label occurences
        labelpops = 0
        do i=1,pop
            labelpops(nint(vals(i))) = labelpops(nint(vals(i)))+1
        end do
        ! find median label
        medlab   = maxloc(labelpops)
        homo_cls = real(count(medlab(1) == nint(vals)))/real(pop)
        homo_avg = homo_avg+homo_cls
        if( debug ) print *, '>>> CLASS HOMOGENEITY: ', homo_cls
        ! apply homogeneity threshold
        if( homo_cls >= thres ) homo_cnt = homo_cnt+1.
        deallocate(labelpops)
    end subroutine homogeneity_1
    
    !>  \brief  for calculating homogeneity of labels within a state group or class
    subroutine homogeneity_2( self, which, minp, thres, homo_cnt, homo_avg )
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        integer, intent(in)          :: minp
        real, intent(in)             :: thres
        real, intent(out)            :: homo_cnt, homo_avg
        integer                      :: k, ncls
        real                         :: cnt
        logical, parameter           :: debug=.false.
        if( which .eq. 'class' )then
            ncls = self%get_ncls()
        else if( which .eq. 'state' )then
            ncls = self%get_nstates()
        else
            stop 'Unsupported which flag; simple_oris::homogeneity'
        endif
        if( debug ) print *, 'number of clusters to analyse: ', ncls
        homo_cnt = 0.
        homo_avg = 0.
        cnt      = 0.
        do k=1,ncls
            call self%homogeneity_1(k, which, minp, thres, cnt, homo_cnt, homo_avg)
        end do
        write(*,'(a,2x,i5)') 'This nr of clusters has pop >= minpop:      ', nint(cnt)
        write(*,'(a,2x,i5)') 'This is the number of homogeneous clusters: ', nint(homo_cnt)
        homo_cnt = homo_cnt/cnt
        homo_avg = homo_avg/cnt
    end subroutine homogeneity_2
    
    !>  \brief  for calculating clustering cohesion (normalised)
    function cohesion_norm( self, which, ncls ) result( coh_norm )
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        integer, intent(in)          :: ncls
        real :: coh_norm
        coh_norm = real(self%cohesion(which))/real(self%cohesion_ideal(ncls))
    end function cohesion_norm
    
    !>  \brief  for calculating clustering cohesion 
    function cohesion( self, which ) result( coh )
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        real, allocatable            :: vals(:)
        integer :: coh, ncls, pop, ival, jval, i, j, k
        if( which .eq. 'class' )then
            ncls = self%get_ncls()
        else if( which .eq. 'state' )then
            ncls = self%get_nstates()
        else
            stop 'Unsupported which flag; simple_oris::cohesion'
        endif
        coh = 0
        do k=1,ncls
            if( which .eq. 'class' )then
                pop = self%get_clspop(k)
            else
                pop = self%get_statepop(k)
            endif
            if( pop == 0 ) cycle
            if( pop == 1 )then
                coh = coh+1
                cycle
            endif
            if( which .eq. 'class' )then
                vals = self%get_arr('label', class=k)
            else
                vals = self%get_arr('label', state=k)
            endif
            do i=1,pop-1
                do j=i+1,pop
                    ival = nint(vals(i))
                    jval = nint(vals(j))
                    if( jval == ival ) coh = coh+1
                end do
            end do
            deallocate(vals)
        end do
    end function cohesion
    
    !>  \brief  for calculating the ideal clustering cohesion
    function cohesion_ideal( self, ncls ) result( coh )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: ncls
        integer :: coh, pop, i, j, k
        pop = nint(real(self%n)/real(ncls))
        coh = 0
        do k=1,ncls
            do i=1,pop-1
                do j=i+1,pop
                    coh = coh+1
                end do
            end do
        end do
    end function cohesion_ideal
    
    !>  \brief  for calculating the degree of separation between all classes based on label (normalised)
    function separation_norm( self, which, ncls ) result( sep_norm )
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        integer, intent(in)          :: ncls
        real :: sep_norm
        sep_norm = real(self%separation(which))/real(self%separation_ideal(ncls))
    end function separation_norm
    
    !>  \brief  for calculating the degree of separation between all classes based on label
    function separation( self, which ) result( sep )
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        integer, allocatable         :: iclsarr(:), jclsarr(:)
        integer :: iclass, jclass, ncls, sep
        integer :: i, j, ilab, jlab, ipop, jpop
        if( which .eq. 'class' )then
            ncls = self%get_ncls()
        else if( which .eq. 'state' )then
            ncls = self%get_nstates()
        else
            stop 'Unsupported which flag; simple_oris::cohesion'
        endif
        if( ncls == 1 )then
            sep = 1
            return
        endif
        sep = 0
        do iclass=1,ncls-1
            if( which .eq. 'class' )then
                ipop = self%get_clspop(iclass)
                if( ipop > 0 ) iclsarr = self%get_cls_pinds(iclass)
            else
                ipop = self%get_statepop(iclass)
                if( ipop > 0 ) iclsarr = self%get_state(iclass)
            endif
            if( ipop > 0 )then
                do jclass=iclass+1,ncls
                    if( which .eq. 'class' )then
                        jpop = self%get_clspop(jclass)
                        if( jpop > 0 ) jclsarr = self%get_cls_pinds(jclass)
                    else
                        jpop = self%get_statepop(jclass)
                        if( jpop > 0 ) jclsarr = self%get_state(jclass)
                    endif
                    if( jpop > 0 )then
                        do i=1,size(iclsarr)
                            do j=1,size(jclsarr)
                                ilab = nint(self%get(iclsarr(i), 'label'))
                                jlab = nint(self%get(jclsarr(j), 'label'))
                                if(jlab /= ilab ) sep = sep+1
                            end do
                        end do
                    endif
                    if( allocated(iclsarr) ) deallocate(iclsarr)
                    if( allocated(jclsarr) ) deallocate(jclsarr)
                end do
            endif
        end do
    end function separation
    
    !>  \brief  for calculating the degree of separation between all classes based on label
    function separation_ideal( self, ncls ) result( sep )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: ncls
        integer :: iclass, jclass, sep
        integer :: i, j, pop
        pop = nint(real(self%n)/real(ncls))
        if( ncls == 1 )then
            sep = 1
            return
        endif
        sep = 0
        do iclass=1,ncls-1
            do jclass=iclass+1,ncls
                do i=1,pop
                    do j=1,pop
                        sep = sep+1
                    end do
                end do
            end do
        end do
    end function separation_ideal

    !>  \brief  plots a histogram using gnuplot
    subroutine histogram( self, which, class, state )
        use simple_stat, only: plot_hist
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which 
        integer, intent(in), optional :: class
        integer, intent(in), optional :: state
        real, allocatable             :: vals(:)
        vals = self%get_arr( which, class, state )
        call plot_hist(vals, 100, 10)
        deallocate(vals)
    end subroutine histogram
    
    !>  \brief  for rounding the origin shifts
    subroutine round_shifts( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%round_shifts
        end do
    end subroutine round_shifts
    
    !>  \brief  for introducing alignment errors
    subroutine introd_alig_err( self, angerr, sherr )
        use simple_rnd, only: ran3
        class(oris), intent(inout) :: self
        real, intent(in) :: angerr, sherr
        real :: x, y, e1, e2, e3
        integer :: i
        do i=1,self%n
            e1 = self%e1get(i)+ran3()*angerr-angerr/2.
            if( e1 > 360. ) e1 = e1-360.
            if( e1 < 0. ) e1 = e1+360.
            call self%o(i)%e1set(e1)
            e2 = self%e2get(i)+ran3()*angerr-angerr/2.
            if( e2 > 180. ) e2 = e2-180.
            if( e2 < 0. ) e2 = e2+180.
            call self%o(i)%e2set(e2)
            e3 = self%e3get(i)+ran3()*angerr-angerr/2.
            if( e3 > 360. ) e3 = e3-360.
            if( e3 < 0. ) e3 = e3+360.
            call self%o(i)%e3set(e3)
            x = self%o(i)%get('x')
            y = self%o(i)%get('y')
            x = x+ran3()*sherr-sherr/2.
            y = y+ran3()*sherr-sherr/2.
            call self%o(i)%set('x', x)
            call self%o(i)%set('y', y)
        end do
    end subroutine introd_alig_err
    
    !>  \brief  for introducing CTF errors
    subroutine introd_ctf_err( self, dferr )
        use simple_rnd, only: ran3
        class(oris), intent(inout) :: self
        real, intent(in) :: dferr
        real    :: dfx, dfy
        integer :: i
        do i=1,self%n
            if( self%o(i)%isthere('dfx') )then
                do
                    dfx = self%o(i)%get('dfx')+ran3()*dferr-dferr/2.
                    if( dfx > 0. ) exit
                end do
                call self%o(i)%set('dfx', dfx)
            endif
            if( self%o(i)%isthere('dfy') )then
                do
                    dfy = self%o(i)%get('dfy')+ran3()*dferr-dferr/2.
                    if( dfy > 0. ) exit
                end do
                call self%o(i)%set('dfy', dfy)
            endif
        end do
    end subroutine introd_ctf_err
    
    !>  \brief  is an Euler angle composer
    subroutine rot_1( self, e )
        class(oris), intent(inout) :: self
        class(ori), intent(in)     :: e
        type(ori)                  :: o_tmp
        integer                    :: i
        do i=1,self%n
             o_tmp = e.compose.self%o(i)
             call self%o(i)%set_euler(o_tmp%get_euler())
        end do
    end subroutine rot_1
    
    !>  \brief  is an Euler angle composer
    subroutine rot_2( self, i, e )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: i
        class(ori), intent(in)     :: e
        type(ori)                  :: o_tmp
        o_tmp = e.compose.self%o(i)
        call self%o(i)%set_euler(o_tmp%get_euler())
    end subroutine rot_2
    
    !>  \brief  for identifying the median value of parameter which within the cluster class
    function median_1( self, which, class ) result( med )
        use simple_math, only: median_nocopy
        class(oris), intent(inout)    :: self
        character(len=*), intent(in)  :: which 
        integer, intent(in), optional :: class
        real, allocatable             :: vals(:)
        integer :: pop, alloc_stat, i
        real :: med
        if( present(class) )then        
            med = 0.
            pop = self%get_clspop(class)
            if( pop == 0 ) return
            vals = self%get_arr(which, class)
            if( pop == 1 )then
                med = vals(1)
                return
            endif
            if( pop == 2 )then
                med = (vals(1)+vals(2))/2.
                return
            endif
        else
            allocate( vals(self%n), stat=alloc_stat )
            call alloc_err('In: median_1, module: simple_oris', alloc_stat)
            vals = 0.
            do i=1,self%n
                vals(i) = self%o(i)%get(which)
            end do
        endif
        ! calculate median
        med = median_nocopy(vals)
        deallocate(vals)
    end function median_1
    
    !>  \brief  is for calculating variable statistics
    subroutine stats( self, which, ave, sdev, var, err )
        use simple_stat, only: moment
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        real, intent(out)            :: ave, sdev, var
        logical, intent(out)         :: err
        real, allocatable            :: vals(:)
        integer                      :: alloc_stat, i, cnt, mystate
        allocate( vals(self%n), stat=alloc_stat )
        call alloc_err('In: stat_1, module: simple_oris', alloc_stat)
        vals = 0.
        ! fish values
        cnt = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate /= 0 )then
                cnt = cnt+1
                vals(cnt) = self%o(i)%get(which)
            endif
        end do
        ! calculate statistics
        call moment(vals(:cnt), ave, sdev, var, err)
        deallocate(vals)
    end subroutine stats
    
    !>  \brief  is for calculating the minimum/maximum values of a variable
    subroutine minmax( self, which, minv, maxv )
        use simple_stat, only: moment
        class(oris), intent(inout)   :: self
        character(len=*), intent(in) :: which
        real, intent(out)            :: minv, maxv
        real, allocatable            :: vals(:)
        integer                      :: alloc_stat, i, cnt, mystate
        allocate( vals(self%n), stat=alloc_stat )
        call alloc_err('In: minmax, module: simple_oris', alloc_stat)
        vals = 0.
        ! fish values
        cnt = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate /= 0 )then
                cnt = cnt+1
                vals(cnt) = self%o(i)%get(which)
            endif
        end do
        minv = minval(vals(:cnt))
        maxv = maxval(vals(:cnt))
        deallocate(vals)
    end subroutine minmax
    
    !>  \brief  is for generating evenly distributed projection directions
    subroutine spiral_1( self )
        use simple_math, only: rad2deg
        class(oris), intent(inout) :: self
        real                       :: h, theta, psi
        integer                    :: k
        if( self%n == 1 )then
            call self%o(1)%set_euler([0.,0.,0.])
        else if( self%n > 1 )then
            do k=1,self%n
                h = -1.+((2.*(real(k)-1.))/(real(self%n)-1.))
                theta = acos(h)
                if( k == 1 .or. k == self%n )then
                    psi = 0.
                else
                    psi = psi+3.6/(sqrt(real(self%n))*sqrt(1.-real(h)**2.))
                endif
                do while( psi > 2.*pi )
                    psi = psi-2.*pi
                end do
                call self%o(k)%set_euler([rad2deg(psi),rad2deg(theta),0.])
            end do
        else
            stop 'object nonexistent; spiral_1; simple_oris' 
        endif 
    end subroutine spiral_1
    
    !>  \brief  is for generating evenly distributed projection directions 
    !!          within the asymetric unit
    subroutine spiral_2( self, nsym, eullims )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: nsym
        real, intent(in)           :: eullims(3,2)
        type(oris) :: tmp
        integer    :: cnt, i
        real       :: e1lim, e2lim
        if( nsym == 1 )then
            call self%spiral_1
            return
        endif
        e1lim = eullims(1,2)
        e2lim = eullims(2,2)
        tmp = oris(self%n*nsym)
        call tmp%spiral_1
        cnt = 0
        do i=1,self%n*nsym
            if( tmp%o(i)%e1get() <= e1lim .and. tmp%o(i)%e2get() <= e2lim )then
                cnt = cnt+1
                self%o(cnt) = tmp%o(i)
                if( cnt == self%n ) exit
            endif
         end do
    end subroutine spiral_2

    !>  \brief  is for generating a quasi-spiral within the asymetric unit
    subroutine qspiral( self, thresh, nsym, eullims )
        use simple_math, only: rad2deg
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsym
        real,        intent(inout) :: eullims(3,2)
        real,        intent(in)    :: thresh
        type(ori) :: o, o_prev
        integer   :: i
        logical   :: found
        call self%spiral_2( nsym, eullims )
        do i=1,self%n
            o_prev = self%o(i)
            o      = o_prev
            found  = .false.
            do while( .not.found )
                call o%rnd_euler( eullims )
                if( rad2deg( o.euldist.o_prev ) > 2.*thresh )cycle
                found = .true.
                call o%e3set( 0.)
                self%o( i ) = o
            end do
        enddo
    end subroutine qspiral

    !>  \brief  orders oris according to quality
    function order( self ) result( arr )
        class(oris), intent(in), target :: self
        integer, allocatable :: arr(:)
        integer :: i, alloc_stat
        op => self%o
        allocate( arr(self%n), stat=alloc_stat )
        call alloc_err('order; simple_oris', alloc_stat)
        arr = (/(i,i=1,self%n)/)
        call hpsort( self%n, arr, o1_gt_o2 )
    end function order
    
    !>  \brief  orders clusters according to population
    function order_cls( self ) result( arr )
        class(oris), intent(inout), target :: self
        integer, allocatable :: arr(:)
        integer :: i, alloc_stat, ncls
        ops => self
        ncls = self%get_ncls()
        allocate(arr(ncls), classpops(ncls), stat=alloc_stat)
        call alloc_err('order_cls; simple_oris', alloc_stat)
        ! calculate class populations
        do i=1,ncls
            classpops(i)  = self%get_clspop(i)
        end do
        arr = (/(i,i=1,ncls)/)
        call hpsort( ncls, arr, class1_gt_class2 )
        deallocate(classpops)
    end function order_cls
    
    !>  \brief  calculates hard weights based on ptcl ranking      
    subroutine calc_hard_ptcl_weights( self, frac, bystate )
        class(oris),       intent(inout) :: self
        real,              intent(in)    :: frac
        logical, optional, intent(in)    :: bystate
        type(oris)            :: os
        integer, allocatable  :: order(:), inds(:)
        integer               :: i, lim, n, nstates, s, pop
        logical               :: ibystate
        ibystate = .false.
        if( present(bystate) )ibystate = bystate
        if( .not.ibystate )then
            ! treated as single state
            order = self%order()
            n = 0
            do i=1,self%n
                if( self%o(i)%get("state") > 0 ) n = n+1
            end do        
            lim = nint(frac*real(n))
            do i=1,self%n
                if( i <= lim )then
                    call self%o(order(i))%set('w', 1.)
                else
                    call self%o(order(i))%set('w', 0.)
                endif           
            end do
            deallocate(order)
        else
            ! per state frac
            nstates = self%get_nstates()
            if( nstates==1 )then
                call self%calc_hard_ptcl_weights( frac, bystate=.false. )
            else
                do s=1,nstates
                    pop = self%get_statepop( s )
                    if( pop==0 )cycle
                    inds = self%get_state( s )
                    os = oris( pop )
                    do i=1,pop
                        call os%set_ori( i, self%get_ori( inds(i) ))
                    enddo
                    call os%calc_hard_ptcl_weights( frac, bystate=.false. )
                    do i=1,pop
                        call self%set_ori( inds(i), os%get_ori( i ))
                    enddo
                    deallocate(inds)
                enddo
            endif
        endif
    end subroutine calc_hard_ptcl_weights

    !>  \brief  to find the closest matching projection direction
    function find_closest_proj( self, o_in, state ) result( closest )
        class(oris),       intent(inout) :: self
        class(ori),        intent(in) :: o_in
        integer, optional, intent(in) :: state
        real     :: dists(self%n), large
        integer  :: loc(1), closest, i
        if( present(state) )then
            if( state<0 )stop 'Invalid state in simple_oris%find_closest_proj'
            dists(:) = huge(large)
            do i=1,self%n
                if( nint(self%o(i)%get('state'))==state )dists(i) = self%o(i).euldist.o_in
            end do
        else
            forall( i=1:self%n )dists(i)=self%o(i).euldist.o_in
        endif          
        loc     = minloc( dists )
        closest = loc(1)
    end function find_closest_proj

    !>  \brief  to find the closest matching projection directions
    subroutine find_closest_projs( self, o_in, pdirinds )
        class(oris), intent(in)  :: self
        class(ori),  intent(in)  :: o_in
        integer,     intent(out) :: pdirinds(:)
        real                     :: dists(self%n)
        integer                  :: inds(self%n), i
        do i=1,self%n
            dists(i) = self%o(i).euldist.o_in
            inds(i)  = i 
        end do
        call hpsort(self%n, dists, inds)
        do i=1,size(pdirinds)
            pdirinds(i) = inds(i)
        end do
    end subroutine find_closest_projs

    !>  \brief  to find the closest matching orientation
    function find_closest_ori( self, o_in, state ) result( closest )
        class(oris),       intent(inout) :: self
        class(ori),        intent(in) :: o_in
        integer, optional, intent(in) :: state
        real     :: dists(self%n), large
        integer  :: loc(1), closest, i
        if( present(state) )then
            if( state<0 )stop 'Invalid state in simple_oris%find_closest_proj'
            dists(:) = huge(large)
            do i=1,self%n
                if( nint(self%o(i)%get('state'))==state )dists(i) = self%o(i).geod.o_in
            end do
        else
            forall( i=1:self%n )dists(i)=self%o(i).geod.o_in
        endif          
        loc     = minloc( dists )
        closest = loc(1)
    end function find_closest_ori

    !>  \brief  to find the closest matching orientations
    subroutine find_closest_oris( self, o_in, oriinds )
        class(oris), intent(in)  :: self
        class(ori),  intent(in)  :: o_in
        integer,     intent(out) :: oriinds(:)
        real                     :: dists(self%n)
        integer                  :: inds(self%n), i
        do i=1,self%n
            dists(i) = self%o(i).geod.o_in
            inds(i)  = i 
        end do
        call hpsort(self%n, dists, inds)
        do i=1,size(oriinds)
            oriinds(i) = inds(i)
        end do
    end subroutine find_closest_oris
    
    !>  \brief  method for discretization of the projection directions
    subroutine discretize( self, n )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: n
        type(oris) :: d
        integer    :: closest, i
        if( n < self%n )then
            call d%new(n)
            call d%spiral
            do i=1,self%n
                closest = d%find_closest_proj(self%o(i))
                call self%o(i)%e1set(d%e1get(closest))
                call self%o(i)%e2set(d%e2get(closest))
                call self%o(i)%set('class', real(closest))
            end do
        else
            stop 'the number of discrete oris is too large; discretize; simple_oris'
        endif
    end subroutine discretize

    
    
    !>  \brief  to identify the indices of the k nearest neighbors (inclusive)
    function nearest_neighbors( self, k ) result( nnmat ) 
        class(oris), intent(inout) :: self
        integer, intent(in)        :: k
        integer, allocatable       :: nnmat(:,:)
        real                       :: dists(self%n)
        integer                    :: inds(self%n), i, j, alloc_stat
        type(ori)                  :: o
        if( k >= self%n ) stop 'need to identify fewer nearest_neighbors; simple_oris'
        if( allocated(nnmat) ) deallocate(nnmat)
        allocate( nnmat(self%n,k), stat=alloc_stat )
        call alloc_err("In: nearest_neighbors; simple_oris", alloc_stat)
        do i=1,self%n
            o = self%get_ori(i)
            do j=1,self%n
                inds(j) = j
                if( i == j )then
                    dists(j) = 0.
                else
                    dists(j) = self%o(j).geod.o
                endif
            end do
            call hpsort(self%n, dists, inds)
            do j=1,k
                nnmat(i,j) = inds(j)
            end do
        end do
    end function nearest_neighbors
    
    !>  \brief  to find angular resolution of an even orientation distribution (in degrees)
    function find_angres( self ) result( res )
        use simple_math, only: rad2deg
        class(oris), intent(in) :: self
        real                    :: dists(self%n), res, x
        integer                 :: i, j
        res = 0.
        do j=1,self%n
            do i=1,self%n
                if( i == j )then
                    dists(i) = huge(x)
                else
                    dists(i) = self%o(i).euldist.self%o(j)
                endif
            end do
            call hpsort(self%n, dists)
            res = res+sum(dists(:3))/3. ! average of three nearest neighbors
        end do
        res = rad2deg(res/real(self%n))  
    end function find_angres

    !>  \brief  to find angular resolution of an even orientation distribution (in degrees)
    function find_angres_geod( self ) result( res )
        use simple_math, only: rad2deg
        class(oris), intent(in) :: self
        real                    :: dists(self%n), res, x
        integer                 :: i, j
        res = 0.
        do j=1,self%n
            do i=1,self%n
                if( i == j )then
                    dists(i) = huge(x)
                else
                    dists(i) = self%o(i).geodsc.self%o(j)
                endif
            end do
            call hpsort(self%n, dists)
            res = res+sum(dists(:3))/3. ! average of three nearest neighbors
        end do
        res = rad2deg(res/real(self%n))  
    end function find_angres_geod
    
    !>  \brief  to find the neighborhood size for weighted orientation assignment
    function find_npeaks( self, res, moldiam ) result( npeaks )
        class(oris), intent(in) :: self
        real, intent(in)        :: res, moldiam
        real                    :: dists(self%n), tres, npeaksum
        integer                 :: i, j, npeaks
        tres = atan(res/(moldiam/2.))
        npeaksum = 0.
        !$omp parallel do schedule(auto) default(shared) private(j,i,dists) reduction(+:npeaksum)
        do j=1,self%n
            do i=1,self%n
                dists(i) = self%o(i).euldist.self%o(j)
            end do
            call hpsort(self%n, dists)
            do i=1,self%n
                if( dists(i) <= tres )then
                    npeaksum = npeaksum+1.
                else
                    exit
                endif 
            end do
        end do
        !$omp end parallel do
        npeaks = nint(npeaksum/real(self%n)) ! nr of peaks is average of peaksum
    end function find_npeaks
    
    !>  \brief  calculate fraction of ptcls in class within lims
    function class_calc_frac( self, class, lims ) result( frac )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: class, lims(2)
        integer, allocatable       :: clsarr(:)
        integer                    :: i, cnt, pop
        real                       :: frac
        clsarr = self%get_cls_pinds(class)
        cnt = 0
        do i=1,size(clsarr)
            if( clsarr(i) >= lims(1) .and. clsarr(i) <= lims(2) ) cnt = cnt+1
        end do
        pop = size(clsarr)
        frac = 0.
        if( pop > 0 ) frac = real(cnt)/real(pop)
        deallocate(clsarr)  
    end function class_calc_frac
    
    !>  \brief  is for calculating the average and variance of euler distances (in degrees) within a class
    subroutine class_dist_stat( self, class, distavg, distsdev )
    ! can envision many other class statistics to calculate that would be interesting
    ! this routine could also be used to reduce bias in refinement (enforce class variance reduction for new orientations)
    ! could also be used to measure convergence (while reducing class variance, continue, else stop)
        use simple_math, only: rad2deg
        use simple_stat, only: moment
        class(oris), intent(inout) :: self
        integer, intent(in)        :: class
        real, intent(out)          :: distavg, distsdev 
        integer, allocatable       :: clsarr(:)
        integer                    :: i, j, cnt, npairs, alloc_stat, pop
        real, allocatable          :: dists(:)
        logical                    :: err
        real                       :: distvar
        clsarr = self%get_cls_pinds(class)
        pop    = size(clsarr)
        cnt    = 0
        npairs = (pop*(pop-1))/2
        allocate( dists(npairs), stat=alloc_stat )
        call alloc_err( 'class_dist_stat; simple_oris', alloc_stat )
        cnt = 0
        distavg = 0.
        do i=1,pop-1
            do j=i+1,pop
                cnt = cnt+1
                dists(cnt) = self%o(clsarr(i)).euldist.self%o(clsarr(j))                
            end do
        end do
        call moment( dists, distavg, distsdev, distvar, err )
        distavg  = rad2deg(distavg)
        distsdev = rad2deg(distsdev)
        deallocate(clsarr)
    end subroutine class_dist_stat
    
    !>  \brief  for calculating the average correlation within a class
    function class_corravg( self, class ) result( avg )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: class
        integer                    :: i, cnt, clsnr
        real                       :: corr, avg
        avg = 0.
        cnt = 0
        do i=1,self%n
            clsnr = nint(self%get(i, 'class'))
            corr  = self%get(i, 'corr')
            if( clsnr == class )then
                avg = avg+corr
                cnt = cnt+1
            endif
        end do
        avg = avg/real(cnt)
    end function class_corravg
    
    !>  \brief  for finding the best (highest correlating) member in a class
    function best_in_class( self, class ) result( best )
        class(oris), intent(inout) :: self
        integer, intent(in)        :: class
        integer                    :: i, best, clsnr
        real                       :: corr, best_corr
        best_corr = -1.
        best = 0
        do i=1,self%n
            clsnr = nint(self%get(i, 'class'))
            corr  = self%get(i, 'corr')
            if( clsnr == class )then
                if( corr > best_corr )then
                    best_corr = corr
                    best = i
                endif
            endif
        end do
    end function best_in_class
    
    !>  \brief  for calculating the average class shifts
    subroutine calc_avg_class_shifts_1( self, avg_shifts )
        class(oris),       intent(inout) :: self
        real, allocatable, intent(out)   :: avg_shifts(:,:)
        real    :: classpop
        integer :: ncls, icls, iptcl, mycls
        ncls = self%get_ncls()
        if( allocated(avg_shifts) ) deallocate(avg_shifts)
        allocate(avg_shifts(ncls,2))
        avg_shifts = 0.
        do icls=1,ncls
            classpop = 0.0
            do iptcl=1,self%n
                mycls = nint(self%o(iptcl)%get('class'))
                if( mycls == icls )then
                    classpop           = classpop           + 1.0
                    avg_shifts(icls,1) = avg_shifts(icls,1) + self%o(iptcl)%get('x')
                    avg_shifts(icls,2) = avg_shifts(icls,2) + self%o(iptcl)%get('y')
                endif
            end do
            if( classpop > 0.5 ) avg_shifts(icls,:) = avg_shifts(icls,:)/classpop
        end do
    end subroutine calc_avg_class_shifts_1
    
    !>  \brief  for calculating the average class shifts
    subroutine calc_avg_class_shifts_2( self, o_avg_shifts )
        class(oris), intent(inout) :: self
        real, allocatable          :: avg_shifts(:,:)
        type(oris)                 :: o_avg_shifts
        integer                    :: ncls, icls
        call self%calc_avg_class_shifts_1(avg_shifts)
        ncls = size(avg_shifts,1)
        call o_avg_shifts%new(ncls)
        do icls=1,ncls
            call o_avg_shifts%set(icls, 'x', avg_shifts(icls,1))
            call o_avg_shifts%set(icls, 'y', avg_shifts(icls,2))
        end do
        deallocate(avg_shifts)
    end subroutine calc_avg_class_shifts_2
    
    !>  \brief  for mapping a 3D shift of volume to 2D shifts of the projections
    subroutine map3dshift22d_1( self, sh3d, state )
        class(oris), intent(inout)    :: self
        real, intent(in)              :: sh3d(3)
        integer, optional, intent(in) :: state
        integer :: i
        do i=1,self%n
            call self%map3dshift22d_2(i, sh3d, state)
        end do
    end subroutine map3dshift22d_1
    
    !>  \brief  for mapping a 3D shift of volume to 2D shifts of the projections
    subroutine map3dshift22d_2( self, i, sh3d, state )
        class(oris), intent(inout)    :: self
        integer, intent(in)           :: i
        real, intent(in)              :: sh3d(3)
        integer, optional, intent(in) :: state
        integer :: mystate
        if( present(state) )then
            mystate = nint(self%o(i)%get('state'))
            if( mystate == state ) call self%o(i)%map3dshift22d(sh3d)
        else
            call self%o(i)%map3dshift22d(sh3d)
        endif
    end subroutine map3dshift22d_2

    !>  \brief  modulates the shifts (additive) within a class
    subroutine add_shift2class( self, class, sh2d )
        use simple_math, only: rotmat2d
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: class
        real,        intent(in)    :: sh2d(2)
        integer :: i
        real    :: sh3d(3)
        sh3d(1:2) = sh2d
        sh3d(3)   = 0.
        do i=1,self%n
            if( nint(self%o(i)%get('class')) == class )then
                call self%o(i)%map3dshift22d(sh3d)
            endif
        end do
    end subroutine add_shift2class

    ! PRIVATE STUFF
    
    !>  \brief  orientation 1 greater than (better) than orientation 2 ?
    function o1_gt_o2( o1, o2 ) result( val )
        integer, intent(in) :: o1, o2
        logical             :: val
        real                :: corr1, corr2
        corr1 = op(o1)%get('corr')
        corr2 = op(o2)%get('corr')
        if( corr1 > corr2 )then
            val = .true.
        else
            val = .false.
        endif
    end function o1_gt_o2
    
    !>  \brief  orientation 1 less than (worse) than orientation 2 ?
    function o1_lt_o2( o1, o2 ) result( val )
        integer, intent(in) :: o1, o2
        logical             :: val
        val = .not. o1_gt_o2( o1, o2 )
    end function o1_lt_o2
    
    !>  \brief  class 1 greater than (better) than class 2 ?
    function class1_gt_class2( class1, class2 ) result( val )
        integer, intent(in) :: class1, class2
        logical             :: val
        integer             :: pop1, pop2
        pop1  = classpops(class1)
        pop2  = classpops(class2)
        if( pop1 > pop2 )then
            val = .true.
        else
            val = .false.
        endif
    end function class1_gt_class2
    
    !>  \brief  class 1 less than (worse) than class 2 ?
    function class1_lt_class2( class1, class2 ) result( val )
        integer, intent(in) :: class1, class2
        logical             :: val
        val = .not. class1_gt_class2(class1, class2)
    end function class1_lt_class2
    
    !>  \brief  for correlating oris objs, just for testing purposes
    function corr_oris( self1, self2 ) result( corr )
        use simple_stat, only: pearsn
        class(oris), intent(inout) :: self1, self2
        real :: arr1(5), arr2(5), corr
        integer :: i
        corr = 0.
        do i=1,self1%n
            arr1(1:3) = self1%get_euler(i) 
            arr1(4)   = self1%get(i,'x')
            arr1(5)   = self1%get(i,'y')
            arr2(1:3) = self2%get_euler(i) 
            arr2(4)   = self2%get(i,'x')
            arr2(5)   = self2%get(i,'y')
            corr = corr+pearsn(arr1,arr2)
        end do
        corr = corr/real(self1%n)
    end function corr_oris

    !>  \brief  uses a greedy approach to generate a maximally diverse set of 
    !!          orientations by maximizing the geodesic distance upon every 
    !!          addition to the growing set
    subroutine gen_diverse( self )
        class(oris), intent(inout) :: self
        logical, allocatable       :: o_is_set(:)
        integer                    :: i
        type(ori)                  :: o1, o2
        ! create the first diverse pair
        call o1%oripair_diverse(o2)
        allocate(o_is_set(self%n))
        o_is_set = .false.
        call self%set_ori(1,o1)
        call self%set_ori(2,o2)
        o_is_set(1) = .true.
        o_is_set(2) = .true.
        ! use a greedy approach to generate the maximally diverse set
        do i=3,self%n
            call progress(i,self%n)
            o1 = self%ori_generator('diverse', o_is_set)
            call self%set_ori(i,o1)
            o_is_set(i) = .true.
        end do
        deallocate(o_is_set)
    end subroutine gen_diverse

    !>  \brief  scores the orientation (o) according to diversity with 
    !!          respect to the orientations in self
    function gen_diversity_score( self, o ) result( diversity_score )
        class(oris), intent(in) :: self
        class(ori),  intent(in) :: o
        real    :: diversity_score
        real    :: dist, dist_min
        integer :: i, j
        dist_min = 2.0*sqrt(2.0)
        do i=1,self%n
            dist = self%o(i).geod.o
            if( dist < dist_min ) dist_min = dist
        end do
        diversity_score = dist_min
    end function gen_diversity_score

    !>  \brief  scores the orientations in self according to diversity with respect to ref_set
    function gen_diversity_scores( self, ref_set ) result( diversity_scores )
        class(oris), intent(in) :: self, ref_set
        real, allocatable :: diversity_scores(:)
        real    :: dist, dist_min
        integer :: i, j
        allocate(diversity_scores(self%n))
        do i=1,self%n
            dist_min = 2.0*sqrt(2.0)
            do j=1,ref_set%n
                dist = self%o(i).geod.ref_set%o(j)
                if( dist < dist_min ) dist_min = dist
            end do
            diversity_scores(i) = dist_min
        end do
    end function gen_diversity_scores

    !>  \brief  if mode='median' this routine creates a spatial median rotation matrix that 
    !!          is representative of the rotation matrices in the instance in a geodesic 
    !!          distance sense. If mode='diverse' it creates the rotation matrix that is 
    !!          maximally diverse to the rotation matrices in the instance 
    function ori_generator( self, mode, part_of_set, weights ) result( oout )
        use simple_simplex_opt, only: simplex_opt
        use simple_opt_spec,    only: opt_spec
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: mode
        logical, optional, intent(in) :: part_of_set(self%n)
        real,    optional, intent(in) :: weights(self%n)
        real, parameter               :: TOL=1e-6
        type(ori)                     :: oout, otst
        type(opt_spec)                :: ospec
        type(simplex_opt)             :: opt
        real                          :: dist, lims(3,2)
        logical                       :: ppart_of_set(self%n)
        real                          :: wweights(self%n)
        if( present(part_of_set) )then
            ppart_of_set = part_of_set
        else
            ppart_of_set = .true.
        endif
        if( present(weights) )then
            wweights = weights
        else
            wweights = 1.0
        endif
        call oout%new
        call otst%new
        call otst%rnd_euler
        lims(1,1) = 0.
        lims(1,2) = 359.99
        lims(2,1) = 0.
        lims(2,2) = 180.
        lims(3,1) = 0.
        lims(3,2) = 359.99
        call ospec%specify('simplex', 3, ftol=TOL, limits=lims)
        ospec%x(1) = otst%e1get()
        ospec%x(2) = otst%e2get()
        ospec%x(3) = otst%e3get()
        call opt%new(ospec)
        select case(mode)
            case('median')
                call ospec%set_costfun(costfun_median)
            case('diverse')
                call ospec%set_costfun(costfun_diverse)
            case DEFAULT
                stop 'unsupported mode inputted; simple_ori :: ori_generator'
        end select
        call opt%minimize(ospec, dist)
        call oout%set_euler(ospec%x)
        call ospec%kill
        call opt%kill

        contains

            function costfun_median( vec, D ) result( dist )
                integer, intent(in) :: D
                real,    intent(in) :: vec(D)
                real :: dist
                call otst%set_euler(vec)
                ! we are minimizing the distance 
                dist = self%geodesic_dist(otst,ppart_of_set,wweights)
            end function costfun_median

            function costfun_diverse( vec, D ) result( dist )
                integer, intent(in) :: D
                real,    intent(in) :: vec(D)
                real :: dist
                call otst%set_euler(vec)
                ! we are maximizing the distance 
                dist = -self%geodesic_dist(otst,ppart_of_set,wweights)
            end function costfun_diverse

    end function ori_generator

    !>  \brief  generates a similarity matrix for the oris in self
    function gen_smat( self ) result( smat )
        class(oris), intent(in) :: self
        real, allocatable :: smat(:,:)
        integer :: alloc_stat, i, j
        allocate( smat(self%n,self%n), stat=alloc_stat )
        call alloc_err("In: simple_oris :: gen_smat", alloc_stat)
        smat = 0.
        !$omp parallel do schedule(auto) default(shared) private(i,j)
        do i=1,self%n-1
            do j=i+1,self%n
                smat(i,j) = -(self%o(i).geod.self%o(j))
                smat(j,i) = smat(i,j)
            end do
        end do
        !$omp end parallel do
    end function gen_smat

    !>  \brief  calculates the average geodesic distance between the input orientation
    !!          and all other orientations in the instance, part_of_set acts as mask
    real function geodesic_dist( self, o, part_of_set, weights )
        class(oris),       intent(in) :: self
        class(ori),        intent(in) :: o
        logical, optional, intent(in) :: part_of_set(self%n)
        real,    optional, intent(in) :: weights(self%n)
        integer :: i
        real    :: dists(self%n)
        logical :: ppart_of_set(self%n)
        real    :: wweights(self%n)
        if( present(part_of_set) )then
            ppart_of_set = part_of_set
        else
            ppart_of_set = .true.
        endif
        if( present(weights) )then
            wweights = weights
        else
            wweights = 1.0
        endif
        dists = 0.
        !$omp parallel do schedule(auto) default(shared) private(i)
        do i=1,self%n
            dists(i) = self%o(i).geod.o
        end do
        !$omp end parallel do
        geodesic_dist = sum(dists*wweights,mask=ppart_of_set)/sum(wweights,mask=ppart_of_set)
    end function geodesic_dist

    !>  \brief  for generation a subset of orientations (spatial medians of all pairs)
    function gen_subset( self ) result( os )
        class(oris), intent(in) :: self
        type(ori)  :: omed
        type(oris) :: os
        integer    :: iref, jref, nsub, cnt
        nsub = (self%n*(self%n-1))/2
        call os%new(nsub)
        cnt = 0
        do iref=1,self%n-1
            do jref=iref+1,self%n
                cnt  = cnt + 1
                omed = self%o(iref)%ori_generator(self%o(jref), 'median')
                call os%set_ori(cnt, omed)
            end do
        end do
    end function gen_subset
    
    !>  \brief  for calculating statistics of distances within a single distribution
    subroutine diststat_1( self, sumd, avgd, sdevd, mind, maxd )
        use simple_stat, only: moment
        class(oris), intent(in) :: self
        real, intent(out)       :: mind, maxd, avgd, sdevd, sumd
        integer :: i, j
        real    :: dists((self%n*(self%n-1))/2), vard, x
        logical :: err
        mind = huge(x)
        maxd = -mind
        sumd = 0.
        do i=1,self%n-1
            do j=i+1,self%n
                 dists(i) = self%o(i).euldist.self%o(j)
                 if( dists(i) < mind ) mind = dists(i)
                 if( dists(i) > maxd ) maxd = dists(i)
                 sumd = sumd+dists(i)
            end do
        end do
        call moment(dists, avgd, sdevd, vard, err )   
    end subroutine diststat_1

    !>  \brief  for calculating statistics of distances between two equally sized distributions
    subroutine diststat_2( self1, self2, sumd, avgd, sdevd, mind, maxd )
        use simple_stat, only: moment
        class(oris), intent(in) :: self1, self2
        real, intent(out)       :: mind, maxd, avgd, sdevd, sumd
        integer :: i
        real    :: dists(self1%n), vard, x
        logical :: err
        if( self1%n /= self2%n )then
            stop 'cannot calculate distance between sets of different size; euldist_2; simple_oris'
        endif
        mind = huge(x)
        maxd = -mind
        sumd = 0.
        do i=1,self1%n
            dists(i) = (self1%o(i).euldist.self2%o(i))
            if( dists(i) < mind ) mind = dists(i)
            if( dists(i) > maxd ) maxd = dists(i)
            sumd = sumd+dists(i)
        end do
        call moment(dists, avgd, sdevd, vard, err )
    end subroutine diststat_2

    !>  \brief  for calculating statistics of geodesic distances within clusters of orientations
    subroutine cluster_diststat( self, avgd, sdevd, maxd, mind )
        use simple_stat, only: moment
        class(oris), intent(inout) :: self
        real,        intent(out)   :: avgd, sdevd, maxd, mind
        integer, allocatable       :: clsarr(:)
        real,    allocatable       :: dists(:)
        integer                    :: ncls, icls, iptcl, jptcl, sz, cnt, nobs, n
        real                       :: avgd_here, sdevd_here, vard_here, mind_here, maxd_here
        type(ori)                  :: iori, jori
        logical                    :: err
        ncls  = self%get_ncls()
        nobs  = 0
        avgd  = 0.
        sdevd = 0.
        maxd  = 0.
        mind  = 0.
        do icls=1,ncls
            clsarr = self%get_cls_pinds(icls)
            if( allocated(clsarr) )then
                sz = size(clsarr)
                n  = (sz*(sz-1))/2
                allocate( dists(n) )
                cnt       = 0
                mind_here = 2.*pi
                maxd_here = 0.
                do iptcl=1,sz-1
                    do jptcl=iptcl+1,sz
                        cnt        = cnt + 1
                        iori       = self%get_ori(clsarr(iptcl))
                        jori       = self%get_ori(clsarr(jptcl))
                        dists(cnt) = iori.geodsc.jori
                        if( dists(cnt) < mind_here ) mind_here = dists(cnt)
                        if( dists(cnt) > maxd_here ) maxd_here = dists(cnt)
                    end do
                end do
                call moment(dists, avgd_here, sdevd_here, vard_here, err)
                avgd  = avgd  + avgd_here
                sdevd = sdevd + sdevd_here
                maxd  = maxd  + maxd_here
                mind  = mind  + mind_here
                nobs  = nobs  + 1
                deallocate( dists, clsarr )
            endif
        end do
        avgd  = avgd/real(nobs)
        sdevd = sdevd/real(nobs)
        maxd  = maxd/real(nobs)
        mind  = mind/real(nobs)
    end subroutine cluster_diststat
    
    !>  \brief  generates the opposite hand of the set of Euler angles
    subroutine mirror3d( self )
        class(oris), intent(inout) :: self
        integer                    :: i
        do i=1,self%n
            call self%o(i)%mirror3d
        end do
    end subroutine mirror3d
    
    !>  \brief  generates the opposite hand of the set of Euler angles
    !!          according to the spider convention
    subroutine mirror2d( self )
        class(oris), intent(inout) :: self
        integer                    :: i
        do i=1,self%n
            call self%o(i)%mirror2d
        end do
    end subroutine mirror2d
    
    !>  \brief  returns cluster population overlap between two already clustered orientations sets
    function cls_overlap( self, os, key )result( overlap )
        class(oris), intent(inout)        :: self ! previous oris
        type(oris), intent(inout)         :: os   ! new oris (eg number of clusters clusters > number of self clusters)
        character(len=*), intent(in)      :: key
        integer, allocatable              :: cls_pop(:), pop_spawn(:)
        real                              :: overlap, ov
        integer                           :: i, iptcl, pop, pop_prev, ncls, ncls_prev
        integer                           :: alloc_stat, loc(1), membership, prev_membership
        overlap   = 0.
        ncls      = os%get_nstates()
        ncls_prev = self%get_nstates()
        if( ncls<ncls_prev)stop 'inconsistent arguments order in simple_oris; overlap'
        if( (key.ne.'class').and.(key.ne.'state') ) stop 'Unsupported field in simple_oris; overlap'
        allocate( cls_pop(ncls), pop_spawn(ncls_prev), stat=alloc_stat)
        call alloc_err('In: simple_oris; overlap', alloc_stat)
        do i=1,ncls_prev
            ! Identifies the newly created cluster (with max pop overlap to current prev cluster)
            cls_pop(:) = 0
            do iptcl=1,self%n
                prev_membership = nint( self%get( iptcl,key ) )
                if( prev_membership.eq.i)then
                    membership          = nint( os%get( iptcl,key ) )
                    cls_pop(membership) = cls_pop(membership)+1
                endif
            enddo
            loc          = maxloc(cls_pop) ! location of new cluster with max pop belonging to prev cluster
            pop_spawn(i) = cls_pop(loc(1))
        enddo
        ! Calculates overlap between prev and new clusters
        do i=1,ncls_prev        
            pop      = pop_spawn(i) ! pop of new cluster spawning from prev cluster
            if(key.eq.'state')then
                pop_prev = self%get_statepop( i )
            else
                pop_prev = self%get_clspop( i )
            endif                
            if( (pop>0).and.(pop_prev>0) )then
                ov       = real(pop)/real(pop_prev)
                overlap  = overlap+ov
            endif
        enddo
        overlap = overlap/real(ncls_prev)
        deallocate(cls_pop,pop_spawn)
    end function cls_overlap
    
    ! UNIT TESTS
    
    !>  \brief  oris class unit test
    subroutine test_oris( doprint )
        logical, intent(in)  :: doprint
        type(oris)           :: os, os2
        real                 :: euls(3), corr, x, x2, y, y2
        integer              :: i
        integer, allocatable :: order(:)
        logical              :: passed
        write(*,'(a)') '**info(simple_oris_unit_test, part1): testing getters/setters'
        os = oris(100)
        os2 = oris(100)
        passed = .false.
        if( os%get_noris() == 100 ) passed = .true.
        if( .not. passed ) stop 'get_noris failed!'
        passed = .false.
        call os%set_euler(1, [1.,2.,3.])
        euls = os%get_euler(1)
        if( abs(euls(1)-1.+euls(2)-2.+euls(3)-3.) < 0.0001 ) passed = .true.
        if( .not. passed ) stop 'get/set eulers failed!'
        passed = .false.
        call os%e1set(1,4.)
        call os%e2set(1,5.)
        call os%e3set(1,6.) 
        euls(1) = os%e1get(1)
        euls(2) = os%e2get(1)
        euls(3) = os%e3get(1)
        if( abs(euls(1)-1.+euls(2)-2.+euls(3)-3.) < 0.0001 ) passed = .true.
        if( doprint )then
            call os%rnd_oris(5.)
            write(*,*) '********'
            do i=1,100
                call os%print(i)
            end do
            call os2%rnd_oris(5.)
            write(*,*) '********'
            do i=1,100
                call os2%print(i)
            end do
            call os%merge(os2)
            write(*,*) '********'
            do i=1,200
                call os%print(i)
            end do
        endif
        write(*,'(a)') '**info(simple_oris_unit_test, part2): testing assignment'
        os = oris(2)
        call os%rnd_oris(5.)
        os2 = os
        do i=1,2
            x  = os%get(i,'x')
            x2 = os2%get(i,'x')
            passed = (abs(x-x2)<TINY)
            if( passed )then
                cycle
            else
                exit
            endif
            y  = os%get(i,'y')
            y2 = os2%get(i,'y')
            passed = (abs(y-y2)<TINY)
            if( passed )then
                cycle
            else
                exit
            endif
            passed = all((os%get_euler(i)-os2%get_euler(i))<TINY)
            if( passed )then
                cycle
            else
                exit
            endif
        end do
        if( .not. passed ) stop 'assignment test failed!'
        write(*,'(a)') '**info(simple_oris_unit_test, part2): testing i/o'
        passed = .false.
        os = oris(100) 
        os2 = oris(100) 
        call os%rnd_oris(5.)
        call os%write('test_oris_rndoris.txt')
        call os2%read('test_oris_rndoris.txt')
        call os2%write('test_oris_rndoris_copy.txt')
        corr = corr_oris(os,os2)
        if( corr > 0.99 ) passed = .true.
        if( .not. passed ) stop 'read/write failed'
        passed = .false.
        call os%rnd_states(5)
        call os%write('test_oris_rndoris_rndstates.txt')
        if( corr_oris(os,os2) > 0.99 ) passed = .true.
        if( .not. passed ) stop 'statedoc read/write failed!'
        write(*,'(a)') '**info(simple_oris_unit_test, part3): testing calculators'
        passed = .false.
        call os%rnd_lps()
        call os%write('test_oris_rndoris_rndstates_rndlps.txt')
        call os%spiral
        call os%write('test_oris_rndoris_rndstates_rndlps_spiral.txt')
        call os%rnd_corrs()
        order = os%order()
        if( doprint )then
            do i=1,100
                call os%print(order(i))
            end do
            write(*,*) 'median:', os%median('lp')
        endif
        write(*,'(a)') '**info(simple_oris_unit_test, part4): testing destructor'
        call os%kill
        call os2%kill
        ! test find_angres vs. find_angres_geod
        os = oris(1000)
        call os%spiral
        print *, 'angres:      ', os%find_angres()
        print *, 'angres_geod: ', os%find_angres_geod()
        write(*,'(a)') 'SIMPLE_ORIS_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_oris
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(oris), intent(inout) :: self
        integer :: i
        if( allocated(self%o) )then
            deallocate( self%o )
        endif
    end subroutine kill
    
end module simple_oris