! an agglomeration of orientations
module simple_oris
!$ use omp_lib
!$ use omp_lib_kinds
#include "simple_lib.f08"

use simple_ran_tabu,    only: ran_tabu
use simple_ori,         only: ori
implicit none

public :: oris, test_oris
private
#include "simple_local_flags.inc"

!>  \brief struct type aggregates ori objects
type :: oris
    private
    type(ori), allocatable :: o(:)
    integer :: n=0
  contains
    ! CONSTRUCTORS
    procedure          :: new
    procedure          :: new_clean
    ! GETTERS
    procedure          :: e1get
    procedure          :: e2get
    procedure          :: e3get
    procedure          :: get_euler
    procedure          :: get_noris
    procedure          :: get_ori
    procedure          :: get
    procedure, private :: getter_1
    procedure, private :: getter_2
    generic            :: getter => getter_1, getter_2
    procedure          :: get_all
    procedure, private :: getter_all_1
    procedure, private :: getter_all_2
    generic            :: getter_all => getter_all_1, getter_all_2
    procedure          :: get_mat
    procedure          :: get_normal
    procedure          :: get_2Dshift
    procedure, private :: isthere_1
    procedure, private :: isthere_2
    generic            :: isthere => isthere_1, isthere_2
    procedure          :: max_hash_size
    procedure          :: get_n
    procedure          :: get_pop
    procedure          :: get_pops
    procedure          :: get_pinds
    procedure          :: gen_mask
    procedure          :: states_exist
    procedure          :: get_arr
    procedure, private :: calc_sum
    procedure          :: get_sum
    procedure          :: get_avg
    procedure, private :: calc_nonzero_sum
    procedure          :: get_nonzero_sum
    procedure          :: get_nonzero_avg
    procedure          :: extract_table
    procedure          :: ang_sdev
    procedure          :: included
    procedure          :: get_nevenodd
    procedure, private :: get_neven
    procedure, private :: get_nodd
    procedure          :: print_
    procedure          :: print_matrices
    procedure          :: sample4update_and_incrcnt
    ! SETTERS
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: reject
    procedure          :: set_euler
    procedure          :: set_shift
    procedure          :: e1set
    procedure          :: e2set
    procedure          :: e3set
    procedure, private :: set_1
    procedure, private :: set_2
    generic            :: set => set_1, set_2
    procedure          :: set_ori
    procedure, private :: set_all_1
    procedure, private :: set_all_2
    generic            :: set_all => set_all_1, set_all_2
    procedure, private :: set_all2single_1
    procedure, private :: set_all2single_2
    generic            :: set_all2single => set_all2single_1, set_all2single_2
    procedure          :: set_projs
    procedure          :: e3swapsgn
    procedure          :: swape1e3
    procedure          :: zero
    procedure          :: zero_projs
    procedure          :: zero_shifts
    procedure          :: mul_shifts
    procedure          :: rnd_oris
    procedure          :: rnd_proj_space
    procedure          :: rnd_neighbors
    procedure          :: rnd_gau_neighbors
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
    procedure          :: ini_tseries
    procedure          :: symmetrize
    procedure          :: merge
    procedure          :: partition_eo
    procedure          :: transf_proj2class
    ! I/O
    procedure          :: read
    procedure          :: read_ctfparams_state_eo
    procedure, private :: write_1
    procedure, private :: write_2
    generic            :: write => write_1, write_2
    ! CALCULATORS
    procedure          :: compress
    procedure          :: split_state
    procedure          :: split_class
    procedure          :: expand_classes
    procedure          :: fill_empty_classes
    procedure          :: remap_classes
    procedure          :: shift_classes
    procedure          :: create_conforming_npeaks_set
    procedure          :: merge_classes
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
    procedure          :: spiral_1
    procedure          :: spiral_2
    generic            :: spiral => spiral_1, spiral_2
    procedure          :: qspiral
    procedure          :: order
    procedure          :: order_corr
    procedure          :: order_cls
    procedure          :: reduce_projs
    procedure, private :: balance_1
    procedure, private :: balance_2
    generic            :: balance => balance_1, balance_2
    procedure          :: calc_hard_weights
    procedure          :: calc_spectral_weights
    procedure          :: reject_above
    procedure          :: find_closest_proj
    procedure          :: find_closest_projs
    procedure          :: find_closest_ori
    procedure          :: find_closest_oris
    procedure          :: calc_euldists
    procedure, private :: create_proj_subspace_1
    procedure, private :: create_proj_subspace_2
    generic            :: create_proj_subspace => create_proj_subspace_1, create_proj_subspace_2
    procedure          :: discretize
    procedure          :: nearest_proj_neighbors
    procedure          :: find_angres
    procedure          :: find_angres_geod
    procedure          :: extremal_bound
    procedure          :: find_npeaks
    procedure          :: find_athres_from_npeaks
    procedure          :: find_npeaks_from_athres
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
    ! DESTRUCTORS
    procedure          :: kill_chash
    procedure          :: kill
end type oris

interface oris
    module procedure constructor
end interface oris

type(oris), pointer  :: ops  =>null()
logical, allocatable :: class_part_of_set(:)
real,    allocatable :: class_weights(:)
type(ori)            :: o_glob

contains

    ! CONSTRUCTORS

    !>  \brief  is an abstract constructor
    function constructor( n ) result( self )
        integer, intent(in) :: n
        type(oris) :: self
        call self%new(n)
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        integer :: i
        call self%kill
        self%n = n
        allocate( self%o(self%n), stat=alloc_stat )
        allocchk('new; simple_oris')
        do i=1,n
            call self%o(i)%new
        end do
    end subroutine new

    !>  \brief  is a constructor
    subroutine new_clean( self, n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        integer :: alloc_stat, i
        call self%kill
        self%n = n
        allocate( self%o(self%n), stat=alloc_stat )
        allocchk('new_clean; simple_oris')
        do i=1,n
            call self%o(i)%new_ori_clean
        end do
    end subroutine new_clean

    ! GETTERS

    !>  \brief  returns the i:th first Euler angle
    pure function e1get( self, i ) result( e1 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e1
        e1 = self%o(i)%e1get()
    end function e1get

    !>  \brief  returns the i:th second Euler angle
    pure function e2get( self, i ) result( e2 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e2
        e2 = self%o(i)%e2get()
    end function e2get

    !>  \brief  returns the i:th third Euler angle
    pure function e3get( self, i ) result( e3 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e3
        e3 = self%o(i)%e3get()
    end function e3get

    !>  \brief  returns the i:th Euler triplet
    pure function get_euler( self, i ) result( euls )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: euls(3)
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
        integer,     intent(in)    :: i
        type(ori) :: o
        if( self%n == 0 ) stop 'oris object does not exist; get_ori; simple_oris'
        if( i > self%n .or. i < 1 )then
            write(*,*) 'trying to get ori: ', i, ' among: ', self%n, ' oris'
            stop 'i out of range; get_ori; simple_oris'
        endif
        o = self%o(i)
    end function get_ori

    !>  \brief  is a multiparameter getter
    function get( self, i, key ) result( val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real :: val
        val = self%o(i)%get(key)
    end function get

    !>  \brief  is a getter
    subroutine getter_1( self, i, key, val )
        class(oris),                   intent(inout) :: self
        integer,                       intent(in)    :: i
        character(len=*),              intent(in)    :: key
        character(len=:), allocatable, intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_1

    ! !>  \brief  is a getter
    subroutine getter_2( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_2

    !>  \brief  is for getting an array of 'key' values
    function get_all( self, key, fromto ) result( arr )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: key
        integer, optional, intent(in)    :: fromto(2)
        real, allocatable :: arr(:)
        integer :: i, ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        allocate( arr(ffromto(1):ffromto(2)), stat=alloc_stat)
        allocchk('get_all; simple_oris')
        do i=ffromto(1),ffromto(2)
            arr(i) = self%o(i)%get(key)
        enddo
    end function get_all

    !>  \brief  is for getting an array of 'key' values
    subroutine getter_all_1( self, key, vals )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: key
        real, allocatable, intent(out)   :: vals(:)
        integer :: i
        if( allocated(vals) ) deallocate(vals)
        allocate( vals(self%n), stat=alloc_stat)
        allocchk('getter_all_1; simple_oris')
        do i=1,self%n
            call self%o(i)%getter(key, vals(i))
        enddo
    end subroutine getter_all_1

    !>  \brief  is for getting an array of 'key' values
    subroutine getter_all_2( self, key, vals )
        class(oris),                        intent(inout) :: self
        character(len=*),                   intent(in)    :: key
        character(len=STDLEN), allocatable, intent(out)   :: vals(:)
        integer :: i
        character(len=:), allocatable :: tmp
        if( allocated(vals) ) deallocate(vals)
        allocate( vals(self%n), stat=alloc_stat)
        allocchk('getter_all_2; simple_oris')
        do i=1,self%n
            call self%o(i)%getter(key, tmp)
            vals(i) = tmp
        enddo
    end subroutine getter_all_2

    !>  \brief  is for getting the i:th rotation matrix
    pure function get_mat( self, i ) result( mat )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: mat(3,3)
        mat = self%o(i)%get_mat()
    end function get_mat

    !>  \brief  is for getting the i:th normal
    pure function get_normal( self, i ) result( normal )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: normal(3)
        normal = self%o(i)%get_normal()
    end function get_normal

    !>  \brief  is a setter
    function get_2Dshift( self, i )  result(shvec)
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real :: shvec(2)
        shvec = self%o(i)%get_2Dshift()
    end function get_2Dshift

    !>  \brief  is for checking if parameter is present
    function isthere_1( self, key ) result( is )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        logical :: is
        integer :: i
        is = .false.
        do i=1,self%n
            is = self%o(i)%isthere(key)
            if( is ) exit
        end do
    end function isthere_1

    !>  \brief  is for checking if parameter is present
    function isthere_2( self, i, key ) result( is )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        logical :: is
        is = self%o(i)%isthere(key)
    end function isthere_2

    !>  \brief  is for getting the maximum hash size
    integer function max_hash_size( self )
        class(oris),      intent(inout) :: self
        integer :: sz, i
        max_hash_size = 0
        do i=1,self%n
            sz = self%o(i)%hash_size()
            if( sz > max_hash_size ) max_hash_size = sz
        end do
    end function max_hash_size

    !>  \brief  is for getting the max val of integer label
    function get_n( self, label ) result( n )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: label
        integer :: i, n, ival
        n = 1
        do i=1,self%n
            ival = nint(self%o(i)%get(trim(label)))
            if( ival > n ) n = ival
        end do
    end function get_n

    !>  \brief  is for checking label population
    function get_pop( self, ind, label, consider_w, eo ) result( pop )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: ind
        character(len=*),  intent(in)    :: label
        logical, optional, intent(in)    :: consider_w
        integer, optional, intent(in)    :: eo
        integer :: mylab, pop, i, mystate, myeo
        logical :: cconsider_w, consider_eo
        real    :: w
        cconsider_w = .false.
        if( present(consider_w) ) cconsider_w = consider_w
        if( cconsider_w )then
            if( .not. self%isthere('w') ) stop 'ERROR, oris :: get_pop with optional consider_w assumes w set'
        endif
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        pop = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            myeo    = nint(self%o(i)%get('eo'))
            if( consider_eo )then
                if( myeo /= eo ) cycle
            endif
            w = 1.0
            if( cconsider_w ) w = self%o(i)%get('w')
            if( mystate > 0 .and. w > TINY )then
                mylab = nint(self%o(i)%get(label))
                if( mylab == ind )  pop = pop + 1
            endif
        end do
    end function get_pop

    ! !>  \brief  is for checking discrete label populations
    ! function get_pops( self, label, consider_w, maxn, eo ) result( pops )
    !     class(oris),       intent(inout) :: self
    !     character(len=*),  intent(in)    :: label
    !     logical, optional, intent(in)    :: consider_w
    !     integer, optional, intent(in)    :: maxn ! max label, for the case where the last class/state is missing
    !     integer, optional, intent(in)    :: eo
    !     integer, allocatable :: pops(:)
    !     integer :: i, mystate, myval, n, myeo
    !     real    :: w
    !     logical :: cconsider_w, consider_eo
    !     cconsider_w = .false.
    !     if( present(consider_w) ) cconsider_w = consider_w
    !     if( cconsider_w )then
    !         if( .not. self%isthere('w') ) stop 'ERROR, oris :: get_pops with optional consider_w assumes w set'
    !     endif
    !     consider_eo = .false.
    !     if( present(eo) ) consider_eo = .true.
    !     n = self%get_n(label)
    !     if( present(maxn) )then
    !         n = max(n, maxn)
    !     endif
    !     if(allocated(pops))deallocate(pops)
    !     allocate(pops(n),source=0,stat=alloc_stat)
    !     allocchk('In: get_pops, module: simple_oris')
    !     do i=1,self%n
    !         mystate = nint(self%o(i)%get('state'))
    !         myeo    = nint(self%o(i)%get('eo'))
    !         if( consider_eo )then
    !             if( myeo /= eo ) cycle
    !         endif
    !         w = 1.0
    !         if( cconsider_w )  w = self%o(i)%get('w')
    !         if( mystate > 0 .and. w > TINY )then
    !             myval = nint(self%o(i)%get(label))
    !             if( myval > 0 ) pops(myval) = pops(myval) + 1
    !         endif
    !     end do
    ! end function get_pops

    !>  \brief  is for checking discrete label populations
    subroutine get_pops( self, pops, label, consider_w, maxn, eo )
        class(oris),          intent(inout) :: self
        integer, allocatable, intent(out)   :: pops(:)
        character(len=*),     intent(in)    :: label
        logical, optional,    intent(in)    :: consider_w
        integer, optional,    intent(in)    :: maxn ! max label, for the case where the last class/state is missing
        integer, optional,    intent(in)    :: eo
        integer :: i, mystate, myval, n, myeo
        real    :: w
        logical :: cconsider_w, consider_eo
        cconsider_w = .false.
        if( present(consider_w) ) cconsider_w = consider_w
        if( cconsider_w )then
            if( .not. self%isthere('w') ) stop 'ERROR, oris :: get_pops with optional consider_w assumes w set'
        endif
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        n = self%get_n(label)
        if( present(maxn) )then
            n = max(n, maxn)
        endif
        if(allocated(pops))deallocate(pops)
        allocate(pops(n),source=0,stat=alloc_stat)
        allocchk('In: get_pops, module: simple_oris')
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            myeo    = nint(self%o(i)%get('eo'))
            if( consider_eo )then
                if( myeo /= eo ) cycle
            endif
            w = 1.0
            if( cconsider_w )  w = self%o(i)%get('w')
            if( mystate > 0 .and. w > TINY )then
                myval = nint(self%o(i)%get(label))
                if( myval > 0 ) pops(myval) = pops(myval) + 1
            endif
        end do
    end subroutine get_pops

    !>  \brief  returns a logical array of state existence
    function states_exist(self, nstates) result(exists)
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nstates
        integer :: i
        logical :: exists(nstates)
        do i=1,nstates
            exists(i) = (self%get_pop(i, 'state') > 0)
        end do
    end function states_exist

    !>  \brief  compresses the aggregated object according to input mask
    subroutine compress( self, mask )
        class(oris), intent(inout) :: self
        logical,     intent(in)    :: mask(:)
        type(oris) :: os_tmp
        integer    :: i, cnt
        if( size(mask) /= self%n )then
            print *, 'self%n:     ', self%n
            print *, 'size(mask): ', size(mask)
            stop 'nonconforming mask size; oris :: compress'
        endif
        call os_tmp%new(count(mask))
        cnt = 0
        do i=1,self%n
            if( mask(i) )then
                cnt = cnt + 1
                os_tmp%o(cnt) = self%o(i)
            endif 
        end do
        self = os_tmp
        call os_tmp%kill
    end subroutine compress

    !>  \brief  for balanced split of a state group
    subroutine split_state( self, which )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: which
        integer, allocatable :: ptcls_in_which(:)
        integer, allocatable :: states(:)
        type(ran_tabu)       :: rt
        integer              ::  n, nstates, iptcl
        nstates = self%get_n('state')
        if( which < 1 .or. which > nstates )then
            stop 'which (state) is out of range; simple_oris::split_state'
        endif
        ptcls_in_which = self%get_pinds(which, 'state')
        n = size(ptcls_in_which)
        allocate(states(n), stat=alloc_stat)
        allocchk("simple_oris::split_state")
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
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: which
        integer, allocatable :: ptcls_in_which(:)
        integer, allocatable :: members(:)
        type(ran_tabu)       :: rt
        integer              ::  n, nmembers, iptcl
        nmembers = self%get_n('class')
        if( which < 1 .or. which > nmembers )then
            stop 'which member is out of range; simple_oris :: split_class'
        endif
        ptcls_in_which = self%get_pinds(which, 'class')
        n = size(ptcls_in_which)
        allocate(members(n), stat=alloc_stat)
        allocchk("simple_oris::split_class")
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
        integer,     intent(in)    :: ncls_target
        integer, allocatable       :: pops(:)
        integer :: ncls, loc(1), myncls, icls
        ncls = self%get_n('class')
        if( ncls_target <= ncls ) stop 'Number of target classes cannot be <= original number'
        ! calculate class populations
        allocate(pops(ncls_target),stat=alloc_stat)
        allocchk('In: expand_classes, module: simple_oris')
        pops = 0
        do icls=1,ncls
            pops(icls) = self%get_pop(icls, 'class')
        end do
        myncls = ncls
        do while( myncls < ncls_target )
            ! split largest
            loc = maxloc(pops)
            call self%split_class(loc(1))
            ! update number of classes
            myncls = myncls+1
            ! update pops
            pops(loc(1)) = self%get_pop(loc(1), 'class')
            pops(myncls) = self%get_pop(myncls, 'class')
        end do
    end subroutine expand_classes

    !>  \brief  is for filling empty classes from the highest populated ones
    subroutine fill_empty_classes( self, ncls, fromtocls )
        class(oris),                    intent(inout) :: self
        integer,                        intent(in)    :: ncls
        integer, allocatable, optional, intent(out)   :: fromtocls(:,:)
        integer, allocatable :: inds2split(:), membership(:), pops(:), fromtoall(:,:)
        type(ran_tabu)       :: rt
        integer              :: iptcl, cls2split(1)
        integer              :: cnt, nempty, n_incl, maxpop, i, icls, ncls_here
        logical              :: consider_w_here
        ncls_here = max(self%get_n('class'), ncls)
        call self%get_pops(pops, 'class', consider_w=.true., maxn=ncls_here)            
        nempty    = count(pops == 0)
        if(nempty == 0)return
        if( present(fromtocls) )then
            if(allocated(fromtocls))deallocate(fromtocls)
            allocate(fromtoall(ncls,2), source=0, stat=alloc_stat)
            allocchk("simple_oris:: fill_empty_classes fromtoall")
        endif
        do icls = 1, ncls
            deallocate(pops, stat=alloc_stat)
            call self%get_pops(pops, 'class', consider_w=.true., maxn=ncls_here)            
            if( pops(icls) /= 0 )cycle
            ! identify class to split
            cls2split  = maxloc(pops)
            maxpop     = pops(cls2split(1))
            if( maxpop <= 2*MINCLSPOPLIM )exit
            ! migration
            inds2split = self%get_pinds(cls2split(1), 'class', consider_w=.false.) ! needs to be false
            rt = ran_tabu(maxpop)
            allocate(membership(maxpop), stat=alloc_stat)
            allocchk("simple_oris:: fill_empty_classes membership")
            call rt%balanced(2, membership)
            do i=1,maxpop
                if(membership(i) == 2)cycle
                iptcl = inds2split(i)
                call self%o(iptcl)%set('class', real(icls))
            enddo
            ! updates populations and migration
            if(present(fromtocls))then
                fromtoall(icls,1) = cls2split(1)
                fromtoall(icls,2) = icls
            endif
            deallocate(membership, inds2split)
        enddo
        if(present(fromtocls))then
            ! updates classes migration
            n_incl = count(fromtoall(:,1)>0)
            if(n_incl>0)then
                allocate(fromtocls(n_incl,2), stat=alloc_stat)
                allocchk("simple_oris:: fill_empty_classes fromtocls")
                cnt = 0
                do icls = 1,ncls
                    if( fromtoall(icls,1) > 0 )then
                        cnt = cnt + 1
                        fromtocls(cnt,:) = fromtoall(icls,:)
                    endif
                enddo
            endif
        endif
    end subroutine fill_empty_classes

    !>  \brief  for remapping classes after exclusion
    subroutine remap_classes( self )
        class(oris), intent(inout) :: self
        integer :: ncls, clsind_remap, pop, icls, iptcl, old_cls
        integer , allocatable :: clspops(:)
        ncls = self%get_n('class')
        allocate(clspops(ncls),stat=alloc_stat)
        allocchk('In: remap_classes, module: simple_oris')
        do icls=1,ncls
            clspops(icls) = self%get_pop(icls, 'class')
        end do
        if( any(clspops == 0) )then
            clsind_remap = ncls
            do icls=1,ncls
                pop = clspops(icls)
                if( pop > 1 )then
                    clsind_remap = clsind_remap + 1
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

    !>  \brief  for shifting class indices after chunk-based prime2D exec
    !!             1 .. 10
    !!    + 10 => 11 .. 20
    !!    + 10 => 21 .. 31 etc.
    subroutine shift_classes( self, ishift )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ishift
        integer :: iptcl, new_cls, old_cls
        do iptcl=1,self%n
            old_cls = nint(self%o(iptcl)%get('class'))
            new_cls = old_cls + ishift
            call self%o(iptcl)%set('class', real(new_cls))
        end do
    end subroutine shift_classes

    !>  \brief  is for contracting/expanding the ori set in self to make 
    !!          a set of conforming size
    function create_conforming_npeaks_set( self, npeaks ) result( os_conf )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: npeaks
        real, allocatable :: ows(:)
        type(ori)         :: o_template
        type(oris)        :: os_conf
        integer           :: inds(self%n), i, cnt
        if( .not. self%isthere('ow') )&
        &stop 'orientation weights (ow) needs to be set; oris :: create_conforming_npeaks_set'
        if( self%n == npeaks )then
            ! set is conforming
            os_conf = self
        else if( self%n > npeaks )then
            ! shrink set
            ows  = self%get_all('ow')
            inds = (/(i,i=1,self%n)/)
            call hpsort(ows, inds)
            cnt = 0
            do i=self%n,self%n - npeaks + 1,-1
                cnt = cnt + 1
                os_conf%o(cnt) = self%o(inds(i))
            end do
        else
            ! expand set
            o_template = self%o(1)
            call o_template%set_euler([0.,0.,0.])
            call o_template%set_shift([0.,0.])
            call o_template%set('state', 0.0)
            call o_template%set('ow', 0.0)
            call os_conf%new(npeaks)
            do i=1,npeaks
                if( i <= self%n )then
                    os_conf%o(i) = self%o(i)
                else
                    os_conf%o(i) = o_template
                endif
            end do
        endif
    end function create_conforming_npeaks_set

    !>  \brief  is for getting an allocatable array with ptcl indices of the label 'label'
    function get_pinds( self, ind, label, consider_w ) result( indices )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: label
        integer,           intent(in)    :: ind
        logical, optional, intent(in)    :: consider_w
        integer, allocatable :: indices(:)
        integer :: pop, mystate, cnt, myval, i
        logical :: cconsider_w
        real    :: w
        cconsider_w = .false.
        if( present(consider_w) ) cconsider_w = consider_w
        if( cconsider_w )then
            if( .not. self%isthere('w') ) stop 'ERROR, oris :: get_pinds with optional consider_w assumes w set'
        endif
        pop = self%get_pop(ind, label, cconsider_w ) 
        if( pop > 0 )then
            allocate( indices(pop), stat=alloc_stat )
            allocchk('get_pinds; simple_oris')
            cnt = 0
            do i=1,self%n
                mystate = nint(self%o(i)%get('state'))
                w = 1.0
                if( cconsider_w ) w = self%o(i)%get('w')
                if( mystate > 0 .and. w > TINY )then
                    myval = nint(self%get(i, trim(label)))
                    if( myval == ind )then
                        cnt = cnt + 1
                        indices(cnt) = i
                    endif
                endif
            end do
        endif
    end function get_pinds

    !>  \brief  generate a mask with the oris with mystate == state/ind == get(label)
    subroutine gen_mask( self, state, ind, label, l_mask, consider_w, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state, ind
        character(len=*),     intent(in)    :: label
        logical, allocatable, intent(out)   :: l_mask(:)
        logical, optional,    intent(in)    :: consider_w
        integer, optional,    intent(in)    :: fromto(2)
        logical :: cconsider_w
        real    :: w
        integer :: i, mystate, myval, ffromto(2)
        cconsider_w = .false.
        if( present(consider_w) ) cconsider_w = consider_w
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        if( allocated(l_mask) ) deallocate(l_mask)
        allocate(l_mask(ffromto(1):ffromto(2)))
        l_mask = .false.
        do i=ffromto(1),ffromto(2)
            w = 1.0
            if( cconsider_w ) w = self%o(i)%get('w')
            mystate = nint(self%o(i)%get('state'))
            myval   = nint(self%o(i)%get(trim(label)))
            if( mystate == state .and. (myval == ind .and. w > TINY) ) l_mask(i) = .true.
        end do
    end subroutine gen_mask

    !>  \brief  is for getting an array of 'which' variables with
    !!          filtering based on class/state/proj
    function get_arr( self, which, class, state ) result( vals )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        real, allocatable :: vals(:)
        integer :: pop, cnt, clsnr, i, mystate
        real    :: val
        logical :: class_present, state_present
        class_present = class_present
        state_present = present(state)
        if( class_present )then
            pop = self%get_pop(class, 'class')
        else if( state_present )then
            pop = self%get_pop(state, 'state')
        else
            pop = self%n
        endif
        if( pop > 0 )then
            allocate( vals(pop), stat=alloc_stat )
            allocchk('get_arr; simple_oris')
            cnt = 0
            do i=1,self%n
                val = self%get(i, which)
                if( class_present )then
                    mystate = nint(self%get( i, 'state'))
                    if( mystate > 0 )then
                        clsnr   = nint(self%get( i, 'class'))
                        if( clsnr == class )then
                            cnt = cnt+1
                            vals(cnt) = val
                        endif
                    endif
                else if( state_present )then
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
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: sum
        integer,           intent(out)   :: cnt
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        integer :: clsnr, i, mystate, istart, istop
        real    :: val
        logical :: proceed, class_present, state_present, mask_present
        class_present = present(class)
        state_present = present(state)
        mask_present  = present(mask)
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
            if( mask_present )then
                proceed = mask(i)
            else
                proceed = .true.
            endif
            if( proceed )then
                val = self%get(i, which)
                if( .not. is_a_number(val) ) val=0.
                if( class_present )then
                    clsnr = nint(self%get( i, 'class'))
                    if( clsnr == class )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else if( state_present )then
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
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        integer :: cnt
        real    :: sum
        call self%calc_sum(which, sum, cnt, class, state, fromto, mask)
    end function get_sum

    !>  \brief  is for getting the average of 'which' variables with
    !!          filtering based on class/state/fromto
    function get_avg( self, which, class, state, fromto, mask ) result( avg )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        integer :: cnt
        real    :: avg, sum
        call self%calc_sum(which, sum, cnt, class, state, fromto, mask)
        avg = sum/real(cnt)
    end function get_avg

    !>  \brief  is for calculating the nonzero sum of 'which' variables with
    !!          filtering based on class/state/fromto
    subroutine calc_nonzero_sum( self, which, sum, cnt, class, state, fromto )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: sum
        integer,           intent(out)   :: cnt
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        integer :: clsnr, i, mystate, istart, istop
        real    :: val
        logical :: class_present, state_present
        class_present = present(class)
        state_present = present(state)

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
            if( .not. is_a_number(val) ) val = 0.
            if( val > 0. )then
                if( class_present )then
                    clsnr = nint(self%get( i, 'class'))
                    if( clsnr == class )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else if( state_present )then
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
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        integer :: cnt
        real    :: sum
        call self%calc_nonzero_sum(which, sum, cnt, class, state, fromto)
    end function get_nonzero_sum

    !>  \brief  is for getting the average of 'which' variables with
    !!          filtering based on class/state/fromto
    function get_nonzero_avg( self, which, class, state, fromto ) result( avg )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        integer :: cnt
        real    :: avg, sum
        call self%calc_nonzero_sum(which, sum, cnt, class, state, fromto)
        avg = sum/real(cnt)
    end function get_nonzero_avg

    !>  \brief  is for extracting a table from the character hash
    function extract_table( self, which ) result( table )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        character(len=STDLEN), allocatable :: table(:)
        character(len=:),      allocatable :: str
        integer :: i
        allocate(table(self%n), stat=alloc_stat)
        allocchk('In: extract_table, module: simple_oris')
        do i=1,self%n
            call self%getter(i, which, str)
            table(i) = str
        end do
    end function extract_table

    !>  \brief  angular standard deviation
    real function ang_sdev( self, refine, nstates, npeaks )
        class(oris),           intent(inout) :: self
        character(len=STDLEN), intent(in)    :: refine
        integer,               intent(in)    :: nstates, npeaks
        integer :: instates, state, pop
        ang_sdev = 0.
        if(npeaks < 3 .or. trim(refine).eq.'shc' .or. trim(refine).eq.'shcneigh')return
        instates = 0
        do state=1,nstates
            pop = self%get_pop( state, 'state' )
            if( pop > 0 )then
                ang_sdev = ang_sdev + ang_sdev_state( state )
                instates = instates + 1
            endif
        enddo
        ang_sdev = ang_sdev / real(instates)

        contains

            function ang_sdev_state( istate )result( isdev )
                integer, intent(in)  :: istate
                type(ori)            :: o_best, o
                type(oris)           :: os
                real,    allocatable :: dists(:), ws(:)
                integer, allocatable :: inds(:)
                real                 :: ave, isdev, var
                integer              :: loc(1), i, ind, n, cnt
                logical              :: err
                isdev = 0.
                inds = self%get_pinds( istate, 'state' )
                n = size(inds)
                if( n < 3 )return ! because one is excluded in the next step & moment needs at least 2 objs
                call os%new( n )
                allocate( ws(n), dists(n-1), stat=alloc_stat )
                allocchk( 'ang_sdev_state; simple_oris')
                ws    = 0.
                dists = 0.
                ! get best ori
                do i=1,n
                    ind   = inds(i)
                    ws(i) = self%get(ind,'ow')
                    call os%set_ori( i, self%o(ind) )
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

    !>  \brief  for getting a logical mask of the included particles
    function included( self, consider_w ) result( incl )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: consider_w
        logical, allocatable :: incl(:)
        integer :: i, istate
        logical :: cconsider_w
        real    :: w
        cconsider_w = .false.
        if( present(consider_w) ) cconsider_w = consider_w
        if( cconsider_w )then
            if( .not. self%isthere('w') ) stop 'ERROR, oris :: included with optional consider_w assumes w set'
        endif
        allocate(incl(self%n), stat=alloc_stat)
        allocchk('In: included, module: simple_oris')
        incl = .false.
        do i=1,self%n
            istate = nint(self%o(i)%get('state'))
            w = 1.0
            if( cconsider_w ) w = self%o(i)%get('w')
            if( istate > 0 .and. w > TINY ) incl(i) = .true.
        end do
    end function included

    !>  \brief  is getting the number of oris assigned to the even partion
    ! integer function get_neven( self, fromto )
    !     class(oris),       intent(inout) :: self
    !     integer, optional, intent(in)    :: fromto(2)
    !     integer :: i, from, to
    !     from = 1
    !     to   = self%n
    !     if(present(fromto))then
    !         from = fromto(1)
    !         to   = fromto(2)
    !     endif
    !     get_neven = 0
    !     do i = from, to
    !         if( self%o(i)%isthere('eo') )then
    !             if( self%o(i)%get('state') < 0.5 )cycle
    !             if( self%o(i)%get('weights') < TINY )cycle
    !             if( self%o(i)%iseven() )get_neven = get_neven + 1
    !         endif
    !     enddo
    ! end function get_neven

    integer function get_neven( self )
        class(oris), intent(inout) :: self
        integer, allocatable :: eopart(:)
        eopart = nint(self%get_all('eo'))
        get_neven = count(eopart == 0)
    end function get_neven

    !>  \brief  is getting the number of oris assigned to the odd partion
    ! integer function get_nodd( self, fromto )
    !     class(oris),       intent(inout) :: self
    !     integer, optional, intent(in)    :: fromto(2)
    !     integer :: i, from, to
    !     from = 1
    !     to   = self%n
    !     if(present(fromto))then
    !         from = fromto(1)
    !         to   = fromto(2)
    !     endif
    !     get_nodd = 0
    !     do i = from, to
    !         if( self%o(i)%isthere('eo') )then
    !             if( self%o(i)%get('state') < 0.5 )cycle
    !             if( self%o(i)%get('weights') < TINY )cycle
    !             if( self%o(i)%isodd() )get_nodd = get_nodd + 1
    !         endif
    !     enddo
    ! end function get_nodd

    !>  \brief  is getting the number of oris assigned to the odd partion
    integer function get_nodd( self )
        class(oris), intent(inout) :: self
        integer, allocatable :: eopart(:)
        eopart = nint(self%get_all('eo'))
        get_nodd = count(eopart == 1)
    end function get_nodd

    !>  \brief  is getting the number of oris assigned to the odd partion
    integer function get_nevenodd( self, fromto )
        class(oris),       intent(inout) :: self
        integer, optional, intent(in)    :: fromto(2)
        get_nevenodd = self%get_neven() + self%get_nodd()
    end function get_nevenodd

    !>  \brief  is for printing
    subroutine print_( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        call self%o(i)%print_ori()
    end subroutine print_

    !>  \brief  is for printing
    subroutine print_matrices( self )
        class(oris), intent(inout) :: self
        integer :: i
        write(*,*) 'ORDER OF ROTATION MATRIX ELEMENTS: (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)'
        do i=1,self%n
            call self%o(i)%print_mat()
        end do
    end subroutine print_matrices

    subroutine sample4update_and_incrcnt( self, fromto, update_frac, nsamples, inds, mask )
        use simple_rnd, only: multinomal_many
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(out)   :: nsamples
        integer, allocatable, intent(out)   :: inds(:)
        logical,              intent(out)   :: mask(fromto(1):fromto(2))
        real,    allocatable :: counts(:), states(:)
        integer, allocatable :: inds_here(:)
        integer :: i, cnt, n_incl
        real    :: val
        ! gather info
        allocate( states(fromto(1):fromto(2)), counts(fromto(1):fromto(2)), inds_here(fromto(1):fromto(2)))
        do i=fromto(1),fromto(2)
            if( self%o(i)%isthere('updatecnt') )then
                counts(i) = self%o(i)%get('updatecnt')
            else
                counts(i) = 0
            endif
            states(i)    = self%o(i)%get('state')
            inds_here(i) = i
        end do
        ! order counts
        call hpsort(counts, inds_here)
        ! figure out how many samples
        nsamples = nint(update_frac * real(count(states > 0.5)))
        ! allocate output index array
        if( allocated(inds) ) deallocate(inds)
        allocate( inds(nsamples) )
        ! sample
        cnt = 0
        do i=fromto(1),fromto(2)
            if( states(i) > 0.5 )then
                cnt = cnt + 1
                inds(cnt) = inds_here(i)
            endif
            if( cnt == nsamples ) exit
        end do
        ! increment update counter and set mask
        mask = .false.
        do i=1,nsamples
            val = self%o(inds(i))%get('updatecnt')
            val = val + 1.0
            call self%o(inds(i))%set('updatecnt', val)
            mask(inds(i)) = .true.
        end do 
    end subroutine sample4update_and_incrcnt

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
    subroutine reject( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        call self%o(i)%reject
    end subroutine reject

    !>  \brief  is a setter
    subroutine set_euler( self, i, euls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: euls(3)
        call self%o(i)%set_euler(euls)
    end subroutine set_euler

    !>  \brief  is a setter
    subroutine set_shift( self, i, vec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: vec(2)
        call self%o(i)%set_shift(vec)
    end subroutine set_shift

    !>  \brief  is a setter
    subroutine e1set( self, i, e1 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e1
        call self%o(i)%e1set(e1)
    end subroutine e1set

    !>  \brief  is a setter
    subroutine e2set( self, i, e2 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e2
        call self%o(i)%e2set(e2)
    end subroutine e2set

    !>  \brief  is a setter
    subroutine e3set( self, i, e3 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e3
        call self%o(i)%e3set(e3)
    end subroutine e3set

    !>  \brief  is a setter
    subroutine set_1( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_1

    !>  \brief  is a setter
    subroutine set_2( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_2

    !>  \brief  is a setter
    subroutine set_ori( self, i, o )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: o
        self%o(i) = o
    end subroutine set_ori

    !>  \brief  is a setter
    subroutine set_all_1( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_1

    !>  \brief  is a setter
    subroutine set_all_2( self, which, vals )
        class(oris),           intent(inout) :: self
        character(len=*),      intent(in)    :: which
        character(len=STDLEN), intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_2

    !>  \brief  is a setter
    subroutine set_all2single_1( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, val)
        end do
    end subroutine set_all2single_1

    !>  \brief  is a setter
    subroutine set_all2single_2( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, val)
        end do
    end subroutine set_all2single_2

    !>  \brief  is a setter
    subroutine set_projs( self, e_space )
        class(oris), intent(inout) :: self
        class(oris), intent(inout) :: e_space
        integer    :: i
        type(ori)  :: o
        do i=1,self%n
            o = self%get_ori(i)
            call self%set(i, 'proj', real(e_space%find_closest_proj(o)))
        end do
    end subroutine set_projs

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
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
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
        real,        intent(in)    :: mul
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', mul*self%o(i)%get('x'))
            call self%o(i)%set('y', mul*self%o(i)%get('y'))
        end do
    end subroutine mul_shifts

    !>  \brief  randomizes eulers in oris
    subroutine rnd_oris( self, trs, eullims )
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
        integer :: i
        do i=1,self%n
            call self%o(i)%rnd_ori(trs, eullims)
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

    !>  \brief  generates nnn stochastic projection direction neighbors to o_prev
    !!          with a sigma of the random Gaussian tilt angle of asig
    subroutine rnd_gau_neighbors( self, nnn, o_prev, asig, eullims )
        class(oris),    intent(inout) :: self
        integer,        intent(in)    :: nnn
        class(ori),     intent(in)    :: o_prev
        real,           intent(in)    :: asig
        real, optional, intent(inout) :: eullims(3,2)
        type(ori) :: o_transform, o_rnd
        real      :: val
        integer   :: i
        call o_transform%new
        call self%new(nnn)
        do i=1,nnn
            call o_transform%rnd_euler(eullims)
            val = gasdev(0., asig)
            do while(abs(val) > eullims(2,2))
                val = gasdev(0., asig)
            enddo
            call o_transform%e2set(val)
            call o_transform%e3set(0.)
            o_rnd = o_prev.compose.o_transform
            call self%set_ori(i, o_rnd)
        end do
    end subroutine rnd_gau_neighbors

    !>  \brief  generate random projection direction space around a given one
    subroutine rnd_proj_space( self, nsample, o_prev, thres, eullims )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: nsample      !< # samples
        class(ori), optional, intent(inout) :: o_prev       !< orientation
        real,       optional, intent(inout) :: eullims(3,2) !< Euler limits
        real,       optional, intent(in)    :: thres        !< half-angle of spherical cap
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
                    if( rad2deg( o_stoch.euldist.o_prev ) > thres )cycle
                    found = .true.
                    call o_stoch%e3set( 0.)
                    self%o( i ) = o_stoch
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
        class(oris),    intent(inout) :: self
        integer,        intent(in)    :: i
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
        call self%o(i)%rnd_ori( trs, eullims )
    end subroutine rnd_ori

    !>  \brief  for generating random clustering
    subroutine rnd_cls( self, ncls, srch_inpl )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: ncls
        logical, optional, intent(in)    :: srch_inpl
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

    !>  \brief  for generating an initial clustering of time series
    subroutine ini_tseries( self, nsplit, state_or_class )
        use simple_map_reduce, only: split_nobjs_even
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: nsplit
        character(len=*), intent(in)    :: state_or_class
        integer, allocatable :: parts(:,:)
        integer :: ipart, iptcl
        parts = split_nobjs_even(self%n, nsplit)
        do ipart=1,nsplit
            do iptcl=parts(ipart,1),parts(ipart,2)
                call self%o(iptcl)%set(trim(state_or_class), real(ipart))
            end do
        end do
        deallocate(parts)
    end subroutine ini_tseries

    !>  \brief  randomizes eulers in oris
    subroutine rnd_oris_discrete( self, discrete, nsym, eullims )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: discrete, nsym
        real,        intent(in)    :: eullims(3,2)
        type(oris) :: tmp
        real       :: euls(3)
        integer    :: i, irnd
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
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
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
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nstates
        integer, allocatable       :: states(:)
        type(ran_tabu)             :: rt
        integer :: i, state
        if( nstates > 1 )then
            allocate( states(self%n), stat=alloc_stat )
            allocchk("simple_oris::rnd_states")
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
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls
        integer, allocatable       :: classes(:)
        type(ran_tabu)             :: rt
        integer :: i
        if( ncls > 1 )then
            allocate(classes(self%n), stat=alloc_stat)
            allocchk("simple_oris::rnd_classes")
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
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('lp', ran3()*100.)
        end do
    end subroutine rnd_lps

    !>  \brief  randomizes correlations in oris
    subroutine rnd_corrs( self )
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
        integer,     intent(in)    :: class_merged, class
        integer                    :: i, clsnr
        do i=1,self%n
            clsnr = nint(self%get(i, 'class'))
            if(clsnr == class) call self%set(i, 'class', real(class_merged))
        end do
    end subroutine merge_classes

    !>  \brief  for extending algndoc according to nr of symmetry ops
    subroutine symmetrize( self, nsym )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsym
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
        type(oris) :: self
        integer    :: ntot, cnt, i
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

    !>  \brief  for balanced assignment of even/odd partitions
    subroutine partition_eo( self, tseries )
        class(oris),       intent(inout) :: self    !< instance
        logical, optional, intent(in)    :: tseries !< logical tseries flag
        type(ran_tabu)       :: rt
        integer, allocatable :: eopart(:)
        logical              :: ttseries
        integer              :: i, i_incr
        ttseries = .false.
        if( present(tseries) ) ttseries = tseries
        allocate(eopart(self%n))
        if( ttseries )then
            do i=1,self%n,2
                eopart(i) = 1
                i_incr = i + 1 
                if(i_incr > self%n) exit
                eopart(i_incr) = 2
            end do
        else
            rt = ran_tabu(self%n)
            call rt%balanced(2, eopart)
            call rt%kill
        endif
        eopart = eopart - 1 ! 0 is even 1 is odd
        call self%set_all('eo', real(eopart))
        deallocate(eopart)
    end subroutine partition_eo

    !>  \brief  transfers proj indices to class indices in self
    subroutine transf_proj2class( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('class', self%o(i)%get('proj'))
        enddo
    end subroutine transf_proj2class

    ! I/O

    !>  \brief  reads orientation info from file
    subroutine read( self, orifile, fromto, nst )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: orifile
        integer, optional, intent(in)    :: fromto(2)
        integer, optional, intent(out)   :: nst
        character(len=100) :: io_message
        integer :: file_stat, i, fnr, state, istart, iend
        if( .not. file_exists(orifile) )then
            call simple_stop("oris ; read; The file you are trying to read: "//trim(orifile)//' does not exist in cwd' )
        endif
        if( str_has_substr(orifile,'.bin') )then
            call simple_stop('this method does not support binary files; simple_oris :: read')
        endif
        io_message='No error'
        call fopen(fnr, FILE=orifile, STATUS='OLD', action='READ', iostat=file_stat,iomsg=io_message)
        call fileio_errmsg("oris ; read ,Error when opening file for reading: "//trim(orifile)//':'//trim(io_message), file_stat)
        if( present(nst) ) nst = 0
        if( present(fromto) )then
            istart = fromto(1)
            iend   = fromto(2)
            if(istart < 1) stop 'Invalid index; simple_oris%read'
            if(iend > self%n) stop 'Invalid index; simple_oris%read'
        else
            istart = 1
            iend   = self%n 
        endif
        do i = istart, iend
            call self%o(i)%read(fnr)
            if( present(nst) )then
                state = int(self%o(i)%get('state'))
                nst = max(1,max(state,nst))
            endif
        end do
        call fclose(fnr,errmsg="oris ; read ,Error when closing file")
    end subroutine read

    !>  \brief  reads CTF parameters and state info from file
    subroutine read_ctfparams_state_eo( self, ctfparamfile )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: ctfparamfile
        logical    :: params_are_there(11)
        integer    :: i
        type(oris) :: os_tmp
        if( .not. file_exists(ctfparamfile) )then
            call simple_stop ("oris ; read_ctfparams_state_eo; The file you are trying to read: "&
                &//trim(ctfparamfile)//' does not exist in cwd' )
        endif
        if( str_has_substr(ctfparamfile,'.bin') )then
            call simple_stop('this method does not support binary files; simple_oris :: read_ctfparams_state_eo')
        endif
        call os_tmp%new(self%n)
        call os_tmp%read(ctfparamfile)
        params_are_there(1)  = os_tmp%isthere('smpd')
        params_are_there(2)  = os_tmp%isthere('kv')
        params_are_there(3)  = os_tmp%isthere('cs')
        params_are_there(4)  = os_tmp%isthere('fraca')
        params_are_there(5)  = os_tmp%isthere('phshift')
        params_are_there(6)  = os_tmp%isthere('dfx')
        params_are_there(7)  = os_tmp%isthere('dfy')
        params_are_there(8)  = os_tmp%isthere('angast')
        params_are_there(9)  = os_tmp%isthere('bfac')
        params_are_there(10) = os_tmp%isthere('state')
        params_are_there(11) = os_tmp%isthere('eo')
        do i=1,self%n
            if( params_are_there(1) )  call self%set(i, 'smpd',    os_tmp%get(i, 'smpd')   )
            if( params_are_there(2) )  call self%set(i, 'kv',      os_tmp%get(i, 'kv')     )
            if( params_are_there(3) )  call self%set(i, 'cs',      os_tmp%get(i, 'cs')     )
            if( params_are_there(4) )  call self%set(i, 'fraca',   os_tmp%get(i, 'fraca')  )
            if( params_are_there(5) )  call self%set(i, 'phshift', os_tmp%get(i, 'phshift'))
            if( params_are_there(6) )  call self%set(i, 'dfx',     os_tmp%get(i, 'dfx')    )
            if( params_are_there(7) )  call self%set(i, 'dfy',     os_tmp%get(i, 'dfy')    )
            if( params_are_there(8) )  call self%set(i, 'angast',  os_tmp%get(i, 'angast') )
            if( params_are_there(9) )  call self%set(i, 'bfac',    os_tmp%get(i, 'bfac')   )
            if( params_are_there(10) ) call self%set(i, 'state',   os_tmp%get(i, 'state')  )
            if( params_are_there(11) ) call self%set(i, 'eo',      os_tmp%get(i, 'eo')     )
        end do
        call os_tmp%kill
    end subroutine read_ctfparams_state_eo

    !>  \brief  writes orientation info to file
    subroutine write_1( self, orifile, fromto )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: orifile
        integer, optional, intent(in)    :: fromto(2)
        character(len=100) :: io_message
        integer            :: file_stat, fnr, i, ffromto(2), cnt
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        call fopen(fnr, orifile, status='REPLACE', action='WRITE',&
        &iostat=file_stat, iomsg=io_message)
        call fileio_errmsg(' Error opening file for writing: '//trim(orifile)//' ; '//trim(io_message), file_stat)
        cnt = 0
        do i=ffromto(1),ffromto(2)
            cnt = cnt + 1
            call self%o(i)%write(fnr)
        end do
        call fclose(fnr, errmsg=' Error closing file for writing: '//trim(orifile))
    end subroutine write_1

    !>  \brief  writes orientation info to file
    subroutine write_2( self, i, orifile  )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: orifile
        integer,          intent(in)    :: i
        integer :: fnr, file_stat
        call fopen(fnr, orifile, status='UNKNOWN', action='WRITE', position='APPEND', iostat=file_stat)
        call fileio_errmsg( 'In: write_2, module: simple_oris.f90  opening '//trim(orifile), file_stat )
        call self%o(i)%write(fnr)
        call fclose(fnr, errmsg=' Error closing file for writing: '//trim(orifile))
    end subroutine write_2

    ! CALCULATORS

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
        class(oris), intent(inout) :: self
        real,        intent(in)    :: angerr, sherr
        real    :: x, y, e1, e2, e3
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
        class(oris), intent(inout) :: self
        real,        intent(in)    :: dferr
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
        class(ori),  intent(in)    :: e
        type(ori) :: o_tmp
        integer   :: i
        do i=1,self%n
             o_tmp = e.compose.self%o(i)
             call self%o(i)%set_euler(o_tmp%get_euler())
        end do
    end subroutine rot_1

    !>  \brief  is an Euler angle composer
    subroutine rot_2( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
        type(ori)                  :: o_tmp
        o_tmp = e.compose.self%o(i)
        call self%o(i)%set_euler(o_tmp%get_euler())
    end subroutine rot_2

    !>  \brief  for identifying the median value of parameter which within the cluster class
    function median_1( self, which, class ) result( med )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        real, allocatable :: vals(:)
        integer :: pop, i
        real :: med
        if( present(class) )then
            med = 0.
            pop = self%get_pop(class, 'class')
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
            allocchk('In: median_1, module: simple_oris')
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
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        real, allocatable :: vals(:), all_vals(:)
        real, allocatable :: states(:)
        states   = self%get_all('state')
        all_vals = self%get_all(which)
        vals     = pack(all_vals, mask=(states > 0.5))
        call moment(vals, ave, sdev, var, err)
        deallocate(vals, all_vals, states)
    end subroutine stats

    !>  \brief  is for calculating the minimum/maximum values of a variable
    subroutine minmax( self, which, minv, maxv )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: minv, maxv
        real, allocatable :: vals(:)
        integer           :: i, cnt, mystate
        allocate( vals(self%n), stat=alloc_stat )
        allocchk('In: minmax, module: simple_oris')
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
        class(oris),    intent(inout) :: self
        real    :: h, theta, psi
        integer :: k
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
        integer,     intent(in)    :: nsym
        real,        intent(in)    :: eullims(3,2)
        logical, allocatable :: avail(:)
        type(oris) :: tmp
        integer    :: cnt, i, n, nprojs, lim
        real       :: e1lim, e2lim, frac, frac1, frac2
        if( nsym == 1 )then
            call self%spiral_1
            return
        endif
        e1lim  = eullims(1,2)
        e2lim  = eullims(2,2)
        frac1  = 360./e1lim
        frac2  = 1. / ( (1.-cos(deg2rad(e2lim))) /2. )
        frac   = frac1 * frac2 ! area sphere / area asu
        n      = ceiling(real(self%n) * frac) ! was n = nsym * self%n
        call gen_c1
        nprojs = count(avail)
        if( nprojs < self%n )then
            ! under sampling
            n = n + self%n/2
            call gen_c1
            nprojs = count(avail)
        endif
        if( nprojs > self%n )then
            ! over sampling
            lim = max(2, floor(real(nprojs)/real(nprojs-self%n)))
            cnt  = 0
            do i = 1, n
                if(.not.avail(i))cycle
                cnt  = cnt + 1
                if(cnt == lim) then
                    avail(i) = .false.
                    cnt      = 0
                    nprojs   = nprojs-1
                    if( nprojs == self%n )exit
                endif
            enddo
        endif
        ! copy asu
        cnt = 0
        do i = 1, n
            if( avail(i) )then
                cnt = cnt + 1
                if(cnt > self%n)exit
                self%o(cnt) = tmp%o(i)
            endif
        enddo

        contains

            subroutine gen_c1
                integer :: i
                if( allocated(avail) )deallocate(avail, stat=alloc_stat)
                allocate(avail(n), source=.false., stat=alloc_stat)
                call tmp%new(n)
                call tmp%spiral_1
                do i = 1, n
                    if( tmp%o(i)%e1get() <= e1lim .and. tmp%o(i)%e2get() <= e2lim )&
                    &avail(i) = .true.
                end do
            end subroutine gen_c1
    end subroutine spiral_2

    !>  \brief  is for generating a quasi-spiral within the asymetric unit
    subroutine qspiral( self, thresh, nsym, eullims )
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

    !>  \brief  orders oris according to specscore
    function order( self ) result( inds )
        class(oris), intent(inout) :: self
        real,    allocatable :: specscores(:)
        integer, allocatable :: inds(:)
        integer :: i
        allocate( inds(self%n), stat=alloc_stat )
        allocchk('order; simple_oris')
        specscores = self%get_all('specscore')
        inds = (/(i,i=1,self%n)/)
        call hpsort(specscores, inds)
        call reverse(inds)
    end function order

    !>  \brief  orders oris according to corr
    function order_corr( self ) result( inds )
        class(oris), intent(inout) :: self
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        integer :: i
        allocate( inds(self%n), stat=alloc_stat )
        allocchk('order; simple_oris')
        corrs = self%get_all('corr')
        inds = (/(i,i=1,self%n)/)
        call hpsort(corrs, inds)
    end function order_corr

    !>  \brief  orders clusters according to population
    function order_cls( self, ncls ) result( inds )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls
        integer, allocatable :: inds(:)
        real    :: classpops(ncls)
        integer :: i
        if(ncls <= 0) stop 'invalid number of classes; simple_oris%order_cls'
        allocate(inds(ncls), stat=alloc_stat)
        allocchk('order_cls; simple_oris')
        classpops = 0.0
        ! calculate class populations
        do i=1,ncls
            classpops(i) = real(self%get_pop(i, 'class'))
        end do
        inds = (/(i,i=1,ncls)/)
        call hpsort(classpops, inds)
        call reverse(inds)
    end function order_cls

    !>  \brief  reduces the number of projection directions, useful for balancing 
    !!          operations
    subroutine reduce_projs( self, nprojs_reduced, nsym, eullims )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nprojs_reduced, nsym
        real,        intent(in)    :: eullims(3,2)
        type(oris) :: osubspace
        type(ori)  :: o_single
        integer    :: iptcl
        ! generate discrete projection direction space
        call osubspace%new( nprojs_reduced )
        call osubspace%spiral( nsym, eullims )
        ! reduction
        do iptcl=1,self%n
            o_single = self%get_ori(iptcl)
            call self%set(iptcl, 'proj', real(osubspace%find_closest_proj(o_single)))
        end do
        call osubspace%kill
    end subroutine reduce_projs

    !>  \brief  applies a one-sided balance restraint on the number of particles
    !!          in 2D classes based on corr/specscore order
    subroutine balance_1( self, popmax, skewness, use_specscore )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: popmax
        real,              intent(out)   :: skewness
        logical, optional, intent(in)    :: use_specscore
        integer, allocatable :: inds(:), inds_tmp(:)
        logical, allocatable :: included(:)
        real,    allocatable :: scores(:)
        integer :: i, j, n, pop
        logical :: uuse_specscore
        if( .not. self%isthere('w') ) stop 'ERROR, oris :: balance_1 assumes w  set'
        uuse_specscore = .false.
        if( present(use_specscore) )then
            if( use_specscore )then
                uuse_specscore = self%isthere('specscore')
            endif
        endif
        if( .not. self%isthere('class') ) stop 'class label must be set; oris :: balance_1'
        n = self%get_n('class')
        allocate(included(self%n),stat=alloc_stat)
        allocchk('In: balance, module: simple_oris')
        included = .false.
        do i=1,n ! n is number of classes
            ! get class population
            pop = self%get_pop(i, 'class', consider_w=.true.)
            ! zero case
            if( pop < 1 ) cycle
            ! get indices of particles in class subject to weight not zero
            inds = self%get_pinds(i, 'class', consider_w=.true.)
            if( pop <= popmax )then ! all inclded
                do j=1,pop
                    included(inds(j)) = .true.
                end do
            else ! popmax threshold applied
                allocate(scores(pop), inds_tmp(pop))
                do j=1,pop
                    if( uuse_specscore )then
                        scores(j) = self%o(inds(j))%get('specscore')
                    else
                        scores(j) = self%o(inds(j))%get('corr')
                    endif
                    inds_tmp(j) = inds(j)
                end do
                call hpsort(scores, inds_tmp)
                do j=pop,pop - popmax + 1,-1
                    included(inds_tmp(j)) = .true. 
                end do
                deallocate(scores, inds_tmp)
            endif
        end do
        ! communicate selection to instance
        do i=1,self%n
            if( included(i) )then
                call self%o(i)%set('state_balance', 1.0)
            else
                call self%o(i)%set('state_balance', 0.0)
            endif
        end do
        ! communicate skewness as fraction of excluded
        skewness = real(self%n - count(included))/real(self%n)
    end subroutine balance_1

    !>  \brief  applies a one-sided balance restraint on the number of particles
    !!          in projection groups based on corr/specscore order
    subroutine balance_2( self, popmax, nspace_bal, nsym, eullims, skewness, use_specscore )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: popmax, nspace_bal, nsym
        real,              intent(in)    :: eullims(3,2)
        real,              intent(out)   :: skewness
        logical, optional, intent(in)    :: use_specscore
        type(oris)           :: osubspace
        type(ori)            :: o_single
        integer, allocatable :: inds(:), inds_tmp(:), clustering(:), clustszs(:)
        real,    allocatable :: scores(:)
        logical, allocatable :: ptcl_mask(:), included(:)
        integer              :: i, j, n, pop, iptcl, icls, cnt
        logical              :: uuse_specscore
        if( .not. self%isthere('w') ) stop 'ERROR, oris :: balance assumes w  set'
        uuse_specscore = .false.
        if( present(use_specscore) )then
            if( use_specscore )then
                uuse_specscore = self%isthere('specscore')
            endif
        endif
        if( .not. self%isthere('proj') ) stop 'proj label must be set; oris :: balance'
        n = self%get_n('proj')
        ptcl_mask = self%included(consider_w=.true.)
        allocate(clustering(self%n), clustszs(nspace_bal), included(self%n), stat=alloc_stat)
        allocchk('In oris :: balance')
        clustering = 0
        clustszs   = 0
        included   = .false.
        ! generate discrete projection direction space
        call osubspace%new( nspace_bal )
        call osubspace%spiral( nsym, eullims )
        ! cluster the projection directions
        do iptcl=1,self%n
            if( ptcl_mask(iptcl) )then
                o_single = self%get_ori(iptcl)
                clustering(iptcl) = osubspace%find_closest_proj(o_single)
            else
                clustering(iptcl) = 0
            endif
        end do
        ! loop over clusters
        do icls=1,nspace_bal
            pop = count(clustering == icls)
            ! zero case
            if( pop < 1 ) cycle
            ! get indices of particles in cluster
            allocate(inds(pop))
            cnt = 0 
            do iptcl=1,self%n
                if( clustering(iptcl) == icls )then
                    cnt = cnt + 1
                    inds(cnt) = iptcl
                endif
            end do
            if( pop <= popmax )then ! all inclded
                do j=1,pop
                    included(inds(j)) = .true.
                end do
            else ! popmax threshold applied
                allocate(scores(pop), inds_tmp(pop))
                do j=1,pop
                    if( uuse_specscore )then
                        scores(j) = self%o(inds(j))%get('specscore')
                    else
                        scores(j) = self%o(inds(j))%get('corr')
                    endif
                    inds_tmp(j) = inds(j)
                end do
                call hpsort(scores, inds_tmp)
                do j=pop,pop - popmax + 1,-1
                    included(inds_tmp(j)) = .true. 
                end do
                deallocate(scores, inds_tmp)
            endif
            deallocate(inds)
        end do
        ! communicate selection to instance
        do i=1,self%n
            if( included(i) )then
                call self%o(i)%set('state_balance', 1.0)
            else
                call self%o(i)%set('state_balance', 0.0)
            endif
        end do
        ! communicate skewness as fraction of excluded
        skewness = real(self%n - count(included))/real(self%n)
    end subroutine balance_2

    !>  \brief  calculates spectral particle weights
    subroutine calc_spectral_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        integer           :: i, nstates, istate, cnt, mystate
        real, allocatable :: weights(:), specscores(:)
        call self%calc_hard_weights( frac )
        if( self%isthere('specscore') )then
            nstates = self%get_n('state')
            if( nstates > 1 )then
                do istate=1,nstates
                    specscores = self%get_arr('specscore', state=istate)
                    weights    = self%get_arr('w',         state=istate)
                    where( weights < 0.5 )
                        specscores = 0.
                    end where
                    weights = corrs2weights(specscores)
                    cnt = 0
                    do i=1,self%n
                        mystate = nint(self%o(i)%get('state'))
                        if( mystate == istate )then
                            cnt = cnt + 1
                            call self%o(i)%set('w', weights(cnt))
                        else if( mystate == 0 )then
                            call self%o(i)%set('w', 0.0)
                        endif
                    enddo
                    deallocate(specscores,weights)
                enddo
            else
                specscores = self%get_all('specscore')
                weights    = self%get_all('w')
                where( weights < 0.5 )
                    specscores = 0.
                end where
                weights = corrs2weights(specscores)
                do i=1,self%n
                    call self%o(i)%set('w', weights(i))
                end do
                deallocate(specscores,weights)
            endif
        else
            !stop 'specscore not part of oris; simple_oris :: calc_spectral_weights'
        endif
    end subroutine calc_spectral_weights

    subroutine reject_above( self, which, thres )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: thres
        integer :: i
        real    :: val
        if( self%isthere(which) )then
            do i=1,self%n
                val = self%o(i)%get(which)
                if( val > thres ) call self%o(i)%set('state_balance', 0.)
            end do
        else
            print *, 'which: ', trim(which)
            stop 'variable not part of oris; simple_oris :: reject_above'
        endif
    end subroutine reject_above

    !>  \brief  calculates hard weights based on ptcl ranking
    subroutine calc_hard_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        integer, allocatable :: order(:)
        integer :: i, lim, ind
        if( frac < 0.99 )then
            lim   = nint(frac*real(self%n))
            order = self%order() ! specscore ranking
            do i=1,self%n
                ind = order(i)
                if( i <= lim )then
                    call self%o(ind)%set('w', 1.)
                else
                    call self%o(ind)%set('w', 0.)
                endif
            end do
        else
            call self%set_all2single('w', 1.)
        endif
    end subroutine calc_hard_weights

    !>  \brief  to find the closest matching projection direction
    function find_closest_proj( self, o_in, state ) result( closest )
        class(oris),       intent(inout) :: self
        class(ori),        intent(in) :: o_in
        integer, optional, intent(in) :: state
        real    :: dists(self%n), large
        integer :: loc(1), closest, i
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
        real    :: dists(self%n)
        integer :: inds(self%n), i
        do i=1,self%n
            dists(i) = self%o(i).euldist.o_in
            inds(i)  = i
        end do
        call hpsort(dists, inds)
        do i=1,size(pdirinds)
            pdirinds(i) = inds(i)
        end do
    end subroutine find_closest_projs

    !>  \brief  to find the closest matching orientation with respect to state
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
        real    :: dists(self%n)
        integer :: inds(self%n), i
        do i=1,self%n
            dists(i) = self%o(i).geod.o_in
            inds(i)  = i
        end do
        call hpsort(dists, inds)
        do i=1,size(oriinds)
            oriinds(i) = inds(i)
        end do
    end subroutine find_closest_oris

    !>  \brief  to calculate all euler distances to one orientation
    subroutine calc_euldists( self, o_in, dists )
        class(oris),       intent(in)  :: self
        class(ori),        intent(in)  :: o_in
        real,              intent(out) :: dists(self%n)
        integer :: i
        do i=1,self%n
            dists(i) = self%o(i).euldist.o_in
        end do
    end subroutine calc_euldists

    !>  \brief  to identify a subspace of projection directions
    function create_proj_subspace_1( self, nsub ) result( subspace_projs )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsub
        type(ori)            :: o
        type(oris)           :: suboris
        integer, allocatable :: subspace_projs(:)
        integer :: isub
        call suboris%new(nsub)
        call suboris%spiral
        allocate( subspace_projs(nsub) ,stat=alloc_stat)
        allocchk('In: create_proj_subspace_1, module: simple_oris')
        do isub=1,nsub
            o = suboris%get_ori(isub)
            subspace_projs(isub) = self%find_closest_proj(o)
        end do
    end function create_proj_subspace_1

    !>  \brief  to identify a subspace of projection directions
    function create_proj_subspace_2( self, nsub, nsym, eullims ) result( subspace_projs )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsub
        integer,     intent(in)    :: nsym
        real,        intent(in)    :: eullims(3,2)
        type(ori)            :: o
        type(oris)           :: suboris
        integer, allocatable :: subspace_projs(:)
        integer :: isub
        call suboris%new(nsub)
        call suboris%spiral(nsym, eullims)
        allocate( subspace_projs(nsub) , stat=alloc_stat)
        allocchk('In: create_proj_subspace_2, module: simple_oris')
        do isub=1,nsub
            o = suboris%get_ori(isub)
            subspace_projs(isub) = self%find_closest_proj(o)
        end do
    end function create_proj_subspace_2

    !>  \brief  method for discretization of the projection directions
    subroutine discretize( self, n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
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

    !>  \brief  to identify the indices of the k nearest projection neighbors (inclusive)
    subroutine nearest_proj_neighbors( self, k, nnmat )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: k
        integer, allocatable, intent(inout) :: nnmat(:,:)
        real      :: dists(self%n)
        integer   :: inds(self%n), i, j
        type(ori) :: o
        if( k >= self%n ) stop 'need to identify fewer nearest_neighbors; simple_oris'
        if( allocated(nnmat) ) deallocate(nnmat)
        allocate( nnmat(self%n,k), stat=alloc_stat )
        allocchk("In: nearest_neighbors; simple_oris")
        do i=1,self%n
            o = self%get_ori(i)
            do j=1,self%n
                inds(j) = j
                if( i == j )then
                    dists(j) = 0.
                else
                    dists(j) = self%o(j).euldist.o
                endif
            end do
            call hpsort(dists, inds)
            do j=1,k
                nnmat(i,j) = inds(j)
            end do
        end do
    end subroutine nearest_proj_neighbors

    !>  \brief  to find angular resolution of an even orientation distribution (in degrees)
    function find_angres( self ) result( res )
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
            call hpsort(dists)
            res = res+sum(dists(:3))/3. ! average of three nearest neighbors
        end do
        res = rad2deg(res/real(self%n))
    end function find_angres

    !>  \brief  to find angular resolution of an even orientation distribution (in degrees)
    function find_angres_geod( self ) result( res )
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
            call hpsort(dists)
            res = res+sum(dists(:3))/3. ! average of three nearest neighbors
        end do
        res = rad2deg(res/real(self%n))
    end function find_angres_geod

    !>  \brief  to find the correlation bound in extremal search
    function extremal_bound( self, thresh ) result( corr_bound )
        class(oris),       intent(inout) :: self
        real,              intent(in)    :: thresh
        real,    allocatable       :: corrs(:), corrs_incl(:)
        logical, allocatable       :: incl(:)
        integer :: n_incl, thresh_ind
        real    :: corr_bound
        ! grab relevant correlations
        corrs      = self%get_all('corr')
        incl       = self%included()
        corrs_incl = pack(corrs, mask=incl)
        ! sort correlations & determine threshold
        n_incl     = size(corrs_incl)
        call hpsort(corrs_incl)
        thresh_ind  = nint(real(n_incl) * thresh)
        corr_bound  = corrs_incl(thresh_ind)
        deallocate(corrs, incl, corrs_incl)
    end function extremal_bound

    !>  \brief  to find the neighborhood size for weighted orientation assignment
    function find_npeaks( self, res, moldiam ) result( npeaks )
        class(oris), intent(in) :: self
        real,        intent(in) :: res, moldiam
        real                    :: dists(self%n), tres, npeaksum
        integer                 :: i, j, npeaks
        tres = atan(res/(moldiam/2.))
        npeaksum = 0.
        !$omp parallel do schedule(static) default(shared) proc_bind(close)&
        !$omp private(j,i,dists) reduction(+:npeaksum)
        do j=1,self%n
            do i=1,self%n
                dists(i) = self%o(i).euldist.self%o(j)
            end do
            call hpsort(dists)
            do i=1,self%n
                if( dists(i) <= tres )then
                    npeaksum = npeaksum + 1.
                else
                    exit
                endif
            end do
        end do
        !$omp end parallel do
        npeaks = nint(npeaksum/real(self%n)) ! nr of peaks is average of peaksum
    end function find_npeaks

    !>  \brief  find angular threshold from number of peaks
    function find_athres_from_npeaks( self, npeaks ) result( athres )
        class(oris), intent(in) :: self
        integer,     intent(in) :: npeaks
        real :: athres, dist, maxdist
        integer :: pdirinds(npeaks), i, j
        athres = 0.
        do i=1,self%n
            call self%find_closest_projs(self%o(i), pdirinds)
            maxdist = 0.
            do j=1,npeaks
                if( i /= pdirinds(j) )then
                    dist = self%o(i).euldist.self%o(pdirinds(j))
                    if( dist > maxdist )then
                        maxdist = dist
                    endif
                endif
            end do
            athres = athres + maxdist
        end do
        athres = rad2deg(athres / real(self%n))
    end function find_athres_from_npeaks
    
    !>  \brief  find number of peaks from angular threshold
    function find_npeaks_from_athres( self, athres ) result( npeaks )
        class(oris), intent(in) :: self
        real,        intent(in) :: athres
        integer :: i, j, npeaks
        real :: dists(self%n), rnpeaks
        rnpeaks = 0.
        do i=1,self%n
            do j=1,self%n
                dists(j) = rad2deg( self%o(i).euldist.self%o(j) )
            end do
            rnpeaks = rnpeaks + real(count(dists <= athres))
        end do
        npeaks = nint(rnpeaks/real(self%n))
    end function find_npeaks_from_athres

    !>  \brief  for mapping a 3D shift of volume to 2D shifts of the projections
    subroutine map3dshift22d_1( self, sh3d, state )
        class(oris),       intent(inout) :: self
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
        integer :: i
        do i=1,self%n
            call self%map3dshift22d_2(i, sh3d, state)
        end do
    end subroutine map3dshift22d_1

    !>  \brief  for mapping a 3D shift of volume to 2D shifts of the projections
    subroutine map3dshift22d_2( self, i, sh3d, state )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
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

    !>  \brief  for correlating oris objs, just for testing purposes
    function corr_oris( self1, self2 ) result( corr )
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
        use simple_jiffys, only: progress
        class(oris), intent(inout) :: self
        logical, allocatable       :: o_is_set(:)
        integer                    :: i
        type(ori)                  :: o1, o2
        ! create the first diverse pair
        call o1%oripair_diverse(o2)
        allocate(o_is_set(self%n), stat=alloc_stat)
        allocchk('In: gen_diverse, module: simple_oris')
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
        integer :: i
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
        class(oris), target, intent(in) :: self
        character(len=*),    intent(in) :: mode
        logical, optional,   intent(in) :: part_of_set(self%n)
        real,    optional,   intent(in) :: weights(self%n)
        real, parameter                 :: TOL=1e-6
        type(ori)                       :: oout
        type(opt_spec)                  :: ospec
        type(simplex_opt)               :: opt
        real                            :: dist, lims(3,2)
        class(*), pointer               :: fun_self => null()
        ! manage class vars
        if( allocated(class_part_of_set) ) deallocate(class_part_of_set)
        if( allocated(class_weights)     ) deallocate(class_weights)
        if( present(part_of_set) )then
            allocate(class_part_of_set(self%n), source=part_of_set, stat=alloc_stat)
        else
            allocate(class_part_of_set(self%n), stat=alloc_stat)
            class_part_of_set = .true.
        endif
        allocchk('In: ori_generator, module: simple_oris part_of_set')
        if( present(weights) )then
            allocate(class_weights(self%n), source=weights, stat=alloc_stat)
        else
            allocate(class_weights(self%n), stat=alloc_stat)
            class_weights = 1.0
        endif
        allocchk('In: ori_generator, module: simple_oris weights')
        ! set globals
        ops => self
        call o_glob%new
        call o_glob%rnd_euler
        ! init
        call oout%new
        lims(1,1) = 0.
        lims(1,2) = 359.99
        lims(2,1) = 0.
        lims(2,2) = 180.
        lims(3,1) = 0.
        lims(3,2) = 359.99
        call ospec%specify('simplex', 3, ftol=TOL, limits=lims)
        ospec%x(1) = o_glob%e1get()
        ospec%x(2) = o_glob%e2get()
        ospec%x(3) = o_glob%e3get()
        call opt%new(ospec)
        select case(mode)
            case('median')
                call ospec%set_costfun(costfun_median)
            case('diverse')
                call ospec%set_costfun(costfun_diverse)
            case DEFAULT
                stop 'unsupported mode inputted; simple_ori :: ori_generator'
        end select
        call opt%minimize(ospec, fun_self, dist)
        call oout%set_euler(ospec%x)
        call ospec%kill
        call opt%kill
    end function ori_generator

    !>  \brief  generates a similarity matrix for the oris in self
    function gen_smat( self ) result( smat )
        class(oris), intent(in) :: self
        real, allocatable :: smat(:,:)
        integer :: i, j
        allocate( smat(self%n,self%n), stat=alloc_stat )
        allocchk("In: simple_oris :: gen_smat")
        smat = 0.
        !$omp parallel do schedule(guided) default(shared) private(i,j) proc_bind(close)
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
        if( allocated(class_part_of_set) ) deallocate(class_part_of_set)
        if( allocated(class_weights)     ) deallocate(class_weights)
        if( present(part_of_set) )then
            allocate(class_part_of_set(self%n), source=part_of_set)
        else
            allocate(class_part_of_set(self%n))
            class_part_of_set = .true.
        endif
        if( present(weights) )then
            allocate(class_weights(self%n), source=weights)
        else
            allocate(class_weights(self%n))
            class_weights = 1.0
        endif
        dists = 0.
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i=1,self%n
            dists(i) = self%o(i).geod.o
        end do
        !$omp end parallel do
        geodesic_dist = sum(dists*class_weights,mask=class_part_of_set)/sum(class_weights,mask=class_part_of_set)
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
        class(oris), intent(in)  :: self
        real,        intent(out) :: mind, maxd, avgd, sdevd, sumd
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
        class(oris), intent(in)  :: self1, self2
        real,        intent(out) :: mind, maxd, avgd, sdevd, sumd
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
        class(oris), intent(inout) :: self
        real,        intent(out)   :: avgd, sdevd, maxd, mind
        integer, allocatable       :: clsarr(:)
        real,    allocatable       :: dists(:)
        integer                    :: ncls, icls, iptcl, jptcl, sz, cnt, nobs, n
        real                       :: avgd_here, sdevd_here, vard_here, mind_here, maxd_here
        type(ori)                  :: iori, jori
        logical                    :: err
        ncls  = self%get_n('class')
        nobs  = 0
        avgd  = 0.
        sdevd = 0.
        maxd  = 0.
        mind  = 0.
        do icls=1,ncls
            clsarr = self%get_pinds(icls, 'class')
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

    ! PRIVATE ROUTINES

    function costfun_median( fun_self, vec, D ) result( dist )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: dist
        call o_glob%set_euler(vec)
        ! we are minimizing the distance
        dist = ops%geodesic_dist(o_glob,class_part_of_set,class_weights)
    end function costfun_median

    function costfun_diverse( fun_self, vec, D ) result( dist )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: dist
        call o_glob%set_euler(vec)
        ! we are maximizing the distance
        dist = -ops%geodesic_dist(o_glob,class_part_of_set,class_weights)
    end function costfun_diverse

    ! UNIT TEST

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
                call os%print_(i)
            end do
            call os2%rnd_oris(5.)
            write(*,*) '********'
            do i=1,100
                call os2%print_(i)
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
        call os%write('test_oris_rndoris'//METADATEXT)
        call os2%read('test_oris_rndoris'//METADATEXT)
        call os2%write('test_oris_rndoris_copy'//METADATEXT)
        corr = corr_oris(os,os2)
        if( corr > 0.99 ) passed = .true.
        if( .not. passed ) stop 'read/write failed'
        passed = .false.
        call os%rnd_states(5)
        call os%write('test_oris_rndoris_rndstates'//METADATEXT)
        if( corr_oris(os,os2) > 0.99 ) passed = .true.
        if( .not. passed ) stop 'statedoc read/write failed!'
        write(*,'(a)') '**info(simple_oris_unit_test, part3): testing calculators'
        passed = .false.
        call os%rnd_lps()
        call os%write('test_oris_rndoris_rndstates_rndlps'//METADATEXT)
        call os%spiral
        call os%write('test_oris_rndoris_rndstates_rndlps_spiral'//METADATEXT)
        call os%rnd_corrs()
        order = os%order()
        if( doprint )then
            do i=1,100
                call os%print_(order(i))
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

    ! DESTRUCTORS

    !>  \brief  is a destructor
    subroutine kill_chash( self )
        class(oris), intent(inout) :: self
        integer :: i
        if( allocated(self%o) )then
            do i=1,self%n
                call self%o(i)%kill_chash
            end do
        endif
    end subroutine kill_chash

    !>  \brief  is a destructor
    subroutine kill( self )
        class(oris), intent(inout) :: self
        integer :: i
        if( allocated(self%o) )then
            do i=1,self%n
                call self%o(i)%kill
            end do
            deallocate( self%o , stat=alloc_stat)
            allocchk('In: kill, module: simple_oris o')
            self%n = 0
        endif
    end subroutine kill

end module simple_oris
