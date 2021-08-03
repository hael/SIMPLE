! an agglomeration of orientations
module simple_oris
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_ori, only: ori
use simple_defs_ori
implicit none

public :: oris, test_oris
private
#include "simple_local_flags.inc"

!>  \brief struct type aggregates ori objects
type :: oris
    private
    type(ori), allocatable :: o(:)
    integer :: n = 0
  contains
    ! CONSTRUCTORS
    procedure          :: new
    procedure          :: reallocate
    ! GETTERS
    procedure          :: exists
    procedure          :: e1get
    procedure          :: e2get
    procedure          :: e3get
    procedure          :: get_euler
    procedure          :: get_noris
    procedure          :: get_ori
    procedure          :: get
    procedure          :: get_static
    procedure, private :: getter_1
    procedure, private :: getter_2
    generic            :: getter => getter_1, getter_2
    procedure          :: get_all
    procedure          :: get_all_rmats
    procedure          :: get_mat
    procedure          :: get_normal
    procedure          :: get_2Dshift
    procedure          :: get_2Dshift_incr
    procedure          :: get_state
    procedure          :: get_class
    procedure          :: get_tseries_neighs
    procedure, private :: isthere_1
    procedure, private :: isthere_2
    generic            :: isthere => isthere_1, isthere_2
    procedure          :: ischar
    procedure          :: is_particle
    procedure          :: max_ori_strlen_trim
    procedure          :: get_n
    procedure          :: get_pop
    procedure          :: get_pops
    procedure          :: get_pinds
    procedure          :: gen_mask
    procedure          :: mask_from_state
    procedure, private :: get_all_normals
    procedure          :: states_exist
    procedure          :: get_arr
    procedure, private :: calc_sum
    procedure          :: get_sum
    procedure          :: get_avg
    procedure          :: included
    procedure          :: get_nevenodd
    procedure, private :: get_neven
    procedure, private :: get_nodd
    procedure          :: print_
    procedure          :: print_matrices
    procedure          :: sample4update_and_incrcnt_nofrac
    procedure          :: sample4update_and_incrcnt
    procedure          :: sample4update_and_incrcnt2D
    procedure          :: sample_rnd_subset
    procedure          :: incr_updatecnt
    procedure          :: has_been_searched
    procedure          :: any_state_zero
    procedure          :: ori2str
    procedure          :: ori2prec
    procedure          :: prec2ori
    procedure          :: get_ctfvars
    ! SETTERS
    procedure          :: copy
    procedure          :: reject
    generic            :: delete_entry => delete_entry_1, delete_entry_2
    procedure          :: delete_entry_1
    procedure          :: delete_entry_2
    procedure          :: delete_2Dclustering
    procedure          :: delete_3Dalignment
    procedure, private :: transfer_2Dparams_1, transfer_2Dparams_2
    generic            :: transfer_2Dparams => transfer_2Dparams_1, transfer_2Dparams_2
    procedure          :: transfer_3Dparams
    procedure          :: set_euler
    procedure          :: set_shift
    procedure          :: set_shift_incr
    procedure          :: e1set
    procedure          :: e2set
    procedure          :: e3set
    procedure, private :: set_1
    procedure, private :: set_2
    generic            :: set => set_1, set_2
    procedure          :: set_ori
    procedure          :: transfer_ori
    procedure, private :: set_all_1
    procedure, private :: set_all_2
    generic            :: set_all => set_all_1, set_all_2
    procedure, private :: set_all2single_1
    procedure, private :: set_all2single_2
    generic            :: set_all2single => set_all2single_1, set_all2single_2
    procedure          :: set_projs
    procedure          :: remap_projs
    procedure          :: proj2class
    procedure          :: e3swapsgn
    procedure          :: swape1e3
    procedure          :: zero
    procedure          :: zero_projs
    procedure          :: zero_shifts
    procedure          :: mul_shifts
    procedure          :: rnd_oris
    procedure          :: rnd_ori
    procedure          :: rnd_oris_discrete_from
    procedure          :: rnd_inpls
    procedure          :: rnd_ctf
    procedure          :: rnd_states
    procedure          :: rnd_lps
    procedure          :: rnd_corrs
    procedure          :: rnd_proj_space
    procedure          :: revshsgn
    procedure          :: revorisgn
    procedure          :: ini_tseries
    procedure          :: symmetrize
    procedure          :: merge
    procedure          :: clean_updatecnt
    procedure          :: partition_eo
    procedure          :: str2ori
    procedure          :: str2ori_ctfparams_state_eo
    procedure          :: set_ctfvars
    procedure          :: set_boxfile
    procedure          :: reset
    ! I/O
    procedure          :: read
    procedure          :: read_ctfparams_state_eo
    procedure, private :: write_1
    procedure, private :: write_2
    generic            :: write => write_1, write_2
    procedure          :: write2bild
    ! CALCULATORS
    procedure          :: compress
    procedure          :: split_state
    procedure          :: split_class
    procedure          :: expand_classes
    procedure          :: fill_empty_classes
    procedure          :: remap_cls
    procedure          :: merge_classes
    procedure          :: round_shifts
    procedure          :: introd_alig_err
    procedure          :: introd_ctf_err
    procedure, private :: rot_1
    procedure, private :: rot_2
    generic            :: rot => rot_1, rot_2
    procedure, private :: rot_transp_1
    procedure, private :: rot_transp_2
    generic            :: rot_transp => rot_transp_1, rot_transp_2
    procedure, private :: median_1
    generic            :: median => median_1
    procedure, private :: stats_1
    procedure, private :: stats_2
    generic            :: stats => stats_1, stats_2
    procedure          :: minmax
    procedure          :: spiral_1
    procedure          :: spiral_2
    generic            :: spiral => spiral_1, spiral_2
    procedure          :: order
    procedure          :: order_corr
    procedure          :: order_cls
    procedure          :: calc_hard_weights
    procedure          :: calc_soft_weights
    procedure          :: calc_hard_weights2D
    procedure          :: calc_soft_weights2D
    procedure          :: find_best_classes
    procedure          :: find_closest_proj
    procedure          :: discretize
    procedure, private :: nearest_proj_neighbors_1
    procedure, private :: nearest_proj_neighbors_2
    procedure, private :: nearest_proj_neighbors_3
    procedure, private :: nearest_proj_neighbors_4
    generic            :: nearest_proj_neighbors => nearest_proj_neighbors_1, nearest_proj_neighbors_2, nearest_proj_neighbors_3, nearest_proj_neighbors_4
    procedure          :: min_euldist
    procedure          :: find_angres
    procedure          :: extremal_bound
    procedure          :: set_extremal_vars
    procedure, private :: map3dshift22d_1
    procedure, private :: map3dshift22d_2
    generic            :: map3dshift22d => map3dshift22d_1, map3dshift22d_2
    procedure          :: mirror2d
    procedure          :: mirror3d
    procedure          :: add_shift2class
    procedure          :: corr_oris
    procedure, private :: diststat_1
    procedure, private :: diststat_2
    generic            :: diststat => diststat_1, diststat_2
    procedure          :: overlap
    ! DESTRUCTORS
    procedure          :: kill_chash
    procedure          :: kill
end type oris

interface oris
    module procedure constructor
end interface oris

contains

    ! CONSTRUCTORS

    function constructor( n, is_ptcl ) result( self )
        integer, intent(in) :: n
        logical, intent(in) :: is_ptcl
        type(oris) :: self
        call self%new(n, is_ptcl)
    end function constructor

    subroutine new( self, n, is_ptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        logical,     intent(in)    :: is_ptcl
        integer :: i
        call self%kill
        self%n = n
        allocate( self%o(self%n), stat=alloc_stat )
        if( alloc_stat.ne.0 ) call allocchk('new; simple_oris', alloc_stat)
        do i=1,n
            call self%o(i)%new(is_ptcl)
        end do
    end subroutine new

    subroutine reallocate( self, new_n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: new_n
        type(oris) :: tmp
        integer    :: old_n, i
        logical    :: is_ptcl
        if( self%n == 0 )     THROW_HARD('cannot reallocate non-existing oris; reallocate')
        if( new_n <= self%n ) THROW_HARD('reallocation to smaller size not supported; reallocate')
        is_ptcl = self%o(1)%is_particle()
        ! make a copies
        old_n = self%n
        tmp   = self
        ! reallocate
        call self%new(new_n, is_ptcl)
        ! stash back the old data
        do i=1,old_n
            self%o(i) = tmp%o(i)
        end do
    end subroutine reallocate

    ! GETTERS

    pure logical function exists( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        exists = self%o(i)%exists()
    end function exists

    pure function e1get( self, i ) result( e1 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e1
        e1 = self%o(i)%e1get()
    end function e1get

    pure function e2get( self, i ) result( e2 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e2
        e2 = self%o(i)%e2get()
    end function e2get

    pure function e3get( self, i ) result( e3 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e3
        e3 = self%o(i)%e3get()
    end function e3get

    pure function get_euler( self, i ) result( euls )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: euls(3)
        euls = self%o(i)%get_euler()
    end function get_euler

    pure function get_noris( self, consider_state ) result( n )
        class(oris),       intent(in) :: self
        logical, optional, intent(in) :: consider_state
        integer :: i, n
        logical :: consider_state_here
        consider_state_here = .false.
        if(present(consider_state)) consider_state_here = consider_state
        if( consider_state_here )then
            n = 0
            do i=1,self%n
                if( self%o(i)%isthere('state') )then
                    if( self%o(i)%get('state') > 0.5 ) n = n+1
                else
                    n = n+1 ! included by default
                endif
            enddo
        else
            n = self%n
        endif
    end function get_noris

    subroutine get_ori( self, i, o )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(inout) :: o
        if( self%n == 0 ) THROW_HARD('oris object does not exist; get_ori')
        if( i > self%n .or. i < 1 )then
            write(logfhandle,*) 'trying to get ori: ', i, ' among: ', self%n, ' oris'
            THROW_HARD('i out of range; get_ori')
        endif
        o = self%o(i)
    end subroutine get_ori

    pure function get( self, i, key ) result( val )
        class(oris),      intent(in) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real :: val
        val = self%o(i)%get(key)
    end function get

    subroutine getter_1( self, i, key, val )
        class(oris),                   intent(inout) :: self
        integer,                       intent(in)    :: i
        character(len=*),              intent(in)    :: key
        character(len=:), allocatable, intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_1

    subroutine getter_2( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_2

    !>  \brief  is a getter with fixed length return string
    function get_static( self, i, key )result( val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        character(len=STDLEN)           :: val
        val = trim(self%o(i)%get_static(key))
    end function get_static

    !>  \brief  is for getting an array of 'key' values
    function get_all( self, key, fromto ) result( arr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: fromto(2)
        real, allocatable :: arr(:)
        integer :: i, ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        if(allocated(arr))deallocate(arr)
        allocate( arr(ffromto(1):ffromto(2)), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('get_all; simple_oris',alloc_stat)
        do i=ffromto(1),ffromto(2)
            arr(i) = self%o(i)%get(key)
        enddo
    end function get_all

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

    function get_2Dshift( self, i )  result(shvec)
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real :: shvec(2)
        shvec = self%o(i)%get_2Dshift()
    end function get_2Dshift

    function get_2Dshift_incr( self, i )  result(shvec)
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real :: shvec(2)
        shvec = self%o(i)%get_2Dshift_incr()
    end function get_2Dshift_incr

    pure integer function get_state( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_state = self%o(i)%get_state()
    end function get_state

    pure integer function get_class( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_class = self%o(i)%get_class()
    end function get_class

    ! assumes project has been pruned to remove state=0 particles
    subroutine get_tseries_neighs( self, nsz, ptcls2neigh  )
       class(oris),          intent(in)    :: self
       integer,              intent(in)    :: nsz
       integer, allocatable, intent(inout) :: ptcls2neigh(:,:)
       integer :: i, j, cls1, cls2
       if( allocated(ptcls2neigh) ) deallocate(ptcls2neigh)
       allocate(ptcls2neigh(self%n,2), source=0)
       do i=1,self%n-nsz
           cls1 = nint(self%o(i)%get('class'))
           cls2 = nint(self%o(i+1)%get('class'))
           if( cls2 == cls1 )then
               cycle
           else
               do j=max(1,i-nsz+1),min(self%n,i+nsz)
                   ptcls2neigh(j,1) = cls1
                   ptcls2neigh(j,2) = cls2
               end do
           endif
       end do
   end subroutine get_tseries_neighs

    !>  \brief  is for checking if parameter is present
    pure function isthere_1( self, key ) result( is )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: key
        logical :: is
        integer :: i
        is = .false.
        do i=1,self%n
            is = self%o(i)%isthere(key)
            if( is ) exit
        end do
    end function isthere_1

    !>  \brief  is for checking if parameter is present
    pure function isthere_2( self, i, key ) result( is )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: i
        character(len=*),  intent(in) :: key
        logical :: is
        is = self%o(i)%isthere(key)
    end function isthere_2

    function ischar( self, i, key ) result ( is )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: i
        character(len=*),  intent(in) :: key
        logical :: is
        is = self%o(i)%ischar(key)
    end function ischar

    pure function is_particle( self ) result( t )
        class(oris), intent(in) :: self
        logical :: t
        t = self%o(1)%is_particle()
    end function is_particle

    !>  \brief  is for getting the maximum string length of a trimed string ori representation
    integer function max_ori_strlen_trim( self )
        class(oris), intent(in) :: self
        integer :: i
        max_ori_strlen_trim = 0
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)&
        !$omp reduction(max:max_ori_strlen_trim)
        do i=1,self%n
            max_ori_strlen_trim = max(max_ori_strlen_trim,self%o(i)%ori_strlen_trim())
        end do
        !$omp end parallel do
    end function max_ori_strlen_trim

    !>  \brief  is for getting the max val of integer label
    function get_n( self, label ) result( n )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: label
        integer :: i, n, ival
        n = 1
        do i=1,self%n
            ival = nint(self%o(i)%get(trim(label)))
            if( ival > n ) n = ival
        end do
    end function get_n

    !>  \brief  is for checking label population
    function get_pop( self, ind, label, consider_w, eo ) result( pop )
        class(oris),       intent(in)    :: self
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
            if( .not. self%isthere('w') ) THROW_HARD('get_pop with optional consider_w assumes w set')
        endif
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        pop = 0
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( consider_eo )then
                myeo = nint(self%o(i)%get('eo'))
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

    !>  \brief  is for getting all rotation matrices
    function get_all_rmats( self ) result( mat )
        class(oris), intent(in) :: self
        real, allocatable       :: mat(:,:,:)
        integer :: i,n
        n = self%n
        if(allocated(mat))deallocate(mat)
        allocate(mat(n,3,3),source=0.,stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: geteulmats, module: simple_oris')
        do i=1,self%n
            mat(i,:,:) = self%o(i)%get_mat()
        end do
    end function get_all_rmats

    subroutine get_pops( self, pops, label, consider_w, maxn, eo )
        class(oris),          intent(in)    :: self
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
            if( .not. self%isthere('w') ) THROW_HARD('get_pops with optional consider_w assumes w set')
        endif
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        n = self%get_n(label)
        if( present(maxn) )then
            n = max(n, maxn)
        endif
        if(allocated(pops))deallocate(pops)
        allocate(pops(n),source=0,stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: get_pops, module: simple_oris')
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( consider_eo )then
                myeo = nint(self%o(i)%get('eo'))
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

    function get_all_normals(self) result(normals)
        class(oris), intent(inout) :: self
        real, allocatable :: normals(:,:)
        integer :: i
        allocate(normals(self%n,3), stat=alloc_stat)
        do i=1,self%n
            normals(i,:) = self%o(i)%get_normal()
        end do
    end function get_all_normals

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
        logical    :: is_ptcl
        if( size(mask) /= self%n )then
            write(logfhandle,*) 'self%n:     ', self%n
            write(logfhandle,*) 'size(mask): ', size(mask)
            THROW_HARD('nonconforming mask size; compress')
        endif
        is_ptcl = self%o(1)%is_particle()
        call os_tmp%new(count(mask), is_ptcl)
        cnt = 0
        do i=1,self%n
            if( mask(i) )then
                cnt = cnt + 1
                os_tmp%o(cnt) = self%o(i)
            endif
        end do
        call self%copy(os_tmp, is_ptcl)
        call os_tmp%kill
    end subroutine compress

    !>  \brief  for balanced split of a state group
    subroutine split_state( self, which )
        use simple_ran_tabu, only: ran_tabu
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: which
        integer, allocatable :: ptcls_in_which(:)
        integer, allocatable :: states(:)
        type(ran_tabu)       :: rt
        integer              ::  n, nstates, iptcl
        nstates = self%get_n('state')
        if( which < 1 .or. which > nstates )then
            THROW_HARD('which (state) is out of range; split_state')
        endif
        call self%get_pinds(which, 'state', ptcls_in_which)
        n = size(ptcls_in_which)
        allocate(states(n), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("simple_oris::split_state",alloc_stat)
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
            THROW_HARD('which member is out of range; split_class')
        endif
        call self%get_pinds(which, 'class', ptcls_in_which)
        n = size(ptcls_in_which)
        allocate(members(n), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("simple_oris::split_class",alloc_stat)
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
        if( ncls_target <= ncls ) THROW_HARD('nr of target classes cannot be <= original number')
        ! calculate class populations
        allocate(pops(ncls_target),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: expand_classes, module: simple_oris',alloc_stat)
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

    !>  \brief  is for filling empty classes from a stochastically selected highly populated one
    subroutine fill_empty_classes( self, ncls, chunk, fromtocls)
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: ncls, chunk
        integer, allocatable, intent(out)   :: fromtocls(:,:)
        integer, allocatable :: inds2split(:), pops(:), fromtoall(:,:)
        real,    allocatable :: probs(:), corrs(:)
        logical, allocatable :: chunk_mask(:), ptcl_mask(:)
        integer              :: cnt, n_incl, pop, i, icls, ncls_here, iptcl
        integer              :: cls2split, halfpop, fromto(2)
        if(allocated(fromtocls))deallocate(fromtocls)
        ncls_here = max(self%get_n('class'), ncls)
        allocate(pops(ncls_here),fromtoall(ncls,2),source=0)
        ! chunks & populations
        allocate(ptcl_mask(self%n),chunk_mask(ncls_here),source=.false.)
        fromto(1) = huge(fromto(1))
        fromto(2) = -huge(fromto(2))
        do iptcl=1,self%n
            if( self%get(iptcl,'state') < 0.5 ) cycle
            if( nint(self%get(iptcl,'chunk')) == chunk )then
                icls = nint(self%get(iptcl,'class'))
                if( icls<1 .or. icls>ncls_here ) cycle
                ptcl_mask(iptcl) = .true.
                pops(icls)       = pops(icls)+1
                chunk_mask(icls) = .true.
                fromto(1)        = min(fromto(1),icls)
                fromto(2)        = max(fromto(2),icls)
            endif
        enddo
        if(count(pops==0) == 0)then
            deallocate(pops,fromtoall,chunk_mask,ptcl_mask)
            return
        endif
        ! temptatively assign empty classes to chunks
        chunk_mask(fromto(1):fromto(2)) = .true.
        ! randomly reassign lowly populated classes
        do icls=fromto(1),fromto(2)
            if( maxval(pops) <= 2*MINCLSPOPLIM ) exit
            if( .not.chunk_mask(icls) ) cycle
            if( pops(icls)>2 .or. pops(icls)==0 ) cycle
            call self%get_pinds(icls, 'class', inds2split, consider_w=.false.)
            do i=1,size(inds2split)
                iptcl = inds2split(i)
                cnt   = irnd_uni(ncls)
                do while( .not.chunk_mask(cnt) )
                    cnt = irnd_uni(ncls)
                enddo
                pops(cnt)  = pops(cnt)+1
                call self%set(iptcl,'class',real(cnt))
            enddo
            pops(icls) = 0
            deallocate(inds2split)
        enddo
        ! splitting
        do icls = 1,ncls_here
            if( maxval(pops) <= 2*MINCLSPOPLIM ) exit
            if( .not.chunk_mask(icls) ) cycle
            if( pops(icls) > 0 )cycle
            ! full remapping for empty classes
            chunk_mask(icls) = .false. ! exclude
            ! split class preferentially with high population
            allocate(probs(ncls))
            probs = real(pops - 2*MINCLSPOPLIM)
            where( probs<0. .or. .not.chunk_mask ) probs=0.
            probs     = probs/sum(probs)
            cls2split = multinomal(probs)
            pop       = pops(cls2split)
            if( pop <= 2*MINCLSPOPLIM )exit
            ! migration: the worst moves
            call self%get_pinds(cls2split, 'class', inds2split, consider_w=.false.)
            allocate(corrs(pop),source=-1.)
            do i=1,pop
                corrs(i) = self%get(inds2split(i),'corr')
            enddo
            call hpsort(corrs,inds2split)
            halfpop = floor(real(pop)/2.)
            do i=1,halfpop
                iptcl = inds2split(i)
                if(.not.ptcl_mask(iptcl))cycle
                ! updates populations
                call self%o(iptcl)%set('class', real(icls))
                pops(icls)      = pops(icls) + 1
                pops(cls2split) = pops(cls2split) - 1
            enddo
            ! updates populations and migration
            fromtoall(icls,1) = cls2split
            fromtoall(icls,2) = icls
            ! cleanup
            deallocate(corrs,inds2split,probs)
        enddo
        ! updates classes migration
        n_incl = count(fromtoall(:,1)>0)
        if(n_incl>0)then
            allocate(fromtocls(n_incl,2))
            cnt = 0
            do icls = 1,ncls
                if( fromtoall(icls,1) > 0 )then
                    cnt = cnt + 1
                    fromtocls(cnt,:) = fromtoall(icls,:)
                endif
            enddo
        endif
        ! cleanup
        deallocate(fromtoall,pops,ptcl_mask,chunk_mask)
    end subroutine fill_empty_classes

    !>  \brief  for remapping clusters
    subroutine remap_cls( self )
        class(oris), intent(inout) :: self
        integer :: ncls, clsind_remap, pop, icls, iptcl, old_cls
        integer , allocatable :: clspops(:)
        ncls = self%get_n('class')
        allocate(clspops(ncls),stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: remap_cls, module: simple_oris',alloc_stat)
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
    end subroutine remap_cls

    !>  \brief  is for getting an allocatable array with ptcl indices of the label 'label'
    subroutine get_pinds( self, ind, label, indices, consider_w )
        class(oris),          intent(inout) :: self
        character(len=*),     intent(in)    :: label
        integer,              intent(in)    :: ind
        integer, allocatable, intent(out)   :: indices(:)
        logical, optional,    intent(in)    :: consider_w
        integer :: pop, mystate, cnt, myval, i
        logical :: cconsider_w
        real    :: w
        cconsider_w = .false.
        if( present(consider_w) ) cconsider_w = consider_w
        if( cconsider_w )then
            if( .not. self%isthere('w') ) THROW_HARD('get_pinds with optional consider_w assumes w set')
        endif
        if( allocated(indices) )deallocate(indices)
        pop = self%get_pop(ind, label, cconsider_w )
        if( pop > 0 )then
            allocate( indices(pop), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('get_pinds; simple_oris',alloc_stat)
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
    end subroutine get_pinds

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

    !>  \brief  generate a mask for a single value of label state
    subroutine mask_from_state( self, state, l_mask, pinds, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state
        logical, allocatable, intent(out)   :: l_mask(:)
        integer, allocatable, intent(out)   :: pinds(:)
        integer, optional,    intent(in)    :: fromto(2)
        integer :: i, cnt, ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        if( allocated(l_mask) ) deallocate(l_mask)
        if( allocated(pinds)  ) deallocate(pinds)
        allocate(l_mask(ffromto(1):ffromto(2)))
        l_mask = .false.
        do i=ffromto(1),ffromto(2)
            if( self%o(i)%get_state() == state ) l_mask(i)=.true.
        end do
        allocate(pinds(1:count(l_mask)))
        cnt = 0
        do i=ffromto(1),ffromto(2)
            if( .not.l_mask(i) ) cycle
            cnt = cnt + 1
            pinds(cnt) = i
        end do
    end subroutine mask_from_state

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
        class_present = present(class)
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
            if(alloc_stat.ne.0)call allocchk('get_arr; simple_oris',alloc_stat)
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
        if( mask_present )then
            if( count(mask)==0 )then
                write(logfhandle,*)'Empty mask; simple_oris :: clac_sum'
                return
            endif
        endif
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
            if( .not. self%isthere('w') ) THROW_HARD('included with optional consider_w assumes w set')
        endif
        if(.not.allocated(incl))allocate(incl(self%n), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('In: included, module: simple_oris',alloc_stat)
        incl = .false.
        do i=1,self%n
            istate = nint(self%o(i)%get('state'))
            w = 1.0
            if( cconsider_w ) w = self%o(i)%get('w')
            if( istate > 0 .and. w > TINY ) incl(i) = .true.
        end do
    end function included

    integer function get_neven( self )
        class(oris), intent(inout) :: self
        integer :: i
        get_neven = 0
        do i=1,self%n
            if(self%o(i)%isthere('eo'))then
                if( nint(self%o(i)%get('eo'))==0) get_neven = get_neven + 1
            endif
        enddo
    end function get_neven

    !>  \brief  is getting the number of oris assigned to the odd partion
    integer function get_nodd( self )
        class(oris), intent(inout) :: self
        integer, allocatable :: eopart(:)
        eopart = nint(self%get_all('eo'))
        get_nodd = count(eopart == 1)
    end function get_nodd

    !>  \brief  is getting the number of oris assigned to the odd partion
    integer function get_nevenodd( self )
        class(oris),       intent(inout) :: self
        get_nevenodd = self%get_neven() + self%get_nodd()
    end function get_nevenodd

    subroutine print_( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        call self%o(i)%print_ori()
    end subroutine print_

    subroutine print_matrices( self )
        class(oris), intent(inout) :: self
        integer :: i
        write(logfhandle,*) 'ORDER OF ROTATION MATRIX ELEMENTS: (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)'
        do i=1,self%n
            call self%o(i)%print_mat()
        end do
    end subroutine print_matrices

    subroutine sample4update_and_incrcnt_nofrac( self, fromto, nsamples, inds, mask )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(out)   :: nsamples
        integer, allocatable, intent(out)   :: inds(:)
        logical,              intent(out)   :: mask(fromto(1):fromto(2))
        real, allocatable :: states(:)
        real    :: val
        integer :: iptcl, cnt
        ! gather info
        allocate( states(fromto(1):fromto(2)) )
        do iptcl=fromto(1),fromto(2)
            states(iptcl) = self%o(iptcl)%get('state')
        end do
        ! figure out how many samples
        nsamples = count(states > 0.5)
        ! allocate output index array
        if( allocated(inds) ) deallocate(inds)
        allocate( inds(nsamples) )
        ! update mask & count
        mask = .false.
        cnt  = 0
        do iptcl = fromto(1),fromto(2)
            if( states(iptcl) > 0.5 )then
                cnt         = cnt + 1
                inds(cnt)   = iptcl
                mask(iptcl) = .true.
                val         = self%o(iptcl)%get('updatecnt')
                call self%o(iptcl)%set('updatecnt', val+1.0)
            endif
        end do
    end subroutine sample4update_and_incrcnt_nofrac

    subroutine sample4update_and_incrcnt( self, fromto, update_frac, nsamples, inds, mask )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(out)   :: nsamples
        integer, allocatable, intent(out)   :: inds(:)
        logical,              intent(out)   :: mask(fromto(1):fromto(2))
        real,    allocatable :: counts(:), states(:)
        integer, allocatable :: inds_here(:)
        integer :: i, iptcl, cnt
        real    :: val
        ! gather info
        allocate( states(fromto(1):fromto(2)), counts(fromto(1):fromto(2)), inds_here(fromto(1):fromto(2)))
        do iptcl=fromto(1),fromto(2)
            if( self%o(iptcl)%isthere('updatecnt') )then
                counts(iptcl) = self%o(iptcl)%get('updatecnt')
            else
                counts(iptcl) = 0
            endif
            states(iptcl)    = self%o(iptcl)%get('state')
            inds_here(iptcl) = iptcl
        end do
        ! order counts
        call hpsort(counts, inds_here)
        ! figure out how many samples
        nsamples = nint(update_frac * real(count(states > 0.5)))
        ! allocate output index array
        if( allocated(inds) ) deallocate(inds)
        allocate( inds(nsamples) )
        ! update mask, index range & count
        mask = .false.
        cnt  = 0
        do i = fromto(1),fromto(2)
            iptcl = inds_here(i)
            if( states(iptcl) > 0.5 )then
                cnt         = cnt + 1
                inds(cnt)   = iptcl
                mask(iptcl) = .true.
                val         = self%o(iptcl)%get('updatecnt')
                call self%o(iptcl)%set('updatecnt', val+1.0)
            endif
            if( cnt == nsamples ) exit
        end do
    end subroutine sample4update_and_incrcnt

    subroutine sample4update_and_incrcnt2D( self, ncls, fromto, update_frac, nsamples, inds, mask )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: ncls, fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(out)   :: nsamples
        integer, allocatable, intent(out)   :: inds(:)
        logical,              intent(out)   :: mask(fromto(1):fromto(2))
        real,    allocatable :: counts(:)
        integer, allocatable :: clsarr(:)
        logical :: state_mask(fromto(1):fromto(2))
        integer :: i, cnt, icls, sz
        real    :: val, update_frac_here
        state_mask = .false.
        mask       = .false.
        do i=fromto(1),fromto(2)
            if(self%o(i)%get_state() > 0) state_mask(i) = .true.
        enddo
        update_frac_here = update_frac * real(count(state_mask)) / real(self%n)
        do icls=1,ncls
            call self%get_pinds(icls, 'class', clsarr)
            if( allocated(clsarr) )then
                sz       = size(clsarr)
                nsamples = max(1,nint(real(sz) * update_frac_here ))
                allocate(counts(sz), source=0.)
                do i=1,sz
                    if( self%o(clsarr(i))%isthere('updatecnt'))then
                        counts(i) = self%o(clsarr(i))%get('updatecnt')
                    endif
                end do
                ! order counts
                call hpsort(counts, clsarr)
                cnt = 0
                do i=1,sz
                    if( clsarr(i) >= fromto(1) .and. clsarr(i) <= fromto(2) )then
                        if( .not.state_mask(clsarr(i)) )cycle
                        cnt = cnt +1
                        if( cnt > nsamples )exit
                        mask(clsarr(i)) = .true.
                    endif
                end do
                deallocate(counts,clsarr)
            endif
        end do
        nsamples = count(mask)
        if( nsamples > 0 )then
            allocate(inds(nsamples), source=0)
            ! increment update counter and set mask
            cnt  = 0
            do i=fromto(1),fromto(2)
                if( mask(i) )then
                    cnt       = cnt + 1
                    inds(cnt) = i
                    val       = self%o(i)%get('updatecnt')
                    val       = val + 1.0
                    call self%o(i)%set('updatecnt', val)
                endif
            end do
        else
            nsamples = count(state_mask)
            allocate(inds(nsamples), source=0)
            cnt = 0
            do i=fromto(1),fromto(2)
                mask(i) = state_mask(i)
                if(mask(i))then
                    cnt       = cnt+1
                    inds(cnt) = i
                    call self%o(i)%set('updatecnt', 1.)
                endif
            end do
        endif
    end subroutine sample4update_and_incrcnt2D

    subroutine sample_rnd_subset( self, ncls, fromto, min_nsamples, nptcls_per_cls, nparts, mask, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: ncls, fromto(2), min_nsamples, nptcls_per_cls, nparts
        logical,              intent(out)   :: mask(fromto(1):fromto(2))
        integer, allocatable, intent(out)   :: inds(:)
        real    :: counts(fromto(1):fromto(2)), a
        integer :: order(1:fromto(2)-fromto(1)+1)
        integer :: i,j, nptcls_mask, nsamples, cnt
        nptcls_mask = 0
        !$omp parallel do default(shared) private(i,j) schedule(static) proc_bind(close) reduction(+:nptcls_mask)
        do i=fromto(1),fromto(2)
            j = i - fromto(1) + 1
            order(j) = j
            mask(i)  = self%o(i)%get_state() > 0
            if( mask(i) )then
                nptcls_mask = nptcls_mask + 1
                if( self%o(i)%isthere('updatecnt') )then
                    counts(i) = self%o(i)%get('updatecnt')
                else
                    counts(i) = 0.
                endif
            else
                counts(i) = huge(a)
            endif
        enddo
        !$omp end parallel do
        nsamples    = nint(max(real(min_nsamples)/real(nparts), real(nptcls_per_cls*ncls)/real(nparts)))
        if( nsamples >= nptcls_mask )then
            ! all included
        else
            call hpsort(counts,order)
            do i = nsamples+1,nptcls_mask
                j = fromto(1) + order(i) - 1
                if( mask(j) ) mask(j) = .false.
            enddo
        endif
        nsamples = count(mask)
        allocate(inds(nsamples), source=0)
        cnt = 0
        do i=fromto(1),fromto(2)
            if(mask(i))then
                cnt       = cnt+1
                inds(cnt) = i
            endif
        end do
    end subroutine sample_rnd_subset

    subroutine incr_updatecnt( self, fromto, mask )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: fromto(2)
        logical,     intent(in)    :: mask(fromto(1):fromto(2))
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=fromto(1),fromto(2)
            if( .not.mask(i) )cycle
            call self%o(i)%set('updatecnt', self%o(i)%get('updatecnt') + 1.0)
        end do
        !$omp end parallel do
    end subroutine incr_updatecnt

    !>  \brief  check wether the orientation has any typical search parameter
    logical function has_been_searched( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        has_been_searched = self%o(i)%has_been_searched()
    end function has_been_searched

    logical function any_state_zero( self )
        class(oris), intent(in) :: self
        integer :: i
        any_state_zero = .false.
        do i=1,self%n
            if( self%o(i)%get_state() == 0 )then
                any_state_zero = .true.
                return
            endif
        end do
    end function any_state_zero

    !>  \brief  joins the hashes into a string that represent the ith ori
    function ori2str( self, i ) result( str )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        character(len=:), allocatable :: str
        str = self%o(i)%ori2str()
    end function ori2str

    subroutine ori2prec( self, i, prec )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: i
        real,        intent(inout) :: prec(N_PTCL_ORIPARAMS)
        call self%o(i)%ori2prec(prec)
    end subroutine ori2prec

    subroutine prec2ori( self, i, prec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: prec(N_PTCL_ORIPARAMS)
        call self%o(i)%prec2ori(prec)
    end subroutine prec2ori

    function get_ctfvars( self, i ) result( ctfvars )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        type(ctfparams) :: ctfvars
        ctfvars = self%o(i)%get_ctfvars()
    end function get_ctfvars

    ! SETTERS

    subroutine copy( self_out, self_in, is_ptcl )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        logical,     intent(in)    :: is_ptcl
        integer   :: i
        call self_out%new(self_in%n, is_ptcl)
        do i=1,self_in%n
            self_out%o(i) = self_in%o(i)
        end do
    end subroutine copy

    subroutine reject( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        call self%o(i)%reject
    end subroutine reject

    subroutine delete_entry_1( self, key )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: key
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_entry(key)
        end do
    end subroutine delete_entry_1

    subroutine delete_entry_2( self, ind, key )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: key
        if( ind < 0 .or. ind > self%n )then
            THROW_WARN('index out of range; simple_oris % delete_entry_2')
            return
        endif
        call self%o(ind)%delete_entry(key)
    end subroutine delete_entry_2

    subroutine delete_2Dclustering( self, keepshifts )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_2Dclustering(keepshifts)
        end do
    end subroutine delete_2Dclustering

    subroutine delete_3Dalignment( self, keepshifts )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_3Dalignment(keepshifts)
        end do
    end subroutine delete_3Dalignment

    subroutine set_euler( self, i, euls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: euls(3)
        call self%o(i)%set_euler(euls)
    end subroutine set_euler

    subroutine set_shift( self, i, vec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: vec(2)
        call self%o(i)%set_shift(vec)
    end subroutine set_shift

    subroutine set_shift_incr( self, i, vec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: vec(2)
        call self%o(i)%set_shift_incr(vec)
    end subroutine set_shift_incr

    subroutine e1set( self, i, e1 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e1
        call self%o(i)%e1set(e1)
    end subroutine e1set

    subroutine e2set( self, i, e2 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e2
        call self%o(i)%e2set(e2)
    end subroutine e2set

    subroutine e3set( self, i, e3 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e3
        call self%o(i)%e3set(e3)
    end subroutine e3set

    subroutine set_1( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_1

    subroutine set_2( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_2

    subroutine set_ori( self, i, o )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: o
        self%o(i) = o
    end subroutine set_ori

    !> transfer external j-th ori j to self at i-th position
    subroutine transfer_ori( self, i, self2transfer, i2transfer )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: self2transfer
        integer,     intent(in)    :: i, i2transfer
        self%o(i) = self2transfer%o(i2transfer)
    end subroutine transfer_ori

    subroutine set_all_1( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_1

    subroutine set_all_2( self, which, vals )
        class(oris),           intent(inout) :: self
        character(len=*),      intent(in)    :: which
        character(len=STDLEN), intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_2

    subroutine set_all2single_1( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, val)
        end do
    end subroutine set_all2single_1

    subroutine set_all2single_2( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, val)
        end do
    end subroutine set_all2single_2

    subroutine set_projs( self, e_space )
        class(oris), intent(inout) :: self
        class(oris), intent(inout) :: e_space
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            call self%set(i, 'proj', real(e_space%find_closest_proj(self%o(i))))
        end do
        !$omp end parallel do
    end subroutine set_projs

    subroutine remap_projs( self, e_space, mapped_projs )
        class(oris), intent(inout) :: self
        class(oris), intent(inout) :: e_space
        integer,     intent(out)   :: mapped_projs(self%n)
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            mapped_projs(i) = e_space%find_closest_proj(self%o(i))
        end do
        !$omp end parallel do
    end subroutine remap_projs

    subroutine proj2class( self )
        class(oris), intent(inout) :: self
        integer :: i
        if( .not. self%isthere('proj') ) THROW_HARD('No proj indices to turn into class indices; proj2class')
        do i=1,self%n
            call self%o(i)%set('class', self%o(i)%get('proj'))
        end do
    end subroutine proj2class

    subroutine e3swapsgn( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%e3set(i,360.-self%e3get(i))
        end do
    end subroutine e3swapsgn

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

    subroutine zero_shifts( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', 0.)
            call self%o(i)%set('y', 0.)
        end do
    end subroutine zero_shifts

    subroutine mul_shifts( self, mul )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: mul
        integer :: i
        do i=1,self%n
            call self%o(i)%set('x', mul*self%o(i)%get('x'))
            call self%o(i)%set('y', mul*self%o(i)%get('y'))
        end do
    end subroutine mul_shifts

    subroutine rnd_oris( self, trs, eullims )
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
        integer :: i
        do i=1,self%n
            call self%o(i)%rnd_ori(trs, eullims)
        end do
    end subroutine rnd_oris

    subroutine transfer_2Dparams_1( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
        call self%o(i)%transfer_2Dparams(o_in)
    end subroutine transfer_2Dparams_1

    subroutine transfer_2Dparams_2( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
        call self%o(i)%transfer_2Dparams(os_in%o(i_in))
    end subroutine transfer_2Dparams_2

    subroutine transfer_3Dparams( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
        call self%o(i)%transfer_3Dparams(o_in)
    end subroutine transfer_3Dparams

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
            & THROW_HARD('missing angular threshold in rnd_proj_space')
        if( .not.present(o_prev) .and. present(thres) ) &
            & THROW_HARD('missing orientation in rnd_proj_space')
        within_lims = .false.
        if( present(eullims) )within_lims = .true.
        call self%new(nsample, is_ptcl=.false.)
        if( present(o_prev).and.present(thres) )then
            do i=1,self%n
                found = .false.
                do while( .not.found )
                    call o_stoch%new(is_ptcl=.false.)
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

    subroutine rnd_ori( self, i, trs, eullims )
        class(oris),    intent(inout) :: self
        integer,        intent(in)    :: i
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
        call self%o(i)%rnd_ori( trs, eullims )
    end subroutine rnd_ori

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
        if(allocated(parts))deallocate(parts)
    end subroutine ini_tseries

    !>  \brief  randomizes eulers in oris from an oris object that is assumed non-redundant
    subroutine rnd_oris_discrete_from( self, os_discrete )
        class(oris), intent(inout) :: self, os_discrete
        real       :: euls(3)
        integer    :: i, irnd, ndiscrete
        euls(3)   = 0.
        ndiscrete = os_discrete%get_noris()
        do i=1,self%n
            irnd    = irnd_uni(ndiscrete)
            euls    = os_discrete%get_euler( irnd )
            euls(3) = self%o( i )%e3get()
            call self%o(i)%set_euler( euls )
            call self%o(i)%set( 'proj', real(irnd) )
        end do
    end subroutine rnd_oris_discrete_from

    !>  \brief  randomizes the in-plane degrees of freedom
    subroutine rnd_inpls( self, trs )
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
        integer :: i
        real :: x, y
        do i=1,self%n
            if( present(trs) )then
                if( abs(trs) < TINY )then
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
            if(alloc_stat.ne.0)call allocchk("simple_oris::rnd_states",alloc_stat)
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
            THROW_HARD('invalid value for nstates; rnd_states')
        else
            ! nstates = 1; zero-preserving
            do i=1,self%n
                state = nint(self%o(i)%get('state'))
                if(  state /= 0 )call self%o(i)%set('state', 1.)
            end do
        endif
    end subroutine rnd_states

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
        type(oris) :: tmp
        type(ori)  :: o
        integer    :: cnt, i, j
        logical    :: is_ptcl
        is_ptcl = self%o(1)%is_particle()
        tmp = oris(self%get_noris()*nsym, is_ptcl)
        cnt = 0
        do i=1,self%get_noris()
            do j=1,nsym
                cnt = cnt+1
                call self%get_ori(i, o)
                call tmp%set_ori(cnt, o)
            end do
        end do
        call self%copy(tmp, is_ptcl)
        call tmp%kill
        call o%kill
    end subroutine symmetrize

    !>  \brief  for merging two oris objects into one
    subroutine merge( self1, self2add )
        class(oris), intent(inout) :: self1, self2add
        type(oris) :: self
        integer    :: ntot, cnt, i
        logical    :: is_ptcl
        is_ptcl = self1%o(1)%is_particle()
        ntot    = self1%n+self2add%n
        self    = oris(ntot, is_ptcl)
        cnt     = 0
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
        call self1%copy(self, is_ptcl)
        call self%kill
        call self2add%kill
    end subroutine merge

    subroutine clean_updatecnt( self )
        class(oris),       intent(inout) :: self
        integer :: i
        do i = 1,self%n
            call self%o(i)%delete_entry('updatecnt')
        enddo
    end subroutine clean_updatecnt

    !>  \brief  for balanced assignment of even/odd partitions
    subroutine partition_eo( self )
        class(oris),       intent(inout) :: self    !< instance
        integer :: i
        do i = 1, self%n-1, 2
            call self%set(i,  'eo', 0.)
            call self%set(i+1,'eo', 1.)
        end do
        if(is_even(self%n))then
            call self%set(self%n,'eo',1.)
        else
            call self%set(self%n,'eo',0.)
        endif
    end subroutine partition_eo

    subroutine str2ori( self, i, line )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(inout) :: line
        call self%o(i)%str2ori(line, self%o(1)%is_particle())
    end subroutine str2ori

    subroutine str2ori_ctfparams_state_eo( self, i, line )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(inout) :: line
        type(ori) :: o_tmp
        call o_tmp%str2ori(line, self%o(1)%is_particle())
        if( o_tmp%isthere('smpd')    ) call self%o(i)%set('smpd',    o_tmp%get('smpd'))
        if( o_tmp%isthere('kv')      ) call self%o(i)%set('kv',      o_tmp%get('kv'))
        if( o_tmp%isthere('cs')      ) call self%o(i)%set('cs',      o_tmp%get('cs'))
        if( o_tmp%isthere('fraca')   ) call self%o(i)%set('fraca',   o_tmp%get('fraca'))
        if( o_tmp%isthere('phshift') ) call self%o(i)%set('phshift', o_tmp%get('phshift'))
        if( o_tmp%isthere('dfx')     ) call self%o(i)%set('dfx',     o_tmp%get('dfx'))
        if( o_tmp%isthere('dfy')     ) call self%o(i)%set('dfy',     o_tmp%get('dfy'))
        if( o_tmp%isthere('angast')  ) call self%o(i)%set('angast',  o_tmp%get('angast'))
        if( o_tmp%isthere('state')   )then
            call self%o(i)%set('state', o_tmp%get('state'))
        else
            call self%o(i)%set('state', 1.0)
        endif
        if( o_tmp%isthere('eo') ) call self%o(i)%set('eo',      o_tmp%get('eo'))
        call o_tmp%kill
    end subroutine str2ori_ctfparams_state_eo

    subroutine set_ctfvars( self, i, ctfvars )
        class(oris),     intent(inout) :: self
        integer,         intent(in)    :: i
        type(ctfparams), intent(in)    :: ctfvars
        call self%o(i)%set_ctfvars(ctfvars)
    end subroutine set_ctfvars

    subroutine set_boxfile( self, i, boxfname, nptcls )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        character(len=*),  intent(in)    :: boxfname
        integer, optional, intent(in)    :: nptcls
        call self%o(i)%set_boxfile(boxfname, nptcls)
    end subroutine set_boxfile

    subroutine reset( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i = 1,self%n
            call self%o(i)%kill_hash
            call self%o(i)%kill_chash
            call self%o(i)%reset_pparms
        enddo
    end subroutine reset

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
            THROW_HARD("the file you are trying to read: "//trim(orifile)//' does not exist in cwd' )
        endif
        if( trim(fname2ext(orifile)) == 'bin' )then
            THROW_HARD('this method does not support binary files; read')
        endif
        io_message='No error'
        call fopen(fnr, FILE=orifile, STATUS='OLD', action='READ', iostat=file_stat,iomsg=io_message)
        call fileiochk("oris ; read ,Error when opening file for reading: "//trim(orifile)//':'//trim(io_message), file_stat)
        if( present(nst) ) nst = 0
        if( present(fromto) )then
            istart = fromto(1)
            iend   = fromto(2)
            if(istart <      1) THROW_HARD('Invalid index; read')
            if(iend   > self%n) THROW_HARD('Invalid index; read')
        else
            istart = 1
            iend   = self%n
        endif
        do i = istart, iend
            call self%o(i)%read(fnr)
            if( present(nst) )then
                state = self%o(i)%get_state()
                nst   = max(1,max(state,nst))
            endif
        end do
        call fclose(fnr,errmsg="oris ; read ,Error when closing file")
    end subroutine read

    !>  \brief  reads CTF parameters and state info from file
    subroutine read_ctfparams_state_eo( self, ctfparamfile )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: ctfparamfile
        logical    :: params_are_there(10)
        integer    :: i
        type(oris) :: os_tmp
        if( .not. file_exists(ctfparamfile) )then
            THROW_HARD ("read_ctfparams_state_eo; The file you are trying to read: "//trim(ctfparamfile)//' does not exist')
        endif
        if( str_has_substr(ctfparamfile,'.bin') )then
            THROW_HARD('this method does not support binary files; read_ctfparams_state_eo')
        endif
        call os_tmp%new(self%n, self%o(1)%is_particle())
        call os_tmp%read(ctfparamfile)
        params_are_there(1)  = os_tmp%isthere('smpd')
        params_are_there(2)  = os_tmp%isthere('kv')
        params_are_there(3)  = os_tmp%isthere('cs')
        params_are_there(4)  = os_tmp%isthere('fraca')
        params_are_there(5)  = os_tmp%isthere('phshift')
        params_are_there(6)  = os_tmp%isthere('dfx')
        params_are_there(7)  = os_tmp%isthere('dfy')
        params_are_there(8)  = os_tmp%isthere('angast')
        params_are_there(9)  = os_tmp%isthere('state')
        params_are_there(10) = os_tmp%isthere('eo')
        do i=1,self%n
            if( params_are_there(1) )  call self%set(i, 'smpd',    os_tmp%get(i, 'smpd')   )
            if( params_are_there(2) )  call self%set(i, 'kv',      os_tmp%get(i, 'kv')     )
            if( params_are_there(3) )  call self%set(i, 'cs',      os_tmp%get(i, 'cs')     )
            if( params_are_there(4) )  call self%set(i, 'fraca',   os_tmp%get(i, 'fraca')  )
            if( params_are_there(5) )  call self%set(i, 'phshift', os_tmp%get(i, 'phshift'))
            if( params_are_there(6) )  call self%set(i, 'dfx',     os_tmp%get(i, 'dfx')    )
            if( params_are_there(7) )  call self%set(i, 'dfy',     os_tmp%get(i, 'dfy')    )
            if( params_are_there(8) )  call self%set(i, 'angast',  os_tmp%get(i, 'angast') )
            if( params_are_there(9) )  call self%set(i, 'state',   os_tmp%get(i, 'state')  )
            if( params_are_there(10) ) call self%set(i, 'eo',      os_tmp%get(i, 'eo')     )
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
        call fileiochk(' Error opening file for writing: '//trim(orifile)//' ; '//trim(io_message), file_stat)
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
        call fileiochk( 'In: write_2, module: simple_oris.f90  opening '//trim(orifile), file_stat )
        call self%o(i)%write(fnr)
        call fclose(fnr, errmsg=' Error closing file for writing: '//trim(orifile))
    end subroutine write_2

    !>  \brief  writes object to BILD Chimera readable format
    subroutine write2bild( self, file  )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: file
        integer :: i,funit, file_stat
        call fopen(funit, file, status='REPLACE', action='WRITE',iostat=file_stat)
        call fileiochk( 'In: write2bild, module: simple_oris.f90  opening '//trim(file), file_stat )
        ! header
        write(funit,'(A)')".translate 0.0 0.0 0.0"
        write(funit,'(A)')".scale 10"
        write(funit,'(A)')".comment -- unit sphere --"
        write(funit,'(A)')".color 0.8 0.8 0.8"
        write(funit,'(A)')".sphere 0 0 0 1.0"
        write(funit,'(A)')".comment -- planes --"
        write(funit,'(A)')".color 0.3 0.3 0.3"
        write(funit,'(A)')".cylinder -0.02 0 0 0.02 0 0 1.02"
        write(funit,'(A)')".cylinder 0 -0.02 0 0 0.02 0 1.02"
        write(funit,'(A)')".cylinder 0 0 -0.02 0 0 0.02 1.02"
        write(funit,'(A)')".comment -- x-axis --"
        write(funit,'(A)')".color 1 0 0"
        write(funit,'(A)')".cylinder -1.5 0 0 1.5 0 0 0.02"
        write(funit,'(A)')".comment -- y-axis --"
        write(funit,'(A)')".color 0 1 0"
        write(funit,'(A)')".cylinder 0 -1.5 0 0 1.5 0 0.02"
        write(funit,'(A)')".comment -- z-axis --"
        write(funit,'(A)')".color 0 0 1"
        write(funit,'(A)')".cylinder 0 0 -1.5 0 0 1.5 0.02"
        write(funit,'(A)')".comment -- north pole --"
        write(funit,'(A)')".color 0 0 1"
        write(funit,'(A)')".sphere 0 0 1.5 0.1"
        ! body
        write(funit,'(A)')".color 0.4 0.4 0.4"
        do i=1,self%n
            call self%o(i)%write2bild(funit)
        enddo
        call fclose(funit, errmsg=' Error closing file for writing: '//trim(file))
    end subroutine write2bild

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
             call e%compose(self%o(i), o_tmp)
             call self%o(i)%set_euler(o_tmp%get_euler())
         end do
         call o_tmp%kill
    end subroutine rot_1

    !>  \brief  is an Euler angle composer
    subroutine rot_2( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
        type(ori)                  :: o_tmp
        call e%compose(self%o(i), o_tmp)
        call self%o(i)%set_euler(o_tmp%get_euler())
        call o_tmp%kill
    end subroutine rot_2

    !>  \brief  is an Euler angle composer
    subroutine rot_transp_1( self, e )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: e
        type(ori) :: o_tmp, e_transp
        integer   :: i
        e_transp = e
        call e_transp%transp
        do i=1,self%n
             call e_transp%compose(self%o(i), o_tmp)
             call self%o(i)%set_euler(o_tmp%get_euler())
         end do
         call o_tmp%kill
         call e_transp%kill
    end subroutine rot_transp_1

    !>  \brief  is an Euler angle composer
    subroutine rot_transp_2( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
        type(ori) :: o_tmp, e_transp
        e_transp = e
        call e_transp%transp
        call e_transp%compose(self%o(i), o_tmp)
        call self%o(i)%set_euler(o_tmp%get_euler())
        call o_tmp%kill
        call e_transp%kill
    end subroutine rot_transp_2

    !>  \brief  for identifying the median value of parameter which
    function median_1( self, which ) result( med )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,    allocatable :: vals(:), vals4med(:)
        logical, allocatable :: incl(:)
        real :: med
        incl     = self%included()
        vals     = self%get_all(which)
        vals4med = pack(vals, incl)
        med      = median_nocopy(vals)
    end function median_1

    !>  \brief  is for calculating variable statistics
    subroutine stats_1( self, which, ave, sdev, var, err )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        real, allocatable :: vals(:), states(:)
        states = self%get_all('state')
        vals   = self%get_all(which)
        call moment(vals, ave, sdev, var, err, states > 0.5)
        deallocate(vals, states)
    end subroutine stats_1

    !>  \brief  is for calculating variable statistics
    subroutine stats_2( self, which, statvars, mask )
        class(oris),        intent(inout) :: self
        character(len=*),   intent(in)    :: which
        type(stats_struct), intent(out)   :: statvars
        logical, optional,  intent(in)    :: mask(self%n)
        real, allocatable :: vals(:), states(:)
        logical :: err, mmask(self%n)
        real    :: var
        if( present(mask) )then
            mmask = mask
        else
            mmask = .true.
        endif
        states = self%get_all('state')
        mmask  = mmask .and. states > 0.5
        vals   = self%get_all(which)
        call moment(vals, statvars%avg, statvars%sdev, var, err, mmask)
        statvars%minv = minval(vals, mask=mmask)
        statvars%maxv = maxval(vals, mask=mmask)
        deallocate(states, vals)
    end subroutine stats_2

    !>  \brief  is for calculating the minimum/maximum values of a variable
    subroutine minmax( self, which, minv, maxv )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: minv, maxv
        real    :: val, x
        integer :: i, mystate
        minv = huge(x)
        maxv = -huge(x)
        do i=1,self%n
            mystate = nint(self%o(i)%get('state'))
            if( mystate /= 0 )then
                val = self%o(i)%get(which)
                if( val < minv ) minv = val
                if( val > maxv ) maxv = val
            endif
        end do
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
            THROW_HARD('object nonexistent; spiral_1')
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
        deallocate(avail)

        contains

            subroutine gen_c1
                integer :: i
                if( allocated(avail) )deallocate(avail, stat=alloc_stat)
                allocate(avail(n), source=.false., stat=alloc_stat)
                call tmp%new(n, self%o(1)%is_particle())
                call tmp%spiral_1
                do i = 1, n
                    if( tmp%o(i)%e1get() <= e1lim .and. tmp%o(i)%e2get() <= e2lim )&
                    &avail(i) = .true.
                end do
            end subroutine gen_c1

    end subroutine spiral_2

    !>  \brief  orders oris according to specscore
    function order( self ) result( inds )
        class(oris), intent(inout) :: self
        real,    allocatable :: specscores(:)
        integer, allocatable :: inds(:)
        integer :: i
        allocate( inds(self%n), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('order; simple_oris',alloc_stat)
        specscores = self%get_all('specscore')
        inds = (/(i,i=1,self%n)/)
        call hpsort(specscores, inds)
        call reverse(inds)
        if(allocated(specscores))then
            deallocate( specscores, stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('order; simple_oris dealloc ',alloc_stat)
        end if
    end function order

    !>  \brief  orders oris according to corr
    function order_corr( self ) result( inds )
        class(oris), intent(inout) :: self
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        integer :: i
        allocate( inds(self%n), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('order; simple_oris',alloc_stat)
        corrs = self%get_all('corr')
        inds = (/(i,i=1,self%n)/)
        call hpsort(corrs, inds)
        call reverse(inds)
    end function order_corr

    !>  \brief  orders clusters according to population
    function order_cls( self, ncls ) result( inds )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls
        integer, allocatable :: inds(:)
        real    :: classpops(ncls)
        integer :: i
        if(ncls <= 0) THROW_HARD('invalid number of classes; order_cls')
        allocate(inds(ncls), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('order_cls; simple_oris',alloc_stat)
        classpops = 0.0
        ! calculate class populations
        do i=1,ncls
            classpops(i) = real(self%get_pop(i, 'class'))
        end do
        inds = (/(i,i=1,ncls)/)
        call hpsort(classpops, inds)
        call reverse(inds)
    end function order_cls

    subroutine find_best_classes( self, box, smpd, res_thresh,cls_mask, ndev )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: box
        real,        intent(in)    :: smpd, res_thresh, ndev
        logical,     intent(inout) :: cls_mask(1:self%n)
        real,    allocatable :: rfinds(:), corrs(:)
        logical, allocatable :: msk(:)
        real    :: ave, sdev, res, res_threshold, corr_threshold
        integer :: icls, nincl
        logical :: has_res, has_corr
        allocate(msk(self%n), source=.true.)
        if( self%isthere('pop') )then
            do icls=1,self%n
                if(self%get(icls,'pop')<0.5) msk(icls) = .false.
            enddo
        endif
        has_corr = self%isthere('corr')
        if( has_corr )then
            do icls=1,self%n
                if(self%get(icls,'corr') < 0.0001) msk(icls) = .false.
            enddo
        endif
        has_res = self%isthere('res')
        if( has_res )then
            do icls=1,self%n
                if(self%get(icls,'res') > res_thresh) msk(icls) = .false.
            enddo
        endif
        nincl    = count(msk)
        cls_mask = msk
        if( has_res )then
            allocate(rfinds(self%n), source=0.)
            do icls=1,self%n
                res = self%get(icls, 'res')
                rfinds(icls) = real(calc_fourier_index(res,box,smpd))
                call self%set(icls,'find',rfinds(icls))
            enddo
            ave  = sum(rfinds,mask=msk)/real(nincl)
            sdev = sqrt(sum((rfinds-ave)**2.,mask=msk)/real(nincl))
            res_threshold = max(ave-ndev*sdev,2.)
        else
            allocate(rfinds(self%n), source=huge(res_threshold))
            res_threshold = 0.
        endif
        if( has_corr )then
            corrs = self%get_all('corr')
            ave   = sum(corrs,mask=msk) / real(nincl)
            sdev  = sqrt(sum((corrs-ave)**2., mask=msk)/real(nincl))
            corr_threshold = ave-ndev*sdev
        else
            allocate(corrs(self%n), source=huge(corr_threshold))
            corr_threshold = 0.
        endif
        do icls=1,self%n
            if( cls_mask(icls) )then
                if(rfinds(icls)<res_threshold .and. corrs(icls)<corr_threshold)then
                    cls_mask(icls) = .false.
                endif
            endif
        enddo
        deallocate(msk,rfinds,corrs)
    end subroutine find_best_classes

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

    !>  \brief  calculates soft weights based on specscore
    subroutine calc_soft_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        real,    allocatable :: specscores(:), states(:), weights(:), weights_glob(:)
        integer, allocatable :: order(:), states_int(:)
        real    :: minw, mins
        integer :: i, lim, state, nstates
        if( self%isthere('specscore') )then
            specscores = self%get_all('specscore')
            mins       = minval(specscores, mask=specscores > TINY)
            if( mins >= 0.8 )then
                call self%set_all2single('w', 1.0)
                return
            endif
            if( self%isthere('states') )then
                states = self%get_all('states')
            else
                allocate(states(self%n), source=1.0)
            endif
            nstates = nint(maxval(states))
            if( nstates == 1 )then
                weights = z_scores(specscores, mask=specscores > TINY .and. states > 0.5)
                minw    = minval(weights,      mask=specscores > TINY .and. states > 0.5)
                where( specscores > TINY .and. states > 0.5 )
                    weights = weights + abs(minw)
                elsewhere
                    weights = 0. ! nuke
                endwhere
                call self%set_all('w', weights)
                deallocate(weights)
            else
                allocate(states_int(self%n),   source=nint(states))
                allocate(weights_glob(self%n), source=0.)
                do state=1,nstates
                    weights = z_scores(specscores, mask=specscores > TINY .and. states_int == state)
                    minw    = minval(weights,      mask=specscores > TINY .and. states_int == state)
                    where( specscores > TINY .and. states_int == state ) weights_glob = weights + abs(minw)
                    deallocate(weights)
                end do
                call self%set_all('w', weights_glob)
                deallocate(states_int, weights_glob)
            endif
            if( frac < 0.99 )then
                ! in 3D frac operates globally, independent of state
                lim   = nint(frac*real(self%n))
                order = self%order() ! specscore ranking
                do i=1,self%n
                    if( i > lim ) call self%o(order(i))%set('w', 0.) ! nuke
                end do
                deallocate(order)
            endif
            deallocate(specscores, states)
        else
            call self%set_all2single('w', 1.0)
        endif
    end subroutine calc_soft_weights

    !>  \brief  calculates hard weights based on ptcl ranking
    subroutine calc_hard_weights2D( self, frac, ncls )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        integer,     intent(in)    :: ncls
        type(oris)           :: os
        type(ori)            :: o
        integer, allocatable :: pinds(:)
        integer :: i, icls, pop
        if( frac < 0.99 )then
            do icls=1,ncls
                call self%get_pinds(icls, 'class', pinds, consider_w=.false.)
                if(.not.allocated(pinds)) cycle
                pop = size(pinds)
                call os%new(pop, self%o(1)%is_particle())
                do i=1,pop
                    call self%get_ori(pinds(i), o)
                    call os%set_ori(i, o)
                enddo
                call os%calc_hard_weights(frac)
                do i=1,pop
                    call self%set(pinds(i), 'w', os%get(i,'w'))
                enddo
                deallocate(pinds)
            enddo
        else
            call self%set_all2single('w', 1.)
        endif
        call o%kill
        call os%kill
    end subroutine calc_hard_weights2D

    !>  \brief  calculates soft weights based on specscore
    subroutine calc_soft_weights2D( self )
        class(oris), intent(inout) :: self
        real,    allocatable :: specscores(:), states(:)
        real,    allocatable :: weights(:), weights_glob(:)
        integer, allocatable :: classes(:)
        integer :: icls, pop, ncls
        real    :: minw
        if( self%isthere('specscore') )then
            specscores = self%get_all('specscore')
            classes    = nint(self%get_all('class'))
            if( self%isthere('states') )then
                states = self%get_all('states')
            else
                allocate(states(self%n), source=1.0)
            endif
            allocate(weights_glob(self%n), source=0.)
            ncls = maxval(classes)
            do icls=1,ncls
                pop = count(classes == icls .and. states > 0.5)
                if( pop == 0 )then
                    cycle
                else if( pop <= MINCLSPOPLIM )then
                    where( specscores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = 1.0
                else
                    if( count(specscores > TINY .and. (classes == icls .and. states > 0.5)) <= MINCLSPOPLIM )then
                        where( specscores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = 1.0
                    else
                        weights = z_scores(specscores, mask=specscores > TINY .and. (classes == icls .and. states > 0.5))
                        minw    = minval(weights,      mask=specscores > TINY .and. (classes == icls .and. states > 0.5))
                        where( specscores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = weights + abs(minw)
                        deallocate(weights)
                    endif
                endif
            end do
            call self%set_all('w', weights_glob)
            deallocate(specscores, classes, states, weights_glob)
        else
            call self%set_all2single('w', 1.)
        endif
    end subroutine calc_soft_weights2D

    !>  \brief  to find the closest matching projection direction
    !! KEEP THIS ROUTINE SERIAL
    function find_closest_proj( self, o_in ) result( closest )
        class(oris), intent(inout) :: self
        class(ori),  intent(in) :: o_in
        real    :: dists(self%n)
        integer :: closest, i
        do i=1,self%n
            dists(i)=self%o(i).euldist.o_in
        end do
        closest = minloc( dists, dim=1 )
    end function find_closest_proj

    !>  \brief  method for discretization of the projection directions
    subroutine discretize( self, n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        type(oris) :: d
        integer    :: closest, i
        if( n < self%n )then
            call d%new(n, self%o(1)%is_particle())
            call d%spiral
            do i=1,self%n
                closest = d%find_closest_proj(self%o(i))
                call self%o(i)%e1set(d%e1get(closest))
                call self%o(i)%e2set(d%e2get(closest))
                call self%o(i)%set('class', real(closest))
            end do
        else
            THROW_HARD('the number of discrete oris is too large; discretize')
        endif
    end subroutine discretize

    !>  \brief  to identify the indices of the k nearest projection neighbors (inclusive)
    subroutine nearest_proj_neighbors_1( self, k, nnmat )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: k
        integer,     intent(inout) :: nnmat(self%n,k)
        real      :: dists(self%n)
        integer   :: inds(self%n), i, j
        type(ori) :: o
        if( k >= self%n ) THROW_HARD('need to identify fewer nearest_proj_neighbors')
        do i=1,self%n
            call self%get_ori(i, o)
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
        call o%kill
    end subroutine nearest_proj_neighbors_1

    !>  \brief  to identify the k nearest projection neighbors (exclusive), returned as logical array
    !!          self is search space with finer angular resolution
    subroutine nearest_proj_neighbors_2( self, o_in, k, lnns )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_in
        integer,     intent(in)    :: k
        logical,     intent(inout) :: lnns(self%n)
        real      :: dists(self%n)
        integer   :: inds(self%n), i, j
        type(ori) :: o
        lnns  = .false.
        do i=1,self%n
            inds(i) = i
            call self%get_ori(i, o)
            dists(i) = o.euldist.o_in
        end do
        call hpsort(dists, inds)
        do j=1,k
            lnns(inds(j)) = .true.
        end do
    end subroutine nearest_proj_neighbors_2

    !>  \brief  to identify the nearest projection neighbors based on euldist threshold
    subroutine nearest_proj_neighbors_3( self, o, euldist_thres, lnns )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o
        real,        intent(in)    :: euldist_thres ! in degrees
        logical,     intent(inout) :: lnns(self%n)
        real      :: dists(self%n), euldist_thres_rad
        integer   :: inds(self%n), j
        euldist_thres_rad = deg2rad(euldist_thres)
        lnns  = .false.
        do j=1,self%n
            inds(j)  = j
            dists(j) = self%o(j).euldist.o
        end do
        call hpsort(dists, inds)
        do j=1,self%n
            if( dists(j) <= euldist_thres_rad ) lnns(inds(j)) = .true.
        end do
    end subroutine nearest_proj_neighbors_3

    !>  \brief  to identify the nearest projection neighbors based on nearest neigh matrix and euldist threshold
    subroutine nearest_proj_neighbors_4( self, os_cls, icls, nnn, nnmat, euldist_thres, lnns )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: os_cls
        integer,     intent(in)    :: icls, nnn, nnmat(os_cls%n,nnn)
        real,        intent(in)    :: euldist_thres ! in degrees
        logical,     intent(inout) :: lnns(self%n)
        real      :: dists(self%n), euldist_thres_rad
        integer   :: inds(self%n), j, k
        euldist_thres_rad = deg2rad(euldist_thres)
        lnns  = .false.
        do k=1,nnn
            do j=1,self%n
                inds(j)  = j
                dists(j) = self%o(j).euldist.os_cls%o(nnmat(icls,k))
            end do
            call hpsort(dists, inds)
            do j=1,self%n
                if( dists(j) <= euldist_thres_rad ) lnns(inds(j)) = .true.
            end do
        end do
    end subroutine nearest_proj_neighbors_4

    subroutine min_euldist( self, o_in, mindist )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_in
        real,        intent(inout) :: mindist
        real, allocatable :: ows(:)
        real      :: dists(self%n), x
        integer   :: inds(self%n), i, loc(1)
        type(ori) :: o
        dists = huge(x)
        do i=1,self%n
            inds(i) = i
            call self%get_ori(i, o)
            dists(i) = o.euldist.o_in
        end do
        loc = minloc(dists)
        mindist = rad2deg(dists(loc(1)))
    end subroutine min_euldist

    !>  \brief  to find angular resolution of an even orientation distribution (in degrees)
    function find_angres( self ) result( res )
        class(oris), intent(in) :: self
        real    :: dists(self%n), res, x, nearest3(3)
        integer :: i, j
        res = 0.
        !$omp parallel do default(shared) proc_bind(close)&
        !$omp private(j,i,dists,nearest3) reduction(+:res)
        do j=1,self%n
            do i=1,self%n
                if( i == j )then
                    dists(i) = huge(x)
                else
                    dists(i) = self%o(i).euldist.self%o(j)
                endif
            end do
            nearest3 = min3(dists)
            res = res + sum(nearest3) / 3. ! average of three nearest neighbors
        end do
        !$omp end parallel do
        res = rad2deg(res / real(self%n))
    end function find_angres

    !>  \brief  to find the correlation bound in extremal search
    function extremal_bound( self, thresh, which ) result( score_bound )
        class(oris),             intent(inout) :: self
        real,                    intent(in)    :: thresh ! is a fraction
        character(len=*), optional, intent(in) :: which
        real,    allocatable  :: scores(:), scores_incl(:)
        logical, allocatable  :: incl(:)
        character(len=KEYLEN) :: which_here
        integer :: n_incl, thresh_ind
        real    :: score_bound
        if( present(which) )then
            which_here = trim(which)
        else
            which_here = 'corr'
        endif
        select case(trim(which_here))
            case('corr')
                ! to use in conjunction with objfun=cc
            case('specscore')
                ! to be tested, should be amenable to any objective function
            case DEFAULT
                write(logfhandle,*)'Invalid metric: ', trim(which_here)
                THROW_HARD('extremal_bound')
        end select
        if( .not.self%isthere(which_here) )then
            THROW_HARD('Metric: '//trim(which_here)//' is unpopulated; extremal_bound')
        endif
        ! fetch scores
        scores      = self%get_all(which_here)
        incl        = self%included()
        scores_incl = pack(scores, mask=incl)
        ! sort & determine threshold
        n_incl      = size(scores_incl)
        call hpsort(scores_incl)
        thresh_ind  = nint(real(n_incl) * thresh)
        score_bound = scores_incl(thresh_ind)
        deallocate(scores, incl, scores_incl)
    end function extremal_bound

    !>  \brief  utility function for setting extremal optimization parameters
    subroutine set_extremal_vars(self, extr_init, extr_iter, iter, frac_srch_space,&
            &do_extr, iextr_lim, update_frac)
        class(oris),       intent(in)  :: self
        real,              intent(in)  :: extr_init
        integer,           intent(in)  :: extr_iter, iter
        real,              intent(in)  :: frac_srch_space
        logical,           intent(out) :: do_extr
        integer,           intent(out) :: iextr_lim
        real,    optional, intent(in)  :: update_frac
        integer :: zero_pop
        logical :: l_update_frac
        do_extr           = .false.
        l_update_frac     = .false.
        if( present(update_frac) )then
            if( update_frac > 0.001 ) l_update_frac = .true.
        endif
        zero_pop  = self%n - self%get_noris(consider_state=.true.)
        iextr_lim = ceiling(2.*log(real(self%n-zero_pop)) * (2.-update_frac))
        if( l_update_frac )then
            iextr_lim = ceiling(2.*log(real(self%n-zero_pop)) * (2.-update_frac))
            if(iter==1 .or.(frac_srch_space <= 99. .and. extr_iter <= iextr_lim)) do_extr = .true.
        else
            iextr_lim = ceiling(2.*log(real(self%n-zero_pop)))
            if(iter==1 .or.(frac_srch_space <= 98. .and. extr_iter <= iextr_lim)) do_extr = .true.
        endif
    end subroutine set_extremal_vars

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

    !>  \brief  generates the mirror of the projection
    !!          can be equivalently accomplished by mirror('y')
    !!          the image after projection
    subroutine mirror2d( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i = 1, self%n
            call self%o(i)%mirror2d()
        enddo
    end subroutine mirror2d

    !>  \brief  generates the opposite hand of an Euler angle
    !!          so that a set of Euler angles transformed by
    !!          this operation changes the handedness of the volume
    subroutine mirror3d( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i = 1, self%n
            call self%o(i)%mirror3d()
        enddo
    end subroutine mirror3d

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

    !>  \brief  for calculating statistics of distances within a single distribution
    subroutine diststat_1( self, sumd, avgd, sdevd, mind, maxd )
        class(oris), intent(in)  :: self
        real,        intent(out) :: mind, maxd, avgd, sdevd, sumd
        integer :: i, j, cnt
        real    :: dists((self%n*(self%n-1))/2), vard
        logical :: err
        cnt  = 0
        do i=1,self%n-1
           do j=i+1,self%n
              cnt = cnt+1
              dists(cnt) = self%o(i).euldist.self%o(j)
           end do
        end do
        mind = minval(dists)
        maxd = maxval(dists)
        sumd = sum(dists)
        call moment(dists, avgd, sdevd, vard, err )
    end subroutine diststat_1

    !>  \brief  for calculating statistics of distances between two equally sized distributions
    subroutine diststat_2( self1, self2, sumd, avgd, sdevd, mind, maxd )
        class(oris), intent(inout)  :: self1, self2
        real,        intent(out) :: mind, maxd, avgd, sdevd, sumd
        real, allocatable :: onormals1(:,:),onormals2(:,:)
        real    :: dists(self1%n), vard, x
        integer :: i
        logical :: err
        if( self1%n /= self2%n )then
            THROW_HARD('cannot calculate distance between sets of different size; euldist_2')
        endif
        mind = huge(x)
        maxd = -mind
        sumd = 0.
        onormals1 = self1%get_all_normals()
        onormals2 = self2%get_all_normals()
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)&
        !$omp reduction(+:sumd) reduction(min:mind) reduction(max:maxd)
         do i=1,self1%n
            dists(i) = vector_angle_norm(onormals1(i,:),onormals2(i,:)) ! (self1%o(i).euldist.self2%o(i))
            if( dists(i) < mind ) mind = dists(i)
            if( dists(i) > maxd ) maxd = dists(i)
            sumd = sumd+dists(i)
        end do
        !$omp end parallel do
        call moment(dists, avgd, sdevd, vard, err )
        deallocate(onormals1,onormals2)
    end subroutine diststat_2

    real function overlap( self1, self2, which, state )
        class(oris),      intent(inout) :: self1, self2
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: state
        real,    allocatable :: arr1(:), arr2(:), tmp(:), ows(:), states(:)
        logical, allocatable :: mask1(:), mask2(:)
        integer :: n1, n2, sz, n_min, i
        overlap = 0.
        if( self1%n == 0 .or. self2%n == 0 ) return
        if( .not. self1%isthere(trim(which)) )then
            THROW_WARN('key: '//trim(which)//' not present in self1; overlap')
            return
        endif
        if( .not. self2%isthere(trim(which)) )then
            THROW_WARN('key: '//trim(which)//' not present in self2; overlap')
            return
        endif
        if( .not. self1%isthere('state') )     THROW_HARD('key: state not present in self1; overlap')
        if( .not. self2%isthere('state') )     THROW_HARD('key: state not present in self2; overlap')
        if( .not. self1%isthere('ow') )        THROW_HARD('key: ow not present in self1; overlap')
        if( .not. self2%isthere('ow') )        THROW_HARD('key: ow not present in self2; overlap')
        ! extract self1 which labels
        ows    = self1%get_all('ow')
        states = self1%get_all('state')
        where(nint(states) .ne. state) ows = 0.
        if( .not. any(ows > TINY) ) return
        tmp    = self1%get_all(trim(which))
        arr1   = pack(tmp, mask=ows > TINY)
        n1     = size(arr1)
        ! extract self2 which labels
        ows    = self2%get_all('ow')
        states = self2%get_all('state')
        where(nint(states) .ne. state) ows = 0.
        if( .not. any(ows > TINY) ) return
        tmp    = self2%get_all(trim(which))
        arr2   = pack(tmp, mask=ows > TINY)
        n2     = size(arr2)
        ! translate to masks (this prevents counting duplicates)
        sz = nint(max(maxval(arr1), maxval(arr2)))
        allocate(mask1(sz), mask2(sz), source=.false.)
        forall(i=1:n1) mask1(nint(arr1(i))) = .true.
        forall(i=1:n2) mask2(nint(arr2(i))) = .true.
        ! compare and normalize
        n_min   = min(count(mask1),count(mask2))
        overlap = real(count(mask1 .and. mask2)) / real(n_min)
    end function overlap

    pure real function geodesic_distance( rmat1, rmat2 )
        real, intent(in) :: rmat1(3,3), rmat2(3,3)
        real :: Imat(3,3), sumsq, diffmat(3,3)
        Imat      = 0.
        Imat(1,1) = 1.
        Imat(2,2) = 1.
        Imat(3,3) = 1.
        diffmat = Imat - matmul(rmat1,transpose(rmat2))
        sumsq   = sum(diffmat*diffmat)
        if( sumsq > 0.0001 )then
            geodesic_distance = sqrt(sumsq)
        else
            geodesic_distance = 0.
        endif
    end function geodesic_distance

    pure real function geodesic_scaled_dist( rmat1, rmat2 )
        real, intent(in) :: rmat1(3,3), rmat2(3,3)
        real, parameter :: old_max = 2.*sqrt(2.)
        geodesic_scaled_dist = geodesic_distance(rmat1,rmat2)*(pi/old_max)
    end function geodesic_scaled_dist

    ! UNIT TEST

    !>  \brief  oris class unit test
    subroutine test_oris( doprint )
        logical, intent(in)  :: doprint
        type(oris)           :: os, os2
        real                 :: euls(3), corr, x, x2, y, y2
        integer              :: i
        integer, allocatable :: order(:)
        logical              :: passed
        write(logfhandle,'(a)') '**info(simple_oris_unit_test, part1): testing getters/setters'
        os  = oris(100, is_ptcl=.false.)
        os2 = oris(100, is_ptcl=.false.)
        passed = .false.
        if( os%get_noris() == 100 ) passed = .true.
        if( .not. passed ) THROW_HARD('get_noris failed!')
        passed = .false.
        call os%set_euler(1, [1.,2.,3.])
        euls = os%get_euler(1)
        if( abs(euls(1)-1.+euls(2)-2.+euls(3)-3.) < 0.0001 ) passed = .true.
        if( .not. passed ) THROW_HARD('get/set eulers failed!')
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
            write(logfhandle,*) '********'
            do i=1,100
                call os%print_(i)
            end do
            call os2%rnd_oris(5.)
            write(logfhandle,*) '********'
            do i=1,100
                call os2%print_(i)
            end do
        endif
        write(logfhandle,'(a)') '**info(simple_oris_unit_test, part2): testing assignment'
        os = oris(2, is_ptcl=.false.)
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
        if( .not. passed ) THROW_HARD('assignment test failed!')
        write(logfhandle,'(a)') '**info(simple_oris_unit_test, part2): testing i/o'
        passed = .false.
        os  = oris(100, is_ptcl=.false.)
        os2 = oris(100, is_ptcl=.false.)
        call os%rnd_oris(5.)
        call os%write('test_oris_rndoris.txt')
        call os2%read('test_oris_rndoris.txt')
        call os2%write('test_oris_rndoris_copy.txt')
        corr = os%corr_oris(os2)
        if( corr > 0.99 ) passed = .true.
        if( .not. passed ) THROW_HARD('read/write failed')
        passed = .false.
        call os%rnd_states(5)
        call os%write('test_oris_rndoris_rndstates.txt')
        if( os%corr_oris(os2) > 0.99 ) passed = .true.
        if( .not. passed ) THROW_HARD('statedoc read/write failed!')
        write(logfhandle,'(a)') '**info(simple_oris_unit_test, part3): testing calculators'
        passed = .false.
        call os%rnd_lps()
        call os%write('test_oris_rndoris_rndstates_rndlps.txt')
        call os%spiral
        call os%write('test_oris_rndoris_rndstates_rndlps_spiral.txt')
        call os%rnd_corrs()
        order = os%order()
        if( doprint )then
            do i=1,100
                call os%print_(order(i))
            end do
            write(logfhandle,*) 'median:', os%median('lp')
        endif
        write(logfhandle,'(a)') '**info(simple_oris_unit_test, part4): testing destructor'
        call os%kill
        call os2%kill
        ! test find_angres
        os = oris(1000, is_ptcl=.false.)
        call os%spiral
        write(logfhandle,*) 'angres:      ', os%find_angres()
        write(logfhandle,'(a)') 'SIMPLE_ORIS_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
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
            if(alloc_stat.ne.0)call allocchk('In: kill, module: simple_oris o')
            self%n = 0
        endif
    end subroutine kill

end module simple_oris
