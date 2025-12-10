! an agglomeration of orientations
module simple_oris
!$ use omp_lib
!$ use omp_lib_kinds
use json_kinds
use json_module
use simple_defs
use simple_defs_ori
use simple_fileio
use simple_fileio
use simple_is_check_assert
use simple_math
use simple_math_ft
use simple_ori
use simple_ran_tabu
use simple_rnd
use simple_srch_sort_loc
use simple_stat
use simple_string
use simple_string_utils
use simple_syslib
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
    generic            :: new => new_1, new_2
    procedure, private :: new_1, new_2
    procedure          :: reallocate
    procedure, private :: extract_subset_1, extract_subset_2
    generic            :: extract_subset => extract_subset_1, extract_subset_2
    ! GETTERS
    procedure          :: exists
    procedure          :: e1get
    procedure          :: e2get
    procedure          :: e3get
    procedure          :: get_euler
    procedure          :: get_noris
    procedure          :: get_ori
    procedure          :: get
    procedure          :: get_int
    procedure          :: get_str
    procedure          :: get_static
    procedure, private :: getter_1, getter_2, getter_3
    generic            :: getter => getter_1, getter_2, getter_3
    procedure          :: get_all, get_all_asint
    procedure          :: gen_ptcl_mask
    procedure          :: get_all_sampled
    procedure          :: get_all_rmats
    procedure          :: get_mat
    procedure          :: get_normal
    procedure          :: get_2Dshift
    procedure          :: get_dfx, get_dfy
    procedure          :: get_state
    procedure          :: get_class
    procedure          :: get_proj
    procedure          :: get_label_inds
    procedure          :: get_eo
    procedure          :: get_fromp, get_top
    procedure          :: get_sampled
    procedure          :: get_updatecnt
    procedure          :: get_tseries_neighs
    procedure, private :: isthere_1, isthere_2
    generic            :: isthere => isthere_1, isthere_2
    procedure          :: ischar
    procedure          :: is_particle
    procedure          :: max_ori_strlen_trim
    procedure          :: get_n
    procedure, private :: get_pop_1, get_pop_2
    generic            :: get_pop => get_pop_1, get_pop_2
    procedure          :: get_pops
    procedure          :: get_pinds
    procedure          :: gen_mask
    procedure          :: mask_from_state
    procedure, private :: get_all_normals
    procedure          :: states_exist
    procedure          :: projs_exist
    procedure          :: get_arr
    procedure, private :: calc_sum
    procedure          :: get_sum
    procedure          :: get_avg
    procedure          :: included
    procedure          :: get_nevenodd
    procedure          :: get_neven
    procedure          :: get_nodd
    procedure          :: print
    procedure          :: print_matrices
    procedure          :: select_particles_set
    procedure          :: sample4rec
    procedure          :: sample4update_all
    procedure          :: sample4update_rnd
    procedure          :: sample4update_cnt
    procedure          :: sample4update_class
    procedure          :: sample4update_reprod
    procedure          :: sample4update_updated
    procedure          :: sample4update_fillin
    procedure          :: get_update_frac
    procedure          :: get_class_sample_stats
    procedure          :: get_proj_sample_stats
    procedure, private :: sample_balanced_1, sample_balanced_2
    generic            :: sample_balanced => sample_balanced_1, sample_balanced_2
    procedure          :: sample_balanced_inv
    procedure          :: sample_balanced_parts
    procedure          :: sample_ranked_parts
    procedure          :: balance_ptcls_within_cls
    procedure, private :: get_sample_ind
    procedure, private :: incr_sampled_updatecnt
    procedure          :: is_first_update
    procedure          :: set_nonzero_updatecnt
    procedure          :: set_updatecnt
    procedure          :: clean_entry
    procedure          :: has_been_sampled
    procedure          :: has_been_searched
    procedure          :: any_state_zero
    procedure          :: count_state_gt_zero
    procedure          :: ori2str
    procedure          :: ori2json
    procedure          :: ori2prec
    procedure          :: prec2ori
    procedure          :: get_ctfvars
    ! SETTERS
    procedure, private :: append_1, append_2
    generic            :: append => append_1, append_2
    procedure, private :: copy_1, copy_2
    generic            :: copy => copy_1, copy_2
    procedure          :: reject
    procedure, private :: delete_entry_1, delete_entry_2
    generic            :: delete_entry => delete_entry_1, delete_entry_2
    procedure, private :: delete_2Dclustering_1, delete_2Dclustering_2
    generic            :: delete_2Dclustering => delete_2Dclustering_1, delete_2Dclustering_2
    procedure          :: delete_3Dalignment
    procedure          :: delete
    procedure          :: transfer_2Dshifts
    procedure, private :: transfer_2Dparams_1, transfer_2Dparams_2
    generic            :: transfer_2Dparams => transfer_2Dparams_1, transfer_2Dparams_2
    procedure, private :: transfer_3Dparams_1, transfer_3Dparams_2
    generic            :: transfer_3Dparams => transfer_3Dparams_1, transfer_3Dparams_2
    procedure          :: transfer_class_assignment
    procedure          :: set_euler
    procedure          :: set_shift
    procedure          :: set_state
    procedure          :: set_class
    procedure          :: set_stkind
    procedure          :: set_ogid
    procedure          :: e1set
    procedure          :: e2set
    procedure          :: e3set
    procedure, private :: set_1, set_2, set_3, set_4, set_5
    generic            :: set => set_1, set_2, set_3, set_4, set_5
    procedure          :: set_dfx, set_dfy
    procedure          :: set_ori
    procedure          :: transfer_ori
    procedure, private :: set_all_1, set_all_2, set_all_3
    generic            :: set_all => set_all_1, set_all_2, set_all_3
    procedure, private :: set_all2single_1, set_all2single_2, set_all2single_3
    generic            :: set_all2single => set_all2single_1, set_all2single_2, set_all2single_3
    procedure, private :: set_field2single_1, set_field2single_2, set_field2single_3
    generic            :: set_field2single => set_field2single_1, set_field2single_2, set_field2single_3
    procedure          :: set_projs
    procedure          :: remap_projs
    procedure          :: proj2class
    procedure          :: e3swapsgn
    procedure          :: swape1e3
    procedure          :: zero
    procedure          :: zero_projs
    procedure          :: zero_shifts
    procedure          :: zero_inpl
    procedure          :: mul_shifts
    procedure          :: rnd_oris
    procedure          :: gau_rnd_shifts
    procedure          :: rnd_ori
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
    procedure          :: partition_eo
    procedure          :: str2ori
    procedure          :: str2ori_ctfparams_state_eo
    procedure          :: set_ctfvars
    procedure          :: gen_balanced_partitions
    ! I/O
    procedure          :: read
    procedure          :: read_ctfparams_state_eo
    procedure, private :: write_1, write_2
    generic            :: write => write_1, write_2
    procedure          :: write2bild
    ! CALCULATORS
    procedure          :: compress
    procedure          :: split_state
    procedure          :: split_class
    procedure          :: expand_classes
    procedure          :: remap_cls
    procedure          :: merge_classes
    procedure          :: round_shifts
    procedure          :: introd_alig_err
    procedure          :: introd_ctf_err
    procedure, private :: rot_1, rot_2
    generic            :: rot => rot_1, rot_2
    procedure, private :: rot_transp_1, rot_transp_2
    generic            :: rot_transp => rot_transp_1, rot_transp_2
    procedure, private :: median_1
    generic            :: median => median_1
    procedure, private :: stats_1, stats_2, stats_3
    generic            :: stats => stats_1, stats_2, stats_3
    procedure          :: minmax
    procedure, private :: spiral_1, spiral_2
    generic            :: spiral => spiral_1, spiral_2
    procedure          :: order
    procedure          :: order_cls
    procedure          :: calc_hard_weights
    procedure          :: calc_soft_weights, calc_cavg_soft_weights
    procedure          :: calc_hard_weights2D
    procedure          :: calc_soft_weights2D
    procedure          :: find_best_classes
    procedure, private :: euldist_1, euldist_2
    generic            :: euldist => euldist_1, euldist_2
    procedure          :: find_closest_proj
    procedure          :: discretize
    procedure, private :: nearest_proj_neighbors_1, nearest_proj_neighbors_2, nearest_proj_neighbors_3
    generic            :: nearest_proj_neighbors => nearest_proj_neighbors_1, nearest_proj_neighbors_2, nearest_proj_neighbors_3
    procedure          :: extract_subspace
    procedure          :: detect_peaks
    procedure          :: min_euldist
    procedure          :: find_angres
    procedure          :: extremal_bound
    procedure          :: set_extremal_vars
    procedure, private :: map3dshift22d_1, map3dshift22d_2
    generic            :: map3dshift22d => map3dshift22d_1, map3dshift22d_2
    procedure          :: calc_avg_offset2D, calc_avg_offset3D
    procedure          :: mirror2d
    procedure          :: mirror3d
    procedure          :: add_shift2class
    procedure          :: corr_oris
    procedure, private :: diststat_1, diststat_2
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

    subroutine new_1( self, n, is_ptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        logical,     intent(in)    :: is_ptcl
        integer :: i
        call self%kill
        self%n = n
        allocate( self%o(self%n) )
        do i=1,n
            call self%o(i)%new_ori(is_ptcl)
        end do
    end subroutine new_1

    subroutine new_2( self, o_arr )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_arr(:)
        integer :: i
        call self%kill
        self%n = size(o_arr)
        allocate( self%o(self%n) )
        do i=1,self%n
            call self%o(i)%copy(o_arr(i))
        end do
    end subroutine new_2

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
        call tmp%kill
    end subroutine reallocate

    function extract_subset_1( self, from, to ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: from, to
        type(oris) :: self_sub
        integer    :: n, cnt, i
        n  = to - from + 1
        call self_sub%new(n, self%is_particle())
        cnt = 0
        do i = from, to
            cnt = cnt + 1
            call self_sub%o(cnt)%copy(self%o(i))
        end do
    end function extract_subset_1

    function extract_subset_2( self, inds ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: inds(:)
        type(oris) :: self_sub
        integer    :: n, i
        n  = size(inds)
        call self_sub%new(n, self%is_particle())
        do i = 1,n
            call self_sub%o(i)%copy(self%o(inds(i)))
        end do
    end function extract_subset_2

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
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: i
        type(ori),   intent(inout) :: o
        if( self%n == 0 ) THROW_HARD('oris object does not exist; get_ori')
        if( i > self%n .or. i < 1 )then
            write(logfhandle,*) 'trying to get ori: ', i, ' among: ', self%n, ' oris'
            THROW_HARD('i out of range; get_ori')
        endif
        o = self%o(i)
    end subroutine get_ori

    pure real function get( self, i, key )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        get = self%o(i)%get(key)
    end function get

    pure integer function get_int( self, i, key )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        call self%o(i)%getter(key, get_int)
    end function get_int

    function get_str( self, i, key ) result( val )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        type(string) :: val
        call self%getter_1(i, key, val)
    end function get_str

    pure subroutine get_static( self, i, key, val )
        class(oris),      intent(in)  :: self
        integer,          intent(in)  :: i
        character(len=*), intent(in)  :: key
        character(len=*), intent(out) :: val
        call self%o(i)%get_static(key, val)
    end subroutine get_static

    subroutine getter_1( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        type(string),     intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_1

    subroutine getter_2( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_2

    subroutine getter_3( self, i, key, ival )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(inout) :: ival
        call self%o(i)%getter(key, ival)
    end subroutine getter_3

    !>  \brief  is for getting an array of 'key' values
    function get_all( self, key, fromto, nonzero ) result( arr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: fromto(2)
        logical, optional, intent(in) :: nonzero
        real,    allocatable :: arr(:)
        integer :: ffromto(2)
        logical :: nnonzero
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        nnonzero = .false.
        if( present(nonzero) ) nnonzero = nonzero
        allocate( arr(ffromto(1):ffromto(2)), source=self%o(ffromto(1):ffromto(2))%get(key) )
        if( nnonzero )then
            arr = pack(arr, mask=self%o(ffromto(1):ffromto(2))%get_state()>0)
        endif
    end function get_all

    !>  \brief  is for getting an array of 'key' values cast to integer
    function get_all_asint( self, key, fromto ) result( iarr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: fromto(2)
        integer, allocatable :: iarr(:)
        integer :: ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        allocate(iarr(ffromto(1):ffromto(2)), source=nint(self%o(ffromto(1):ffromto(2))%get(key)))
    end function get_all_asint

    function gen_ptcl_mask( self, key, ival, frac_best ) result( mask )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer,           intent(in) :: ival
        real,    optional, intent(in) :: frac_best
        logical, allocatable :: mask(:)
        integer, allocatable :: ivals(:)
        real,    allocatable :: corrs(:), tmp(:)
        real    :: corr_t
        integer :: ncorrs, nsample
        allocate(mask(self%n))
        allocate(ivals(self%n), source=self%o(1:self%n)%get_int(key))
        if( present(frac_best) )then
            allocate(corrs(self%n), source=self%o(self%n)%get('corr'))
            tmp     = pack(corrs, mask=ivals == ival)
            ncorrs  = size(tmp)
            nsample = ceiling(frac_best * real(ncorrs))
            call hpsort(tmp)
            corr_t  = tmp(ncorrs - nsample)
            deallocate(tmp)
            where(ivals == ival .and. corrs >= corr_t )
                mask = .true.
            elsewhere
                mask = .false.
            endwhere
        else
            where(ivals == ival )
                mask = .true.
            elsewhere
                mask = .false.
            endwhere
        endif
    end function gen_ptcl_mask

    !>  \brief  is for getting an array of 'key' values
    function get_all_sampled( self, key, state, lowerbound ) result( arr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: state
        real,    optional, intent(in) :: lowerbound
        real,    allocatable :: arr(:), sampled(:)
        integer, allocatable :: states(:)
        integer :: i
        real    :: lb
        allocate(arr(self%n), sampled(self%n), source=0.)
        if( present(state) ) states = nint(self%get_all('state'))
        do i=1,self%n
            arr(i)     = self%o(i)%get(key)
            sampled(i) = self%o(i)%get('sampled')
        enddo
        if( present(lowerbound) )then
            lb = lowerbound
        else
            lb = maxval(sampled) - 0.5
        endif
        if( present(state) )then
            arr = pack(arr, mask=(sampled > lb .and. states == state))
        else
            arr = pack(arr, mask=sampled > lb)
        endif
    end function get_all_sampled

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

    pure function get_2Dshift( self, i )  result( shvec )
        class(oris), intent(in) :: self
        integer,     intent(in)    :: i
        real :: shvec(2)
        shvec = self%o(i)%get_2Dshift()
    end function get_2Dshift

    pure real function get_dfx( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_dfx = self%o(i)%get_dfx()
    end function get_dfx

    pure real function get_dfy( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_dfy = self%o(i)%get_dfy()
    end function get_dfy

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

    pure integer function get_proj( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_proj = self%o(i)%get_proj()
    end function get_proj

    function get_label_inds( self, label ) result( inds )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: label
        integer :: i, icls, icls_max
        integer, allocatable :: inds(:), inds_all(:)
        logical, allocatable :: isthere(:)
        icls_max = self%get_n(label)
        inds_all = (/(i,i=1,icls_max)/)
        allocate(isthere(icls_max), source=.false.)
        do i = 1, self%n
            icls = self%get_class(i)
            isthere(icls) = .true.
        end do
        inds = pack(inds_all, mask=isthere)
    end function get_label_inds

    pure integer function get_eo( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_eo = self%o(i)%get_eo()
    end function get_eo

    pure integer function get_fromp( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_fromp = self%o(i)%get_fromp()
    end function get_fromp

    pure integer function get_top( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_top = self%o(i)%get_top()
    end function get_top

     pure integer function get_sampled( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_sampled = self%o(i)%get_sampled()
    end function get_sampled

    pure integer function get_updatecnt( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_updatecnt = self%o(i)%get_updatecnt()
    end function get_updatecnt

    ! assumes project has been pruned to remove state=0 particles
    subroutine get_tseries_neighs( self, nsz, ptcls2neigh )
       class(oris),          intent(in)    :: self
       integer,              intent(in)    :: nsz
       integer, allocatable, intent(inout) :: ptcls2neigh(:,:)
       integer :: i, j, cls1, cls2
       if( allocated(ptcls2neigh) ) deallocate(ptcls2neigh)
       allocate(ptcls2neigh(self%n,2), source=0)
       do i=1,self%n-nsz
           cls1 = self%o(i)%get_class()
           cls2 = self%o(i+1)%get_class()
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

    function ischar( self, i, key ) result( is )
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
            ival = self%o(i)%get_int(label)
            if( ival > n ) n = ival
        end do
    end function get_n

    !>  \brief  is for checking label population
    function get_pop_1( self, ind, label, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: ind
        character(len=*),  intent(in) :: label
        integer, optional, intent(in) :: eo
        integer :: pop, i
        logical :: consider_eo
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        pop = 0
        if( consider_eo )then
            !$omp parallel do private(i) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero()  ) cycle
                if( self%o(i)%get_eo() /= eo ) cycle
                if( self%o(i)%get_int(label) == ind )  pop = pop + 1
            end do
            !$omp end parallel do
        else
            !$omp parallel do private(i) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero() ) cycle
                if( self%o(i)%get_int(label) == ind )  pop = pop + 1
            end do
            !$omp end parallel do
        endif
    end function get_pop_1

    function get_pop_2( self, inds, labels, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: inds(:)
        character(len=*),  intent(in) :: labels(:)
        integer, optional, intent(in) :: eo
        integer :: pop, i, j, n
        logical :: consider_eo, l_okay
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        n   = size(inds)
        pop = 0
        if( consider_eo )then
            !$omp parallel do private(i,j,l_okay) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero()  ) cycle
                if( self%o(i)%get_eo() /= eo ) cycle
                l_okay = .true.
                do j=1,n
                    if( self%o(i)%get_int(labels(j)) /= inds(j) )then
                        l_okay = .false.
                        exit
                    endif
                enddo
                if( l_okay )  pop = pop + 1
            end do
            !$omp end parallel do
        else
            !$omp parallel do private(i,j,l_okay) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero() ) cycle
                l_okay = .true.
                do j=1,n
                    if( self%o(i)%get_int(labels(j)) /= inds(j) )then
                        l_okay = .false.
                        exit
                    endif
                enddo
                if( l_okay )  pop = pop + 1
            end do
            !$omp end parallel do
        endif
    end function get_pop_2

    !>  \brief  is for getting all rotation matrices
    function get_all_rmats( self ) result( mat )
        class(oris), intent(in) :: self
        real, allocatable       :: mat(:,:,:)
        integer :: i,n
        n = self%n
        if(allocated(mat))deallocate(mat)
        allocate(mat(n,3,3),source=0.)
        do i=1,self%n
            mat(i,:,:) = self%o(i)%get_mat()
        end do
    end function get_all_rmats

    subroutine get_pops( self, pops, label, maxn, weight )
        class(oris),          intent(in)    :: self
        integer, allocatable, intent(out)   :: pops(:)
        character(len=*),     intent(in)    :: label
        integer, optional,    intent(in)    :: maxn   ! max label, for the case where the last class/state is missing
        logical, optional,    intent(in)    :: weight ! whether to consider weights
        integer :: i, val, n
        logical :: consider_w
        consider_w = .false.
        if( present(weight) ) consider_w = weight
        n = self%get_n(label)
        if( present(maxn) ) n = max(n, maxn)
        if(allocated(pops)) deallocate(pops)
        allocate(pops(n),source=0)
        !$omp parallel do private(i,val) default(shared) proc_bind(close) reduction(+:pops)
        do i = 1,self%n
            if( self%o(i)%isstatezero() ) cycle
            if( consider_w )then
                if( self%o(i)%get('w') < 1.e-6 ) cycle
            endif
            val = self%o(i)%get_int(label)
            if( val > 0 ) pops(val) = pops(val) + 1
        end do
        !$omp end parallel do
    end subroutine get_pops

    function get_all_normals( self ) result( normals )
        class(oris), intent(inout) :: self
        real, allocatable :: normals(:,:)
        integer :: i
        allocate(normals(self%n,3))
        do i=1,self%n
            normals(i,:) = self%o(i)%get_normal()
        end do
    end function get_all_normals

    !>  \brief  returns a logical array of state existence
    function states_exist( self, nstates, thres ) result( state_exists )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: nstates
        integer, optional, intent(in)    :: thres
        integer :: i, min_pop
        logical :: state_exists(nstates)
        min_pop = 0
        if( present(thres) ) min_pop = thres
        do i=1,nstates
            state_exists(i) = (self%get_pop(i, 'state') > min_pop)
        end do
    end function states_exist

    !>  \brief  returns a logical array of state existence
    function projs_exist( self, nstates, nprojs, thres ) result( proj_exists )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: nstates
        integer,           intent(in)    :: nprojs
        integer, optional, intent(in)    :: thres
        integer :: i, j, min_pop
        logical :: proj_exists(nprojs,nstates)
        min_pop = 0
        if( present(thres) ) min_pop = thres
        do i=1,nstates
            do j=1,nprojs
                proj_exists(j,i) = (self%get_pop([i,j], ['state', 'proj ']) > min_pop)
            enddo
        end do
    end function projs_exist

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
        allocate(states(n))
        rt = ran_tabu(n)
        call rt%balanced(2, states)
        do iptcl=1,n
            if( states(iptcl) == 1 )then
                ! do nothing, leave this state as is
            else
                call self%o(ptcls_in_which(iptcl))%set_state(nstates+1)
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
        allocate(members(n))
        rt = ran_tabu(n)
        call rt%balanced(2, members)
        do iptcl=1,n
            if( members(iptcl) == 1 )then
                ! do nothing, leave this state as is
            else
                call self%o(ptcls_in_which(iptcl))%set_class(nmembers+1)
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
        allocate(pops(ncls_target))
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

    !>  \brief  for remapping clusters
    subroutine remap_cls( self )
        class(oris), intent(inout) :: self
        integer :: ncls, clsind_remap, pop, icls, iptcl, old_cls
        integer , allocatable :: clspops(:)
        ncls = self%get_n('class')
        allocate(clspops(ncls))
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
                        old_cls = self%o(iptcl)%get_class()
                        if( old_cls == icls ) call self%o(iptcl)%set_class(clsind_remap)
                    end do
                else
                    do iptcl=1,self%n
                        old_cls = self%o(iptcl)%get_class()
                        if( old_cls == icls )then
                            call self%o(iptcl)%set_class(0)
                            call self%o(iptcl)%set_state(0)
                        endif
                    end do
                endif
            end do
            do iptcl=1,self%n
                old_cls = self%o(iptcl)%get_class()
                if( old_cls /= 0 ) call self%o(iptcl)%set_class(old_cls-ncls)
            end do
        endif
        deallocate(clspops)
    end subroutine remap_cls

    !>  \brief  is for getting an allocatable array with ptcl indices of the label 'label'
    subroutine get_pinds( self, ind, label, indices, l_shuffle, l_require_updated )
        class(oris),          intent(in)  :: self
        character(len=*),     intent(in)  :: label
        integer,              intent(in)  :: ind
        integer, allocatable, intent(out) :: indices(:)
        logical, optional,    intent(in)  :: l_shuffle, l_require_updated
        type(ran_tabu)       :: rt
        logical, allocatable :: mask(:)
        integer :: pop, i
        logical :: ll_shuffle, ll_require_updated
        ll_shuffle = .false.
        if( present(l_shuffle) ) ll_shuffle = l_shuffle
        ll_require_updated = .false.
        if( present(l_require_updated) ) ll_require_updated = l_require_updated
        if( allocated(indices) )deallocate(indices)
        allocate(indices(self%n),mask(self%n))
        !$omp parallel do private(i) default(shared) proc_bind(close)
        do i = 1,self%n
            if( self%o(i)%isstatezero() )then
                mask(i) = .false.
            else
                if( ll_require_updated )then
                    if( self%o(i)%get_int('updatecnt') == 0 )then
                        mask(i) = .false.
                    else
                        mask(i) = self%o(i)%get_int(label) == ind
                    endif
                else
                    mask(i) = self%o(i)%get_int(label) == ind
                endif
            endif
            if( mask(i) ) indices(i) = i
        end do
        !$omp end parallel do
        pop = count(mask)
        if( pop > 0 )then
            indices = pack(indices, mask=mask)
            deallocate(mask)
            if( ll_shuffle )then
                rt  = ran_tabu(pop)
                call rt%shuffle(indices)
                call rt%kill
            endif
        else
            deallocate(indices,mask)
        endif
    end subroutine get_pinds

    !>  \brief  generate a mask with the oris with mystate == state/ind == get(label)
    subroutine gen_mask( self, state, ind, label, l_mask, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state, ind
        character(len=*),     intent(in)    :: label
        logical, allocatable, intent(out)   :: l_mask(:)
        integer, optional,    intent(in)    :: fromto(2)
        integer :: i, mystate, myval, ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        if( allocated(l_mask) ) deallocate(l_mask)
        allocate(l_mask(ffromto(1):ffromto(2)))
        l_mask = .false.
        do i=ffromto(1),ffromto(2)
            mystate = self%o(i)%get_state()
            myval   = self%o(i)%get_int(label)
            if( mystate == state .and. myval == ind ) l_mask(i) = .true.
        end do
    end subroutine gen_mask

    !>  \brief  generate a mask for a single value of label state
    subroutine mask_from_state( self, state, l_mask, pinds, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state
        logical, allocatable, intent(inout) :: l_mask(:)
        integer, allocatable, intent(inout) :: pinds(:)
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
            allocate( vals(pop) )
            cnt = 0
            do i=1,self%n
                val = self%get(i, which)
                if( class_present )then
                    mystate = self%o(i)%get_state()
                    if( mystate > 0 )then
                        clsnr   = self%o(i)%get_class()
                        if( clsnr == class )then
                            cnt = cnt+1
                            vals(cnt) = val
                        endif
                    endif
                else if( state_present )then
                    mystate = self%o(i)%get_state()
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
            mystate = self%o(i)%get_state()
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
                    clsnr = self%o(i)%get_class()
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
    function included( self ) result( incl )
        class(oris),       intent(inout) :: self
        logical, allocatable :: incl(:)
        integer :: i
        if(.not.allocated(incl))allocate(incl(self%n))
        incl = .false.
        do i=1,self%n
            incl(i) = self%o(i)%get_state() > 0
        end do
    end function included

    integer function get_neven( self )
        class(oris), intent(inout) :: self
        integer :: i
        get_neven = 0
        do i=1,self%n
            if(self%o(i)%isthere('eo'))then
                if( self%o(i)%get_eo()==0) get_neven = get_neven + 1
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

    subroutine print( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        call self%o(i)%print_ori()
    end subroutine print

    subroutine print_matrices( self )
        class(oris), intent(inout) :: self
        integer :: i
        write(logfhandle,*) 'ORDER OF ROTATION MATRIX ELEMENTS: (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)'
        do i=1,self%n
            call self%o(i)%print_mat()
        end do
    end subroutine print_matrices

    subroutine select_particles_set( self, nptcls, inds )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nptcls
        integer, allocatable, intent(inout) :: inds(:)
        type(ran_tabu) :: rt
        integer        :: i,n
        if( allocated(inds) ) deallocate(inds)
        inds = (/(i,i=1,self%n)/)
        inds = pack(inds, mask=self%get_all('state') > 0.5)
        n    = size(inds)
        if( n < nptcls )then
            ! if less than desired particles select all
        else
            rt   = ran_tabu(n)
            call rt%shuffle(inds)
            call rt%kill
            inds = inds(1:nptcls)
            call hpsort(inds) !!
        endif
    end subroutine select_particles_set

    subroutine sample4rec( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        integer, allocatable :: states(:), updatecnts(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), updatecnts(nptcls), inds(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            states(cnt)     = self%o(i)%get_state()
            updatecnts(cnt) = self%o(i)%get_int('updatecnt')
            inds(cnt)       = i
        end do
        if( any(updatecnts > 0) )then
            nsamples = count(states > 0 .and. updatecnts > 0)
            inds     = pack(inds, mask=states > 0 .and. updatecnts > 0)
        else
            nsamples = count(states > 0)
            inds     = pack(inds, mask=states > 0)
        endif
    end subroutine sample4rec

    subroutine sample4update_all( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: states(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), inds(nptcls), source=0)
        cnt      = 0
        nsamples = 0
        do i = fromto(1), fromto(2)
            cnt         = cnt + 1
            states(cnt) = self%o(i)%get_state()
            inds(cnt)   = i
            if( states(cnt) > 0 ) nsamples = nsamples + 1
        end do
        inds = pack(inds, mask=states > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_all

    subroutine sample4update_rnd( self, fromto, update_frac, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        type(ran_tabu) :: rt
        integer, allocatable :: states(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), inds(nptcls), source=0)
        cnt    = 0
        nptcls = 0
        do i = fromto(1), fromto(2)
            cnt         = cnt + 1
            states(cnt) = self%o(i)%get_state()
            inds(cnt)   = i
            if( states(cnt) > 0 ) nptcls = nptcls + 1
        end do
        inds     = pack(inds,       mask=states > 0)
        nsamples = min(nptcls, nint(update_frac * real(nptcls)))
        rt       = ran_tabu(nptcls)
        call rt%shuffle(inds)
        call rt%kill
        inds = inds(1:nsamples)
        call hpsort(inds) ! indices in increasing order
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_rnd

    subroutine sample4update_cnt( self, fromto, update_frac, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        type(ran_tabu)       :: rt
        integer, allocatable :: states(:), updatecnts(:)
        integer              :: i, cnt, nptcls, ucnt, ucnt_min, ucnt_max, ucnt_lim
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), updatecnts(nptcls), inds(nptcls), source=0)
        cnt    = 0
        nptcls = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            states(cnt)     = self%o(i)%get_state()
            updatecnts(cnt) = self%o(i)%get_int('updatecnt')
            inds(cnt)       = i
            if( states(cnt) > 0 ) nptcls = nptcls + 1
        end do
        inds        = pack(inds,       mask=states > 0)
        updatecnts  = pack(updatecnts, mask=states > 0)
        ucnt_min    = minval(updatecnts)
        ucnt_max    = maxval(updatecnts)
        nsamples    = min(nptcls, nint(update_frac * real(nptcls)))
        ucnt_lim    = 0
        if( ucnt_max > ucnt_min )then
            do ucnt = ucnt_min,ucnt_max
                if( count(updatecnts < ucnt) >= nsamples )then
                    ucnt_lim = ucnt
                    exit
                endif
            end do
        endif
        if( ucnt_lim > 0 )then
            inds   = pack(inds, mask=updatecnts < ucnt)
            nptcls = size(inds)
        endif
        rt = ran_tabu(nptcls)
        call rt%shuffle(inds)
        call rt%kill
        inds = inds(1:nsamples)
        call hpsort(inds) ! indices in increasing order
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_cnt

    subroutine sample4update_class( self, clssmp, fromto, update_frac, nsamples, inds, incr_sampled, l_greedy, frac_best )
        class(oris),          intent(inout) :: self
        type(class_sample),   intent(inout) :: clssmp(:) ! data structure for balanced sampling
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled, l_greedy
        real,    optional,    intent(in)    :: frac_best
        integer, allocatable :: states(:)
        real,    allocatable :: rstates(:)
        integer :: i, cnt, nptcls, nsamples_class, states_bal(self%n)
        ! balanced sampling is global
        rstates        = self%get_all('state')
        nsamples_class = nint(update_frac * real(count(rstates > 0.5)))
        deallocate(rstates)
        ! class-biased selection
        if( present(frac_best) )then
            call self%sample_balanced(clssmp, nsamples_class, frac_best,  states_bal) ! stochastic sampling from frac_best fraction
        else
            call self%sample_balanced(clssmp, nsamples_class, l_greedy, states_bal)
        endif
        ! now, we deal with the partition
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), inds(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt          = cnt + 1
            states(cnt)  = states_bal(i)
            inds(cnt)    = i
        end do
        nsamples = count(states > 0)
        inds     = pack(inds, mask=states > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_class

    subroutine sample4update_reprod( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        integer, allocatable :: sampled(:)
        integer :: i, cnt, nptcls, sample_ind
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), sampled(nptcls), source=0)
        cnt        = 0
        sample_ind = self%get_sample_ind(.false.)
        do i = fromto(1), fromto(2)
            cnt          = cnt + 1
            inds(cnt)    = i
            sampled(cnt) = self%o(i)%get_sampled()
        end do
        if( sample_ind == 0 ) THROW_HARD('requires previous sampling')
        nsamples = count(sampled == sample_ind)
        inds     = pack(inds, mask=sampled == sample_ind)
    end subroutine sample4update_reprod

    subroutine sample4update_updated( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: updatecnts(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), updatecnts(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            inds(cnt)       = i
            updatecnts(cnt) = self%o(i)%get_updatecnt()
        end do
        if( .not. any(updatecnts > 0) ) THROW_HARD('requires previous update')
        nsamples = count(updatecnts > 0)
        inds     = pack(inds, mask=updatecnts > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_updated

    subroutine sample4update_fillin( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: updatecnts(:), states(:)
        integer :: i, cnt, nptcls, minucnt
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), states(nptcls), updatecnts(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            inds(cnt)       = i
            states(cnt)     = self%o(i)%get_state()
            updatecnts(cnt) = self%o(i)%get_updatecnt()
        end do
        updatecnts = pack(updatecnts, mask=states > 0)
        inds       = pack(inds,       mask=states > 0)
        minucnt    = minval(updatecnts)
        nsamples   = count(updatecnts == minucnt)
        inds       = pack(inds, mask=updatecnts == minucnt)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_fillin

    function get_update_frac( self ) result( update_frac )
        class(oris), intent(inout) :: self
        integer :: updatecnts(self%n), sampled(self%n), states(self%n), updatecnt_max, sampled_max, i
        real    :: update_frac
        sampled_max   = 0
        updatecnt_max = 0
        do i = 1, self%n
            updatecnts(i) = self%o(i)%get_updatecnt()
            sampled(i)    = self%o(i)%get_sampled()
            states(i)     = self%o(i)%get_state()
            sampled_max   = max(sampled_max,sampled(i))
            updatecnt_max = max(updatecnt_max,updatecnts(i))
        end do
        if( sampled_max   == 0 ) THROW_HARD('requires previous sampling')
        if( updatecnt_max == 0 ) THROW_HARD('requires previous update')
        update_frac = real(count(sampled == sampled_max .and. states > 0)) / real(count(updatecnts > 0 .and. states > 0))
    end function get_update_frac

    subroutine get_class_sample_stats( self, clsinds, clssmp, label )
        class(oris),                     intent(inout) :: self
        integer,                         intent(in)    :: clsinds(:) ! class indices to sample from
        type(class_sample), allocatable, intent(inout) :: clssmp(:)  ! data structure for balanced samplign
        character(len=*),      optional, intent(in)    :: label
        character(len=:), allocatable :: flag
        integer :: n, i, j, nc
        if( present(label) )then
            flag = trim(label)
        else
            flag = 'class'
        endif
        ! init data structure
        n = size(clsinds)
        if( allocated(clssmp) )then
            nc = size(clssmp)
            do i = 1, nc
                if( allocated(clssmp(i)%ccs)   ) deallocate(clssmp(i)%ccs)
                if( allocated(clssmp(i)%pinds) ) deallocate(clssmp(i)%pinds)
            end do
            deallocate(clssmp)
        endif
        allocate(clssmp(n))
        ! fetch information necessary for balanced sampling
        do i = 1, n
            call self%get_pinds(clsinds(i), flag, clssmp(i)%pinds)
            if( allocated(clssmp(i)%pinds) )then
                clssmp(i)%clsind = clsinds(i)
                clssmp(i)%pop    = size(clssmp(i)%pinds)
                allocate(clssmp(i)%ccs(clssmp(i)%pop), source=0.)
                do j = 1, clssmp(i)%pop
                    clssmp(i)%ccs(j) = self%o(clssmp(i)%pinds(j))%get('corr')
                end do
                call hpsort(clssmp(i)%ccs, clssmp(i)%pinds)
                call reverse(clssmp(i)%ccs)   ! best first
                call reverse(clssmp(i)%pinds) ! best first
            endif
        end do
    end subroutine get_class_sample_stats

    ! for balanced sampling in refine3D_auto, after abinitio3D
    ! use the abinitio3D alignment to create a class assignment
    ! use the re-projection correlations for ranking
    ! while leaving the original instance untouched
    subroutine get_proj_sample_stats( self, eulspace, clssmp )
        class(oris),                     intent(inout) :: self
        class(oris),                     intent(in)    :: eulspace
        type(class_sample), allocatable, intent(inout) :: clssmp(:)  ! data structure for balanced samplign
        integer, allocatable :: tmpinds(:), clsinds(:), clspops(:)
        integer              :: ncls, icls
        type(oris)           :: self_copy
        ncls = eulspace%get_noris()
        call self_copy%copy(self)
        call self_copy%set_projs(eulspace)
        call self_copy%proj2class
        allocate(clspops(ncls))
        do icls=1,ncls
            clspops(icls) = self_copy%get_pop(icls, 'class')
        end do
        tmpinds = (/(icls,icls=1,ncls)/)
        clsinds = pack(tmpinds, mask= clspops > 0)
        call self_copy%get_class_sample_stats(clsinds, clssmp)
        call self_copy%kill
    end subroutine get_proj_sample_stats

    subroutine sample_balanced_1( self, clssmp, nptcls, l_greedy, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)  ! data structure for balanced sampling
        integer,            intent(in)    :: nptcls     ! # particles to sample in total
        logical,            intent(in)    :: l_greedy   ! greedy or stochastic sampling
        integer,            intent(inout) :: states(self%n)
        integer,            allocatable   :: pinds_left(:)
        type(ran_tabu) :: rt
        integer        :: i, j
        ! calculate sampling size for each class
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0
        if( l_greedy )then
            ! completely greedy selection based on objective function value
            do i = 1, size(clssmp)
                do j = 1, clssmp(i)%nsample
                    states(clssmp(i)%pinds(j)) = 1
                end do
            end do
        else
            ! completely random selection (only class assignment matters)
            do i = 1, size(clssmp)
                allocate(pinds_left(clssmp(i)%pop), source=clssmp(i)%pinds)
                rt = ran_tabu(clssmp(i)%pop)
                call rt%shuffle(pinds_left)
                call rt%kill
                do j = 1, clssmp(i)%nsample
                    states(pinds_left(j)) = 1
                end do
                deallocate(pinds_left)
            end do
        endif
    end subroutine sample_balanced_1

    subroutine sample_balanced_2( self, clssmp, nptcls, frac_best, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)  ! data structure for balanced sampling
        integer,            intent(in)    :: nptcls     ! # particles to sample in total
        real,               intent(in)    :: frac_best  ! fraction of best scoring particles to sample from
        integer,            intent(inout) :: states(self%n)
        integer,            allocatable   :: pinds2sample(:)
        type(ran_tabu) :: rt
        integer        :: i, j, nbest
        ! calculate sampling size for each class
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0    
        do i = 1, size(clssmp)
            nbest = max(clssmp(i)%nsample, nint(frac_best * real(clssmp(i)%pop)))
            if( nbest == clssmp(i)%nsample )then
                ! completely greedy selection based on objective function value
                do j = 1, clssmp(i)%nsample
                    states(clssmp(i)%pinds(j)) = 1
                end do
            else
                ! random selection from the frac_best scoring ones
                allocate(pinds2sample(nbest), source=clssmp(i)%pinds(:nbest))
                rt = ran_tabu(nbest)
                call rt%shuffle(pinds2sample)
                call rt%kill
                do j = 1, clssmp(i)%nsample
                    states(pinds2sample(j)) = 1
                end do
                deallocate(pinds2sample)
            endif
        end do
    end subroutine sample_balanced_2

    subroutine sample_balanced_inv( self, clssmp, nptcls, frac_worst, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)  ! data structure for balanced sampling
        integer,            intent(in)    :: nptcls     ! # particles to sample in total
        real,               intent(in)    :: frac_worst ! fraction of worst scoring particles to sample from
        integer,            intent(inout) :: states(self%n)
        integer,            allocatable   :: pinds_rev(:), pinds2sample(:)
        type(ran_tabu) :: rt
        integer        :: i, j, nworst
        ! calculate sampling size for each class
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0    
        do i = 1, size(clssmp)
            nworst = max(clssmp(i)%nsample, nint(frac_worst * real(clssmp(i)%pop)))
            if( nworst == clssmp(i)%nsample )then
                ! completely deterministic selection based on objective function value
                do j = clssmp(i)%nsample,1,-1
                    states(clssmp(i)%pinds(j)) = 1
                end do
            else
                ! random selection from the frac_worst scoring ones
                allocate(pinds_rev(clssmp(i)%pop), source=clssmp(i)%pinds)
                call reverse(pinds_rev)
                allocate(pinds2sample(nworst), source=pinds_rev(:nworst))
                rt = ran_tabu(nworst)
                call rt%shuffle(pinds2sample)
                call rt%kill
                do j = 1, clssmp(i)%nsample
                    states(pinds2sample(j)) = 1
                end do
                deallocate(pinds_rev, pinds2sample)
            endif
        end do
    end subroutine sample_balanced_inv

    subroutine sample_balanced_parts( self, clssmp, nparts, states, nptcls_per_part )
        class(oris),        intent(inout) :: self
        type(class_sample), intent(inout) :: clssmp(:) ! data structure for balanced sampling
        integer,            intent(in)    :: nparts
        integer,            intent(inout) :: states(self%n)
        integer, optional,  intent(in)    :: nptcls_per_part
        integer :: nptcls, i, j, k, nptcls_eff
        ! calculate total # particles to sample
        nptcls_eff = self%count_state_gt_zero() 
        if( present(nptcls_per_part) )then
            nptcls = min(nptcls_eff, nparts * nptcls_per_part)
        else
            nptcls = nptcls_eff
        endif
        ! calculate sampling size for each class
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        ! create nparts balanced partitions of unique particles, state=1 is part 1, state=2 is part 2, etc.
        states = 0
        do i = 1, size(clssmp)
            j = 1
            do while( j < clssmp(i)%nsample )
                if( j + nparts - 1 > clssmp(i)%pop ) exit
                do k = 1, nparts
                    states(clssmp(i)%pinds(j)) = k
                    j = j + 1
                end do
            end do
        end do
    end subroutine sample_balanced_parts

    subroutine sample_ranked_parts( self, clssmp, nparts, states, nptcls_per_part )
        use simple_map_reduce, only: split_nobjs_even
        class(oris),        intent(inout) :: self
        type(class_sample), intent(inout) :: clssmp(:) ! data structure for balanced sampling
        integer,            intent(in)    :: nparts
        integer,            intent(inout) :: states(self%n)
        integer, optional,  intent(in)    :: nptcls_per_part
        integer, allocatable :: parts(:,:)
        integer :: nptcls, i, j, nptcls_eff, ipart
        ! calculate total # particles to sample
        nptcls_eff = self%count_state_gt_zero() 
        if( present(nptcls_per_part) )then
            nptcls = min(nptcls_eff, nparts * nptcls_per_part)
        else
            nptcls = nptcls_eff
        endif
        ! calculate sampling size for each class
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        ! create nparts balanced partitions of ranked particles, state=1 is rank 1, state=2 is rank 2, etc.
        states = 0
        do i = 1, size(clssmp)
            if( clssmp(i)%nsample <= nparts )then
                do j = 1, clssmp(i)%nsample
                    states(clssmp(i)%pinds(j)) = j
                end do
            else
                parts = split_nobjs_even(clssmp(i)%nsample, nparts)
                do ipart = 1, nparts
                    do j = parts(ipart,1), parts(ipart,2)
                        states(clssmp(i)%pinds(j)) = ipart
                    end do
                end do
                deallocate(parts)
            endif
        enddo
    end subroutine sample_ranked_parts

    subroutine balance_ptcls_within_cls( self, nptcls, pinds, maxpop, nparts )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nptcls, maxpop, nparts
        integer,     intent(in)    :: pinds(1:nptcls)
        real,    allocatable :: cls_scores(:)
        integer, allocatable :: cls_inds(:)
        integer :: classes(nptcls), i,icls,pop,n2reject,ncls,mpop
        if( maxpop < 1 )return
        classes = self%o(pinds(:))%get_class()
        ncls    = maxval(classes)
        mpop    = floor(real(maxpop)/real(nparts))
        do icls = 1,ncls
            cls_inds = pack((/(i,i=1,nptcls)/), mask=classes==icls)
            if( .not.allocated(cls_inds) ) cycle
            pop = size(cls_inds)
            if( pop <= mpop )cycle
            cls_scores = self%o(pinds(cls_inds(:)))%get('corr')
            call hpsort(cls_scores, cls_inds)
            n2reject = pop - mpop
            do i = 1,n2reject
                call self%o(pinds(cls_inds(i)))%set('w',0.)
            enddo
            deallocate(cls_inds,cls_scores)
        enddo
    end subroutine balance_ptcls_within_cls

    function get_sample_ind( self, incr_sampled ) result( sample_ind )
        class(oris), intent(in) :: self
        logical,     intent(in) :: incr_sampled
        integer :: i, sample_ind
        sample_ind = 0
        do i = 1, self%n
            if( self%o(i)%get_state() > 0 ) sample_ind = max(sample_ind, self%o(i)%get_sampled())
        end do
        if( incr_sampled ) sample_ind = sample_ind + 1
    end function get_sample_ind

    subroutine incr_sampled_updatecnt( self, inds, incr_sampled )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: inds(:)
        logical,     intent(in)    :: incr_sampled
        integer :: i, iptcl, sample_ind
        real    :: val
        sample_ind = self%get_sample_ind(incr_sampled)
        do i = 1, size(inds)
            iptcl = inds(i)
            val   = self%o(iptcl)%get('updatecnt')
            call self%o(iptcl)%set('updatecnt', val + 1.0)
            call self%o(iptcl)%set('sampled',   sample_ind)
        end do
    end subroutine incr_sampled_updatecnt

    logical function is_first_update( self, iter, iptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: iter, iptcl
        is_first_update = (self%o(iptcl)%get_int('updatecnt') == 1) .and. (iter > 1)
    end function is_first_update

    subroutine set_nonzero_updatecnt( self, updatecnt  )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: updatecnt
        integer :: i
        do i = 1,self%n
            if( self%o(i)%get('updatecnt') > 0 )then
                call self%o(i)%set('updatecnt', updatecnt)
            endif
        enddo
    end subroutine set_nonzero_updatecnt

    subroutine set_updatecnt( self, updatecnt, pinds )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: updatecnt, pinds(:)
        integer :: i, n
        ! zero them all
        do i = 1, self%n
            call self%o(i)%set('updatecnt', 0)
        enddo
        ! set the pinds to inputted value
        n = size(pinds)
        do i = 1, n
            call self%o(pinds(i))%set('updatecnt', updatecnt)
        enddo
    end subroutine set_updatecnt

    subroutine clean_entry( self, varflag1, varflag2 )
        class(oris),                 intent(inout) :: self
        character (len=*),           intent(in)    :: varflag1
        character (len=*), optional, intent(in)    :: varflag2
        logical :: varflag2_present
        integer :: i
        varflag2_present = present(varflag2)
        do i = 1,self%n
            call self%o(i)%delete_entry(varflag1)
            if( varflag2_present ) call self%o(i)%delete_entry(varflag2)
        enddo
    end subroutine clean_entry

    logical function has_been_sampled( self )
        class(oris), intent(inout) :: self
        integer :: i
        has_been_sampled =.false.
        do i = 1,self%n
            if( self%o(i)%get_state() > 0 )then
                if( self%o(i)%get_sampled() > 0 )then
                    has_been_sampled = .true.
                    exit
                endif
            endif
        end do
    end function has_been_sampled

    !>  \brief  check wether the orientation has any typical search parameter
    logical function has_been_searched( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
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

    function count_state_gt_zero( self ) result( cnt )
        class(oris), intent(in) :: self
        integer :: i, cnt
        cnt = 0
        do i = 1, self%n
            if( self%o(i)%get_state() > 0 ) cnt = cnt + 1
        end do
    end function count_state_gt_zero

    !>  \brief  joins the hashes into a string that represent the ith ori
    function ori2str( self, i ) result( str )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        type(string) :: str
        str = self%o(i)%ori2str()
    end function ori2str

    subroutine ori2json( self, i, json_ori, boxes)
        class(oris),                  intent(in)    :: self
        integer,                      intent(in)    :: i
        logical,          optional,   intent(in)    :: boxes
        type(json_value), pointer,    intent(inout) :: json_ori
        type(json_core)                             :: json
        logical :: l_boxes = .false.
        if(present(boxes)) l_boxes = boxes
        call self%o(i)%ori2json(json_ori, boxes=l_boxes)
        call json%add(json_ori, 'n', i) 
    end subroutine ori2json

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

    subroutine append_1( self, i, ori_in )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: ori_in
        integer,     intent(in)    :: i
        if( i < 0 .or. i > self%n )then
            THROW_WARN('index out of range; simple_oris % append')
            return
        endif
        call self%o(i)%append_ori(ori_in)
    end subroutine append_1

    subroutine append_2( self1, self2 )
        class(oris), intent(inout) :: self1
        class(oris), intent(in)    :: self2
        integer :: i,nprev
        logical :: self1_isptcl, self2_isptcl
        if( self2%n == 0 ) return
        self2_isptcl = self2%o(1)%is_particle()
        if( self1%n == 0 )then
            call self1%copy(self2, self2_isptcl)
        else
            self1_isptcl = self1%o(1)%is_particle()
            if( self1_isptcl.eqv.self2_isptcl )then
                nprev = self1%n
                call self1%reallocate(self1%n+self2%n)
                do i = 1,self2%n
                    self1%o(nprev+i) = self2%o(i)
                end do
            else
                THROW_HARD('self1 and self2 do not have equivalent is_ptcl status')
            endif
        endif
    end subroutine append_2

    subroutine copy_1( self_out, self_in, is_ptcl )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        logical,     intent(in)    :: is_ptcl
        integer :: i
        call self_out%new(self_in%n, is_ptcl)
        do i=1,self_in%n
            self_out%o(i) = self_in%o(i)
        end do
    end subroutine copy_1

    subroutine copy_2( self_out, self_in )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        logical :: is_ptcl
        integer :: i
        if(self_in%get_noris() == 0) then
            call self_out%kill()
        else
            is_ptcl = self_in%is_particle()
            call self_out%new(self_in%n, is_ptcl)
            do i=1,self_in%n
                self_out%o(i) = self_in%o(i)
            end do
        end if
    end subroutine copy_2

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

    subroutine delete_2Dclustering_1( self, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts, keepcls
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_2Dclustering(keepshifts, keepcls)
        end do
    end subroutine delete_2Dclustering_1

    subroutine delete_2Dclustering_2( self, i, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        logical, optional, intent(in)    :: keepshifts, keepcls
        if( i < 0 .or. i > self%n )then
            THROW_WARN('index out of range; simple_oris % delete_2Dclustering_2')
            return
        endif
        call self%o(i)%delete_2Dclustering(keepshifts, keepcls)
    end subroutine delete_2Dclustering_2

    subroutine delete_3Dalignment( self, keepshifts )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts
        integer :: i
        do i=1,self%n
            call self%o(i)%delete_3Dalignment(keepshifts)
        end do
    end subroutine delete_3Dalignment

    subroutine delete( self, ind )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: ind
        type(ori),        allocatable   :: tmporis(:)
        if( ind < 0 .or. ind > self%n )then
            THROW_WARN('index out of range; simple_oris % delete')
            return
        endif
        allocate( tmporis(self%n - 1) )
        tmporis=[self%o(1:ind-1), self%o(ind+1:self%n)]
        call self%new(self%n - 1, .false.)
        self%o=tmporis
        if(allocated(tmporis)) deallocate(tmporis)
    end subroutine delete

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

    subroutine set_state( self, i, state )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, state
        call self%o(i)%set_state(state)
    end subroutine set_state

    subroutine set_class( self, i, cls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, cls
        call self%o(i)%set_class(cls)
    end subroutine set_class

    subroutine set_stkind( self, i, stkind )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, stkind
        call self%o(i)%set_stkind(stkind)
    end subroutine set_stkind

    subroutine set_ogid( self, i, ogid )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, ogid
        call self%o(i)%set_ogid(ogid)
    end subroutine set_ogid

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

    subroutine set_3( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_3

    subroutine set_4( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real(dp),         intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_4
    
    subroutine set_5( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        class(string),    intent(in)    :: val
        call self%o(i)%set(key, val)
    end subroutine set_5

    subroutine set_dfx( self, i, dfx )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: dfx
        call self%o(i)%set_dfx(dfx)
    end subroutine set_dfx

    subroutine set_dfy( self, i, dfy )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: dfy
        call self%o(i)%set_dfy(dfy)
    end subroutine set_dfy

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

    subroutine set_all_3( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: vals(self%n)
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, vals(i))
        enddo
    end subroutine set_all_3

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

    subroutine set_all2single_3( self, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
        integer :: i
        do i=1,self%n
            call self%o(i)%set(which, ival)
        end do
    end subroutine set_all2single_3

    subroutine set_field2single_1( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            if( self%o(i)%get(trim(field)) == ind ) call self%o(i)%set(which, val)
        end do
        !$omp end parallel do 
    end subroutine set_field2single_1

    subroutine set_field2single_2( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            if( self%o(i)%get(trim(field)) == ind ) call self%o(i)%set(which, val)
        end do
        !$omp end parallel do 
    end subroutine set_field2single_2

    subroutine set_field2single_3( self, field, ind, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            if( self%o(i)%get(trim(field)) == ind ) call self%o(i)%set(which, ival)
        end do
        !$omp end parallel do 
    end subroutine set_field2single_3

    subroutine set_projs( self, e_space )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: e_space
        integer :: i
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%n
            call self%set(i, 'proj', e_space%find_closest_proj(self%o(i)))
        end do
        !$omp end parallel do
    end subroutine set_projs

    subroutine remap_projs( self, e_space, mapped_projs )
        class(oris), intent(in)  :: self
        class(oris), intent(in)  :: e_space
        integer,     intent(out) :: mapped_projs(self%n)
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

    subroutine zero_inpl( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%e3set(0.)
            call self%o(i)%set('x', 0.)
            call self%o(i)%set('y', 0.)
        end do
    end subroutine zero_inpl

    subroutine mul_shifts( self, mul )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: mul
        integer :: i
        do i=1,self%n
            call self%o(i)%set_shift(mul*self%o(i)%get_2Dshift())
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

    subroutine gau_rnd_shifts( self, std )
        class(oris),    intent(inout) :: self
        real,           intent(in)    :: std
        integer :: i
        do i=1,self%n
            call self%o(i)%gau_rnd_shift(std)
        end do
    end subroutine gau_rnd_shifts

    subroutine transfer_2Dshifts( self_out, self_in )
        class(oris), intent(inout) :: self_out
        type(oris),   intent(in)   :: self_in
        integer :: i
        do i = 1,self_out%n
            call self_out%o(i)%set_shift(self_in%o(i)%get_2Dshift())
        enddo
    end subroutine transfer_2Dshifts

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

    subroutine transfer_3Dparams_1( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
        call self%o(i)%transfer_3Dparams(o_in)
    end subroutine transfer_3Dparams_1

    subroutine transfer_3Dparams_2( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
        call self%o(i)%transfer_3Dparams(os_in%o(i_in))
    end subroutine transfer_3Dparams_2

    subroutine transfer_class_assignment( self_from, self_to )
        class(oris), intent(in)    :: self_from
        class(oris), intent(inout) :: self_to
        integer :: i
        if( self_from%n /= self_to%n ) THROW_HARD('Incongruent object instances')
        do i = 1,self_from%n
            call self_to%o(i)%set('class', self_from%o(i)%get_class())
        end do
    end subroutine transfer_class_assignment

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
                call self%o(iptcl)%set(trim(state_or_class), ipart)
            end do
        end do
        if(allocated(parts))deallocate(parts)
    end subroutine ini_tseries

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
            call self%o(i)%set_dfx(dfx)
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
                call self%o(i)%set_dfy(dfy)
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
            allocate( states(self%n) )
            rt = ran_tabu(self%n)
            call rt%balanced(nstates, states)
            do i=1,self%n
                state = self%o(i)%get_state()
                if( state /= 0 )then
                    call self%o(i)%set_state(states(i))
                endif
            end do
            call rt%kill
            deallocate(states)
        else if( nstates<=0)then
            THROW_HARD('invalid value for nstates; rnd_states')
        else
            ! nstates = 1; zero-preserving
            do i=1,self%n
                state = self%o(i)%get_state()
                if(  state /= 0 ) call self%o(i)%set_state(1)
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
            clsnr = self%get_class(i)
            if(clsnr == class) call self%set_class(i, class_merged)
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
        character(len=*), intent(in)    :: line
        call self%o(i)%str2ori(line, self%o(1)%is_particle())
    end subroutine str2ori

    subroutine str2ori_ctfparams_state_eo( self, i, line )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: line
        type(ori) :: o_tmp
        call o_tmp%str2ori(line, self%o(1)%is_particle())
        if( o_tmp%isthere('smpd')    ) call self%o(i)%set('smpd',    o_tmp%get('smpd'))
        if( o_tmp%isthere('kv')      ) call self%o(i)%set('kv',      o_tmp%get('kv'))
        if( o_tmp%isthere('cs')      ) call self%o(i)%set('cs',      o_tmp%get('cs'))
        if( o_tmp%isthere('fraca')   ) call self%o(i)%set('fraca',   o_tmp%get('fraca'))
        if( o_tmp%isthere('phshift') ) call self%o(i)%set('phshift', o_tmp%get('phshift'))
        if( o_tmp%isthere('dfx')     ) call self%o(i)%set_dfx(       o_tmp%get_dfx())
        if( o_tmp%isthere('dfy')     ) call self%o(i)%set_dfy(       o_tmp%get_dfy())
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

    subroutine gen_balanced_partitions( self, nparts, parts, err )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nparts
        integer, allocatable, intent(inout) :: parts(:,:)
        logical,              intent(out)   :: err
        real, allocatable :: rstates(:)
        integer :: nobjs_per_part, i, ipart, m
        err = .false.
        if( allocated(parts) ) deallocate(parts)
        if( .not.self%isthere('state') )then
            ! requires state field to operate )
            err = .true.
            return
        endif
        allocate(parts(nparts,2),source=0)
        nobjs_per_part = ceiling(real(self%n)/real(nparts))
        rstates    = self%get_all('state')
        parts(1,1) = 1
        ipart      = 1
        m          = 0
        do i = 1,self%n
            if( rstates(i) < 0.5 )cycle
            m = m+1
            if( m == nobjs_per_part )then
                m = 0
                parts(ipart,2) = i
                if( ipart < nparts ) parts(ipart+1,1) = i+1
                ipart = ipart + 1
                if( ipart == nparts )exit
            endif
        enddo
        parts(nparts,2) = self%n
        deallocate(rstates)
    end subroutine gen_balanced_partitions

    ! I/O

    !>  \brief  reads orientation info from file
    subroutine read( self, orifile, fromto, nst )
        class(oris),       intent(inout) :: self
        class(string),     intent(in)    :: orifile
        integer, optional, intent(in)    :: fromto(2)
        integer, optional, intent(out)   :: nst
        character(len=100) :: io_message
        integer :: file_stat, i, fnr, state, istart, iend
        if( .not. file_exists(orifile) )then
            THROW_HARD("the file you are trying to read: "//orifile%to_char()//' does not exist in cwd' )
        endif
        if( fname2ext(orifile) == 'bin' )then
            THROW_HARD('this method does not support binary files; read')
        endif
        io_message='No error'
        call fopen(fnr, FILE=orifile, STATUS='OLD', action='READ', iostat=file_stat,iomsg=io_message)
        call fileiochk("oris ; read ,Error when opening file for reading: "//orifile%to_char()//':'//trim(io_message), file_stat)
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
        call fclose(fnr)
    end subroutine read

    !>  \brief  reads CTF parameters and state info from file
    subroutine read_ctfparams_state_eo( self, ctfparamfile )
        class(oris),    intent(inout) :: self
        class(string),  intent(in)    :: ctfparamfile
        logical    :: params_are_there(10)
        integer    :: i
        type(oris) :: os_tmp
        if( .not. file_exists(ctfparamfile) )then
            THROW_HARD ("read_ctfparams_state_eo; The file you are trying to read: "//ctfparamfile%to_char()//' does not exist')
        endif
        if( ctfparamfile%has_substr('.bin') )then
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
            if( params_are_there(1)  ) call self%set(i, 'smpd',    os_tmp%get(i, 'smpd')   )
            if( params_are_there(2)  ) call self%set(i, 'kv',      os_tmp%get(i, 'kv')     )
            if( params_are_there(3)  ) call self%set(i, 'cs',      os_tmp%get(i, 'cs')     )
            if( params_are_there(4)  ) call self%set(i, 'fraca',   os_tmp%get(i, 'fraca')  )
            if( params_are_there(5)  ) call self%set(i, 'phshift', os_tmp%get(i, 'phshift'))
            if( params_are_there(6)  ) call self%set_dfx(i,        os_tmp%get_dfx(i)       )
            if( params_are_there(7)  ) call self%set_dfy(i,        os_tmp%get_dfy(i)       )
            if( params_are_there(8)  ) call self%set(i, 'angast',  os_tmp%get(i, 'angast') )
            if( params_are_there(9)  ) call self%set(i, 'state',   os_tmp%get(i, 'state')  )
            if( params_are_there(10) ) call self%set(i, 'eo',      os_tmp%get(i, 'eo')     )
        end do
        call os_tmp%kill
    end subroutine read_ctfparams_state_eo

    !>  \brief  writes orientation info to file
    subroutine write_1( self, orifile, fromto )
        class(oris),       intent(in) :: self
        class(string),     intent(in) :: orifile
        integer, optional, intent(in) :: fromto(2)
        character(len=100) :: io_message
        integer            :: file_stat, fnr, i, ffromto(2), cnt
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        call fopen(fnr, orifile, status='REPLACE', action='WRITE', iostat=file_stat, iomsg=io_message)
        call fileiochk(' Error opening file for writing: '//orifile%to_char()//' ; '//trim(io_message), file_stat)
        cnt = 0
        do i=ffromto(1),ffromto(2)
            cnt = cnt + 1
            call self%o(i)%write(fnr)
        end do
        call fclose(fnr)
    end subroutine write_1

    !>  \brief  writes orientation info to file
    subroutine write_2( self, i, orifile )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: orifile
        integer,       intent(in)    :: i
        integer :: fnr, file_stat
        call fopen(fnr, orifile, status='UNKNOWN', action='WRITE', position='APPEND', iostat=file_stat)
        call fileiochk( 'In: write_2, module: simple_oris.f90  opening '//orifile%to_char(), file_stat )
        call self%o(i)%write(fnr)
        call fclose(fnr)
    end subroutine write_2

    !>  \brief  writes object to BILD Chimera readable format
    subroutine write2bild( self, file )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: file
        integer :: i,funit, file_stat
        call fopen(funit, file, status='REPLACE', action='WRITE',iostat=file_stat)
        call fileiochk( 'In: write2bild, module: simple_oris.f90  opening '//file%to_char(), file_stat )
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
        call fclose(funit)
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
                    dfx = self%o(i)%get_dfx()+ran3()*dferr-dferr/2.
                    if( dfx > 0. ) exit
                end do
                call self%o(i)%set_dfx(dfx)
            endif
            if( self%o(i)%isthere('dfy') )then
                do
                    dfy = self%o(i)%get_dfy()+ran3()*dferr-dferr/2.
                    if( dfy > 0. ) exit
                end do
                call self%o(i)%set_dfy(dfy)
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
        med      = median_nocopy(vals4med)
    end function median_1

    !>  \brief  is for calculating variable statistics
    subroutine stats_1( self, which, ave, sdev, var, err )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        real,    allocatable :: vals(:)
        integer, allocatable :: states(:)
        states        = self%get_all_asint('state')
        vals          = self%get_all(which)
        vals          = pack(vals, mask=states > 0)
        call moment(vals, ave, sdev, var, err)
        deallocate(vals, states)
    end subroutine stats_1

    !>  \brief  is for calculating variable statistics
    subroutine stats_2( self, which, statvars, mask, nozero )
        class(oris),        intent(inout) :: self
        character(len=*),   intent(in)    :: which
        type(stats_struct), intent(inout) :: statvars
        logical,            intent(in)    :: mask(self%n)
        logical, optional,  intent(in)    :: nozero
        real,    allocatable :: vals(:)
        integer, allocatable :: states(:)
        logical :: err, nnozero
        real    :: var
        nnozero = .false.
        if( present(nozero) ) nnozero = nozero
        vals          = self%get_all(which)
        states        = self%get_all_asint('state')
        vals          = pack(vals, states > 0 .and. mask)
        if( nnozero ) vals = pack(vals, mask=vals(:) > TINY)
        call moment(vals, statvars%avg, statvars%sdev, var, err)
        statvars%minv = minval(vals)
        statvars%maxv = maxval(vals)
        statvars%med  = median_nocopy(vals)
        deallocate(vals)
    end subroutine stats_2

    !>  \brief  is for calculating variable statistics
    subroutine stats_3( self, which, ave, sdev, var, err, fromto )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        integer,          intent(in)    :: fromto(2)
        real,    allocatable            :: vals(:)
        integer, allocatable            :: states(:)
        states        = self%get_all_asint('state', fromto=fromto)
        vals          = self%get_all(which, fromto=fromto)
        vals          = pack(vals, mask=states > 0)
        call moment(vals, ave, sdev, var, err)
        deallocate(vals, states)
    end subroutine stats_3

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
            mystate = self%o(i)%get_state()
            if( mystate /= 0 )then
                val = self%o(i)%get(which)
                if( val < minv ) minv = val
                if( val > maxv ) maxv = val
            endif
        end do
    end subroutine minmax

    !>  \brief  is for generating evenly distributed projection directions
    subroutine spiral_1( self )
        class(oris), intent(inout) :: self
        real    :: h, theta, psi
        integer :: k
        if( self%n == 1 )then
            call self%o(1)%set_euler([0.,0.,0.])
        else if( self%n > 1 )then
            do k=1,self%n
                h     = -1.+2.*real(k-1)/real(self%n-1)
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
                if( allocated(avail) )deallocate(avail)
                allocate(avail(n), source=.false.)
                call tmp%new(n, self%o(1)%is_particle())
                call tmp%spiral_1
                do i = 1, n
                    if( tmp%o(i)%e1get() <= e1lim .and. tmp%o(i)%e2get() <= e2lim )&
                    &avail(i) = .true.
                end do
            end subroutine gen_c1

    end subroutine spiral_2

    !>  \brief  orders oris according to alignment score
    function order( self ) result( inds )
        class(oris), intent(inout) :: self
        real,    allocatable :: scores(:)
        integer, allocatable :: inds(:)
        integer :: i
        inds   = (/(i,i=1,self%n)/)
        scores = self%get_all('corr')
        call hpsort(scores, inds)
        call reverse(inds)
        deallocate( scores )
    end function order

    !>  \brief  orders clusters according to population
    function order_cls( self, ncls ) result( inds )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls
        integer, allocatable :: inds(:)
        real    :: classpops(ncls)
        integer :: i
        if(ncls <= 0) THROW_HARD('invalid number of classes; order_cls')
        allocate(inds(ncls))
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
        msk = cls_mask
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
        integer, allocatable :: order(:), states(:)
        integer :: i, j, lim, ind, n
        if( frac < 0.99 )then
            states = nint(self%get_all('state'))
            n      = count(states>0)
            lim    = nint(frac*real(n))
            order  = self%order() ! score ranking
            j = 0
            do i = 1,self%n
                ind = order(i)
                if( states(ind) == 0 )then
                    call self%o(ind)%set('w', 0.)
                else
                    j = j+1
                    if( j <= lim )then
                        call self%o(ind)%set('w', 1.)
                    else
                        call self%o(ind)%set('w', 0.)
                    endif
                endif
            end do
            deallocate(states,order)
        else
            call self%set_all2single('w', 1.)
        endif
    end subroutine calc_hard_weights

    !>  \brief  calculates soft weights based on score
    subroutine calc_soft_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        real,    allocatable :: scores(:), states(:), weights(:), weights_glob(:)
        integer, allocatable :: order(:), states_int(:)
        real    :: minw, mins
        integer :: i, lim, state, nstates
        if( self%isthere('corr') )then
            scores = self%get_all('corr')
            mins   = minval(scores, mask=scores > TINY)
            if( mins >= 0.8 )then
                call self%set_all2single('w', 1.0)
                return
            endif
            if( self%isthere('state') )then
                states = self%get_all('state')
            else
                allocate(states(self%n), source=1.0)
            endif
            nstates = nint(maxval(states))
            if( nstates == 1 )then
                weights = z_scores(scores, mask=scores > TINY .and. states > 0.5)
                minw    = minval(weights,  mask=scores > TINY .and. states > 0.5)
                where( scores > TINY .and. states > 0.5 )
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
                    weights = z_scores(scores, mask=scores > TINY .and. states_int == state)
                    minw    = minval(weights,  mask=scores > TINY .and. states_int == state)
                    where( scores > TINY .and. states_int == state ) weights_glob = weights + abs(minw)
                    deallocate(weights)
                end do
                call self%set_all('w', weights_glob)
                deallocate(states_int, weights_glob)
            endif
            if( frac < 0.99 )then
                ! in 3D frac operates globally, independent of state
                lim   = nint(frac*real(self%n))
                order = self%order() ! score ranking
                do i=1,self%n
                    if( i > lim ) call self%o(order(i))%set('w', 0.) ! nuke
                end do
                deallocate(order)
            endif
            deallocate(scores, states)
        else
            call self%set_all2single('w', 1.0)
        endif
    end subroutine calc_soft_weights

    !>  \brief  calculates soft weights based on score
    subroutine calc_cavg_soft_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        real,    allocatable :: scores(:)
        integer, allocatable :: order(:), states(:)
        integer :: i, lim, state, nstates, nnonzero
        if( .not.self%isthere('corr') )then
            call self%set_all2single('w', 1.0)
            return
        endif
        scores = self%get_all('corr')
        if( self%isthere('state') )then
            states = self%get_all_asint('state')
            nstates = maxval(states)
        else
            allocate(states(self%n), source=1)
            nstates = 1
        endif
        ! Weighing
        if( minval(scores, mask=(scores>TINY).and.states>0) >= 0.85 )then
            call self%set_all2single('w', 1.0)
            return
        endif
        do state = 1,nstates
            call calc_weights( state )
        enddo
        ! Thresholding a faction of particles
        if( frac < 0.999 )then
            order    = pack((/(i,i=1,self%n)/), mask=states > 0)
            nnonzero = size(order)
            lim      = nint((1.-frac)*real(nnonzero))
            scores   = scores(order(:))
            call hpsort(scores, order)
            do i = 1,lim,1
                call self%o(order(i))%set('w', 0.)
            end do
        endif
      contains

        subroutine calc_weights( s )
            integer,  intent(in) :: s
            real,    parameter   :: LOWER_BOUND_THRESHOLD = -4.0
            real,    allocatable :: weights(:), tmp(:)
            integer, allocatable :: inds(:), inds2(:)
            real    :: minw
            integer :: i, n
            n = count(states==s)
            if( n == 0 )then
                return  ! empty state
            else if( n < 5 )then
                ! not populated enough for weights to be calculated
                inds = pack((/(i,i=1,self%n)/), mask=states==s)
                allocate(weights(n),source=1.0)
            else
                inds    = pack((/(i,i=1,self%n)/), mask=states==s)
                weights = scores(inds(:))
                inds2 = pack((/(i,i=1,n)/), mask=weights > TINY)
                tmp   = weights(inds2)
                where(weights <= TINY) weights = 0.
                ! robust IQR-based normalization
                call robust_scaling(tmp)
                minw = max(LOWER_BOUND_THRESHOLD,minval(tmp))
               ! scale such that minimum=0 and 1 remains 1
                tmp = (tmp - minw) / (1.0 - minw)
                ! Values superior to 1 are set to 1
                do i = 1,size(inds2)
                    weights(inds2(i)) = min(1.0,max(0.0,tmp(i)))
                enddo
            endif
            do i = 1,n
                call self%o(inds(i))%set('w', weights(i))
            enddo
        end subroutine calc_weights

    end subroutine calc_cavg_soft_weights

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
                call self%get_pinds(icls, 'class', pinds)
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

    !>  \brief  calculates soft weights based on score
    subroutine calc_soft_weights2D( self )
        class(oris), intent(inout) :: self
        real,    allocatable :: scores(:), states(:)
        real,    allocatable :: weights(:), weights_glob(:)
        integer, allocatable :: classes(:)
        integer :: icls, pop, ncls
        real    :: minw
        if( self%isthere('corr') )then
            scores = self%get_all('corr')
            classes    = nint(self%get_all('class'))
            if( self%isthere('state') )then
                states = self%get_all('state')
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
                    where( scores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = 1.0
                else
                    if( count(scores > TINY .and. (classes == icls .and. states > 0.5)) <= MINCLSPOPLIM )then
                        where( scores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = 1.0
                    else
                        weights = z_scores(scores, mask=scores > TINY .and. (classes == icls .and. states > 0.5))
                        minw    = minval(weights,      mask=scores > TINY .and. (classes == icls .and. states > 0.5))
                        where( scores > TINY .and. (classes == icls .and. states > 0.5) ) weights_glob = weights + abs(minw)
                        deallocate(weights)
                    endif
                endif
            end do
            call self%set_all('w', weights_glob)
            deallocate(scores, classes, states, weights_glob)
        else
            call self%set_all2single('w', 1.)
        endif
    end subroutine calc_soft_weights2D

    pure real function euldist_1( self, i, j )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i, j
        euldist_1 = self%o(i).euldist.self%o(j)
    end function euldist_1

    pure real function euldist_2( self, i, o )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        class(ori),  intent(in) :: o
        euldist_2 = self%o(i).euldist.o
    end function euldist_2

    !>  \brief  to find the closest matching projection direction
    !! KEEP THIS ROUTINE SERIAL
    function find_closest_proj( self, o_in ) result( closest )
        class(oris), intent(in) :: self
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
                call self%o(i)%set_class(closest)
            end do
        else
            THROW_HARD('the number of discrete oris is too large; discretize')
        endif
    end subroutine discretize

    !>  \brief  to identify the indices of the k nearest projection neighbors (including self)
    subroutine nearest_proj_neighbors_1( self, k, nnmat )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: k
        integer,     intent(inout) :: nnmat(k,self%n)
        real      :: dists(self%n)
        integer   :: inds(self%n), i, j
        if( k >= self%n ) THROW_HARD('need to identify fewer nearest_proj_neighbors')
        !$omp parallel do default(shared) proc_bind(close) private(i,j,inds,dists)
        do i=1,self%n
            do j=1,self%n
                inds(j)  = j
                dists(j) = self%o(j).euldist.self%o(i)
            end do
            call hpsort(dists, inds)
            do j=1,k
                nnmat(j,i) = inds(j)
            end do
        end do
        !$omp end parallel do
    end subroutine nearest_proj_neighbors_1

    !>  \brief  to identify the nearest projection neighbors based on euldist threshold
    !! the policy here is based solely on angular distance and initialization of lnns is
    !! deferred to the calling unit, so that we can add additional neighborhoods on top of
    !! of each other to create more complex search spaces
    subroutine nearest_proj_neighbors_2( self, o, euldist_thres, lnns )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o
        real,        intent(in)    :: euldist_thres ! in degrees
        logical,     intent(inout) :: lnns(self%n)
        real    :: dists(self%n), euldist_thres_rad
        integer :: j
        euldist_thres_rad = deg2rad(euldist_thres)
        do j=1,self%n
            dists(j) = self%o(j).euldist.o
        end do
        where( dists <= euldist_thres_rad ) lnns = .true.
    end subroutine nearest_proj_neighbors_2

    !>  \brief  to identify the indices of the k nearest projection neighbors (including self)
    !! initialization of lnns is
    !! deferred to the calling unit, so that we can add additional neighborhoods on top of
    !! of each other to create more complex search spaces
    subroutine nearest_proj_neighbors_3( self, o, k, lnns )
        class(oris), intent(in)    :: self
        class(ori),  intent(in)    :: o
        integer,     intent(in)    :: k
        logical,     intent(inout) :: lnns(self%n)
        real    :: dists(self%n)
        integer :: inds(self%n), i, j
        if( k >= self%n ) THROW_HARD('need to identify fewer nearest_proj_neighbors')
        do i=1,self%n
            do j=1,self%n
                inds(j)  = j
                dists(j) = self%o(j).euldist.o
            end do
            call hpsort(dists, inds)
            do j=1,k
                lnns(inds(j)) = .true.
            end do
        end do
    end subroutine nearest_proj_neighbors_3

    subroutine extract_subspace( self, lnns, subself )
        class(oris), intent(in)    :: self
        logical,     intent(in)    :: lnns(self%n)
        class(oris), intent(inout) :: subself
        integer :: n, cnt, i
        n = count(lnns)
        if( n < 1 ) THROW_HARD('logical array for subspace generation empty')
        call subself%new(n, is_ptcl=self%o(1)%is_particle())
        cnt = 0
        do i = 1, self%n
            if( lnns(i) )then
                cnt            = cnt + 1
                subself%o(cnt) = self%o(i)
            endif
        enddo
    end subroutine extract_subspace

    subroutine detect_peaks( self, nnmat, corrs, peaks )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: nnmat(4,self%n) ! 4 because "self" is included
        real,        intent(in)    :: corrs(self%n)
        logical,     intent(inout) :: peaks(self%n)
        real, allocatable :: corrs_packed(:)
        integer :: i, npeaks
        real    :: corr_t
        corr_t = 0.0
        do i = 1,self%n
            if( i /= nnmat(1,i) ) THROW_HARD('self is not set to the first entry of the 2nd dimension')
        end do
        ! search for peaks
        do i = 1,self%n
            if( corrs(nnmat(1,i)) > corr_t )then
                peaks(i) = ( corrs(nnmat(1,i)) > corrs(nnmat(2,i)) .and.&
                &corrs(nnmat(1,i)) > corrs(nnmat(3,i)) .and.&
                &corrs(nnmat(1,i)) > corrs(nnmat(4,i)) )
            else
                peaks(i) = .false.
            endif
        end do
        npeaks = count(peaks .and. corrs > 0.)
        if( npeaks > 0 )then
            ! good/bad binning with Otsu's algorithm
            corrs_packed = pack(corrs, mask=peaks .and. corrs > 0.)
            call otsu(npeaks, corrs_packed, corr_t)
            deallocate(corrs_packed)
            where( corrs <= corr_t ) peaks = .false.
        endif
    end subroutine detect_peaks

    subroutine min_euldist( self, o_in, mindist )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_in
        real,        intent(inout) :: mindist
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
        real    :: dists(self%n), dists_max(self%n), x, nearest3(3), res
        integer :: i, j
        !$omp parallel do default(shared) proc_bind(close) private(j,i,dists,nearest3)
        do j=1,self%n
            do i=1,self%n
                if( i == j )then
                    dists(i) = huge(x)
                else
                    dists(i) = self%o(i).euldist.self%o(j)
                endif
            end do
            nearest3     = min3(dists)
            dists_max(j) = maxval(nearest3)
        end do
        !$omp end parallel do
        res = rad2deg(maxval(dists_max))
    end function find_angres

    !>  \brief  to find the correlation bound in extremal search
    function extremal_bound( self, thresh ) result( score_bound )
        class(oris),             intent(inout) :: self
        real,                    intent(in)    :: thresh ! is a fraction
        real,    allocatable  :: scores(:), scores_incl(:)
        logical, allocatable  :: incl(:)
        integer :: n_incl, thresh_ind
        real    :: score_bound
        if( .not.self%isthere('corr') )then
            THROW_HARD('Metric: corr is unpopulated; extremal_bound')
        endif
        ! fetch scores
        scores      = self%get_all('corr')
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
    subroutine set_extremal_vars( self, extr_init, extr_iter, iter, frac_srch_space, do_extr, iextr_lim, update_frac )
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
            mystate = self%o(i)%get_state()
            if( mystate == state ) call self%o(i)%map3dshift22d(sh3d)
        else
            call self%o(i)%map3dshift22d(sh3d)
        endif
    end subroutine map3dshift22d_2

    !>  \brief  Calculates the average particles offset within a class
    subroutine calc_avg_offset2D( self, class, offset_avg )
        class(oris),    intent(in)    :: self
        integer,        intent(in)    :: class
        real,           intent(inout) :: offset_avg(2)
        real(dp) :: avg(2)
        real     :: sh(2)
        integer  :: i, n
        n   = 0
        avg = 0.d0
        do i=1,self%n
            if( self%o(i)%isstatezero() ) cycle
            if( self%o(i)%get_class() == class )then
                call self%o(i)%calc_offset2D(sh)
                n   = n+1
                avg = avg + real(sh,dp)
            endif
        end do
        if( n > 0 )then
            avg        = avg / real(n,dp)
            offset_avg = real(avg)
        else
            offset_avg = 0.0
        endif
    end subroutine calc_avg_offset2D

    subroutine calc_avg_offset3D( self, offset_avg, state )
        class(oris),       intent(inout)    :: self
        real,              intent(inout) :: offset_avg(3)
        integer, optional, intent(in)    :: state
        integer,    parameter :: N = 300
        integer,    parameter :: PROJDIRMINPOP = 10
        type(oris)            :: spiral
        type(ori)             :: o, oxy, oxz, oyz
        integer,  allocatable :: closest_proj(:)
        real(dp) :: avg(3), offset(3), w(3), sumw(3)
        real     :: sh3d(3)
        integer  :: i,j, istate, pop, npop
        istate = 1
        if( present(state) ) istate = state
        call spiral%new(N,.true.)
        call spiral%spiral
        call spiral%set_all2single('state', 1.)
        allocate(closest_proj(self%n),source=0)
        call self%remap_projs(spiral, closest_proj)
        avg  = 0.d0
        npop = 0
        o    = ori(.true.)
        oxy  = ori(.true.)
        oxz  = ori(.true.)
        oyz  = ori(.true.)
        call oxy%set_euler([ 0., 0.,0.])
        call oxz%set_euler([90.,90.,0.])
        call oyz%set_euler([ 0.,90.,0.])
        !$omp parallel do default(shared) private(i,j,pop,offset,o,sh3d,w,sumw) &
        !$omp schedule(static) proc_bind(close) reduction(+:npop,avg)
        do i = 1,N
            pop    = 0
            offset = 0.d0
            sumw   = 0.d0
            call spiral%get_ori(i, o)
            do j = 1,self%n
                if( self%o(j)%get_state() /= istate ) cycle
                if( closest_proj(j) /= i ) cycle
                call self%o(j)%compose2dshift3d(sh3d)
                w(1)   = abs(sin(oyz.euldist.self%o(j)))
                w(2)   = abs(sin(oxz.euldist.self%o(j)))
                w(3)   = abs(sin(oxy.euldist.self%o(j)))
                offset = offset + w*real(sh3d,dp)
                sumw   = sumw   + w**2
                pop    = pop    + 1
            enddo
            if( pop < PROJDIRMINPOP ) cycle
            where( sumw < 1.d-6 )
                offset = 0.d0
            else where
                offset = offset / sumw
            end where
            avg  = avg  + offset
            npop = npop + 1
        enddo
        !$omp end parallel do
        avg        = avg / real(npop,dp)
        offset_avg = real(avg)
        ! cleanup
        call spiral%kill
        call o%kill
        call oxy%kill
        call oxz%kill
        call oyz%kill
        deallocate(closest_proj)
    end subroutine calc_avg_offset3D

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
            if( self%o(i)%get_class() == class )then
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
        use simple_linalg, only: vector_angle_norm
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
        os     = oris(100, is_ptcl=.false.)
        os2    = oris(100, is_ptcl=.false.)
        passed = .false.
        if( os%get_noris() == 100 ) passed = .true.
        if( .not. passed ) THROW_HARD('get_noris failed!')
        passed = .false.
        call os%set_euler(1, [1.,2.,3.])
        euls   = os%get_euler(1)
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
                call os%print(i)
            end do
            call os2%rnd_oris(5.)
            write(logfhandle,*) '********'
            do i=1,100
                call os2%print(i)
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
            y      = os%get(i,'y')
            y2     = os2%get(i,'y')
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
        call os%write(string('test_oris_rndoris.txt'))
        call os2%read(string('test_oris_rndoris.txt'))
        call os2%write(string('test_oris_rndoris_copy.txt'))
        corr = os%corr_oris(os2)
        if( corr > 0.99 ) passed = .true.
        if( .not. passed ) THROW_HARD('read/write failed')
        passed = .false.
        call os%rnd_states(5)
        call os%write(string('test_oris_rndoris_rndstates.txt'))
        if( os%corr_oris(os2) > 0.99 ) passed = .true.
        if( .not. passed ) THROW_HARD('statedoc read/write failed!')
        write(logfhandle,'(a)') '**info(simple_oris_unit_test, part3): testing calculators'
        passed = .false.
        call os%rnd_lps()
        call os%write(string('test_oris_rndoris_rndstates_rndlps.txt'))
        call os%spiral
        call os%write(string('test_oris_rndoris_rndstates_rndlps_spiral.txt'))
        call os%rnd_corrs()
        order = os%order()
        if( doprint )then
            do i=1,100
                call os%print(order(i))
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
            deallocate( self%o )
        endif
        self%n = 0
    end subroutine kill

end module simple_oris
