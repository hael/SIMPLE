!@descr: an agglomeration of orientations
module simple_oris
use simple_ori_api
use simple_ori, only: ori
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
    !======================================================================
    ! LIFE CYCLE / ALLOCATION (simple_oris_life.f90)
    !======================================================================
    generic            :: new => new_1, new_2
    procedure, private :: new_1, new_2
    procedure          :: reallocate
    procedure, private :: extract_subset_1, extract_subset_2
    generic            :: extract_subset => extract_subset_1, extract_subset_2
    procedure          :: kill_chash
    procedure          :: kill
    !======================================================================
    ! GETTERS / QUERIES (simple_oris_getters.f90)
    !======================================================================
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
    procedure          :: get_neven
    procedure          :: get_nodd
    procedure          :: get_nevenodd
    procedure          :: get_all_rmats
    procedure          :: count_state_gt_zero
    procedure          :: ori2str
    procedure          :: ori2json
    procedure          :: ori2prec
    procedure          :: prec2ori
    procedure          :: get_ctfvars
    procedure          :: has_been_sampled
    procedure          :: has_been_searched
    procedure          :: any_state_zero
    procedure          :: is_first_update
    procedure          :: get_update_frac
    procedure          :: get_class_sample_stats
    procedure          :: get_proj_sample_stats
    !======================================================================
    ! SETTERS / MUTATORS (simple_oris_setters.f90)
    !======================================================================
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
    !======================================================================
    ! SAMPLING / UPDATECOUNT (simple_oris_sampling.f90)
    !======================================================================
    procedure          :: select_particles_set
    procedure          :: sample4rec
    procedure          :: sample4update_all
    procedure          :: sample4update_rnd
    procedure          :: sample4update_cnt
    procedure          :: sample4update_class
    procedure          :: sample4update_reprod
    procedure          :: sample4update_updated
    procedure          :: sample4update_fillin
    procedure, private :: sample_balanced_1, sample_balanced_2
    generic            :: sample_balanced => sample_balanced_1, sample_balanced_2
    procedure          :: sample_balanced_inv
    procedure          :: sample_balanced_parts
    procedure          :: sample_ranked_parts
    procedure          :: balance_ptcls_within_cls
    procedure, private :: get_sample_ind
    procedure, private :: incr_sampled_updatecnt
    procedure          :: set_nonzero_updatecnt
    procedure          :: set_updatecnt
    procedure          :: clean_entry
    !======================================================================
    ! DISTANCES (simple_oris_dist.f90)
    !======================================================================
    procedure, private :: euldist_1, euldist_2
    generic            :: euldist => euldist_1, euldist_2
    procedure          :: min_euldist
    procedure          :: find_angres
    !======================================================================
    ! PRINT & FILE I/O (simple_oris_io.f90)
    !======================================================================
    procedure          :: print
    procedure          :: print_matrices
    procedure          :: read
    procedure          :: read_ctfparams_state_eo
    procedure, private :: write_1, write_2
    generic            :: write => write_1, write_2
    procedure          :: write2bild
    !======================================================================
    ! RESHAPE / PARTITION / REMAP (simple_oris_reshape.f90)
    !======================================================================
    procedure          :: compress
    procedure          :: split_state
    procedure          :: split_class
    procedure          :: expand_classes
    procedure          :: remap_cls
    procedure          :: merge_classes
    procedure          :: discretize
    procedure          :: extract_subspace
    !======================================================================
    ! TRANSFORMS / OFFSETS / ROTATIONS (simple_oris_transform.f90)
    !======================================================================
    procedure          :: round_shifts
    procedure          :: introd_alig_err
    procedure          :: introd_ctf_err
    procedure, private :: rot_1, rot_2
    generic            :: rot => rot_1, rot_2
    procedure, private :: rot_transp_1, rot_transp_2
    generic            :: rot_transp => rot_transp_1, rot_transp_2
    procedure, private :: map3dshift22d_1, map3dshift22d_2
    generic            :: map3dshift22d => map3dshift22d_1, map3dshift22d_2
    procedure          :: calc_avg_offset2D, calc_avg_offset3D
    procedure          :: mirror2d
    procedure          :: mirror3d
    procedure          :: add_shift2class
    !======================================================================
    ! STATS / ORDERING / PEAKS (simple_oris_stats.f90)
    !======================================================================
    procedure, private :: median_1
    generic            :: median => median_1
    procedure, private :: stats_1, stats_2, stats_3
    generic            :: stats => stats_1, stats_2, stats_3
    procedure          :: minmax
    procedure, private :: spiral_1, spiral_2
    generic            :: spiral => spiral_1, spiral_2
    procedure          :: order
    procedure          :: order_cls
    procedure          :: detect_peaks
    procedure          :: extremal_bound
    procedure          :: set_extremal_vars
    !======================================================================
    ! WEIGHTS / SCORING (simple_oris_weights.f90)
    !======================================================================
    procedure          :: calc_hard_weights
    procedure          :: calc_soft_weights, calc_cavg_soft_weights
    procedure          :: calc_hard_weights2D
    procedure          :: calc_soft_weights2D
    procedure          :: find_best_classes
    !======================================================================
    ! NEIGHBORS / CORRELATION / OVERLAP (simple_oris_neigh.f90)
    !======================================================================
    procedure          :: find_closest_proj
    procedure, private :: nearest_proj_neighbors_1, nearest_proj_neighbors_2, nearest_proj_neighbors_3
    generic            :: nearest_proj_neighbors => nearest_proj_neighbors_1, nearest_proj_neighbors_2, nearest_proj_neighbors_3
    procedure          :: corr_oris
    procedure, private :: diststat_1, diststat_2
    generic            :: diststat => diststat_1, diststat_2
    procedure          :: overlap
end type oris

interface oris
     module procedure constructor
end interface oris

interface

    !======================================================================
    ! LIFE CYCLE / ALLOCATION (simple_oris_life.f90)
    !======================================================================

    module function constructor( n, is_ptcl ) result( self )
        integer, intent(in) :: n
        logical, intent(in) :: is_ptcl
        type(oris) :: self
    end function constructor

    module subroutine new_1( self, n, is_ptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        logical,     intent(in)    :: is_ptcl
    end subroutine new_1

    module subroutine new_2( self, o_arr )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_arr(:)
    end subroutine new_2

    module subroutine reallocate( self, new_n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: new_n
    end subroutine reallocate

    module function extract_subset_1( self, from, to ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: from, to
        type(oris) :: self_sub
    end function extract_subset_1

    module function extract_subset_2( self, inds ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: inds(:)
        type(oris) :: self_sub
    end function extract_subset_2

    module subroutine kill_chash( self )
        class(oris), intent(inout) :: self
    end subroutine kill_chash

    module subroutine kill( self )
        class(oris), intent(inout) :: self
    end subroutine kill

    !======================================================================
    ! GETTERS / QUERIES (simple_oris_getters.f90)
    !======================================================================

    module pure function exists( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        logical :: exists
    end function exists

    module pure function e1get( self, i ) result( e1 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e1
    end function e1get

    module pure function e2get( self, i ) result( e2 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e2
    end function e2get

    module pure function e3get( self, i ) result( e3 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e3
    end function e3get

    module pure function get_euler( self, i ) result( euls )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: euls(3)
    end function get_euler

    module pure function get_noris( self, consider_state ) result( n )
        class(oris),       intent(in) :: self
        logical, optional, intent(in) :: consider_state
        integer :: n
    end function get_noris

    module subroutine get_ori( self, i, o )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: i
        type(ori),   intent(inout) :: o
    end subroutine get_ori

    module pure function get( self, i, key )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        real :: get
    end function get

    module pure function get_int( self, i, key )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        integer :: get_int
    end function get_int

    module function get_str( self, i, key ) result( val )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        type(string) :: val
    end function get_str

    module pure subroutine get_static( self, i, key, val )
        class(oris),      intent(in)  :: self
        integer,          intent(in)  :: i
        character(len=*), intent(in)  :: key
        character(len=*), intent(out) :: val
    end subroutine get_static

    module subroutine getter_1( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        type(string),     intent(inout) :: val
    end subroutine getter_1

    module subroutine getter_2( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
    end subroutine getter_2

    module subroutine getter_3( self, i, key, ival )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(inout) :: ival
    end subroutine getter_3

    module function get_all( self, key, fromto, nonzero ) result( arr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: fromto(2)
        logical, optional, intent(in) :: nonzero
        real, allocatable :: arr(:)
    end function get_all

    module function get_all_asint( self, key, fromto ) result( iarr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: fromto(2)
        integer, allocatable :: iarr(:)
    end function get_all_asint

    module function gen_ptcl_mask( self, key, ival, frac_best ) result( mask )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer,           intent(in) :: ival
        real,    optional, intent(in) :: frac_best
        logical, allocatable :: mask(:)
    end function gen_ptcl_mask

    module function get_all_sampled( self, key, state, lowerbound ) result( arr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: state
        real,    optional, intent(in) :: lowerbound
        real, allocatable :: arr(:)
    end function get_all_sampled

    module pure function get_mat( self, i ) result( mat )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: mat(3,3)
    end function get_mat

    module pure function get_normal( self, i ) result( normal )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: normal(3)
    end function get_normal

    module pure function get_2Dshift( self, i )  result( shvec )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: shvec(2)
    end function get_2Dshift

    module pure function get_dfx( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: get_dfx
    end function get_dfx

    module pure function get_dfy( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: get_dfy
    end function get_dfy

    module pure function get_state( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_state
    end function get_state

    module pure function get_class( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_class
    end function get_class

    module pure function get_proj( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_proj
    end function get_proj

    module function get_label_inds( self, label ) result( inds )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: label
        integer, allocatable :: inds(:)
    end function get_label_inds

    module pure function get_eo( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_eo
    end function get_eo

    module pure function get_fromp( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_fromp
    end function get_fromp

    module pure function get_top( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_top
    end function get_top

    module pure function get_sampled( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_sampled
    end function get_sampled

    module pure function get_updatecnt( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_updatecnt
    end function get_updatecnt

    module subroutine get_tseries_neighs( self, nsz, ptcls2neigh )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nsz
        integer, allocatable, intent(inout) :: ptcls2neigh(:,:)
    end subroutine get_tseries_neighs

    module pure function isthere_1( self, key ) result( is )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: key
        logical :: is
    end function isthere_1

    module pure function isthere_2( self, i, key ) result( is )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: i
        character(len=*),  intent(in) :: key
        logical :: is
    end function isthere_2

    module function ischar( self, i, key ) result( is )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: i
        character(len=*),  intent(in) :: key
        logical :: is
    end function ischar

    module pure function is_particle( self ) result( t )
        class(oris), intent(in) :: self
        logical :: t
    end function is_particle

    module function max_ori_strlen_trim( self )
        class(oris), intent(in) :: self
        integer :: max_ori_strlen_trim
    end function max_ori_strlen_trim

    module function get_n( self, label ) result( n )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: label
        integer :: n
    end function get_n

    module function get_pop_1( self, ind, label, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: ind
        character(len=*),  intent(in) :: label
        integer, optional, intent(in) :: eo
        integer :: pop
    end function get_pop_1

    module function get_pop_2( self, inds, labels, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: inds(:)
        character(len=*),  intent(in) :: labels(:)
        integer, optional, intent(in) :: eo
        integer :: pop
    end function get_pop_2

    module subroutine get_pops( self, pops, label, maxn, weight )
        class(oris),          intent(in)    :: self
        integer, allocatable, intent(out)   :: pops(:)
        character(len=*),     intent(in)    :: label
        integer, optional,    intent(in)    :: maxn
        logical, optional,    intent(in)    :: weight
    end subroutine get_pops

    module subroutine get_pinds( self, ind, label, indices, l_shuffle, l_require_updated )
        class(oris),          intent(in)  :: self
        character(len=*),     intent(in)  :: label
        integer,              intent(in)  :: ind
        integer, allocatable, intent(out) :: indices(:)
        logical, optional,    intent(in)  :: l_shuffle, l_require_updated
    end subroutine get_pinds

    module subroutine gen_mask( self, state, ind, label, l_mask, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state, ind
        character(len=*),     intent(in)    :: label
        logical, allocatable, intent(out)   :: l_mask(:)
        integer, optional,    intent(in)    :: fromto(2)
    end subroutine gen_mask

    module subroutine mask_from_state( self, state, l_mask, pinds, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state
        logical, allocatable, intent(inout) :: l_mask(:)
        integer, allocatable, intent(inout) :: pinds(:)
        integer, optional,    intent(in)    :: fromto(2)
    end subroutine mask_from_state

    module function get_all_normals( self ) result( normals )
        class(oris), intent(inout) :: self
        real, allocatable :: normals(:,:)
    end function get_all_normals

    module function states_exist( self, nstates, thres ) result( state_exists )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: nstates
        integer, optional, intent(in)    :: thres
        logical :: state_exists(nstates)
    end function states_exist

    module function projs_exist( self, nstates, nprojs, thres ) result( proj_exists )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: nstates
        integer,           intent(in)    :: nprojs
        integer, optional, intent(in)    :: thres
        logical :: proj_exists(nprojs,nstates)
    end function projs_exist

    module function get_arr( self, which, class, state ) result( vals )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        real, allocatable :: vals(:)
    end function get_arr

    module subroutine calc_sum( self, which, sum, cnt, class, state, fromto, mask )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: sum
        integer,           intent(out)   :: cnt
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
    end subroutine calc_sum

    module function get_sum( self, which, class, state, fromto, mask) result( sum )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        real :: sum
    end function get_sum

    module function get_avg( self, which, class, state, fromto, mask ) result( avg )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        real :: avg
    end function get_avg

    module function included( self ) result( incl )
        class(oris), intent(inout) :: self
        logical, allocatable :: incl(:)
    end function included

    module function get_neven( self )
        class(oris), intent(inout) :: self
        integer :: get_neven
    end function get_neven

    module function get_nodd( self )
        class(oris), intent(inout) :: self
        integer :: get_nodd
    end function get_nodd

    module function get_nevenodd( self )
        class(oris), intent(inout) :: self
        integer :: get_nevenodd
    end function get_nevenodd

    module function get_all_rmats( self ) result( mat )
        class(oris), intent(in) :: self
        real, allocatable :: mat(:,:,:)
    end function get_all_rmats

    module function count_state_gt_zero( self ) result( cnt )
        class(oris), intent(in) :: self
        integer :: cnt
    end function count_state_gt_zero

    module function ori2str( self, i ) result( str )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        type(string) :: str
    end function ori2str

    module subroutine ori2json( self, i, json_ori, boxes)
        class(oris),               intent(in)    :: self
        integer,                   intent(in)    :: i
        logical,       optional,   intent(in)    :: boxes
        type(json_value), pointer, intent(inout) :: json_ori
    end subroutine ori2json

    module subroutine ori2prec( self, i, prec )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: i
        real,        intent(inout) :: prec(N_PTCL_ORIPARAMS)
    end subroutine ori2prec

    module subroutine prec2ori( self, i, prec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: prec(N_PTCL_ORIPARAMS)
    end subroutine prec2ori

    module function get_ctfvars( self, i ) result( ctfvars )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        type(ctfparams) :: ctfvars
    end function get_ctfvars

    module function has_been_sampled( self )
        class(oris), intent(inout) :: self
        logical :: has_been_sampled
    end function has_been_sampled

    module function has_been_searched( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        logical :: has_been_searched
    end function has_been_searched

    module function any_state_zero( self )
        class(oris), intent(in) :: self
        logical :: any_state_zero
    end function any_state_zero

    module function is_first_update( self, iter, iptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: iter, iptcl
        logical :: is_first_update
    end function is_first_update

    module function get_update_frac( self ) result( update_frac )
        class(oris), intent(inout) :: self
        real :: update_frac
    end function get_update_frac

    module subroutine get_class_sample_stats( self, clsinds, clssmp, label )
        class(oris),                     intent(inout) :: self
        integer,                         intent(in)    :: clsinds(:)
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        character(len=*),      optional, intent(in)    :: label
    end subroutine get_class_sample_stats

    module subroutine get_proj_sample_stats( self, eulspace, clssmp )
        class(oris),                     intent(inout) :: self
        class(oris),                     intent(in)    :: eulspace
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
    end subroutine get_proj_sample_stats

    !======================================================================
    ! SETTERS / MUTATORS (simple_oris_setters.f90)
    !======================================================================

    module subroutine append_1( self, i, ori_in )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: ori_in
        integer,     intent(in)    :: i
    end subroutine append_1

    module subroutine append_2( self1, self2 )
        class(oris), intent(inout) :: self1
        class(oris), intent(in)    :: self2
    end subroutine append_2

    module subroutine copy_1( self_out, self_in, is_ptcl )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        logical,     intent(in)    :: is_ptcl
    end subroutine copy_1

    module subroutine copy_2( self_out, self_in )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
    end subroutine copy_2

    module subroutine reject( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
    end subroutine reject

    module subroutine delete_entry_1( self, key )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: key
    end subroutine delete_entry_1

    module subroutine delete_entry_2( self, ind, key )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: key
    end subroutine delete_entry_2

    module subroutine delete_2Dclustering_1( self, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts, keepcls
    end subroutine delete_2Dclustering_1

    module subroutine delete_2Dclustering_2( self, i, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        logical, optional, intent(in)    :: keepshifts, keepcls
    end subroutine delete_2Dclustering_2

    module subroutine delete_3Dalignment( self, keepshifts )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts
    end subroutine delete_3Dalignment

    module subroutine delete( self, ind )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ind
    end subroutine delete

    module subroutine transfer_2Dshifts( self_out, self_in )
        class(oris), intent(inout) :: self_out
        type(oris),  intent(in)    :: self_in
    end subroutine transfer_2Dshifts

    module subroutine transfer_2Dparams_1( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
    end subroutine transfer_2Dparams_1

    module subroutine transfer_2Dparams_2( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
    end subroutine transfer_2Dparams_2

    module subroutine transfer_3Dparams_1( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
    end subroutine transfer_3Dparams_1

    module subroutine transfer_3Dparams_2( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
    end subroutine transfer_3Dparams_2

    module subroutine transfer_class_assignment( self_from, self_to )
        class(oris), intent(in)    :: self_from
        class(oris), intent(inout) :: self_to
    end subroutine transfer_class_assignment

    module subroutine set_euler( self, i, euls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: euls(3)
    end subroutine set_euler

    module subroutine set_shift( self, i, vec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: vec(2)
    end subroutine set_shift

    module subroutine set_state( self, i, state )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, state
    end subroutine set_state

    module subroutine set_class( self, i, cls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, cls
    end subroutine set_class

    module subroutine set_stkind( self, i, stkind )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, stkind
    end subroutine set_stkind

    module subroutine set_ogid( self, i, ogid )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, ogid
    end subroutine set_ogid

    module subroutine e1set( self, i, e1 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e1
    end subroutine e1set

    module subroutine e2set( self, i, e2 )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: e2
    end subroutine e2set

    module subroutine e3set( self, i, e3 )
        class(oris), intent(inout) :: self
       integer,      intent(in)    :: i
        real,         intent(in)   :: e3
    end subroutine e3set

    module subroutine set_1( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
    end subroutine set_1

    module subroutine set_2( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
    end subroutine set_2

    module subroutine set_3( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: val
    end subroutine set_3

    module subroutine set_4( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real(dp),         intent(in)    :: val
    end subroutine set_4

    module subroutine set_5( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        class(string),    intent(in)    :: val
    end subroutine set_5

    module subroutine set_dfx( self, i, dfx )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: dfx
    end subroutine set_dfx

    module subroutine set_dfy( self, i, dfy )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: dfy
    end subroutine set_dfy

    module subroutine set_ori( self, i, o )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: o
    end subroutine set_ori

    module subroutine transfer_ori( self, i, self2transfer, i2transfer )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: self2transfer
        integer,     intent(in)    :: i, i2transfer
    end subroutine transfer_ori

    module subroutine set_all_1( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: vals(self%n)
    end subroutine set_all_1

    module subroutine set_all_2( self, which, vals )
        class(oris),           intent(inout) :: self
        character(len=*),      intent(in)    :: which
        character(len=STDLEN), intent(in)    :: vals(self%n)
    end subroutine set_all_2

    module subroutine set_all_3( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: vals(self%n)
    end subroutine set_all_3

    module subroutine set_all2single_1( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
    end subroutine set_all2single_1

    module subroutine set_all2single_2( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
    end subroutine set_all2single_2

    module subroutine set_all2single_3( self, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
    end subroutine set_all2single_3

    module subroutine set_field2single_1( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
    end subroutine set_field2single_1

    module subroutine set_field2single_2( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
    end subroutine set_field2single_2

    module subroutine set_field2single_3( self, field, ind, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
    end subroutine set_field2single_3

    module subroutine set_projs( self, e_space )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: e_space
    end subroutine set_projs

    module subroutine remap_projs( self, e_space, mapped_projs )
        class(oris), intent(in)  :: self
        class(oris), intent(in)  :: e_space
        integer,     intent(out) :: mapped_projs(self%n)
    end subroutine remap_projs

    module subroutine proj2class( self )
        class(oris), intent(inout) :: self
    end subroutine proj2class

    module subroutine e3swapsgn( self )
        class(oris), intent(inout) :: self
    end subroutine e3swapsgn

    module subroutine swape1e3( self )
        class(oris), intent(inout) :: self
    end subroutine swape1e3

    module subroutine zero( self, which )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
    end subroutine zero

    module subroutine zero_projs( self )
        class(oris), intent(inout) :: self
    end subroutine zero_projs

    module subroutine zero_shifts( self )
        class(oris), intent(inout) :: self
    end subroutine zero_shifts

    module subroutine zero_inpl( self )
        class(oris), intent(inout) :: self
    end subroutine zero_inpl

    module subroutine mul_shifts( self, mul )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: mul
    end subroutine mul_shifts

    module subroutine rnd_oris( self, trs, eullims )
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
    end subroutine rnd_oris

    module subroutine gau_rnd_shifts( self, std )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: std
    end subroutine gau_rnd_shifts

    module subroutine rnd_ori( self, i, trs, eullims )
        class(oris),    intent(inout) :: self
        integer,        intent(in)    :: i
        real, optional, intent(in)    :: trs
        real, optional, intent(inout) :: eullims(3,2)
    end subroutine rnd_ori

    module subroutine rnd_inpls( self, trs )
        class(oris),    intent(inout) :: self
        real, optional, intent(in)    :: trs
    end subroutine rnd_inpls

    module subroutine rnd_ctf( self, kv, cs, fraca, defocus, deferr, astigerr )
        class(oris),    intent(inout) :: self
        real,           intent(in)    :: kv, cs, fraca, defocus, deferr
        real, optional, intent(in)    :: astigerr
    end subroutine rnd_ctf

    module subroutine rnd_states( self, nstates )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nstates
    end subroutine rnd_states

    module subroutine rnd_lps( self )
        class(oris), intent(inout) :: self
    end subroutine rnd_lps

    module subroutine rnd_corrs( self )
        class(oris), intent(inout) :: self
    end subroutine rnd_corrs

    module subroutine rnd_proj_space( self, nsample, o_prev, thres, eullims )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: nsample
        class(ori), optional, intent(inout) :: o_prev
        real,       optional, intent(inout) :: eullims(3,2)
        real,       optional, intent(in)    :: thres
    end subroutine rnd_proj_space

    module subroutine revshsgn( self )
        class(oris), intent(inout) :: self
    end subroutine revshsgn

    module subroutine revorisgn( self )
        class(oris), intent(inout) :: self
    end subroutine revorisgn

    module subroutine ini_tseries( self, nsplit, state_or_class )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: nsplit
        character(len=*), intent(in)    :: state_or_class
    end subroutine ini_tseries

    module subroutine symmetrize( self, nsym )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsym
    end subroutine symmetrize

    module subroutine merge( self1, self2add )
        class(oris), intent(inout) :: self1, self2add
    end subroutine merge

    module subroutine partition_eo( self )
        class(oris), intent(inout) :: self
    end subroutine partition_eo

    module subroutine str2ori( self, i, line )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: line
    end subroutine str2ori

    module subroutine str2ori_ctfparams_state_eo( self, i, line )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: line
    end subroutine str2ori_ctfparams_state_eo

    module subroutine set_ctfvars( self, i, ctfvars )
        class(oris),     intent(inout) :: self
        integer,         intent(in)    :: i
        type(ctfparams), intent(in)    :: ctfvars
    end subroutine set_ctfvars

    module subroutine gen_balanced_partitions( self, nparts, parts, err )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nparts
        integer, allocatable, intent(inout) :: parts(:,:)
        logical,              intent(out)   :: err
    end subroutine gen_balanced_partitions

    !======================================================================
    ! SAMPLING / UPDATECOUNT (simple_oris_sampling.f90)
    !======================================================================

    module subroutine select_particles_set( self, nptcls, inds )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nptcls
        integer, allocatable, intent(inout) :: inds(:)
    end subroutine select_particles_set

    module subroutine sample4rec( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
    end subroutine sample4rec

    module subroutine sample4update_all( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_all

    module subroutine sample4update_rnd( self, fromto, update_frac, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_rnd

    module subroutine sample4update_cnt( self, fromto, update_frac, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_cnt

    module subroutine sample4update_class( self, clssmp, fromto, update_frac, nsamples, inds, &
                                        incr_sampled, l_greedy, frac_best )
        class(oris),          intent(inout) :: self
        type(class_sample),   intent(inout) :: clssmp(:)
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled, l_greedy
        real,    optional,    intent(in)    :: frac_best
    end subroutine sample4update_class

    module subroutine sample4update_reprod( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
    end subroutine sample4update_reprod

    module subroutine sample4update_updated( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_updated

    module subroutine sample4update_fillin( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_fillin

    module subroutine sample_balanced_1( self, clssmp, nptcls, l_greedy, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        logical,            intent(in)    :: l_greedy
        integer,            intent(inout) :: states(self%n)
    end subroutine sample_balanced_1

    module subroutine sample_balanced_2( self, clssmp, nptcls, frac_best, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        real,               intent(in)    :: frac_best
        integer,            intent(inout) :: states(self%n)
    end subroutine sample_balanced_2

    module subroutine sample_balanced_inv( self, clssmp, nptcls, frac_worst, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        real,               intent(in)    :: frac_worst
        integer,            intent(inout) :: states(self%n)
    end subroutine sample_balanced_inv

    module subroutine sample_balanced_parts( self, clssmp, nparts, states, nptcls_per_part )
        class(oris),        intent(inout) :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nparts
        integer,            intent(inout) :: states(self%n)
        integer, optional,  intent(in)    :: nptcls_per_part
    end subroutine sample_balanced_parts

    module subroutine sample_ranked_parts( self, clssmp, nparts, states, nptcls_per_part )
        class(oris),        intent(inout) :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nparts
        integer,            intent(inout) :: states(self%n)
        integer, optional,  intent(in)    :: nptcls_per_part
    end subroutine sample_ranked_parts

    module subroutine balance_ptcls_within_cls( self, nptcls, pinds, maxpop, nparts )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nptcls, maxpop, nparts
        integer,     intent(in)    :: pinds(1:nptcls)
    end subroutine balance_ptcls_within_cls

    module function get_sample_ind( self, incr_sampled ) result( sample_ind )
        class(oris), intent(in) :: self
        logical,     intent(in) :: incr_sampled
        integer :: sample_ind
    end function get_sample_ind

    module subroutine incr_sampled_updatecnt( self, inds, incr_sampled )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: inds(:)
        logical,     intent(in)    :: incr_sampled
    end subroutine incr_sampled_updatecnt

    module subroutine set_nonzero_updatecnt( self, updatecnt )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: updatecnt
    end subroutine set_nonzero_updatecnt

    module subroutine set_updatecnt( self, updatecnt, pinds )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: updatecnt
        integer,     intent(in)    :: pinds(:)
    end subroutine set_updatecnt

    module subroutine clean_entry( self, varflag1, varflag2 )
        class(oris),                 intent(inout) :: self
        character(len=*),            intent(in)    :: varflag1
        character(len=*),  optional, intent(in)    :: varflag2
    end subroutine clean_entry




    ! module subroutine sample4update_all( self, fromto, nsampled, sampled_inds, l_update_state )
    !     class(oris),       intent(inout) :: self
    !     integer,           intent(in)    :: fromto(2)
    !     integer,           intent(out)   :: nsampled
    !     integer, allocatable, intent(inout) :: sampled_inds(:)
    !     logical, optional, intent(in)    :: l_update_state
    ! end subroutine sample4update_all

    ! module subroutine sample4update_rnd( self, fromto, frac, nsampled, sampled_inds, l_update_state )
    !     class(oris),       intent(inout) :: self
    !     integer,           intent(in)    :: fromto(2)
    !     real,              intent(in)    :: frac
    !     integer,           intent(out)   :: nsampled
    !     integer, allocatable, intent(inout) :: sampled_inds(:)
    !     logical, optional, intent(in)    :: l_update_state
    ! end subroutine sample4update_rnd

    ! module subroutine sample4update_cnt( self, fromto, cnt, nsampled, sampled_inds, l_update_state )
    !     class(oris),       intent(inout) :: self
    !     integer,           intent(in)    :: fromto(2), cnt
    !     integer,           intent(out)   :: nsampled
    !     integer, allocatable, intent(inout) :: sampled_inds(:)
    !     logical, optional, intent(in)    :: l_update_state
    ! end subroutine sample4update_cnt

    ! module subroutine sample4update_class( self, class, nsampled, sampled_inds, l_update_state )
    !     class(oris),       intent(inout) :: self
    !     integer,           intent(in)    :: class
    !     integer,           intent(out)   :: nsampled
    !     integer, allocatable, intent(inout) :: sampled_inds(:)
    !     logical, optional, intent(in)    :: l_update_state
    ! end subroutine sample4update_class

    ! module subroutine sample4update_reprod( self, nsampled, sampled_inds, l_update_state )
    !     class(oris),       intent(inout) :: self
    !     integer,           intent(out)   :: nsampled
    !     integer, allocatable, intent(inout) :: sampled_inds(:)
    !     logical, optional, intent(in)    :: l_update_state
    ! end subroutine sample4update_reprod

    ! module subroutine sample4update_updated( self, nsampled, sampled_inds, l_update_state )
    !     class(oris),       intent(inout) :: self
    !     integer,           intent(out)   :: nsampled
    !     integer, allocatable, intent(inout) :: sampled_inds(:)
    !     logical, optional, intent(in)    :: l_update_state
    ! end subroutine sample4update_updated

    ! module subroutine sample4update_fillin( self, nsampled, sampled_inds, l_update_state )
    !     class(oris),       intent(inout) :: self
    !     integer,           intent(out)   :: nsampled
    !     integer, allocatable, intent(inout) :: sampled_inds(:)
    !     logical, optional, intent(in)    :: l_update_state
    ! end subroutine sample4update_fillin

    ! module subroutine sample_balanced_1( self, label, nsample, greedy, state )
    !     class(oris),       intent(inout) :: self
    !     character(len=*),  intent(in)    :: label
    !     integer,           intent(in)    :: nsample
    !     logical,           intent(in)    :: greedy
    !     integer, optional, intent(out)   :: state(:)
    ! end subroutine sample_balanced_1

    ! module subroutine sample_balanced_2( self, label, nsample, greedy, state )
    !     class(oris),       intent(inout) :: self
    !     character(len=*),  intent(in)    :: label
    !     integer,           intent(in)    :: nsample
    !     logical,           intent(in)    :: greedy
    !     integer, optional, intent(out)   :: state(:)
    ! end subroutine sample_balanced_2

    ! module subroutine sample_balanced_inv( self, label, nsample, frac_worst, state )
    !     class(oris),       intent(inout) :: self
    !     character(len=*),  intent(in)    :: label
    !     integer,           intent(in)    :: nsample
    !     real,              intent(in)    :: frac_worst
    !     integer, optional, intent(out)   :: state(:)
    ! end subroutine sample_balanced_inv

    ! module subroutine sample_balanced_parts( self, label, nparts, state, nptcls_per_part )
    !     class(oris),       intent(inout) :: self
    !     character(len=*),  intent(in)    :: label
    !     integer,           intent(in)    :: nparts
    !     integer, optional, intent(out)   :: state(:)
    !     integer, optional, intent(in)    :: nptcls_per_part
    ! end subroutine sample_balanced_parts

    ! module subroutine sample_ranked_parts( self, label, nparts, state, nptcls_per_part )
    !     class(oris),       intent(inout) :: self
    !     character(len=*),  intent(in)    :: label
    !     integer,           intent(in)    :: nparts
    !     integer, optional, intent(out)   :: state(:)
    !     integer, optional, intent(in)    :: nptcls_per_part
    ! end subroutine sample_ranked_parts

    ! module subroutine balance_ptcls_within_cls( self, label )
    !     class(oris),      intent(inout) :: self
    !     character(len=*), intent(in)    :: label
    ! end subroutine balance_ptcls_within_cls

    ! module subroutine get_sample_ind( self, i, ind )
    !     class(oris), intent(in)  :: self
    !     integer,     intent(in)  :: i
    !     integer,     intent(out) :: ind
    ! end subroutine get_sample_ind

    ! module subroutine incr_sampled_updatecnt( self, i )
    !     class(oris), intent(inout) :: self
    !     integer,     intent(in)    :: i
    ! end subroutine incr_sampled_updatecnt

    !======================================================================
    ! DISTANCES (simple_oris_dist.f90)
    !======================================================================

    module pure function euldist_1( self, i, j )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i, j
        real :: euldist_1
    end function euldist_1

    module pure function euldist_2( self, i, o )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        class(ori),  intent(in) :: o
        real :: euldist_2
    end function euldist_2

    module subroutine min_euldist( self, o_in, mindist )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_in
        real,        intent(inout) :: mindist
    end subroutine min_euldist

    module function find_angres( self ) result( res )
        class(oris), intent(in) :: self
        real :: res
    end function find_angres

    !======================================================================
    ! PRINT & FILE I/O (simple_oris_io.f90)
    !======================================================================

    module subroutine print( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
    end subroutine print

    module subroutine print_matrices( self )
        class(oris), intent(inout) :: self
    end subroutine print_matrices

    module subroutine read( self, orifile, fromto, nst )
        class(oris),       intent(inout) :: self
        class(string),     intent(in)    :: orifile
        integer, optional, intent(in)    :: fromto(2)
        integer, optional, intent(out)   :: nst
    end subroutine read

    module subroutine read_ctfparams_state_eo( self, ctfparamfile )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: ctfparamfile
    end subroutine read_ctfparams_state_eo

    module subroutine write_1( self, orifile, fromto )
        class(oris),       intent(in) :: self
        class(string),     intent(in) :: orifile
        integer, optional, intent(in) :: fromto(2)
    end subroutine write_1

    module subroutine write_2( self, i, orifile )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: orifile
        integer,       intent(in)    :: i
    end subroutine write_2

    module subroutine write2bild( self, file )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: file
    end subroutine write2bild

    !======================================================================
    ! RESHAPE / PARTITION / REMAP (simple_oris_reshape.f90)
    !======================================================================

    module subroutine compress( self, mask )
        class(oris), intent(inout) :: self
        logical,     intent(in)    :: mask(:)
    end subroutine compress

    module subroutine split_state( self, which )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: which
    end subroutine split_state

    module subroutine split_class( self, which )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: which
    end subroutine split_class

    module subroutine expand_classes( self, ncls_target )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls_target
    end subroutine expand_classes

    module subroutine remap_cls( self )
        class(oris), intent(inout) :: self
    end subroutine remap_cls

    module subroutine merge_classes( self, class_merged, class )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: class_merged, class
    end subroutine merge_classes

    module subroutine discretize( self, n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
    end subroutine discretize

    module subroutine extract_subspace( self, lnns, subself )
        class(oris), intent(in)    :: self
        logical,     intent(in)    :: lnns(self%n)
        class(oris), intent(inout) :: subself
    end subroutine extract_subspace

    !======================================================================
    ! TRANSFORMS / OFFSETS / ROTATIONS (simple_oris_transform.f90)
    !======================================================================

    module subroutine round_shifts( self )
        class(oris), intent(inout) :: self
    end subroutine round_shifts

    module subroutine introd_alig_err( self, angerr, sherr )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: angerr, sherr
    end subroutine introd_alig_err

    module subroutine introd_ctf_err( self, dferr )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: dferr
    end subroutine introd_ctf_err

    module subroutine rot_1( self, e )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: e
    end subroutine rot_1

    module subroutine rot_2( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
    end subroutine rot_2

    module subroutine rot_transp_1( self, e )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: e
    end subroutine rot_transp_1

    module subroutine rot_transp_2( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
    end subroutine rot_transp_2

    module subroutine map3dshift22d_1( self, sh3d, state )
        class(oris),       intent(inout) :: self
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
    end subroutine map3dshift22d_1

    module subroutine map3dshift22d_2( self, i, sh3d, state )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
    end subroutine map3dshift22d_2

    module subroutine calc_avg_offset2D( self, class, offset_avg )
        class(oris),    intent(in)    :: self
        integer,        intent(in)    :: class
        real,           intent(inout) :: offset_avg(2)
    end subroutine calc_avg_offset2D

    module subroutine calc_avg_offset3D( self, offset_avg, state )
        class(oris),       intent(inout) :: self
        real,              intent(inout) :: offset_avg(3)
        integer, optional, intent(in)    :: state
    end subroutine calc_avg_offset3D

    module subroutine mirror2d( self )
        class(oris), intent(inout) :: self
    end subroutine mirror2d

    module subroutine mirror3d( self )
        class(oris), intent(inout) :: self
    end subroutine mirror3d

    module subroutine add_shift2class( self, class, sh2d )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: class
        real,        intent(in)    :: sh2d(2)
    end subroutine add_shift2class

    !======================================================================
    ! STATS / ORDERING / PEAKS (simple_oris_stats.f90)
    !======================================================================

    module function median_1( self, which ) result( med )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real :: med
    end function median_1

    module subroutine stats_1( self, which, ave, sdev, var, err )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
    end subroutine stats_1

    module subroutine stats_2( self, which, statvars, mask, nozero )
        class(oris),        intent(inout) :: self
        character(len=*),   intent(in)    :: which
        type(stats_struct), intent(inout) :: statvars
        logical,            intent(in)    :: mask(self%n)
        logical, optional,  intent(in)    :: nozero
    end subroutine stats_2

    module subroutine stats_3( self, which, ave, sdev, var, err, fromto )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        integer,          intent(in)    :: fromto(2)
    end subroutine stats_3

    module subroutine minmax( self, which, minv, maxv )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: minv, maxv
    end subroutine minmax

    module subroutine spiral_1( self )
        class(oris), intent(inout) :: self
    end subroutine spiral_1

    module subroutine spiral_2( self, nsym, eullims )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsym
        real,        intent(in)    :: eullims(3,2)
    end subroutine spiral_2

    module function order( self ) result( inds )
        class(oris), intent(inout) :: self
        integer, allocatable :: inds(:)
    end function order

    module function order_cls( self, ncls ) result( inds )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls
        integer, allocatable :: inds(:)
    end function order_cls

    module subroutine detect_peaks( self, nnmat, corrs, peaks )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: nnmat(4,self%n)
        real,        intent(in)    :: corrs(self%n)
        logical,     intent(inout) :: peaks(self%n)
    end subroutine detect_peaks

    module function extremal_bound( self, thresh ) result( score_bound )
        class(oris),             intent(inout) :: self
        real,                    intent(in)    :: thresh
        real :: score_bound
    end function extremal_bound

    module subroutine set_extremal_vars( self, extr_init, extr_iter, iter, frac_srch_space, do_extr, iextr_lim, update_frac )
        class(oris),       intent(in)  :: self
        real,              intent(in)  :: extr_init
        integer,           intent(in)  :: extr_iter, iter
        real,              intent(in)  :: frac_srch_space
        logical,           intent(out) :: do_extr
        integer,           intent(out) :: iextr_lim
        real,    optional, intent(in)  :: update_frac
    end subroutine set_extremal_vars

    !======================================================================
    ! WEIGHTS / SCORING (simple_oris_weights.f90)
    !======================================================================

    module subroutine calc_hard_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
    end subroutine calc_hard_weights

    module subroutine calc_soft_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
    end subroutine calc_soft_weights

    module subroutine calc_cavg_soft_weights( self, frac )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
    end subroutine calc_cavg_soft_weights

    module subroutine calc_hard_weights2D( self, frac, ncls )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: frac
        integer,     intent(in)    :: ncls
    end subroutine calc_hard_weights2D

    module subroutine calc_soft_weights2D( self )
        class(oris), intent(inout) :: self
    end subroutine calc_soft_weights2D

    module subroutine find_best_classes( self, box, smpd, res_thresh, cls_mask, ndev )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: box
        real,        intent(in)    :: smpd, res_thresh, ndev
        logical,     intent(inout) :: cls_mask(1:self%n)
    end subroutine find_best_classes

    !======================================================================
    ! NEIGHBORS / CORRELATION / OVERLAP (simple_oris_neigh.f90)
    !======================================================================

    module function find_closest_proj( self, o_in ) result( closest )
        class(oris), intent(in) :: self
        class(ori),  intent(in) :: o_in
        integer :: closest
    end function find_closest_proj

    module subroutine nearest_proj_neighbors_1( self, k, nnmat )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: k
        integer,     intent(inout) :: nnmat(k,self%n)
    end subroutine nearest_proj_neighbors_1

    module subroutine nearest_proj_neighbors_2( self, o, euldist_thres, lnns )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o
        real,        intent(in)    :: euldist_thres
        logical,     intent(inout) :: lnns(self%n)
    end subroutine nearest_proj_neighbors_2

    module subroutine nearest_proj_neighbors_3( self, o, k, lnns )
        class(oris), intent(in)    :: self
        class(ori),  intent(in)    :: o
        integer,     intent(in)    :: k
        logical,     intent(inout) :: lnns(self%n)
    end subroutine nearest_proj_neighbors_3

    module function corr_oris( self1, self2 ) result( corr )
        class(oris), intent(inout) :: self1, self2
        real :: corr
    end function corr_oris

    module subroutine diststat_1( self, sumd, avgd, sdevd, mind, maxd )
        class(oris), intent(in)  :: self
        real,        intent(out) :: sumd, avgd, sdevd, mind, maxd
    end subroutine diststat_1

    module subroutine diststat_2( self1, self2, sumd, avgd, sdevd, mind, maxd )
        class(oris), intent(inout) :: self1, self2
        real,        intent(out)   :: sumd, avgd, sdevd, mind, maxd
    end subroutine diststat_2

    module real function overlap( self1, self2, which, state )
        class(oris),      intent(inout) :: self1, self2
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: state
    end function overlap

end interface

contains

    ! Test subroutine remains in main module
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
        os = oris(1000, is_ptcl=.false.)
        call os%spiral
        write(logfhandle,*) 'angres:      ', os%find_angres()
        write(logfhandle,'(a)') 'SIMPLE_ORIS_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_oris

end module simple_oris
