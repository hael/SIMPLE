!@descr: an agglomeration of orientations
module simple_oris
use simple_ori_api
use simple_ori, only: ori
use simple_rnd, only: rnd_key
implicit none

public :: oris, population_blend_weights, class_sample_quotas, class_sample_sweep
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
    generic            :: new => new_by_n, new_from_arr
    procedure, private :: new_by_n, new_from_arr
    procedure          :: reallocate
    procedure, private :: extract_subset_range, extract_subset_inds
    generic            :: extract_subset => extract_subset_range, extract_subset_inds
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
    procedure, private :: getter_str, getter_real, getter_int
    generic            :: getter => getter_str, getter_real, getter_int
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
    procedure          :: get_eo
    procedure          :: get_half
    procedure          :: nonslot_key
    procedure          :: strip_nonslot_keys
    procedure          :: get_fromp, get_top
    procedure          :: get_sampled
    procedure          :: get_updatecnt
    procedure          :: get_tseries_neighs
    procedure, private :: isthere_any, isthere_by_ind
    generic            :: isthere => isthere_any, isthere_by_ind
    procedure          :: ischar
    procedure          :: is_particle
    procedure          :: max_ori_strlen_trim
    procedure          :: get_n
    procedure, private :: get_pop_scalar, get_pop_array
    generic            :: get_pop => get_pop_scalar, get_pop_array
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
    procedure          :: get_update_frac
    procedure          :: get_group_update_counts
    procedure          :: get_update_rows
    procedure          :: get_state_update_fracs
    procedure          :: get_class_sample_stats
    !======================================================================
    ! SETTERS / MUTATORS (simple_oris_setters.f90)
    !======================================================================
    procedure, private :: append_at, append_oris
    generic            :: append => append_at, append_oris
    procedure, private :: copy_as_kind, copy_same_kind
    generic            :: copy => copy_as_kind, copy_same_kind
    procedure          :: reject
    procedure, private :: delete_entry_all, delete_entry_by_ind
    generic            :: delete_entry => delete_entry_all, delete_entry_by_ind
    procedure, private :: delete_2Dclustering_all, delete_2Dclustering_by_ind
    generic            :: delete_2Dclustering => delete_2Dclustering_all, delete_2Dclustering_by_ind
    procedure          :: delete_3Dalignment
    procedure          :: delete
    procedure          :: transfer_2Dshifts
    procedure, private :: transfer_2Dparams_ori, transfer_2Dparams_oris
    generic            :: transfer_2Dparams => transfer_2Dparams_ori, transfer_2Dparams_oris
    procedure, private :: transfer_3Dparams_ori, transfer_3Dparams_oris
    generic            :: transfer_3Dparams => transfer_3Dparams_ori, transfer_3Dparams_oris
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
    procedure, private :: set_real, set_char, set_int, set_dp, set_str
    generic            :: set => set_real, set_char, set_int, set_dp, set_str
    procedure          :: set_dfx, set_dfy
    procedure          :: set_ori
    procedure          :: transfer_ori
    procedure, private :: set_all_real, set_all_char, set_all_int
    generic            :: set_all => set_all_real, set_all_char, set_all_int
    procedure, private :: set_all2single_real, set_all2single_char, set_all2single_int
    generic            :: set_all2single => set_all2single_real, set_all2single_char, set_all2single_int
    procedure, private :: set_field2single_real, set_field2single_char, set_field2single_int
    generic            :: set_field2single => set_field2single_real, set_field2single_char, set_field2single_int
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
    procedure          :: gau_rnd_shifts
    procedure          :: rnd_ori
    procedure          :: rnd_inpls
    procedure          :: rnd_ctf
    procedure          :: rnd_states
    procedure          :: rnd_lps
    procedure          :: rnd_corrs
    procedure          :: rnd_proj_space
    procedure          :: rnd_cls
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
    procedure          :: sample4rec
    procedure          :: get_state_rec_pops
    procedure          :: sample4update_all
    procedure          :: sample4update_rnd
    procedure          :: sample4update_cnt
    procedure          :: sample4update_class
    procedure          :: sample4update_reprod
    procedure          :: sample4update_rescore
    procedure          :: sample4update_updated
    procedure          :: sample4update_fillin
    procedure          :: sample4update_missing
    procedure, private :: sample_balanced_mode, sample_balanced_frac_best
    generic            :: sample_balanced => sample_balanced_mode, sample_balanced_frac_best
    procedure          :: sample_balanced_inv
    procedure          :: sample_balanced_parts
    procedure          :: sample_ranked_parts
    procedure, private :: get_sample_ind
    procedure, private :: incr_sampled_updatecnt
    procedure          :: set_nonzero_updatecnt
    procedure          :: set_updatecnt
    procedure          :: clean_entry
    !======================================================================
    ! DISTANCES (simple_oris_dists.f90)
    !======================================================================
    procedure, private :: euldist_inds, euldist_ori
    generic            :: euldist => euldist_inds, euldist_ori
    procedure          :: min_euldist
    procedure          :: find_angres
    !======================================================================
    ! PRINT & FILE I/O (simple_oris_io.f90)
    !======================================================================
    procedure          :: print
    procedure          :: print_matrices
    procedure          :: read
    procedure          :: read_ctfparams_state_eo
    procedure, private :: write_range, write_by_ind
    generic            :: write => write_range, write_by_ind
    procedure          :: write2bild
    procedure          :: write_projdir_heatmap
    !======================================================================
    ! RESHAPE / PARTITION / REMAP (simple_oris_reshape.f90)
    !======================================================================
    procedure          :: compress
    procedure          :: split_state
    procedure          :: split_class
    procedure          :: expand_classes
    procedure          :: remap_cls
    procedure          :: reseed_classes
    procedure          :: merge_classes
    procedure          :: discretize
    procedure          :: extract_subspace
    !======================================================================
    ! TRANSFORMS / OFFSETS / ROTATIONS (simple_oris_transform.f90)
    !======================================================================
    procedure          :: round_shifts
    procedure          :: introd_alig_err
    procedure          :: introd_ctf_err
    procedure, private :: rot_all, rot_by_ind
    generic            :: rot => rot_all, rot_by_ind
    procedure, private :: rot_transp_all, rot_transp_by_ind
    generic            :: rot_transp => rot_transp_all, rot_transp_by_ind
    procedure, private :: map3dshift22d_all, map3dshift22d_by_ind
    generic            :: map3dshift22d => map3dshift22d_all, map3dshift22d_by_ind
    procedure          :: calc_avg_offset2D, calc_avg_offset3D
    procedure          :: mirror2d
    procedure          :: mirror3d
    procedure          :: add_shift2class
    !======================================================================
    ! STATS / ORDERING / PEAKS (simple_oris_stats.f90)
    !======================================================================
    procedure, private :: median_key
    generic            :: median => median_key
    procedure, private :: stats_all, stats_masked, stats_range
    generic            :: stats => stats_all, stats_masked, stats_range
    procedure          :: minmax
    procedure, private :: spiral_full, spiral_lims
    generic            :: spiral => spiral_full, spiral_lims
    procedure          :: order
    procedure          :: order_cls
    procedure          :: detect_peaks
    procedure          :: extremal_bound
    procedure          :: set_extremal_vars
    !======================================================================
    ! NEIGHBORS / CORRELATION / OVERLAP (simple_oris_neigh.f90)
    !======================================================================
    procedure          :: find_closest_proj
    procedure, private :: nearest_proj_neighbors_mat, nearest_proj_neighbors_thres, nearest_proj_neighbors_k
    generic            :: nearest_proj_neighbors => nearest_proj_neighbors_mat, nearest_proj_neighbors_thres, nearest_proj_neighbors_k
    procedure          :: replace_with_closest
    procedure, private :: diststat_within, diststat_between
    generic            :: diststat => diststat_within, diststat_between
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

    module subroutine new_by_n( self, n, is_ptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        logical,     intent(in)    :: is_ptcl
    end subroutine new_by_n

    module subroutine new_from_arr( self, o_arr )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_arr(:)
    end subroutine new_from_arr

    module subroutine reallocate( self, new_n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: new_n
    end subroutine reallocate

    module function extract_subset_range( self, from, to ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: from, to
        type(oris) :: self_sub
    end function extract_subset_range

    module function extract_subset_inds( self, inds ) result( self_sub )
        class(oris), intent(in) :: self
        integer,     intent(in) :: inds(:)
        type(oris) :: self_sub
    end function extract_subset_inds

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

    module subroutine getter_str( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        type(string),     intent(inout) :: val
    end subroutine getter_str

    module subroutine getter_real( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
    end subroutine getter_real

    module subroutine getter_int( self, i, key, ival )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(inout) :: ival
    end subroutine getter_int

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

    module pure function get_eo( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: get_eo
    end function get_eo

    module function get_half( self, i ) result( eo )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        integer :: eo
    end function get_half

    module function nonslot_key( self, i ) result( key )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        type(string) :: key
    end function nonslot_key

    module subroutine strip_nonslot_keys( self )
        class(oris), intent(inout) :: self
    end subroutine strip_nonslot_keys

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

    module pure function isthere_any( self, key ) result( is )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: key
        logical :: is
    end function isthere_any

    module pure function isthere_by_ind( self, i, key ) result( is )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: i
        character(len=*),  intent(in) :: key
        logical :: is
    end function isthere_by_ind

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

    module function get_pop_scalar( self, ind, label, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: ind
        character(len=*),  intent(in) :: label
        integer, optional, intent(in) :: eo
        integer :: pop
    end function get_pop_scalar

    module function get_pop_array( self, inds, labels, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: inds(:)
        character(len=*),  intent(in) :: labels(:)
        integer, optional, intent(in) :: eo
        integer :: pop
    end function get_pop_array

    module subroutine get_pops( self, pops, label, maxn, eo )
        class(oris),          intent(in)    :: self
        integer, allocatable, intent(out)   :: pops(:)
        character(len=*),     intent(in)    :: label
        integer, optional,    intent(in)    :: maxn, eo
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

    module function get_update_frac( self ) result( update_frac )
        class(oris), intent(inout) :: self
        real :: update_frac
    end function get_update_frac

    module subroutine get_update_rows( self, rep_rows, smp_rows, rec_rows )
        class(oris),          intent(inout) :: self
        integer, allocatable, intent(inout) :: rep_rows(:), smp_rows(:), rec_rows(:)
    end subroutine get_update_rows

    module subroutine get_group_update_counts( self, label, ngroups, nrep, nsmp )
        class(oris),          intent(inout) :: self
        character(len=*),     intent(in)    :: label
        integer,              intent(in)    :: ngroups
        integer, allocatable, intent(inout) :: nrep(:), nsmp(:)
    end subroutine get_group_update_counts

    module subroutine get_state_update_fracs( self, nstates, rho )
        class(oris),                     intent(inout) :: self
        integer,                         intent(in)    :: nstates
        real, allocatable,               intent(inout) :: rho(:)
    end subroutine get_state_update_fracs

    module subroutine get_class_sample_stats( self, clsinds, clssmp, label, score_key )
        class(oris),                     intent(inout) :: self
        integer,                         intent(in)    :: clsinds(:)
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        character(len=*),      optional, intent(in)    :: label, score_key
    end subroutine get_class_sample_stats

    !======================================================================
    ! SETTERS / MUTATORS (simple_oris_setters.f90)
    !======================================================================

    module subroutine append_at( self, i, ori_in )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: ori_in
        integer,     intent(in)    :: i
    end subroutine append_at

    module subroutine append_oris( self1, self2 )
        class(oris), intent(inout) :: self1
        class(oris), intent(in)    :: self2
    end subroutine append_oris

    module subroutine copy_as_kind( self_out, self_in, is_ptcl )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
        logical,     intent(in)    :: is_ptcl
    end subroutine copy_as_kind

    module subroutine copy_same_kind( self_out, self_in )
        class(oris), intent(inout) :: self_out
        class(oris), intent(in)    :: self_in
    end subroutine copy_same_kind

    module subroutine reject( self, i )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
    end subroutine reject

    module subroutine delete_entry_all( self, key )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: key
    end subroutine delete_entry_all

    module subroutine delete_entry_by_ind( self, ind, key )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: key
    end subroutine delete_entry_by_ind

    module subroutine delete_2Dclustering_all( self, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        logical, optional, intent(in)    :: keepshifts, keepcls
    end subroutine delete_2Dclustering_all

    module subroutine delete_2Dclustering_by_ind( self, i, keepshifts, keepcls )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        logical, optional, intent(in)    :: keepshifts, keepcls
    end subroutine delete_2Dclustering_by_ind

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

    module subroutine transfer_2Dparams_ori( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
    end subroutine transfer_2Dparams_ori

    module subroutine transfer_2Dparams_oris( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
    end subroutine transfer_2Dparams_oris

    module subroutine transfer_3Dparams_ori( self, i, o_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        type(ori),   intent(in)    :: o_in
    end subroutine transfer_3Dparams_ori

    module subroutine transfer_3Dparams_oris( self, i, os_in, i_in )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i, i_in
        class(oris), intent(in)    :: os_in
    end subroutine transfer_3Dparams_oris

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

    module subroutine set_real( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(in)    :: val
    end subroutine set_real

    module subroutine set_char( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        character(len=*), intent(in)    :: val
    end subroutine set_char

    module subroutine set_int( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(in)    :: val
    end subroutine set_int

    module subroutine set_dp( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real(dp),         intent(in)    :: val
    end subroutine set_dp

    module subroutine set_str( self, i, key, val )
        class(oris),      intent(inout) :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        class(string),    intent(in)    :: val
    end subroutine set_str

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

    module subroutine set_all_real( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: vals(self%n)
    end subroutine set_all_real

    module subroutine set_all_char( self, which, vals )
        class(oris),           intent(inout) :: self
        character(len=*),      intent(in)    :: which
        character(len=STDLEN), intent(in)    :: vals(self%n)
    end subroutine set_all_char

    module subroutine set_all_int( self, which, vals )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: vals(self%n)
    end subroutine set_all_int

    module subroutine set_all2single_real( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
    end subroutine set_all2single_real

    module subroutine set_all2single_char( self, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
    end subroutine set_all2single_char

    module subroutine set_all2single_int( self, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
    end subroutine set_all2single_int

    module subroutine set_field2single_real( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        real,             intent(in)    :: val
    end subroutine set_field2single_real

    module subroutine set_field2single_char( self, field, ind, which, val )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        character(len=*), intent(in)    :: val
    end subroutine set_field2single_char

    module subroutine set_field2single_int( self, field, ind, which, ival )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: field
        integer,          intent(in)    :: ind
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: ival
    end subroutine set_field2single_int

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

    module subroutine rnd_cls( self, ncls )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls
    end subroutine rnd_cls

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

    module subroutine sample4rec( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
    end subroutine sample4rec

    module subroutine get_state_rec_pops( self, nstates, pops )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: nstates
        integer, allocatable, intent(inout) :: pops(:)
    end subroutine get_state_rec_pops

    module subroutine sample4update_all( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_all

    module subroutine sample4update_rnd( self, fromto, update_frac, nsamples, inds, incr_sampled, key )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        type(rnd_key), optional, intent(in) :: key
    end subroutine sample4update_rnd

    module subroutine sample4update_cnt( self, fromto, update_frac, nsamples, inds, incr_sampled, allow_empty, key )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        logical, optional,    intent(in)    :: allow_empty
        type(rnd_key), optional, intent(in) :: key
    end subroutine sample4update_cnt

    module subroutine sample4update_class( self, clssmp, fromto, update_frac, nsamples, inds, &
                                        incr_sampled, l_greedy, frac_best, allow_empty, key )
        class(oris),          intent(inout) :: self
        type(class_sample),   intent(inout) :: clssmp(:)
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled, l_greedy
        real,    optional,    intent(in)    :: frac_best
        logical, optional,    intent(in)    :: allow_empty
        type(rnd_key), optional, intent(in) :: key
    end subroutine sample4update_class

    module subroutine sample4update_reprod( self, fromto, nsamples, inds, allow_empty )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical, optional,    intent(in)    :: allow_empty
    end subroutine sample4update_reprod

    module subroutine sample4update_rescore( self, fromto, nsamples, inds, allow_empty )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical, optional,    intent(in)    :: allow_empty
    end subroutine sample4update_rescore

    module subroutine sample4update_updated( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_updated

    module subroutine sample4update_fillin( self, fromto, update_frac, nsamples, inds, incr_sampled, allow_empty, key )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        logical, optional,    intent(in)    :: allow_empty
        type(rnd_key), optional, intent(in) :: key
    end subroutine sample4update_fillin

    module subroutine sample4update_missing( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
    end subroutine sample4update_missing

    module subroutine sample_balanced_mode( self, clssmp, nptcls, l_greedy, states, exact )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        logical,            intent(in)    :: l_greedy
        integer,            intent(inout) :: states(self%n)
        logical, optional,  intent(in)    :: exact
    end subroutine sample_balanced_mode

    module subroutine sample_balanced_frac_best( self, clssmp, nptcls, frac_best, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        real,               intent(in)    :: frac_best
        integer,            intent(inout) :: states(self%n)
    end subroutine sample_balanced_frac_best

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

    ! Population rule of the fractional blends (2D class sums, 3D trailing chains). Not
    ! type-bound: it acts on the counts of get_group_update_counts (or on applied weight masses),
    ! not on rows.
    elemental module subroutine population_blend_weights_counts( nrep, nsmp, mrep, s, w, mnew, ufrac )
        integer,        intent(in)  :: nrep, nsmp
        real,           intent(in)  :: mrep
        real,           intent(out) :: s, w, mnew
        real, optional, intent(in)  :: ufrac
    end subroutine population_blend_weights_counts

    elemental module subroutine population_blend_weights_mass( nrep, nsmp, mrep, s, w, mnew, ufrac )
        real,           intent(in)  :: nrep, nsmp
        real,           intent(in)  :: mrep
        real,           intent(out) :: s, w, mnew
        real, optional, intent(in)  :: ufrac
    end subroutine population_blend_weights_mass

    ! Expected particles per iteration of every sampling unit under the nested equal quota of
    ! sample4update_class (equal over groups, then over the units of a group, capped at populations)
    module subroutine class_sample_quotas( clssmp, ntarget, quotas )
        type(class_sample), intent(in)  :: clssmp(:)
        integer,            intent(in)  :: ntarget
        real,               intent(out) :: quotas(size(clssmp))
    end subroutine class_sample_quotas

    ! Iterations that visit every particle of every unit once: max over units of ceil(pop/quota)
    module function class_sample_sweep( clssmp, ntarget ) result( sweep )
        type(class_sample), intent(in) :: clssmp(:)
        integer,            intent(in) :: ntarget
        integer :: sweep
    end function class_sample_sweep

    !======================================================================
    ! DISTANCES (simple_oris_dist.f90)
    !======================================================================

    module pure function euldist_inds( self, i, j )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i, j
        real :: euldist_inds
    end function euldist_inds

    module pure function euldist_ori( self, i, o )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        class(ori),  intent(in) :: o
        real :: euldist_ori
    end function euldist_ori

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

    module subroutine write_range( self, orifile, fromto )
        class(oris),       intent(in) :: self
        class(string),     intent(in) :: orifile
        integer, optional, intent(in) :: fromto(2)
    end subroutine write_range

    module subroutine write_by_ind( self, i, orifile )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: orifile
        integer,       intent(in)    :: i
    end subroutine write_by_ind

    module subroutine write2bild( self, file )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: file
    end subroutine write2bild

    module subroutine write_projdir_heatmap( self, state, nspace, fname, legend )
        class(oris),            intent(in) :: self
        integer,                intent(in) :: state, nspace
        class(string),          intent(in) :: fname
        type(string), optional, intent(in) :: legend
    end subroutine write_projdir_heatmap

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

    module subroutine reseed_classes( self, clsinds, ncls_target, parent_of_seed, seed_pops, ndropped )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: clsinds(:)
        integer,              intent(in)    :: ncls_target
        integer, allocatable, intent(inout) :: parent_of_seed(:), seed_pops(:)
        integer,              intent(out)   :: ndropped
    end subroutine reseed_classes

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

    module subroutine rot_all( self, e )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: e
    end subroutine rot_all

    module subroutine rot_by_ind( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
    end subroutine rot_by_ind

    module subroutine rot_transp_all( self, e )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: e
    end subroutine rot_transp_all

    module subroutine rot_transp_by_ind( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
    end subroutine rot_transp_by_ind

    module subroutine map3dshift22d_all( self, sh3d, state )
        class(oris),       intent(inout) :: self
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
    end subroutine map3dshift22d_all

    module subroutine map3dshift22d_by_ind( self, i, sh3d, state )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
    end subroutine map3dshift22d_by_ind

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

    module function median_key( self, which ) result( med )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real :: med
    end function median_key

    module subroutine stats_all( self, which, ave, sdev, var, err )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
    end subroutine stats_all

    module subroutine stats_masked( self, which, statvars, mask, nozero )
        class(oris),        intent(inout) :: self
        character(len=*),   intent(in)    :: which
        type(stats_struct), intent(inout) :: statvars
        logical,            intent(in)    :: mask(self%n)
        logical, optional,  intent(in)    :: nozero
    end subroutine stats_masked

    module subroutine stats_range( self, which, ave, sdev, var, err, fromto )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: ave, sdev, var
        logical,          intent(out)   :: err
        integer,          intent(in)    :: fromto(2)
    end subroutine stats_range

    module subroutine minmax( self, which, minv, maxv )
        class(oris),      intent(inout) :: self
        character(len=*), intent(in)    :: which
        real,             intent(out)   :: minv, maxv
    end subroutine minmax

    module subroutine spiral_full( self )
        class(oris), intent(inout) :: self
    end subroutine spiral_full

    module subroutine spiral_lims( self, nsym, eullims )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nsym
        real,        intent(in)    :: eullims(3,2)
    end subroutine spiral_lims

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
    ! NEIGHBORS / CORRELATION / OVERLAP (simple_oris_neigh.f90)
    !======================================================================

    module function find_closest_proj( self, o_in ) result( closest )
        class(oris), intent(in) :: self
        class(ori),  intent(in) :: o_in
        integer :: closest
    end function find_closest_proj

    module subroutine nearest_proj_neighbors_mat( self, k, nnmat )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: k
        integer,     intent(inout) :: nnmat(k,self%n)
    end subroutine nearest_proj_neighbors_mat

    module subroutine nearest_proj_neighbors_thres( self, o, euldist_thres, lnns )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o
        real,        intent(in)    :: euldist_thres
        logical,     intent(inout) :: lnns(self%n)
    end subroutine nearest_proj_neighbors_thres

    module subroutine nearest_proj_neighbors_k( self, o, k, lnns )
        class(oris), intent(in)    :: self
        class(ori),  intent(in)    :: o
        integer,     intent(in)    :: k
        logical,     intent(inout) :: lnns(self%n)
    end subroutine nearest_proj_neighbors_k

    module subroutine replace_with_closest( self, other )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: other
    end subroutine replace_with_closest

    module subroutine diststat_within( self, sumd, avgd, sdevd, mind, maxd )
        class(oris), intent(in)  :: self
        real,        intent(out) :: sumd, avgd, sdevd, mind, maxd
    end subroutine diststat_within

    module subroutine diststat_between( self1, self2, sumd, avgd, sdevd, mind, maxd )
        class(oris), intent(inout) :: self1, self2
        real,        intent(out)   :: sumd, avgd, sdevd, mind, maxd
    end subroutine diststat_between

    module real function overlap( self1, self2, which, state )
        class(oris),      intent(inout) :: self1, self2
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: state
    end function overlap

end interface

! the population rule on row counts or on applied weight masses (same result for integer masses)
interface population_blend_weights
    module procedure population_blend_weights_counts, population_blend_weights_mass
end interface population_blend_weights

end module simple_oris
