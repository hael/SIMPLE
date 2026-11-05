!@descr: shared utility routines for probabilistic alignment tables
module simple_eul_prob_tab_utils
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,      only: TINY, dp, DPI, PI, logfhandle
use simple_linalg,    only: rad2deg
use simple_error,     only: simple_exception
use simple_math,      only: hpsort
use simple_oris,      only: oris
use simple_sym,       only: sym
use simple_ori_utils, only: euler2m
use simple_rnd,       only: ran3, rnd_key, rnd_scope_key, rnd_draw_key, rnd_run_seed, rnd_keyed_partial_shuffle, &
    &RND_PROG_REFINE3D, RND_PROG_REFINE2D, RND_PASS_ALIGN, RND_FAM_CAND
use simple_type_defs, only: OBJFUN_CC, OBJFUN_EUCLID, ptcl_ref
use simple_string,    only: string
use simple_string_utils, only: int2str_pad
use simple_fileio,    only: fopen, fclose, fileiochk, del_file, del_files, file2rarr
use simple_syslib,    only: file_exists
use simple_math_ft,   only: get_resarr
use simple_estimate_ssnr,   only: get_resolution
use simple_refine3D_fnames, only: refine3D_fsc_fname
implicit none

public :: angle_sampling, calc_athres, calc_num2sample
public :: prob_rnd_scope, keyed_candidate_order
public :: eulprob_corr_switch, eulprob_dist_switch
public :: record_seed_shift, write_fill_profile_bench
public :: sample_bounded_dists, sample_power_dists
public :: prob_candidate, prob_candidate_buffer, prob_candidate_store, prob_fill_profile
public :: prob_dense_entry, prob_refined_entry
public :: prob_posterior_stat, posterior_temperature, posterior_row_stats, draw_posterior_row, top_posterior_extras
public :: temperature_or_argmin, prob_beta_at, shift_prior_sigma2, SHIFT_PRIOR_SIGMA_MIN, SHIFT_PRIOR_LP_FRAC
public :: prob_neff_readout, PROB_NEFF_ROWS, PROB_NEFF_NMIN, PROB_NEFF_MIN_ROWS
public :: prob_readout, prob_readout_fname, write_prob_readout, read_prob_readout, add_prob_readout, report_prob_readout
public :: merge_prob_readouts
public :: PROB_READOUT_3D_DENSE, PROB_READOUT_3D_NEIGH, PROB_READOUT_2D_DENSE, PROB_READOUT_2D_SPARSE
public :: prob_insert_pose, prob_insert_table, prob_insert_poses, prob_insert_cone
public :: screen_pool_size, SCREEN_QUANTILE
public :: shift_needs_disk, disk_radius_px, shift_box_radius, next_shift_width, record_shift_state
public :: map_prior_weight, insert_pose_loss, update_screen_best, print_disk_set, insert_cone_from_fsc
private
#include "simple_local_flags.inc"

!> The floor of the shift prior's width: one pixel, or SHIFT_PRIOR_LP_FRAC of the matching wavelength when that
!! is wider, both in the pixels the shifts are measured in: the original sampling (smpd), whatever the crop, as
!! the polar searches and the stored shifts and increments are. Below a pixel the increments are the
!! refinement's sub-pixel moves, which a prior that tight would only fight; and increments measured under the
!! prior can never exceed it, so without a band-aware floor a coarse band's prior could pin the width wherever it
!! started
real, parameter :: SHIFT_PRIOR_SIGMA_MIN = 1.0
real, parameter :: SHIFT_PRIOR_LP_FRAC   = 0.25

!> the screen pool (shift_sigma > 0): the quantile of a particle's origin posterior the coarse-mapped references hold
real, parameter :: SCREEN_QUANTILE = 0.99

type :: prob_candidate
    integer :: iref   = 0
    integer :: inpl   = 0
    real    :: dist   = huge(1.0)
    real    :: x      = 0.
    real    :: y      = 0.
    logical :: has_sh = .false.
    real    :: sdev   = -1.     !< the shift's posterior width from the coarse map (px; below zero: none)
end type prob_candidate

!> The n_eff readout of the tables' posterior draws: the mean effective number of candidates exp(entropy) over a
!! subsample of rows, at the draw's tempering factor and at beta 1 (the physical posterior), logged every iteration
!! as POSTERIOR DRAW N_EFF. The rows are chosen from the global sample by rank, the same rows whatever the
!! partitioning (plan section 5): the particle's rank in the sample, with the stride the global sample size over
!! PROB_NEFF_ROWS, every stride-th rank among the first PROB_NEFF_ROWS strides is wanted; a wanted row is kept when
!! it holds at least PROB_NEFF_NMIN valid distances and a positive temperature, packed as distances above the row
!! minimum with the particle's temperature (a row of two candidates lands on one of its best two at any factor and
!! would pull the level of the rows where the choice is open). A part offers its rows with their global rank (its
!! rank offset in the sample); the master merges the parts' sums in part order. A readout, not a control: a factor
!! solved to pin n_eff cancels the data's own annealing at the band entries (the schedule is prob_beta flat with
!! the tail)
integer, parameter :: PROB_NEFF_ROWS     = 1000  !< rows of the subsample
integer, parameter :: PROB_NEFF_NMIN     = 6     !< valid candidates a row needs to be kept
integer, parameter :: PROB_NEFF_MIN_ROWS = 50    !< below which nothing is reported

type :: prob_neff_readout
    private
    integer :: stride = 1, nmax = 0, nmin = 2, nrows = 0, nvals = 0, rank0 = 0
    integer, allocatable :: offs(:)
    real,    allocatable :: dists(:), temps(:)
contains
    procedure :: new       => neff_readout_new
    procedure :: wanted    => neff_readout_wanted
    procedure :: offer     => neff_readout_offer
    procedure :: nkept     => neff_readout_nkept
    procedure :: mean_neff => neff_readout_mean_neff
    procedure :: sum_neff  => neff_readout_sum_neff
    procedure :: report    => neff_readout_report
    procedure :: kill      => neff_readout_kill
end type prob_neff_readout

!> the per-part publication of the readouts an assignment prints (plan section 5): the counts, the insert report's
!! sums and the n_eff readout's kept rows and sums, written by a distributed part as prob_readout_part<NN>.txt,
!! read and summed by the master in part order, then reported as the owner reports its own
integer, parameter :: PROB_READOUT_3D_DENSE = 1, PROB_READOUT_3D_NEIGH = 2, PROB_READOUT_2D_DENSE = 3, PROB_READOUT_2D_SPARSE = 4
type :: prob_readout
    integer        :: route = 0            !< PROB_READOUT_*
    integer        :: nptcls = 0           !< particles assigned (3D: the part's sampled particles; 2D: assigned)
    integer        :: nrefs  = 0           !< references or classes (the same on every part)
    integer        :: nmode  = 0           !< draws that landed on the row mode
    integer        :: ntemp_missing = 0    !< particles without a temperature
    integer        :: nbad   = 0           !< rows without a valid distance (3D dense)
    integer        :: nflat  = 0           !< rows with a flat distance (3D dense)
    integer        :: nkept  = 0           !< particles without a valid candidate, orientation kept (neighbourhood)
    integer        :: kmax   = 1           !< topk
    logical        :: l_extras = .false.   !< the insertion table was filled (a cone was open)
    integer        :: ins_nptcls = 0, ins_npose = 0, ins_ncapped = 0
    integer(int64) :: ins_nelig = 0_int64
    real(dp)       :: ins_wc_sum = 0.d0    !< the committed weights' sum
    integer        :: neff_nrows = 0
    real(dp)       :: neff_sum_draw = 0.d0 !< sum over the kept rows of exp(entropy) at the draw's factor
    real(dp)       :: neff_sum_phys = 0.d0 !< the same at factor 1
end type prob_readout

!> One extra pose of a particle for the weighted insertion (topk > 1): the reference (the full index over the states), its in-plane index, its
!! refined shift in the stored frame and its weight. The committed pose is the project's, with the weight the
!! table keeps per particle
type :: prob_insert_pose
    integer :: iref = 0, inpl = 0
    real    :: x = 0., y = 0., w = 0.
end type prob_insert_pose

!> The weighted insertion's table of one iteration: per particle (the table's order) the committed pose's weight
!! and its extra poses, selected by the assignment (top_posterior_extras); the search inserts each particle at its
!! committed pose with its weight and at its extras with theirs, every particle's weights summing to 1. The table
!! also counts the cone's occupancy (the eligible entries per particle and the particles with more of them than
!! the cap kmax - 1 admits) for the PROB_TAB_INSERT line, published by a distributed part through its counts
type :: prob_insert_table
    private
    integer :: nptcls = 0, kmax = 1, npose = 0, ncapped = 0
    integer(int64) :: nelig_sum = 0_int64
    integer,                allocatable :: pinds(:), nextra(:), first(:)
    real,                   allocatable :: wc(:)
    type(prob_insert_pose), allocatable :: poses(:)
contains
    procedure :: new    => insert_table_new
    procedure :: add    => insert_table_add
    procedure :: report => insert_table_report
    procedure :: counts => insert_table_counts
    procedure :: get_nptcls   => insert_table_get_nptcls
    procedure :: get_kmax     => insert_table_get_kmax
    procedure :: get_pind     => insert_table_get_pind
    procedure :: get_committed_weight => insert_table_get_committed_weight
    procedure :: get_nextra   => insert_table_get_nextra
    procedure :: get_extra    => insert_table_get_extra
    procedure :: get_mean_eligible => insert_table_get_mean_eligible
    procedure :: get_frac_capped   => insert_table_get_frac_capped
    procedure :: kill   => insert_table_kill
end type prob_insert_table

!> The weighted insertion's poses of one part's particles, resolved: per particle (by global index through loc)
!! the committed pose's weight and its extra poses as Euler triplets, absolute stored shifts and weights. setup
!! lays the structure out from the part's table (counts, weights, slots); the search fills a particle's slots
!! when it commits the particle's pose (strategy3D_prob, through candidate_pose, the rule the committed pose is
!! converted by); calc_3Drec inserts. The weights of a particle sum to 1, so the reconstruction's mass is the
!! particle count as without the extras
type :: prob_insert_poses
    integer              :: nptcls = 0, npose = 0
    integer, allocatable :: loc(:)       !< global particle index -> local index (0: no entry)
    real,    allocatable :: wc(:)        !< (nptcls) the committed pose's weight
    integer, allocatable :: nextra(:)    !< (nptcls)
    integer, allocatable :: first(:)     !< (nptcls) first extra in euler/shift/w
    real,    allocatable :: euler(:,:)   !< (3,npose)
    real,    allocatable :: shift(:,:)   !< (2,npose) absolute stored shifts
    real,    allocatable :: w(:)         !< (npose), 0 until the search fills the slot
contains
    procedure :: setup_range => insert_poses_setup_range
    procedure :: has         => insert_poses_has
    procedure :: kill        => insert_poses_kill
end type prob_insert_poses

!> The weighted insertion's cone: per state the resolution (the state map's FSC=0.5, Angstrom) that bounds how far
!! a pose may move the particle's rim from the committed pose, with the geometry the test needs (set_space, once
!! per assignment: the rotation matrices of the directions, of the rotation grid and of the point group). A pose is
!! within the cone of the committed one when the rotation between the two, the smallest over the point group's
!! equivalents of the candidate, taken as a displacement over the mask's radius, plus the distance between the two
!! stored shifts, is at most the resolution: a particle averaged over such poses is moved by no more than the map
!! resolves at its rim (a heuristic scale: at the cone's edge the displacement is one period of the resolution's
!! frequency; within it the weights at the particle's own temperature decide). Beyond the cone a pose is another
!! hypothesis, which the search's draw decides. A state without a resolution has no cone and its particles enter
!! at the committed pose alone
type :: prob_insert_cone
    private
    real, allocatable :: res(:)       !< (nstates) the resolution, Angstrom; not positive: no cone
    real, allocatable :: dmat(:,:,:)  !< (3,3,nspace) the directions' rotation matrices, in-plane angle zero
    real, allocatable :: rz(:,:,:)    !< (3,3,nrots) the rotation grid's in-plane matrices
    real, allocatable :: rsym(:,:,:)  !< (3,3,nsym) the point group's matrices
    real    :: diam = 0.              !< the mask diameter, Angstrom
    real    :: smpd = 0.              !< the sampling distance the stored shifts are measured in, Angstrom per pixel
    integer :: nspace = 0, nrots = 0, nsym = 0
contains
    procedure :: new       => insert_cone_new
    procedure :: set_space => insert_cone_set_space
    procedure :: is_open   => insert_cone_is_open
    procedure :: has_space => insert_cone_has_space
    procedure :: get_angle => insert_cone_get_angle
    procedure :: within    => insert_cone_within
    procedure :: kill      => insert_cone_kill
end type prob_insert_cone

! One entry of the dense table, a particle and a reference: the distance and the in-plane index. No default values:
! the fill and the read write every entry.
type :: prob_dense_entry
    real    :: dist
    integer :: inpl
end type prob_dense_entry

! A refined entry of the dense table: the particle (the table's local index),
! the compact reference, the shift of its search and the loss at that pose (the row's distance under the window
! search; under the polar shift map the drawn cell's posterior loss, its data term with the map's prior on the increment at the
! physical temperature, where the row holds the reference's evidence)
type :: prob_refined_entry
    integer :: i
    integer :: ri
    real    :: x
    real    :: y
    real    :: loss
    real    :: sdev   !< the shift's posterior width from the coarse map (px; below zero: none)
end type prob_refined_entry

!> Per-particle statistics of the pose posterior over a dense table row (prob_post_diag=yes): the mode
!! (state, direction, loss), the temperature, the mode's posterior mass, the effective number of
!! references exp(entropy), the number holding 99% of the mass and their posterior-weighted mean angular
!! distance from the mode in degrees (the mode's state only, symmetry-aware)
type :: prob_posterior_stat
    ! the production row (posterior-drawn in-plane, the best projs_ns shift-refined): mode, loss,
    ! temperature, mass, breadth and spread
    integer :: iptcl = 0, istate = 0, iproj = 0, k99 = 0, nvalid = 0
    real    :: lmin = 0., temp = 0., pmax = 0., neff = 0., spread = 0.
    ! the identical-pose null test, in nats: the kernel re-scored at the production mode's exact (reference,
    ! shift, grid rotation) minus the stored loss. mode_src says where the stored value came from and what the
    ! two columns are: 1 the shift search (refined): null_a the shifted-reference kernel at the stored
    ! particle-frame shift converted to the kernel's frame (shift_ptcl2ref; expected zero), null_b the stored
    ! shift taken as it is (the wrong frame, for contrast); 2 the seed bank (unrefined, seeded): null_a the bank
    ! kernel at the stored rotation (expected zero), null_b the shifted-reference kernel at the frame-converted
    ! seed (the two kernels differ by the angular Nyquist bin); 3 unshifted: both the kernel at zero shift
    integer :: mode_src = 0
    real    :: null_a = 0., null_b = 0.
    ! the deterministic row on the matcher's references (best grid rotation at the particle-frame seed the
    ! production search scored at, through the seed bank when there is one, no sampling, no refinement): the
    ! search test against the production row; rank_d the production mode's rank in it, dist_d the angle
    ! between the two modes
    integer :: iproj_d = 0, k99_d = 0, rank_d = 0
    real    :: lmin_d = 0., pmax_d = 0., neff_d = 0., spread_d = 0., dist_d = 0.
    ! the deterministic row on the particle's own real half (the diagnostic model of the independent
    ! even/odd maps, loaded beside the matcher's references): rank_h the matcher deterministic mode's rank
    ! in it, dist_h the angle between the two modes
    integer :: iproj_h = 0, k99_h = 0, rank_h = 0
    real    :: lmin_h = 0., pmax_h = 0., neff_h = 0., spread_h = 0., dist_h = 0.
    ! the same deterministic row on the other real half: the reference test against the own-half row;
    ! rank_x the own-half mode's rank in it, dist_x the angle between the two modes
    integer :: iproj_x = 0, k99_x = 0, rank_x = 0
    real    :: lmin_x = 0., pmax_x = 0., neff_x = 0., spread_x = 0., dist_x = 0.
    ! the self-reinforcement test, in nats, at the pose the particle went into its half with (its previous
    ! direction on the current grid, iproj_p, its previous in-plane rotation, zero residual shift): the
    ! other half's loss minus the own half's, T (L_other - L_own); positive when the own half, which holds
    ! the particle, favours that pose. ctrl_nats the same difference a half turn away in-plane, where the
    ! particle's own contribution does not correlate with it: the noise scale of self_nats. iproj_p 0: no
    ! previous pose (a virgin field, an invalid state)
    integer :: iproj_p = 0
    real    :: self_nats = 0., ctrl_nats = 0.
end type prob_posterior_stat

! Cost record of a partition's table fill for the bench file: wall time outside the parallel region,
! thread-seconds summed over threads.
type :: prob_fill_profile
    real(dp)       :: wall_region    = 0.d0   !< wall time of the parallel region, summed over batches
    real(dp)       :: ts_seed        = 0.d0   !< thread-seconds in the shift-seed estimation
    real(dp)       :: ts_score       = 0.d0   !< thread-seconds in the all-reference score loop
    real(dp)       :: ts_refine      = 0.d0   !< thread-seconds in the refinement of the best references
    integer        :: nthr           = 1      !< threads of the parallel region
    integer        :: nrefs          = 0      !< references scored per particle
    integer        :: projs_ns       = 0      !< references refined per particle (largest over batches)
    integer(int64) :: nptcls         = 0_int64 !< particles filled
    integer(int64) :: nrefined       = 0_int64 !< references handed to the shift search
    integer(int64) :: objective_evals = 0_int64 !< objective evaluations of the search objects (rotations count nrots)
    integer(int64) :: gradient_evals  = 0_int64 !< L-BFGS-B gradient evaluations, E
    integer(int64) :: alternations    = 0_int64 !< L-BFGS-B / in-plane alternations
    integer(int64) :: grid_points     = 0_int64 !< coarse-grid points scored
    integer(int64) :: minimize_calls  = 0_int64 !< shift searches started
    integer(int64) :: bank_builds     = 0_int64 !< shifted particle memos built
    integer(int64) :: nscored         = 0_int64 !< references scored, summed over particles (neighbourhood tables; 0 in the full table, which scores nrefs per particle)
    integer(int64) :: screens         = 0_int64 !< coarse maps over either search domain
    integer(int64) :: disk_ptcls      = 0_int64 !< particles that searched the disk (shift_sigma > 0, the centring stage)
    real(dp)       :: map_times(6)    = 0.d0   !< construction, two reserved slots, selection, posterior reduction, setup wall time
    integer(int64) :: map_counts(8)   = 0_int64
    integer        :: map_geometry(3,3) = 0   !< radius/direction/order counts of full/disk/window maps
    integer(int64) :: fit_counts(0:5) = 0_int64 !< quadratic grid step: (0) fits kept, (1:5) fallbacks per reason; all zero in a single-state table
    !> Is the refinement budget binding? Per particle of the dense fill, the rank in the unrefined row of the entry
    !! that is lowest after the shift refinement: (0) rank 1, the unrefined argmin kept; (1:10) the rank by decile
    !! of the refined count projs_ns (rank 1 is in decile 1); (11) an unrefined entry below every refined one. Mass
    !! in the top deciles says the truth lands near the cap and more refinements would help; mass at rank 1 says
    !! the cap is slack
    integer(int64) :: rank_counts(0:11) = 0_int64
end type prob_fill_profile

! Thread-owned append buffer used while candidate counts are not known in advance.
! Keeping this common and mode-independent avoids per-search table variants.
type :: prob_candidate_buffer
    integer, allocatable              :: particle_indices(:)
    type(prob_candidate), allocatable :: candidates(:)
    integer                           :: nused = 0
contains
    procedure :: append            => append_candidate
    procedure :: append_or_replace => append_or_replace_candidate
    procedure :: get_inpl          => get_buffer_candidate_inpl
    procedure :: kill              => kill_candidate_buffer
end type prob_candidate_buffer

type :: prob_candidate_store
    integer,                 allocatable :: pinds(:)
    integer(int64),          allocatable :: offsets(:)
    type(prob_candidate),    allocatable :: candidates(:)
    real,                    allocatable :: seed_shifts(:,:)
    logical,                 allocatable :: seed_has_sh(:)
contains
    procedure :: new_fixed => new_fixed_candidate_store
    procedure :: new_ragged => new_ragged_candidate_store
    procedure :: kill      => kill_candidate_store
end type prob_candidate_store

interface angle_sampling
    module procedure angle_sampling_1
    module procedure angle_sampling_2
end interface

contains

    !> The best losses bound width readout only; final candidate ordering is unchanged.
    pure subroutine update_screen_best( loss, best_losses )
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real, intent(in) :: loss
        real, intent(inout) :: best_losses(:)
        integer :: worst
        if( size(best_losses) < 1 ) return
        if( .not. ieee_is_finite(loss) ) return
        worst = maxloc(best_losses,dim=1)
        if( loss < best_losses(worst) ) best_losses(worst) = loss
    end subroutine update_screen_best

    subroutine append_candidate( self, particle_index, candidate )
        class(prob_candidate_buffer), intent(inout) :: self
        integer,                      intent(in)    :: particle_index
        type(prob_candidate),         intent(in)    :: candidate
        integer, allocatable :: particle_indices_tmp(:)
        type(prob_candidate), allocatable :: candidates_tmp(:)
        integer :: old_capacity, new_capacity
        if( .not. allocated(self%candidates) )then
            allocate(self%particle_indices(256), self%candidates(256))
        else if( self%nused == size(self%candidates) )then
            old_capacity = size(self%candidates)
            new_capacity = max(old_capacity + 1, 2 * old_capacity)
            allocate(particle_indices_tmp(new_capacity), candidates_tmp(new_capacity))
            particle_indices_tmp(1:old_capacity) = self%particle_indices
            candidates_tmp(1:old_capacity)       = self%candidates
            call move_alloc(particle_indices_tmp, self%particle_indices)
            call move_alloc(candidates_tmp,       self%candidates)
        endif
        self%nused = self%nused + 1
        self%particle_indices(self%nused) = particle_index
        self%candidates(self%nused)       = candidate
    end subroutine append_candidate

    subroutine append_or_replace_candidate( self, particle_index, candidate )
        class(prob_candidate_buffer), intent(inout) :: self
        integer,                      intent(in)    :: particle_index
        type(prob_candidate),         intent(in)    :: candidate
        integer :: i
        do i = self%nused, 1, -1
            if( self%particle_indices(i) /= particle_index ) exit
            if( self%candidates(i)%iref == candidate%iref )then
                self%candidates(i) = candidate
                return
            endif
        enddo
        call self%append(particle_index,candidate)
    end subroutine append_or_replace_candidate

    integer function get_buffer_candidate_inpl( self, particle_index, full_ref ) result(inpl)
        class(prob_candidate_buffer), intent(in) :: self
        integer, intent(in) :: particle_index, full_ref
        integer :: i
        inpl = 0
        do i = self%nused, 1, -1
            if( self%particle_indices(i) /= particle_index ) exit
            if( self%candidates(i)%iref == full_ref )then
                inpl = self%candidates(i)%inpl
                return
            endif
        enddo
    end function get_buffer_candidate_inpl

    subroutine kill_candidate_buffer( self )
        class(prob_candidate_buffer), intent(inout) :: self
        if( allocated(self%particle_indices) ) deallocate(self%particle_indices)
        if( allocated(self%candidates)        ) deallocate(self%candidates)
        self%nused = 0
    end subroutine kill_candidate_buffer

    subroutine new_fixed_candidate_store( self, pinds, ncandidates_per_particle )
        class(prob_candidate_store), intent(inout) :: self
        integer,                     intent(in)    :: pinds(:)
        integer,                     intent(in)    :: ncandidates_per_particle
        integer :: i, nptcls
        integer(int64) :: ntot
        call self%kill
        if( ncandidates_per_particle < 0 ) THROW_HARD('negative candidate count in new_fixed_candidate_store')
        nptcls = size(pinds)
        ntot   = int(nptcls,int64) * int(ncandidates_per_particle,int64)
        if( ntot > int(huge(1),int64) ) THROW_HARD('candidate store exceeds default-integer array bounds')
        allocate(self%pinds(nptcls), source=pinds)
        allocate(self%offsets(nptcls+1))
        allocate(self%candidates(int(ntot)))
        allocate(self%seed_shifts(2,nptcls), source=0.)
        allocate(self%seed_has_sh(nptcls), source=.false.)
        do i = 1,nptcls+1
            self%offsets(i) = 1_int64 + int(i-1,int64) * int(ncandidates_per_particle,int64)
        enddo
    end subroutine new_fixed_candidate_store

    subroutine new_ragged_candidate_store( self, pinds, candidate_counts )
        class(prob_candidate_store), intent(inout) :: self
        integer,                     intent(in)    :: pinds(:), candidate_counts(:)
        integer :: i, nptcls
        integer(int64) :: ntot
        call self%kill
        nptcls = size(pinds)
        if( size(candidate_counts) /= nptcls ) THROW_HARD('candidate-count size mismatch in new_ragged_candidate_store')
        if( any(candidate_counts < 0) ) THROW_HARD('negative candidate count in new_ragged_candidate_store')
        allocate(self%pinds(nptcls), source=pinds)
        allocate(self%offsets(nptcls+1))
        self%offsets(1) = 1_int64
        do i = 1,nptcls
            self%offsets(i+1) = self%offsets(i) + int(candidate_counts(i),int64)
        enddo
        ntot = self%offsets(nptcls+1) - 1_int64
        if( ntot > int(huge(1),int64) ) THROW_HARD('candidate store exceeds default-integer array bounds')
        allocate(self%candidates(int(ntot)))
        allocate(self%seed_shifts(2,nptcls), source=0.)
        allocate(self%seed_has_sh(nptcls), source=.false.)
    end subroutine new_ragged_candidate_store

    subroutine kill_candidate_store( self )
        class(prob_candidate_store), intent(inout) :: self
        if( allocated(self%pinds)        ) deallocate(self%pinds)
        if( allocated(self%offsets)      ) deallocate(self%offsets)
        if( allocated(self%candidates)   ) deallocate(self%candidates)
        if( allocated(self%seed_shifts)  ) deallocate(self%seed_shifts)
        if( allocated(self%seed_has_sh)  ) deallocate(self%seed_has_sh)
    end subroutine kill_candidate_store

    !> The shift recorded for an assigned candidate that was not refined: the particle-frame seed every candidate of
    !! the particle was scored with, as it is (every probability table seeds in the particle frame).
    subroutine record_seed_shift( assgn, seed_shift, seed_has_sh, l_doshift )
        type(ptcl_ref), intent(inout) :: assgn
        real,           intent(in)    :: seed_shift(2)
        logical,        intent(in)    :: seed_has_sh, l_doshift
        if( .not. l_doshift .or. assgn%has_sh .or. .not. seed_has_sh ) return
        assgn%x      = seed_shift(1)
        assgn%y      = seed_shift(2)
        assgn%has_sh = .true.
    end subroutine record_seed_shift

    subroutine calc_num2sample( os, num_all, field_str, num_smpl, prob_athres, state)
        class(oris),       intent(in)  :: os
        integer,           intent(in)  :: num_all
        character(len=*),  intent(in)  :: field_str
        integer,           intent(out) :: num_smpl
        real,              intent(in)  :: prob_athres
        integer, optional, intent(in)  :: state
        real :: athres
        athres   = calc_athres(os, field_str, prob_athres, state=state)
        num_smpl = min(num_all,max(1,int(athres * real(num_all) / 180.)))
    end subroutine calc_num2sample

    !> The width of the shift prior (shift_prior=yes): the variance per axis of the shift increments the sampled
    !! particles made in their last search (shincarg, the increment's length, zero when the search kept the shift
    !! and zero for a particle never searched, which a particle record cannot tell apart: for a 2D Gaussian the
    !! mean square length is twice the variance per axis), floored at the square of max(SHIFT_PRIOR_SIGMA_MIN,
    !! SHIFT_PRIOR_LP_FRAC x lp_pix), lp_pix the matching wavelength in the shifts' pixels: lp / smpd, the
    !! original sampling, not lp / smpd_crop (the crop is tied to lp, so in its pixels the band term is a fixed
    !! fraction of a pixel and the floor collapses to the one-pixel guard). A fresh field (no increments yet, or
    !! only the zero increments of iterations that searched no shifts) sits at the floor, so the prior holds a
    !! search within a quarter wavelength or a pixel of its seed
    real function shift_prior_sigma2( os, lp_pix ) result( sigma2 )
        class(oris), intent(in) :: os
        real,        intent(in) :: lp_pix
        real, allocatable :: vals(:)
        real :: sigma_min
        sigma_min = max(SHIFT_PRIOR_SIGMA_MIN, SHIFT_PRIOR_LP_FRAC * lp_pix)
        sigma2    = sigma_min**2
        vals      = os%get_all_sampled('shincarg')
        if( size(vals) > 0 ) sigma2 = max(sigma2, real(0.5d0 * sum(real(vals,dp)**2) / real(size(vals),dp)))
    end function shift_prior_sigma2

    function calc_athres( os, field_str, prob_athres, state ) result( athres )
        class(oris),       intent(in) :: os
        character(len=*),  intent(in) :: field_str
        real,              intent(in) :: prob_athres
        integer, optional, intent(in) :: state
        real, allocatable :: vals(:)
        real :: athres, dist_thres
        vals = os%get_all_sampled(trim(field_str), state=state)
        dist_thres = sum(vals) / real(size(vals))
        athres     = prob_athres
        if( dist_thres > TINY ) athres = min(athres, dist_thres)
    end function calc_athres

    !> The temperature of the pose posterior over a particle's table row. The row holds the normalized
    !! Euclidean loss L_r = sum_k (k/sigma2_k) sum_p |X_kp - CTF R_kp|^2 / wsqsum(X), summed over the pftsz
    !! polar samples of each ring (gen_raw_euclid_vals). The Gaussian negative log-likelihood runs over the
    !! Cartesian half-plane components, pi*k of them in ring k, so each polar sample stands for pi*k/pftsz
    !! of them, and sigma2 is the variance per real or imaginary part (gen_sigma_contrib):
    !!   NLL_r = sum_k (pi k/pftsz) sum_p |r_kp|^2 / (2 sigma2_k) = pi wsqsum(X) L_r / (2 pftsz).
    !! The posterior is p_r proportional to exp(-T (L_r - L_min)) with T = pi wsqsum(X) / (2 pftsz).
    !! wsqsum is the calculator's memoized weighted power of the particle (wsqsum_ptcl_now).
    pure real function posterior_temperature( wsqsum, pftsz ) result( temp )
        real(dp), intent(in) :: wsqsum
        integer,  intent(in) :: pftsz
        temp = real(DPI * wsqsum / (2.d0 * real(pftsz,dp)))
    end function posterior_temperature

    !> The tempering factor of iteration iter of a run of niters (the tables' posterior draws): a geometric ramp from
    !! beta0 at the first iteration to beta1 at the last, so every iteration cools by the same factor (0.1 to 1
    !! over 25 iterations is x1.1 per iteration). beta0 = beta1 is a flat
    !! run; a run of one iteration is at beta1, the factor the map is weighted at
    pure real function prob_beta_at( iter, niters, beta0, beta1 ) result( beta )
        integer, intent(in) :: iter, niters
        real,    intent(in) :: beta0, beta1
        real :: frac
        if( niters <= 1 )then
            beta = beta1
            return
        endif
        frac = real(max(1, min(iter, niters)) - 1) / real(niters - 1)
        beta = beta0 * (beta1 / beta0)**frac
    end function prob_beta_at

    ! ---- the n_eff readout -----------------------------------------------------------------------------------

    !> nptcls rows will be offered; keep up to nmax of them, each with at least nmin valid distances
    !> nptcls is the global sample's size (the stride's), rank0 the rank in it of the row before this table's first
    !! (a part's offset; 0 for the owner's whole sample)
    subroutine neff_readout_new( self, nptcls, nmax, nmin, rank0 )
        class(prob_neff_readout), intent(inout) :: self
        integer,                  intent(in)    :: nptcls, nmax, nmin
        integer,        optional, intent(in)    :: rank0
        call self%kill
        self%nmax   = max(1, nmax)
        self%nmin   = max(2, nmin)
        self%stride = max(1, nptcls / self%nmax)
        self%rank0  = 0
        if( present(rank0) ) self%rank0 = rank0
        self%nrows  = 0
        self%nvals  = 0
        allocate(self%offs(self%nmax + 1), self%temps(self%nmax), self%dists(64 * self%nmax))
        self%offs(1) = 1
    end subroutine neff_readout_new

    !> is the row of the table's i-th particle wanted (its global rank's turn in the stride, within the first nmax
    !! strides)? the caller skips building the rest
    pure logical function neff_readout_wanted( self, i ) result( l_wanted )
        class(prob_neff_readout), intent(in) :: self
        integer,                  intent(in) :: i
        integer :: rank
        rank     = self%rank0 + i
        l_wanted = mod(rank - 1, self%stride) == 0 .and. (rank - 1) / self%stride < self%nmax
    end function neff_readout_wanted

    !> offer a row: its valid distances (finite, below the invalid marker) are packed above their minimum when
    !! there are at least nmin of them and the temperature is positive
    subroutine neff_readout_offer( self, n, dists, temp )
        class(prob_neff_readout), intent(inout) :: self
        integer,                  intent(in)    :: n
        real,                     intent(in)    :: dists(n), temp
        real, allocatable :: tmp(:)
        real    :: dmin, invalid
        integer :: r, nvalid, need
        if( .not. (temp > 0.) ) return
        invalid = 0.1 * huge(invalid)
        nvalid  = 0
        dmin    = huge(dmin)
        do r = 1,n
            if( ieee_is_finite(dists(r)) .and. dists(r) < invalid )then
                nvalid = nvalid + 1
                dmin   = min(dmin, dists(r))
            endif
        enddo
        if( nvalid < self%nmin ) return
        need = self%nvals + nvalid
        if( need > size(self%dists) )then
            allocate(tmp(max(need, 2 * size(self%dists))))
            tmp(1:self%nvals) = self%dists(1:self%nvals)
            call move_alloc(tmp, self%dists)
        endif
        do r = 1,n
            if( ieee_is_finite(dists(r)) .and. dists(r) < invalid )then
                self%nvals             = self%nvals + 1
                self%dists(self%nvals) = dists(r) - dmin
            endif
        enddo
        self%nrows             = self%nrows + 1
        self%temps(self%nrows) = temp
        self%offs(self%nrows + 1) = self%nvals + 1
    end subroutine neff_readout_offer

    pure integer function neff_readout_nkept( self ) result( nrows )
        class(prob_neff_readout), intent(in) :: self
        nrows = self%nrows
    end function neff_readout_nkept

    !> the mean over the kept rows of exp(entropy) of the row posterior at tempering factor beta
    real function neff_readout_mean_neff( self, beta ) result( neff )
        class(prob_neff_readout), intent(in) :: self
        real,                     intent(in) :: beta
        neff = 0.
        if( self%nrows > 0 ) neff = real(self%sum_neff(beta) / real(self%nrows,dp))
    end function neff_readout_mean_neff

    !> the sum over the kept rows of exp(entropy) of the row posterior at tempering factor beta:
    !! p_r proportional to exp(-beta T d_r), d_r the distance above the row minimum; in one pass per row,
    !! H = log z + beta T <d>_p with z the partition sum
    real(dp) function neff_readout_sum_neff( self, beta ) result( acc )
        class(prob_neff_readout), intent(in) :: self
        real,                     intent(in) :: beta
        real(dp) :: z, s, a, e
        integer  :: irow, k
        acc = 0.d0
        !$omp parallel do default(shared) private(irow,k,z,s,a,e) schedule(static) reduction(+:acc)
        do irow = 1,self%nrows
            z = 0.d0
            s = 0.d0
            do k = self%offs(irow), self%offs(irow + 1) - 1
                a = real(beta,dp) * real(self%temps(irow),dp) * real(self%dists(k),dp)
                if( a < 80.d0 )then
                    e = exp(-a)
                    z = z + e
                    s = s + a * e
                endif
            enddo
            acc = acc + exp(log(z) + s / z)
        enddo
        !$omp end parallel do
    end function neff_readout_sum_neff

    !> the readout of one assignment: the mean n_eff at the draw's factor beta_draw and at beta 1 (the physical
    !! posterior), over the kept rows; nothing below PROB_NEFF_MIN_ROWS rows
    subroutine neff_readout_report( self, beta_draw )
        class(prob_neff_readout), intent(in) :: self
        real,                     intent(in) :: beta_draw
        call neff_report_sums(self%nrows, self%sum_neff(beta_draw), self%sum_neff(1.), beta_draw)
    end subroutine neff_readout_report

    !> the readout line from the kept rows and their sums (one table's, or the parts' summed in part order)
    subroutine neff_report_sums( nrows, sum_draw, sum_phys, beta_draw )
        integer,  intent(in) :: nrows
        real(dp), intent(in) :: sum_draw, sum_phys
        real,     intent(in) :: beta_draw
        character(len=32) :: draw_neff, beta, physical_neff
        if( nrows < PROB_NEFF_MIN_ROWS ) return
        write(draw_neff,'(F32.2)') real(sum_draw / real(nrows,dp))
        write(beta,'(F32.2)') beta_draw
        write(physical_neff,'(F32.2)') real(sum_phys / real(nrows,dp))
        write(logfhandle,'(A,I0,A)') '>>> POSTERIOR DRAW N_EFF: '//trim(adjustl(draw_neff))//&
            &' at beta '//trim(adjustl(beta))//', '//trim(adjustl(physical_neff))//' at beta 1.00 (', nrows, ' rows)'
        call flush(logfhandle)
    end subroutine neff_report_sums

    subroutine neff_readout_kill( self )
        class(prob_neff_readout), intent(inout) :: self
        if( allocated(self%offs)  ) deallocate(self%offs)
        if( allocated(self%dists) ) deallocate(self%dists)
        if( allocated(self%temps) ) deallocate(self%temps)
        self%nrows = 0
        self%nvals = 0
    end subroutine neff_readout_kill

    ! ---- the weighted insertion's table ------------------------------------------------------------------------

    !> nptcls particles (the global table's), every one at its committed pose with weight 1 until add says otherwise
    subroutine insert_table_new( self, nptcls, pinds, kmax )
        class(prob_insert_table), intent(inout) :: self
        integer,                  intent(in)    :: nptcls, pinds(nptcls), kmax
        call self%kill
        self%nptcls = nptcls
        self%kmax   = max(1, kmax)
        self%npose  = 0
        allocate(self%pinds(nptcls), source=pinds)
        allocate(self%nextra(nptcls), self%first(nptcls), source=0)
        allocate(self%wc(nptcls), source=1.)
        allocate(self%poses(max(1, nptcls * (self%kmax - 1))))
    end subroutine insert_table_new

    !> particle i's weights: the committed pose's wc and its nextra extra poses, the candidates cands with weights w
    !! (reference, in-plane index and refined shift from the candidate), nelig the eligible entries its extras were
    !! selected among (the cone's occupancy); a particle's extras are added once
    subroutine insert_table_add( self, i, wc, nextra, cands, w, nelig )
        class(prob_insert_table), intent(inout) :: self
        integer,                  intent(in)    :: i, nextra, nelig
        real,                     intent(in)    :: wc
        type(prob_candidate),     intent(in)    :: cands(:)
        real,                     intent(in)    :: w(:)
        type(prob_insert_pose), allocatable :: tmp(:)
        integer :: k
        if( i < 1 .or. i > self%nptcls ) THROW_HARD('particle index out of range; insert_table_add')
        if( self%nextra(i) > 0 ) THROW_HARD('a particle''s extra poses are added once; insert_table_add')
        self%wc(i)     = wc
        self%nelig_sum = self%nelig_sum + int(max(0, nelig),int64)
        if( nelig > self%kmax - 1 ) self%ncapped = self%ncapped + 1
        if( nextra < 1 ) return
        if( self%npose + nextra > size(self%poses) )then
            allocate(tmp(2 * size(self%poses) + nextra))
            tmp(1:self%npose) = self%poses(1:self%npose)
            call move_alloc(tmp, self%poses)
        endif
        self%first(i)  = self%npose + 1
        self%nextra(i) = nextra
        do k = 1,nextra
            self%npose = self%npose + 1
            self%poses(self%npose)%iref = cands(k)%iref
            self%poses(self%npose)%inpl = cands(k)%inpl
            self%poses(self%npose)%x    = cands(k)%x
            self%poses(self%npose)%y    = cands(k)%y
            self%poses(self%npose)%w    = w(k)
        enddo
    end subroutine insert_table_add

    !> the log line: insertions per particle, the committed pose's mean weight, and the cone's occupancy (the
    !! eligible entries per particle and the fraction of particles with more of them than the cap admits: at a
    !! fraction near zero the insertion is the full average within the cone, above it topk is binding)
    subroutine insert_table_report( self )
        class(prob_insert_table), intent(in) :: self
        integer        :: nptcls, npose, ncapped, kmax
        integer(int64) :: nelig
        real(dp)       :: wc_sum
        call self%counts(nptcls, npose, wc_sum, nelig, ncapped, kmax)
        call insert_report_counts(nptcls, npose, wc_sum, nelig, ncapped, kmax)
    end subroutine insert_table_report

    !> the counts the report is made of (a part publishes them, the master sums them in part order); the committed
    !! weights' sum is the single-precision sum the report always took
    subroutine insert_table_counts( self, nptcls, npose, wc_sum, nelig, ncapped, kmax )
        class(prob_insert_table), intent(in)  :: self
        integer,                  intent(out) :: nptcls, npose, ncapped, kmax
        integer(int64),           intent(out) :: nelig
        real(dp),                 intent(out) :: wc_sum
        nptcls  = self%nptcls
        npose   = self%npose
        ncapped = self%ncapped
        kmax    = self%kmax
        nelig   = self%nelig_sum
        wc_sum  = 0.d0
        if( self%nptcls > 0 ) wc_sum = real(sum(self%wc(1:self%nptcls)),dp)
    end subroutine insert_table_counts

    !> the fill's cost record in a matcher's bench report: the counts of the fill (references scored and refined,
    !! the shift searches, the maps, the optimiser's evaluations, the memos) and the thread-seconds of its regions;
    !! the table workers' bench families carried them, the matcher's part record carries them now
    subroutine write_fill_profile_bench( fnr, label, fill, nstates )
        integer,                 intent(in) :: fnr
        character(len=*),        intent(in) :: label
        type(prob_fill_profile), intent(in) :: fill
        integer,                 intent(in) :: nstates
        real(dp) :: ts_sum
        if( fill%nptcls < 1 ) return
        write(fnr,'(a,i0)') key('fill particles'), fill%nptcls
        if( fill%nscored > 0 )then
            ! a neighbourhood table (and the 2D prob_neigh table) scores a different set of references per particle
            write(fnr,'(a,i0)') key('fill scored references (all)'), fill%nscored
            if( fill%nrefined > 0 ) write(fnr,'(a,i0)') key('fill refined refs (all)'), fill%nrefined
        else
            write(fnr,'(a,i0)') key('fill references'), fill%nrefs
            write(fnr,'(a,i0)') key('fill refined refs per particle'), fill%projs_ns
            write(fnr,'(a,i0)') key('fill refined refs (all)'), fill%nrefined
            if( fill%nrefined > 0 )then
                ! where the entry lowest after refinement sat in the unrefined ranking: is the budget binding?
                write(fnr,'(a,i0)')        key('fill refined argmin at unrefined rank 1'), fill%rank_counts(0)
                write(fnr,'(a,10(1x,i0))') key('fill refined argmin by decile of the cap'), fill%rank_counts(1:10)
                write(fnr,'(a,i0)')        key('fill refined argmin outside refined set'), fill%rank_counts(11)
            endif
        endif
        write(fnr,'(a,i0)') key('fill shift searches'), fill%minimize_calls
        if( fill%screens > 0 )    write(fnr,'(a,i0)') key('fill coarse maps (all)'), fill%screens
        if( fill%disk_ptcls > 0 ) write(fnr,'(a,i0)') key('fill particles on the disk'), fill%disk_ptcls
        if( sum(fill%map_counts(1:3)) > 0 )then
            write(fnr,'(a,i0)') key('fill disk maps'), fill%map_counts(1)
            write(fnr,'(a,i0)') key('fill window maps'), fill%map_counts(2)
            write(fnr,'(a,i0)') key('fill full-band maps'), fill%map_counts(3)
            write(fnr,'(a,i0)') key('fill map selections'), fill%map_counts(4)
            write(fnr,'(a,i0)') key('fill map widths'), fill%map_counts(5)
            write(fnr,'(a,i0)') key('fill map particle expansions'), fill%map_counts(6)
            write(fnr,'(a,i0)') key('fill map geometry builds'), fill%map_counts(7)
            write(fnr,'(a,i0)') key('fill map geometry reuses'), fill%map_counts(8)
            write(fnr,'(a,3(1x,i0))') key('fill disk grid radii dirs order'), fill%map_geometry(:,2)
            write(fnr,'(a,3(1x,i0))') key('fill window grid radii dirs order'), fill%map_geometry(:,3)
        endif
        write(fnr,'(a,i0)') key('fill alternations'), fill%alternations
        write(fnr,'(a,i0)') key('fill grid evaluations'), fill%grid_points
        write(fnr,'(a,i0)') key('fill objective evaluations'), fill%objective_evals
        write(fnr,'(a,i0)') key('fill gradient evaluations (E)'), fill%gradient_evals
        if( nstates > 1 )then
            write(fnr,'(a,i0)') key('fill grid fits kept'), fill%fit_counts(0)
            write(fnr,'(a,i0)') key('fill fit fallback, edge'), fill%fit_counts(1)
            write(fnr,'(a,i0)') key('fill fit fallback, curvature'), fill%fit_counts(2)
            write(fnr,'(a,i0)') key('fill fit fallback, step'), fill%fit_counts(3)
            write(fnr,'(a,i0)') key('fill fit fallback, search box'), fill%fit_counts(4)
            write(fnr,'(a,i0)') key('fill fit fallback, value'), fill%fit_counts(5)
        endif
        write(fnr,'(a,i0)') key('fill shifted particle memos'), fill%bank_builds
        if( fill%minimize_calls > 0 )then
            write(fnr,'(a,1x,f0.2)') key('fill alternations per search'), real(fill%alternations,dp) / real(fill%minimize_calls,dp)
            write(fnr,'(a,1x,f0.2)') key('fill gradient evals per search'), real(fill%gradient_evals,dp) / real(fill%minimize_calls,dp)
        endif
        write(fnr,'(a,1x,f0.2)') key('fill map setup (within fill, s)'), fill%map_times(6)
        ! thread-seconds are never part of a report's '% accounted for'
        ts_sum = fill%ts_seed + fill%ts_score + fill%ts_refine
        write(fnr,'(a,1x,f0.2)') key('fill seed estimation (thread-s)'), fill%ts_seed
        write(fnr,'(a,1x,f0.2)') key('fill score loop (thread-s)'), fill%ts_score
        write(fnr,'(a,1x,f0.2)') key('fill refinement (thread-s)'), fill%ts_refine
        if( sum(fill%map_counts(1:3)) > 0 )then
            write(fnr,'(a,1x,f0.4)') key('fill map construction (thread-s)'), fill%map_times(1)
            write(fnr,'(a,1x,f0.4)') key('fill map cell selection (thread-s)'), fill%map_times(4)
            write(fnr,'(a,1x,f0.4)') key('fill map posterior reduction (thread-s)'), fill%map_times(5)
        endif
        write(fnr,'(a,1x,f0.2)') key('fill thread-seconds, sum'), ts_sum
        write(fnr,'(a,1x,f0.2)') key('fill region wall time (s)'), fill%wall_region
        if( fill%wall_region > 0.d0 .and. fill%nthr > 0 )then
            write(fnr,'(a,1x,f0.3)') key('fill load balance'), ts_sum / (fill%wall_region * real(fill%nthr,dp))
        endif
        if( fill%nscored > 0 )then
            write(fnr,'(a,1x,es10.3)') key('fill score loop s/(particle*ref)'), fill%ts_score / real(fill%nscored,dp)
        else if( fill%nptcls > 0 .and. fill%nrefs > 0 )then
            write(fnr,'(a,1x,es10.3)') key('fill score loop s/(particle*ref)'), &
                &fill%ts_score / (real(fill%nptcls,dp) * real(fill%nrefs,dp))
        endif

    contains

        !> the line's key: the report's label, the counter's name, padded to the report's column
        function key( name ) result( k )
            character(len=*), intent(in) :: name
            character(len=44) :: k
            k = label//' '//name
            k(43:44) = ': '
        end function key

    end subroutine write_fill_profile_bench

    !> the insertion report from counts
    subroutine insert_report_counts( nptcls, npose, wc_sum, nelig, ncapped, kmax )
        integer,        intent(in) :: nptcls, npose, ncapped, kmax
        integer(int64), intent(in) :: nelig
        real(dp),       intent(in) :: wc_sum
        character(len=32) :: ninsert, wcommit, nelig_str, pct_capped
        if( nptcls < 1 ) return
        write(ninsert,'(F32.2)') 1. + real(npose) / real(nptcls)
        write(wcommit,'(F32.2)') real(wc_sum) / real(nptcls)
        write(nelig_str,'(F32.2)') real(nelig) / real(nptcls)
        write(pct_capped,'(F32.1)') 100. * real(ncapped) / real(nptcls)
        write(logfhandle,'(A,I0,A)') '>>> PROB_TAB_INSERT (topk ', kmax, '): '//trim(adjustl(ninsert))//&
            &' poses per particle, committed weight '//trim(adjustl(wcommit))//', '//trim(adjustl(nelig_str))//&
            &' eligible per particle, '//trim(adjustl(pct_capped))//'% capped'
        call flush(logfhandle)
    end subroutine insert_report_counts

    !> the readout record of part ipart: prob_readout_part<NN>.txt
    function prob_readout_fname( ipart, numlen ) result( fname )
        integer, intent(in) :: ipart, numlen
        type(string) :: fname
        fname = string('prob_readout_part')//int2str_pad(ipart, numlen)//'.txt'
    end function prob_readout_fname

    !> a part's record, one key and value per line (reals to the last digit)
    subroutine write_prob_readout( rec, fname )
        type(prob_readout), intent(in) :: rec
        class(string),      intent(in) :: fname
        integer :: funit, io_stat
        call fopen(funit, fname, status='replace', action='write', iostat=io_stat)
        call fileiochk('write_prob_readout; fopen '//fname%to_char(), io_stat)
        write(funit,'(a,1x,i0)')      'route',         rec%route
        write(funit,'(a,1x,i0)')      'nptcls',        rec%nptcls
        write(funit,'(a,1x,i0)')      'nrefs',         rec%nrefs
        write(funit,'(a,1x,i0)')      'nmode',         rec%nmode
        write(funit,'(a,1x,i0)')      'ntemp_missing', rec%ntemp_missing
        write(funit,'(a,1x,i0)')      'nbad',          rec%nbad
        write(funit,'(a,1x,i0)')      'nflat',         rec%nflat
        write(funit,'(a,1x,i0)')      'nkept',         rec%nkept
        write(funit,'(a,1x,i0)')      'kmax',          rec%kmax
        write(funit,'(a,1x,l1)')      'l_extras',      rec%l_extras
        write(funit,'(a,1x,i0)')      'ins_nptcls',    rec%ins_nptcls
        write(funit,'(a,1x,i0)')      'ins_npose',     rec%ins_npose
        write(funit,'(a,1x,i0)')      'ins_ncapped',   rec%ins_ncapped
        write(funit,'(a,1x,i0)')      'ins_nelig',     rec%ins_nelig
        write(funit,'(a,1x,es24.16)') 'ins_wc_sum',    rec%ins_wc_sum
        write(funit,'(a,1x,i0)')      'neff_nrows',    rec%neff_nrows
        write(funit,'(a,1x,es24.16)') 'neff_sum_draw', rec%neff_sum_draw
        write(funit,'(a,1x,es24.16)') 'neff_sum_phys', rec%neff_sum_phys
        call fclose(funit)
    end subroutine write_prob_readout

    subroutine read_prob_readout( rec, fname )
        type(prob_readout), intent(out) :: rec
        class(string),      intent(in)  :: fname
        character(len=64) :: line_value
        integer :: funit, io_stat
        if( .not. file_exists(fname) ) THROW_HARD('a part published no readout record: '//fname%to_char())
        call fopen(funit, fname, status='old', action='read', iostat=io_stat)
        call fileiochk('read_prob_readout; fopen '//fname%to_char(), io_stat)
        call read_int('route',         rec%route)
        call read_int('nptcls',        rec%nptcls)
        call read_int('nrefs',         rec%nrefs)
        call read_int('nmode',         rec%nmode)
        call read_int('ntemp_missing', rec%ntemp_missing)
        call read_int('nbad',          rec%nbad)
        call read_int('nflat',         rec%nflat)
        call read_int('nkept',         rec%nkept)
        call read_int('kmax',          rec%kmax)
        call read_log('l_extras',      rec%l_extras)
        call read_int('ins_nptcls',    rec%ins_nptcls)
        call read_int('ins_npose',     rec%ins_npose)
        call read_int('ins_ncapped',   rec%ins_ncapped)
        call read_i8 ('ins_nelig',     rec%ins_nelig)
        call read_dp ('ins_wc_sum',    rec%ins_wc_sum)
        call read_int('neff_nrows',    rec%neff_nrows)
        call read_dp ('neff_sum_draw', rec%neff_sum_draw)
        call read_dp ('neff_sum_phys', rec%neff_sum_phys)
        call fclose(funit)

        contains

            !> one line: its key checked against the one expected (a stale or foreign file is refused)
            subroutine read_key( expected )
                character(len=*), intent(in) :: expected
                character(len=32) :: key
                read(funit,*,iostat=io_stat) key, line_value
                if( io_stat /= 0 ) THROW_HARD('readout record '//fname%to_char()//' ends before '//trim(expected))
                if( trim(key) /= trim(expected) ) THROW_HARD('readout record '//fname%to_char()//' holds '//trim(key)//' where '//trim(expected)//' is expected')
            end subroutine read_key

            subroutine read_int( expected, val )
                character(len=*), intent(in)  :: expected
                integer,          intent(out) :: val
                call read_key(expected)
                read(line_value,*) val
            end subroutine read_int

            subroutine read_i8( expected, val )
                character(len=*), intent(in)  :: expected
                integer(int64),   intent(out) :: val
                call read_key(expected)
                read(line_value,*) val
            end subroutine read_i8

            subroutine read_dp( expected, val )
                character(len=*), intent(in)  :: expected
                real(dp),         intent(out) :: val
                call read_key(expected)
                read(line_value,*) val
            end subroutine read_dp

            subroutine read_log( expected, val )
                character(len=*), intent(in)  :: expected
                logical,          intent(out) :: val
                call read_key(expected)
                read(line_value,*) val
            end subroutine read_log

    end subroutine read_prob_readout

    !> the parts' records summed in part order; the route, the reference count and topk are the same on every part,
    !! the extras flag is any part's (a part without a sampled particle has none)
    subroutine add_prob_readout( acc, rec )
        type(prob_readout), intent(inout) :: acc
        type(prob_readout), intent(in)    :: rec
        if( rec%route > 0 ) acc%route = rec%route
        acc%nptcls        = acc%nptcls        + rec%nptcls
        acc%nrefs         = max(acc%nrefs, rec%nrefs)
        acc%nmode         = acc%nmode         + rec%nmode
        acc%ntemp_missing = acc%ntemp_missing + rec%ntemp_missing
        acc%nbad          = acc%nbad          + rec%nbad
        acc%nflat         = acc%nflat         + rec%nflat
        acc%nkept         = acc%nkept         + rec%nkept
        acc%kmax          = max(acc%kmax, rec%kmax)
        acc%l_extras      = acc%l_extras .or. rec%l_extras
        acc%ins_nptcls    = acc%ins_nptcls    + rec%ins_nptcls
        acc%ins_npose     = acc%ins_npose     + rec%ins_npose
        acc%ins_ncapped   = acc%ins_ncapped   + rec%ins_ncapped
        acc%ins_nelig     = acc%ins_nelig     + rec%ins_nelig
        acc%ins_wc_sum    = acc%ins_wc_sum    + rec%ins_wc_sum
        acc%neff_nrows    = acc%neff_nrows    + rec%neff_nrows
        acc%neff_sum_draw = acc%neff_sum_draw + rec%neff_sum_draw
        acc%neff_sum_phys = acc%neff_sum_phys + rec%neff_sum_phys
    end subroutine add_prob_readout

    !> the distributed master's merge of its parts' records in part order, reported as the owner reports its own
    !! (plan section 5: sums and counts, never means of means); the records are deleted after the merge
    subroutine merge_prob_readouts( nparts, numlen, prob_beta )
        integer, intent(in) :: nparts, numlen
        real,    intent(in) :: prob_beta
        type(prob_readout) :: acc, rec
        type(string)       :: fname
        integer :: ipart
        do ipart = 1,nparts
            fname = prob_readout_fname(ipart, numlen)
            call read_prob_readout(rec, fname)
            call add_prob_readout(acc, rec)
            call del_file(fname)
        enddo
        call report_prob_readout(acc, prob_beta)
    end subroutine merge_prob_readouts

    !> the lines an assignment prints, from a record (the owner's own, or the parts' sum): the n_eff readout, the
    !! insert report or the "no FSC yet" line (3D), the row-mode share, the warnings; the dense 3D stops when every
    !! row is invalid or flat
    subroutine report_prob_readout( rec, prob_beta )
        type(prob_readout), intent(in) :: rec
        real,               intent(in) :: prob_beta
        character(len=32) :: mode_pct
        integer :: nassigned
        call neff_report_sums(rec%neff_nrows, rec%neff_sum_draw, rec%neff_sum_phys, prob_beta)
        select case(rec%route)
            case(PROB_READOUT_3D_DENSE, PROB_READOUT_3D_NEIGH)
                if( rec%l_extras )then
                    call insert_report_counts(rec%ins_nptcls, rec%ins_npose, rec%ins_wc_sum, rec%ins_nelig, &
                        &rec%ins_ncapped, rec%kmax)
                else if( rec%kmax > 1 )then
                    write(logfhandle,'(A,I0,A)') '>>> PROB_TAB_INSERT (topk ', rec%kmax, '): no FSC yet, committed pose only'
                endif
                nassigned = rec%nptcls
                if( rec%route == PROB_READOUT_3D_NEIGH ) nassigned = rec%nptcls - rec%nkept
                write(mode_pct,'(F32.1)') 100. * real(rec%nmode) / real(max(1,nassigned))
                write(logfhandle,'(A,I0,A,I0,A)') '>>> PROB_TAB_ASSIGN: ', nassigned, ' particles x ', rec%nrefs,&
                    &' references, '//trim(adjustl(mode_pct))//'% at the row mode'
                if( rec%route == PROB_READOUT_3D_NEIGH )then
                    if( rec%nkept > 0 ) write(logfhandle,'(A,I0)') &
                        &'WARNING: particles without a valid neighbourhood candidate, orientation kept: ', rec%nkept
                    if( rec%ntemp_missing > 0 ) write(logfhandle,'(A,I0)') &
                        &'WARNING: particles without a temperature in the neighbourhood rows, assigned to their argmin: ', &
                        &rec%ntemp_missing
                else
                    if( rec%ntemp_missing > 0 ) write(logfhandle,'(A,I0)') &
                        &'WARNING: particles without a temperature in the dense rows, assigned to their argmin: ', &
                        &rec%ntemp_missing
                    if( rec%nptcls > 0 )then
                        if( rec%nbad  == rec%nptcls ) THROW_HARD('all probability-table reference distances are invalid')
                        if( rec%nflat == rec%nptcls ) THROW_HARD('all probability-table reference distances are flat')
                    endif
                    if( rec%nbad  > 0 ) write(logfhandle,*) 'WARNING: particles with invalid probability-table distances: ', rec%nbad
                    if( rec%nflat > 0 ) write(logfhandle,*) 'WARNING: particles with flat probability-table reference distances: ', rec%nflat
                endif
            case(PROB_READOUT_2D_DENSE, PROB_READOUT_2D_SPARSE)
                write(mode_pct,'(F32.1)') 100. * real(rec%nmode) / real(max(1,rec%nptcls))
                write(logfhandle,'(A,I0,A,I0,A)') '>>> PROB_TAB2D_ASSIGN: ', rec%nptcls, ' particles x ', rec%nrefs,&
                    &' classes, '//trim(adjustl(mode_pct))//'% at the row mode'
                if( rec%ntemp_missing > 0 )then
                    if( rec%route == PROB_READOUT_2D_SPARSE )then
                        write(logfhandle,'(A,I0)') &
                            &'WARNING: particles without a temperature in the 2D candidate rows, assigned to their best candidate: ', &
                            &rec%ntemp_missing
                    else
                        write(logfhandle,'(A,I0)') &
                            &'WARNING: particles without a temperature in the 2D rows, assigned to their argmin: ', rec%ntemp_missing
                    endif
                endif
            case DEFAULT
                THROW_HARD('unknown readout route; report_prob_readout')
        end select
        call flush(logfhandle)
    end subroutine report_prob_readout

    !> the eligible entries per particle, over the particles added
    pure real function insert_table_get_mean_eligible( self )
        class(prob_insert_table), intent(in) :: self
        insert_table_get_mean_eligible = 0.
        if( self%nptcls > 0 ) insert_table_get_mean_eligible = real(self%nelig_sum) / real(self%nptcls)
    end function insert_table_get_mean_eligible

    !> the fraction of particles with more eligible entries than kmax - 1
    pure real function insert_table_get_frac_capped( self )
        class(prob_insert_table), intent(in) :: self
        insert_table_get_frac_capped = 0.
        if( self%nptcls > 0 ) insert_table_get_frac_capped = real(self%ncapped) / real(self%nptcls)
    end function insert_table_get_frac_capped

    pure integer function insert_table_get_nptcls( self )
        class(prob_insert_table), intent(in) :: self
        insert_table_get_nptcls = self%nptcls
    end function insert_table_get_nptcls

    pure integer function insert_table_get_kmax( self )
        class(prob_insert_table), intent(in) :: self
        insert_table_get_kmax = self%kmax
    end function insert_table_get_kmax

    pure integer function insert_table_get_pind( self, i )
        class(prob_insert_table), intent(in) :: self
        integer,                  intent(in) :: i
        insert_table_get_pind = self%pinds(i)
    end function insert_table_get_pind

    pure real function insert_table_get_committed_weight( self, i )
        class(prob_insert_table), intent(in) :: self
        integer,                  intent(in) :: i
        insert_table_get_committed_weight = self%wc(i)
    end function insert_table_get_committed_weight

    pure integer function insert_table_get_nextra( self, i )
        class(prob_insert_table), intent(in) :: self
        integer,                  intent(in) :: i
        insert_table_get_nextra = self%nextra(i)
    end function insert_table_get_nextra

    !> particle i's k-th extra pose
    pure function insert_table_get_extra( self, i, k ) result( pose )
        class(prob_insert_table), intent(in) :: self
        integer,                  intent(in) :: i, k
        type(prob_insert_pose) :: pose
        pose = self%poses(self%first(i) + k - 1)
    end function insert_table_get_extra

    subroutine insert_table_kill( self )
        class(prob_insert_table), intent(inout) :: self
        if( allocated(self%pinds)  ) deallocate(self%pinds)
        if( allocated(self%nextra) ) deallocate(self%nextra)
        if( allocated(self%first)  ) deallocate(self%first)
        if( allocated(self%wc)     ) deallocate(self%wc)
        if( allocated(self%poses)  ) deallocate(self%poses)
        self%nptcls    = 0
        self%kmax      = 1
        self%npose     = 0
        self%ncapped   = 0
        self%nelig_sum = 0_int64
    end subroutine insert_table_kill

    ! ---- the resolved poses of a part -----------------------------------------------------------------------------

    !> the pipeline's layout, a batch at a time in particle order: the first call (i_first = 1) lays the structure out
    !! for the part with room for every particle's kmax - 1 extras, each later call takes the counts and weights of
    !! the batch's particles from the table as their extras are selected, so the slots are those setup gives when the
    !! whole table is known
    subroutine insert_poses_setup_range( self, tab, pinds, i_first, i_last )
        class(prob_insert_poses), intent(inout) :: self
        class(prob_insert_table), intent(in)    :: tab
        integer,                  intent(in)    :: pinds(:), i_first, i_last
        integer :: i, max_pind, capacity
        if( i_first == 1 )then
            call self%kill
            if( tab%nptcls /= size(pinds) ) THROW_HARD('the table is not the part''s; insert_poses_setup_range')
            self%nptcls = tab%nptcls
            self%npose  = 0
            capacity    = max(1, self%nptcls * (tab%kmax - 1))
            max_pind    = 0
            if( size(pinds) > 0 ) max_pind = maxval(pinds)
            allocate(self%loc(max(1,max_pind)), self%nextra(max(1,self%nptcls)), self%first(max(1,self%nptcls)), source=0)
            allocate(self%wc(max(1,self%nptcls)), source=1.)
            allocate(self%euler(3,capacity), self%shift(2,capacity), source=0.)
            allocate(self%w(capacity), source=0.)
            do i = 1,self%nptcls
                if( tab%pinds(i) /= pinds(i) ) THROW_HARD('the table''s particles are not the part''s; insert_poses_setup_range')
                if( pinds(i) > 0 ) self%loc(pinds(i)) = i
            enddo
        endif
        if( .not. allocated(self%loc) ) THROW_HARD('the layout starts with the first batch; insert_poses_setup_range')
        if( i_first < 1 .or. i_last > self%nptcls .or. i_last < i_first ) THROW_HARD('invalid range; insert_poses_setup_range')
        do i = i_first,i_last
            self%wc(i)     = tab%wc(i)
            self%nextra(i) = tab%nextra(i)
            self%first(i)  = self%npose + 1
            self%npose     = self%npose + tab%nextra(i)
        enddo
        if( self%npose > size(self%w) ) THROW_HARD('more extras than the layout holds; insert_poses_setup_range')
    end subroutine insert_poses_setup_range

    !> whether particle iptcl has an entry
    pure logical function insert_poses_has( self, iptcl )
        class(prob_insert_poses), intent(in) :: self
        integer,                  intent(in) :: iptcl
        insert_poses_has = .false.
        if( .not. allocated(self%loc) ) return
        if( iptcl < 1 .or. iptcl > size(self%loc) ) return
        insert_poses_has = self%loc(iptcl) > 0
    end function insert_poses_has

    subroutine insert_poses_kill( self )
        class(prob_insert_poses), intent(inout) :: self
        if( allocated(self%loc)    ) deallocate(self%loc)
        if( allocated(self%wc)     ) deallocate(self%wc)
        if( allocated(self%nextra) ) deallocate(self%nextra)
        if( allocated(self%first)  ) deallocate(self%first)
        if( allocated(self%euler)  ) deallocate(self%euler)
        if( allocated(self%shift)  ) deallocate(self%shift)
        if( allocated(self%w)      ) deallocate(self%w)
        self%nptcls = 0
        self%npose  = 0
    end subroutine insert_poses_kill

    !> One draw from the posterior over a table row at temperature temp (the tables' assignment): p_r proportional to
    !! exp(-temp (dist_r - dist_min)) over the valid distances (finite, below the invalid marker), zero
    !! elsewhere. which is the drawn index and dist its distance; which = 0 when no distance is valid.
    !! temp = 0 draws uniformly over the valid distances; a huge temp returns the argmin when the minimum
    !! is unique and one of the equal minima with equal probability otherwise (a flat row is a uniform draw
    !! at any temperature: the posterior, not the balanced rule's lowest-index convention). w is a work
    !! array of size >= n; u is the draw's unit variate in [0,1), the caller's keyed draw (rnd_keyed), so
    !! the selection depends on the row and the key alone.
    subroutine draw_posterior_row( n, dists, temp, w, u, which, dist )
        integer,  intent(in)    :: n
        real,     intent(in)    :: dists(n), temp
        real,     intent(inout) :: w(:)
        real(dp), intent(in)    :: u
        integer,  intent(out)   :: which
        real,     intent(out)   :: dist
        real(dp) :: wsum, csum, tgt, arg
        real     :: dmin, invalid
        integer  :: r, imin
        invalid = 0.1 * huge(invalid)
        which   = 0
        dist    = huge(dist)
        imin    = 0
        dmin    = huge(dmin)
        do r = 1,n
            if( ieee_is_finite(dists(r)) .and. dists(r) < invalid )then
                if( dists(r) < dmin )then
                    dmin = dists(r)
                    imin = r
                endif
            endif
        enddo
        if( imin == 0 ) return
        wsum = 0.d0
        do r = 1,n
            w(r) = 0.
            if( ieee_is_finite(dists(r)) .and. dists(r) < invalid )then
                ! the mode's weight is exp(0) = 1 whatever temp, so wsum >= 1; weights below exp(-80) are 0
                ! (not exp'd: a single-precision underflow would flag)
                arg = real(max(0.,temp),dp) * real(dists(r) - dmin,dp)
                if( arg < 80.d0 )then
                    w(r) = real(exp(-arg))
                    wsum = wsum + real(w(r),dp)
                endif
            endif
        enddo
        tgt    = u * wsum
        csum   = 0.d0
        do r = 1,n
            if( w(r) <= 0. ) cycle
            csum = csum + real(w(r),dp)
            if( csum >= tgt )then
                which = r
                exit
            endif
        enddo
        if( which == 0 ) which = imin ! rounding at the top of the cumulative sum
        dist = dists(which)
    end subroutine draw_posterior_row

    ! ---- the weighted insertion's cone ------------------------------------------------------------------------

    !> per state the resolution res(s) (FSC 0.5, Angstrom) the rim displacement is bounded by, for a particle of
    !! mask diameter diam whose stored shifts are in pixels of smpd; a state with res <= 0 has no cone, and without
    !! a diameter or a sampling distance no state has
    subroutine insert_cone_new( self, res, diam, smpd )
        class(prob_insert_cone), intent(inout) :: self
        real,                    intent(in)    :: res(:), diam, smpd
        call self%kill
        allocate(self%res(size(res)), source=res)
        if( diam > 0. .and. smpd > 0. )then
            self%diam = diam
            self%smpd = smpd
        else
            self%res = -1.
        endif
    end subroutine insert_cone_new

    !> the geometry of the test, as rotation matrices: the directions' (eulspace, in-plane angle zero), the rotation
    !! grid's (the in-plane index i is the Euler angle 360 - (i - 1) 360 / nrots, as the committed pose is converted)
    !! and the point group's
    subroutine insert_cone_set_space( self, eulspace, se, nrots )
        class(prob_insert_cone), intent(inout) :: self
        class(oris),             intent(in)    :: eulspace
        class(sym),              intent(in)    :: se
        integer,                 intent(in)    :: nrots
        real    :: e(3)
        integer :: i
        if( allocated(self%dmat) ) deallocate(self%dmat)
        if( allocated(self%rz)   ) deallocate(self%rz)
        if( allocated(self%rsym) ) deallocate(self%rsym)
        self%nspace = eulspace%get_noris()
        self%nrots  = nrots
        self%nsym   = se%get_nsym()
        allocate(self%dmat(3,3,max(1,self%nspace)), self%rz(3,3,max(1,self%nrots)), self%rsym(3,3,max(1,self%nsym)), source=0.)
        do i = 1,self%nspace
            e    = eulspace%get_euler(i)
            e(3) = 0.
            self%dmat(:,:,i) = euler2m(e)
        enddo
        do i = 1,self%nrots
            self%rz(:,:,i) = euler2m([0., 0., 360. - real(i - 1) * 360. / real(self%nrots)])
        enddo
        do i = 1,self%nsym
            call se%get_sym_rmat(i, self%rsym(:,:,i))
        enddo
    end subroutine insert_cone_set_space

    !> whether any state has a cone
    pure logical function insert_cone_is_open( self )
        class(prob_insert_cone), intent(in) :: self
        insert_cone_is_open = .false.
        if( allocated(self%res) ) insert_cone_is_open = any(self%res > 0.)
    end function insert_cone_is_open

    pure logical function insert_cone_has_space( self )
        class(prob_insert_cone), intent(in) :: self
        insert_cone_has_space = self%nspace > 0 .and. self%nrots > 0 .and. self%nsym > 0
    end function insert_cone_has_space

    !> the half-angle of state istate's cone at equal shifts, degrees (the rotation whose displacement over the
    !! mask's radius is the resolution, capped at the half-turn; below zero: no cone)
    pure real function insert_cone_get_angle( self, istate )
        class(prob_insert_cone), intent(in) :: self
        integer,                 intent(in) :: istate
        insert_cone_get_angle = -1.
        if( .not. allocated(self%res) ) return
        if( istate < 1 .or. istate > size(self%res) ) return
        if( .not. (self%res(istate) > 0.) ) return
        insert_cone_get_angle = rad2deg(min(2. * self%res(istate) / self%diam, PI))
    end function insert_cone_get_angle

    !> whether the pose (direction proj_r, in-plane index inpl_r, stored-frame shift sh_r) lies within state istate's
    !! cone of the committed pose (proj_c, inpl_c, sh_c): the angle of the rotation between the two poses, the
    !! smallest over the point group's equivalents of the candidate (R_r R_g, as the group is applied to a pose),
    !! times the mask's radius, plus the distance between the shifts, at most the resolution. A rotation by alpha
    !! moves a point at the rim by at most 2 r sin(alpha/2) <= r alpha, a shift moves every point by its length
    logical function insert_cone_within( self, istate, proj_c, inpl_c, sh_c, proj_r, inpl_r, sh_r ) result( l_in )
        class(prob_insert_cone), intent(in) :: self
        integer,                 intent(in) :: istate, proj_c, inpl_c, proj_r, inpl_r
        real,                    intent(in) :: sh_c(2), sh_r(2)
        real    :: rc(3,3), rr(3,3), m(3,3), ds, trmax, tr, alpha
        integer :: g
        l_in = .false.
        if( .not. self%has_space() ) return
        if( istate < 1 .or. istate > size(self%res) ) return
        if( .not. (self%res(istate) > 0.) ) return
        if( proj_c < 1 .or. proj_c > self%nspace .or. proj_r < 1 .or. proj_r > self%nspace ) return
        if( inpl_c < 1 .or. inpl_c > self%nrots  .or. inpl_r < 1 .or. inpl_r > self%nrots  ) return
        ds = self%smpd * hypot(sh_r(1) - sh_c(1), sh_r(2) - sh_c(2))
        if( ds > self%res(istate) ) return
        rc    = matmul(self%rz(:,:,inpl_c), self%dmat(:,:,proj_c))
        rr    = matmul(self%rz(:,:,inpl_r), self%dmat(:,:,proj_r))
        m     = matmul(transpose(rc), rr)
        trmax = -huge(trmax)
        do g = 1,self%nsym
            ! the trace of R_g M: the cosine of the angle between the committed pose and the candidate's equivalent
            tr    = sum(self%rsym(:,:,g) * transpose(m))
            trmax = max(trmax, tr)
        enddo
        alpha = acos(max(-1., min(1., 0.5 * (trmax - 1.))))
        l_in  = 0.5 * self%diam * alpha + ds <= self%res(istate)
    end function insert_cone_within

    subroutine insert_cone_kill( self )
        class(prob_insert_cone), intent(inout) :: self
        if( allocated(self%res)  ) deallocate(self%res)
        if( allocated(self%dmat) ) deallocate(self%dmat)
        if( allocated(self%rz)   ) deallocate(self%rz)
        if( allocated(self%rsym) ) deallocate(self%rsym)
        self%diam   = 0.
        self%smpd   = 0.
        self%nspace = 0
        self%nrots  = 0
        self%nsym   = 0
    end subroutine insert_cone_kill

    !> The extra poses of the weighted insertion (topk > 1): the kmax - 1 eligible entries of a row (mask: the refined
    !! ones within the cone) of highest posterior p_r proportional to exp(-temp (dist_r - dist_min)), other than the
    !! committed entry ifirst, with the weights that posterior renormalised over the committed entry and the extras
    !! (the committed entry's weight is its own share, whatever its rank). temp is the particle's own temperature,
    !! untempered: the committed pose is the search's draw from the tempered row, hot so that the search explores;
    !! the map is weighted by what the data say of each pose, so the reconstruction carries the data's pose
    !! uncertainty and not the search's. An extra whose weight falls below wmin is dropped (the rest renormalised
    !! again), so no pose enters the reconstruction with less than that share of its particle; the committed entry
    !! is always kept, so kmax = 1 keeps it alone at weight 1. The selection is deterministic (the top entries, no
    !! draw), so topk touches no random stream. On return nkept entries, the committed first and the extras in
    !! descending weight, with their row indices inds and weights w summing to 1, and nelig the eligible entries
    !! with a valid distance other than the committed one: the cone's occupancy, which above kmax - 1 the cap cuts
    subroutine top_posterior_extras( n, dists, temp, mask, ifirst, kmax, wmin, nkept, inds, w, nelig )
        integer, intent(in)  :: n, ifirst, kmax
        real,    intent(in)  :: dists(n), temp, wmin
        logical, intent(in)  :: mask(n)
        integer, intent(out) :: nkept, inds(:), nelig
        real,    intent(out) :: w(:)
        real(dp) :: arg, tot
        real     :: dmin, invalid, wr
        integer  :: r, k, m, nsel
        nkept   = 1
        inds(1) = ifirst
        w(1)    = 1.
        nelig   = 0
        invalid = 0.1 * huge(invalid)
        dmin    = huge(dmin)
        ! the local normaliser over the committed entry and the eligible ones: an excluded row never changes their weights
        do r = 1,n
            if( r /= ifirst .and. .not. mask(r) ) cycle
            if( ieee_is_finite(dists(r)) .and. dists(r) < invalid ) dmin = min(dmin, dists(r))
        enddo
        if( dmin >= invalid ) return
        do r = 1,n
            if( r == ifirst .or. .not. mask(r) ) cycle
            if( ieee_is_finite(dists(r)) .and. dists(r) < invalid ) nelig = nelig + 1
        enddo
        if( kmax <= 1 .or. nelig < 1 ) return
        ! the unnormalised posterior of the committed entry
        w(1) = 0.
        if( ieee_is_finite(dists(ifirst)) .and. dists(ifirst) < invalid )then
            arg = real(max(0.,temp),dp) * real(dists(ifirst) - dmin,dp)
            if( arg < 80.d0 ) w(1) = real(exp(-arg))
        endif
        ! the kmax - 1 eligible entries of highest weight, kept in descending order in inds(2:) and w(2:); weights
        ! below exp(-80) are zero (not exp'd: a single-precision underflow would flag) and never selected; on a tie
        ! the earlier entry stays
        nsel = 0
        do r = 1,n
            if( r == ifirst .or. .not. mask(r) ) cycle
            if( .not. (ieee_is_finite(dists(r)) .and. dists(r) < invalid) ) cycle
            arg = real(max(0.,temp),dp) * real(dists(r) - dmin,dp)
            if( .not. (arg < 80.d0) ) cycle
            wr = real(exp(-arg))
            if( nsel == kmax - 1 )then
                if( .not. (wr > w(nsel + 1)) ) cycle
                k = nsel + 1
            else
                nsel = nsel + 1
                k    = nsel + 1
            endif
            inds(k) = r
            w(k)    = wr
            do while( k > 2 )
                if( .not. (w(k) > w(k - 1)) ) exit
                m         = inds(k - 1)
                inds(k-1) = inds(k)
                inds(k)   = m
                wr        = w(k - 1)
                w(k-1)    = w(k)
                w(k)      = wr
                k         = k - 1
            enddo
        enddo
        nkept = 1 + nsel
        ! the weights over the kept set; extras below wmin dropped (the tail of the descending order), the committed
        ! entry kept whatever its share
        tot = sum(real(w(1:nkept),dp))
        if( .not. (tot > 0.d0) )then
            nkept = 1
            w(1)  = 1.
            return
        endif
        w(1:nkept) = real(real(w(1:nkept),dp) / tot)
        m = 1
        do k = 2,nkept
            if( w(k) < wmin ) exit
            m = k
        enddo
        nkept = m
        tot   = sum(real(w(1:nkept),dp))
        w(1:nkept) = real(real(w(1:nkept),dp) / tot)
    end subroutine top_posterior_extras

    !> the temperature a row is weighted at: the particle's own, or, without one from the fill, a huge one so that
    !! only the row's minimum carries weight (the argmin: the particle has no scale)
    pure real function temperature_or_argmin( temp ) result( t )
        real, intent(in) :: temp
        t = merge(temp, huge(temp), temp > 0.)
    end function temperature_or_argmin

    !> The posterior over one table row at temperature temp: p_r proportional to exp(-temp (dist_r - dist_min))
    !! over the valid distances (finite and below the invalid marker, as ref_assign counts them), zero
    !! elsewhere. On return w holds the posterior sorted ascending with the references in inds (both of
    !! size >= n; the mass sits at the top), imin the mode, pmax its mass, neff the effective number of
    !! references exp(entropy), k99 the number of references holding 99% of the mass, nvalid the number of
    !! valid distances. No valid distance: every output zero. temp = 0 gives the uniform posterior.
    subroutine posterior_row_stats( n, dists, temp, w, inds, imin, pmax, neff, k99, nvalid )
        integer, intent(in)    :: n
        real,    intent(in)    :: dists(n), temp
        real,    intent(inout) :: w(:)
        integer, intent(inout) :: inds(:)
        integer, intent(out)   :: imin, k99, nvalid
        real,    intent(out)   :: pmax, neff
        real(dp) :: z, h, csum
        real     :: dmin, invalid_dist, arg
        integer  :: i
        imin   = 0
        nvalid = 0
        pmax   = 0.
        neff   = 0.
        k99    = 0
        w(1:n) = 0.
        do i = 1,n
            inds(i) = i
        end do
        invalid_dist = 0.1 * huge(invalid_dist)
        dmin = huge(dmin)
        do i = 1,n
            if( ieee_is_finite(dists(i)) .and. dists(i) < invalid_dist )then
                nvalid = nvalid + 1
                if( dists(i) < dmin )then
                    dmin = dists(i)
                    imin = i
                endif
            endif
        end do
        if( nvalid == 0 ) return
        z = 0.d0
        do i = 1,n
            if( ieee_is_finite(dists(i)) .and. dists(i) < invalid_dist )then
                arg = temp * (dists(i) - dmin)
                if( arg < 80. ) w(i) = exp(-arg)
                z = z + real(w(i),dp)
            endif
        end do
        w(1:n) = real(w(1:n) / z)
        h = 0.d0
        do i = 1,n
            if( w(i) > 0. ) h = h - real(w(i),dp) * log(real(w(i),dp))
        end do
        neff = real(exp(h))
        pmax = w(imin)
        call hpsort(w(1:n), inds(1:n))
        csum = 0.d0
        do i = n,1,-1
            csum = csum + real(w(i),dp)
            k99  = k99 + 1
            if( csum >= 0.99d0 ) exit
        end do
    end subroutine posterior_row_stats

    elemental function eulprob_dist_switch( corr, cc_objfun ) result(dist)
        real,    intent(in) :: corr
        integer, intent(in) :: cc_objfun
        real :: dist
        dist = corr
        select case(cc_objfun)
            case(OBJFUN_CC)
                if( corr < 0. )then
                    dist = 0.
                else
                    dist = corr
                endif
                dist = 1. - dist
            case(OBJFUN_EUCLID)
                if( corr < TINY )then
                    dist = huge(dist)
                else
                    dist = - log(corr)
                endif
        end select
    end function eulprob_dist_switch

    elemental function eulprob_corr_switch( dist, cc_objfun ) result(corr)
        real,    intent(in) :: dist
        integer, intent(in) :: cc_objfun
        real :: corr
        corr = dist
        select case(cc_objfun)
            case(OBJFUN_CC)
                corr = 1 - dist
            case(OBJFUN_EUCLID)
                corr = exp(-dist)
        end select
    end function eulprob_corr_switch

    function angle_sampling_1( pvec, athres_ub_in, prob_athres ) result( which )
        real,    intent(in)  :: pvec(:)        !< probabilities
        real,    intent(in)  :: athres_ub_in
        real,    intent(in)  :: prob_athres
        real,    allocatable :: pvec_sorted(:)
        integer, allocatable :: sorted_inds(:)
        integer :: which, n
        n = size(pvec)
        allocate(pvec_sorted(n),sorted_inds(n))
        which = angle_sampling_2(pvec, pvec_sorted, sorted_inds, athres_ub_in, prob_athres)
    end function angle_sampling_1

    function angle_sampling_2( pvec, pvec_sorted, sorted_inds, athres_ub_in, prob_athres ) result( which )
        real,    intent(in)    :: pvec(:)        !< probabilities/distances, lower is better
        real,    intent(inout) :: pvec_sorted(:) !< work buffer, size >= size(pvec)
        integer, intent(inout) :: sorted_inds(:)
        real,    intent(in)    :: athres_ub_in
        real,    intent(in)    :: prob_athres
        integer :: which, n, num_lb, num_ub, i, pick
        real    :: athres_ub, athres_lb, sum_pvec, rnd
        n         = size(pvec)
        athres_ub = min(prob_athres, athres_ub_in)
        athres_lb = min(athres_ub / 10., 1.)
        num_ub    = min(n,max(1,int(athres_ub * real(n) / 180.)))
        num_lb    = 1 + floor(athres_lb / athres_ub * num_ub)
        pvec_sorted(1:num_ub) = pvec(1:num_ub)
        sorted_inds(1:num_ub) = (/(i,i=1,num_ub)/)
        do i = num_ub / 2, 1, -1
            call maxheap_sift_down(pvec_sorted, sorted_inds, i, num_ub)
        enddo
        do i = num_ub + 1, n
            if( pvec(i) < pvec_sorted(1) )then
                pvec_sorted(1) = pvec(i)
                sorted_inds(1) = i
                call maxheap_sift_down(pvec_sorted, sorted_inds, 1, num_ub)
            endif
        enddo
        call hpsort(pvec_sorted(1:num_ub), sorted_inds(1:num_ub))
        rnd      = ran3()
        sum_pvec = sum(pvec_sorted(1:num_ub))
        if( sum_pvec < TINY )then
            pick  = min(num_ub, 1 + floor(real(num_ub) * rnd))
            which = sorted_inds(pick)
        else
            pvec_sorted(1:num_ub) = pvec_sorted(1:num_ub) / sum_pvec
            if( rnd > sum(pvec_sorted(1:num_lb)) )then
                which = sorted_inds(1)
            else
                which = sorted_inds(num_ub)
            endif
        endif

    end function angle_sampling_2

    !> the bounded draw over distances already evaluated, dists(1:n): among the num_ub best (athres_ub of 180 degrees)
    !! the best or, with the probability of the num_lb best's share, the worst of them
    subroutine sample_bounded_dists( n, dists, athres_ub_in, prob_athres, dist, which, pvec_sorted, sorted_inds )
        integer, intent(in)    :: n
        real,    intent(in)    :: dists(:)
        real,    intent(in)    :: athres_ub_in, prob_athres
        real,    intent(out)   :: dist
        integer, intent(out)   :: which
        real,    intent(inout) :: pvec_sorted(:)
        integer, intent(inout) :: sorted_inds(:)
        integer :: num_lb, num_ub, pick
        real    :: athres_ub, athres_lb, sum_pvec, rnd
        athres_ub = min(prob_athres, athres_ub_in)
        athres_lb = min(athres_ub / 10., 1.)
        num_ub    = min(n,max(1,int(athres_ub * real(n) / 180.)))
        num_lb    = 1 + floor(athres_lb / athres_ub * num_ub)
        call select_lowest_dists(n, dists, num_ub, pvec_sorted, sorted_inds)
        rnd      = ran3()
        sum_pvec = sum(pvec_sorted(1:num_ub))
        if( sum_pvec < TINY )then
            pick  = min(num_ub, 1 + floor(real(num_ub) * rnd))
            which = sorted_inds(pick)
        else
            pvec_sorted(1:num_ub) = pvec_sorted(1:num_ub) / sum_pvec
            if( rnd > sum(pvec_sorted(1:num_lb)) )then
                which = sorted_inds(1)
            else
                which = sorted_inds(num_ub)
            endif
        endif
        dist = dists(which)
    end subroutine sample_bounded_dists

    !> the power draw over distances already evaluated, dists(1:n), with the caller's work buffer cdf (size >= nsample):
    !! a rank among the nsample best, the uniform draw raised to power biasing it toward the best
    subroutine sample_power_dists( n, dists, power, nsample, dist, corr, which, pvec_sorted, sorted_inds, cdf )
        integer, intent(in)    :: n, nsample
        real,    intent(in)    :: dists(:)
        real,    intent(in)    :: power
        real,    intent(out)   :: dist, corr
        integer, intent(out)   :: which
        real,    intent(inout) :: pvec_sorted(:)
        integer, intent(inout) :: sorted_inds(:)
        real,    intent(inout) :: cdf(:)
        real    :: rnd, csum
        integer :: num_ub, i, rank, pick
        num_ub = min(n,max(1,nsample))
        call select_lowest_dists(n, dists, num_ub, pvec_sorted, sorted_inds)
        do i = 1,num_ub
            cdf(i) = exp(-pvec_sorted(num_ub-i+1))
        enddo
        if( all(cdf(1:num_ub) < TINY) )then
            pick  = 1
            which = sorted_inds(pick)
            dist  = pvec_sorted(pick)
            corr  = cdf(num_ub)
            return
        endif
        where( cdf(1:num_ub) < TINY ) cdf(1:num_ub) = 0.
        do i = 2,num_ub
            cdf(i) = cdf(i) + cdf(i-1)
        enddo
        csum = cdf(num_ub)
        cdf(1:num_ub) = cdf(1:num_ub) / csum
        rnd = ran3()
        rnd = min(1.,max(0.,1.-rnd**power))
        rank = 0
        do i = 1,num_ub
            if( cdf(i) > rnd )then
                rank = i
                exit
            endif
        enddo
        if( rank == 0 ) rank = num_ub
        pick  = num_ub - rank + 1
        which = sorted_inds(pick)
        dist  = pvec_sorted(pick)
        corr  = exp(-dist)
    end subroutine sample_power_dists

    !> the num_ub lowest of dists(1:n), sorted ascending into pvec_sorted with their indices: a bounded max-heap,
    !! then a sort
    subroutine select_lowest_dists( n, dists, num_ub, pvec_sorted, sorted_inds )
        integer, intent(in)    :: n, num_ub
        real,    intent(in)    :: dists(:)
        real,    intent(inout) :: pvec_sorted(:)
        integer, intent(inout) :: sorted_inds(:)
        integer :: i
        do i = 1,num_ub
            pvec_sorted(i) = dists(i)
            sorted_inds(i) = i
        enddo
        do i = num_ub / 2, 1, -1
            call maxheap_sift_down(pvec_sorted, sorted_inds, i, num_ub)
        enddo
        do i = num_ub + 1,n
            if( dists(i) < pvec_sorted(1) )then
                pvec_sorted(1) = dists(i)
                sorted_inds(1) = i
                call maxheap_sift_down(pvec_sorted, sorted_inds, 1, num_ub)
            endif
        enddo
        call hpsort(pvec_sorted(1:num_ub), sorted_inds(1:num_ub))
    end subroutine select_lowest_dists

    subroutine maxheap_sift_down( pvec_sorted, sorted_inds, root_in, heap_size )
        real,    intent(inout) :: pvec_sorted(:)
        integer, intent(inout) :: sorted_inds(:)
        integer, intent(in)    :: root_in, heap_size
        integer :: root, child, swap_i, tmp_ind
        real    :: tmp_val
        root = root_in
        do
            child = 2 * root
            if( child > heap_size ) exit
            swap_i = root
            if( pvec_sorted(swap_i) < pvec_sorted(child) ) swap_i = child
            if( child + 1 <= heap_size )then
                if( pvec_sorted(swap_i) < pvec_sorted(child + 1) ) swap_i = child + 1
            endif
            if( swap_i == root ) exit
            tmp_val             = pvec_sorted(root)
            pvec_sorted(root)   = pvec_sorted(swap_i)
            pvec_sorted(swap_i) = tmp_val
            tmp_ind             = sorted_inds(root)
            sorted_inds(root)   = sorted_inds(swap_i)
            sorted_inds(swap_i) = tmp_ind
            root = swap_i
        enddo
    end subroutine maxheap_sift_down

    !> the size of a particle's own screen pool (shift_sigma > 0): the smallest leading set of the ascending origin
    !! distances holding SCREEN_QUANTILE of the posterior at temperature temp (weights exp(-temp (d - d_min))); the
    !! breadth the particle's own ambiguity at the origin asks for
    pure integer function screen_pool_size( dists_sorted, n, temp )
        real,     intent(in) :: dists_sorted(:)
        integer,  intent(in) :: n
        real(dp), intent(in) :: temp
        real(dp) :: w(n), tot, acc
        integer  :: k
        screen_pool_size = n
        if( n < 1 ) return
        do k = 1, n
            w(k) = exp(-max(0.d0, temp) * real(dists_sorted(k) - dists_sorted(1), dp))
        enddo
        tot = sum(w)
        if( .not. (tot > 0.d0) ) return
        acc = 0.d0
        do k = 1, n
            acc = acc + w(k)
            if( acc >= real(SCREEN_QUANTILE, dp) * tot )then
                screen_pool_size = k
                return
            endif
        enddo
    end function screen_pool_size

    !> The dynamic centring's rule (shift_sigma > 0): a particle searches the disk when it has no history in this
    !! workflow (the sentinel) or when three of its prior widths exceed the stage's window
    pure logical function shift_needs_disk( shsig, trs )
        real, intent(in) :: shsig, trs
        shift_needs_disk = shsig < 0. .or. 3. * shsig > trs
    end function shift_needs_disk

    !> the disk's radius: three sigma of the map's prior in the calculator's pixels, bounded by the box
    pure real function disk_radius_px( shift_sigma, smpd, box )
        real,    intent(in) :: shift_sigma, smpd
        integer, intent(in) :: box
        disk_radius_px = min(3. * shift_sigma / smpd, 0.5 * real(box) - 1.)
    end function disk_radius_px

    !> the map's prior's weight 1 / (2 sigma^2) on the squared shift from the stored one, sigma in the calculator's
    !! pixels
    pure real(dp) function map_prior_weight( shift_sigma, smpd )
        real, intent(in) :: shift_sigma, smpd
        map_prior_weight = 1.d0 / (2.d0 * (real(shift_sigma, dp) / real(smpd, dp))**2)
    end function map_prior_weight

    !> A pose's loss for the weighted insertion under the polar shift map, the dense and the neighbourhood tables
    !! alike: the data term dist with the prior (weight prior_w) on the pose's increment du from the stored shift
    !! (the prior moves with the particle) at the particle's physical temperature temp, so exp(-temp loss) is the
    !! pose's posterior up to a constant
    pure real function insert_pose_loss( dist, du, prior_w, temp )
        real,     intent(in) :: dist, temp
        real(dp), intent(in) :: du(2), prior_w
        insert_pose_loss = dist
        if( temp > 0. ) insert_pose_loss = dist + real(prior_w * sum(du**2) / real(temp, dp))
    end function insert_pose_loss

    !> the radius of the box a particle's search covers: the disk when it needs one, else the window
    pure real function shift_box_radius( shsig, shift_sigma, smpd, box, trs )
        real,    intent(in) :: shsig, shift_sigma, smpd, trs
        integer, intent(in) :: box
        shift_box_radius = trs
        if( shift_needs_disk(shsig, trs) ) shift_box_radius = disk_radius_px(shift_sigma, smpd, box)
    end function shift_box_radius

    !> The particle's prior width for its next search: the posterior width sdev of this one (the previous width
    !! when the assigned pose carries none), at least the radius r searched when the move d reached that radius's
    !! bound (within a tenth of it, or half a map step where the map is finer): a posterior truncated at the bound
    !! says nothing of the mass beyond it
    pure real function next_shift_width( sdev, sig_prev, d, r, step ) result( sig )
        real, intent(in) :: sdev, sig_prev, d, r, step
        sig = sig_prev
        if( sdev >= 0. ) sig = sdev
        if( d >= max(0.9 * r, r - 0.5 * step) ) sig = max(sig, r)
    end function next_shift_width

    !> The dynamic centring's state after a search, written by the consumer: shincr the move d and shsig the prior
    !! width of the next search (next_shift_width); the radius searched follows from the state before the search
    subroutine record_shift_state( os, iptcl, d, sdev, shift_sigma, smpd, box, trs, kfromto )
        class(oris), intent(inout) :: os
        integer,     intent(in)    :: iptcl, box, kfromto(2)
        real,        intent(in)    :: d, sdev, shift_sigma, smpd, trs
        real :: sig_prev, r, step
        sig_prev = os%get(iptcl, 'shsig')
        r        = shift_box_radius(sig_prev, shift_sigma, smpd, box, trs)
        ! the coarse map's radial step (setup_shiftmaps: the box over four times the coarse band's last shell)
        step = real(box) / (4. * real(max(kfromto(1) + 1, kfromto(2) / 2)))
        call os%set(iptcl, 'shincr', d)
        call os%set(iptcl, 'shsig',  next_shift_width(sdev, sig_prev, d, r, step))
    end subroutine record_shift_state

    !> The scope of the probabilistic alignment's keyed draws of one iteration: the run's seed (SIMPLE_SEED), the
    !! program (refine3D or refine2D), the discrete alignment pass and the iteration. The tables take it at every
    !! fill and assignment and set it on the calculator, so that a draw depends on the key alone.
    function prob_rnd_scope( which_iter, l_2D ) result( scope )
        integer, intent(in) :: which_iter
        logical, intent(in) :: l_2D
        type(rnd_key) :: scope
        scope = rnd_scope_key(rnd_run_seed(), merge(RND_PROG_REFINE2D, RND_PROG_REFINE3D, l_2D), RND_PASS_ALIGN, which_iter)
    end function prob_rnd_scope

    !> The sparse fills' candidate order of a particle: its previous reference (or class) last, the others in a
    !! keyed partial Fisher-Yates permutation (family CAND, one keyed uniform per swap), so the order depends on
    !! the particle and the scope alone, never on the thread
    subroutine keyed_candidate_order( scope, iptcl, n, iprev, order )
        type(rnd_key), intent(in)    :: scope
        integer,       intent(in)    :: iptcl, n, iprev
        integer,       intent(inout) :: order(n)
        integer :: j
        do j = 1, n
            order(j) = j
        enddo
        if( iprev >= 1 .and. iprev <= n )then
            order(iprev) = n
            order(n)     = iprev
            if( n > 1 ) call rnd_keyed_partial_shuffle(rnd_draw_key(scope, iptcl, RND_FAM_CAND), order(1:n-1), n-1)
        else
            call rnd_keyed_partial_shuffle(rnd_draw_key(scope, iptcl, RND_FAM_CAND), order, n)
        endif
    end subroutine keyed_candidate_order

    !> The weighted insertion's cone from the current FSCs: per state the rim limit is the FSC=0.5 resolution,
    !! bounded by the matching low-pass limit when one is set (lp > 0); a state without an FSC file gets no cone.
    !! The rim limits are reported unless l_report is false (a distributed part, whose log the next schedule wipes;
    !! its master reports the same cone from the same files)
    function insert_cone_from_fsc( nstates, box, smpd, mskdiam, lp, l_report ) result( cone )
        integer,           intent(in) :: nstates, box
        real,              intent(in) :: smpd, mskdiam, lp
        logical, optional, intent(in) :: l_report
        type(prob_insert_cone) :: cone
        real, allocatable :: res05(:), fsc(:), res(:)
        type(string) :: fname
        real    :: fsc05, fsc0143
        integer :: s
        logical :: ll_report
        character(len=32) :: limit, angle, shift
        ll_report = .true.
        if( present(l_report) ) ll_report = l_report
        allocate(res05(nstates), source=0.)
        res = get_resarr(box, smpd)
        do s = 1,nstates
            fname = refine3D_fsc_fname(s)
            if( .not. file_exists(fname) ) cycle
            fsc = file2rarr(fname)
            if( size(fsc) < 1 .or. size(fsc) > size(res) ) cycle
            call get_resolution(fsc, res(:size(fsc)), fsc05, fsc0143)
            res05(s) = fsc05
            if( lp > 0. ) res05(s) = min(res05(s), lp)
            deallocate(fsc)
        enddo
        call cone%new(res05, mskdiam, smpd)
        do s = 1,nstates
            if( .not. ll_report ) exit
            if( cone%get_angle(s) < 0. ) cycle
            write(limit,'(F32.2)') res05(s)
            write(angle,'(F32.2)') cone%get_angle(s)
            write(shift,'(F32.2)') res05(s) / smpd
            write(logfhandle,'(A,I0,A)') '>>> PROB_TAB_INSERT CONE (state ', s, &
                &'): rim limit '//trim(adjustl(limit))//' A, i.e. '//trim(adjustl(angle))//' deg of rotation or '//&
                &trim(adjustl(shift))//' px of shift'
        enddo
        call fname%kill
    end function insert_cone_from_fsc

    !> the dynamic centring's disk set of the iteration (shift_sigma > 0): the sampled particles not yet centred in
    !! this workflow and those whose prior width exceeds a third of the window, and the widths' median and 90th
    !! percentile over the centred ones; printed by the project owner after the sample draw
    subroutine print_disk_set( os, l_centring, trs, pinds, nptcls )
        class(oris), intent(inout) :: os
        logical,     intent(in)    :: l_centring
        real,        intent(in)    :: trs
        integer,     intent(in)    :: pinds(:), nptcls
        real, allocatable :: shsig(:), widths(:)
        integer :: i, nnew, nwide, nw
        if( .not. l_centring ) return
        if( nptcls < 1 ) return
        shsig = os%get_all('shsig')
        nnew  = 0
        nwide = 0
        do i = 1, nptcls
            if( shsig(pinds(i)) < 0. )then
                nnew = nnew + 1
            else if( shift_needs_disk(shsig(pinds(i)), trs) )then
                nwide = nwide + 1
            endif
        enddo
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A)') '>>> DISK SEARCHED BY ', nnew + nwide, ' OF ', nptcls,&
            &' SAMPLED PARTICLES (', nnew, ' NOT YET CENTRED, ', nwide, ' WIDE)'
        nw = count(shsig(pinds(1:nptcls)) >= 0.)
        if( nw > 0 )then
            widths = pack(shsig(pinds(1:nptcls)), shsig(pinds(1:nptcls)) >= 0.)
            call hpsort(widths)
            write(logfhandle,'(A,F0.2,A,F0.2,A)') '>>> SHIFT PRIOR WIDTH (PX) MEDIAN ', widths(max(1, nw / 2)),&
                &' 90TH PCT ', widths(max(1, (9 * nw) / 10)), ' OVER THE CENTRED PARTICLES'
            deallocate(widths)
        endif
        deallocate(shsig)
    end subroutine print_disk_set

end module simple_eul_prob_tab_utils
