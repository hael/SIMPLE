!@descr: shared utility routines for probabilistic alignment tables
module simple_eul_prob_tab_utils
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,      only: TINY, dp, DPI, logfhandle
use simple_error,     only: simple_exception
use simple_math,      only: hpsort
use simple_oris,      only: oris
use simple_rnd,       only: ran3
use simple_type_defs, only: OBJFUN_CC, OBJFUN_EUCLID, ptcl_ref
use simple_string,    only: string
use simple_fileio,    only: fopen, fclose, fileiochk, del_file, del_files
use simple_syslib,    only: file_exists
use simple_defs_fname, only: DIST_FBODY, ASSIGNMENT_FBODY, INSERT_FBODY
implicit none

public :: angle_sampling, build_pind_lookup, calc_athres, calc_num2sample
public :: eulprob_corr_switch, eulprob_dist_switch
public :: record_seed_shift, read_seed_shift_table, write_seed_shift_table
public :: sample_bounded_dists, sample_power_dists
public :: prob_candidate, prob_candidate_buffer, prob_candidate_store, prob_fill_profile
public :: prob_dense_entry, prob_refined_entry
public :: prob_posterior_stat, posterior_temperature, posterior_row_stats, draw_posterior_row, draw_posterior_extras
public :: prob_beta_at, shift_prior_sigma2, SHIFT_PRIOR_SIGMA_MIN, SHIFT_PRIOR_LP_FRAC
public :: prob_neff_readout, PROB_NEFF_ROWS, PROB_NEFF_NMIN, PROB_NEFF_MIN_ROWS
public :: prob_insert_pose, prob_insert_table, prob_insert_poses, del_prob_align_files
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

type :: prob_candidate
    integer :: iref   = 0
    integer :: inpl   = 0
    real    :: dist   = huge(1.0)
    real    :: x      = 0.
    real    :: y      = 0.
    logical :: has_sh = .false.
end type prob_candidate

!> The n_eff readout of the tables' posterior draws: the mean effective number of candidates exp(entropy) over a subsample of rows, at the draw's
!! tempering factor and at beta 1 (the physical posterior), logged every iteration as POSTERIOR DRAW N_EFF. The
!! tables offer their rows in particle order before the draws: every stride-th particle whose row holds at least
!! PROB_NEFF_NMIN valid distances is kept, up to PROB_NEFF_ROWS rows, packed as distances above the row minimum
!! with the particle's temperature; a row of two candidates lands on one of its best two at any factor and would
!! pull the level of the rows where the choice is open. Until 2026-10-29 the same object solved the factor so that
!! the mean n_eff met a target (prob_neff, the control); pinning n_eff cancelled the data's own annealing at the
!! band entries and failed PorF 7 of 10, so the schedule is prob_beta flat with the tail, and n_eff is a readout
integer, parameter :: PROB_NEFF_ROWS     = 1000  !< rows of the subsample
integer, parameter :: PROB_NEFF_NMIN     = 6     !< valid candidates a row needs to be kept
integer, parameter :: PROB_NEFF_MIN_ROWS = 50    !< below which nothing is reported

type :: prob_neff_readout
    private
    integer :: stride = 1, nmax = 0, nmin = 2, nrows = 0, nvals = 0
    integer, allocatable :: offs(:)
    real,    allocatable :: dists(:), temps(:)
contains
    procedure :: new       => neff_readout_new
    procedure :: wanted    => neff_readout_wanted
    procedure :: offer     => neff_readout_offer
    procedure :: nkept     => neff_readout_nkept
    procedure :: mean_neff => neff_readout_mean_neff
    procedure :: report    => neff_readout_report
    procedure :: kill      => neff_readout_kill
end type prob_neff_readout

!> One extra pose of a particle for the weighted insertion (topk > 1): the reference (the full index over the states), its in-plane index, its
!! refined shift in the stored frame and its weight. The committed pose is the project's, with the weight the
!! table keeps per particle
type :: prob_insert_pose
    integer :: iref = 0, inpl = 0
    real    :: x = 0., y = 0., w = 0.
end type prob_insert_pose

!> The weighted insertion's table of one iteration: per particle (the global table's order) the committed pose's
!! weight and its extra poses, drawn by the assignment (draw_posterior_extras) and written by the master beside
!! the assignment (INSERT_FBODY); the workers read it with the particle indices and insert each particle at its
!! committed pose with its weight and at its extras with theirs, every particle's weights summing to 1
type :: prob_insert_table
    private
    integer :: nptcls = 0, kmax = 1, npose = 0
    integer,                allocatable :: pinds(:), nextra(:), first(:)
    real,                   allocatable :: wc(:)
    type(prob_insert_pose), allocatable :: poses(:)
contains
    procedure :: new    => insert_table_new
    procedure :: add    => insert_table_add
    procedure :: subset => insert_table_subset
    procedure :: write  => insert_table_write
    procedure :: read   => insert_table_read
    procedure :: report => insert_table_report
    procedure :: get_nptcls   => insert_table_get_nptcls
    procedure :: get_kmax     => insert_table_get_kmax
    procedure :: get_pind     => insert_table_get_pind
    procedure :: get_committed_weight => insert_table_get_committed_weight
    procedure :: get_nextra   => insert_table_get_nextra
    procedure :: get_extra    => insert_table_get_extra
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
    procedure :: setup => insert_poses_setup
    procedure :: has   => insert_poses_has
    procedure :: kill  => insert_poses_kill
end type prob_insert_poses

! One entry of the dense table, a particle and a reference: the distance and the in-plane index. No default values:
! the fill and the read write every entry.
type :: prob_dense_entry
    real    :: dist
    integer :: inpl
end type prob_dense_entry

! A refined entry of the dense table: the particle (part-local in a part file, the master's index on the master),
! the compact reference and the shift of its search
type :: prob_refined_entry
    integer :: i
    integer :: ri
    real    :: x
    real    :: y
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

    subroutine build_pind_lookup( glob_pinds, loc_pinds, pind2glob, max_pind )
        integer,              intent(in)  :: glob_pinds(:), loc_pinds(:)
        integer, allocatable, intent(out) :: pind2glob(:)
        integer,              intent(out) :: max_pind
        integer :: i, pind
        max_pind = max(maxval(glob_pinds), maxval(loc_pinds))
        if( max_pind < 1 )then
            allocate(pind2glob(0))
            return
        endif
        allocate(pind2glob(max_pind), source=0)
        do i = 1, size(glob_pinds)
            pind = glob_pinds(i)
            if( pind > 0 .and. pind <= max_pind ) pind2glob(pind) = i
        enddo
    end subroutine build_pind_lookup

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

    subroutine write_seed_shift_table( funit, addr, seed_nrots, seed_shifts, seed_has_sh )
        integer, intent(in)    :: funit
        integer(int64), intent(inout) :: addr
        integer, intent(in)    :: seed_nrots
        real,    intent(in)    :: seed_shifts(:,:)
        logical, intent(in)    :: seed_has_sh(:)
        write(funit, pos=addr) seed_nrots
        addr = addr + sizeof(seed_nrots)
        write(funit, pos=addr) seed_shifts
        addr = addr + sizeof(seed_shifts)
        write(funit, pos=addr) seed_has_sh
        addr = addr + sizeof(seed_has_sh)
    end subroutine write_seed_shift_table

    subroutine read_seed_shift_table( funit, addr, seed_nrots, seed_shifts, seed_has_sh )
        integer, intent(in)    :: funit
        integer(int64), intent(inout) :: addr
        integer, intent(out)   :: seed_nrots
        real,    intent(out)   :: seed_shifts(:,:)
        logical, intent(out)   :: seed_has_sh(:)
        read(funit, pos=addr) seed_nrots
        addr = addr + sizeof(seed_nrots)
        read(funit, pos=addr) seed_shifts
        addr = addr + sizeof(seed_shifts)
        read(funit, pos=addr) seed_has_sh
        addr = addr + sizeof(seed_has_sh)
    end subroutine read_seed_shift_table

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
    subroutine neff_readout_new( self, nptcls, nmax, nmin )
        class(prob_neff_readout), intent(inout) :: self
        integer,                  intent(in)    :: nptcls, nmax, nmin
        call self%kill
        self%nmax   = max(1, nmax)
        self%nmin   = max(2, nmin)
        self%stride = max(1, nptcls / self%nmax)
        self%nrows  = 0
        self%nvals  = 0
        allocate(self%offs(self%nmax + 1), self%temps(self%nmax), self%dists(64 * self%nmax))
        self%offs(1) = 1
    end subroutine neff_readout_new

    !> is particle i's row wanted (its turn in the stride, room left)? the caller skips building the rest
    pure logical function neff_readout_wanted( self, i ) result( l_wanted )
        class(prob_neff_readout), intent(in) :: self
        integer,                  intent(in) :: i
        l_wanted = self%nrows < self%nmax .and. mod(i - 1, self%stride) == 0
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
        if( self%nrows >= self%nmax ) return
        if( .not. (temp > 0.) )       return
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

    !> the mean over the kept rows of exp(entropy) of the row posterior at tempering factor beta:
    !! p_r proportional to exp(-beta T d_r), d_r the distance above the row minimum; in one pass per row,
    !! H = log z + beta T <d>_p with z the partition sum
    real function neff_readout_mean_neff( self, beta ) result( neff )
        class(prob_neff_readout), intent(in) :: self
        real,                     intent(in) :: beta
        real(dp) :: z, s, acc, a, e
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
        neff = 0.
        if( self%nrows > 0 ) neff = real(acc / real(self%nrows,dp))
    end function neff_readout_mean_neff

    !> the readout of one assignment: the mean n_eff at the draw's factor beta_draw and at beta 1 (the physical
    !! posterior), over the kept rows; nothing below PROB_NEFF_MIN_ROWS rows
    subroutine neff_readout_report( self, beta_draw )
        class(prob_neff_readout), intent(in) :: self
        real,                     intent(in) :: beta_draw
        if( self%nrows < PROB_NEFF_MIN_ROWS ) return
        write(logfhandle,'(A,F8.2,A,F7.4,A,F8.2,A,I0,A)') '>>> POSTERIOR DRAW N_EFF: ', self%mean_neff(beta_draw), &
            &' at beta ', beta_draw, ';', self%mean_neff(1.), ' at beta 1 (over ', self%nrows, ' rows)'
        call flush(logfhandle)
    end subroutine neff_readout_report

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
    !! (reference, in-plane index and refined shift from the candidate); a particle's extras are added once
    subroutine insert_table_add( self, i, wc, nextra, cands, w )
        class(prob_insert_table), intent(inout) :: self
        integer,                  intent(in)    :: i, nextra
        real,                     intent(in)    :: wc
        type(prob_candidate),     intent(in)    :: cands(:)
        real,                     intent(in)    :: w(:)
        type(prob_insert_pose), allocatable :: tmp(:)
        integer :: k
        if( i < 1 .or. i > self%nptcls ) THROW_HARD('particle index out of range; insert_table_add')
        if( self%nextra(i) > 0 ) THROW_HARD('a particle''s extra poses are added once; insert_table_add')
        self%wc(i) = wc
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

    !> this table as the part of glob that holds the particles loc_pinds, in their order: a particle glob has no
    !! entry for keeps weight 1 and no extras
    subroutine insert_table_subset( self, glob, loc_pinds )
        class(prob_insert_table), intent(inout) :: self
        class(prob_insert_table), intent(in)    :: glob
        integer,                  intent(in)    :: loc_pinds(:)
        integer, allocatable :: pind2glob(:)
        integer :: max_pind, i_loc, i_glob, pind, k, n
        call self%new(size(loc_pinds), loc_pinds, glob%kmax)
        if( glob%nptcls < 1 ) return
        call build_pind_lookup(glob%pinds(1:glob%nptcls), loc_pinds, pind2glob, max_pind)
        if( max_pind < 1 ) return
        do i_loc = 1,size(loc_pinds)
            pind = loc_pinds(i_loc)
            if( pind < 1 .or. pind > max_pind ) cycle
            i_glob = pind2glob(pind)
            if( i_glob < 1 ) cycle
            n = glob%nextra(i_glob)
            self%wc(i_loc) = glob%wc(i_glob)
            if( n < 1 ) cycle
            if( self%npose + n > size(self%poses) )then
                block
                    type(prob_insert_pose), allocatable :: tmp(:)
                    allocate(tmp(2 * size(self%poses) + n))
                    tmp(1:self%npose) = self%poses(1:self%npose)
                    call move_alloc(tmp, self%poses)
                end block
            endif
            self%first(i_loc)  = self%npose + 1
            self%nextra(i_loc) = n
            do k = 1,n
                self%npose = self%npose + 1
                self%poses(self%npose) = glob%poses(glob%first(i_glob) + k - 1)
            enddo
        enddo
        deallocate(pind2glob)
    end subroutine insert_table_subset

    !> the table to a binary file: nptcls, kmax, npose; pinds, wc, nextra, first; the poses
    subroutine insert_table_write( self, fname )
        class(prob_insert_table), intent(in) :: self
        class(string),            intent(in) :: fname
        integer :: funit, io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('insert_table_write; file: '//fname%to_char(), io_stat)
        write(unit=funit) self%nptcls, self%kmax, self%npose
        if( self%nptcls > 0 )then
            write(unit=funit) self%pinds(1:self%nptcls), self%wc(1:self%nptcls), self%nextra(1:self%nptcls), &
                &self%first(1:self%nptcls)
        endif
        if( self%npose > 0 ) write(unit=funit) self%poses(1:self%npose)
        call fclose(funit)
    end subroutine insert_table_write

    subroutine insert_table_read( self, fname )
        class(prob_insert_table), intent(inout) :: self
        class(string),            intent(in)    :: fname
        integer :: funit, io_stat, nptcls, kmax, npose
        call self%kill
        if( .not. file_exists(fname) ) THROW_HARD('file '//fname%to_char()//' does not exist; insert_table_read')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('insert_table_read; file: '//fname%to_char(), io_stat)
        read(unit=funit) nptcls, kmax, npose
        self%nptcls = nptcls
        self%kmax   = kmax
        self%npose  = npose
        allocate(self%pinds(max(1,nptcls)), self%nextra(max(1,nptcls)), self%first(max(1,nptcls)), source=0)
        allocate(self%wc(max(1,nptcls)), source=1.)
        allocate(self%poses(max(1,npose)))
        if( nptcls > 0 ) read(unit=funit) self%pinds(1:nptcls), self%wc(1:nptcls), self%nextra(1:nptcls), self%first(1:nptcls)
        if( npose  > 0 ) read(unit=funit) self%poses(1:npose)
        call fclose(funit)
    end subroutine insert_table_read

    !> the log line: insertions per particle and the committed pose's mean weight
    subroutine insert_table_report( self )
        class(prob_insert_table), intent(in) :: self
        if( self%nptcls < 1 ) return
        write(logfhandle,'(A,I0,A,F6.3,A,F6.3)') '>>> PROB_TAB_INSERT: topk ', self%kmax, '; insertions per particle ', &
            &1. + real(self%npose) / real(self%nptcls), '; mean weight of the committed pose ', &
            &sum(self%wc(1:self%nptcls)) / real(self%nptcls)
        call flush(logfhandle)
    end subroutine insert_table_report

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
        self%nptcls = 0
        self%kmax   = 1
        self%npose  = 0
    end subroutine insert_table_kill

    ! ---- the resolved poses of a part -----------------------------------------------------------------------------

    !> the layout from the part's table tab (its particles in the order of pinds, tab's order): counts, committed
    !! weights and slots; the slots' poses are zero (weight 0, which the insertion skips) until the search fills them
    subroutine insert_poses_setup( self, tab, pinds )
        class(prob_insert_poses), intent(inout) :: self
        class(prob_insert_table), intent(in)    :: tab
        integer,                  intent(in)    :: pinds(:)
        integer :: i, ipose, max_pind
        call self%kill
        if( tab%nptcls /= size(pinds) ) THROW_HARD('the table is not the part''s; insert_poses_setup')
        self%nptcls = tab%nptcls
        self%npose  = tab%npose
        max_pind    = 0
        if( size(pinds) > 0 ) max_pind = maxval(pinds)
        allocate(self%loc(max(1,max_pind)), self%nextra(max(1,self%nptcls)), self%first(max(1,self%nptcls)), source=0)
        allocate(self%wc(max(1,self%nptcls)), source=1.)
        allocate(self%euler(3,max(1,self%npose)), self%shift(2,max(1,self%npose)), source=0.)
        allocate(self%w(max(1,self%npose)), source=0.)
        ipose = 0
        do i = 1,self%nptcls
            if( tab%pinds(i) /= pinds(i) ) THROW_HARD('the table''s particles are not the part''s; insert_poses_setup')
            if( pinds(i) > 0 ) self%loc(pinds(i)) = i
            self%wc(i)     = tab%wc(i)
            self%nextra(i) = tab%nextra(i)
            self%first(i)  = ipose + 1
            ipose = ipose + tab%nextra(i)
        enddo
    end subroutine insert_poses_setup

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

    !> the files the probabilistic alignment of an iteration leaves in the run directory: the tables' parts, the
    !! assignment and the weighted insertion's extra poses; one place for the four workflows that remove them
    subroutine del_prob_align_files( nparts )
        integer, intent(in) :: nparts
        call del_files(DIST_FBODY,       nparts, ext='.dat')
        call del_files(ASSIGNMENT_FBODY, nparts, ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        call del_file(INSERT_FBODY//'.dat')
    end subroutine del_prob_align_files

    !> One draw from the posterior over a table row at temperature temp (the tables' assignment): p_r proportional to
    !! exp(-temp (dist_r - dist_min)) over the valid distances (finite, below the invalid marker), zero
    !! elsewhere. which is the drawn index and dist its distance; which = 0 when no distance is valid.
    !! temp = 0 draws uniformly over the valid distances; a huge temp returns the argmin when the minimum
    !! is unique and one of the equal minima with equal probability otherwise (a flat row is a uniform draw
    !! at any temperature: the posterior, not the balanced rule's lowest-index convention). w is a work
    !! array of size >= n; the one random number is ran3's, so the stream is the caller's.
    subroutine draw_posterior_row( n, dists, temp, w, which, dist )
        integer, intent(in)    :: n
        real,    intent(in)    :: dists(n), temp
        real,    intent(inout) :: w(:)
        integer, intent(out)   :: which
        real,    intent(out)   :: dist
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
        tgt    = real(ran3(),dp) * wsum
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

    !> Extra poses for the weighted insertion (topk > 1): up to kmax - 1 further entries of a row, drawn WITHOUT
    !! replacement from the posterior p_r proportional to exp(-temp (dist_r - dist_min)) over the eligible entries
    !! (mask: the refined ones) other than the committed draw ifirst, each draw from the remainder renormalised,
    !! with the weights that posterior renormalised over the committed entry and the extras. temp is the
    !! particle's own temperature, untempered: the committed pose is the search's draw from the tempered row, hot
    !! so that the search explores, but the map is weighted by what the data say of each pose, so the
    !! reconstruction carries the data's pose uncertainty and not the search's (on PfCRT, extras drawn and
    !! weighted at the search's factor 0.3 left the committed pose 0.5-0.65 of its particle where the untempered
    !! posterior held about 1.4 effective poses: a blur that capped the resolution a stage early, 2026-10-29). An
    !! extra whose weight falls below wmin is dropped (the rest renormalised again), so no pose enters the
    !! reconstruction with less than that share of its particle; the committed entry is always kept, so kmax = 1
    !! keeps it alone at weight 1 and draws no random number: the single draw's stream and reconstruction, bit for
    !! bit. On return nkept entries, the committed first, with their row indices inds and weights w summing to 1.
    !! work is a work array of size >= n
    subroutine draw_posterior_extras( n, dists, temp, mask, ifirst, kmax, wmin, work, nkept, inds, w )
        integer, intent(in)    :: n, ifirst, kmax
        real,    intent(in)    :: dists(n), temp, wmin
        logical, intent(in)    :: mask(n)
        real,    intent(inout) :: work(:)
        integer, intent(out)   :: nkept
        integer, intent(out)   :: inds(:)
        real,    intent(out)   :: w(:)
        real(dp) :: wsum, csum, tgt, arg, tot
        real     :: dmin, invalid
        integer  :: r, k, m
        nkept   = 1
        inds(1) = ifirst
        w(1)    = 1.
        if( kmax <= 1 ) return
        invalid = 0.1 * huge(invalid)
        dmin    = huge(dmin)
        do r = 1,n
            if( ieee_is_finite(dists(r)) .and. dists(r) < invalid ) dmin = min(dmin, dists(r))
        enddo
        if( dmin >= invalid ) return
        ! the unnormalised posterior of the committed entry and of the eligible others
        w(1) = 0.
        if( ieee_is_finite(dists(ifirst)) .and. dists(ifirst) < invalid )then
            arg = real(max(0.,temp),dp) * real(dists(ifirst) - dmin,dp)
            if( arg < 80.d0 ) w(1) = real(exp(-arg))
        endif
        wsum = 0.d0
        do r = 1,n
            work(r) = 0.
            if( r == ifirst .or. .not. mask(r) ) cycle
            if( .not. (ieee_is_finite(dists(r)) .and. dists(r) < invalid) ) cycle
            arg = real(max(0.,temp),dp) * real(dists(r) - dmin,dp)
            if( arg < 80.d0 )then
                work(r) = real(exp(-arg))
                wsum    = wsum + real(work(r),dp)
            endif
        enddo
        ! the extras, without replacement, each from the remainder
        do k = 2,kmax
            if( .not. (wsum > 0.d0) ) exit
            tgt  = real(ran3(),dp) * wsum
            csum = 0.d0
            m    = 0
            do r = 1,n
                if( work(r) <= 0. ) cycle
                csum = csum + real(work(r),dp)
                m    = r
                if( csum >= tgt ) exit
            enddo
            if( m == 0 ) exit
            nkept       = nkept + 1
            inds(nkept) = m
            w(nkept)    = work(m)
            wsum        = wsum - real(work(m),dp)
            work(m)     = 0.
        enddo
        ! the weights over the kept set; extras below wmin dropped, the committed entry kept whatever its share
        tot = sum(real(w(1:nkept),dp))
        if( .not. (tot > 0.d0) )then
            nkept = 1
            w(1)  = 1.
            return
        endif
        w(1:nkept) = real(real(w(1:nkept),dp) / tot)
        m = 1
        do k = 2,nkept
            if( w(k) < wmin ) cycle
            m       = m + 1
            inds(m) = inds(k)
            w(m)    = w(k)
        enddo
        nkept = m
        tot   = sum(real(w(1:nkept),dp))
        w(1:nkept) = real(real(w(1:nkept),dp) / tot)
    end subroutine draw_posterior_extras

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

end module simple_eul_prob_tab_utils
