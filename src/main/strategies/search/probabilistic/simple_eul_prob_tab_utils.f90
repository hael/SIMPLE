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
implicit none

public :: angle_sampling, build_pind_lookup, calc_athres, calc_num2sample
public :: eulprob_corr_switch, eulprob_dist_switch
public :: record_seed_shift, read_seed_shift_table, write_seed_shift_table
public :: sample_bounded_dists, sample_power_dists
public :: prob_candidate, prob_candidate_buffer, prob_candidate_store, prob_fill_profile
public :: prob_dense_entry, prob_refined_entry
public :: prob_posterior_stat, posterior_temperature, posterior_row_stats, draw_posterior_row
public :: prob_beta_at, shift_prior_sigma2, SHIFT_PRIOR_SIGMA_MIN, SHIFT_PRIOR_LP_FRAC
public :: prob_neff_control, prob_beta_control_last, PROB_NEFF_ROWS, PROB_NEFF_MIN_ROWS
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

!> The n_eff control of the tables' posterior draws (prob_neff > 0; doc/implementation_notes/planned/
!! pose_marginalization.md section 20.2): one tempering factor per iteration, solved so that the mean effective
!! number of candidates exp(entropy) over a subsample of rows equals the target. The tables offer their rows in
!! particle order before the draws: every stride-th particle whose row holds at least nmin valid distances is
!! kept, up to nmax rows, packed as distances above the row minimum with the particle's temperature. Rows with
!! fewer candidates than nmin are drawn at the solved factor but do not set it: a row of two candidates lands
!! on one of its best two at any factor, and would pull the level of the rows where the choice is open. Below
!! PROB_NEFF_MIN_ROWS kept rows the control falls back to the flat factor. The last factor applied is kept for
!! the commander (prob_beta_control_last), which hands it to the next iteration's fills as prob_beta.
integer, parameter :: PROB_NEFF_ROWS     = 1000  !< rows of the subsample
integer, parameter :: PROB_NEFF_MIN_ROWS = 50    !< below which the control falls back to prob_beta
real,    parameter :: PROB_NEFF_BETA_LO  = 1.e-3 !< the bracket of the solve, in beta
real,    parameter :: PROB_NEFF_BETA_HI  = 1.e3
integer, parameter :: PROB_NEFF_BISECT   = 40    !< bisections in log beta (a factor 1e6 bracket to 1e-6)
real               :: prob_beta_last_solved = 0. !< the factor the last control applied (0: none this process)

type :: prob_neff_control
    private
    integer :: stride = 1, nmax = 0, nmin = 2, nrows = 0, nvals = 0
    integer, allocatable :: offs(:)
    real,    allocatable :: dists(:), temps(:)
contains
    procedure :: new       => neff_control_new
    procedure :: wanted    => neff_control_wanted
    procedure :: offer     => neff_control_offer
    procedure :: nkept     => neff_control_nkept
    procedure :: mean_neff => neff_control_mean_neff
    procedure :: solve     => neff_control_solve
    procedure :: apply     => neff_control_apply
    procedure :: kill      => neff_control_kill
end type prob_neff_control

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
    !! the cap is slack (pose_marginalization.md section 19, item 10)
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
    !! over 25 iterations is x1.1 per iteration; pose_marginalization.md section 17). beta0 = beta1 is a flat
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

    ! ---- the n_eff control ------------------------------------------------------------------------------------

    !> nptcls rows will be offered; keep up to nmax of them, each with at least nmin valid distances
    subroutine neff_control_new( self, nptcls, nmax, nmin )
        class(prob_neff_control), intent(inout) :: self
        integer,                  intent(in)    :: nptcls, nmax, nmin
        call self%kill
        self%nmax   = max(1, nmax)
        self%nmin   = max(2, nmin)
        self%stride = max(1, nptcls / self%nmax)
        self%nrows  = 0
        self%nvals  = 0
        allocate(self%offs(self%nmax + 1), self%temps(self%nmax), self%dists(64 * self%nmax))
        self%offs(1) = 1
    end subroutine neff_control_new

    !> is particle i's row wanted (its turn in the stride, room left)? the caller skips building the rest
    pure logical function neff_control_wanted( self, i ) result( l_wanted )
        class(prob_neff_control), intent(in) :: self
        integer,                  intent(in) :: i
        l_wanted = self%nrows < self%nmax .and. mod(i - 1, self%stride) == 0
    end function neff_control_wanted

    !> offer a row: its valid distances (finite, below the invalid marker) are packed above their minimum when
    !! there are at least nmin of them and the temperature is positive
    subroutine neff_control_offer( self, n, dists, temp )
        class(prob_neff_control), intent(inout) :: self
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
    end subroutine neff_control_offer

    pure integer function neff_control_nkept( self ) result( nrows )
        class(prob_neff_control), intent(in) :: self
        nrows = self%nrows
    end function neff_control_nkept

    !> the mean over the kept rows of exp(entropy) of the row posterior at tempering factor beta:
    !! p_r proportional to exp(-beta T d_r), d_r the distance above the row minimum; in one pass per row,
    !! H = log z + beta T <d>_p with z the partition sum
    real function neff_control_mean_neff( self, beta ) result( neff )
        class(prob_neff_control), intent(in) :: self
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
    end function neff_control_mean_neff

    !> the factor at which the mean n_eff of the kept rows is the target: a bisection in log beta over
    !! [PROB_NEFF_BETA_LO, PROB_NEFF_BETA_HI] (the mean n_eff falls monotonically with beta: the entropy of
    !! each tempered row does, dH/dbeta = -beta Var(T d)); an unreachable target returns the bracket's end.
    !! neff_at is the mean at the solution, neff_phys the mean at beta 1 (the physical posterior)
    subroutine neff_control_solve( self, target, beta, neff_at, neff_phys )
        class(prob_neff_control), intent(in)  :: self
        real,                     intent(in)  :: target
        real,                     intent(out) :: beta, neff_at, neff_phys
        real    :: lo, hi, mid, neff_lo, neff_hi
        integer :: it
        neff_phys = self%mean_neff(1.)
        lo        = log(PROB_NEFF_BETA_LO)
        hi        = log(PROB_NEFF_BETA_HI)
        neff_lo   = self%mean_neff(PROB_NEFF_BETA_LO)
        neff_hi   = self%mean_neff(PROB_NEFF_BETA_HI)
        if( neff_lo <= target )then
            beta    = PROB_NEFF_BETA_LO
            neff_at = neff_lo
            return
        endif
        if( neff_hi >= target )then
            beta    = PROB_NEFF_BETA_HI
            neff_at = neff_hi
            return
        endif
        do it = 1,PROB_NEFF_BISECT
            mid     = 0.5 * (lo + hi)
            neff_at = self%mean_neff(exp(mid))
            if( neff_at > target )then
                lo = mid
            else
                hi = mid
            endif
        enddo
        beta    = exp(0.5 * (lo + hi))
        neff_at = self%mean_neff(beta)
    end subroutine neff_control_solve

    !> the control of one assignment: the factor the draws use. The solved factor when enough rows were kept,
    !! the flat factor beta_flat otherwise; in the annealing tail (anneal_frac in (0,1]) the factor ramps
    !! geometrically from the solved value to 1, the physical posterior, so the tail of section 19.8 keeps its
    !! shape from wherever the control sits. Logs one line; keeps the factor for the commander
    subroutine neff_control_apply( self, target, beta_flat, anneal_frac, beta )
        class(prob_neff_control), intent(in)  :: self
        real,                     intent(in)  :: target, beta_flat, anneal_frac
        real,                     intent(out) :: beta
        real :: neff_at, neff_phys, frac
        if( self%nrows < PROB_NEFF_MIN_ROWS )then
            beta = beta_flat
            write(logfhandle,'(A,I0,A,I0,A,F9.4)') '>>> POSTERIOR DRAW CONTROL: ', self%nrows, &
                &' rows with at least ', self%nmin, ' candidates, fewer than needed; the flat factor ', beta
        else
            call self%solve(target, beta, neff_at, neff_phys)
            frac = min(1., max(0., anneal_frac))
            if( frac > 0. .and. beta > 0. )then
                beta = beta**(1. - frac)
                write(logfhandle,'(A,F9.4,A,F6.2,A,F8.2,A,I0,A,F8.2,A,F4.2,A)') '>>> POSTERIOR DRAW CONTROL: BETA ', beta, &
                    &' (solved for n_eff ', target, ': ', neff_at, ' over ', self%nrows, ' rows; at beta 1: ', neff_phys, &
                    &'; annealing tail, frac ', frac, ' to 1)'
            else
                write(logfhandle,'(A,F9.4,A,F6.2,A,F8.2,A,I0,A,F8.2,A)') '>>> POSTERIOR DRAW CONTROL: BETA ', beta, &
                    &' FOR N_EFF ', target, ' (', neff_at, ' over ', self%nrows, ' rows; n_eff at beta 1: ', neff_phys, ')'
            endif
        endif
        prob_beta_last_solved = beta
        call flush(logfhandle)
    end subroutine neff_control_apply

    subroutine neff_control_kill( self )
        class(prob_neff_control), intent(inout) :: self
        if( allocated(self%offs)  ) deallocate(self%offs)
        if( allocated(self%dists) ) deallocate(self%dists)
        if( allocated(self%temps) ) deallocate(self%temps)
        self%nrows = 0
        self%nvals = 0
    end subroutine neff_control_kill

    !> the factor the last n_eff control of this process applied; 0 when none ran
    real function prob_beta_control_last() result( beta )
        beta = prob_beta_last_solved
    end function prob_beta_control_last

    !> One draw from the posterior over a table row at temperature temp (the tables' assignment,
    !! doc/implementation_notes/planned/pose_marginalization.md section 4.1): p_r proportional to
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
