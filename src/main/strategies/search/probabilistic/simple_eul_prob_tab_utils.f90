!@descr: shared utility routines for probabilistic alignment tables
module simple_eul_prob_tab_utils
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,      only: TINY, dp, DPI
use simple_error,     only: simple_exception
use simple_math,      only: hpsort
use simple_oris,      only: oris
use simple_rnd,       only: ran3
use simple_type_defs, only: OBJFUN_CC, OBJFUN_EUCLID, ptcl_ref
implicit none

public :: angle_sampling, build_pind_lookup, calc_athres, calc_num2sample
public :: eulprob_corr_switch, eulprob_dist_switch
public :: record_seed_shift, read_seed_shift_table, write_seed_shift_table
public :: sample_bounded_dists, sample_likelihood_dist, sample_likelihood_dists, sample_power_dists, sample_likelihood_index
public :: prob_candidate, prob_candidate_buffer, prob_candidate_store, prob_fill_profile
public :: prob_dense_entry, prob_refined_entry
public :: prob_posterior_stat, posterior_temperature, posterior_row_stats
private
#include "simple_local_flags.inc"

type :: prob_candidate
    integer :: iref   = 0
    integer :: inpl   = 0
    real    :: dist   = huge(1.0)
    real    :: x      = 0.
    real    :: y      = 0.
    logical :: has_sh = .false.
end type prob_candidate

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
    ! the production row (likelihood-sampled in-plane, the best projs_ns shift-refined): mode, loss,
    ! temperature, mass, breadth and spread
    integer :: iptcl = 0, istate = 0, iproj = 0, k99 = 0, nvalid = 0
    real    :: lmin = 0., temp = 0., pmax = 0., neff = 0., spread = 0.
    ! the identical-pose null test, in nats: the kernel re-scored at the production mode's exact (reference,
    ! shift, grid rotation) minus the stored loss, with the stored shift rotated back into the kernel's frame
    ! (null_a, expected zero) and taken as it is (null_b, the alternative convention)
    real    :: null_a = 0., null_b = 0.
    ! the deterministic row on the own slot (best grid rotation at the seed shift, no sampling, no
    ! refinement): the search test against the production row; rank_d the production mode's rank in it,
    ! dist_d the angle between the two modes
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

abstract interface
    function dist_eval_fun( ind ) result( dist )
        integer, intent(in) :: ind
        real :: dist
    end function dist_eval_fun
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

    ! Truncated, numerically-stabilized softmax over a lazily-evaluated distance function.
    ! The distances are noise-variance-normalized negative-log-likelihoods (whitened per
    ! Fourier shell by the estimated sigma2_noise upstream), so exp(-dist) is a properly
    ! calibrated likelihood weight and no external temperature/tau factor is required. The
    ! per-candidate min-shift only guards against overflow; it cancels in the normalized
    ! weights. Only the best nsample candidates form the sampling support (top-K sampling).
    subroutine sample_likelihood_dist( n, dist_fun, nsample, dist, corr, which, pvec_sorted, sorted_inds )
        integer,   intent(in)    :: n, nsample
        procedure(dist_eval_fun)  :: dist_fun
        real,      intent(out)   :: dist, corr
        integer,   intent(out)   :: which
        real,      intent(inout) :: pvec_sorted(:)
        integer,   intent(inout) :: sorted_inds(:)
        real    :: cdf(max(1,nsample)), dist_tmp
        integer :: num_ub, i
        num_ub = min(n,max(1,nsample))
        do i = 1,num_ub
            pvec_sorted(i) = dist_fun(i)
            sorted_inds(i) = i
        enddo
        do i = num_ub / 2, 1, -1
            call maxheap_sift_down(pvec_sorted, sorted_inds, i, num_ub)
        enddo
        do i = num_ub + 1,n
            dist_tmp = dist_fun(i)
            if( dist_tmp < pvec_sorted(1) )then
                pvec_sorted(1) = dist_tmp
                sorted_inds(1) = i
                call maxheap_sift_down(pvec_sorted, sorted_inds, 1, num_ub)
            endif
        enddo
        call hpsort(pvec_sorted(1:num_ub), sorted_inds(1:num_ub))
        call draw_likelihood(num_ub, pvec_sorted, sorted_inds, cdf, dist, corr, which)
    end subroutine sample_likelihood_dist

    !> sample_likelihood_dist over distances already evaluated, dists(1:n), with the caller's work buffer cdf
    !! (size >= nsample): the same draw without an internal procedure or an allocation per call
    subroutine sample_likelihood_dists( n, dists, nsample, dist, corr, which, pvec_sorted, sorted_inds, cdf )
        integer, intent(in)    :: n, nsample
        real,    intent(in)    :: dists(:)
        real,    intent(out)   :: dist, corr
        integer, intent(out)   :: which
        real,    intent(inout) :: pvec_sorted(:)
        integer, intent(inout) :: sorted_inds(:)
        real,    intent(inout) :: cdf(:)
        integer :: num_ub
        num_ub = min(n,max(1,nsample))
        call select_lowest_dists(n, dists, num_ub, pvec_sorted, sorted_inds)
        call draw_likelihood(num_ub, pvec_sorted, sorted_inds, cdf, dist, corr, which)
    end subroutine sample_likelihood_dists

    !> the draw of the likelihood samplers from the num_ub best distances, sorted ascending
    subroutine draw_likelihood( num_ub, pvec_sorted, sorted_inds, cdf, dist, corr, which )
        integer, intent(in)    :: num_ub
        real,    intent(in)    :: pvec_sorted(:)
        integer, intent(in)    :: sorted_inds(:)
        real,    intent(inout) :: cdf(:)
        real,    intent(out)   :: dist, corr
        integer, intent(out)   :: which
        real    :: rnd, csum, min_dist
        integer :: i, rank, pick
        min_dist = pvec_sorted(1)
        do i = 1,num_ub
            cdf(i) = exp(-(pvec_sorted(num_ub-i+1) - min_dist))
        enddo
        if( all(cdf(1:num_ub) < TINY) )then
            pick  = 1
            which = sorted_inds(pick)
            dist  = pvec_sorted(pick)
            corr  = exp(-dist)
            return
        endif
        where( cdf(1:num_ub) < TINY ) cdf(1:num_ub) = 0.
        do i = 2,num_ub
            cdf(i) = cdf(i) + cdf(i-1)
        enddo
        csum = cdf(num_ub)
        cdf(1:num_ub) = cdf(1:num_ub) / csum
        rnd = ran3()
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
    end subroutine draw_likelihood

    ! Array-based counterpart of sample_likelihood_dist for the global assignment stage.
    ! Given the already-evaluated distances of n candidates, draw one index with the same
    ! truncated, min-shifted softmax exp(-(dist-min_dist)) over the best ntop candidates.
    ! Candidates with dist >= huge_val (exhausted/sentinel entries) are excluded. sorted and
    ! order are caller-provided work buffers of size >= n. Shared by the 2D sparse and dense
    ! class-assignment frontiers to keep a single softmax implementation.
    subroutine sample_likelihood_index( n, dists, ntop, huge_val, sorted, order, which )
        integer, intent(in)    :: n, ntop
        real,    intent(in)    :: dists(:)
        real,    intent(in)    :: huge_val
        real,    intent(inout) :: sorted(:)
        integer, intent(inout) :: order(:)
        integer, intent(out)   :: which
        integer :: j, ntop_eff
        real    :: min_dist, csum, rnd, acc, weight
        if( n < 1 )then
            which = 1
            return
        endif
        do j = 1,n
            sorted(j) = dists(j)
            order(j)  = j
        enddo
        call hpsort(sorted(1:n), order(1:n))
        min_dist = sorted(1)
        ntop_eff = min(max(1,ntop), n)
        csum = 0.
        do j = 1,ntop_eff
            if( sorted(j) >= huge_val ) cycle
            weight = exp(-(sorted(j) - min_dist))
            if( weight > TINY ) csum = csum + weight
        enddo
        if( csum < TINY )then
            which = order(1)
            return
        endif
        rnd = ran3() * csum
        acc = 0.
        do j = 1,ntop_eff
            if( sorted(j) >= huge_val ) cycle
            weight = exp(-(sorted(j) - min_dist))
            if( weight <= TINY ) cycle
            acc = acc + weight
            if( acc >= rnd )then
                which = order(j)
                return
            endif
        enddo
        which = order(1)
    end subroutine sample_likelihood_index

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
