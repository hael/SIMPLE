!@descr: shared utility routines for probabilistic alignment tables
module simple_eul_prob_tab_utils
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,      only: TINY
use simple_error,     only: simple_exception
use simple_math,      only: hpsort, rotmat2d
use simple_oris,      only: oris
use simple_rnd,       only: ran3
use simple_type_defs, only: OBJFUN_CC, OBJFUN_EUCLID, ptcl_ref
implicit none

public :: angle_sampling, build_pind_lookup, calc_athres, calc_num2sample
public :: eulprob_corr_switch, eulprob_dist_switch
public :: materialize_seed_shift, read_seed_shift_table, write_seed_shift_table
public :: sample_bounded_dist, sample_likelihood_dist, sample_power_dist, sample_likelihood_index
public :: prob_candidate, prob_candidate_buffer, prob_candidate_store
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

    subroutine materialize_seed_shift( assgn, seed_shift, seed_has_sh, l_doshift, seed_nrots )
        type(ptcl_ref),       intent(inout) :: assgn
        real,                 intent(in)    :: seed_shift(2)
        logical,              intent(in)    :: seed_has_sh, l_doshift
        integer,              intent(in)    :: seed_nrots
        real    :: rotmat(2,2), rot_xy(2)
        integer :: irot
        if( .not. l_doshift ) return
        if( assgn%has_sh    ) return
        if( .not. seed_has_sh ) return
        irot = assgn%inpl
        if( irot < 1 .or. irot > seed_nrots .or. seed_nrots < 1 )then
            assgn%x      = 0.
            assgn%y      = 0.
            assgn%has_sh = .true.
            return
        endif
        call rotmat2d(real(irot - 1) * 360. / real(seed_nrots), rotmat)
        rot_xy(1) = seed_shift(1) * rotmat(1,1) + seed_shift(2) * rotmat(2,1)
        rot_xy(2) = seed_shift(1) * rotmat(1,2) + seed_shift(2) * rotmat(2,2)
        assgn%x      = rot_xy(1)
        assgn%y      = rot_xy(2)
        assgn%has_sh = .true.
    end subroutine materialize_seed_shift

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

    subroutine sample_bounded_dist( n, dist_fun, athres_ub_in, prob_athres, dist, which, pvec_sorted, sorted_inds )
        integer,   intent(in)    :: n
        procedure(dist_eval_fun)  :: dist_fun
        real,      intent(in)    :: athres_ub_in, prob_athres
        real,      intent(out)   :: dist
        integer,   intent(out)   :: which
        real,      intent(inout) :: pvec_sorted(:)
        integer,   intent(inout) :: sorted_inds(:)
        integer :: num_lb, num_ub, i, pick
        real    :: athres_ub, athres_lb, sum_pvec, rnd, dist_tmp
        athres_ub = min(prob_athres, athres_ub_in)
        athres_lb = min(athres_ub / 10., 1.)
        num_ub    = min(n,max(1,int(athres_ub * real(n) / 180.)))
        num_lb    = 1 + floor(athres_lb / athres_ub * num_ub)
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
        dist = dist_fun(which)
    end subroutine sample_bounded_dist

    subroutine sample_power_dist( n, dist_fun, power, nsample, dist, corr, which, pvec_sorted, sorted_inds )
        integer,   intent(in)    :: n, nsample
        procedure(dist_eval_fun)  :: dist_fun
        real,      intent(in)    :: power
        real,      intent(out)   :: dist, corr
        integer,   intent(out)   :: which
        real,      intent(inout) :: pvec_sorted(:)
        integer,   intent(inout) :: sorted_inds(:)
        real    :: cdf(max(1,nsample)), rnd, dist_tmp, csum
        integer :: num_ub, i, rank, pick
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
    end subroutine sample_power_dist

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
        real    :: cdf(max(1,nsample)), rnd, dist_tmp, csum, min_dist
        integer :: num_ub, i, rank, pick
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
    end subroutine sample_likelihood_dist

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
