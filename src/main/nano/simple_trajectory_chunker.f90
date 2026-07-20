!@descr: time-constrained flex-latent segmentation for nanoparticle trajectories
module simple_trajectory_chunker
use simple_core_module_api
use simple_projected_latent_result, only: projected_latent_fit_result
use simple_srch_sort_loc,          only: hpsort
implicit none

public :: trajectory_chunk, trajectory_chunk_plan
public :: make_trajectory_chunk_plan, select_trajectory_chunk_plan
public :: trajectory_chunks_to_parts, write_trajectory_chunks_csv
private
#include "simple_local_flags.inc"

type :: trajectory_chunk
    integer  :: fromp = 0
    integer  :: top = 0
    integer  :: frame_start = 0
    integer  :: frame_end = 0
    integer  :: nframes = 0
    integer  :: medoid_pind = 0
    integer  :: medoid_frame = 0
    real(dp) :: cost = 0.d0
end type trajectory_chunk

type :: trajectory_chunk_plan
    type(trajectory_chunk), allocatable :: chunks(:)
    real(dp), allocatable :: mode_weights(:)
    real(dp) :: total_cost = 0.d0
    real(dp) :: temporal_silhouette = 0.d0
    real(dp) :: selection_score = 0.d0
contains
    procedure :: kill => kill_trajectory_chunk_plan
end type trajectory_chunk_plan

contains

    subroutine make_trajectory_chunk_plan( fit, frame_inds, nchunks_in, target_len_in, min_len_in, &
        &max_len_in, max_shift_in, plan )
        type(projected_latent_fit_result), intent(in) :: fit
        integer, intent(in) :: frame_inds(:)
        integer, intent(in) :: nchunks_in, target_len_in, min_len_in, max_len_in, max_shift_in
        type(trajectory_chunk_plan), intent(inout) :: plan
        real(dp), allocatable :: x(:,:), weights(:,:), prefix_w(:,:), prefix_wx(:,:), prefix_wx2(:,:)
        real(dp), allocatable :: dp_cost(:,:), mode_evidence(:)
        integer, allocatable :: sorted_frames(:), perm(:), pinds(:), backptr(:,:), bound_lo(:), bound_hi(:)
        integer :: n, ncomp, nchunks, target_len, min_len, max_len, max_shift
        integer :: c, i, q, last, prev, first, best_prev, medoid
        real(dp) :: val, best_val
        call plan%kill
        n     = fit%nptcls
        ncomp = fit%ncomp
        if( n <= 0 .or. ncomp <= 0 ) THROW_HARD('trajectory chunking requires a nonempty projected latent fit')
        if( size(frame_inds) /= n ) THROW_HARD('trajectory chunking frame index count does not match latent fit')
        if( .not. allocated(fit%pinds) .or. .not. allocated(fit%z) .or. &
            &.not. allocated(fit%z_postcov) .or. .not. allocated(fit%basis_energy) )then
            THROW_HARD('trajectory chunking received an incomplete projected latent fit')
        endif
        nchunks = max(1, min(nchunks_in, n))
        target_len = target_len_in
        if( target_len <= 0 ) target_len = max(1, nint(real(n) / real(nchunks)))
        min_len = min_len_in
        if( min_len <= 0 ) min_len = max(2, target_len / 2)
        max_len = max_len_in
        if( max_len <= 0 ) max_len = max(min_len, 2 * target_len)
        max_shift = max_shift_in
        if( max_shift <= 0 ) max_shift = max(1, target_len / 2)
        if( min_len * nchunks > n ) THROW_HARD('trajectory chunk_min_len is incompatible with the requested chunk count')
        if( max_len * nchunks < n ) THROW_HARD('trajectory chunk_max_len is incompatible with the requested chunk count')

        allocate(sorted_frames(n), source=frame_inds)
        allocate(perm(n), source=[(i, i=1,n)])
        call hpsort(sorted_frames, perm)
        if( n > 1 )then
            if( any(sorted_frames(2:n) <= sorted_frames(1:n-1)) )then
                THROW_HARD('trajectory frame identifiers must be strictly increasing after sorting')
            endif
        endif
        allocate(pinds(n), source=fit%pinds(perm))
        if( n > 1 )then
            if( any(pinds(2:n) /= pinds(1:n-1) + 1) )then
                THROW_HARD('time-sorted trajectory particles are not a contiguous project-row range')
            endif
        endif

        allocate(x(n,ncomp), weights(n,ncomp), mode_evidence(ncomp), source=0.d0)
        call prepare_weighted_features(fit, perm, x, weights, mode_evidence)
        allocate(plan%mode_weights(ncomp), source=mode_evidence)
        allocate(prefix_w(ncomp,0:n), prefix_wx(ncomp,0:n), prefix_wx2(ncomp,0:n), source=0.d0)
        do i = 1, n
            do q = 1, ncomp
                prefix_w(q,i)   = prefix_w(q,i-1)   + weights(i,q)
                prefix_wx(q,i)  = prefix_wx(q,i-1)  + weights(i,q) * x(i,q)
                prefix_wx2(q,i) = prefix_wx2(q,i-1) + weights(i,q) * x(i,q) * x(i,q)
            end do
        end do

        allocate(bound_lo(0:nchunks), bound_hi(0:nchunks), source=0)
        call make_boundary_bands(n, nchunks, min_len, max_len, max_shift, bound_lo, bound_hi)
        allocate(dp_cost(nchunks,0:n), source=huge(1.d0))
        allocate(backptr(nchunks,0:n), source=-1)
        do last = bound_lo(1), bound_hi(1)
            if( last < min_len .or. last > max_len ) cycle
            dp_cost(1,last) = interval_cost(1, last, prefix_w, prefix_wx, prefix_wx2)
            backptr(1,last) = 0
        end do
        do c = 2, nchunks
            do last = bound_lo(c), bound_hi(c)
                best_val  = huge(1.d0)
                best_prev = -1
                do prev = bound_lo(c-1), bound_hi(c-1)
                    if( dp_cost(c-1,prev) >= huge(1.d0) * 0.5d0 ) cycle
                    if( last - prev < min_len .or. last - prev > max_len ) cycle
                    val = dp_cost(c-1,prev) + interval_cost(prev + 1, last, prefix_w, prefix_wx, prefix_wx2)
                    if( val < best_val )then
                        best_val  = val
                        best_prev = prev
                    endif
                end do
                if( best_prev >= 0 )then
                    dp_cost(c,last) = best_val
                    backptr(c,last) = best_prev
                endif
            end do
        end do
        if( dp_cost(nchunks,n) >= huge(1.d0) * 0.5d0 )then
            THROW_HARD('time-constrained latent chunking found no feasible contiguous partition')
        endif

        allocate(plan%chunks(nchunks))
        plan%total_cost = dp_cost(nchunks,n)
        last = n
        do c = nchunks, 1, -1
            prev  = backptr(c,last)
            first = prev + 1
            medoid = interval_medoid(first, last, x, prefix_w, prefix_wx)
            plan%chunks(c)%fromp        = pinds(first)
            plan%chunks(c)%top          = pinds(last)
            plan%chunks(c)%frame_start  = sorted_frames(first)
            plan%chunks(c)%frame_end    = sorted_frames(last)
            plan%chunks(c)%nframes      = last - first + 1
            plan%chunks(c)%medoid_pind  = pinds(medoid)
            plan%chunks(c)%medoid_frame = sorted_frames(medoid)
            plan%chunks(c)%cost = interval_cost(first, last, prefix_w, prefix_wx, prefix_wx2)
            last = prev
        end do
        write(logfhandle,'(A,I0,A,ES12.4)') '>>> TRAJECTORY_CHUNK TIME-CONSTRAINED CHUNKS: ', &
            &nchunks, ' total_cost=', plan%total_cost
        do q = 1, ncomp
            write(logfhandle,'(A,I0,A,ES12.4)') '>>> TRAJECTORY_CHUNK MODE WEIGHT ', q, ': ', plan%mode_weights(q)
        end do
        call flush(logfhandle)
        deallocate(sorted_frames, perm, pinds, x, weights, mode_evidence, prefix_w, prefix_wx, prefix_wx2, &
            &dp_cost, backptr, bound_lo, bound_hi)
    end subroutine make_trajectory_chunk_plan

    subroutine select_trajectory_chunk_plan( fit, frame_inds, nchunks_min_in, nchunks_max_in, &
        &min_len_in, max_len_in, max_shift_in, count_penalty, plan, scan_fname )
        type(projected_latent_fit_result), intent(in) :: fit
        integer, intent(in) :: frame_inds(:)
        integer, intent(in) :: nchunks_min_in, nchunks_max_in, min_len_in, max_len_in, max_shift_in
        real, intent(in) :: count_penalty
        type(trajectory_chunk_plan), intent(inout) :: plan
        character(len=*), intent(in) :: scan_fname
        type(trajectory_chunk_plan) :: candidate
        integer :: nchunks_min, nchunks_max, nchunks, target_len, funit, io_stat
        real(dp) :: silhouette, complexity_penalty, score, best_score
        logical :: found
        call plan%kill
        if( nchunks_min_in < 1 .or. nchunks_max_in < nchunks_min_in )then
            THROW_HARD('automatic trajectory chunking requires 1 <= nchunks_min <= nchunks_max')
        endif
        nchunks_min = max(1, nchunks_min_in)
        nchunks_max = min(fit%nptcls, nchunks_max_in)
        if( nchunks_min > nchunks_max )then
            THROW_HARD('automatic trajectory chunk-count range contains no candidates')
        endif
        open(newunit=funit, file=trim(scan_fname), status='replace', action='write', iostat=io_stat)
        call fileiochk('opening '//trim(scan_fname), io_stat)
        write(funit,'(A)') 'NCHUNKS,VALID,TOTAL_COST,ADJACENT_CENTROID_SILHOUETTE,COMPLEXITY_PENALTY,SCORE'
        best_score = -huge(1.d0)
        found = .false.
        do nchunks = nchunks_min, nchunks_max
            target_len = max(1, nint(real(fit%nptcls) / real(nchunks)))
            if( .not. chunk_count_is_feasible(fit%nptcls, nchunks, target_len, min_len_in, &
                &max_len_in, max_shift_in) )then
                write(funit,'(I0,A,A,A)') nchunks, ',', 'no', ',,,,'
                cycle
            endif
            call make_trajectory_chunk_plan(fit, frame_inds, nchunks, target_len, min_len_in, &
                &max_len_in, max_shift_in, candidate)
            silhouette = temporal_partition_silhouette(fit, frame_inds, candidate)
            complexity_penalty = max(0.d0, real(count_penalty,dp)) * real(nchunks - 1,dp)
            score = silhouette - complexity_penalty
            candidate%temporal_silhouette = silhouette
            candidate%selection_score = score
            write(funit,'(I0,A,A,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') nchunks, ',', 'yes', ',', &
                &candidate%total_cost, ',', silhouette, ',', complexity_penalty, ',', score
            if( .not. found .or. score > best_score )then
                plan = candidate
                best_score = score
                found = .true.
            endif
            call candidate%kill
        end do
        close(funit)
        if( .not. found )then
            THROW_HARD('automatic trajectory chunking found no feasible chunk count')
        endif
        write(logfhandle,'(A,I0,A,F8.4,A,F8.4)') '>>> TRAJECTORY_CHUNK SELECTED NCHUNKS: ', &
            &size(plan%chunks), ' silhouette=', plan%temporal_silhouette, ' score=', plan%selection_score
        write(logfhandle,'(A,A)') '>>> TRAJECTORY_CHUNK WROTE COUNT SCAN: ', trim(scan_fname)
        call flush(logfhandle)
    end subroutine select_trajectory_chunk_plan

    subroutine prepare_weighted_features( fit, perm, x, weights, evidence )
        type(projected_latent_fit_result), intent(in) :: fit
        integer, intent(in) :: perm(:)
        real(dp), intent(out) :: x(:,:), weights(:,:), evidence(:)
        real(dp), allocatable :: raw_evidence(:)
        real(dp) :: avg, var, post_avg, signal_var, denom, sumw
        integer :: n, ncomp, i, q, row
        n = size(perm)
        ncomp = fit%ncomp
        allocate(raw_evidence(ncomp), source=0.d0)
        do q = 1, ncomp
            avg = sum(fit%z(:,q)) / real(n,dp)
            var = sum((fit%z(:,q) - avg) ** 2) / real(max(1,n-1),dp)
            post_avg = 0.d0
            do i = 1, n
                post_avg = post_avg + max(0.d0, fit%z_postcov(i,q,q))
            end do
            post_avg = post_avg / real(n,dp)
            signal_var = max(0.d0, var - post_avg)
            raw_evidence(q) = max(0.d0, fit%basis_energy(q)) * signal_var
            denom = sqrt(max(var, DTINY))
            do i = 1, n
                row = perm(i)
                x(i,q) = (fit%z(row,q) - avg) / denom
                weights(i,q) = max(0.d0, fit%z_postcov(row,q,q)) / max(var, DTINY)
            end do
        end do
        if( sum(raw_evidence) <= DTINY )then
            raw_evidence = max(0.d0, fit%basis_energy)
            if( sum(raw_evidence) <= DTINY ) raw_evidence(1) = 1.d0
        endif
        evidence = raw_evidence / sum(raw_evidence)
        do q = 1, ncomp
            sumw = 0.d0
            do i = 1, n
                weights(i,q) = 1.d0 / (1.d0 + weights(i,q))
                sumw = sumw + weights(i,q)
            end do
            if( sumw > DTINY )then
                weights(:,q) = weights(:,q) * (evidence(q) * real(n,dp) / sumw)
            else
                weights(:,q) = evidence(q)
            endif
        end do
        deallocate(raw_evidence)
    end subroutine prepare_weighted_features

    subroutine make_boundary_bands( n, nchunks, min_len, max_len, max_shift, bound_lo, bound_hi )
        integer, intent(in) :: n, nchunks, min_len, max_len, max_shift
        integer, intent(out) :: bound_lo(0:nchunks), bound_hi(0:nchunks)
        integer :: c, center
        bound_lo = 0
        bound_hi = 0
        do c = 1, nchunks - 1
            center = nint(real(c*n) / real(nchunks))
            bound_lo(c) = max(c * min_len, n - (nchunks-c) * max_len, center - max_shift)
            bound_hi(c) = min(c * max_len, n - (nchunks-c) * min_len, center + max_shift)
            if( bound_lo(c) > bound_hi(c) )then
                THROW_HARD('trajectory chunk boundary band is empty; relax length or shift constraints')
            endif
        end do
        bound_lo(nchunks) = n
        bound_hi(nchunks) = n
    end subroutine make_boundary_bands

    logical function chunk_count_is_feasible( n, nchunks, target_len, min_len_in, max_len_in, &
        &max_shift_in ) result(feasible)
        integer, intent(in) :: n, nchunks, target_len, min_len_in, max_len_in, max_shift_in
        integer :: min_len, max_len, max_shift, c, center, bound_lo, bound_hi
        min_len = min_len_in
        if( min_len <= 0 ) min_len = max(2, target_len / 2)
        max_len = max_len_in
        if( max_len <= 0 ) max_len = max(min_len, 2 * target_len)
        max_shift = max_shift_in
        if( max_shift <= 0 ) max_shift = max(1, target_len / 2)
        feasible = min_len * nchunks <= n .and. max_len * nchunks >= n
        if( .not. feasible ) return
        do c = 1, nchunks - 1
            center = nint(real(c*n) / real(nchunks))
            bound_lo = max(c * min_len, n - (nchunks-c) * max_len, center - max_shift)
            bound_hi = min(c * max_len, n - (nchunks-c) * min_len, center + max_shift)
            if( bound_lo > bound_hi )then
                feasible = .false.
                return
            endif
        end do
    end function chunk_count_is_feasible

    real(dp) function temporal_partition_silhouette( fit, frame_inds, plan ) result(score)
        type(projected_latent_fit_result), intent(in) :: fit
        integer, intent(in) :: frame_inds(:)
        type(trajectory_chunk_plan), intent(in) :: plan
        real(dp), allocatable :: x(:,:), weights(:,:), evidence(:), centers(:,:), center_weights(:,:)
        integer, allocatable :: sorted_frames(:), perm(:), chunk_ids(:)
        real(dp) :: own_dist, adjacent_dist, denom
        integer :: n, ncomp, nchunks, i, q, c, first, last
        n = fit%nptcls
        ncomp = fit%ncomp
        nchunks = size(plan%chunks)
        if( nchunks <= 1 )then
            score = 0.d0
            return
        endif
        allocate(sorted_frames(n), source=frame_inds)
        allocate(perm(n), source=[(i, i=1,n)])
        call hpsort(sorted_frames, perm)
        allocate(x(n,ncomp), weights(n,ncomp), evidence(ncomp), source=0.d0)
        call prepare_weighted_features(fit, perm, x, weights, evidence)
        allocate(centers(nchunks,ncomp), center_weights(nchunks,ncomp), source=0.d0)
        allocate(chunk_ids(n), source=0)
        first = 1
        do c = 1, nchunks
            last = first + plan%chunks(c)%nframes - 1
            if( last > n ) THROW_HARD('trajectory chunk plan exceeds latent frame count')
            chunk_ids(first:last) = c
            do i = first, last
                do q = 1, ncomp
                    centers(c,q) = centers(c,q) + weights(i,q) * x(i,q)
                    center_weights(c,q) = center_weights(c,q) + weights(i,q)
                end do
            end do
            first = last + 1
        end do
        if( first /= n + 1 ) THROW_HARD('trajectory chunk plan does not cover all latent frames')
        do c = 1, nchunks
            do q = 1, ncomp
                if( center_weights(c,q) > DTINY ) centers(c,q) = centers(c,q) / center_weights(c,q)
            end do
        end do
        score = 0.d0
        do i = 1, n
            c = chunk_ids(i)
            own_dist = weighted_center_distance(x(i,:), weights(i,:), centers(c,:))
            adjacent_dist = huge(1.d0)
            if( c > 1 ) adjacent_dist = min(adjacent_dist, &
                &weighted_center_distance(x(i,:), weights(i,:), centers(c-1,:)))
            if( c < nchunks ) adjacent_dist = min(adjacent_dist, &
                &weighted_center_distance(x(i,:), weights(i,:), centers(c+1,:)))
            denom = max(own_dist, adjacent_dist)
            if( denom > DTINY ) score = score + (adjacent_dist - own_dist) / denom
        end do
        score = score / real(n,dp)
        deallocate(sorted_frames, perm, x, weights, evidence, centers, center_weights, chunk_ids)
    end function temporal_partition_silhouette

    pure real(dp) function weighted_center_distance( xrow, wrow, center ) result(dist)
        real(dp), intent(in) :: xrow(:), wrow(:), center(:)
        dist = sum(wrow * (xrow - center) ** 2)
    end function weighted_center_distance

    real(dp) function interval_cost( first, last, prefix_w, prefix_wx, prefix_wx2 ) result( cost )
        integer, intent(in) :: first, last
        real(dp), intent(in) :: prefix_w(:,0:), prefix_wx(:,0:), prefix_wx2(:,0:)
        real(dp) :: sw, swx, swx2
        integer :: q
        cost = 0.d0
        do q = 1, size(prefix_w,dim=1)
            sw   = prefix_w(q,last)   - prefix_w(q,first-1)
            swx  = prefix_wx(q,last)  - prefix_wx(q,first-1)
            swx2 = prefix_wx2(q,last) - prefix_wx2(q,first-1)
            if( sw > DTINY ) cost = cost + max(0.d0, swx2 - swx * swx / sw)
        end do
    end function interval_cost

    integer function interval_medoid( first, last, x, prefix_w, prefix_wx ) result( medoid )
        integer, intent(in) :: first, last
        real(dp), intent(in) :: x(:,:), prefix_w(:,0:), prefix_wx(:,0:)
        real(dp), allocatable :: center(:)
        real(dp) :: sw, dist, best_dist
        integer :: i, q
        allocate(center(size(x,dim=2)), source=0.d0)
        do q = 1, size(x,dim=2)
            sw = prefix_w(q,last) - prefix_w(q,first-1)
            if( sw > DTINY ) center(q) = (prefix_wx(q,last) - prefix_wx(q,first-1)) / sw
        end do
        medoid = first
        best_dist = huge(1.d0)
        do i = first, last
            dist = 0.d0
            do q = 1, size(x,dim=2)
                sw = prefix_w(q,last) - prefix_w(q,first-1)
                dist = dist + sw * (x(i,q) - center(q)) ** 2
            end do
            if( dist < best_dist )then
                best_dist = dist
                medoid = i
            endif
        end do
        deallocate(center)
    end function interval_medoid

    subroutine trajectory_chunks_to_parts( plan, parts )
        type(trajectory_chunk_plan), intent(in) :: plan
        integer, allocatable, intent(inout) :: parts(:,:)
        integer :: i
        if( allocated(parts) ) deallocate(parts)
        if( .not. allocated(plan%chunks) ) THROW_HARD('trajectory chunk plan is empty')
        allocate(parts(size(plan%chunks),2), source=0)
        do i = 1, size(plan%chunks)
            parts(i,:) = [plan%chunks(i)%fromp, plan%chunks(i)%top]
        end do
    end subroutine trajectory_chunks_to_parts

    subroutine write_trajectory_chunks_csv( plan, fname )
        type(trajectory_chunk_plan), intent(in) :: plan
        character(len=*), intent(in) :: fname
        integer :: funit, io_stat, i, lifetime
        if( .not. allocated(plan%chunks) ) return
        open(newunit=funit, file=trim(fname), status='replace', action='write', iostat=io_stat)
        call fileiochk('opening '//trim(fname), io_stat)
        write(funit,'(A)') 'PARTITION,FROMP,TOP,FRAME_START,FRAME_END,NFRAMES,LIFETIME,MEDOID_PIND,MEDOID_FRAME,COST'
        do i = 1, size(plan%chunks)
            lifetime = plan%chunks(i)%frame_end - plan%chunks(i)%frame_start + 1
            write(funit,'(I0,8(A,I0),A,ES14.6)') i, ',', plan%chunks(i)%fromp, ',', plan%chunks(i)%top, ',', &
                &plan%chunks(i)%frame_start, ',', plan%chunks(i)%frame_end, ',', plan%chunks(i)%nframes, ',', &
                &lifetime, ',', plan%chunks(i)%medoid_pind, ',', plan%chunks(i)%medoid_frame, ',', plan%chunks(i)%cost
        end do
        close(funit)
    end subroutine write_trajectory_chunks_csv

    subroutine kill_trajectory_chunk_plan( self )
        class(trajectory_chunk_plan), intent(inout) :: self
        if( allocated(self%chunks) ) deallocate(self%chunks)
        if( allocated(self%mode_weights) ) deallocate(self%mode_weights)
        self%total_cost = 0.d0
        self%temporal_silhouette = 0.d0
        self%selection_score = 0.d0
    end subroutine kill_trajectory_chunk_plan

end module simple_trajectory_chunker
