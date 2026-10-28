!@descr: binary classification statistics for class-average quality analysis
module simple_cavg_quality_stats
use simple_error,              only: simple_exception
use simple_defs,               only: dp
use simple_stat,               only: median, mad_gau
use simple_srch_sort_loc,      only: hpsort
use simple_cavg_quality_types, only: EPS
implicit none
private
#include "simple_local_flags.inc"

public :: calc_confusion
public :: calc_binary_metrics
public :: auc_for_values
public :: median_by_state
public :: mad_by_state
public :: safe_div

contains

    subroutine calc_confusion( pred, ref, tp, fp, tn, fn )
        logical, intent(in)  :: pred(:), ref(:)
        integer, intent(out) :: tp, fp, tn, fn
        if( size(pred) /= size(ref) ) THROW_HARD('calc_confusion: size mismatch')
        tp = count( pred .and.  ref)
        fp = count( pred .and. .not. ref)
        tn = count(.not. pred .and. .not. ref)
        fn = count(.not. pred .and.  ref)
    end subroutine calc_confusion

    subroutine calc_binary_metrics( tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy )
        integer, intent(in)  :: tp, fp, tn, fn
        real,    intent(out) :: precision, recall, specificity, f1, balacc, accuracy
        precision   = safe_div(real(tp), real(tp + fp))
        recall      = safe_div(real(tp), real(tp + fn))
        specificity = safe_div(real(tn), real(tn + fp))
        f1          = safe_div(2.0 * precision * recall, precision + recall)
        balacc      = 0.5 * (recall + specificity)
        accuracy    = safe_div(real(tp + tn), real(tp + fp + tn + fn))
    end subroutine calc_binary_metrics

    ! 0.9-1.0: Excellent
    ! 0.8-0.9: Good
    ! 0.7-0.8: Fair
    ! 0.6-0.7: Poor
    ! 0.5:     No better than chance
    !> Mann-Whitney AUC: a good value above a bad one wins 1, one at most EPS below it 1/2. Counted against the
    !! sorted bad values (two bisections per good value) in double precision with a 64-bit-safe denominator, so it
    !! holds at trajectory scale (tens of millions of pairs) and equals the pairwise count where that is exact.
    real function auc_for_values( vals, reference_states )
        real,    intent(in) :: vals(:)
        integer, intent(in) :: reference_states(:)
        real,    allocatable :: bad(:)
        real(dp) :: wins
        integer  :: i, ngood, nbad, nbelow, nnear
        if( size(vals) /= size(reference_states) ) THROW_HARD('auc_for_values: size mismatch')
        ngood = count(reference_states > 0)
        nbad  = size(reference_states) - ngood
        if( ngood == 0 .or. nbad == 0 ) then
            auc_for_values = 0.5
            return
        end if
        bad = pack(vals, reference_states <= 0)
        call hpsort(bad)
        wins = 0._dp
        do i = 1, size(vals)
            if( reference_states(i) <= 0 ) cycle
            nbelow = count_below(vals(i))
            nnear  = count_within_eps(vals(i)) - nbelow
            wins   = wins + real(nbelow, dp) + 0.5_dp * real(nnear, dp)
        end do
        auc_for_values = real(wins / (real(ngood, dp) * real(nbad, dp)))
        deallocate(bad)

        contains

            ! number of bad values strictly below v
            integer function count_below( v )
                real, intent(in) :: v
                integer :: lo, hi, mid
                lo = 0
                hi = nbad
                do while( lo < hi )
                    mid = lo + (hi - lo) / 2 + 1
                    if( bad(mid) < v )then
                        lo = mid
                    else
                        hi = mid - 1
                    endif
                end do
                count_below = lo
            end function count_below

            ! number of bad values below v or at most EPS above it
            integer function count_within_eps( v )
                real, intent(in) :: v
                integer :: lo, hi, mid
                lo = 0
                hi = nbad
                do while( lo < hi )
                    mid = lo + (hi - lo) / 2 + 1
                    if( bad(mid) - v <= EPS )then
                        lo = mid
                    else
                        hi = mid - 1
                    endif
                end do
                count_within_eps = lo
            end function count_within_eps

    end function auc_for_values

    real function median_by_state( vals, mask )
        real,    intent(in) :: vals(:)
        logical, intent(in) :: mask(:)
        real, allocatable :: packed(:)
        if( size(vals) /= size(mask) ) THROW_HARD('median_by_state: size mismatch')
        if( count(mask) == 0 ) then
            median_by_state = 0.0
            return
        end if
        packed = pack(vals, mask)
        median_by_state = median(packed)
        deallocate(packed)
    end function median_by_state

    real function mad_by_state( vals, mask, med )
        real,    intent(in) :: vals(:), med
        logical, intent(in) :: mask(:)
        real, allocatable :: packed(:)
        if( size(vals) /= size(mask) ) THROW_HARD('mad_by_state: size mismatch')
        if( count(mask) < 2 ) then
            mad_by_state = 0.0
            return
        end if
        packed = pack(vals, mask)
        mad_by_state = mad_gau(packed, med)
        deallocate(packed)
    end function mad_by_state

    real function safe_div( num, den )
        real, intent(in) :: num, den
        if( abs(den) <= EPS ) then
            safe_div = 0.0
        else
            safe_div = num / den
        end if
    end function safe_div

end module simple_cavg_quality_stats
