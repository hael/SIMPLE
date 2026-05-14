!@descr: binary classification statistics for class-average quality analysis
module simple_cavg_quality_stats
use simple_error,              only: simple_exception
use simple_stat,               only: median, mad_gau
use simple_cavg_quality_feats, only: EPS
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

    real function auc_for_values( vals, reference_states )
        real,    intent(in) :: vals(:)
        integer, intent(in) :: reference_states(:)
        integer :: i, j, ngood, nbad
        real    :: wins
        if( size(vals) /= size(reference_states) ) THROW_HARD('auc_for_values: size mismatch')
        ngood = count(reference_states > 0)
        nbad  = size(reference_states) - ngood
        if( ngood == 0 .or. nbad == 0 ) then
            auc_for_values = 0.5
            return
        end if
        wins = 0.0
        do i = 1, size(vals)
            if( reference_states(i) <= 0 ) cycle
            do j = 1, size(vals)
                if( reference_states(j) > 0 ) cycle
                if( vals(i) > vals(j) ) then
                    wins = wins + 1.0
                else if( abs(vals(i) - vals(j)) <= EPS ) then
                    wins = wins + 0.5
                end if
            end do
        end do
        auc_for_values = wins / real(ngood * nbad)
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
