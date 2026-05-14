!@descr: evaluation and analyze-mode reporting for class-average quality analysis
module simple_cavg_quality_analysis
use simple_defs,               only: logfhandle, LONGSTRLEN
use simple_error,              only: simple_exception
use simple_image,              only: image
use simple_oris,               only: oris
use simple_cavg_quality_feats, only: cavg_quality_feature_name, &
    calc_histogram_quality_signal, calc_power_spectrum_quality_signal, extract_cavg_quality_features, &
    normalize_cavg_quality_features, write_cavg_quality_feature_inventory, I_HIST_ENTROPY, I_LOG_HIST_VARIANCE
use simple_cavg_quality_model, only: cavg_quality_model
use simple_cavg_quality_stats, only: calc_confusion, calc_binary_metrics, auc_for_values, &
    median_by_state, mad_by_state, safe_div
use simple_cavg_quality_types, only: CAVG_QUALITY_HIST_DMAT_METRIC, CAVG_QUALITY_NFEATS, &
    CAVG_QUALITY_SPEC_DMAT_METRIC, EPS, cavg_quality_result, reset_cavg_quality_result
implicit none
private
#include "simple_local_flags.inc"

public :: evaluate_cavg_quality
public :: write_cavg_quality_analysis
public :: write_cavg_quality_feature_table

contains

    subroutine evaluate_cavg_quality( imgs, cls_oris, mskdiam, quality, model )
        class(image),              intent(inout) :: imgs(:)
        type(oris),                intent(in)    :: cls_oris
        real,                      intent(in)    :: mskdiam
        type(cavg_quality_result), intent(inout) :: quality
        type(cavg_quality_model),  intent(in)    :: model
        real, allocatable :: hist_entropy(:), hist_variance(:)
        call reset_cavg_quality_result(quality)
        call extract_cavg_quality_features(imgs, cls_oris, mskdiam, quality%raw, quality%hard_reject)
        call calc_histogram_quality_signal(imgs, mskdiam, quality%hard_reject, quality%hist_dmat, &
                                           hist_entropy, hist_variance)
        if( allocated(hist_entropy) )then
            if( size(hist_entropy) /= size(quality%raw, dim=1) ) &
                THROW_HARD('evaluate_cavg_quality: invalid histogram entropy size')
            quality%raw(:, I_HIST_ENTROPY) = hist_entropy
            deallocate(hist_entropy)
        endif
        if( allocated(hist_variance) )then
            if( size(hist_variance) /= size(quality%raw, dim=1) ) &
                THROW_HARD('evaluate_cavg_quality: invalid histogram variance size')
            quality%raw(:, I_LOG_HIST_VARIANCE) = hist_variance
            deallocate(hist_variance)
        endif
        call calc_power_spectrum_quality_signal(imgs, mskdiam, quality%hard_reject, quality%spec_dmat)
        call normalize_cavg_quality_features(quality%raw, quality%hard_reject, quality%features)
        call model%classify(quality)
    end subroutine evaluate_cavg_quality

    subroutine write_cavg_quality_analysis( quality, reference_states, model, fname, dataset_id )
        type(cavg_quality_result), intent(in) :: quality
        integer,                   intent(in) :: reference_states(:)
        type(cavg_quality_model),  intent(in) :: model
        character(len=*),          intent(in) :: fname
        character(len=*),          intent(in) :: dataset_id
        logical, allocatable :: auto(:), ref(:)
        real,    allocatable :: suggested_weights(:)
        character(len=LONGSTRLEN) :: dataset
        integer :: ncls, ifeat, icls, tp, fp, tn, fn, ngood, nbad, funit
        real    :: precision, recall, specificity, f1, balacc, accuracy
        real    :: best_bal_thr, best_balacc, best_f1_thr, best_f1
        real    :: auc, med_good, med_bad, mad_good, mad_bad, sep
        ncls = size(reference_states)
        if( size(quality%states) /= ncls ) THROW_HARD('write_cavg_quality_analysis: state size mismatch')
        if( size(quality%scores) /= ncls ) THROW_HARD('write_cavg_quality_analysis: score size mismatch')
        if( size(quality%features, dim=1) /= ncls ) THROW_HARD('write_cavg_quality_analysis: feature size mismatch')
        dataset = trim(dataset_id)
        ref  = reference_states > 0
        auto = quality%states > 0
        ngood = count(ref)
        nbad  = ncls - ngood
        call calc_confusion(auto, ref, tp, fp, tn, fn)
        call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
        call calc_best_thresholds(quality%scores, reference_states, best_bal_thr, best_balacc, best_f1_thr, best_f1)
        allocate(suggested_weights(CAVG_QUALITY_NFEATS), source=0.0)
        do ifeat = 1, CAVG_QUALITY_NFEATS
            suggested_weights(ifeat) = max(0.0, auc_for_values(quality%features(:,ifeat), reference_states) - 0.5)
        end do
        if( sum(suggested_weights) > EPS ) then
            suggested_weights = suggested_weights / sum(suggested_weights)
        else
            suggested_weights = model%weights
        end if
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# cluster_cavgs_quality_analysis_version=1'
        write(funit,'(A,A)') '# dataset_id=', trim(dataset)
        write(funit,'(A,A)') '# model_name=', trim(model%name)
        write(funit,'(A,A)') '# model_family=', trim(model%family)
        write(funit,'(A,A)') '# model_context=', trim(model%context)
        write(funit,'(A,I0)') '# n_classes=', ncls
        write(funit,'(A,I0)') '# manual_selected=', ngood
        write(funit,'(A,I0)') '# manual_rejected=', nbad
        write(funit,'(A,I0)') '# auto_selected=', count(auto)
        write(funit,'(A,I0)') '# true_positive=', tp
        write(funit,'(A,I0)') '# false_positive=', fp
        write(funit,'(A,I0)') '# true_negative=', tn
        write(funit,'(A,I0)') '# false_negative=', fn
        write(funit,'(A,F10.5)') '# precision=', precision
        write(funit,'(A,F10.5)') '# recall=', recall
        write(funit,'(A,F10.5)') '# specificity=', specificity
        write(funit,'(A,F10.5)') '# f1=', f1
        write(funit,'(A,F10.5)') '# balanced_accuracy=', balacc
        write(funit,'(A,F10.5)') '# accuracy=', accuracy
        write(funit,'(A,F10.5)') '# score_auc=', auc_for_values(quality%scores, reference_states)
        write(funit,'(A,F10.5)') '# raw_score_threshold=', quality%raw_threshold
        write(funit,'(A,F10.5)') '# threshold_boundary_margin=', quality%threshold_margin
        write(funit,'(A,F10.5)') '# current_score_threshold=', quality%threshold
        write(funit,'(A,F10.5)') '# best_balacc_threshold=', best_bal_thr
        write(funit,'(A,F10.5)') '# best_balacc=', best_balacc
        write(funit,'(A,F10.5)') '# best_f1_threshold=', best_f1_thr
        write(funit,'(A,F10.5)') '# best_f1=', best_f1
        call write_model_as_analysis_comments(funit, model)
        call write_cavg_quality_feature_inventory(funit)
        call write_hist_dmat_comments(funit, quality)
        call write_spec_dmat_comments(funit, quality)
        write(funit,'(A)') '# feature_summary_header=feature,auc,median_manual_good,median_manual_bad,robust_separation,current_weight,suggested_weight'
        do ifeat = 1, CAVG_QUALITY_NFEATS
            auc      = auc_for_values(quality%features(:,ifeat), reference_states)
            med_good = median_by_state(quality%features(:,ifeat), ref)
            med_bad  = median_by_state(quality%features(:,ifeat), .not. ref)
            mad_good = mad_by_state(quality%features(:,ifeat), ref, med_good)
            mad_bad  = mad_by_state(quality%features(:,ifeat), .not. ref, med_bad)
            sep      = safe_div(med_good - med_bad, 0.5 * (mad_good + mad_bad) + EPS)
            write(funit,'(A,A,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5)') &
                '# feature_summary,', trim(cavg_quality_feature_name(ifeat)), ',', auc, ',', med_good, ',', med_bad, &
                ',', sep, ',', model%weights(ifeat), ',', suggested_weights(ifeat)
        end do
        call write_threshold_scan_comments(funit, quality%scores, reference_states)
        call write_analysis_class_header(funit)
        do icls = 1, ncls
            call write_analysis_class_row(funit, trim(dataset), model, quality, reference_states, icls)
        end do
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
        write(logfhandle,'(A,F8.3,A,F8.3)') '>>> BEST QUALITY THRESHOLD BALACC/F1: ', best_bal_thr, ' / ', best_f1_thr
        deallocate(auto, ref, suggested_weights)
    end subroutine write_cavg_quality_analysis

    subroutine write_cavg_quality_feature_table( quality, model, fname, dataset_id, manual_states )
        type(cavg_quality_result), intent(in) :: quality
        type(cavg_quality_model),  intent(in) :: model
        character(len=*),          intent(in) :: fname
        character(len=*),          intent(in) :: dataset_id
        integer, optional,         intent(in) :: manual_states(:)
        integer :: funit, icls, ncls, manual_state, auto_match
        logical :: have_manual
        ncls = size(quality%states)
        if( size(quality%scores) /= ncls ) THROW_HARD('write_cavg_quality_feature_table: score size mismatch')
        if( size(quality%features, dim=1) /= ncls ) THROW_HARD('write_cavg_quality_feature_table: feature size mismatch')
        have_manual = present(manual_states)
        if( have_manual )then
            if( size(manual_states) /= ncls ) THROW_HARD('write_cavg_quality_feature_table: manual state size mismatch')
        endif
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        call write_analysis_class_header(funit)
        do icls = 1, ncls
            if( have_manual )then
                manual_state = manual_states(icls)
                auto_match   = merge(1, 0, (manual_state > 0) .eqv. (quality%states(icls) > 0))
            else
                manual_state = 0
                auto_match   = 0
            endif
            call write_feature_table_class_row(funit, dataset_id, model, quality, icls, manual_state, auto_match)
        end do
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
    end subroutine write_cavg_quality_feature_table

    subroutine write_model_as_analysis_comments( funit, model )
        integer,                  intent(in) :: funit
        type(cavg_quality_model), intent(in) :: model
        integer :: i
        write(funit,'(A)', advance='no') '# model_feature_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') model%weights(i)
        end do
        write(funit,*)
        write(funit,'(A,ES14.6)') '# model_boundary_margin=', model%boundary_margin
        write(funit,'(A,ES14.6)') '# model_min_score_separation=', model%min_score_separation
        write(funit,'(A,ES14.6)') '# model_hist_dmat_weight=', model%hist_dmat_weight
        write(funit,'(A,ES14.6)') '# model_spec_dmat_weight=', model%spec_dmat_weight
        write(funit,'(A,ES14.6)') '# model_otsu_min_offset=', model%otsu_min_offset
        write(funit,'(A,ES14.6)') '# model_otsu_max_offset=', model%otsu_max_offset
        write(funit,'(A,ES14.6)') '# model_cluster_rescue_margin=', model%cluster_rescue_margin
        write(funit,'(A,ES14.6)') '# model_min_accept_frac=', model%min_accept_frac
        write(funit,'(A,L1)') '# model_use_lowsep_otsu=', model%use_lowsep_otsu
        write(funit,'(A,L1)') '# model_use_otsu_window=', model%use_otsu_window
        write(funit,'(A,L1)') '# model_use_cluster_rescue=', model%use_cluster_rescue
        write(funit,'(A,L1)') '# model_enforce_min_accept_frac=', model%enforce_min_accept_frac
    end subroutine write_model_as_analysis_comments

    subroutine write_hist_dmat_comments( funit, quality )
        integer,                   intent(in) :: funit
        type(cavg_quality_result), intent(in) :: quality
        integer :: i, j, ncls
        if( .not. allocated(quality%hist_dmat) ) return
        ncls = size(quality%hist_dmat, dim=1)
        if( size(quality%hist_dmat, dim=2) /= ncls ) &
            THROW_HARD('write_hist_dmat_comments: histogram distance matrix is not square')
        write(funit,'(A,I0)') '# hist_dmat_nclasses=', ncls
        write(funit,'(A,A)') '# hist_dmat_metric=', CAVG_QUALITY_HIST_DMAT_METRIC
        write(funit,'(A)') '# hist_dmat_header=i,j,distance'
        do i = 1, ncls - 1
            do j = i + 1, ncls
                write(funit,'(A,I0,A,I0,A,ES14.6)') '# hist_dmat,', i, ',', j, ',', quality%hist_dmat(i,j)
            end do
        end do
    end subroutine write_hist_dmat_comments

    subroutine write_spec_dmat_comments( funit, quality )
        integer,                   intent(in) :: funit
        type(cavg_quality_result), intent(in) :: quality
        integer :: i, j, ncls
        if( .not. allocated(quality%spec_dmat) ) return
        ncls = size(quality%spec_dmat, dim=1)
        if( size(quality%spec_dmat, dim=2) /= ncls ) &
            THROW_HARD('write_spec_dmat_comments: spectrum distance matrix is not square')
        write(funit,'(A,I0)') '# spec_dmat_nclasses=', ncls
        write(funit,'(A,A)') '# spec_dmat_metric=', CAVG_QUALITY_SPEC_DMAT_METRIC
        write(funit,'(A)') '# spec_dmat_header=i,j,distance'
        do i = 1, ncls - 1
            do j = i + 1, ncls
                write(funit,'(A,I0,A,I0,A,ES14.6)') '# spec_dmat,', i, ',', j, ',', quality%spec_dmat(i,j)
            end do
        end do
    end subroutine write_spec_dmat_comments

    subroutine write_analysis_class_header( funit )
        integer, intent(in) :: funit
        integer :: ifeat
        write(funit,'(A)', advance='no') 'dataset_id,model_context,model_name,class,state,hard_reject,quality_cluster,quality_score,manual_state,auto_matches_manual'
        do ifeat = 1, CAVG_QUALITY_NFEATS
            write(funit,'(A)', advance='no') ',raw_'//trim(cavg_quality_feature_name(ifeat))// &
                ',z_'//trim(cavg_quality_feature_name(ifeat))
        enddo
        write(funit,*)
    end subroutine write_analysis_class_header

    subroutine write_analysis_class_row( funit, dataset_id, model, quality, manual_states, icls )
        integer,                   intent(in) :: funit, icls
        character(len=*),          intent(in) :: dataset_id
        type(cavg_quality_model),  intent(in) :: model
        type(cavg_quality_result), intent(in) :: quality
        integer,                   intent(in) :: manual_states(:)
        call write_feature_table_class_row(funit, dataset_id, model, quality, icls, manual_states(icls), &
            merge(1, 0, (manual_states(icls) > 0) .eqv. (quality%states(icls) > 0)))
    end subroutine write_analysis_class_row

    subroutine write_feature_table_class_row( funit, dataset_id, model, quality, icls, manual_state, auto_match )
        integer,                   intent(in) :: funit, icls, manual_state, auto_match
        character(len=*),          intent(in) :: dataset_id
        type(cavg_quality_model),  intent(in) :: model
        type(cavg_quality_result), intent(in) :: quality
        integer :: ifeat
        write(funit,'(A,A,A,A,A,A,I0,A,I0,A,L1,A,I0,A,ES14.6,A,I0,A,I0)', advance='no') &
            trim(dataset_id), ',', trim(model%context), ',', trim(model%name), ',', icls, ',', quality%states(icls), ',', &
            quality%hard_reject(icls), ',', quality%labels(icls), ',', quality%scores(icls), ',', manual_state, ',', &
            auto_match
        do ifeat = 1, CAVG_QUALITY_NFEATS
            write(funit,'(A,ES14.6,A,ES14.6)', advance='no') ',', quality%raw(icls,ifeat), ',', quality%features(icls,ifeat)
        enddo
        write(funit,*)
    end subroutine write_feature_table_class_row

    subroutine write_threshold_scan_comments( funit, scores, reference_states )
        integer, intent(in) :: funit
        real,    intent(in) :: scores(:)
        integer, intent(in) :: reference_states(:)
        logical, allocatable :: ref(:), pred(:)
        integer :: i, tp, fp, tn, fn
        real    :: threshold, precision, recall, specificity, f1, balacc, accuracy
        allocate(ref(size(reference_states)), source=(reference_states > 0))
        allocate(pred(size(scores)), source=.false.)
        write(funit,'(A)') '# threshold_scan_header=score_threshold,selected,tp,fp,tn,fn,precision,recall,specificity,f1,balanced_accuracy,accuracy'
        call write_one_threshold(minval(scores) - EPS)
        do i = 1, size(scores)
            call write_one_threshold(scores(i))
        end do
        call write_one_threshold(maxval(scores) + EPS)
        deallocate(ref, pred)

    contains

        subroutine write_one_threshold( thresh )
            real, intent(in) :: thresh
            threshold = thresh
            pred = scores >= threshold
            call calc_confusion(pred, ref, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
            write(funit,'(A,F12.6,A,I0,A,I0,A,I0,A,I0,A,I0,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5)') &
                '# threshold_scan,', threshold, ',', count(pred), ',', tp, ',', fp, ',', tn, ',', fn, ',', precision, ',', &
                recall, ',', specificity, ',', f1, ',', balacc, ',', accuracy
        end subroutine write_one_threshold

    end subroutine write_threshold_scan_comments

    subroutine calc_best_thresholds( scores, reference_states, best_bal_thr, best_balacc, best_f1_thr, best_f1 )
        real,    intent(in)  :: scores(:)
        integer, intent(in)  :: reference_states(:)
        real,    intent(out) :: best_bal_thr, best_balacc, best_f1_thr, best_f1
        logical, allocatable :: ref(:), pred(:)
        integer :: i, tp, fp, tn, fn
        real    :: threshold, precision, recall, specificity, f1, balacc, accuracy
        allocate(ref(size(reference_states)), source=(reference_states > 0))
        allocate(pred(size(scores)), source=.false.)
        best_bal_thr = 0.0
        best_f1_thr  = 0.0
        best_balacc  = -huge(1.0)
        best_f1      = -huge(1.0)
        call eval_threshold(minval(scores) - EPS)
        do i = 1, size(scores)
            call eval_threshold(scores(i))
        end do
        call eval_threshold(maxval(scores) + EPS)
        deallocate(ref, pred)

    contains

        subroutine eval_threshold( thresh )
            real, intent(in) :: thresh
            threshold = thresh
            pred = scores >= threshold
            call calc_confusion(pred, ref, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
            if( balacc > best_balacc ) then
                best_balacc  = balacc
                best_bal_thr = threshold
            end if
            if( f1 > best_f1 ) then
                best_f1     = f1
                best_f1_thr = threshold
            end if
        end subroutine eval_threshold

    end subroutine calc_best_thresholds

end module simple_cavg_quality_analysis
