!@descr: evaluation and analyze-mode reporting for class-average quality analysis
module simple_cavg_quality_analysis
use simple_defs,               only: logfhandle, LONGSTRLEN
use simple_error,              only: simple_exception
use simple_image,              only: image
use simple_oris,               only: oris
use simple_parameters,         only: parameters
use simple_cavg_quality_feats, only: cavg_quality_feature_name, &
    extract_cavg_quality_features, normalize_cavg_quality_features, &
    write_cavg_quality_feature_inventory
use simple_cavg_quality_model, only: cavg_quality_model
use simple_cavg_quality_relations, only: CAVG_RELATIONAL_FEATURE_NAME, cavg_quality_relation_analysis
use simple_cavg_quality_stats, only: calc_confusion, calc_binary_metrics, auc_for_values, &
    median_by_state, mad_by_state, safe_div
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, &
    CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1, cavg_quality_result
implicit none
private
#include "simple_local_flags.inc"

public :: evaluate_cavg_quality
public :: evaluate_cavg_quality_hard_reject
public :: prepare_cavg_quality
public :: write_cavg_quality_training_table
public :: write_cavg_quality_feature_table

contains

    subroutine prepare_cavg_quality( imgs, cls_oris, mskdiam, quality, quality_context )
        class(image),              intent(inout) :: imgs(:)
        type(oris),                intent(in)    :: cls_oris
        real,                      intent(in)    :: mskdiam
        type(cavg_quality_result), intent(inout) :: quality
        character(len=*), optional,intent(in)    :: quality_context
        character(len=32) :: context
        call quality%kill()
        context = 'chunk'
        if( present(quality_context) ) context = trim(quality_context)
        call extract_cavg_quality_features(imgs, cls_oris, mskdiam, quality, trim(context))
        call normalize_cavg_quality_features(quality%raw, quality%hard_reject, quality%features)
    end subroutine prepare_cavg_quality

    subroutine evaluate_cavg_quality( imgs, cls_oris, mskdiam, quality, model, quality_context, relation_params, &
                                      relation_result )
        class(image),              intent(inout) :: imgs(:)
        type(oris),                intent(in)    :: cls_oris
        real,                      intent(in)    :: mskdiam
        type(cavg_quality_result), intent(inout) :: quality
        type(cavg_quality_model),  intent(in)    :: model
        character(len=*), optional,intent(in)    :: quality_context
        type(parameters), optional,intent(in)    :: relation_params
        type(cavg_quality_relation_analysis), optional, intent(inout) :: relation_result
        type(cavg_quality_relation_analysis) :: local_relation
        character(len=32) :: context
        context = trim(model%context)
        if( present(quality_context) ) context = trim(quality_context)
        call prepare_cavg_quality(imgs, cls_oris, mskdiam, quality, trim(context))
        if( model%supports_relational() )then
            if( .not. present(relation_params) ) &
                THROW_HARD('evaluate_cavg_quality: relational model requires parameters')
            if( present(relation_result) )then
                call relation_result%calculate(imgs, quality, relation_params, model%relational_knn, &
                    model%relational_corr_hp, model%relational_corr_lp, model%relational_corr_trs)
                call model%classify(quality, relation_result%feature)
            else
                call local_relation%calculate(imgs, quality, relation_params, model%relational_knn, &
                    model%relational_corr_hp, model%relational_corr_lp, model%relational_corr_trs)
                call model%classify(quality, local_relation%feature)
                call local_relation%kill()
            end if
        else
            call model%classify(quality)
        end if
    end subroutine evaluate_cavg_quality

    subroutine evaluate_cavg_quality_hard_reject( imgs, cls_oris, mskdiam, quality, quality_context )
        class(image),              intent(inout) :: imgs(:)
        type(oris),                intent(in)    :: cls_oris
        real,                      intent(in)    :: mskdiam
        type(cavg_quality_result), intent(inout) :: quality
        character(len=*), optional,intent(in)    :: quality_context
        logical,                     allocatable :: standard_hard_reject(:)
        integer                                  :: icls, ncls
        call quality%kill()
        call extract_cavg_quality_features(imgs, cls_oris, mskdiam, quality, quality_context)
        call normalize_cavg_quality_features(quality%raw, quality%hard_reject, quality%features)
        if( .not. allocated(quality%features)    ) THROW_HARD('evaluate_cavg_quality_hard_reject: missing features')
        if( .not. allocated(quality%hard_reject) ) &
            THROW_HARD('evaluate_cavg_quality_hard_reject: missing hard-reject mask')
        ncls = size(quality%features, dim=1)
        if( size(quality%features, dim=2) /= CAVG_QUALITY_NFEATS ) &
            THROW_HARD('evaluate_cavg_quality_hard_reject: invalid feature count')
        if( size(quality%hard_reject) /= ncls ) THROW_HARD('evaluate_cavg_quality_hard_reject: invalid mask size')
        allocate(standard_hard_reject(ncls), source=quality%hard_reject)
        if( allocated(quality%states)  ) deallocate(quality%states)
        if( allocated(quality%labels)  ) deallocate(quality%labels)
        if( allocated(quality%medoids) ) deallocate(quality%medoids)
        if( allocated(quality%scores)  ) deallocate(quality%scores)
        allocate(quality%states(ncls), quality%labels(ncls), source=0)
        allocate(quality%scores(ncls), source=-CLIP_Z)
        quality%threshold        = 0.0
        quality%raw_threshold    = 0.0
        quality%threshold_offset = 0.0
        quality%separation       = 0.0
        quality%nclust           = 2
        quality%good_label       = 1
        quality%used_threshold   = .false.
        quality%model_name       = 'standard_hard_gates'
        quality%soft_decision    = 'hard_only'
        quality%soft_reason      = 'standard_hard_gates_only'
        do icls = 1, ncls
            if( standard_hard_reject(icls) ) cycle
            quality%states(icls) = 1
            quality%labels(icls) = 1
            quality%scores(icls) = 1.0
        end do
        deallocate(standard_hard_reject)
    end subroutine evaluate_cavg_quality_hard_reject

    subroutine write_cavg_quality_training_table( quality, reference_states, model, fname, dataset_id, relation )
        type(cavg_quality_result), intent(in) :: quality
        integer,                   intent(in) :: reference_states(:)
        type(cavg_quality_model),  intent(in) :: model
        character(len=*),          intent(in) :: fname
        character(len=*),          intent(in) :: dataset_id
        type(cavg_quality_relation_analysis), intent(in) :: relation
        logical, allocatable :: auto(:), ref(:)
        real,    allocatable :: suggested_weights(:)
        character(len=LONGSTRLEN) :: dataset
        integer :: ncls, ifeat, icls, tp, fp, tn, fn, ngood, nbad, funit
        real    :: precision, recall, specificity, f1, balacc, accuracy
        real    :: best_bal_thr, best_balacc, best_f1_thr, best_f1
        real    :: auc, med_good, med_bad, mad_good, mad_bad, sep
        ncls = size(reference_states)
        if( size(quality%states) /= ncls ) THROW_HARD('write_cavg_quality_training_table: state size mismatch')
        if( size(quality%scores) /= ncls ) THROW_HARD('write_cavg_quality_training_table: score size mismatch')
        if( size(quality%features, dim=1) /= ncls ) &
            THROW_HARD('write_cavg_quality_training_table: feature size mismatch')
        if( .not. allocated(relation%raw) .or. .not. allocated(relation%feature) ) &
            THROW_HARD('write_cavg_quality_training_table: missing relational feature')
        if( size(relation%raw) /= ncls .or. size(relation%feature) /= ncls ) &
            THROW_HARD('write_cavg_quality_training_table: relational feature size mismatch')
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
        write(funit,'(A)') '# cavg_quality_training_version=1'
        write(funit,'(A,A)') '# dataset_id=', trim(dataset)
        write(funit,'(A,A)') '# model_name=', trim(model%name)
        write(funit,'(A,A)') '# relational_feature_schema=', CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1
        write(funit,'(A)') '# relation_anchor=rotmax_cc'
        write(funit,'(A,A)') '# relational_feature_name=', CAVG_RELATIONAL_FEATURE_NAME
        write(funit,'(A,I0)') '# relational_knn=', relation%k_requested
        write(funit,'(A,F10.4)') '# relational_corr_hp=', relation%hp
        write(funit,'(A,F10.4)') '# relational_corr_lp=', relation%lp
        write(funit,'(A,F10.4)') '# relational_corr_trs=', relation%trs
        write(funit,'(A)') '# signal_stats_provider=cluster_cavgs_signal_stats_v1'
        write(funit,'(A,I0)') '# signal_stats_histogram_bins=', relation%signal%histogram_bins
        write(funit,'(A,F10.4)') '# signal_stats_power_hp=', relation%signal%power_hp
        write(funit,'(A,F10.4)') '# signal_stats_power_lp=', relation%signal%power_lp
        write(funit,'(A,F10.4)') '# signal_stats_mask_radius_pixels=', relation%signal_msk
        write(funit,'(A,F14.6)') '# signal_stats_oa_min=', relation%signal_oa_minmax(1)
        write(funit,'(A,F14.6)') '# signal_stats_oa_max=', relation%signal_oa_minmax(2)
        write(funit,'(A,I0)') '# n_classes=', ncls
        write(funit,'(A,I0)') '# manual_selected=', ngood
        write(funit,'(A,I0)') '# manual_rejected=', nbad
        write(funit,'(A,I0)') '# auto_selected=', count(auto)
        write(funit,'(A,A)') '# soft_decision=', trim(quality%soft_decision)
        write(funit,'(A,A)') '# soft_reason=', trim(quality%soft_reason)
        write(funit,'(A,I0)') '# soft_trainable_rejected=', count(quality%states <= 0 .and. .not. quality%hard_reject)
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
        write(funit,'(A,F10.5)') '# threshold_offset=', quality%threshold_offset
        write(funit,'(A,F10.5)') '# current_score_threshold=', quality%threshold
        write(funit,'(A,F10.5)') '# best_balacc_threshold=', best_bal_thr
        write(funit,'(A,F10.5)') '# best_balacc=', best_balacc
        write(funit,'(A,F10.5)') '# best_f1_threshold=', best_f1_thr
        write(funit,'(A,F10.5)') '# best_f1=', best_f1
        call write_model_as_analysis_comments(funit, model)
        call write_cavg_quality_feature_inventory(funit)
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
        call write_analysis_class_header(funit, .true.)
        do icls = 1, ncls
            call write_analysis_class_row(funit, trim(dataset), model, quality, reference_states, icls, &
                relation%raw(icls), relation%feature(icls))
        end do
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
        write(logfhandle,'(A,F8.3,A,F8.3)') '>>> BEST QUALITY THRESHOLD BALACC/F1: ', best_bal_thr, ' / ', best_f1_thr
        deallocate(auto, ref, suggested_weights)
    end subroutine write_cavg_quality_training_table

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
                ! No reference is available; keep auto_matches_manual from
                ! masquerading as a real mismatch.
                auto_match   = -1
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
        write(funit,'(A,A)') '# model_context=', trim(model%context)
        write(funit,'(A,A)') '# model_feature_policy=', trim(model%feature_policy)
        write(funit,'(A)', advance='no') '# model_feature_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') model%weights(i)
        end do
        write(funit,*)
        write(funit,'(A,ES14.6)') '# model_boundary_margin=', model%boundary_margin
        write(funit,'(A,ES14.6)') '# model_min_score_separation=', model%min_score_separation
        write(funit,'(A,ES14.6)') '# model_otsu_min_offset=', model%otsu_min_offset
        write(funit,'(A,ES14.6)') '# model_otsu_max_offset=', model%otsu_max_offset
        write(funit,'(A,ES14.6)') '# model_cluster_rescue_margin=', model%cluster_rescue_margin
        write(funit,'(A,ES14.6)') '# model_min_accept_frac=', model%min_accept_frac
        write(funit,'(A,L1)') '# model_use_lowsep_otsu=', model%use_lowsep_otsu
        write(funit,'(A,L1)') '# model_use_otsu_window=', model%use_otsu_window
        write(funit,'(A,L1)') '# model_use_cluster_rescue=', model%use_cluster_rescue
        write(funit,'(A,L1)') '# model_enforce_min_accept_frac=', model%enforce_min_accept_frac
        write(funit,'(A,ES14.6)') '# model_intercept=', model%intercept
        write(funit,'(A)', advance='no') '# model_linear_coefficients='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') model%linear_coefficients(i)
        end do
        write(funit,*)
        write(funit,'(A)', advance='no') '# model_interaction_terms='
        do i = 1, model%n_interactions
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(A,A,A)', advance='no') &
                trim(cavg_quality_feature_name(model%interaction_terms(i,1))), ':', &
                trim(cavg_quality_feature_name(model%interaction_terms(i,2)))
        end do
        write(funit,*)
        write(funit,'(A)', advance='no') '# model_interaction_coefficients='
        do i = 1, model%n_interactions
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') model%interaction_coefficients(i)
        end do
        write(funit,*)
        write(funit,'(A,ES14.6)') '# model_prob_threshold=', model%prob_threshold
        write(funit,'(A,ES14.6)') '# model_regularization_lambda=', model%regularization_lambda
        write(funit,'(A,ES14.6)') '# model_calibration_temperature=', model%calibration_temperature
        write(funit,'(A,A)') '# model_relational_feature_schema=', trim(model%relational_feature_schema)
        write(funit,'(A,ES14.6)') '# model_relational_coefficient=', model%relational_coefficient
    end subroutine write_model_as_analysis_comments

    subroutine write_analysis_class_header( funit, include_relation )
        integer, intent(in) :: funit
        logical, optional, intent(in) :: include_relation
        integer :: ifeat
        logical :: with_relation
        with_relation = .false.
        if( present(include_relation) ) with_relation = include_relation
        write(funit,'(A)', advance='no') &
            'dataset_id,model_name,class,state,hard_reject,quality_cluster,quality_score,'//&
            'manual_state,auto_matches_manual'
        do ifeat = 1, CAVG_QUALITY_NFEATS
            write(funit,'(A)', advance='no') ',raw_'//trim(cavg_quality_feature_name(ifeat))// &
                ',z_'//trim(cavg_quality_feature_name(ifeat))
        enddo
        if( with_relation ) write(funit,'(A)', advance='no') ',raw_'//CAVG_RELATIONAL_FEATURE_NAME// &
            ',z_'//CAVG_RELATIONAL_FEATURE_NAME
        write(funit,*)
    end subroutine write_analysis_class_header

    subroutine write_analysis_class_row( funit, dataset_id, model, quality, manual_states, icls, &
                                         relational_raw, relational_feature )
        integer,                   intent(in) :: funit, icls
        character(len=*),          intent(in) :: dataset_id
        type(cavg_quality_model),  intent(in) :: model
        type(cavg_quality_result), intent(in) :: quality
        integer,                   intent(in) :: manual_states(:)
        real, optional,            intent(in) :: relational_raw, relational_feature
        if( present(relational_raw) .and. present(relational_feature) )then
            call write_feature_table_class_row(funit, dataset_id, model, quality, icls, manual_states(icls), &
                merge(1, 0, (manual_states(icls) > 0) .eqv. (quality%states(icls) > 0)), &
                relational_raw, relational_feature)
        else
            call write_feature_table_class_row(funit, dataset_id, model, quality, icls, manual_states(icls), &
                merge(1, 0, (manual_states(icls) > 0) .eqv. (quality%states(icls) > 0)))
        end if
    end subroutine write_analysis_class_row

    subroutine write_feature_table_class_row( funit, dataset_id, model, quality, icls, manual_state, auto_match, &
                                              relational_raw, relational_feature )
        integer,                   intent(in) :: funit, icls, manual_state, auto_match
        character(len=*),          intent(in) :: dataset_id
        type(cavg_quality_model),  intent(in) :: model
        type(cavg_quality_result), intent(in) :: quality
        real, optional,            intent(in) :: relational_raw, relational_feature
        integer :: ifeat
        write(funit,'(A,A,A,A,I0,A,I0,A,L1,A,I0,A,ES14.6,A,I0,A,I0)', advance='no') &
            trim(dataset_id), ',', trim(model%name), ',', &
            icls, ',', quality%states(icls), ',', &
            quality%hard_reject(icls), ',', quality%labels(icls), ',', quality%scores(icls), ',', manual_state, ',', &
            auto_match
        do ifeat = 1, CAVG_QUALITY_NFEATS
            write(funit,'(A,ES14.6,A,ES14.6)', advance='no') ',', quality%raw(icls,ifeat), ',', quality%features(icls,ifeat)
        enddo
        if( present(relational_raw) .and. present(relational_feature) ) &
            write(funit,'(A,ES14.6,A,ES14.6)', advance='no') ',', relational_raw, ',', relational_feature
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
