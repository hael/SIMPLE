!@descr: learn-mode training-table reader and model search for class-average quality
module simple_cavg_quality_learn
use simple_defs,               only: logfhandle, LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string,             only: string
use simple_string_utils,       only: str2int, str2real, str_is_true, csv_field
use simple_cavg_quality_feats, only: cavg_quality_feature_name, CAVG_RES_HARD_REJECT_A, I_LOG_POP, I_NEG_LOG_RES, I_MASK_INSIDE, &
    I_LOCVAR_FG, I_LOCVAR_BG, I_CC_SINGLE, I_CORR_FRC, I_HIST_ENTROPY, I_CC_AREA_FRAC, I_PRESENCE, &
    I_LOG_CONTRAST, I_LOG_HIST_VARIANCE
use simple_cavg_quality_model, only: cavg_quality_model
use simple_cavg_quality_stats, only: calc_confusion, calc_binary_metrics, auc_for_values
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, cavg_quality_model_spec, &
    cavg_quality_result, cavg_quality_training_dataset, cavg_quality_learn_diagnostics
implicit none
private
#include "simple_local_flags.inc"

public :: learn_cavg_quality_model

integer, parameter :: CAVG_QUALITY_LEARN_TOP_K = 10
integer, parameter :: CAVG_QUALITY_LEARN_N_POLICIES = 9
real, parameter :: LEARN_MINSEPS(5)       = [0.05, 0.10, 0.15, 0.20, 0.30]
real, parameter :: LEARN_MARGINS(11)      = [-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, &
                                              -0.05, 0.0, 0.05, 0.10, 0.20]
real, parameter :: LEARN_POOL_FRACS(5)    = [0.50, 0.60, 0.65, 0.70, 0.80]

contains

    subroutine learn_cavg_quality_model( analysis_files, learned_model, model_fname, report_fname )
        class(string),             intent(in)    :: analysis_files(:)
        type(cavg_quality_model),  intent(inout) :: learned_model
        character(len=*),          intent(in)    :: model_fname, report_fname
        type(cavg_quality_training_dataset), allocatable :: dsets(:)
        type(cavg_quality_model)      :: base_model, candidate
        type(cavg_quality_model_spec) :: base_spec, candidate_spec, best_spec
        type(cavg_quality_model_spec) :: top_specs(CAVG_QUALITY_LEARN_TOP_K)
        type(cavg_quality_model_spec), allocatable :: best_tie_specs(:)
        real :: suggested_weights(CAVG_QUALITY_NFEATS)
        real :: top_scores(CAVG_QUALITY_LEARN_TOP_K)
        real :: score, best_score
        integer :: i, ipol, im, isep, ifrac, max_grid, n_grid, n_top, n_best_ties
        if( size(analysis_files) == 0 ) THROW_HARD('learn_cavg_quality_model: empty analysis file table')
        allocate(dsets(size(analysis_files)))
        do i = 1, size(analysis_files)
            call read_quality_training_dataset(analysis_files(i)%to_char(), dsets(i))
        end do
        base_spec = learned_model%get_spec()
        call base_model%init_spec(base_spec)
        call calc_suggested_training_weights(dsets, suggested_weights)
        best_spec  = base_spec
        best_score = -huge(1.0)
        top_scores = -huge(1.0)
        n_top       = 0
        n_grid      = 0
        n_best_ties = 0
        max_grid = CAVG_QUALITY_LEARN_N_POLICIES * size(LEARN_MINSEPS) * size(LEARN_MARGINS)
        if( model_is_pool_context(learned_model) ) max_grid = max_grid * size(LEARN_POOL_FRACS)
        allocate(best_tie_specs(max_grid))
        do ipol = 1, CAVG_QUALITY_LEARN_N_POLICIES
            candidate_spec = base_spec
            candidate_spec%feature_policy = feature_policy_name(ipol)
            candidate_spec%weights = suggested_weights
            call apply_feature_policy(ipol, candidate_spec%weights)
            do isep = 1, size(LEARN_MINSEPS)
                candidate_spec%min_score_separation = LEARN_MINSEPS(isep)
                do im = 1, size(LEARN_MARGINS)
                    candidate_spec%boundary_margin = LEARN_MARGINS(im)
                    if( model_is_pool_context(learned_model) )then
                        do ifrac = 1, size(LEARN_POOL_FRACS)
                            candidate_spec%min_accept_frac = LEARN_POOL_FRACS(ifrac)
                            call candidate%init_spec(candidate_spec)
                            score = macro_balacc_for_model(dsets, candidate)
                            n_grid = n_grid + 1
                            call consider_model_candidate(candidate%get_spec(), score, best_spec, best_score, &
                                best_tie_specs, n_best_ties, top_specs, top_scores, n_top)
                        end do
                    else
                        call candidate%init_spec(candidate_spec)
                        score = macro_balacc_for_model(dsets, candidate)
                        n_grid = n_grid + 1
                        call consider_model_candidate(candidate%get_spec(), score, best_spec, best_score, &
                            best_tie_specs, n_best_ties, top_specs, top_scores, n_top)
                    endif
                end do
            end do
        end do
        best_spec%name = trim(best_spec%context)//'_learned_v1'
        call learned_model%init_spec(best_spec)
        call learned_model%write(model_fname)
        call write_cavg_quality_learn_report(report_fname, dsets, base_model, suggested_weights, learned_model, best_score, &
            n_grid, top_specs, top_scores, n_top, best_tie_specs, n_best_ties)
        call kill_training_datasets(dsets)
        deallocate(best_tie_specs)
        deallocate(dsets)
    end subroutine learn_cavg_quality_model

    logical function model_is_pool_context( model )
        type(cavg_quality_model), intent(in) :: model
        model_is_pool_context = trim(model%context) == 'pool'
    end function model_is_pool_context

    function feature_policy_name( ipolicy ) result( name )
        integer, intent(in) :: ipolicy
        character(len=32) :: name
        select case(ipolicy)
            case(1)
                name = 'base12'
            case(2)
                name = 'base12_no_geom_softs'
            case(3)
                name = 'base12_pruned_plus_histvar'
            case(4)
                name = 'base12_no_cc_area_plus_histvar'
            case(5)
                name = 'base12_no_cc_area'
            case(6)
                name = 'base12_plus_histvar'
            case(7)
                name = 'all15_no_geom_softs'
            case(8)
                name = 'all_features_no_mask_single'
            case(9)
                name = 'all_features'
            case default
                THROW_HARD('feature_policy_name: invalid feature policy')
        end select
    end function feature_policy_name

    subroutine feature_policy_mask( ipolicy, mask )
        integer, intent(in)  :: ipolicy
        logical, intent(out) :: mask(CAVG_QUALITY_NFEATS)
        mask = .false.
        select case(ipolicy)
            case(1)
                mask(1:12) = .true.
            case(2)
                mask(1:12) = .true.
                mask(I_MASK_INSIDE) = .false.
                mask(I_CC_SINGLE)   = .false.
            case(3)
                mask(1:12) = .true.
                mask(I_MASK_INSIDE) = .false.
                mask(I_CC_SINGLE)   = .false.
                mask(I_CC_AREA_FRAC) = .false.
                mask(I_LOG_HIST_VARIANCE) = .true.
            case(4)
                mask(1:12) = .true.
                mask(I_CC_AREA_FRAC) = .false.
                mask(I_LOG_HIST_VARIANCE) = .true.
            case(5)
                mask(1:12) = .true.
                mask(I_CC_AREA_FRAC) = .false.
            case(6)
                mask(1:12) = .true.
                mask(I_LOG_HIST_VARIANCE) = .true.
            case(7)
                mask = .true.
                mask(I_MASK_INSIDE) = .false.
                mask(I_CC_SINGLE)   = .false.
                mask(I_CC_AREA_FRAC) = .false.
            case(8)
                mask = .true.
                mask(I_MASK_INSIDE) = .false.
                mask(I_CC_SINGLE)   = .false.
            case(9)
                mask = .true.
            case default
                THROW_HARD('feature_policy_mask: invalid feature policy')
        end select
    end subroutine feature_policy_mask

    subroutine apply_feature_policy( ipolicy, weights )
        integer, intent(in)    :: ipolicy
        real,    intent(inout) :: weights(CAVG_QUALITY_NFEATS)
        logical :: mask(CAVG_QUALITY_NFEATS)
        integer :: n_active
        call feature_policy_mask(ipolicy, mask)
        where( .not. mask ) weights = 0.0
        if( sum(weights) > EPS )then
            weights = weights / sum(weights)
        else
            n_active = count(mask)
            if( n_active < 1 ) THROW_HARD('apply_feature_policy: empty feature policy')
            where( mask )
                weights = 1.0 / real(n_active)
            elsewhere
                weights = 0.0
            end where
        endif
    end subroutine apply_feature_policy

    subroutine feature_policy_indices( ipolicy, inds, ninds )
        integer, intent(in)  :: ipolicy
        integer, intent(out) :: inds(CAVG_QUALITY_NFEATS)
        integer, intent(out) :: ninds
        logical :: mask(CAVG_QUALITY_NFEATS)
        integer :: ifeat
        call feature_policy_mask(ipolicy, mask)
        inds = 0
        ninds = 0
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( .not. mask(ifeat) ) cycle
            ninds = ninds + 1
            inds(ninds) = ifeat
        end do
    end subroutine feature_policy_indices

    subroutine read_quality_training_dataset( fname, dset )
        character(len=*),                    intent(in)    :: fname
        type(cavg_quality_training_dataset), intent(inout) :: dset
        character(len=XLONGSTRLEN) :: line
        character(len=LONGSTRLEN)  :: field
        integer :: funit, ios, nrows, irow, ifeat, feat_col(CAVG_QUALITY_NFEATS)
        integer :: raw_mask_inside_col, raw_centered_col, raw_neg_log_res_col
        logical :: have_feature_header
        dset%fname = trim(fname)
        open(newunit=funit, file=trim(fname), status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('read_quality_training_dataset: failed to open '//trim(fname))
        nrows = 0
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( is_analysis_data_line(line) ) nrows = nrows + 1
        end do
        if( nrows == 0 ) THROW_HARD('read_quality_training_dataset: no class rows in '//trim(fname))
        allocate(dset%features(nrows, CAVG_QUALITY_NFEATS), source=0.0)
        allocate(dset%manual_states(nrows), source=0)
        allocate(dset%hard_reject(nrows), source=.false.)
        rewind(funit)
        irow = 0
        feat_col = 0
        raw_mask_inside_col = 0
        raw_centered_col    = 0
        raw_neg_log_res_col = 0
        have_feature_header = .false.
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( is_analysis_header_line(line) )then
                call map_analysis_feature_columns(line, feat_col)
                call map_analysis_raw_hard_reject_columns(line, raw_mask_inside_col, raw_centered_col, raw_neg_log_res_col)
                call require_analysis_feature_columns(feat_col, trim(fname))
                have_feature_header = .true.
                cycle
            endif
            if( .not. is_analysis_data_line(line) ) cycle
            irow = irow + 1
            field = csv_field(line, 1)
            if( irow == 1 ) dset%dataset_id = trim(field)
            field = csv_field(line, 6)
            dset%hard_reject(irow) = str_is_true(field)
            if( analysis_row_legacy_hard_reject(line, raw_mask_inside_col, raw_centered_col, raw_neg_log_res_col) ) &
                dset%hard_reject(irow) = .true.
            field = csv_field(line, 9)
            dset%manual_states(irow) = str2int(trim(field))
            do ifeat = 1, CAVG_QUALITY_NFEATS
                field = ''
                if( feat_col(ifeat) > 0 )then
                    field = csv_field(line, feat_col(ifeat))
                else if( .not. have_feature_header )then
                    field = csv_field(line, 10 + 2 * ifeat)
                endif
                if( len_trim(field) > 0 ) dset%features(irow, ifeat) = str2real(trim(field))
            end do
        end do
        dset%ncls = irow
        close(funit)
    end subroutine read_quality_training_dataset

    logical function is_analysis_data_line( line )
        character(len=*), intent(in) :: line
        character(len=XLONGSTRLEN) :: tmp
        tmp = adjustl(line)
        is_analysis_data_line = .false.
        if( len_trim(tmp) == 0 ) return
        if( tmp(1:1) == '#' ) return
        if( index(tmp, 'dataset_id,') == 1 ) return
        is_analysis_data_line = .true.
    end function is_analysis_data_line

    logical function is_analysis_header_line( line )
        character(len=*), intent(in) :: line
        character(len=XLONGSTRLEN) :: tmp
        tmp = adjustl(line)
        is_analysis_header_line = index(tmp, 'dataset_id,') == 1
    end function is_analysis_header_line

    subroutine map_analysis_feature_columns( line, feat_col )
        character(len=*), intent(in)  :: line
        integer,          intent(out) :: feat_col(CAVG_QUALITY_NFEATS)
        character(len=LONGSTRLEN) :: field, name
        integer :: icol, ifeat
        feat_col = 0
        do icol = 1, 512
            field = csv_field(line, icol)
            if( len_trim(field) == 0 ) exit
            if( len_trim(field) < 3 ) cycle
            if( field(1:2) /= 'z_' ) cycle
            name = adjustl(trim(field(3:)))
            do ifeat = 1, CAVG_QUALITY_NFEATS
                if( trim(name) == trim(cavg_quality_feature_name(ifeat)) )then
                    feat_col(ifeat) = icol
                    exit
                endif
            end do
        end do
    end subroutine map_analysis_feature_columns

    subroutine map_analysis_raw_hard_reject_columns( line, raw_mask_inside_col, raw_centered_col, raw_neg_log_res_col )
        character(len=*), intent(in)  :: line
        integer,          intent(out) :: raw_mask_inside_col, raw_centered_col, raw_neg_log_res_col
        character(len=LONGSTRLEN) :: field
        integer :: icol
        raw_mask_inside_col = 0
        raw_centered_col    = 0
        raw_neg_log_res_col = 0
        do icol = 1, 512
            field = csv_field(line, icol)
            if( len_trim(field) == 0 ) exit
            select case(trim(field))
                case('raw_neg_log_res')
                    raw_neg_log_res_col = icol
                case('raw_mask_inside')
                    raw_mask_inside_col = icol
                case('raw_centered')
                    raw_centered_col = icol
            end select
        end do
    end subroutine map_analysis_raw_hard_reject_columns

    logical function analysis_row_legacy_hard_reject( line, raw_mask_inside_col, raw_centered_col, raw_neg_log_res_col )
        character(len=*), intent(in) :: line
        integer,          intent(in) :: raw_mask_inside_col, raw_centered_col, raw_neg_log_res_col
        character(len=LONGSTRLEN) :: field
        real :: raw_mask_inside, raw_centered, raw_neg_log_res
        ! Legacy analysis records lack the exact largest-CC outside-pixel count,
        ! but they do carry the raw mask/centroid diagnostics from the same
        ! segmentation pass. They also predate the hard resolution veto. Use
        ! these raw diagnostics to keep obvious validity failures out of the
        ! learned quality model when retraining from existing analysis files.
        analysis_row_legacy_hard_reject = .false.
        if( raw_neg_log_res_col > 0 )then
            field = csv_field(line, raw_neg_log_res_col)
            if( len_trim(field) > 0 )then
                raw_neg_log_res = str2real(trim(field))
                if( raw_neg_log_res < -log(CAVG_RES_HARD_REJECT_A) - EPS ) &
                    analysis_row_legacy_hard_reject = .true.
            endif
        endif
        if( raw_centered_col > 0 )then
            field = csv_field(line, raw_centered_col)
            if( len_trim(field) > 0 )then
                raw_centered = str2real(trim(field))
                if( raw_centered < -1.0 - EPS ) analysis_row_legacy_hard_reject = .true.
            endif
        endif
        if( raw_mask_inside_col > 0 )then
            field = csv_field(line, raw_mask_inside_col)
            if( len_trim(field) > 0 )then
                raw_mask_inside = str2real(trim(field))
                if( raw_mask_inside < -EPS ) analysis_row_legacy_hard_reject = .true.
            endif
        endif
    end function analysis_row_legacy_hard_reject

    subroutine require_analysis_feature_columns( feat_col, fname )
        integer,          intent(in) :: feat_col(CAVG_QUALITY_NFEATS)
        character(len=*), intent(in) :: fname
        character(len=LONGSTRLEN) :: errmsg
        integer :: ifeat
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( feat_col(ifeat) == 0 )then
                errmsg = 'read_quality_training_dataset: missing z_'//trim(cavg_quality_feature_name(ifeat))//&
                    ' column in '//trim(fname)
                THROW_HARD(trim(errmsg))
            endif
        end do
    end subroutine require_analysis_feature_columns

    subroutine calc_suggested_training_weights( dsets, weights )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        real,                                intent(out) :: weights(CAVG_QUALITY_NFEATS)
        real,    allocatable :: vals(:)
        integer, allocatable :: refs(:)
        integer :: nall, ids, j, off, nfit
        nall = 0
        do ids = 1, size(dsets)
            nall = nall + count_trainable_classes(dsets(ids))
        end do
        if( nall == 0 ) THROW_HARD('calc_suggested_training_weights: no non-hard-rejected training rows')
        allocate(vals(nall), source=0.0)
        allocate(refs(nall), source=0)
        weights = 0.0
        do j = 1, CAVG_QUALITY_NFEATS
            off = 0
            do ids = 1, size(dsets)
                nfit = count_trainable_classes(dsets(ids))
                if( nfit == 0 ) cycle
                vals(off+1:off+nfit) = pack(dsets(ids)%features(:,j), .not. dsets(ids)%hard_reject)
                refs(off+1:off+nfit) = pack(dsets(ids)%manual_states, .not. dsets(ids)%hard_reject)
                off = off + nfit
            end do
            weights(j) = max(0.0, auc_for_values(vals, refs) - 0.5)
        end do
        if( sum(weights) > EPS ) weights = weights / sum(weights)
        deallocate(vals, refs)
    end subroutine calc_suggested_training_weights

    integer function count_trainable_classes( dset )
        type(cavg_quality_training_dataset), intent(in) :: dset
        ! Hard rejects are upstream validity-gate failures. They stay visible
        ! in diagnostics, but the learned boundary is fit only on rescuable rows.
        count_trainable_classes = count(.not. dset%hard_reject)
    end function count_trainable_classes

    real function macro_balacc_for_model( dsets, model )
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: model
        integer :: ids, tp, fp, tn, fn, nused
        real :: precision, recall, specificity, f1, balacc, accuracy
        macro_balacc_for_model = 0.0
        nused = 0
        do ids = 1, size(dsets)
            if( count_trainable_classes(dsets(ids)) == 0 ) cycle
            call classify_training_dataset(dsets(ids), model, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
            macro_balacc_for_model = macro_balacc_for_model + balacc
            nused = nused + 1
        end do
        if( nused == 0 ) THROW_HARD('macro_balacc_for_model: no non-hard-rejected training rows')
        macro_balacc_for_model = macro_balacc_for_model / real(nused)
    end function macro_balacc_for_model

    subroutine consider_model_candidate( spec, score, best_spec, best_score, best_tie_specs, n_best_ties, &
                                         top_specs, top_scores, n_top )
        type(cavg_quality_model_spec), intent(in)    :: spec
        real,                          intent(in)    :: score
        type(cavg_quality_model_spec), intent(inout) :: best_spec
        real,                          intent(inout) :: best_score
        type(cavg_quality_model_spec), intent(inout) :: best_tie_specs(:)
        integer,                       intent(inout) :: n_best_ties
        type(cavg_quality_model_spec), intent(inout) :: top_specs(:)
        real,                          intent(inout) :: top_scores(:)
        integer,                       intent(inout) :: n_top
        if( score > best_score + EPS )then
            best_score  = score
            best_spec   = spec
            n_best_ties = 1
            best_tie_specs(1) = spec
        else if( abs(score - best_score) <= EPS )then
            n_best_ties = n_best_ties + 1
            if( n_best_ties <= size(best_tie_specs) ) best_tie_specs(n_best_ties) = spec
        endif
        call record_top_model_candidate(spec, score, top_specs, top_scores, n_top)
    end subroutine consider_model_candidate

    subroutine record_top_model_candidate( spec, score, top_specs, top_scores, n_top )
        type(cavg_quality_model_spec), intent(in)    :: spec
        real,                          intent(in)    :: score
        type(cavg_quality_model_spec), intent(inout) :: top_specs(:)
        real,                          intent(inout) :: top_scores(:)
        integer,                       intent(inout) :: n_top
        integer :: i, pos, max_top
        max_top = min(size(top_specs), size(top_scores))
        if( max_top <= 0 ) return
        if( n_top < max_top )then
            n_top = n_top + 1
            pos   = n_top
        else
            if( score <= top_scores(max_top) + EPS ) return
            pos = max_top
        endif
        do i = pos, 2, -1
            if( score > top_scores(i-1) + EPS )then
                top_scores(i) = top_scores(i-1)
                top_specs(i)  = top_specs(i-1)
                pos = i - 1
            else
                exit
            endif
        end do
        top_scores(pos) = score
        top_specs(pos)  = spec
    end subroutine record_top_model_candidate

    subroutine classify_training_dataset( dset, model, tp, fp, tn, fn )
        type(cavg_quality_training_dataset), intent(in) :: dset
        type(cavg_quality_model),            intent(in) :: model
        integer,                             intent(out):: tp, fp, tn, fn
        type(cavg_quality_result) :: quality
        call classify_training_dataset_detail(dset, model, quality, tp, fp, tn, fn)
        call quality%kill()
    end subroutine classify_training_dataset

    subroutine classify_training_dataset_detail( dset, model, quality, tp, fp, tn, fn )
        type(cavg_quality_training_dataset), intent(in)    :: dset
        type(cavg_quality_model),            intent(in)    :: model
        type(cavg_quality_result),           intent(inout) :: quality
        integer,                             intent(out)   :: tp, fp, tn, fn
        logical, allocatable :: pred(:), ref(:)
        integer :: nfit
        call quality%kill()
        quality%features    = dset%features
        quality%hard_reject = dset%hard_reject
        call model%classify(quality)
        nfit = count_trainable_classes(dset)
        if( nfit == 0 )then
            tp = 0
            fp = 0
            tn = 0
            fn = 0
            return
        endif
        allocate(pred(nfit), ref(nfit))
        pred = pack(quality%states > 0, .not. dset%hard_reject)
        ref  = pack(dset%manual_states > 0, .not. dset%hard_reject)
        call calc_confusion(pred, ref, tp, fp, tn, fn)
        deallocate(pred, ref)
    end subroutine classify_training_dataset_detail

    subroutine write_cavg_quality_learn_report( fname, dsets, base_model, suggested_weights, learned_model, best_score, &
                                                n_grid, top_specs, top_scores, n_top, best_tie_specs, n_best_ties )
        character(len=*),                    intent(in) :: fname
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: base_model, learned_model
        real,                                intent(in) :: suggested_weights(:)
        real,                                intent(in) :: best_score
        integer,                             intent(in) :: n_grid, n_top, n_best_ties
        type(cavg_quality_model_spec),       intent(in) :: top_specs(:), best_tie_specs(:)
        real,                                intent(in) :: top_scores(:)
        type(cavg_quality_learn_diagnostics) :: diag
        integer :: funit, ids, i, tp, fp, tn, fn
        real :: precision, recall, specificity, f1, balacc, accuracy
        call collect_learn_diagnostics(dsets, learned_model, diag)
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# cluster_cavgs_quality learn report'
        write(funit,'(A,A)') 'context=', trim(learned_model%context)
        write(funit,'(A,A)') 'base_model=', trim(base_model%name)
        write(funit,'(A,A)') 'learned_model=', trim(learned_model%name)
        write(funit,'(A,F10.5)') 'macro_balanced_accuracy=', best_score
        write(funit,'(A,I0)') 'n_datasets=', size(dsets)
        write(funit,'(A,I0)') 'model_search_grid_n=', n_grid
        write(funit,'(A,I0)') 'best_tie_count=', n_best_ties
        write(funit,'(A,I0)') 'top_candidates_reported=', n_top
        write(funit,'(A)') 'note=scalar_feature_space_only_pairwise_distance_matrices_removed'
        write(funit,'(A)') 'note=legacy_raw_resolution_and_geometry_failures_are_honored_as_hard_rejects_when_present'
        write(funit,'(A)') 'note=hard_rejected_rows_are_reported_but_excluded_from_model_fit_and_scoring'
        write(funit,'(A)') 'note=feature_weights_derived_from_training_data_no_base_weight_blending'
        call write_feature_policy_grid(funit)
        call write_real_list(funit, 'grid_min_score_separations=', LEARN_MINSEPS)
        call write_real_list(funit, 'grid_boundary_margins=', LEARN_MARGINS)
        if( model_is_pool_context(learned_model) ) &
            call write_real_list(funit, 'grid_pool_min_accept_fracs=', LEARN_POOL_FRACS)
        write(funit,'(A)', advance='no') 'suggested_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') suggested_weights(i)
        end do
        write(funit,*)
        write(funit,'(A)') ''
        call write_learn_search_diagnostics(funit, learned_model, diag)
        call write_feature_screen_diagnostics(funit, dsets, base_model, suggested_weights, learned_model, best_score)
        write(funit,'(A)') ''
        call write_candidate_table_header(funit, 'top_candidate_header=')
        do i = 1, n_top
            call write_candidate_row(funit, 'top_candidate', i, top_scores(i), top_specs(i))
        end do
        write(funit,'(A)') ''
        call write_candidate_table_header(funit, 'best_tie_header=')
        do i = 1, n_best_ties
            call write_candidate_row(funit, 'best_tie', i, best_score, best_tie_specs(i))
        end do
        write(funit,'(A)') ''
        write(funit,'(A)') 'dataset,n_classes,n_trainable,tp,fp,tn,fn,precision,recall,specificity,f1,'//&
            'balanced_accuracy,accuracy,hard_rejected_manual_good'
        do ids = 1, size(dsets)
            call classify_training_dataset(dsets(ids), learned_model, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
            write(funit,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,I0)') &
                trim(dsets(ids)%dataset_id), ',', dsets(ids)%ncls, ',', count_trainable_classes(dsets(ids)), ',', &
                tp, ',', fp, ',', tn, ',', fn, ',', precision, ',', recall, ',', specificity, ',', f1, ',', &
                balacc, ',', accuracy, ',', count(dsets(ids)%hard_reject .and. (dsets(ids)%manual_states > 0))
        end do
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
    end subroutine write_cavg_quality_learn_report

    subroutine write_feature_screen_diagnostics( funit, dsets, base_model, suggested_weights, learned_model, best_score )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: base_model, learned_model
        real,                                intent(in) :: suggested_weights(:)
        real,                                intent(in) :: best_score
        write(funit,'(A)') ''
        write(funit,'(A)') '# feature-screen diagnostics'
        write(funit,'(A)') &
            'note=feature_group_lodo_rows_measure_holdout_separability_with_a_best_holdout_threshold_not_a_promoted_model'
        call write_feature_signal_diagnostics(funit, dsets, base_model, suggested_weights, learned_model)
        write(funit,'(A)') ''
        call write_feature_drop_diagnostics(funit, dsets, learned_model, best_score)
        write(funit,'(A)') ''
        call write_lodo_group_screen(funit, dsets)
    end subroutine write_feature_screen_diagnostics

    subroutine write_feature_signal_diagnostics( funit, dsets, base_model, suggested_weights, learned_model )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: base_model, learned_model
        real,                                intent(in) :: suggested_weights(:)
        integer :: ifeat, inverted
        real    :: pooled_auc, mean_auc, min_auc, max_auc
        write(funit,'(A)') 'feature_signal_header=feature,pooled_auc,mean_dataset_auc,min_dataset_auc,'//&
            'max_dataset_auc,inverted_datasets,base_weight,suggested_weight,learned_weight'
        do ifeat = 1, CAVG_QUALITY_NFEATS
            call feature_auc_summary(dsets, ifeat, pooled_auc, mean_auc, min_auc, max_auc, inverted)
            write(funit,'(A,A,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,I0,A,ES14.6,A,ES14.6,A,ES14.6)') &
                'feature_signal,', trim(cavg_quality_feature_name(ifeat)), ',', pooled_auc, ',', mean_auc, ',', &
                min_auc, ',', max_auc, ',', inverted, ',', base_model%weights(ifeat), ',', suggested_weights(ifeat), &
                ',', learned_model%weights(ifeat)
        end do
    end subroutine write_feature_signal_diagnostics

    subroutine feature_auc_summary( dsets, ifeat, pooled_auc, mean_auc, min_auc, max_auc, inverted_datasets )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        integer,                             intent(in)  :: ifeat
        real,                                intent(out) :: pooled_auc, mean_auc, min_auc, max_auc
        integer,                             intent(out) :: inverted_datasets
        integer :: ids, nused
        real    :: auc
        pooled_auc        = pooled_feature_auc(dsets, ifeat, 0)
        mean_auc          = 0.0
        min_auc           = huge(1.0)
        max_auc           = -huge(1.0)
        inverted_datasets = 0
        nused             = 0
        do ids = 1, size(dsets)
            if( count_trainable_classes(dsets(ids)) == 0 ) cycle
            auc = dataset_feature_auc(dsets(ids), ifeat)
            mean_auc = mean_auc + auc
            min_auc  = min(min_auc, auc)
            max_auc  = max(max_auc, auc)
            if( auc < 0.5 - EPS ) inverted_datasets = inverted_datasets + 1
            nused = nused + 1
        end do
        if( nused > 0 ) mean_auc = mean_auc / real(nused)
        if( nused == 0 )then
            min_auc = 0.5
            max_auc = 0.5
        endif
    end subroutine feature_auc_summary

    real function dataset_feature_auc( dset, ifeat )
        type(cavg_quality_training_dataset), intent(in) :: dset
        integer,                             intent(in) :: ifeat
        real,    allocatable :: vals(:)
        integer, allocatable :: refs(:)
        integer :: nfit
        nfit = count_trainable_classes(dset)
        if( nfit == 0 )then
            dataset_feature_auc = 0.5
            return
        endif
        allocate(vals(nfit), refs(nfit))
        vals = pack(dset%features(:,ifeat), .not. dset%hard_reject)
        refs = pack(dset%manual_states, .not. dset%hard_reject)
        dataset_feature_auc = auc_for_values(vals, refs)
        deallocate(vals, refs)
    end function dataset_feature_auc

    real function pooled_feature_auc( dsets, ifeat, skip_dataset )
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer,                             intent(in) :: ifeat, skip_dataset
        real,    allocatable :: vals(:)
        integer, allocatable :: refs(:)
        integer :: ids, nall, off, nfit
        nall = count_nclasses_except(dsets, skip_dataset)
        if( nall == 0 )then
            pooled_feature_auc = 0.5
            return
        endif
        allocate(vals(nall), refs(nall))
        off = 0
        do ids = 1, size(dsets)
            if( ids == skip_dataset ) cycle
            nfit = count_trainable_classes(dsets(ids))
            if( nfit == 0 ) cycle
            vals(off+1:off+nfit) = pack(dsets(ids)%features(:,ifeat), .not. dsets(ids)%hard_reject)
            refs(off+1:off+nfit) = pack(dsets(ids)%manual_states, .not. dsets(ids)%hard_reject)
            off = off + nfit
        end do
        pooled_feature_auc = auc_for_values(vals, refs)
        deallocate(vals, refs)
    end function pooled_feature_auc

    integer function count_nclasses_except( dsets, skip_dataset )
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer,                             intent(in) :: skip_dataset
        integer :: ids
        count_nclasses_except = 0
        do ids = 1, size(dsets)
            if( ids == skip_dataset ) cycle
            count_nclasses_except = count_nclasses_except + count_trainable_classes(dsets(ids))
        end do
    end function count_nclasses_except

    subroutine write_feature_drop_diagnostics( funit, dsets, learned_model, best_score )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: best_score
        type(cavg_quality_model)      :: candidate
        type(cavg_quality_model_spec) :: spec
        integer :: ifeat
        real    :: score
        write(funit,'(A)') 'feature_drop_header=feature,learned_weight,macro_balanced_accuracy,delta_vs_learned'
        do ifeat = 1, CAVG_QUALITY_NFEATS
            spec = learned_model%get_spec()
            spec%weights(ifeat) = 0.0
            call candidate%init_spec(spec)
            score = macro_balacc_for_model(dsets, candidate)
            write(funit,'(A,A,A,ES14.6,A,F10.5,A,F10.5)') 'feature_drop,', &
                trim(cavg_quality_feature_name(ifeat)), ',', learned_model%weights(ifeat), ',', score, ',', &
                score - best_score
        end do
    end subroutine write_feature_drop_diagnostics

    subroutine write_lodo_group_screen( funit, dsets )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer :: group_inds(CAVG_QUALITY_NFEATS)
        integer :: i, ipol, ninds
        write(funit,'(A)') 'feature_group_lodo_header=group,n_features,mean_auc,min_auc,min_auc_dataset,'//&
            'mean_oracle_balacc,min_oracle_balacc,min_balacc_dataset,total_tp,total_fp,total_tn,total_fn'
        do ipol = 1, CAVG_QUALITY_LEARN_N_POLICIES
            call feature_policy_indices(ipol, group_inds, ninds)
            call write_lodo_group_row(funit, trim(feature_policy_name(ipol)), dsets, group_inds(1:ninds))
        end do
        ninds = 0
        call append_feature_except(group_inds, ninds, I_CORR_FRC)
        call write_lodo_group_row(funit, 'all_without_corr_frc', dsets, group_inds(1:ninds))
        ninds = 0
        do i = 1, CAVG_QUALITY_NFEATS
            if( i == I_LOG_POP .or. i == I_NEG_LOG_RES .or. i == I_CORR_FRC ) cycle
            ninds = ninds + 1
            group_inds(ninds) = i
        end do
        call write_lodo_group_row(funit, 'general_image_no_pop_res_corr', dsets, group_inds(1:ninds))
        group_inds(1:3) = [I_PRESENCE, I_LOG_CONTRAST, I_LOG_HIST_VARIANCE]
        call write_lodo_group_row(funit, 'new_stable3', dsets, group_inds(1:3))
        group_inds(1:6) = [I_LOCVAR_FG, I_LOCVAR_BG, I_HIST_ENTROPY, I_PRESENCE, &
                           I_LOG_CONTRAST, I_LOG_HIST_VARIANCE]
        call write_lodo_group_row(funit, 'cheap_scalar6', dsets, group_inds(1:6))
    end subroutine write_lodo_group_screen

    subroutine append_feature_except( inds, ninds, excluded_feature )
        integer, intent(inout) :: inds(:)
        integer, intent(inout) :: ninds
        integer, intent(in)    :: excluded_feature
        integer :: ifeat
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( ifeat == excluded_feature ) cycle
            ninds = ninds + 1
            inds(ninds) = ifeat
        end do
    end subroutine append_feature_except

    subroutine write_lodo_group_row( funit, group_name, dsets, feat_inds )
        integer,                             intent(in) :: funit
        character(len=*),                    intent(in) :: group_name
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer,                             intent(in) :: feat_inds(:)
        real, allocatable :: weights(:), scores(:)
        integer, allocatable :: refs(:)
        real :: auc, balacc, sum_auc, sum_balacc, min_auc, min_balacc
        integer :: holdout, tp, fp, tn, fn, total_tp, total_fp, total_tn, total_fn, nused
        integer :: low_auc_id, low_balacc_id
        character(len=LONGSTRLEN) :: low_auc_name, low_balacc_name
        if( size(dsets) == 0 .or. size(feat_inds) == 0 ) return
        allocate(weights(size(feat_inds)), source=0.0)
        sum_auc = 0.0
        sum_balacc = 0.0
        min_auc = huge(1.0)
        min_balacc = huge(1.0)
        nused = 0
        low_auc_id = 1
        low_balacc_id = 1
        total_tp = 0
        total_fp = 0
        total_tn = 0
        total_fn = 0
        do holdout = 1, size(dsets)
            if( count_trainable_classes(dsets(holdout)) == 0 ) cycle
            call lodo_group_weights(dsets, holdout, feat_inds, weights)
            call score_dataset_feature_group(dsets(holdout), feat_inds, weights, scores, refs)
            auc = auc_for_values(scores, refs)
            call best_score_threshold_balacc(scores, refs, balacc, tp, fp, tn, fn)
            sum_auc = sum_auc + auc
            sum_balacc = sum_balacc + balacc
            if( auc < min_auc )then
                min_auc = auc
                low_auc_id = holdout
            endif
            if( balacc < min_balacc )then
                min_balacc = balacc
                low_balacc_id = holdout
            endif
            total_tp = total_tp + tp
            total_fp = total_fp + fp
            total_tn = total_tn + tn
            total_fn = total_fn + fn
            nused = nused + 1
            if( allocated(scores) ) deallocate(scores)
            if( allocated(refs)   ) deallocate(refs)
        end do
        if( nused == 0 )then
            deallocate(weights)
            return
        endif
        low_auc_name = dataset_short_name(dsets(low_auc_id))
        low_balacc_name = dataset_short_name(dsets(low_balacc_id))
        write(funit,'(A,A,A,I0,A,F10.5,A,F10.5,A,A,A,F10.5,A,F10.5,A,A,A,I0,A,I0,A,I0,A,I0)') &
            'feature_group_lodo,', trim(group_name), ',', size(feat_inds), ',', sum_auc / real(nused), ',', &
            min_auc, ',', trim(low_auc_name), ',', sum_balacc / real(nused), ',', min_balacc, ',', &
            trim(low_balacc_name), ',', total_tp, ',', total_fp, ',', total_tn, ',', total_fn
        deallocate(weights)
        if( allocated(scores) ) deallocate(scores)
        if( allocated(refs)   ) deallocate(refs)
    end subroutine write_lodo_group_row

    subroutine lodo_group_weights( dsets, holdout, feat_inds, weights )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        integer,                             intent(in)  :: holdout
        integer,                             intent(in)  :: feat_inds(:)
        real,                                intent(out) :: weights(:)
        integer :: i
        real    :: auc
        if( size(weights) /= size(feat_inds) ) THROW_HARD('lodo_group_weights: invalid weight size')
        do i = 1, size(feat_inds)
            auc = pooled_feature_auc(dsets, feat_inds(i), holdout)
            weights(i) = max(0.0, auc - 0.5)
        end do
        if( sum(weights) > EPS )then
            weights = weights / sum(weights)
        else
            weights = 1.0 / real(size(weights))
        endif
    end subroutine lodo_group_weights

    subroutine score_dataset_feature_group( dset, feat_inds, weights, scores, refs )
        type(cavg_quality_training_dataset), intent(in)  :: dset
        integer,                             intent(in)  :: feat_inds(:)
        real,                                intent(in)  :: weights(:)
        real, allocatable,                   intent(out) :: scores(:)
        integer, allocatable,                intent(out) :: refs(:)
        integer :: i, j, nfit, ifit
        if( size(feat_inds) /= size(weights) ) THROW_HARD('score_dataset_feature_group: invalid weight size')
        nfit = count_trainable_classes(dset)
        allocate(scores(nfit), refs(nfit))
        scores = 0.0
        refs = 0
        ifit = 0
        do i = 1, dset%ncls
            if( dset%hard_reject(i) ) cycle
            ifit = ifit + 1
            refs(ifit) = dset%manual_states(i)
            do j = 1, size(feat_inds)
                scores(ifit) = scores(ifit) + weights(j) * dset%features(i, feat_inds(j))
            end do
        end do
    end subroutine score_dataset_feature_group

    subroutine best_score_threshold_balacc( scores, refs, best_balacc, tp, fp, tn, fn )
        real,    intent(in)  :: scores(:)
        integer, intent(in)  :: refs(:)
        real,    intent(out) :: best_balacc
        integer, intent(out) :: tp, fp, tn, fn
        logical, allocatable :: pred(:), ref(:)
        integer :: i, cand_tp, cand_fp, cand_tn, cand_fn
        real :: precision, recall, specificity, f1, balacc, accuracy
        if( size(scores) /= size(refs) ) THROW_HARD('best_score_threshold_balacc: size mismatch')
        best_balacc = -huge(1.0)
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        if( size(scores) == 0 )then
            best_balacc = 0.5
            return
        endif
        allocate(pred(size(scores)), ref(size(scores)))
        ref = refs > 0
        call test_threshold(minval(scores) - EPS)
        do i = 1, size(scores)
            call test_threshold(scores(i))
        end do
        call test_threshold(maxval(scores) + EPS)
    contains
        subroutine test_threshold( threshold )
            real, intent(in) :: threshold
            pred = scores >= threshold
            call calc_confusion(pred, ref, cand_tp, cand_fp, cand_tn, cand_fn)
            call calc_binary_metrics(cand_tp, cand_fp, cand_tn, cand_fn, precision, recall, specificity, f1, balacc, &
                accuracy)
            if( balacc > best_balacc + EPS )then
                best_balacc = balacc
                tp = cand_tp
                fp = cand_fp
                tn = cand_tn
                fn = cand_fn
            endif
        end subroutine test_threshold
    end subroutine best_score_threshold_balacc

    function dataset_short_name( dset ) result( name )
        type(cavg_quality_training_dataset), intent(in) :: dset
        character(len=LONGSTRLEN) :: name
        character(len=LONGSTRLEN) :: tmp
        integer :: i, last_slash, txt_pos
        tmp = trim(dset%fname)
        last_slash = 0
        do i = 1, len_trim(tmp)
            if( tmp(i:i) == '/' ) last_slash = i
        end do
        name = adjustl(tmp(last_slash+1:))
        if( index(name, 'cavgs_quality_analysis_') == 1 ) &
            name = adjustl(name(len('cavgs_quality_analysis_')+1:))
        txt_pos = index(name, '.txt')
        if( txt_pos > 1 ) name = name(:txt_pos-1)
    end function dataset_short_name

    subroutine collect_learn_diagnostics( dsets, model, diag )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        type(cavg_quality_model),            intent(in)  :: model
        type(cavg_quality_learn_diagnostics),intent(out) :: diag
        type(cavg_quality_result) :: quality
        integer :: ids, tp, fp, tn, fn
        logical :: lowsep, single_cluster, otsu_like, rescue_like, min_accept_like
        diag%n_datasets = size(dsets)
        do ids = 1, size(dsets)
            call classify_training_dataset_detail(dsets(ids), model, quality, tp, fp, tn, fn)
            diag%total_fp = diag%total_fp + fp
            diag%total_fn = diag%total_fn + fn
            lowsep = quality%separation < model%min_score_separation
            single_cluster = quality%nclust <= 1 .or. .not. quality%used_threshold
            otsu_like = threshold_policy_looks_otsu(quality, model)
            rescue_like = cluster_rescue_looks_active(quality, model)
            min_accept_like = min_accept_fraction_looks_active(quality, model)
            if( lowsep )then
                diag%n_lowsep  = diag%n_lowsep + 1
                diag%lowsep_fp = diag%lowsep_fp + fp
                diag%lowsep_fn = diag%lowsep_fn + fn
            endif
            if( single_cluster )then
                diag%n_single_cluster  = diag%n_single_cluster + 1
                diag%single_cluster_fp = diag%single_cluster_fp + fp
                diag%single_cluster_fn = diag%single_cluster_fn + fn
            endif
            if( otsu_like )then
                diag%n_otsu_like  = diag%n_otsu_like + 1
                diag%otsu_like_fp = diag%otsu_like_fp + fp
                diag%otsu_like_fn = diag%otsu_like_fn + fn
            endif
            if( rescue_like )then
                diag%n_rescue_like  = diag%n_rescue_like + 1
                diag%rescue_like_fp = diag%rescue_like_fp + fp
                diag%rescue_like_fn = diag%rescue_like_fn + fn
            endif
            if( min_accept_like )then
                diag%n_min_accept_like  = diag%n_min_accept_like + 1
                diag%min_accept_like_fp = diag%min_accept_like_fp + fp
                diag%min_accept_like_fn = diag%min_accept_like_fn + fn
            endif
            call quality%kill()
        end do
    end subroutine collect_learn_diagnostics

    logical function threshold_policy_looks_otsu( quality, model )
        type(cavg_quality_result), intent(in) :: quality
        type(cavg_quality_model),  intent(in) :: model
        threshold_policy_looks_otsu = .false.
        if( .not. quality%used_threshold ) return
        if( quality%nclust /= 2 ) return
        if( quality%separation < model%min_score_separation )then
            threshold_policy_looks_otsu = model%use_lowsep_otsu
        else if( model%use_otsu_window )then
            threshold_policy_looks_otsu = abs(quality%threshold_margin - model%boundary_margin) > 1.0e-4
        endif
    end function threshold_policy_looks_otsu

    logical function cluster_rescue_looks_active( quality, model )
        type(cavg_quality_result), intent(in) :: quality
        type(cavg_quality_model),  intent(in) :: model
        cluster_rescue_looks_active = .false.
        if( .not. model%use_cluster_rescue ) return
        if( .not. allocated(quality%states) ) return
        if( .not. allocated(quality%scores) ) return
        if( .not. allocated(quality%labels) ) return
        if( quality%good_label <= 0 ) return
        cluster_rescue_looks_active = any(quality%states > 0 .and. &
            quality%scores < quality%threshold - EPS .and. quality%labels == quality%good_label)
    end function cluster_rescue_looks_active

    logical function min_accept_fraction_looks_active( quality, model )
        type(cavg_quality_result), intent(in) :: quality
        type(cavg_quality_model),  intent(in) :: model
        min_accept_fraction_looks_active = .false.
        if( .not. model%enforce_min_accept_frac ) return
        if( .not. quality%used_threshold ) return
        min_accept_fraction_looks_active = quality%threshold_margin > model%boundary_margin + 1.0e-4
    end function min_accept_fraction_looks_active

    subroutine write_learn_search_diagnostics( funit, learned_model, diag )
        integer,                             intent(in) :: funit
        type(cavg_quality_model),            intent(in) :: learned_model
        type(cavg_quality_learn_diagnostics),intent(in) :: diag
        write(funit,'(A)') 'search_diagnostic_header=level,parameter,status,detail'
        call write_search_diagnostic(funit, 'note', 'feature_policy', trim(learned_model%feature_policy), &
            'selected_policy_encoded_by_zeroed_model_weights')
        call write_search_diagnostic(funit, 'note', 'feature_weights', 'data_auc_no_base_blending', &
            'learned_weights_are_derived_abinitio_from_non_hard_rejected_training_rows')
        call write_grid_position_diagnostic(funit, 'min_score_separation', learned_model%min_score_separation, &
            LEARN_MINSEPS, 'best_at_lowest_value_consider_lower_if_accept_all_remains_too_common', &
            'best_at_highest_value_consider_higher_if_unstable_splits_remain')
        call write_grid_position_diagnostic(funit, 'boundary_margin', learned_model%boundary_margin, &
            LEARN_MARGINS, 'best_at_lowest_value_consider_more_negative_if_junk_leaks_after_validation', &
            'best_at_highest_value_consider_more_positive_if_good_classes_are_rejected')
        if( model_is_pool_context(learned_model) )then
            call write_grid_position_diagnostic(funit, 'min_accept_frac', learned_model%min_accept_frac, &
                LEARN_POOL_FRACS, 'best_at_lowest_value_consider_lower_if_pool_model_keeps_too_much_junk', &
                'best_at_highest_value_consider_higher_if_pool_model_rejects_good_classes')
        else
            call write_search_diagnostic(funit, 'note', 'min_accept_frac', 'not_searched_in_chunk_context', &
                'enforce_min_accept_frac is false for the current learned model')
        endif
        call write_policy_parameter_diagnostics(funit, learned_model, diag)
    end subroutine write_learn_search_diagnostics

    subroutine write_grid_position_diagnostic( funit, param, val, grid, low_status, high_status )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: param, low_status, high_status
        real,             intent(in) :: val, grid(:)
        if( real_close(val, minval(grid)) )then
            call write_search_diagnostic(funit, 'warning', param, low_status, &
                'learned_value_is_on_lower_edge_of_current_grid')
        else if( real_close(val, maxval(grid)) )then
            call write_search_diagnostic(funit, 'warning', param, high_status, &
                'learned_value_is_on_upper_edge_of_current_grid')
        else
            call write_search_diagnostic(funit, 'note', param, 'interior', &
                'learned_value_is_not_on_a_search_grid_edge')
        endif
    end subroutine write_grid_position_diagnostic

    logical function real_close( a, b )
        real, intent(in) :: a, b
        real_close = abs(a - b) <= 1.0e-5 * max(1.0, abs(a), abs(b))
    end function real_close

    subroutine write_policy_parameter_diagnostics( funit, model, diag )
        integer,                              intent(in) :: funit
        type(cavg_quality_model),             intent(in) :: model
        type(cavg_quality_learn_diagnostics), intent(in) :: diag
        character(len=LONGSTRLEN) :: detail
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'inherited=', model%use_lowsep_otsu, ';active_datasets=', &
            diag%n_lowsep, ';fp=', diag%lowsep_fp, ';fn=', diag%lowsep_fn
        call write_search_diagnostic(funit, policy_level(diag%n_lowsep, diag%lowsep_fp, diag%lowsep_fn), &
            'use_lowsep_otsu', 'not_searched', trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'inherited=', model%use_otsu_window, ';otsu_like_datasets=', &
            diag%n_otsu_like, ';fp=', diag%otsu_like_fp, ';fn=', diag%otsu_like_fn
        call write_search_diagnostic(funit, policy_level(diag%n_otsu_like, diag%otsu_like_fp, diag%otsu_like_fn), &
            'use_otsu_window', 'not_searched', trim(detail))
        write(detail,'(A,F8.4,A,I0,A,I0,A,I0)') 'inherited=', model%otsu_min_offset, ';otsu_like_datasets=', &
            diag%n_otsu_like, ';fp=', diag%otsu_like_fp, ';fn=', diag%otsu_like_fn
        call write_search_diagnostic(funit, policy_level(diag%n_otsu_like, diag%otsu_like_fp, diag%otsu_like_fn), &
            'otsu_min_offset', 'not_searched', trim(detail))
        write(detail,'(A,F8.4,A,I0,A,I0,A,I0)') 'inherited=', model%otsu_max_offset, ';otsu_like_datasets=', &
            diag%n_otsu_like, ';fp=', diag%otsu_like_fp, ';fn=', diag%otsu_like_fn
        call write_search_diagnostic(funit, policy_level(diag%n_otsu_like, diag%otsu_like_fp, diag%otsu_like_fn), &
            'otsu_max_offset', 'not_searched', trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'inherited=', model%use_cluster_rescue, ';rescue_like_datasets=', &
            diag%n_rescue_like, ';fp=', diag%rescue_like_fp, ';fn=', diag%rescue_like_fn
        call write_search_diagnostic(funit, rescue_policy_level(model, diag), 'use_cluster_rescue', 'not_searched', &
            trim(detail))
        write(detail,'(A,F8.4,A,I0,A,I0,A,I0)') 'inherited=', model%cluster_rescue_margin, ';rescue_like_datasets=', &
            diag%n_rescue_like, ';fp=', diag%rescue_like_fp, ';fn=', diag%rescue_like_fn
        call write_search_diagnostic(funit, rescue_policy_level(model, diag), 'cluster_rescue_margin', 'not_searched', &
            trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'inherited=', model%enforce_min_accept_frac, ';min_accept_datasets=', &
            diag%n_min_accept_like, ';fp=', diag%min_accept_like_fp, ';fn=', diag%min_accept_like_fn
        call write_search_diagnostic(funit, min_accept_policy_level(model, diag), 'enforce_min_accept_frac', &
            'not_searched', trim(detail))
        write(detail,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'single_cluster_datasets=', diag%n_single_cluster, &
            ';fp=', diag%single_cluster_fp, ';fn=', diag%single_cluster_fn, ';total_fp=', diag%total_fp, &
            ';total_fn=', diag%total_fn
        call write_search_diagnostic(funit, policy_level(diag%n_single_cluster, diag%single_cluster_fp, &
            diag%single_cluster_fn), 'accept_all_fallback', 'implicit_not_searched', trim(detail))
    end subroutine write_policy_parameter_diagnostics

    function policy_level( nactive, fp, fn ) result( level )
        integer, intent(in) :: nactive, fp, fn
        character(len=8) :: level
        if( nactive > 0 .and. fp + fn > 0 )then
            level = 'warning'
        else
            level = 'note'
        endif
    end function policy_level

    function rescue_policy_level( model, diag ) result( level )
        type(cavg_quality_model),             intent(in) :: model
        type(cavg_quality_learn_diagnostics), intent(in) :: diag
        character(len=8) :: level
        if( model%use_cluster_rescue )then
            level = policy_level(diag%n_rescue_like, diag%rescue_like_fp, diag%rescue_like_fn)
        else if( diag%total_fn > 0 )then
            level = 'warning'
        else
            level = 'note'
        endif
    end function rescue_policy_level

    function min_accept_policy_level( model, diag ) result( level )
        type(cavg_quality_model),             intent(in) :: model
        type(cavg_quality_learn_diagnostics), intent(in) :: diag
        character(len=8) :: level
        if( model%enforce_min_accept_frac )then
            level = policy_level(diag%n_min_accept_like, diag%min_accept_like_fp, diag%min_accept_like_fn)
        else if( diag%total_fn > 0 .and. model_is_pool_context(model) )then
            level = 'warning'
        else
            level = 'note'
        endif
    end function min_accept_policy_level

    subroutine write_search_diagnostic( funit, level, param, status, detail )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: level, param, status, detail
        write(funit,'(8A)') 'search_diagnostic,', trim(level), ',', trim(param), ',', trim(status), ',', trim(detail)
    end subroutine write_search_diagnostic

    subroutine write_real_list( funit, key, vals )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: key
        real,             intent(in) :: vals(:)
        integer :: i
        write(funit,'(A)', advance='no') trim(key)
        do i = 1, size(vals)
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') vals(i)
        end do
        write(funit,*)
    end subroutine write_real_list

    subroutine write_feature_policy_grid( funit )
        integer, intent(in) :: funit
        integer :: ipol
        write(funit,'(A)', advance='no') 'grid_feature_policies='
        do ipol = 1, CAVG_QUALITY_LEARN_N_POLICIES
            if( ipol > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(A)', advance='no') trim(feature_policy_name(ipol))
        end do
        write(funit,*)
    end subroutine write_feature_policy_grid

    subroutine write_candidate_table_header( funit, key )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: key
        write(funit,'(A)') trim(key)//&
            'rank,balanced_accuracy,feature_policy,boundary_margin,min_score_separation,min_accept_frac,'//&
            'use_lowsep_otsu,use_otsu_window,use_cluster_rescue,enforce_min_accept_frac,feature_weights_semicolon'
    end subroutine write_candidate_table_header

    subroutine write_candidate_row( funit, tag, irank, score, spec )
        integer,                       intent(in) :: funit, irank
        character(len=*),              intent(in) :: tag
        real,                          intent(in) :: score
        type(cavg_quality_model_spec), intent(in) :: spec
        integer :: i
        write(funit,'(A,A,I0,A,F10.5,A,A,A,ES14.6,A,ES14.6,A,ES14.6,A,L1,A,L1,A,L1,A,L1,A)', advance='no') &
            trim(tag), ',', irank, ',', score, ',', trim(spec%feature_policy), ',', spec%boundary_margin, ',', &
            spec%min_score_separation, ',', spec%min_accept_frac, ',', spec%use_lowsep_otsu, ',', &
            spec%use_otsu_window, ',', spec%use_cluster_rescue, ',', spec%enforce_min_accept_frac, ','
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ';'
            write(funit,'(ES14.6)', advance='no') spec%weights(i)
        end do
        write(funit,*)
    end subroutine write_candidate_row

    subroutine kill_training_datasets( dsets )
        type(cavg_quality_training_dataset), intent(inout) :: dsets(:)
        integer :: i
        do i = 1, size(dsets)
            if( allocated(dsets(i)%features)      ) deallocate(dsets(i)%features)
            if( allocated(dsets(i)%manual_states) ) deallocate(dsets(i)%manual_states)
            if( allocated(dsets(i)%hard_reject)   ) deallocate(dsets(i)%hard_reject)
        end do
    end subroutine kill_training_datasets

end module simple_cavg_quality_learn
