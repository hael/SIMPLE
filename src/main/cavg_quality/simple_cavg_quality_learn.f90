!@descr: learn-mode training-table reader and model search for class-average quality
module simple_cavg_quality_learn
use simple_defs,               only: logfhandle, LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string,             only: string
use simple_string_utils,       only: str2int, str2real, str_is_true, csv_field
use simple_cavg_quality_feats, only: cavg_quality_feature_name
use simple_cavg_quality_model, only: CAVG_QUALITY_DEFAULT_WEIGHTS, &
    CAVG_QUALITY_DEFAULT_HIST_DMAT_WEIGHT, cavg_quality_model
use simple_cavg_quality_stats, only: calc_confusion, calc_binary_metrics, auc_for_values
use simple_cavg_quality_types, only: CAVG_REJECTION_POOL, CAVG_QUALITY_NFEATS, EPS, cavg_quality_model_spec, &
    cavg_quality_result, cavg_quality_training_dataset, cavg_quality_learn_diagnostics
implicit none
private
#include "simple_local_flags.inc"

public :: learn_cavg_quality_model

integer, parameter :: CAVG_QUALITY_LEARN_TOP_K = 10
real, parameter :: LEARN_WEIGHT_ALPHAS(5) = [0.0, 0.25, 0.50, 0.75, 1.0]
real, parameter :: LEARN_MINSEPS(5)       = [0.05, 0.10, 0.15, 0.20, 0.30]
real, parameter :: LEARN_MARGINS(11)      = [-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, &
                                              -0.05, 0.0, 0.05, 0.10, 0.20]
real, parameter :: LEARN_POOL_FRACS(5)    = [0.50, 0.60, 0.65, 0.70, 0.80]
real, parameter :: LEARN_HIST_WEIGHTS(5)  = [CAVG_QUALITY_DEFAULT_HIST_DMAT_WEIGHT, 0.0, 0.25, 0.75, 1.0]

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
        real :: score, best_score, alpha
        integer :: i, ia, ih, im, isep, ifrac, max_grid, n_grid, n_top, n_best_ties
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
        max_grid = size(LEARN_WEIGHT_ALPHAS) * size(LEARN_HIST_WEIGHTS) * size(LEARN_MINSEPS) * size(LEARN_MARGINS)
        if( learned_model%rejection_type == CAVG_REJECTION_POOL ) max_grid = max_grid * size(LEARN_POOL_FRACS)
        allocate(best_tie_specs(max_grid))
        do ia = 1, size(LEARN_WEIGHT_ALPHAS)
            alpha = LEARN_WEIGHT_ALPHAS(ia)
            candidate_spec = base_spec
            candidate_spec%weights = (1.0 - alpha) * base_spec%weights + alpha * suggested_weights
            do ih = 1, size(LEARN_HIST_WEIGHTS)
                candidate_spec%hist_dmat_weight = LEARN_HIST_WEIGHTS(ih)
                do isep = 1, size(LEARN_MINSEPS)
                    candidate_spec%min_score_separation = LEARN_MINSEPS(isep)
                    do im = 1, size(LEARN_MARGINS)
                        candidate_spec%boundary_margin = LEARN_MARGINS(im)
                        if( learned_model%rejection_type == CAVG_REJECTION_POOL )then
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

    subroutine read_quality_training_dataset( fname, dset )
        character(len=*),                    intent(in)    :: fname
        type(cavg_quality_training_dataset), intent(inout) :: dset
        character(len=XLONGSTRLEN) :: line
        character(len=LONGSTRLEN)  :: field
        integer :: funit, ios, nrows, irow, ifeat, feat_col(CAVG_QUALITY_NFEATS)
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
        allocate(dset%hist_dmat(nrows, nrows), source=0.0)
        allocate(dset%manual_states(nrows), source=0)
        allocate(dset%hard_reject(nrows), source=.false.)
        rewind(funit)
        irow = 0
        feat_col = 0
        have_feature_header = .false.
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( is_hist_dmat_line(line) )then
                call read_hist_dmat_line(line, dset)
                cycle
            endif
            if( is_analysis_header_line(line) )then
                call map_analysis_feature_columns(line, feat_col)
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

    subroutine read_hist_dmat_line( line, dset )
        character(len=*),                    intent(in)    :: line
        type(cavg_quality_training_dataset), intent(inout) :: dset
        integer :: i, j
        real    :: dist
        i    = str2int(trim(csv_field(line, 2)))
        j    = str2int(trim(csv_field(line, 3)))
        dist = str2real(trim(csv_field(line, 4)))
        if( .not. allocated(dset%hist_dmat) ) &
            THROW_HARD('read_hist_dmat_line: histogram distance matrix is not allocated')
        if( i < 1 .or. j < 1 .or. i > size(dset%hist_dmat, dim=1) .or. j > size(dset%hist_dmat, dim=2) ) &
            THROW_HARD('read_hist_dmat_line: histogram distance index out of range')
        dset%hist_dmat(i,j) = dist
        dset%hist_dmat(j,i) = dist
    end subroutine read_hist_dmat_line

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

    logical function is_hist_dmat_line( line )
        character(len=*), intent(in) :: line
        character(len=XLONGSTRLEN) :: tmp
        tmp = adjustl(line)
        is_hist_dmat_line = index(tmp, '# hist_dmat,') == 1
    end function is_hist_dmat_line

    subroutine calc_suggested_training_weights( dsets, weights )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        real,                                intent(out) :: weights(CAVG_QUALITY_NFEATS)
        real,    allocatable :: vals(:)
        integer, allocatable :: refs(:)
        integer :: nall, ids, i, j, off
        nall = 0
        do ids = 1, size(dsets)
            nall = nall + dsets(ids)%ncls
        end do
        allocate(vals(nall), source=0.0)
        allocate(refs(nall), source=0)
        weights = 0.0
        do j = 1, CAVG_QUALITY_NFEATS
            off = 0
            do ids = 1, size(dsets)
                do i = 1, dsets(ids)%ncls
                    vals(off+i) = dsets(ids)%features(i,j)
                    refs(off+i) = dsets(ids)%manual_states(i)
                end do
                off = off + dsets(ids)%ncls
            end do
            weights(j) = max(0.0, auc_for_values(vals, refs) - 0.5)
        end do
        if( sum(weights) > EPS )then
            weights = weights / sum(weights)
        else
            weights = CAVG_QUALITY_DEFAULT_WEIGHTS
        endif
        deallocate(vals, refs)
    end subroutine calc_suggested_training_weights

    real function macro_balacc_for_model( dsets, model )
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: model
        integer :: ids, tp, fp, tn, fn
        real :: precision, recall, specificity, f1, balacc, accuracy
        macro_balacc_for_model = 0.0
        do ids = 1, size(dsets)
            call classify_training_dataset(dsets(ids), model, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
            macro_balacc_for_model = macro_balacc_for_model + balacc
        end do
        macro_balacc_for_model = macro_balacc_for_model / real(size(dsets))
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
        call quality%kill()
        quality%features    = dset%features
        quality%hard_reject = dset%hard_reject
        quality%hist_dmat   = dset%hist_dmat
        call model%classify(quality)
        allocate(pred(dset%ncls), ref(dset%ncls))
        pred = quality%states > 0
        ref  = dset%manual_states > 0
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
        write(funit,'(A,F10.5)') 'learned_hist_dmat_weight=', learned_model%hist_dmat_weight
        write(funit,'(A,I0)') 'n_datasets=', size(dsets)
        write(funit,'(A,I0)') 'model_search_grid_n=', n_grid
        write(funit,'(A,I0)') 'best_tie_count=', n_best_ties
        write(funit,'(A,I0)') 'top_candidates_reported=', n_top
        write(funit,'(A)') 'note=histogram distance matrix weight was learned from analysis records'
        call write_real_list(funit, 'grid_weight_alphas=', LEARN_WEIGHT_ALPHAS)
        call write_real_list(funit, 'grid_hist_dmat_weights=', LEARN_HIST_WEIGHTS)
        call write_real_list(funit, 'grid_min_score_separations=', LEARN_MINSEPS)
        call write_real_list(funit, 'grid_boundary_margins=', LEARN_MARGINS)
        if( learned_model%rejection_type == CAVG_REJECTION_POOL ) &
            call write_real_list(funit, 'grid_pool_min_accept_fracs=', LEARN_POOL_FRACS)
        write(funit,'(A)', advance='no') 'suggested_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') suggested_weights(i)
        end do
        write(funit,*)
        write(funit,'(A)') ''
        call write_learn_search_diagnostics(funit, base_model, suggested_weights, learned_model, diag)
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
        write(funit,'(A)') 'dataset,n_classes,tp,fp,tn,fn,precision,recall,specificity,f1,'//&
            'balanced_accuracy,accuracy,hard_rejected_manual_good'
        do ids = 1, size(dsets)
            call classify_training_dataset(dsets(ids), learned_model, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
            write(funit,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,I0)') &
                trim(dsets(ids)%dataset_id), ',', dsets(ids)%ncls, ',', tp, ',', fp, ',', tn, ',', fn, ',', &
                precision, ',', recall, ',', specificity, ',', f1, ',', balacc, ',', accuracy, ',', &
                count(dsets(ids)%hard_reject .and. (dsets(ids)%manual_states > 0))
        end do
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
    end subroutine write_cavg_quality_learn_report

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

    subroutine write_learn_search_diagnostics( funit, base_model, suggested_weights, learned_model, diag )
        integer,                             intent(in) :: funit
        type(cavg_quality_model),            intent(in) :: base_model, learned_model
        real,                                intent(in) :: suggested_weights(:)
        type(cavg_quality_learn_diagnostics),intent(in) :: diag
        write(funit,'(A)') 'search_diagnostic_header=level,parameter,status,detail'
        call write_weight_alpha_diagnostic(funit, base_model%weights, suggested_weights, learned_model%weights)
        call write_grid_position_diagnostic(funit, 'hist_dmat_weight', learned_model%hist_dmat_weight, &
            LEARN_HIST_WEIGHTS, 'best_at_zero_hist_distance_disabled', 'best_at_one_hist_distance_dominant')
        call write_grid_position_diagnostic(funit, 'min_score_separation', learned_model%min_score_separation, &
            LEARN_MINSEPS, 'best_at_lowest_value_consider_lower_if_accept_all_remains_too_common', &
            'best_at_highest_value_consider_higher_if_unstable_splits_remain')
        call write_grid_position_diagnostic(funit, 'boundary_margin', learned_model%boundary_margin, &
            LEARN_MARGINS, 'best_at_lowest_value_consider_more_negative_if_junk_leaks_after_validation', &
            'best_at_highest_value_consider_more_positive_if_good_classes_are_rejected')
        if( learned_model%rejection_type == CAVG_REJECTION_POOL )then
            call write_grid_position_diagnostic(funit, 'min_accept_frac', learned_model%min_accept_frac, &
                LEARN_POOL_FRACS, 'best_at_lowest_value_consider_lower_if_pool_model_keeps_too_much_junk', &
                'best_at_highest_value_consider_higher_if_pool_model_rejects_good_classes')
        else
            call write_search_diagnostic(funit, 'note', 'min_accept_frac', 'not_searched_in_chunk_context', &
                'enforce_min_accept_frac is false for the current learned model')
        endif
        call write_policy_parameter_diagnostics(funit, learned_model, diag)
    end subroutine write_learn_search_diagnostics

    subroutine write_weight_alpha_diagnostic( funit, base_weights, suggested_weights, learned_weights )
        integer, intent(in) :: funit
        real,    intent(in) :: base_weights(:), suggested_weights(:), learned_weights(:)
        if( weights_close(base_weights, suggested_weights) )then
            call write_search_diagnostic(funit, 'note', 'feature_weight_alpha', 'suggested_equals_base', &
                'AUC_suggested_weights_match_base_weights_no_alpha_range_signal')
        else if( weights_close(learned_weights, base_weights) )then
            call write_search_diagnostic(funit, 'note', 'feature_weight_alpha', 'best_at_base_edge', &
                'best_candidate_used_alpha_zero_suggested_weights_did_not_help')
        else if( weights_close(learned_weights, suggested_weights) )then
            call write_search_diagnostic(funit, 'warning', 'feature_weight_alpha', 'best_at_suggested_edge', &
                'best_candidate_used_alpha_one_consider_alternative_weight_proposal_if_validation_fails')
        else
            call write_search_diagnostic(funit, 'note', 'feature_weight_alpha', 'interior', &
                'best_candidate_used_a_blend_of_base_and_suggested_weights')
        endif
    end subroutine write_weight_alpha_diagnostic

    logical function weights_close( a, b )
        real, intent(in) :: a(:), b(:)
        weights_close = .false.
        if( size(a) /= size(b) ) return
        weights_close = maxval(abs(a - b)) <= 1.0e-5
    end function weights_close

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
        else if( diag%total_fn > 0 .and. model%rejection_type == CAVG_REJECTION_POOL )then
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

    subroutine write_candidate_table_header( funit, key )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: key
        write(funit,'(A)') trim(key)//&
            'rank,balanced_accuracy,boundary_margin,min_score_separation,hist_dmat_weight,min_accept_frac,'//&
            'use_lowsep_otsu,use_otsu_window,use_cluster_rescue,enforce_min_accept_frac,feature_weights_semicolon'
    end subroutine write_candidate_table_header

    subroutine write_candidate_row( funit, tag, irank, score, spec )
        integer,                       intent(in) :: funit, irank
        character(len=*),              intent(in) :: tag
        real,                          intent(in) :: score
        type(cavg_quality_model_spec), intent(in) :: spec
        integer :: i
        write(funit,'(A,A,I0,A,F10.5,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,L1,A,L1,A,L1,A,L1,A)', advance='no') &
            trim(tag), ',', irank, ',', score, ',', spec%boundary_margin, ',', spec%min_score_separation, ',', &
            spec%hist_dmat_weight, ',', spec%min_accept_frac, ',', spec%use_lowsep_otsu, ',', spec%use_otsu_window, &
            ',', spec%use_cluster_rescue, ',', spec%enforce_min_accept_frac, ','
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
            if( allocated(dsets(i)%hist_dmat)     ) deallocate(dsets(i)%hist_dmat)
            if( allocated(dsets(i)%manual_states) ) deallocate(dsets(i)%manual_states)
            if( allocated(dsets(i)%hard_reject)   ) deallocate(dsets(i)%hard_reject)
        end do
    end subroutine kill_training_datasets

end module simple_cavg_quality_learn
