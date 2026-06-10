!@descr: learn-mode training-table reader and model search for class-average quality
module simple_cavg_quality_learn
use simple_defs,               only: logfhandle, LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string,             only: string
use simple_string_utils,       only: str2int, str2real, str_is_true, csv_field
use simple_cavg_quality_feats, only: cavg_quality_feature_name, cavg_quality_feature_family
use simple_cavg_quality_model, only: cavg_quality_model, cavg_quality_classify_cache, &
    build_classify_cache, kill_classify_cache, cached_decision_confusion
use simple_cavg_quality_stats, only: calc_confusion, calc_binary_metrics, auc_for_values
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, cavg_quality_model_spec, &
    cavg_quality_result, cavg_quality_training_dataset, cavg_quality_learn_diagnostics
use simple_srch_sort_loc,      only: hpsort
implicit none
private
#include "simple_local_flags.inc"

public :: learn_cavg_quality_model
public :: evaluate_cavg_quality_model
public :: evaluate_cavg_quality_result

logical, parameter :: LEARN_OTSU_FLAGS(2)           = [.false., .true.]
integer, parameter :: CAVG_QUALITY_LEARN_TOP_K      = 10
integer, parameter :: CAVG_QUALITY_LEARN_N_POLICIES = 4
integer, parameter :: LEARN_ROLE_SKIP               = 0
integer, parameter :: LEARN_ROLE_BALANCED           = 1
integer, parameter :: LEARN_ROLE_RECALL_ONLY        = 2
integer, parameter :: LEARN_ROLE_SPECIFICITY_ONLY   = 3
real,    parameter :: LEARN_RECALL_ONLY_FLOOR        = 0.95
real,    parameter :: LEARN_RECALL_ONLY_PENALTY      = 3.0
real,    parameter :: LEARN_MINSEPS(7)              = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
! Positive margins deliberately over-select relative to the learned boundary.
! Chunk and pool contexts use separate grids so pool learning can test stronger
! recall protection without making stream-chunk models permissive by default.
real,    parameter :: LEARN_CHUNK_MARGINS(19)        = [-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, &
                                                       -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, &
                                                        0.25, 0.30, 0.40, 0.50, 0.60, 0.70, &
                                                        0.80]
real,    parameter :: LEARN_POOL_MARGINS(37)         = [-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, &
                                                       -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, &
                                                        0.25, 0.30, 0.40, 0.50, 0.60, 0.70, &
                                                        0.80, 0.90, 1.00, 1.10, 1.20, 1.30, &
                                                        1.40, 1.50, 1.60, 1.70, 1.80, 1.90, &
                                                        2.00, 2.25, 2.50, 2.75, 3.00, 3.50, &
                                                        4.00]
real,    parameter :: LEARN_OTSU_MIN_OFFSETS(5)     = [0.05, 0.10, 0.15, 0.25, 0.35]
real,    parameter :: LEARN_OTSU_MAX_OFFSETS(6)     = [0.25, 0.30, 0.35, 0.40, 0.50, 0.65]
real,    parameter :: LEARN_POOL_FRACS(11)          = [0.50, 0.60, 0.65, 0.70, 0.80, 0.85, &
                                                       0.90, 0.925, 0.95, 0.975, 1.00]

contains

    subroutine learn_cavg_quality_model( analysis_files, learned_model, model_fname, report_fname )
        class(string),             intent(in)    :: analysis_files(:)
        type(cavg_quality_model),  intent(inout) :: learned_model
        character(len=*),          intent(in)    :: model_fname, report_fname
        type(cavg_quality_training_dataset), allocatable :: dsets(:)
        type(cavg_quality_classify_cache),   allocatable :: caches(:)
        type(cavg_quality_model_spec) :: base_spec, candidate_spec, best_spec
        type(cavg_quality_model_spec) :: top_specs(CAVG_QUALITY_LEARN_TOP_K)
        type(cavg_quality_model_spec), allocatable :: best_tie_specs(:)
        real :: suggested_weights(CAVG_QUALITY_NFEATS)
        real :: top_scores(CAVG_QUALITY_LEARN_TOP_K)
        real :: best_score
        integer :: ipol, im, isep, ilow, iwin, iomin, iomax, max_grid
        integer :: n_grid, n_top, n_best_ties
        logical :: is_pool_context
        call load_quality_training_datasets(analysis_files, dsets)
        base_spec = learned_model%get_spec()
        is_pool_context = model_is_pool_context(learned_model)
        call calc_suggested_training_weights(dsets, suggested_weights)
        best_spec  = base_spec
        best_score = -huge(1.0)
        top_scores = -huge(1.0)
        n_top       = 0
        n_grid      = 0
        n_best_ties = 0
        max_grid = CAVG_QUALITY_LEARN_N_POLICIES * size(LEARN_MINSEPS) * &
            n_learn_margins(is_pool_context) * n_otsu_grid_combinations()
        if( is_pool_context ) max_grid = max_grid * size(LEARN_POOL_FRACS)
        if( max_grid > 10000 ) write(logfhandle,'(A,I0)') &
            '>>> CAVG QUALITY LEARN CANDIDATE GRID SIZE: ', max_grid
        allocate(best_tie_specs(max_grid))
        allocate(caches(size(dsets)))
        do ipol = 1, CAVG_QUALITY_LEARN_N_POLICIES
            candidate_spec = base_spec
            candidate_spec%feature_policy = feature_policy_name(ipol)
            candidate_spec%weights = suggested_weights
            call apply_feature_policy(ipol, candidate_spec%weights)
            ! For a fixed weight vector, the per-dataset scores,
            ! k-medoids partition, raw threshold, and Otsu threshold are
            ! reused across the whole threshold-control sub-grid.
            call build_policy_caches(dsets, candidate_spec%weights, caches)
            do isep = 1, size(LEARN_MINSEPS)
                candidate_spec%min_score_separation = LEARN_MINSEPS(isep)
                do im = 1, n_learn_margins(is_pool_context)
                    candidate_spec%boundary_margin = learn_boundary_margin(im, is_pool_context)
                    do ilow = 1, size(LEARN_OTSU_FLAGS)
                        candidate_spec%use_lowsep_otsu = LEARN_OTSU_FLAGS(ilow)
                        do iwin = 1, size(LEARN_OTSU_FLAGS)
                            candidate_spec%use_otsu_window = LEARN_OTSU_FLAGS(iwin)
                            if( candidate_spec%use_otsu_window )then
                                do iomin = 1, size(LEARN_OTSU_MIN_OFFSETS)
                                    candidate_spec%otsu_min_offset = LEARN_OTSU_MIN_OFFSETS(iomin)
                                    do iomax = 1, size(LEARN_OTSU_MAX_OFFSETS)
                                        candidate_spec%otsu_max_offset = LEARN_OTSU_MAX_OFFSETS(iomax)
                                        if( candidate_spec%otsu_max_offset <= &
                                            candidate_spec%otsu_min_offset + EPS ) cycle
                                        call evaluate_candidate_spec(dsets, caches, candidate_spec, &
                                            is_pool_context, n_grid, best_spec, best_score, &
                                            best_tie_specs, n_best_ties, top_specs, top_scores, n_top)
                                    end do
                                end do
                            else
                                candidate_spec%otsu_min_offset = base_spec%otsu_min_offset
                                candidate_spec%otsu_max_offset = base_spec%otsu_max_offset
                                call evaluate_candidate_spec(dsets, caches, candidate_spec, &
                                    is_pool_context, &
                                    n_grid, best_spec, best_score, best_tie_specs, n_best_ties, &
                                    top_specs, top_scores, n_top)
                            endif
                        end do
                    end do
                end do
            end do
        end do
        call select_preferred_best_tie(base_spec, best_tie_specs, n_best_ties, best_spec, is_pool_context)
        best_spec%name = trim(best_spec%context)//'_learned_v1'
        call learned_model%init_spec(best_spec)
        call learned_model%write(model_fname)
        call write_cavg_quality_learn_report(report_fname, dsets, base_spec, suggested_weights, learned_model, best_score, &
            n_grid, top_specs, top_scores, n_top, best_tie_specs, n_best_ties)
        call kill_policy_caches(caches)
        deallocate(caches)
        call kill_training_datasets(dsets)
        deallocate(best_tie_specs)
        deallocate(dsets)
    end subroutine learn_cavg_quality_model

    subroutine evaluate_cavg_quality_model( analysis_files, model, report_fname )
        class(string),            intent(in) :: analysis_files(:)
        type(cavg_quality_model), intent(in) :: model
        character(len=*),         intent(in) :: report_fname
        type(cavg_quality_training_dataset), allocatable :: dsets(:)
        real :: eval_score
        call load_quality_training_datasets(analysis_files, dsets)
        eval_score = macro_balacc_for_model(dsets, model)
        call write_cavg_quality_evaluate_report(report_fname, dsets, model, eval_score)
        call kill_training_datasets(dsets)
        deallocate(dsets)
    end subroutine evaluate_cavg_quality_model

    subroutine evaluate_cavg_quality_result( quality, reference_states, model, dataset_id, report_fname )
        type(cavg_quality_result), intent(in) :: quality
        integer,                   intent(in) :: reference_states(:)
        type(cavg_quality_model),  intent(in) :: model
        character(len=*),          intent(in) :: dataset_id, report_fname
        type(cavg_quality_training_dataset), allocatable :: dsets(:)
        integer :: ncls
        real :: eval_score
        ncls = size(reference_states)
        if( .not. allocated(quality%features) ) THROW_HARD('evaluate_cavg_quality_result: missing quality features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('evaluate_cavg_quality_result: missing hard-reject mask')
        if( size(quality%features, dim=1) /= ncls ) THROW_HARD('evaluate_cavg_quality_result: feature size mismatch')
        if( size(quality%hard_reject) /= ncls ) THROW_HARD('evaluate_cavg_quality_result: hard-reject size mismatch')
        allocate(dsets(1))
        dsets(1)%fname      = trim(dataset_id)
        dsets(1)%dataset_id = trim(dataset_id)
        dsets(1)%ncls       = ncls
        dsets(1)%features      = quality%features
        dsets(1)%manual_states = reference_states
        dsets(1)%hard_reject   = quality%hard_reject
        eval_score = macro_balacc_for_model(dsets, model)
        call write_cavg_quality_evaluate_report(report_fname, dsets, model, eval_score)
        call kill_training_datasets(dsets)
        deallocate(dsets)
    end subroutine evaluate_cavg_quality_result

    subroutine load_quality_training_datasets( analysis_files, dsets )
        class(string), intent(in) :: analysis_files(:)
        type(cavg_quality_training_dataset), allocatable, intent(inout) :: dsets(:)
        integer :: i
        if( size(analysis_files) == 0 ) THROW_HARD('load_quality_training_datasets: empty analysis file table')
        if( allocated(dsets) )then
            call kill_training_datasets(dsets)
            deallocate(dsets)
        endif
        allocate(dsets(size(analysis_files)))
        do i = 1, size(analysis_files)
            call read_quality_training_dataset(analysis_files(i)%to_char(), dsets(i))
        end do
    end subroutine load_quality_training_datasets

    logical function model_is_pool_context( model )
        type(cavg_quality_model), intent(in) :: model
        model_is_pool_context = trim(model%context) == 'pool'
    end function model_is_pool_context

    integer function n_learn_margins( is_pool_context )
        logical, intent(in) :: is_pool_context
        if( is_pool_context )then
            n_learn_margins = size(LEARN_POOL_MARGINS)
        else
            n_learn_margins = size(LEARN_CHUNK_MARGINS)
        endif
    end function n_learn_margins

    real function learn_boundary_margin( imargin, is_pool_context )
        integer, intent(in) :: imargin
        logical, intent(in) :: is_pool_context
        if( is_pool_context )then
            if( imargin < 1 .or. imargin > size(LEARN_POOL_MARGINS) ) &
                THROW_HARD('learn_boundary_margin: invalid pool margin index')
            learn_boundary_margin = LEARN_POOL_MARGINS(imargin)
        else
            if( imargin < 1 .or. imargin > size(LEARN_CHUNK_MARGINS) ) &
                THROW_HARD('learn_boundary_margin: invalid chunk margin index')
            learn_boundary_margin = LEARN_CHUNK_MARGINS(imargin)
        endif
    end function learn_boundary_margin

    integer function n_otsu_grid_combinations()
        integer :: ilow, iwin, iomin, iomax
        n_otsu_grid_combinations = 0
        do ilow = 1, size(LEARN_OTSU_FLAGS)
            do iwin = 1, size(LEARN_OTSU_FLAGS)
                if( LEARN_OTSU_FLAGS(iwin) )then
                    do iomin = 1, size(LEARN_OTSU_MIN_OFFSETS)
                        do iomax = 1, size(LEARN_OTSU_MAX_OFFSETS)
                            if( LEARN_OTSU_MAX_OFFSETS(iomax) <= LEARN_OTSU_MIN_OFFSETS(iomin) + EPS ) cycle
                            n_otsu_grid_combinations = n_otsu_grid_combinations + 1
                        end do
                    end do
                else
                    n_otsu_grid_combinations = n_otsu_grid_combinations + 1
                endif
            end do
        end do
    end function n_otsu_grid_combinations

    subroutine evaluate_candidate_spec( dsets, caches, candidate_spec, is_pool_context, n_grid, best_spec, best_score, &
                                        best_tie_specs, n_best_ties, top_specs, top_scores, n_top )
        type(cavg_quality_training_dataset), intent(in)    :: dsets(:)
        type(cavg_quality_classify_cache),   intent(in)    :: caches(:)
        type(cavg_quality_model_spec),       intent(in)    :: candidate_spec
        logical,                             intent(in)    :: is_pool_context
        integer,                             intent(inout) :: n_grid, n_best_ties, n_top
        type(cavg_quality_model_spec),       intent(inout) :: best_spec, best_tie_specs(:), top_specs(:)
        real,                                intent(inout) :: best_score, top_scores(:)
        type(cavg_quality_model)      :: candidate
        type(cavg_quality_model_spec) :: scan_spec
        integer :: ifrac
        real    :: score
        if( is_pool_context )then
            do ifrac = 1, size(LEARN_POOL_FRACS)
                scan_spec = candidate_spec
                scan_spec%min_accept_frac = LEARN_POOL_FRACS(ifrac)
                call candidate%init_spec(scan_spec)
                score = macro_balacc_for_model(dsets, candidate, caches)
                n_grid = n_grid + 1
                call consider_model_candidate(candidate%get_spec(), score, best_spec, best_score, &
                    best_tie_specs, n_best_ties, top_specs, top_scores, n_top)
            end do
        else
            call candidate%init_spec(candidate_spec)
            score = macro_balacc_for_model(dsets, candidate, caches)
            n_grid = n_grid + 1
            call consider_model_candidate(candidate%get_spec(), score, best_spec, best_score, &
                best_tie_specs, n_best_ties, top_specs, top_scores, n_top)
        endif
    end subroutine evaluate_candidate_spec

    subroutine build_policy_caches( dsets, weights, caches )
        type(cavg_quality_training_dataset), intent(in)    :: dsets(:)
        real,                                intent(in)    :: weights(:)
        type(cavg_quality_classify_cache),   intent(inout) :: caches(:)
        integer :: ids
        if( size(caches) /= size(dsets) ) THROW_HARD('build_policy_caches: cache/dataset size mismatch')
        do ids = 1, size(dsets)
            call build_classify_cache(dsets(ids)%features, dsets(ids)%hard_reject, weights, caches(ids))
        end do
    end subroutine build_policy_caches

    subroutine kill_policy_caches( caches )
        type(cavg_quality_classify_cache), intent(inout) :: caches(:)
        integer :: ids
        do ids = 1, size(caches)
            call kill_classify_cache(caches(ids))
        end do
    end subroutine kill_policy_caches

    function feature_policy_name( ipolicy ) result( name )
        integer, intent(in) :: ipolicy
        character(len=32) :: name
        select case(ipolicy)
            case(1)
                name = 'microchunk'
            case(2)
                name = 'microchunk_plus_score'
            case(3)
                name = 'microchunk_plus_signal'
            case(4)
                name = 'microchunk_plus_score_signal'
            case default
                THROW_HARD('feature_policy_name: invalid feature policy')
        end select
    end function feature_policy_name

    subroutine feature_policy_mask( ipolicy, mask )
        integer, intent(in)  :: ipolicy
        logical, intent(out) :: mask(CAVG_QUALITY_NFEATS)
        mask = .false.
        call append_feature_family(mask, 'microchunk')
        select case(ipolicy)
            case(1)
                continue
            case(2)
                call append_feature_family(mask, 'score')
            case(3)
                call append_feature_family(mask, 'signal')
            case(4)
                call append_feature_family(mask, 'score')
                call append_feature_family(mask, 'signal')
            case default
                THROW_HARD('feature_policy_mask: invalid feature policy')
        end select
        if( count(mask) < 1 ) THROW_HARD('feature_policy_mask: empty feature policy')
    end subroutine feature_policy_mask

    subroutine append_feature_family( mask, family )
        logical,          intent(inout) :: mask(CAVG_QUALITY_NFEATS)
        character(len=*), intent(in)    :: family
        integer :: ifeat
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( trim(cavg_quality_feature_family(ifeat)) == trim(family) ) mask(ifeat) = .true.
        end do
    end subroutine append_feature_family

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
        integer :: dataset_col, hard_reject_col, manual_state_col
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
        dataset_col     = 0
        hard_reject_col = 0
        manual_state_col = 0
        have_feature_header = .false.
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( is_analysis_header_line(line) )then
                call map_analysis_columns(line, dataset_col, hard_reject_col, manual_state_col, feat_col)
                call require_analysis_columns(dataset_col, hard_reject_col, manual_state_col, feat_col, trim(fname))
                have_feature_header = .true.
                cycle
            endif
            if( .not. is_analysis_data_line(line) ) cycle
            if( .not. have_feature_header ) &
                THROW_HARD('read_quality_training_dataset: missing analysis header in '//trim(fname))
            irow = irow + 1
            field = csv_field(line, dataset_col)
            if( irow == 1 ) dset%dataset_id = trim(field)
            field = csv_field(line, hard_reject_col)
            dset%hard_reject(irow) = str_is_true(field)
            field = csv_field(line, manual_state_col)
            dset%manual_states(irow) = str2int(trim(field))
            do ifeat = 1, CAVG_QUALITY_NFEATS
                field = csv_field(line, feat_col(ifeat))
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

    subroutine map_analysis_columns( line, dataset_col, hard_reject_col, manual_state_col, feat_col )
        character(len=*), intent(in)  :: line
        integer,          intent(out) :: dataset_col, hard_reject_col, manual_state_col
        integer,          intent(out) :: feat_col(CAVG_QUALITY_NFEATS)
        character(len=LONGSTRLEN) :: field, name
        integer :: icol, ifeat
        feat_col = 0
        dataset_col = 0
        hard_reject_col = 0
        manual_state_col = 0
        do icol = 1, 512
            field = csv_field(line, icol)
            if( len_trim(field) == 0 ) exit
            select case(trim(field))
                case('dataset_id')
                    dataset_col = icol
                case('hard_reject')
                    hard_reject_col = icol
                case('manual_state')
                    manual_state_col = icol
            end select
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
    end subroutine map_analysis_columns

    subroutine require_analysis_columns( dataset_col, hard_reject_col, manual_state_col, feat_col, fname )
        integer,          intent(in) :: dataset_col, hard_reject_col, manual_state_col
        integer,          intent(in) :: feat_col(CAVG_QUALITY_NFEATS)
        character(len=*), intent(in) :: fname
        character(len=LONGSTRLEN) :: errmsg
        integer :: ifeat
        if( dataset_col == 0 ) THROW_HARD('read_quality_training_dataset: missing dataset_id column in '//trim(fname))
        if( hard_reject_col == 0 ) &
            THROW_HARD('read_quality_training_dataset: missing hard_reject column in '//trim(fname))
        if( manual_state_col == 0 ) &
            THROW_HARD('read_quality_training_dataset: missing manual_state column in '//trim(fname))
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( feat_col(ifeat) == 0 )then
                errmsg = 'read_quality_training_dataset: missing z_'//trim(cavg_quality_feature_name(ifeat))//&
                    ' column in '//trim(fname)
                THROW_HARD(trim(errmsg))
            endif
        end do
    end subroutine require_analysis_columns

    subroutine calc_suggested_training_weights( dsets, weights )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        real,                                intent(out) :: weights(CAVG_QUALITY_NFEATS)
        real,    allocatable :: vals(:)
        integer, allocatable :: refs(:)
        integer :: nall, ids, j, off, nfit
        nall = 0
        do ids = 1, size(dsets)
            if( dataset_learn_role(dsets(ids)) /= LEARN_ROLE_BALANCED ) cycle
            nall = nall + count_trainable_classes(dsets(ids))
        end do
        if( nall == 0 ) &
            THROW_HARD('calc_suggested_training_weights: no contrast datasets after hard rejects')
        allocate(vals(nall), source=0.0)
        allocate(refs(nall), source=0)
        weights = 0.0
        do j = 1, CAVG_QUALITY_NFEATS
            off = 0
            do ids = 1, size(dsets)
                if( dataset_learn_role(dsets(ids)) /= LEARN_ROLE_BALANCED ) cycle
                nfit = count_trainable_classes(dsets(ids))
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

    integer function count_trainable_manual_good( dset )
        type(cavg_quality_training_dataset), intent(in) :: dset
        count_trainable_manual_good = count((.not. dset%hard_reject) .and. (dset%manual_states > 0))
    end function count_trainable_manual_good

    integer function count_trainable_manual_bad( dset )
        type(cavg_quality_training_dataset), intent(in) :: dset
        count_trainable_manual_bad = count((.not. dset%hard_reject) .and. (dset%manual_states <= 0))
    end function count_trainable_manual_bad

    integer function count_hard_rejected_manual_good( dset )
        type(cavg_quality_training_dataset), intent(in) :: dset
        count_hard_rejected_manual_good = count(dset%hard_reject .and. (dset%manual_states > 0))
    end function count_hard_rejected_manual_good

    integer function dataset_learn_role( dset )
        type(cavg_quality_training_dataset), intent(in) :: dset
        integer :: ngood, nbad
        ngood = count_trainable_manual_good(dset)
        nbad  = count_trainable_manual_bad(dset)
        if( ngood > 0 .and. nbad > 0 )then
            dataset_learn_role = LEARN_ROLE_BALANCED
        else if( ngood > 0 )then
            dataset_learn_role = LEARN_ROLE_RECALL_ONLY
        else if( nbad > 0 .and. count_hard_rejected_manual_good(dset) == 0 )then
            dataset_learn_role = LEARN_ROLE_SPECIFICITY_ONLY
        else
            dataset_learn_role = LEARN_ROLE_SKIP
        endif
    end function dataset_learn_role

    function dataset_learn_role_name( role ) result( name )
        integer, intent(in) :: role
        character(len=32) :: name
        select case(role)
            case(LEARN_ROLE_BALANCED)
                name = 'balanced'
            case(LEARN_ROLE_RECALL_ONLY)
                name = 'trainable_good_only'
            case(LEARN_ROLE_SPECIFICITY_ONLY)
                name = 'trainable_bad_only'
            case default
                name = 'skip_hard_gate_blocked'
        end select
    end function dataset_learn_role_name

    real function macro_balacc_for_model( dsets, model, caches )
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: model
        type(cavg_quality_classify_cache), optional, intent(in) :: caches(:)
        integer :: ids, tp, fp, tn, fn, nused, role
        real :: balacc
        macro_balacc_for_model = 0.0
        if( present(caches) )then
            if( size(caches) /= size(dsets) ) THROW_HARD('macro_balacc_for_model: cache/dataset size mismatch')
        endif
        nused = 0
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            if( role == LEARN_ROLE_SKIP ) cycle
            if( present(caches) )then
                call classify_training_dataset_cached(dsets(ids), caches(ids), model, tp, fp, tn, fn)
            else
                call classify_training_dataset(dsets(ids), model, tp, fp, tn, fn)
            endif
            balacc = learn_balacc_from_confusion(tp, fp, tn, fn, role)
            macro_balacc_for_model = macro_balacc_for_model + balacc
            nused = nused + 1
        end do
        if( nused == 0 ) THROW_HARD('macro_balacc_for_model: no scoreable training datasets')
        macro_balacc_for_model = macro_balacc_for_model / real(nused)
    end function macro_balacc_for_model

    real function learn_balacc_from_confusion( tp, fp, tn, fn, role )
        integer, intent(in) :: tp, fp, tn, fn, role
        real :: precision, recall, specificity, f1, balacc, accuracy
        call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
        select case(role)
            case(LEARN_ROLE_BALANCED)
                learn_balacc_from_confusion = balacc
            case(LEARN_ROLE_RECALL_ONLY)
                learn_balacc_from_confusion = guarded_recall_score(recall)
            case(LEARN_ROLE_SPECIFICITY_ONLY)
                learn_balacc_from_confusion = specificity
            case default
                learn_balacc_from_confusion = 0.5
        end select
    end function learn_balacc_from_confusion

    real function guarded_recall_score( recall )
        real, intent(in) :: recall
        guarded_recall_score = recall - LEARN_RECALL_ONLY_PENALTY * &
            max(0.0, LEARN_RECALL_ONLY_FLOOR - recall)
    end function guarded_recall_score

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

    subroutine select_preferred_best_tie( base_spec, best_tie_specs, n_best_ties, best_spec, is_pool_context )
        type(cavg_quality_model_spec), intent(in)    :: base_spec, best_tie_specs(:)
        integer,                       intent(in)    :: n_best_ties
        type(cavg_quality_model_spec), intent(inout) :: best_spec
        logical,                       intent(in)    :: is_pool_context
        integer :: i, nstored
        real    :: dist, best_dist
        nstored = min(n_best_ties, size(best_tie_specs))
        if( nstored <= 1 ) return
        best_dist = huge(1.0)
        do i = 1, nstored
            dist = threshold_tie_distance(best_tie_specs(i), base_spec, is_pool_context)
            if( dist < best_dist - EPS )then
                best_dist = dist
                best_spec = best_tie_specs(i)
            endif
        end do
    end subroutine select_preferred_best_tie

    real function threshold_tie_distance( spec, base_spec, is_pool_context )
        type(cavg_quality_model_spec), intent(in) :: spec, base_spec
        logical,                       intent(in) :: is_pool_context
        real :: margin_min, margin_max
        if( is_pool_context )then
            margin_min = minval(LEARN_POOL_MARGINS)
            margin_max = maxval(LEARN_POOL_MARGINS)
        else
            margin_min = minval(LEARN_CHUNK_MARGINS)
            margin_max = maxval(LEARN_CHUNK_MARGINS)
        endif
        threshold_tie_distance = scaled_absdiff(spec%min_score_separation, base_spec%min_score_separation, &
            minval(LEARN_MINSEPS), maxval(LEARN_MINSEPS))
        threshold_tie_distance = threshold_tie_distance + scaled_absdiff(spec%boundary_margin, &
            base_spec%boundary_margin, margin_min, margin_max)
        if( is_pool_context ) threshold_tie_distance = threshold_tie_distance + scaled_absdiff( &
            spec%min_accept_frac, base_spec%min_accept_frac, minval(LEARN_POOL_FRACS), maxval(LEARN_POOL_FRACS))
        ! Boolean knobs count as a quarter of one normalized scalar axis:
        ! enough to prefer inherited behavior among exact-score ties without
        ! overwhelming a nearby threshold or margin value.
        if( spec%use_lowsep_otsu .neqv. base_spec%use_lowsep_otsu ) threshold_tie_distance = threshold_tie_distance + 0.25
        if( spec%use_otsu_window .neqv. base_spec%use_otsu_window ) threshold_tie_distance = threshold_tie_distance + 0.25
        if( spec%use_otsu_window )then
            threshold_tie_distance = threshold_tie_distance + scaled_absdiff(spec%otsu_min_offset, &
                base_spec%otsu_min_offset, minval(LEARN_OTSU_MIN_OFFSETS), maxval(LEARN_OTSU_MIN_OFFSETS))
            threshold_tie_distance = threshold_tie_distance + scaled_absdiff(spec%otsu_max_offset, &
                base_spec%otsu_max_offset, minval(LEARN_OTSU_MAX_OFFSETS), maxval(LEARN_OTSU_MAX_OFFSETS))
        endif
    end function threshold_tie_distance

    real function scaled_absdiff( val, ref, grid_min, grid_max )
        real, intent(in) :: val, ref, grid_min, grid_max
        scaled_absdiff = abs(val - ref) / max(EPS, grid_max - grid_min)
    end function scaled_absdiff

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

    ! Fast path: skip the heavy k-medoids / Otsu work and re-use a cache
    ! that was built once per (feature policy, dataset) by the grid driver.
    subroutine classify_training_dataset_cached( dset, cache, model, tp, fp, tn, fn )
        type(cavg_quality_training_dataset), intent(in) :: dset
        type(cavg_quality_classify_cache),   intent(in) :: cache
        type(cavg_quality_model),            intent(in) :: model
        integer,                             intent(out):: tp, fp, tn, fn
        call cached_decision_confusion(cache, model, dset%manual_states, tp, fp, tn, fn)
    end subroutine classify_training_dataset_cached

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

    subroutine write_cavg_quality_learn_report( fname, dsets, base_spec, suggested_weights, learned_model, best_score, &
                                                n_grid, top_specs, top_scores, n_top, best_tie_specs, n_best_ties )
        character(len=*),                    intent(in) :: fname
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model_spec),       intent(in) :: base_spec
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: suggested_weights(:)
        real,                                intent(in) :: best_score
        integer,                             intent(in) :: n_grid, n_top, n_best_ties
        type(cavg_quality_model_spec),       intent(in) :: top_specs(:), best_tie_specs(:)
        real,                                intent(in) :: top_scores(:)
        type(cavg_quality_learn_diagnostics) :: diag
        integer :: funit, i
        call collect_learn_diagnostics(dsets, learned_model, diag)
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# model_cavgs_rejection learn report'
        write(funit,'(A,A)') 'context=', trim(learned_model%context)
        write(funit,'(A,A)') 'base_model=', trim(base_spec%name)
        write(funit,'(A,A)') 'learned_model=', trim(learned_model%name)
        write(funit,'(A,F10.5)') 'macro_learn_score=', best_score
        write(funit,'(A,I0)') 'n_datasets=', size(dsets)
        write(funit,'(A,I0)') 'model_search_grid_n=', n_grid
        write(funit,'(A,I0)') 'best_tie_count=', n_best_ties
        write(funit,'(A,I0)') 'top_candidates_reported=', n_top
        write(funit,'(A)') 'note=scalar_feature_space_only'
        write(funit,'(A)') 'note=feature_policy_scans_cumulative_family_sets_starting_with_microchunk'
        write(funit,'(A)') 'note=hard_rejected_rows_are_reported_but_excluded_from_model_fit_and_scoring'
        write(funit,'(A)') 'note=feature_weights_use_only_datasets_with_both_manual_states_after_hard_rejects'
        write(funit,'(A)') 'note=trainable_good_only_datasets_are_scored_by_guarded_recall'
        write(funit,'(A)') 'note=trainable_bad_only_datasets_are_scored_by_specificity_unless_good_classes_were_hard_rejected'
        write(funit,'(A)') 'note=feature_weights_derived_from_training_data_no_base_weight_blending'
        call write_feature_policy_grid(funit)
        write(funit,'(A,ES14.6)') 'grid_recall_only_floor=', LEARN_RECALL_ONLY_FLOOR
        write(funit,'(A,ES14.6)') 'grid_recall_only_shortfall_penalty=', LEARN_RECALL_ONLY_PENALTY
        call write_real_list(funit, 'grid_min_score_separations=', LEARN_MINSEPS)
        call write_learn_margin_grid(funit, model_is_pool_context(learned_model))
        call write_logical_list(funit, 'grid_use_lowsep_otsu=', LEARN_OTSU_FLAGS)
        call write_logical_list(funit, 'grid_use_otsu_window=', LEARN_OTSU_FLAGS)
        call write_real_list(funit, 'grid_otsu_min_offsets=', LEARN_OTSU_MIN_OFFSETS)
        call write_real_list(funit, 'grid_otsu_max_offsets=', LEARN_OTSU_MAX_OFFSETS)
        if( model_is_pool_context(learned_model) ) &
            call write_real_list(funit, 'grid_pool_min_accept_fracs=', LEARN_POOL_FRACS)
        write(funit,'(A)', advance='no') 'suggested_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') suggested_weights(i)
        end do
        write(funit,*)
        write(funit,'(A)') ''
        call write_learn_search_diagnostics(funit, learned_model, diag, best_tie_specs, n_best_ties)
        call write_otsu_ablation_diagnostics(funit, dsets, learned_model)
        call write_feature_screen_diagnostics(funit, dsets, base_spec, suggested_weights, learned_model, best_score)
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
        call write_dataset_metric_table(funit, dsets, learned_model, 'learn_score')
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
    end subroutine write_cavg_quality_learn_report

    subroutine write_cavg_quality_evaluate_report( fname, dsets, model, eval_score )
        character(len=*),                    intent(in) :: fname
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: model
        real,                                intent(in) :: eval_score
        type(cavg_quality_learn_diagnostics) :: diag
        integer :: funit
        call collect_learn_diagnostics(dsets, model, diag)
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# model_cavgs_rejection evaluate report'
        write(funit,'(A,A)') 'context=', trim(model%context)
        write(funit,'(A,A)') 'model=', trim(model%name)
        write(funit,'(A,F10.5)') 'macro_evaluate_score=', eval_score
        write(funit,'(A,I0)') 'n_datasets=', size(dsets)
        write(funit,'(A)') 'note=fixed_model_no_refit'
        write(funit,'(A)') 'note=analysis_table_rows_reclassified_with_selected_model'
        write(funit,'(A)') 'note=hard_rejected_rows_are_reported_but_excluded_from_model_scoring'
        write(funit,'(A)') 'note=trainable_good_only_datasets_are_scored_by_guarded_recall'
        write(funit,'(A)') 'note=trainable_bad_only_datasets_are_scored_by_specificity_unless_good_classes_were_hard_rejected'
        call write_fixed_model_summary(funit, model)
        write(funit,'(A)') ''
        call write_evaluate_diagnostics(funit, model, diag)
        call write_otsu_ablation_diagnostics(funit, dsets, model)
        call write_dataset_metric_table(funit, dsets, model, 'evaluate_score')
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
    end subroutine write_cavg_quality_evaluate_report

    subroutine write_fixed_model_summary( funit, model )
        integer,                  intent(in) :: funit
        type(cavg_quality_model), intent(in) :: model
        integer :: i
        write(funit,'(A,A)') 'model_feature_policy=', trim(model%feature_policy)
        write(funit,'(A)', advance='no') 'model_feature_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') model%weights(i)
        end do
        write(funit,*)
        write(funit,'(A,ES14.6)') 'model_boundary_margin=', model%boundary_margin
        write(funit,'(A,ES14.6)') 'model_min_score_separation=', model%min_score_separation
        write(funit,'(A,ES14.6)') 'model_otsu_min_offset=', model%otsu_min_offset
        write(funit,'(A,ES14.6)') 'model_otsu_max_offset=', model%otsu_max_offset
        write(funit,'(A,ES14.6)') 'model_cluster_rescue_margin=', model%cluster_rescue_margin
        write(funit,'(A,ES14.6)') 'model_min_accept_frac=', model%min_accept_frac
        write(funit,'(A,L1)') 'model_use_lowsep_otsu=', model%use_lowsep_otsu
        write(funit,'(A,L1)') 'model_use_otsu_window=', model%use_otsu_window
        write(funit,'(A,L1)') 'model_use_cluster_rescue=', model%use_cluster_rescue
        write(funit,'(A,L1)') 'model_enforce_min_accept_frac=', model%enforce_min_accept_frac
    end subroutine write_fixed_model_summary

    subroutine write_dataset_metric_table( funit, dsets, model, score_name )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: model
        character(len=*),                    intent(in) :: score_name
        integer :: ids, tp, fp, tn, fn, role
        real :: precision, recall, specificity, f1, role_score, accuracy
        write(funit,'(A)') ''
        write(funit,'(A,A,A)') 'dataset,n_classes,n_trainable,trainable_manual_good,trainable_manual_bad,learn_role,', &
            'tp,fp,tn,fn,precision,recall,specificity,f1,', trim(score_name)//',accuracy,hard_rejected_manual_good'
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            call classify_training_dataset(dsets(ids), model, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, role_score, accuracy)
            role_score = learn_balacc_from_confusion(tp, fp, tn, fn, role)
            write(funit,'(A,A,I0,A,I0,A,I0,A,I0,A,A,A,I0,A,I0,A,I0,A,I0,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,I0)') &
                trim(dsets(ids)%dataset_id), ',', dsets(ids)%ncls, ',', count_trainable_classes(dsets(ids)), ',', &
                count_trainable_manual_good(dsets(ids)), ',', count_trainable_manual_bad(dsets(ids)), ',', &
                trim(dataset_learn_role_name(role)), ',', &
                tp, ',', fp, ',', tn, ',', fn, ',', precision, ',', recall, ',', specificity, ',', f1, ',', &
                role_score, ',', accuracy, ',', count_hard_rejected_manual_good(dsets(ids))
        end do
    end subroutine write_dataset_metric_table

    subroutine write_otsu_ablation_diagnostics( funit, dsets, learned_model )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: learned_model
        type(cavg_quality_model)      :: no_otsu_model
        type(cavg_quality_model_spec) :: no_otsu_spec
        integer :: ids, tp1, fp1, tn1, fn1, tp0, fp0, tn0, fn0, role
        real :: precision, recall, specificity, f1, balacc1, balacc0, accuracy
        no_otsu_spec = learned_model%get_spec()
        no_otsu_spec%use_lowsep_otsu = .false.
        no_otsu_spec%use_otsu_window = .false.
        call no_otsu_model%init_spec(no_otsu_spec)
        write(funit,'(A)') ''
        write(funit,'(A)') '# otsu-ablation diagnostics'
        write(funit,'(A)') 'otsu_ablation_header=dataset,learn_role,with_otsu_score,no_otsu_score,delta_vs_no_otsu,'//&
            'with_otsu_fp,with_otsu_fn,no_otsu_fp,no_otsu_fn'
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            call classify_training_dataset(dsets(ids), learned_model, tp1, fp1, tn1, fn1)
            call calc_binary_metrics(tp1, fp1, tn1, fn1, precision, recall, specificity, f1, balacc1, accuracy)
            balacc1 = learn_balacc_from_confusion(tp1, fp1, tn1, fn1, role)
            call classify_training_dataset(dsets(ids), no_otsu_model, tp0, fp0, tn0, fn0)
            call calc_binary_metrics(tp0, fp0, tn0, fn0, precision, recall, specificity, f1, balacc0, accuracy)
            balacc0 = learn_balacc_from_confusion(tp0, fp0, tn0, fn0, role)
            write(funit,'(A,A,A,A,A,F10.5,A,F10.5,A,F10.5,A,I0,A,I0,A,I0,A,I0)') 'otsu_ablation,', &
                trim(dataset_short_name(dsets(ids))), ',', trim(dataset_learn_role_name(role)), ',', &
                balacc1, ',', balacc0, ',', balacc1 - balacc0, ',', fp1, ',', fn1, ',', fp0, ',', fn0
        end do
    end subroutine write_otsu_ablation_diagnostics

    subroutine write_feature_screen_diagnostics( funit, dsets, base_spec, suggested_weights, learned_model, best_score )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model_spec),       intent(in) :: base_spec
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: suggested_weights(:)
        real,                                intent(in) :: best_score
        write(funit,'(A)') ''
        write(funit,'(A)') '# feature-screen diagnostics'
        call write_feature_signal_diagnostics(funit, dsets, base_spec, suggested_weights, learned_model)
        write(funit,'(A)') ''
        call write_feature_drop_diagnostics(funit, dsets, learned_model, best_score)
        write(funit,'(A)') ''
        call write_feature_policy_screen(funit, dsets)
    end subroutine write_feature_screen_diagnostics

    subroutine write_feature_signal_diagnostics( funit, dsets, base_spec, suggested_weights, learned_model )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model_spec),       intent(in) :: base_spec
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: suggested_weights(:)
        integer :: ifeat, inverted
        real    :: pooled_auc, mean_auc, min_auc, max_auc
        write(funit,'(A)') 'feature_signal_header=feature,pooled_auc,mean_dataset_auc,min_dataset_auc,'//&
            'max_dataset_auc,inverted_datasets,base_weight,suggested_weight,learned_weight'
        do ifeat = 1, CAVG_QUALITY_NFEATS
            call feature_auc_summary(dsets, ifeat, pooled_auc, mean_auc, min_auc, max_auc, inverted)
            write(funit,'(A,A,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,I0,A,ES14.6,A,ES14.6,A,ES14.6)') &
                'feature_signal,', trim(cavg_quality_feature_name(ifeat)), ',', pooled_auc, ',', mean_auc, ',', &
                min_auc, ',', max_auc, ',', inverted, ',', base_spec%weights(ifeat), ',', suggested_weights(ifeat), &
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
        pooled_auc        = pooled_feature_auc(dsets, ifeat)
        mean_auc          = 0.0
        min_auc           = huge(1.0)
        max_auc           = -huge(1.0)
        inverted_datasets = 0
        nused             = 0
        do ids = 1, size(dsets)
            if( dataset_learn_role(dsets(ids)) /= LEARN_ROLE_BALANCED ) cycle
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

    real function pooled_feature_auc( dsets, ifeat, holdout )
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer,                             intent(in) :: ifeat
        integer, optional,                   intent(in) :: holdout
        real,    allocatable :: vals(:)
        integer, allocatable :: refs(:)
        integer :: ids, nall, off, nfit
        if( present(holdout) )then
            nall = count_trainable_contrast_classes_all(dsets, holdout)
        else
            nall = count_trainable_contrast_classes_all(dsets)
        endif
        if( nall == 0 )then
            pooled_feature_auc = 0.5
            return
        endif
        allocate(vals(nall), refs(nall))
        off = 0
        do ids = 1, size(dsets)
            if( present(holdout) )then
                if( ids == holdout ) cycle
            endif
            if( dataset_learn_role(dsets(ids)) /= LEARN_ROLE_BALANCED ) cycle
            nfit = count_trainable_classes(dsets(ids))
            vals(off+1:off+nfit) = pack(dsets(ids)%features(:,ifeat), .not. dsets(ids)%hard_reject)
            refs(off+1:off+nfit) = pack(dsets(ids)%manual_states, .not. dsets(ids)%hard_reject)
            off = off + nfit
        end do
        pooled_feature_auc = auc_for_values(vals, refs)
        deallocate(vals, refs)
    end function pooled_feature_auc

    integer function count_trainable_contrast_classes_all( dsets, holdout )
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer, optional,                   intent(in) :: holdout
        integer :: ids
        count_trainable_contrast_classes_all = 0
        do ids = 1, size(dsets)
            if( present(holdout) )then
                if( ids == holdout ) cycle
            endif
            if( dataset_learn_role(dsets(ids)) /= LEARN_ROLE_BALANCED ) cycle
            count_trainable_contrast_classes_all = count_trainable_contrast_classes_all + count_trainable_classes(dsets(ids))
        end do
    end function count_trainable_contrast_classes_all

    subroutine write_feature_drop_diagnostics( funit, dsets, learned_model, best_score )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: best_score
        type(cavg_quality_model)      :: candidate
        type(cavg_quality_model_spec) :: spec
        integer :: ifeat
        real    :: score
        write(funit,'(A)') 'feature_drop_header=feature,learned_weight,macro_learn_score,delta_vs_learned'
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

    subroutine write_feature_policy_screen( funit, dsets )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer :: inds(CAVG_QUALITY_NFEATS)
        integer :: ipol, ninds
        write(funit,'(A)') 'feature_policy_lodo_header=feature_policy,n_features,mean_auc,min_auc,min_auc_dataset,'//&
            'mean_oracle_score,min_oracle_score,min_score_dataset,total_tp,total_fp,total_tn,total_fn'
        do ipol = 1, CAVG_QUALITY_LEARN_N_POLICIES
            call feature_policy_indices(ipol, inds, ninds)
            call write_feature_policy_lodo_row(funit, trim(feature_policy_name(ipol)), dsets, inds(1:ninds))
        end do
    end subroutine write_feature_policy_screen

    subroutine write_feature_policy_lodo_row( funit, policy_name, dsets, feat_inds )
        integer,                             intent(in) :: funit
        character(len=*),                    intent(in) :: policy_name
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        integer,                             intent(in) :: feat_inds(:)
        real, allocatable :: weights(:), scores(:)
        integer, allocatable :: refs(:)
        real :: auc, balacc, sum_auc, sum_balacc, min_auc, min_balacc
        integer :: holdout, tp, fp, tn, fn, total_tp, total_fp, total_tn, total_fn, nused, role
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
            role = dataset_learn_role(dsets(holdout))
            if( role == LEARN_ROLE_SKIP ) cycle
            call lodo_feature_policy_weights(dsets, holdout, feat_inds, weights)
            call score_dataset_feature_policy(dsets(holdout), feat_inds, weights, scores, refs)
            auc = auc_for_values(scores, refs)
            call best_score_threshold_balacc(scores, refs, role, balacc, tp, fp, tn, fn)
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
        if( nused > 0 )then
            low_auc_name = dataset_short_name(dsets(low_auc_id))
            low_balacc_name = dataset_short_name(dsets(low_balacc_id))
            write(funit,'(A,A,A,I0,A,F10.5,A,F10.5,A,A,A,F10.5,A,F10.5,A,A,A,I0,A,I0,A,I0,A,I0)') &
                'feature_policy_lodo,', trim(policy_name), ',', size(feat_inds), ',', sum_auc / real(nused), ',', &
                min_auc, ',', trim(low_auc_name), ',', sum_balacc / real(nused), ',', min_balacc, ',', &
                trim(low_balacc_name), ',', total_tp, ',', total_fp, ',', total_tn, ',', total_fn
        endif
        deallocate(weights)
        if( allocated(scores) ) deallocate(scores)
        if( allocated(refs)   ) deallocate(refs)
    end subroutine write_feature_policy_lodo_row

    subroutine lodo_feature_policy_weights( dsets, holdout, feat_inds, weights )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        integer,                             intent(in)  :: holdout
        integer,                             intent(in)  :: feat_inds(:)
        real,                                intent(out) :: weights(:)
        integer :: i
        real    :: auc
        if( size(weights) /= size(feat_inds) ) THROW_HARD('lodo_feature_policy_weights: invalid weight size')
        do i = 1, size(feat_inds)
            auc = pooled_feature_auc(dsets, feat_inds(i), holdout)
            weights(i) = max(0.0, auc - 0.5)
        end do
        if( sum(weights) > EPS )then
            weights = weights / sum(weights)
        else
            weights = 1.0 / real(size(weights))
        endif
    end subroutine lodo_feature_policy_weights

    subroutine score_dataset_feature_policy( dset, feat_inds, weights, scores, refs )
        type(cavg_quality_training_dataset), intent(in)  :: dset
        integer,                             intent(in)  :: feat_inds(:)
        real,                                intent(in)  :: weights(:)
        real, allocatable,                   intent(out) :: scores(:)
        integer, allocatable,                intent(out) :: refs(:)
        integer :: i, j, nfit, ifit
        if( size(feat_inds) /= size(weights) ) THROW_HARD('score_dataset_feature_policy: invalid weight size')
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
    end subroutine score_dataset_feature_policy

    subroutine best_score_threshold_balacc( scores, refs, role, best_balacc, tp, fp, tn, fn )
        real,    intent(in)  :: scores(:)
        integer, intent(in)  :: refs(:)
        integer, intent(in)  :: role
        real,    intent(out) :: best_balacc
        integer, intent(out) :: tp, fp, tn, fn
        real,    allocatable :: sorted_scores(:)
        integer, allocatable :: order(:)
        integer :: n, i, npos, nneg, cur_tp, cur_fp, cur_tn, cur_fn
        real    :: metric
        if( size(scores) /= size(refs) ) THROW_HARD('best_score_threshold_balacc: size mismatch')
        n  = size(scores)
        tp = 0; fp = 0; tn = 0; fn = 0
        if( n == 0 )then
            best_balacc = 0.5
            return
        endif
        npos = count(refs > 0)
        nneg = n - npos
        ! Baseline: predict none positive (threshold above max score)
        cur_tp = 0
        cur_fp = 0
        cur_fn = npos
        cur_tn = nneg
        best_balacc = learn_balacc_from_confusion(cur_tp, cur_fp, cur_tn, cur_fn, role)
        tp = cur_tp; fp = cur_fp; tn = cur_tn; fn = cur_fn
        ! Sort scores ascending, carrying the original indices in 'order'
        allocate(sorted_scores(n), source=scores)
        allocate(order(n))
        do i = 1, n
            order(i) = i
        end do
        call hpsort(sorted_scores, order)
        ! Sweep thresholds from high to low. Each step lowers the threshold to admit
        ! one more (or a tied group of) sample(s) as positive. This visits the same
        ! set of (tp,fp,tn,fn) states as the original O(n^2) loop, but in O(n log n).
        do i = n, 1, -1
            if( refs(order(i)) > 0 )then
                cur_tp = cur_tp + 1
                cur_fn = cur_fn - 1
            else
                cur_fp = cur_fp + 1
                cur_tn = cur_tn - 1
            endif
            ! Evaluate only at a tie-group boundary (last item, or next-lower sample
            ! has a strictly smaller score).
            if( i == 1 )then
                metric = learn_balacc_from_confusion(cur_tp, cur_fp, cur_tn, cur_fn, role)
                if( metric > best_balacc )then
                    best_balacc = metric
                    tp = cur_tp; fp = cur_fp; tn = cur_tn; fn = cur_fn
                endif
            else if( sorted_scores(i-1) < sorted_scores(i) )then
                metric = learn_balacc_from_confusion(cur_tp, cur_fp, cur_tn, cur_fn, role)
                if( metric > best_balacc )then
                    best_balacc = metric
                    tp = cur_tp; fp = cur_fp; tn = cur_tn; fn = cur_fn
                endif
            endif
        end do
        deallocate(sorted_scores, order)
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
        integer :: ids, tp, fp, tn, fn, role
        logical :: lowsep, single_cluster, otsu_like, rescue_like, min_accept_like
        diag = cavg_quality_learn_diagnostics()
        diag%n_datasets = size(dsets)
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            select case(role)
                case(LEARN_ROLE_BALANCED)
                    diag%n_scored_datasets = diag%n_scored_datasets + 1
                    diag%n_weight_datasets = diag%n_weight_datasets + 1
                case(LEARN_ROLE_RECALL_ONLY)
                    diag%n_scored_datasets = diag%n_scored_datasets + 1
                    diag%n_recall_only = diag%n_recall_only + 1
                case(LEARN_ROLE_SPECIFICITY_ONLY)
                    diag%n_scored_datasets = diag%n_scored_datasets + 1
                    diag%n_specificity_only = diag%n_specificity_only + 1
                case default
                    diag%n_skipped = diag%n_skipped + 1
                    cycle
            end select
            call classify_training_dataset_detail(dsets(ids), model, quality, tp, fp, tn, fn)
            diag%total_fp = diag%total_fp + fp
            diag%total_fn = diag%total_fn + fn
            lowsep = quality%separation < model%min_score_separation
            single_cluster = trim(quality%soft_decision) == 'hard_only' .or. &
                quality%nclust <= 1 .or. .not. quality%used_threshold
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
            threshold_policy_looks_otsu = abs(quality%threshold_offset - model%boundary_margin) > 1.0e-4
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
        min_accept_fraction_looks_active = quality%threshold_offset > model%boundary_margin + 1.0e-4
    end function min_accept_fraction_looks_active

    subroutine write_learn_search_diagnostics( funit, learned_model, diag, best_tie_specs, n_best_ties )
        integer,                             intent(in) :: funit
        type(cavg_quality_model),            intent(in) :: learned_model
        type(cavg_quality_learn_diagnostics),intent(in) :: diag
        type(cavg_quality_model_spec),       intent(in) :: best_tie_specs(:)
        integer,                             intent(in) :: n_best_ties
        character(len=LONGSTRLEN) :: detail
        write(funit,'(A)') 'search_diagnostic_header=level,parameter,status,detail'
        write(detail,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'scored=', diag%n_scored_datasets, &
            ';weight_contrast=', diag%n_weight_datasets, ';recall_only=', diag%n_recall_only, &
            ';specificity_only=', diag%n_specificity_only, ';skipped=', diag%n_skipped
        call write_search_diagnostic(funit, 'note', 'dataset_roles', 'automatic', trim(detail))
        write(detail,'(A,F6.3,A,F6.3)') 'floor=', LEARN_RECALL_ONLY_FLOOR, ';penalty=', &
            LEARN_RECALL_ONLY_PENALTY
        call write_search_diagnostic(funit, 'note', 'trainable_good_only_score', 'guarded_recall', trim(detail))
        call write_search_diagnostic(funit, 'note', 'feature_weights', 'auc_no_base_blending', &
            'one_auc_derived_candidate_per_feature_policy')
        call write_search_diagnostic(funit, 'note', 'feature_policy', trim(learned_model%feature_policy), &
            'selected_family_set_encoded_by_zeroed_model_weights')
        call write_minsep_diagnostic(funit, learned_model%min_score_separation, best_tie_specs, n_best_ties)
        if( model_is_pool_context(learned_model) )then
            call write_grid_position_diagnostic(funit, 'boundary_margin', learned_model%boundary_margin, &
                LEARN_POOL_MARGINS, 'best_at_lowest_value_consider_more_negative_if_junk_leaks_after_validation', &
                'best_at_highest_value_consider_more_positive_if_good_classes_are_rejected')
        else
            call write_grid_position_diagnostic(funit, 'boundary_margin', learned_model%boundary_margin, &
                LEARN_CHUNK_MARGINS, 'best_at_lowest_value_consider_more_negative_if_junk_leaks_after_validation', &
                'best_at_highest_value_consider_more_positive_if_good_classes_are_rejected')
        endif
        call write_bool_grid_diagnostic(funit, 'use_lowsep_otsu', learned_model%use_lowsep_otsu, &
            'searched', 'learned_from_training_grid')
        call write_bool_grid_diagnostic(funit, 'use_otsu_window', learned_model%use_otsu_window, &
            'searched', 'learned_from_training_grid')
        if( learned_model%use_otsu_window )then
            call write_grid_position_diagnostic(funit, 'otsu_min_offset', learned_model%otsu_min_offset, &
                LEARN_OTSU_MIN_OFFSETS, 'best_at_lowest_value_consider_smaller_otsu_window_minimum', &
                'best_at_highest_value_consider_larger_otsu_window_minimum')
            call write_grid_position_diagnostic(funit, 'otsu_max_offset', learned_model%otsu_max_offset, &
                LEARN_OTSU_MAX_OFFSETS, 'best_at_lowest_value_consider_smaller_otsu_window_maximum', &
                'best_at_highest_value_consider_larger_otsu_window_maximum')
        else
            call write_search_diagnostic(funit, 'note', 'otsu_min_offset', 'inactive', &
                'use_otsu_window is false in the learned model')
            call write_search_diagnostic(funit, 'note', 'otsu_max_offset', 'inactive', &
                'use_otsu_window is false in the learned model')
        endif
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

    subroutine write_evaluate_diagnostics( funit, model, diag )
        integer,                              intent(in) :: funit
        type(cavg_quality_model),             intent(in) :: model
        type(cavg_quality_learn_diagnostics), intent(in) :: diag
        character(len=LONGSTRLEN) :: detail
        write(funit,'(A)') 'evaluate_diagnostic_header=level,parameter,status,detail'
        write(detail,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'scored=', diag%n_scored_datasets, &
            ';weight_contrast=', diag%n_weight_datasets, ';recall_only=', diag%n_recall_only, &
            ';specificity_only=', diag%n_specificity_only, ';skipped=', diag%n_skipped
        call write_evaluate_diagnostic(funit, 'note', 'dataset_roles', 'automatic', trim(detail))
        write(detail,'(A,F6.3,A,F6.3)') 'floor=', LEARN_RECALL_ONLY_FLOOR, ';penalty=', &
            LEARN_RECALL_ONLY_PENALTY
        call write_evaluate_diagnostic(funit, 'note', 'trainable_good_only_score', 'guarded_recall', trim(detail))
        call write_evaluate_diagnostic(funit, 'note', 'feature_policy', trim(model%feature_policy), 'fixed_model')
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%use_lowsep_otsu, ';active_datasets=', &
            diag%n_lowsep, ';fp=', diag%lowsep_fp, ';fn=', diag%lowsep_fn
        call write_evaluate_diagnostic(funit, policy_level(diag%n_lowsep, diag%lowsep_fp, diag%lowsep_fn), &
            'use_lowsep_otsu_effect', 'fixed_model', trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%use_otsu_window, ';otsu_like_datasets=', &
            diag%n_otsu_like, ';fp=', diag%otsu_like_fp, ';fn=', diag%otsu_like_fn
        call write_evaluate_diagnostic(funit, policy_level(diag%n_otsu_like, diag%otsu_like_fp, diag%otsu_like_fn), &
            'use_otsu_window_effect', 'fixed_model', trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%use_cluster_rescue, ';rescue_like_datasets=', &
            diag%n_rescue_like, ';fp=', diag%rescue_like_fp, ';fn=', diag%rescue_like_fn
        call write_evaluate_diagnostic(funit, rescue_policy_level(model, diag), 'use_cluster_rescue', 'fixed_model', &
            trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%enforce_min_accept_frac, ';min_accept_datasets=', &
            diag%n_min_accept_like, ';fp=', diag%min_accept_like_fp, ';fn=', diag%min_accept_like_fn
        call write_evaluate_diagnostic(funit, min_accept_policy_level(model, diag), 'enforce_min_accept_frac', &
            'fixed_model', trim(detail))
        write(detail,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'single_cluster_datasets=', diag%n_single_cluster, &
            ';fp=', diag%single_cluster_fp, ';fn=', diag%single_cluster_fn, ';total_fp=', diag%total_fp, &
            ';total_fn=', diag%total_fn
        call write_evaluate_diagnostic(funit, policy_level(diag%n_single_cluster, diag%single_cluster_fp, &
            diag%single_cluster_fn), 'accept_all_fallback', 'fixed_model', trim(detail))
    end subroutine write_evaluate_diagnostics

    subroutine write_minsep_diagnostic( funit, val, best_tie_specs, n_best_ties )
        integer,                       intent(in) :: funit, n_best_ties
        real,                          intent(in) :: val
        type(cavg_quality_model_spec), intent(in) :: best_tie_specs(:)
        character(len=LONGSTRLEN) :: detail
        real :: tie_min, tie_max
        call best_tie_minsep_span(best_tie_specs, n_best_ties, tie_min, tie_max)
        if( n_best_ties > 1 .and. real_close(tie_min, minval(LEARN_MINSEPS)) .and. &
            real_close(tie_max, maxval(LEARN_MINSEPS)) )then
            write(detail,'(A,ES14.6,A,ES14.6,A,ES14.6)') 'best_tie_min=', tie_min, ';best_tie_max=', &
                tie_max, ';selected_by_tie_breaker=', val
            call write_search_diagnostic(funit, 'note', 'min_score_separation', 'flat_across_grid', trim(detail))
        else
            call write_grid_position_diagnostic(funit, 'min_score_separation', val, LEARN_MINSEPS, &
                'best_at_lowest_value_consider_lower_if_accept_all_remains_too_common', &
                'best_at_highest_value_consider_higher_if_unstable_splits_remain')
        endif
    end subroutine write_minsep_diagnostic

    subroutine best_tie_minsep_span( best_tie_specs, n_best_ties, tie_min, tie_max )
        type(cavg_quality_model_spec), intent(in)  :: best_tie_specs(:)
        integer,                       intent(in)  :: n_best_ties
        real,                          intent(out) :: tie_min, tie_max
        integer :: i, nstored
        nstored = min(n_best_ties, size(best_tie_specs))
        tie_min = huge(1.0)
        tie_max = -huge(1.0)
        do i = 1, nstored
            tie_min = min(tie_min, best_tie_specs(i)%min_score_separation)
            tie_max = max(tie_max, best_tie_specs(i)%min_score_separation)
        end do
        if( nstored == 0 )then
            tie_min = 0.0
            tie_max = 0.0
        endif
    end subroutine best_tie_minsep_span

    subroutine write_bool_grid_diagnostic( funit, param, val, status, detail )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: param, status, detail
        logical,          intent(in) :: val
        if( val )then
            call write_search_diagnostic(funit, 'note', param, status//'_selected_true', detail)
        else
            call write_search_diagnostic(funit, 'note', param, status//'_selected_false', detail)
        endif
    end subroutine write_bool_grid_diagnostic

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
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%use_lowsep_otsu, ';active_datasets=', &
            diag%n_lowsep, ';fp=', diag%lowsep_fp, ';fn=', diag%lowsep_fn
        call write_search_diagnostic(funit, policy_level(diag%n_lowsep, diag%lowsep_fp, diag%lowsep_fn), &
            'use_lowsep_otsu_effect', 'searched', trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%use_otsu_window, ';otsu_like_datasets=', &
            diag%n_otsu_like, ';fp=', diag%otsu_like_fp, ';fn=', diag%otsu_like_fn
        call write_search_diagnostic(funit, policy_level(diag%n_otsu_like, diag%otsu_like_fp, diag%otsu_like_fn), &
            'use_otsu_window_effect', 'searched', trim(detail))
        write(detail,'(A,F8.4,A,I0,A,I0,A,I0)') 'selected=', model%otsu_min_offset, ';otsu_like_datasets=', &
            diag%n_otsu_like, ';fp=', diag%otsu_like_fp, ';fn=', diag%otsu_like_fn
        call write_search_diagnostic(funit, policy_level(diag%n_otsu_like, diag%otsu_like_fp, diag%otsu_like_fn), &
            'otsu_min_offset_effect', 'searched', trim(detail))
        write(detail,'(A,F8.4,A,I0,A,I0,A,I0)') 'selected=', model%otsu_max_offset, ';otsu_like_datasets=', &
            diag%n_otsu_like, ';fp=', diag%otsu_like_fp, ';fn=', diag%otsu_like_fn
        call write_search_diagnostic(funit, policy_level(diag%n_otsu_like, diag%otsu_like_fp, diag%otsu_like_fn), &
            'otsu_max_offset_effect', 'searched', trim(detail))
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

    subroutine write_evaluate_diagnostic( funit, level, param, status, detail )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: level, param, status, detail
        write(funit,'(8A)') 'evaluate_diagnostic,', trim(level), ',', trim(param), ',', trim(status), ',', trim(detail)
    end subroutine write_evaluate_diagnostic

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

    subroutine write_learn_margin_grid( funit, is_pool_context )
        integer, intent(in) :: funit
        logical, intent(in) :: is_pool_context
        if( is_pool_context )then
            call write_real_list(funit, 'grid_boundary_margins=', LEARN_POOL_MARGINS)
        else
            call write_real_list(funit, 'grid_boundary_margins=', LEARN_CHUNK_MARGINS)
        endif
    end subroutine write_learn_margin_grid

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

    subroutine write_logical_list( funit, key, vals )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: key
        logical,          intent(in) :: vals(:)
        integer :: i
        write(funit,'(A)', advance='no') trim(key)
        do i = 1, size(vals)
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(L1)', advance='no') vals(i)
        end do
        write(funit,*)
    end subroutine write_logical_list

    subroutine write_candidate_table_header( funit, key )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: key
        write(funit,'(A)') trim(key)//&
            'rank,learn_score,feature_policy,boundary_margin,min_score_separation,min_accept_frac,'//&
            'use_lowsep_otsu,use_otsu_window,otsu_min_offset,otsu_max_offset,use_cluster_rescue,'//&
            'enforce_min_accept_frac,feature_weights_semicolon'
    end subroutine write_candidate_table_header

    subroutine write_candidate_row( funit, tag, irank, score, spec )
        integer,                       intent(in) :: funit, irank
        character(len=*),              intent(in) :: tag
        real,                          intent(in) :: score
        type(cavg_quality_model_spec), intent(in) :: spec
        integer :: i
        write(funit,'(A,A,I0,A,F10.5,A,A,A,ES14.6,A,ES14.6,A,ES14.6,A,L1,A,L1,A,ES14.6,A,ES14.6,A,L1,A,L1,A)', &
            advance='no') &
            trim(tag), ',', irank, ',', score, ',', trim(spec%feature_policy), ',', spec%boundary_margin, ',', &
            spec%min_score_separation, ',', spec%min_accept_frac, ',', spec%use_lowsep_otsu, ',', &
            spec%use_otsu_window, ',', spec%otsu_min_offset, ',', spec%otsu_max_offset, ',', &
            spec%use_cluster_rescue, ',', spec%enforce_min_accept_frac, ','
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
