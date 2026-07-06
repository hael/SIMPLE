!@descr: learn-mode training-table reader and model search for class-average quality
module simple_cavg_quality_learn
use simple_defs,               only: logfhandle, LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string,             only: string
use simple_string_utils,       only: str2int, str2real, str_is_true, csv_field
use simple_cavg_quality_feats, only: cavg_quality_feature_name, cavg_quality_feature_family, &
    I_NEG_LOG_RES, I_NEG_LOCVAR_FG, I_CC_AREA_FRAC, I_BP40_100_CENTER_EDGE_VAR
use simple_cavg_quality_model, only: cavg_quality_model, cavg_quality_classify_cache, &
    build_classify_cache, kill_classify_cache, cached_decision_confusion, apply_cached_decision_to_quality
use simple_cavg_quality_stats, only: calc_confusion, calc_binary_metrics, auc_for_values
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, CAVG_QUALITY_MAX_INTERACTIONS, EPS, CAVG_MODEL_FAMILY_LINEAR, &
    CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC, cavg_quality_model_spec, cavg_quality_result, cavg_quality_training_dataset, &
    cavg_quality_learn_diagnostics
use simple_srch_sort_loc,      only: hpsort
use simple_optimizer,          only: optimizer
use simple_opt_factory,        only: opt_factory
use simple_opt_spec,           only: opt_spec
implicit none
private
#include "simple_local_flags.inc"

public :: learn_cavg_quality_model
public :: evaluate_cavg_quality_model
public :: evaluate_cavg_quality_result

logical, parameter :: LEARN_OTSU_FLAGS(2)           = [.false., .true.]
integer, parameter :: CAVG_QUALITY_LEARN_TOP_K      = 10
integer, parameter :: CAVG_QUALITY_LEARN_MAX_TIES   = 20000
integer, parameter :: CAVG_QUALITY_LEARN_N_STANDARD_POLICIES = 4
integer, parameter :: LEARN_ROLE_SKIP               = 0
integer, parameter :: LEARN_ROLE_BALANCED           = 1
integer, parameter :: LEARN_ROLE_RECALL_ONLY        = 2
integer, parameter :: LEARN_ROLE_SPECIFICITY_ONLY   = 3
real,    parameter :: LEARN_RECALL_ONLY_FLOOR        = 0.95
real,    parameter :: LEARN_RECALL_ONLY_PENALTY      = 3.0
real,    parameter :: LEARN_BALANCED_FN_TOLERANCE_FRAC = 0.02
real,    parameter :: LEARN_BALANCED_FN_RATE_PENALTY = 10.0
real,    parameter :: LEARN_ROBUST_TAIL_FRAC         = 0.25
real,    parameter :: LEARN_ROBUST_TAIL_WEIGHT       = 0.50
real,    parameter :: LEARN_BALANCED_GOOD_LOSS_WEIGHT = 0.65
real,    parameter :: LEARN_BALANCED_BAD_LOSS_WEIGHT  = 1.0 - LEARN_BALANCED_GOOD_LOSS_WEIGHT
real,    parameter :: LEARN_BAD_OVERFIT_LOSS_MULT     = 1.50
real,    parameter :: LEARN_OVERFIT_FOCUS_BAD_LOSS_SCALE = 0.00
real,    parameter :: LEARN_BAD_OVERFIT_LOWVAR_FG_MIN = 0.0
real,    parameter :: LEARN_BAD_OVERFIT_SUPPORT_MAX   = 0.5
real,    parameter :: LEARN_OVERFIT_FP_FOCUS_MIN_TRAINABLE_FRAC = 0.20
real,    parameter :: LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET = 0.20
real,    parameter :: LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY = 8.0
real,    parameter :: LEARN_OVERFIT_FP_BP_Z_MAX       = 0.0
real,    parameter :: LEARN_MINSEPS(7)              = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
! Positive margins deliberately over-select relative to the learned boundary.
real,    parameter :: LEARN_BOUNDARY_MARGINS(37)     = [-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, &
                                                       -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, &
                                                        0.25, 0.30, 0.40, 0.50, 0.60, 0.70, &
                                                        0.80, 0.90, 1.00, 1.10, 1.20, 1.30, &
                                                        1.40, 1.50, 1.60, 1.70, 1.80, 1.90, &
                                                        2.00, 2.25, 2.50, 2.75, 3.00, 3.50, &
                                                        4.00]
real,    parameter :: LEARN_OTSU_MIN_OFFSETS(9)     = [0.05, 0.10, 0.15, 0.25, 0.35, 0.40, 0.45, 0.50, 0.60]
real,    parameter :: LEARN_OTSU_MAX_OFFSETS(9)     = [0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.75, 0.90]
real,    parameter :: LEARN_MIN_ACCEPT_FRACS(11)    = [0.50, 0.60, 0.65, 0.70, 0.80, 0.85, &
                                                       0.90, 0.925, 0.95, 0.975, 1.00]
real,    parameter :: LOGISTIC_LAMBDAS(7)           = [1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, &
                                                       1.0e-2, 3.0e-2, 1.0e-1]
real,    parameter :: LOGISTIC_THRESHOLDS(20)       = [0.02, 0.05, 0.075, 0.10, 0.15, 0.20, &
                                                       0.25, 0.30, 0.35, 0.40, 0.45, 0.50, &
                                                       0.55, 0.60, 0.65, 0.70, 0.75, 0.80, &
                                                       0.85, 0.90]
real,    parameter :: LOGISTIC_COEFF_ABS_BOUND      = 20.0

type :: cavg_quality_logistic_problem
    integer :: nobs           = 0
    integer :: ndim           = 0
    integer :: n_linear       = 0
    integer :: n_interactions = 0
    integer :: linear_features(CAVG_QUALITY_NFEATS) = 0
    integer :: interaction_terms(CAVG_QUALITY_MAX_INTERACTIONS,2) = 0
    real(kind=8) :: lambda       = 0.0d0
    real(kind=8) :: total_weight = 0.0d0
    real(kind=8), allocatable :: design(:,:)
    real(kind=8), allocatable :: labels(:)
    real(kind=8), allocatable :: weights(:)
contains
    procedure :: kill => kill_logistic_problem
end type cavg_quality_logistic_problem

contains

    subroutine learn_cavg_quality_model( analysis_files, learned_model, model_fname, report_fname, model_family, &
                                         trust_resolution )
        class(string),             intent(in)    :: analysis_files(:)
        type(cavg_quality_model),  intent(inout) :: learned_model
        character(len=*),          intent(in)    :: model_fname, report_fname
        character(len=*), optional,intent(in)    :: model_family
        logical,          optional,intent(in)    :: trust_resolution
        type(cavg_quality_training_dataset), allocatable :: dsets(:)
        type(cavg_quality_classify_cache),   allocatable :: caches(:)
        type(cavg_quality_model_spec) :: base_spec, candidate_spec, best_spec
        type(cavg_quality_model_spec) :: top_specs(CAVG_QUALITY_LEARN_TOP_K)
        type(cavg_quality_model_spec), allocatable :: best_tie_specs(:)
        real :: suggested_weights(CAVG_QUALITY_NFEATS)
        real :: top_scores(CAVG_QUALITY_LEARN_TOP_K)
        logical :: feature_mask(CAVG_QUALITY_NFEATS)
        real :: best_score, learn_score
        integer :: ipol, im, isep, ilow, iwin, iomin, iomax, max_grid
        integer :: n_grid, n_top, n_best_ties
        character(len=32) :: requested_family
        call init_learn_feature_mask(feature_mask, trust_resolution)
        requested_family = 'logistic'
        if( present(model_family) ) requested_family = trim(model_family)
        select case(trim(requested_family))
        case('linear', CAVG_MODEL_FAMILY_LINEAR)
            continue
        case('logistic', CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC)
            call learn_cavg_quality_pairwise_logistic_model(analysis_files, learned_model, model_fname, report_fname, &
                feature_mask)
            return
        case default
            THROW_HARD('learn_cavg_quality_model: unsupported model family: '//trim(requested_family))
        end select
        call load_quality_training_datasets(analysis_files, dsets)
        base_spec = abinitio_learn_base_spec()
        call calc_suggested_training_weights(dsets, suggested_weights, feature_mask)
        best_spec  = base_spec
        best_score = -huge(1.0)
        top_scores = -huge(1.0)
        n_top       = 0
        n_grid      = 0
        n_best_ties = 0
        max_grid = n_feature_policies() * size(LEARN_MINSEPS) * &
            size(LEARN_BOUNDARY_MARGINS) * n_otsu_grid_combinations() * size(LEARN_OTSU_FLAGS) * &
            (1 + size(LEARN_MIN_ACCEPT_FRACS))
        if( max_grid > 10000 ) write(logfhandle,'(A,I0)') &
            '>>> CAVG QUALITY LEARN CANDIDATE GRID SIZE: ', max_grid
        allocate(best_tie_specs(min(max_grid, CAVG_QUALITY_LEARN_MAX_TIES)))
        allocate(caches(size(dsets)))
        do ipol = 1, n_feature_policies()
            candidate_spec = base_spec
            candidate_spec%feature_policy = feature_policy_name(ipol)
            candidate_spec%weights = suggested_weights
            call apply_feature_policy(ipol, candidate_spec%weights, feature_mask)
            ! For a fixed weight vector, the per-dataset scores,
            ! k-medoids partition, raw threshold, and Otsu threshold are
            ! reused across the whole threshold-control sub-grid.
            call build_policy_caches(dsets, candidate_spec%weights, caches)
            do isep = 1, size(LEARN_MINSEPS)
                candidate_spec%min_score_separation = LEARN_MINSEPS(isep)
                do im = 1, size(LEARN_BOUNDARY_MARGINS)
                    candidate_spec%boundary_margin = LEARN_BOUNDARY_MARGINS(im)
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
                                        call evaluate_policy_grid(dsets, caches, candidate_spec, &
                                            n_grid, best_spec, best_score, best_tie_specs, n_best_ties, &
                                            top_specs, top_scores, n_top)
                                    end do
                                end do
                            else
                                candidate_spec%otsu_min_offset = base_spec%otsu_min_offset
                                candidate_spec%otsu_max_offset = base_spec%otsu_max_offset
                                call evaluate_policy_grid(dsets, caches, candidate_spec, &
                                    n_grid, best_spec, best_score, best_tie_specs, n_best_ties, &
                                    top_specs, top_scores, n_top)
                            endif
                        end do
                    end do
                end do
            end do
        end do
        call select_preferred_best_tie(base_spec, best_tie_specs, n_best_ties, best_spec)
        best_spec%name = 'learned_v1'
        call learned_model%init_spec(best_spec)
        learn_score = macro_balacc_for_model(dsets, learned_model)
        call learned_model%write(model_fname)
        call write_cavg_quality_learn_report(report_fname, dsets, base_spec, suggested_weights, learned_model, &
            learn_score, n_grid, top_specs, top_scores, n_top, best_tie_specs, n_best_ties, feature_mask)
        call kill_policy_caches(caches)
        deallocate(caches)
        call kill_training_datasets(dsets)
        deallocate(best_tie_specs)
        deallocate(dsets)
    end subroutine learn_cavg_quality_model

    subroutine learn_cavg_quality_pairwise_logistic_model( analysis_files, learned_model, model_fname, report_fname, &
                                                           feature_mask )
        class(string),             intent(in)    :: analysis_files(:)
        type(cavg_quality_model),  intent(inout) :: learned_model
        character(len=*),          intent(in)    :: model_fname, report_fname
        logical,                   intent(in)    :: feature_mask(CAVG_QUALITY_NFEATS)
        type(cavg_quality_training_dataset), allocatable :: dsets(:)
        type(cavg_quality_logistic_problem) :: problem
        type(cavg_quality_model) :: candidate, best_model
        type(cavg_quality_model_spec) :: spec
        real(kind=8), allocatable :: solution(:)
        real :: score, best_score, objective, best_objective, learn_score
        integer :: ipol, ilambda, ithresh, n_candidates
        call load_quality_training_datasets(analysis_files, dsets)
        best_score      = -huge(1.0)
        best_objective  = huge(1.0)
        n_candidates    = 0
        do ipol = 1, n_feature_policies()
            do ilambda = 1, size(LOGISTIC_LAMBDAS)
                call build_logistic_problem(dsets, ipol, LOGISTIC_LAMBDAS(ilambda), problem, feature_mask)
                call fit_logistic_problem(problem, solution, objective)
                call logistic_solution_to_model(solution, problem, feature_policy_name(ipol), candidate)
                do ithresh = 1, size(LOGISTIC_THRESHOLDS)
                    spec = candidate%get_spec()
                    spec%prob_threshold = LOGISTIC_THRESHOLDS(ithresh)
                    call candidate%init_spec(spec)
                    score = macro_balacc_for_model(dsets, candidate)
                    n_candidates = n_candidates + 1
                    if( score > best_score + EPS .or. &
                        (abs(score - best_score) <= EPS .and. objective < best_objective) )then
                        best_score     = score
                        best_objective = objective
                        best_model     = candidate
                    endif
                end do
                if( allocated(solution) ) deallocate(solution)
                call problem%kill()
            end do
        end do
        if( n_candidates == 0 ) THROW_HARD('learn_cavg_quality_pairwise_logistic_model: no logistic candidates')
        learned_model = best_model
        learn_score = macro_balacc_for_model(dsets, learned_model)
        call learned_model%write(model_fname)
        call write_cavg_quality_logistic_learn_report(report_fname, dsets, learned_model, learn_score, &
            n_candidates, best_objective, feature_mask)
        call kill_training_datasets(dsets)
        deallocate(dsets)
    end subroutine learn_cavg_quality_pairwise_logistic_model

    subroutine build_logistic_problem( dsets, ipolicy, lambda, problem, feature_mask )
        type(cavg_quality_training_dataset), intent(in)    :: dsets(:)
        integer,                             intent(in)    :: ipolicy
        real,                                intent(in)    :: lambda
        type(cavg_quality_logistic_problem), intent(inout) :: problem
        logical,                             intent(in)    :: feature_mask(CAVG_QUALITY_NFEATS)
        integer :: ids, irow, iobs, ifeat, jfeat, ilinear, jlinear, iterm
        integer :: role, nfit, ngood, nbad
        real(kind=8) :: good_weight, bad_weight, focus_bad_loss_mult
        call problem%kill()
        call feature_policy_indices(ipolicy, problem%linear_features, problem%n_linear, feature_mask)
        if( problem%n_linear < 1 ) THROW_HARD('build_logistic_problem: empty feature policy')
        problem%n_interactions = 0
        do ilinear = 1, problem%n_linear - 1
            do jlinear = ilinear + 1, problem%n_linear
                problem%n_interactions = problem%n_interactions + 1
                problem%interaction_terms(problem%n_interactions,1) = problem%linear_features(ilinear)
                problem%interaction_terms(problem%n_interactions,2) = problem%linear_features(jlinear)
            end do
        end do
        problem%ndim = 1 + problem%n_linear + problem%n_interactions
        problem%lambda = real(lambda, kind=8)
        problem%nobs = 0
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            if( role == LEARN_ROLE_SKIP ) cycle
            problem%nobs = problem%nobs + count_trainable_classes(dsets(ids))
        end do
        if( problem%nobs == 0 ) THROW_HARD('build_logistic_problem: no trainable class averages')
        allocate(problem%design(problem%nobs, problem%ndim), source=0.0d0)
        allocate(problem%labels(problem%nobs), source=0.0d0)
        allocate(problem%weights(problem%nobs), source=0.0d0)
        iobs = 0
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            if( role == LEARN_ROLE_SKIP ) cycle
            ngood = count_trainable_manual_good(dsets(ids))
            nbad  = count_trainable_manual_bad(dsets(ids))
            call logistic_dataset_class_weights(role, ngood, nbad, good_weight, bad_weight)
            focus_bad_loss_mult = logistic_overfit_focus_bad_loss_multiplier(dsets(ids), role)
            do irow = 1, dsets(ids)%ncls
                if( dsets(ids)%hard_reject(irow) ) cycle
                iobs = iobs + 1
                problem%design(iobs,1) = 1.0d0
                do ilinear = 1, problem%n_linear
                    ifeat = problem%linear_features(ilinear)
                    problem%design(iobs, 1 + ilinear) = real(dsets(ids)%features(irow, ifeat), kind=8)
                end do
                do iterm = 1, problem%n_interactions
                    ifeat = problem%interaction_terms(iterm,1)
                    jfeat = problem%interaction_terms(iterm,2)
                    problem%design(iobs, 1 + problem%n_linear + iterm) = &
                        real(dsets(ids)%features(irow, ifeat) * dsets(ids)%features(irow, jfeat), kind=8)
                end do
                if( dsets(ids)%manual_states(irow) > 0 )then
                    problem%labels(iobs)  = 1.0d0
                    problem%weights(iobs) = good_weight
                else
                    problem%labels(iobs)  = 0.0d0
                    problem%weights(iobs) = bad_weight * logistic_bad_example_loss_multiplier( &
                        dsets(ids)%features(irow,:), focus_bad_loss_mult)
                endif
            end do
            nfit = count_trainable_classes(dsets(ids))
            if( nfit > 0 .and. iobs > problem%nobs ) THROW_HARD('build_logistic_problem: row overflow')
        end do
        problem%total_weight = sum(problem%weights)
        if( problem%total_weight <= 0.0d0 ) THROW_HARD('build_logistic_problem: nonpositive total weight')
    end subroutine build_logistic_problem

    subroutine logistic_dataset_class_weights( role, ngood, nbad, good_weight, bad_weight )
        integer,      intent(in)  :: role, ngood, nbad
        real(kind=8), intent(out) :: good_weight, bad_weight
        good_weight = 0.0d0
        bad_weight  = 0.0d0
        select case(role)
        case(LEARN_ROLE_BALANCED)
            if( ngood <= 0 .or. nbad <= 0 ) THROW_HARD('logistic_dataset_class_weights: invalid balanced dataset')
            good_weight = real(LEARN_BALANCED_GOOD_LOSS_WEIGHT, kind=8) / real(ngood, kind=8)
            bad_weight  = real(LEARN_BALANCED_BAD_LOSS_WEIGHT,  kind=8) / real(nbad,  kind=8)
        case(LEARN_ROLE_RECALL_ONLY)
            if( ngood <= 0 ) THROW_HARD('logistic_dataset_class_weights: invalid recall-only dataset')
            good_weight = 1.0d0 / real(ngood, kind=8)
        case(LEARN_ROLE_SPECIFICITY_ONLY)
            if( nbad <= 0 ) THROW_HARD('logistic_dataset_class_weights: invalid specificity-only dataset')
            bad_weight = 1.0d0 / real(nbad, kind=8)
        case default
            THROW_HARD('logistic_dataset_class_weights: unsupported dataset role')
        end select
    end subroutine logistic_dataset_class_weights

    real(kind=8) function logistic_bad_example_loss_multiplier( z_features, focus_bad_loss_mult )
        real,         intent(in) :: z_features(:)
        real(kind=8), intent(in) :: focus_bad_loss_mult
        if( size(z_features) /= CAVG_QUALITY_NFEATS ) &
            THROW_HARD('logistic_bad_example_loss_multiplier: invalid feature count')
        logistic_bad_example_loss_multiplier = 1.0d0
        if( bad_overfit_signature(z_features) )then
            logistic_bad_example_loss_multiplier = real(LEARN_BAD_OVERFIT_LOSS_MULT, kind=8) * &
                focus_bad_loss_mult
        endif
    end function logistic_bad_example_loss_multiplier

    real(kind=8) function logistic_overfit_focus_bad_loss_multiplier( dset, role )
        type(cavg_quality_training_dataset), intent(in) :: dset
        integer,                             intent(in) :: role
        integer :: nfit, nbad_overfit
        real :: signature_frac
        logistic_overfit_focus_bad_loss_multiplier = 1.0d0
        if( role /= LEARN_ROLE_BALANCED .and. role /= LEARN_ROLE_SPECIFICITY_ONLY ) return
        nfit = count_trainable_classes(dset)
        if( nfit <= 0 ) return
        nbad_overfit = count_trainable_bad_overfit_signature(dset)
        signature_frac = real(nbad_overfit) / real(nfit)
        if( signature_frac < LEARN_OVERFIT_FP_FOCUS_MIN_TRAINABLE_FRAC ) return
        logistic_overfit_focus_bad_loss_multiplier = 1.0d0 + &
            real(LEARN_OVERFIT_FOCUS_BAD_LOSS_SCALE * signature_frac, kind=8)
    end function logistic_overfit_focus_bad_loss_multiplier

    subroutine fit_logistic_problem( problem, solution, objective )
        type(cavg_quality_logistic_problem), intent(inout) :: problem
        real(kind=8), allocatable,           intent(inout) :: solution(:)
        real,                                intent(out)   :: objective
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec
        class(optimizer), pointer :: opt_obj
        real, allocatable         :: limits(:,:)
        real(kind=8)              :: pos_weight, neg_weight
        integer :: i
        if( allocated(solution) ) deallocate(solution)
        allocate(solution(problem%ndim), source=0.0d0)
        allocate(limits(problem%ndim,2), source=0.0)
        limits(:,1) = -LOGISTIC_COEFF_ABS_BOUND
        limits(:,2) =  LOGISTIC_COEFF_ABS_BOUND
        pos_weight = sum(problem%weights, mask=problem%labels > 0.5d0)
        neg_weight = sum(problem%weights, mask=problem%labels <= 0.5d0)
        solution(1) = log(max(pos_weight, 1.0d-12) / max(neg_weight, 1.0d-12))
        opt_obj => null()
        call ospec%specify('lbfgsb', problem%ndim, ftol=1e-7, gtol=1e-7, maxits=250, &
            factr=1.0d+7, pgtol=1.0d-5, limits=limits)
        call ofac%new(ospec, opt_obj)
        call ospec%set_costfun_8(logistic_cost_wrapper)
        call ospec%set_gcostfun_8(logistic_gradient_wrapper)
        call ospec%set_fdfcostfun_8(logistic_fdf_wrapper)
        ospec%x   = real(solution)
        ospec%x_8 = solution
        call opt_obj%minimize(ospec, problem, objective)
        do i = 1, problem%ndim
            solution(i) = ospec%x_8(i)
        end do
        call opt_obj%kill()
        deallocate(opt_obj)
        call ospec%kill()
        deallocate(limits)
    end subroutine fit_logistic_problem

    subroutine logistic_solution_to_model( solution, problem, feature_policy, model )
        real(kind=8),                      intent(in)    :: solution(:)
        type(cavg_quality_logistic_problem), intent(in)  :: problem
        character(len=*),                  intent(in)    :: feature_policy
        type(cavg_quality_model),          intent(inout) :: model
        type(cavg_quality_model_spec) :: spec
        integer :: ilinear, iterm, ifeat
        spec = abinitio_learn_base_spec()
        spec%name               = 'learned_pairwise_logistic_v1'
        spec%feature_policy     = trim(feature_policy)
        spec%model_family       = CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC
        spec%weights            = 0.0
        spec%intercept          = real(solution(1))
        spec%linear_coefficients = 0.0
        do ilinear = 1, problem%n_linear
            ifeat = problem%linear_features(ilinear)
            spec%weights(ifeat) = 1.0 / real(problem%n_linear)
            spec%linear_coefficients(ifeat) = real(solution(1 + ilinear))
        end do
        spec%n_interactions = problem%n_interactions
        spec%interaction_terms = 0
        spec%interaction_coefficients = 0.0
        do iterm = 1, problem%n_interactions
            spec%interaction_terms(iterm,:) = problem%interaction_terms(iterm,:)
            spec%interaction_coefficients(iterm) = real(solution(1 + problem%n_linear + iterm))
        end do
        spec%prob_threshold          = 0.5
        spec%regularization_lambda   = real(problem%lambda)
        spec%calibration_temperature = 1.0
        call model%init_spec(spec)
    end subroutine logistic_solution_to_model

    function logistic_cost_wrapper( fun_self, vec, D ) result( cost )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8)                :: cost
        real(kind=8) :: grad(D)
        select type(problem => fun_self)
        type is(cavg_quality_logistic_problem)
            call logistic_problem_fdf(problem, vec, cost, grad, D)
        class default
            THROW_HARD('logistic_cost_wrapper: invalid optimization context')
        end select
    end function logistic_cost_wrapper

    subroutine logistic_gradient_wrapper( fun_self, vec, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: grad(D)
        real(kind=8) :: cost
        select type(problem => fun_self)
        type is(cavg_quality_logistic_problem)
            call logistic_problem_fdf(problem, vec, cost, grad, D)
        class default
            THROW_HARD('logistic_gradient_wrapper: invalid optimization context')
        end select
    end subroutine logistic_gradient_wrapper

    subroutine logistic_fdf_wrapper( fun_self, vec, f, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        select type(problem => fun_self)
        type is(cavg_quality_logistic_problem)
            call logistic_problem_fdf(problem, vec, f, grad, D)
        class default
            THROW_HARD('logistic_fdf_wrapper: invalid optimization context')
        end select
    end subroutine logistic_fdf_wrapper

    subroutine logistic_problem_fdf( problem, vec, f, grad, D )
        type(cavg_quality_logistic_problem), intent(in) :: problem
        integer,                             intent(in) :: D
        real(kind=8),                        intent(in) :: vec(D)
        real(kind=8),                        intent(out):: f, grad(D)
        integer :: iobs, idim
        real(kind=8) :: eta, prob, resid, loss
        if( D /= problem%ndim ) THROW_HARD('logistic_problem_fdf: dimension mismatch')
        if( problem%total_weight <= 0.0d0 ) THROW_HARD('logistic_problem_fdf: nonpositive total weight')
        f    = 0.0d0
        grad = 0.0d0
        do iobs = 1, problem%nobs
            eta  = dot_product(vec, problem%design(iobs,:))
            prob = stable_sigmoid_8(eta)
            if( eta >= 0.0d0 )then
                loss = (1.0d0 - problem%labels(iobs)) * eta + log(1.0d0 + exp(-min(eta, 80.0d0)))
            else
                loss = -problem%labels(iobs) * eta + log(1.0d0 + exp(max(eta, -80.0d0)))
            endif
            f = f + problem%weights(iobs) * loss
            resid = problem%weights(iobs) * (prob - problem%labels(iobs))
            do idim = 1, D
                grad(idim) = grad(idim) + resid * problem%design(iobs,idim)
            end do
        end do
        f    = f    / problem%total_weight
        grad = grad / problem%total_weight
        do idim = 2, D
            f = f + 0.5d0 * problem%lambda * vec(idim)**2
            grad(idim) = grad(idim) + problem%lambda * vec(idim)
        end do
    end subroutine logistic_problem_fdf

    real(kind=8) function stable_sigmoid_8( eta )
        real(kind=8), intent(in) :: eta
        real(kind=8) :: z
        if( eta >= 0.0d0 )then
            z = exp(-min(eta, 80.0d0))
            stable_sigmoid_8 = 1.0d0 / (1.0d0 + z)
        else
            z = exp(max(eta, -80.0d0))
            stable_sigmoid_8 = z / (1.0d0 + z)
        endif
    end function stable_sigmoid_8

    subroutine kill_logistic_problem( self )
        class(cavg_quality_logistic_problem), intent(inout) :: self
        if( allocated(self%design)  ) deallocate(self%design)
        if( allocated(self%labels)  ) deallocate(self%labels)
        if( allocated(self%weights) ) deallocate(self%weights)
        self%nobs              = 0
        self%ndim              = 0
        self%n_linear          = 0
        self%n_interactions    = 0
        self%linear_features   = 0
        self%interaction_terms = 0
        self%lambda            = 0.0d0
        self%total_weight      = 0.0d0
    end subroutine kill_logistic_problem

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

    function abinitio_learn_base_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        type(cavg_quality_model)      :: defaults
        spec = defaults%get_spec()
        spec%name                    = 'abinitio_learn_base'
        spec%feature_policy          = 'microchunk_plus_score_signal'
        spec%weights                 = 0.0
        spec%boundary_margin         = 0.0
        spec%min_score_separation    = minval(LEARN_MINSEPS)
        spec%otsu_min_offset         = 0.0
        spec%otsu_max_offset         = 0.0
        spec%min_accept_frac         = 0.0
        spec%use_lowsep_otsu         = .false.
        spec%use_otsu_window         = .false.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function abinitio_learn_base_spec

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

    subroutine evaluate_policy_grid( dsets, caches, candidate_spec, n_grid, best_spec, best_score, best_tie_specs, &
                                     n_best_ties, top_specs, top_scores, n_top )
        type(cavg_quality_training_dataset), intent(in)    :: dsets(:)
        type(cavg_quality_classify_cache),   intent(in)    :: caches(:)
        type(cavg_quality_model_spec),       intent(in)    :: candidate_spec
        integer,                             intent(inout) :: n_grid, n_best_ties, n_top
        type(cavg_quality_model_spec),       intent(inout) :: best_spec, best_tie_specs(:), top_specs(:)
        real,                                intent(inout) :: best_score, top_scores(:)
        type(cavg_quality_model_spec) :: scan_spec
        integer :: irescue, iaccept
        do irescue = 1, size(LEARN_OTSU_FLAGS)
            scan_spec = candidate_spec
            scan_spec%use_cluster_rescue = LEARN_OTSU_FLAGS(irescue)
            do iaccept = 1, size(LEARN_OTSU_FLAGS)
                scan_spec%enforce_min_accept_frac = LEARN_OTSU_FLAGS(iaccept)
                call evaluate_candidate_spec(dsets, caches, scan_spec, n_grid, best_spec, best_score, &
                    best_tie_specs, n_best_ties, top_specs, top_scores, n_top)
            end do
        end do
    end subroutine evaluate_policy_grid

    subroutine evaluate_candidate_spec( dsets, caches, candidate_spec, n_grid, best_spec, best_score, &
                                        best_tie_specs, n_best_ties, top_specs, top_scores, n_top )
        type(cavg_quality_training_dataset), intent(in)    :: dsets(:)
        type(cavg_quality_classify_cache),   intent(in)    :: caches(:)
        type(cavg_quality_model_spec),       intent(in)    :: candidate_spec
        integer,                             intent(inout) :: n_grid, n_best_ties, n_top
        type(cavg_quality_model_spec),       intent(inout) :: best_spec, best_tie_specs(:), top_specs(:)
        real,                                intent(inout) :: best_score, top_scores(:)
        type(cavg_quality_model)      :: candidate
        type(cavg_quality_model_spec) :: scan_spec
        integer :: ifrac
        real    :: score
        if( candidate_spec%enforce_min_accept_frac )then
            do ifrac = 1, size(LEARN_MIN_ACCEPT_FRACS)
                scan_spec = candidate_spec
                scan_spec%min_accept_frac = LEARN_MIN_ACCEPT_FRACS(ifrac)
                call candidate%init_spec(scan_spec)
                score = macro_balacc_for_model(dsets, candidate, caches)
                n_grid = n_grid + 1
                call consider_model_candidate(candidate%get_spec(), score, best_spec, best_score, &
                    best_tie_specs, n_best_ties, top_specs, top_scores, n_top)
            end do
        else
            scan_spec = candidate_spec
            scan_spec%min_accept_frac = 0.0
            call candidate%init_spec(scan_spec)
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

    integer function n_feature_policies()
        n_feature_policies = CAVG_QUALITY_LEARN_N_STANDARD_POLICIES
    end function n_feature_policies

    subroutine init_learn_feature_mask( feature_mask, trust_resolution )
        logical,          intent(out) :: feature_mask(CAVG_QUALITY_NFEATS)
        logical, optional,intent(in)  :: trust_resolution
        feature_mask = .true.
        if( present(trust_resolution) )then
            if( .not. trust_resolution ) feature_mask(I_NEG_LOG_RES) = .false.
        endif
        if( count(feature_mask) < 1 ) THROW_HARD('init_learn_feature_mask: empty training feature mask')
    end subroutine init_learn_feature_mask

    function feature_policy_name( ipolicy ) result( name )
        integer, intent(in) :: ipolicy
        character(len=64) :: name
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

    subroutine feature_policy_mask( ipolicy, mask, feature_mask )
        integer, intent(in)  :: ipolicy
        logical, intent(out) :: mask(CAVG_QUALITY_NFEATS)
        logical, optional, intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        mask = .false.
        call append_feature_family(mask, 'microchunk')
        call append_feature_family(mask, 'overfit')
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
        if( present(feature_mask) ) mask = mask .and. feature_mask
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

    subroutine apply_feature_policy( ipolicy, weights, feature_mask )
        integer, intent(in)    :: ipolicy
        real,    intent(inout) :: weights(CAVG_QUALITY_NFEATS)
        logical, optional, intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        logical :: mask(CAVG_QUALITY_NFEATS)
        integer :: n_active
        call feature_policy_mask(ipolicy, mask, feature_mask)
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

    subroutine feature_policy_indices( ipolicy, inds, ninds, feature_mask )
        integer, intent(in)  :: ipolicy
        integer, intent(out) :: inds(CAVG_QUALITY_NFEATS)
        integer, intent(out) :: ninds
        logical, optional, intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        logical :: mask(CAVG_QUALITY_NFEATS)
        integer :: ifeat
        call feature_policy_mask(ipolicy, mask, feature_mask)
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
                if( feat_col(ifeat) == 0 ) cycle
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
                if( trim(cavg_quality_feature_family(ifeat)) == 'overfit' ) cycle
                errmsg = 'read_quality_training_dataset: missing z_'//trim(cavg_quality_feature_name(ifeat))//&
                    ' column in '//trim(fname)
                THROW_HARD(trim(errmsg))
            endif
        end do
    end subroutine require_analysis_columns

    subroutine calc_suggested_training_weights( dsets, weights, feature_mask )
        type(cavg_quality_training_dataset), intent(in)  :: dsets(:)
        real,                                intent(out) :: weights(CAVG_QUALITY_NFEATS)
        logical,                             intent(in)  :: feature_mask(CAVG_QUALITY_NFEATS)
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
            if( .not. feature_mask(j) ) cycle
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
        where( .not. feature_mask ) weights = 0.0
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
        type(cavg_quality_result) :: quality
        integer :: ids, tp, fp, tn, fn, nused, role
        integer :: tail_n
        real :: balacc, mean_score, tail_score, role_scores(size(dsets))
        macro_balacc_for_model = 0.0
        if( present(caches) )then
            if( size(caches) /= size(dsets) ) THROW_HARD('macro_balacc_for_model: cache/dataset size mismatch')
        endif
        nused = 0
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            if( role == LEARN_ROLE_SKIP ) cycle
            if( present(caches) )then
                call classify_training_dataset_cached_detail(dsets(ids), caches(ids), model, quality, tp, fp, tn, fn)
            else
                call classify_training_dataset_detail(dsets(ids), model, quality, tp, fp, tn, fn)
            endif
            balacc = learn_balacc_from_confusion(tp, fp, tn, fn, role)
            balacc = balacc - overfit_false_positive_penalty(dsets(ids), quality, role)
            nused = nused + 1
            role_scores(nused) = balacc
            call quality%kill()
        end do
        if( nused == 0 ) THROW_HARD('macro_balacc_for_model: no scoreable training datasets')
        call sort_real_prefix_ascending(role_scores, nused)
        mean_score = sum(role_scores(1:nused)) / real(nused)
        tail_n     = max(1, min(nused, ceiling(LEARN_ROBUST_TAIL_FRAC * real(nused))))
        tail_score = sum(role_scores(1:tail_n)) / real(tail_n)
        macro_balacc_for_model = (1.0 - LEARN_ROBUST_TAIL_WEIGHT) * mean_score + &
            LEARN_ROBUST_TAIL_WEIGHT * tail_score
    end function macro_balacc_for_model

    real function learn_balacc_from_confusion( tp, fp, tn, fn, role )
        integer, intent(in) :: tp, fp, tn, fn, role
        real :: precision, recall, specificity, f1, balacc, accuracy
        call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
        select case(role)
            case(LEARN_ROLE_BALANCED)
                learn_balacc_from_confusion = fn_tolerant_specificity_score(specificity, tp, fn)
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

    real function fn_tolerant_specificity_score( specificity, tp, fn )
        real,    intent(in) :: specificity
        integer, intent(in) :: tp, fn
        fn_tolerant_specificity_score = specificity - LEARN_BALANCED_FN_RATE_PENALTY * &
            balanced_fn_excess_rate(tp, fn)
    end function fn_tolerant_specificity_score

    real function balanced_fn_excess_rate( tp, fn )
        integer, intent(in) :: tp, fn
        integer :: npos
        npos = tp + fn
        if( npos <= 0 )then
            balanced_fn_excess_rate = 0.0
        else
            balanced_fn_excess_rate = max(0.0, real(fn) / real(npos) - LEARN_BALANCED_FN_TOLERANCE_FRAC)
        endif
    end function balanced_fn_excess_rate

    real function overfit_false_positive_penalty( dset, quality, role )
        type(cavg_quality_training_dataset), intent(in) :: dset
        type(cavg_quality_result),           intent(in) :: quality
        integer,                             intent(in) :: role
        integer :: nbad_overfit, nfp_overfit
        real :: fp_rate, excess_rate
        overfit_false_positive_penalty = 0.0
        if( role /= LEARN_ROLE_BALANCED .and. role /= LEARN_ROLE_SPECIFICITY_ONLY ) return
        if( .not. dataset_is_overfit_focus(dset) ) return
        nbad_overfit = count_trainable_bad_overfit_signature(dset)
        if( nbad_overfit <= 0 ) return
        nfp_overfit = count_accepted_bad_overfit_signature(dset, quality)
        fp_rate = real(nfp_overfit) / real(nbad_overfit)
        excess_rate = max(0.0, fp_rate - LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET)
        overfit_false_positive_penalty = LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY * excess_rate * excess_rate
    end function overfit_false_positive_penalty

    logical function dataset_is_overfit_focus( dset )
        type(cavg_quality_training_dataset), intent(in) :: dset
        integer :: nfit, nbad_overfit
        nfit = count_trainable_classes(dset)
        if( nfit <= 0 )then
            dataset_is_overfit_focus = .false.
            return
        endif
        nbad_overfit = count_trainable_bad_overfit_signature(dset)
        dataset_is_overfit_focus = real(nbad_overfit) / real(nfit) >= LEARN_OVERFIT_FP_FOCUS_MIN_TRAINABLE_FRAC
    end function dataset_is_overfit_focus

    integer function count_trainable_bad_overfit_signature( dset )
        type(cavg_quality_training_dataset), intent(in) :: dset
        integer :: irow
        count_trainable_bad_overfit_signature = 0
        do irow = 1, dset%ncls
            if( dset%hard_reject(irow) ) cycle
            if( dset%manual_states(irow) > 0 ) cycle
            if( bad_overfit_signature(dset%features(irow,:)) ) &
                count_trainable_bad_overfit_signature = count_trainable_bad_overfit_signature + 1
        end do
    end function count_trainable_bad_overfit_signature

    integer function count_accepted_bad_overfit_signature( dset, quality )
        type(cavg_quality_training_dataset), intent(in) :: dset
        type(cavg_quality_result),           intent(in) :: quality
        integer :: irow
        if( .not. allocated(quality%states) ) &
            THROW_HARD('count_accepted_bad_overfit_signature: missing quality states')
        if( size(quality%states) /= dset%ncls ) &
            THROW_HARD('count_accepted_bad_overfit_signature: state/dataset size mismatch')
        count_accepted_bad_overfit_signature = 0
        do irow = 1, dset%ncls
            if( dset%hard_reject(irow) ) cycle
            if( dset%manual_states(irow) > 0 ) cycle
            if( quality%states(irow) <= 0 ) cycle
            if( bad_overfit_signature(dset%features(irow,:)) ) &
                count_accepted_bad_overfit_signature = count_accepted_bad_overfit_signature + 1
        end do
    end function count_accepted_bad_overfit_signature

    logical function bad_overfit_signature( z_features )
        real, intent(in) :: z_features(:)
        logical :: low_local_variance, poor_bandpass_localization, poor_support
        if( size(z_features) /= CAVG_QUALITY_NFEATS ) &
            THROW_HARD('bad_overfit_signature: invalid feature count')
        poor_support = z_features(I_CC_AREA_FRAC) < LEARN_BAD_OVERFIT_SUPPORT_MAX
        low_local_variance = z_features(I_NEG_LOCVAR_FG) > LEARN_BAD_OVERFIT_LOWVAR_FG_MIN
        poor_bandpass_localization = z_features(I_BP40_100_CENTER_EDGE_VAR) < LEARN_OVERFIT_FP_BP_Z_MAX
        bad_overfit_signature = poor_support .and. (low_local_variance .or. poor_bandpass_localization)
    end function bad_overfit_signature

    subroutine sort_real_prefix_ascending( vals, nvals )
        real,    intent(inout) :: vals(:)
        integer, intent(in)    :: nvals
        integer :: i, j
        real :: tmp
        if( nvals > size(vals) ) THROW_HARD('sort_real_prefix_ascending: invalid prefix size')
        do i = 1, nvals - 1
            do j = i + 1, nvals
                if( vals(j) < vals(i) )then
                    tmp = vals(i)
                    vals(i) = vals(j)
                    vals(j) = tmp
                endif
            end do
        end do
    end subroutine sort_real_prefix_ascending

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

    subroutine select_preferred_best_tie( base_spec, best_tie_specs, n_best_ties, best_spec )
        type(cavg_quality_model_spec), intent(in)    :: base_spec, best_tie_specs(:)
        integer,                       intent(in)    :: n_best_ties
        type(cavg_quality_model_spec), intent(inout) :: best_spec
        integer :: i, nstored
        real    :: dist, best_dist
        nstored = min(n_best_ties, size(best_tie_specs))
        if( nstored <= 1 ) return
        best_dist = huge(1.0)
        do i = 1, nstored
            dist = threshold_tie_distance(best_tie_specs(i), base_spec)
            if( dist < best_dist - EPS )then
                best_dist = dist
                best_spec = best_tie_specs(i)
            endif
        end do
    end subroutine select_preferred_best_tie

    real function threshold_tie_distance( spec, base_spec )
        type(cavg_quality_model_spec), intent(in) :: spec, base_spec
        threshold_tie_distance = scaled_absdiff(spec%min_score_separation, base_spec%min_score_separation, &
            minval(LEARN_MINSEPS), maxval(LEARN_MINSEPS))
        threshold_tie_distance = threshold_tie_distance + scaled_absdiff(spec%boundary_margin, &
            base_spec%boundary_margin, minval(LEARN_BOUNDARY_MARGINS), maxval(LEARN_BOUNDARY_MARGINS))
        if( spec%enforce_min_accept_frac ) threshold_tie_distance = threshold_tie_distance + scaled_absdiff( &
            spec%min_accept_frac, base_spec%min_accept_frac, minval(LEARN_MIN_ACCEPT_FRACS), &
            maxval(LEARN_MIN_ACCEPT_FRACS))
        ! Boolean knobs count as a quarter of one normalized scalar axis:
        ! enough to prefer the neutral foundation among exact-score ties without
        ! overwhelming a nearby threshold or margin value.
        if( spec%use_lowsep_otsu .neqv. base_spec%use_lowsep_otsu ) threshold_tie_distance = threshold_tie_distance + 0.25
        if( spec%use_otsu_window .neqv. base_spec%use_otsu_window ) threshold_tie_distance = threshold_tie_distance + 0.25
        if( spec%use_cluster_rescue .neqv. base_spec%use_cluster_rescue ) &
            threshold_tie_distance = threshold_tie_distance + 0.25
        if( spec%enforce_min_accept_frac .neqv. base_spec%enforce_min_accept_frac ) &
            threshold_tie_distance = threshold_tie_distance + 0.25
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
        if( trim(model%model_family) == CAVG_MODEL_FAMILY_LINEAR )then
            call cached_decision_confusion(cache, model, dset%manual_states, tp, fp, tn, fn)
        else
            call classify_training_dataset(dset, model, tp, fp, tn, fn)
        endif
    end subroutine classify_training_dataset_cached

    subroutine classify_training_dataset_cached_detail( dset, cache, model, quality, tp, fp, tn, fn )
        type(cavg_quality_training_dataset), intent(in)    :: dset
        type(cavg_quality_classify_cache),   intent(in)    :: cache
        type(cavg_quality_model),            intent(in)    :: model
        type(cavg_quality_result),           intent(inout) :: quality
        integer,                             intent(out)   :: tp, fp, tn, fn
        logical, allocatable :: pred(:), ref(:)
        integer :: nfit
        if( trim(model%model_family) /= CAVG_MODEL_FAMILY_LINEAR )then
            call classify_training_dataset_detail(dset, model, quality, tp, fp, tn, fn)
            return
        endif
        call quality%kill()
        quality%hard_reject = dset%hard_reject
        call apply_cached_decision_to_quality(cache, model, quality)
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
    end subroutine classify_training_dataset_cached_detail

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
                                                n_grid, top_specs, top_scores, n_top, best_tie_specs, n_best_ties, &
                                                feature_mask )
        character(len=*),                    intent(in) :: fname
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model_spec),       intent(in) :: base_spec
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: suggested_weights(:)
        real,                                intent(in) :: best_score
        integer,                             intent(in) :: n_grid, n_top, n_best_ties
        type(cavg_quality_model_spec),       intent(in) :: top_specs(:), best_tie_specs(:)
        real,                                intent(in) :: top_scores(:)
        logical,                             intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        type(cavg_quality_learn_diagnostics) :: diag
        integer :: funit, i
        call collect_learn_diagnostics(dsets, learned_model, diag)
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# model_cavgs_rejection learn report'
        write(funit,'(A,A)') 'foundation_model=', trim(base_spec%name)
        write(funit,'(A,A)') 'learned_model=', trim(learned_model%name)
        write(funit,'(A,F10.5)') 'macro_learn_score=', best_score
        write(funit,'(A,I0)') 'n_datasets=', size(dsets)
        write(funit,'(A,I0)') 'model_search_grid_n=', n_grid
        write(funit,'(A,I0)') 'best_tie_count=', n_best_ties
        write(funit,'(A,I0)') 'top_candidates_reported=', n_top
        write(funit,'(A)') 'note=scalar_feature_space_only'
        write(funit,'(A)') 'note=feature_policy_scans_include_microchunk_and_overfit_features'
        write(funit,'(A)') 'note=hard_rejected_rows_are_reported_but_excluded_from_model_fit_and_scoring'
        write(funit,'(A)') 'note=feature_weights_use_only_datasets_with_both_manual_states_after_hard_rejects'
        write(funit,'(A)') 'note=macro_learn_score_is_mean_plus_lower_tail_robust_score'
        write(funit,'(A)') 'note=balanced_datasets_are_scored_by_specificity_with_small_fn_rate_tolerance'
        write(funit,'(A)') 'note=overfit_focus_datasets_penalize_excess_accepted_bad_fuzzy_ball_signature_rate'
        write(funit,'(A)') 'note=trainable_good_only_datasets_are_scored_by_guarded_recall'
        write(funit,'(A)') 'note=trainable_bad_only_datasets_are_scored_by_specificity_unless_good_classes_were_hard_rejected'
        write(funit,'(A)') 'note=feature_weights_derived_from_training_data_no_base_weight_blending'
        write(funit,'(A)') 'note=learn_mode_uses_neutral_abinitio_foundation_not_quality_model_or_infile_seed'
        call write_feature_policy_grid(funit)
        call write_learn_feature_mask(funit, feature_mask)
        write(funit,'(A,ES14.6)') 'robust_score_tail_frac=', LEARN_ROBUST_TAIL_FRAC
        write(funit,'(A,ES14.6)') 'robust_score_tail_weight=', LEARN_ROBUST_TAIL_WEIGHT
        write(funit,'(A,ES14.6)') 'balanced_score_fn_tolerance_frac=', LEARN_BALANCED_FN_TOLERANCE_FRAC
        write(funit,'(A,ES14.6)') 'balanced_score_fn_rate_penalty=', LEARN_BALANCED_FN_RATE_PENALTY
        write(funit,'(A,ES14.6)') 'overfit_focus_fp_accept_rate_target=', &
            LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET
        write(funit,'(A,ES14.6)') 'overfit_focus_fp_excess_rate_penalty=', &
            LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY
        write(funit,'(A,ES14.6)') 'overfit_focus_min_trainable_signature_frac=', &
            LEARN_OVERFIT_FP_FOCUS_MIN_TRAINABLE_FRAC
        write(funit,'(A,ES14.6)') 'overfit_focus_bp_z_max=', LEARN_OVERFIT_FP_BP_Z_MAX
        write(funit,'(A,ES14.6)') 'grid_recall_only_floor=', LEARN_RECALL_ONLY_FLOOR
        write(funit,'(A,ES14.6)') 'grid_recall_only_shortfall_penalty=', LEARN_RECALL_ONLY_PENALTY
        call write_real_list(funit, 'grid_min_score_separations=', LEARN_MINSEPS)
        call write_real_list(funit, 'grid_boundary_margins=', LEARN_BOUNDARY_MARGINS)
        call write_logical_list(funit, 'grid_use_lowsep_otsu=', LEARN_OTSU_FLAGS)
        call write_logical_list(funit, 'grid_use_otsu_window=', LEARN_OTSU_FLAGS)
        call write_real_list(funit, 'grid_otsu_min_offsets=', LEARN_OTSU_MIN_OFFSETS)
        call write_real_list(funit, 'grid_otsu_max_offsets=', LEARN_OTSU_MAX_OFFSETS)
        call write_logical_list(funit, 'grid_use_cluster_rescue=', LEARN_OTSU_FLAGS)
        call write_logical_list(funit, 'grid_enforce_min_accept_frac=', LEARN_OTSU_FLAGS)
        call write_real_list(funit, 'grid_min_accept_fracs=', LEARN_MIN_ACCEPT_FRACS)
        write(funit,'(A)', advance='no') 'suggested_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') suggested_weights(i)
        end do
        write(funit,*)
        write(funit,'(A)') ''
        call write_learn_search_diagnostics(funit, learned_model, diag, best_tie_specs, n_best_ties)
        call write_otsu_ablation_diagnostics(funit, dsets, learned_model)
        call write_feature_screen_diagnostics(funit, dsets, base_spec, suggested_weights, learned_model, best_score, &
            feature_mask)
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

    subroutine write_cavg_quality_logistic_learn_report( fname, dsets, learned_model, best_score, &
                                                         n_candidates, best_objective, feature_mask )
        character(len=*),                    intent(in) :: fname
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: best_score, best_objective
        integer,                             intent(in) :: n_candidates
        logical,                             intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        type(cavg_quality_learn_diagnostics) :: diag
        integer :: funit
        call collect_learn_diagnostics(dsets, learned_model, diag)
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# model_cavgs_rejection pairwise-logistic learn report'
        write(funit,'(A,A)') 'learned_model=', trim(learned_model%name)
        write(funit,'(A,A)') 'model_family=', trim(learned_model%model_family)
        write(funit,'(A,F10.5)') 'macro_learn_score=', best_score
        write(funit,'(A,ES14.6)') 'best_weighted_logistic_objective=', best_objective
        write(funit,'(A,I0)') 'n_datasets=', size(dsets)
        write(funit,'(A,I0)') 'model_search_grid_n=', n_candidates
        write(funit,'(A)') 'note=pairwise_logistic_uses_existing_normalized_training_table_features'
        write(funit,'(A)') 'note=hard_rejected_rows_are_reported_but_excluded_from_model_fit_and_scoring'
        write(funit,'(A)') 'note=macro_learn_score_is_mean_plus_lower_tail_robust_score'
        write(funit,'(A)') 'note=balanced_datasets_are_scored_by_specificity_with_small_fn_rate_tolerance'
        write(funit,'(A)') 'note=overfit_focus_datasets_penalize_excess_accepted_bad_fuzzy_ball_signature_rate'
        write(funit,'(A)') 'note=balanced_dataset_logistic_loss_uses_moderate_good_class_weight'
        write(funit,'(A)') 'note=manually_bad_overfit_signature_examples_get_extra_logistic_loss_weight'
        write(funit,'(A)') 'note=overfit_focus_signature_fraction_loss_scale_is_reported_but_currently_neutral'
        write(funit,'(A)') 'note=trainable_good_only_datasets_contribute_recall_evidence'
        write(funit,'(A)') 'note=trainable_bad_only_datasets_contribute_specificity_evidence'
        call write_feature_policy_grid(funit)
        call write_learn_feature_mask(funit, feature_mask)
        write(funit,'(A,ES14.6)') 'robust_score_tail_frac=', LEARN_ROBUST_TAIL_FRAC
        write(funit,'(A,ES14.6)') 'robust_score_tail_weight=', LEARN_ROBUST_TAIL_WEIGHT
        write(funit,'(A,ES14.6)') 'balanced_score_fn_tolerance_frac=', LEARN_BALANCED_FN_TOLERANCE_FRAC
        write(funit,'(A,ES14.6)') 'balanced_score_fn_rate_penalty=', LEARN_BALANCED_FN_RATE_PENALTY
        write(funit,'(A,ES14.6)') 'balanced_loss_good_class_weight=', LEARN_BALANCED_GOOD_LOSS_WEIGHT
        write(funit,'(A,ES14.6)') 'balanced_loss_bad_class_weight=', LEARN_BALANCED_BAD_LOSS_WEIGHT
        write(funit,'(A,ES14.6)') 'bad_overfit_loss_multiplier=', LEARN_BAD_OVERFIT_LOSS_MULT
        write(funit,'(A,ES14.6)') 'bad_overfit_focus_loss_scale=', LEARN_OVERFIT_FOCUS_BAD_LOSS_SCALE
        write(funit,'(A,ES14.6)') 'bad_overfit_lowvar_fg_min=', LEARN_BAD_OVERFIT_LOWVAR_FG_MIN
        write(funit,'(A,ES14.6)') 'bad_overfit_support_max=', LEARN_BAD_OVERFIT_SUPPORT_MAX
        write(funit,'(A,ES14.6)') 'overfit_focus_fp_accept_rate_target=', &
            LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET
        write(funit,'(A,ES14.6)') 'overfit_focus_fp_excess_rate_penalty=', &
            LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY
        write(funit,'(A,ES14.6)') 'overfit_focus_min_trainable_signature_frac=', &
            LEARN_OVERFIT_FP_FOCUS_MIN_TRAINABLE_FRAC
        write(funit,'(A,ES14.6)') 'overfit_focus_bp_z_max=', LEARN_OVERFIT_FP_BP_Z_MAX
        call write_real_list(funit, 'grid_logistic_lambdas=', LOGISTIC_LAMBDAS)
        call write_real_list(funit, 'grid_probability_thresholds=', LOGISTIC_THRESHOLDS)
        call write_fixed_model_summary(funit, learned_model)
        call write_logistic_coefficient_table(funit, learned_model)
        write(funit,'(A)') ''
        call write_evaluate_diagnostics(funit, learned_model, diag)
        call write_dataset_metric_table(funit, dsets, learned_model, 'learn_score')
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(fname)
    end subroutine write_cavg_quality_logistic_learn_report

    subroutine write_learn_feature_mask( funit, feature_mask )
        integer, intent(in) :: funit
        logical, intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        integer :: ifeat
        write(funit,'(A,A)') 'trust_resolution=', yes_no_string(feature_mask(I_NEG_LOG_RES))
        if( count(.not. feature_mask) == 0 )then
            write(funit,'(A)') 'inactive_training_features=none'
            return
        endif
        write(funit,'(A)', advance='no') 'inactive_training_features='
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( feature_mask(ifeat) ) cycle
            if( ifeat > 1 )then
                if( count(.not. feature_mask(1:ifeat-1)) > 0 ) write(funit,'(A)', advance='no') ','
            endif
            write(funit,'(A)', advance='no') trim(cavg_quality_feature_name(ifeat))
        end do
        write(funit,*)
    end subroutine write_learn_feature_mask

    character(len=3) function yes_no_string( val )
        logical, intent(in) :: val
        if( val )then
            yes_no_string = 'yes'
        else
            yes_no_string = 'no'
        endif
    end function yes_no_string

    subroutine write_logistic_coefficient_table( funit, model )
        integer,                  intent(in) :: funit
        type(cavg_quality_model), intent(in) :: model
        integer :: ifeat, iterm
        write(funit,'(A)') ''
        write(funit,'(A)') 'logistic_coefficient_header=term,feature_a,feature_b,coefficient'
        write(funit,'(A,ES14.6)') 'logistic_coefficient,intercept,,,', model%intercept
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( abs(model%linear_coefficients(ifeat)) <= EPS ) cycle
            write(funit,'(A,A,A,ES14.6)') 'logistic_coefficient,linear,', &
                trim(cavg_quality_feature_name(ifeat)), ',,', model%linear_coefficients(ifeat)
        end do
        do iterm = 1, model%n_interactions
            if( abs(model%interaction_coefficients(iterm)) <= EPS ) cycle
            write(funit,'(A,A,A,A,A,ES14.6)') 'logistic_coefficient,pairwise,', &
                trim(cavg_quality_feature_name(model%interaction_terms(iterm,1))), ',', &
                trim(cavg_quality_feature_name(model%interaction_terms(iterm,2))), ',', &
                model%interaction_coefficients(iterm)
        end do
    end subroutine write_logistic_coefficient_table

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
        write(funit,'(A,A)') 'model=', trim(model%name)
        write(funit,'(A,F10.5)') 'macro_evaluate_score=', eval_score
        write(funit,'(A,I0)') 'n_datasets=', size(dsets)
        write(funit,'(A)') 'note=fixed_model_no_refit'
        write(funit,'(A)') 'note=analysis_table_rows_reclassified_with_selected_model'
        write(funit,'(A)') 'note=hard_rejected_rows_are_reported_but_excluded_from_model_scoring'
        write(funit,'(A)') 'note=macro_learn_score_is_mean_plus_lower_tail_robust_score'
        write(funit,'(A)') 'note=balanced_datasets_are_scored_by_specificity_with_small_fn_rate_tolerance'
        write(funit,'(A)') 'note=overfit_focus_datasets_penalize_excess_accepted_bad_fuzzy_ball_signature_rate'
        write(funit,'(A)') 'note=trainable_good_only_datasets_are_scored_by_guarded_recall'
        write(funit,'(A)') 'note=trainable_bad_only_datasets_are_scored_by_specificity_unless_good_classes_were_hard_rejected'
        write(funit,'(A,ES14.6)') 'robust_score_tail_frac=', LEARN_ROBUST_TAIL_FRAC
        write(funit,'(A,ES14.6)') 'robust_score_tail_weight=', LEARN_ROBUST_TAIL_WEIGHT
        write(funit,'(A,ES14.6)') 'balanced_score_fn_tolerance_frac=', LEARN_BALANCED_FN_TOLERANCE_FRAC
        write(funit,'(A,ES14.6)') 'balanced_score_fn_rate_penalty=', LEARN_BALANCED_FN_RATE_PENALTY
        write(funit,'(A,ES14.6)') 'overfit_focus_fp_accept_rate_target=', &
            LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET
        write(funit,'(A,ES14.6)') 'overfit_focus_fp_excess_rate_penalty=', &
            LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY
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
        write(funit,'(A,A)') 'model_family=', trim(model%model_family)
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
        if( trim(model%model_family) == CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC )then
            write(funit,'(A,ES14.6)') 'model_intercept=', model%intercept
            write(funit,'(A)', advance='no') 'model_linear_coefficients='
            do i = 1, CAVG_QUALITY_NFEATS
                if( i > 1 ) write(funit,'(A)', advance='no') ','
                write(funit,'(ES14.6)', advance='no') model%linear_coefficients(i)
            end do
            write(funit,*)
            write(funit,'(A)', advance='no') 'model_interaction_terms='
            do i = 1, model%n_interactions
                if( i > 1 ) write(funit,'(A)', advance='no') ','
                write(funit,'(A,A,A)', advance='no') &
                    trim(cavg_quality_feature_name(model%interaction_terms(i,1))), ':', &
                    trim(cavg_quality_feature_name(model%interaction_terms(i,2)))
            end do
            write(funit,*)
            write(funit,'(A)', advance='no') 'model_interaction_coefficients='
            do i = 1, model%n_interactions
                if( i > 1 ) write(funit,'(A)', advance='no') ','
                write(funit,'(ES14.6)', advance='no') model%interaction_coefficients(i)
            end do
            write(funit,*)
            write(funit,'(A,ES14.6)') 'model_prob_threshold=', model%prob_threshold
            write(funit,'(A,ES14.6)') 'model_regularization_lambda=', model%regularization_lambda
            write(funit,'(A,ES14.6)') 'model_calibration_temperature=', model%calibration_temperature
        endif
    end subroutine write_fixed_model_summary

    subroutine write_dataset_metric_table( funit, dsets, model, score_name )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model),            intent(in) :: model
        character(len=*),                    intent(in) :: score_name
        integer :: ids, tp, fp, tn, fn, role
        type(cavg_quality_result) :: quality
        real :: precision, recall, specificity, f1, role_score, accuracy, overfit_penalty
        write(funit,'(A)') ''
        write(funit,'(A,A,A)') 'dataset,n_classes,n_trainable,trainable_manual_good,trainable_manual_bad,learn_role,', &
            'tp,fp,tn,fn,precision,recall,specificity,f1,', trim(score_name)//&
            ',accuracy,hard_rejected_manual_good,overfit_focus_bad,overfit_focus_fp,overfit_penalty'
        do ids = 1, size(dsets)
            role = dataset_learn_role(dsets(ids))
            call classify_training_dataset_detail(dsets(ids), model, quality, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, role_score, accuracy)
            role_score = learn_balacc_from_confusion(tp, fp, tn, fn, role)
            overfit_penalty = overfit_false_positive_penalty(dsets(ids), quality, role)
            role_score = role_score - overfit_penalty
            write(funit,'(A,A,I0,A,I0,A,I0,A,I0,A,A,A,I0,A,I0,A,I0,A,I0,'//&
                'A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,I0,A,I0,A,I0,A,F10.5)') &
                trim(dsets(ids)%dataset_id), ',', dsets(ids)%ncls, ',', count_trainable_classes(dsets(ids)), ',', &
                count_trainable_manual_good(dsets(ids)), ',', count_trainable_manual_bad(dsets(ids)), ',', &
                trim(dataset_learn_role_name(role)), ',', &
                tp, ',', fp, ',', tn, ',', fn, ',', precision, ',', recall, ',', specificity, ',', f1, ',', &
                role_score, ',', accuracy, ',', count_hard_rejected_manual_good(dsets(ids)), ',', &
                count_trainable_bad_overfit_signature(dsets(ids)), ',', &
                count_accepted_bad_overfit_signature(dsets(ids), quality), ',', overfit_penalty
            call quality%kill()
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

    subroutine write_feature_screen_diagnostics( funit, dsets, base_spec, suggested_weights, learned_model, best_score, &
                                                 feature_mask )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        type(cavg_quality_model_spec),       intent(in) :: base_spec
        type(cavg_quality_model),            intent(in) :: learned_model
        real,                                intent(in) :: suggested_weights(:)
        real,                                intent(in) :: best_score
        logical,                             intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        write(funit,'(A)') ''
        write(funit,'(A)') '# feature-screen diagnostics'
        call write_feature_signal_diagnostics(funit, dsets, base_spec, suggested_weights, learned_model)
        write(funit,'(A)') ''
        call write_feature_drop_diagnostics(funit, dsets, learned_model, best_score)
        write(funit,'(A)') ''
        call write_feature_policy_screen(funit, dsets, feature_mask)
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

    subroutine write_feature_policy_screen( funit, dsets, feature_mask )
        integer,                             intent(in) :: funit
        type(cavg_quality_training_dataset), intent(in) :: dsets(:)
        logical,                             intent(in) :: feature_mask(CAVG_QUALITY_NFEATS)
        integer :: inds(CAVG_QUALITY_NFEATS)
        integer :: ipol, ninds
        write(funit,'(A)') 'feature_policy_lodo_header=feature_policy,n_features,mean_auc,min_auc,min_auc_dataset,'//&
            'mean_oracle_score,min_oracle_score,min_score_dataset,total_tp,total_fp,total_tn,total_fn'
        do ipol = 1, n_feature_policies()
            call feature_policy_indices(ipol, inds, ninds, feature_mask)
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
            if( dataset_is_overfit_focus(dsets(ids)) )then
                diag%n_overfit_focus   = diag%n_overfit_focus + 1
                diag%overfit_focus_bad = diag%overfit_focus_bad + &
                    count_trainable_bad_overfit_signature(dsets(ids))
                diag%overfit_focus_fp  = diag%overfit_focus_fp + &
                    count_accepted_bad_overfit_signature(dsets(ids), quality)
            endif
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
        write(detail,'(A,F6.3,A,F6.3)') 'tail_frac=', LEARN_ROBUST_TAIL_FRAC, ';tail_weight=', &
            LEARN_ROBUST_TAIL_WEIGHT
        call write_search_diagnostic(funit, 'note', 'macro_learn_score', &
            'mean_plus_lower_tail_robust_score', trim(detail))
        write(detail,'(A,F6.3,A,F6.3)') 'fn_tolerance_frac=', LEARN_BALANCED_FN_TOLERANCE_FRAC, &
            ';rate_penalty=', LEARN_BALANCED_FN_RATE_PENALTY
        call write_search_diagnostic(funit, 'note', 'balanced_dataset_score', &
            'specificity_with_small_fn_rate_tolerance', trim(detail))
        write(detail,'(A,I0,A,I0,A,I0,A,F6.3,A,F6.3)') 'datasets=', diag%n_overfit_focus, &
            ';bad_signature=', diag%overfit_focus_bad, ';accepted_bad_signature=', diag%overfit_focus_fp, &
            ';target=', LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET, &
            ';excess_penalty=', LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY
        call write_search_diagnostic(funit, 'note', 'overfit_focus_penalty', &
            'accepted_bad_signature_rate_hinge_quadratic', trim(detail))
        write(detail,'(A,F6.3,A,F6.3)') 'floor=', LEARN_RECALL_ONLY_FLOOR, ';penalty=', &
            LEARN_RECALL_ONLY_PENALTY
        call write_search_diagnostic(funit, 'note', 'trainable_good_only_score', 'guarded_recall', trim(detail))
        call write_search_diagnostic(funit, 'note', 'feature_weights', 'auc_no_base_blending', &
            'one_auc_derived_candidate_per_feature_policy')
        call write_search_diagnostic(funit, 'note', 'feature_policy', trim(learned_model%feature_policy), &
            'selected_family_set_encoded_by_zeroed_model_weights')
        call write_minsep_diagnostic(funit, learned_model%min_score_separation, best_tie_specs, n_best_ties)
        call write_grid_position_diagnostic(funit, 'boundary_margin', learned_model%boundary_margin, &
            LEARN_BOUNDARY_MARGINS, 'best_at_lowest_value_consider_more_negative_if_junk_leaks_after_validation', &
            'best_at_highest_value_consider_more_positive_if_good_classes_are_rejected')
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
        if( learned_model%enforce_min_accept_frac )then
            call write_grid_position_diagnostic(funit, 'min_accept_frac', learned_model%min_accept_frac, &
                LEARN_MIN_ACCEPT_FRACS, 'best_at_lowest_value_consider_lower_if_model_keeps_too_much_junk', &
                'best_at_highest_value_consider_higher_if_model_rejects_good_classes')
        else
            call write_search_diagnostic(funit, 'note', 'min_accept_frac', 'inactive', &
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
        write(detail,'(A,F6.3,A,F6.3)') 'tail_frac=', LEARN_ROBUST_TAIL_FRAC, ';tail_weight=', &
            LEARN_ROBUST_TAIL_WEIGHT
        call write_evaluate_diagnostic(funit, 'note', 'macro_evaluate_score', &
            'mean_plus_lower_tail_robust_score', trim(detail))
        write(detail,'(A,F6.3,A,F6.3)') 'fn_tolerance_frac=', LEARN_BALANCED_FN_TOLERANCE_FRAC, &
            ';rate_penalty=', LEARN_BALANCED_FN_RATE_PENALTY
        call write_evaluate_diagnostic(funit, 'note', 'balanced_dataset_score', &
            'specificity_with_small_fn_rate_tolerance', trim(detail))
        write(detail,'(A,I0,A,I0,A,I0,A,F6.3,A,F6.3)') 'datasets=', diag%n_overfit_focus, &
            ';bad_signature=', diag%overfit_focus_bad, ';accepted_bad_signature=', diag%overfit_focus_fp, &
            ';target=', LEARN_OVERFIT_FP_ACCEPT_RATE_TARGET, &
            ';excess_penalty=', LEARN_OVERFIT_FP_EXCESS_RATE_PENALTY
        call write_evaluate_diagnostic(funit, 'note', 'overfit_focus_penalty', &
            'accepted_bad_signature_rate_hinge_quadratic', trim(detail))
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
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%use_cluster_rescue, ';rescue_like_datasets=', &
            diag%n_rescue_like, ';fp=', diag%rescue_like_fp, ';fn=', diag%rescue_like_fn
        call write_search_diagnostic(funit, rescue_policy_level(model, diag), 'use_cluster_rescue', 'searched', &
            trim(detail))
        write(detail,'(A,F8.4,A,I0,A,I0,A,I0)') 'fixed=', model%cluster_rescue_margin, ';rescue_like_datasets=', &
            diag%n_rescue_like, ';fp=', diag%rescue_like_fp, ';fn=', diag%rescue_like_fn
        call write_search_diagnostic(funit, rescue_policy_level(model, diag), 'cluster_rescue_margin', 'fixed', &
            trim(detail))
        write(detail,'(A,L1,A,I0,A,I0,A,I0)') 'selected=', model%enforce_min_accept_frac, ';min_accept_datasets=', &
            diag%n_min_accept_like, ';fp=', diag%min_accept_like_fp, ';fn=', diag%min_accept_like_fn
        call write_search_diagnostic(funit, min_accept_policy_level(model, diag), 'enforce_min_accept_frac', &
            'searched', trim(detail))
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

    subroutine write_feature_policy_grid( funit )
        integer, intent(in) :: funit
        integer :: ipol
        write(funit,'(A)', advance='no') 'grid_feature_policies='
        do ipol = 1, n_feature_policies()
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
