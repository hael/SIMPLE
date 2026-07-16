!@descr: shared types for class-average quality analysis
module simple_cavg_quality_types
use simple_defs,   only: LONGSTRLEN
use simple_string, only: string
implicit none
private

public :: CAVG_QUALITY_NFEATS
public :: CAVG_QUALITY_MAX_INTERACTIONS
public :: CAVG_MODEL_FAMILY_LINEAR
public :: CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC
public :: CAVG_RELATIONAL_SCHEMA_NONE
public :: CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1
public :: CAVG_RELATIONAL_DEFAULT_KNN
public :: CAVG_RELATIONAL_DEFAULT_CORR_HP
public :: CAVG_RELATIONAL_DEFAULT_CORR_LP
public :: CAVG_RELATIONAL_DEFAULT_CORR_TRS
public :: CAVG_QUALITY_CONTEXT_CHUNK
public :: CAVG_QUALITY_CONTEXT_SIEVE
public :: CAVG_QUALITY_CONTEXT_POOL
public :: CAVG_REJECT_REASON_NONE
public :: CAVG_REJECT_REASON_POP
public :: CAVG_REJECT_REASON_BAD_PIXELS
public :: CAVG_REJECT_REASON_NO_COMPONENT
public :: CAVG_REJECT_REASON_MASK_GEOMETRY
public :: CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW
public :: CAVG_REJECT_REASON_LOCVAR_FG_LOW
public :: CAVG_REJECT_REASON_FUZZY_BALL_SIGNAL_LOW
public :: EPS
public :: CLIP_Z
public :: cavg_quality_feature_def
public :: cavg_quality_model_spec
public :: cavg_quality_result
public :: cavg_quality_training_dataset
public :: cavg_quality_learn_diagnostics
public :: reset_cavg_quality_result

integer, parameter :: CAVG_QUALITY_NFEATS  = 14
integer, parameter :: CAVG_QUALITY_MAX_INTERACTIONS = (CAVG_QUALITY_NFEATS * (CAVG_QUALITY_NFEATS - 1)) / 2
real,    parameter :: EPS                  = 1.0e-6
real,    parameter :: CLIP_Z               = 4.0
character(len=*), parameter :: CAVG_MODEL_FAMILY_LINEAR               = 'linear_score'
character(len=*), parameter :: CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC    = 'pairwise_logistic'
character(len=*), parameter :: CAVG_RELATIONAL_SCHEMA_NONE            = 'none'
character(len=*), parameter :: CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1 = 'corr_knn_signal_v1'
integer,          parameter :: CAVG_RELATIONAL_DEFAULT_KNN             = 5
real,             parameter :: CAVG_RELATIONAL_DEFAULT_CORR_HP         = 100.0
real,             parameter :: CAVG_RELATIONAL_DEFAULT_CORR_LP         = 15.0
real,             parameter :: CAVG_RELATIONAL_DEFAULT_CORR_TRS        = 10.0

! Class-average quality contexts form a workflow ladder:
! - sieve: very small 2D chunks; conservative hard gates only, no learned model.
! - chunk: pre-cleaned 10-30k-particle chunks; logistic model rejection.
! - pool : highly cleaned merged chunk output; final model rejection before 3D.
character(len=*), parameter :: CAVG_QUALITY_CONTEXT_CHUNK             = 'chunk'
character(len=*), parameter :: CAVG_QUALITY_CONTEXT_SIEVE             = 'sieve'
character(len=*), parameter :: CAVG_QUALITY_CONTEXT_POOL              = 'pool'

integer, parameter :: CAVG_REJECT_REASON_NONE                  = 0
integer, parameter :: CAVG_REJECT_REASON_POP                   = 1
integer, parameter :: CAVG_REJECT_REASON_BAD_PIXELS            = 2
integer, parameter :: CAVG_REJECT_REASON_NO_COMPONENT          = 3
integer, parameter :: CAVG_REJECT_REASON_MASK_GEOMETRY         = 4
integer, parameter :: CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW    = 5
integer, parameter :: CAVG_REJECT_REASON_LOCVAR_FG_LOW         = 6
integer, parameter :: CAVG_REJECT_REASON_FUZZY_BALL_SIGNAL_LOW = 7

type :: cavg_quality_feature_def
    ! Keep these fixed-width fields within the current inventory limits:
    ! names/families are compact identifiers, descriptions are one-line labels.
    character(len=32)  :: name        = ''
    character(len=32)  :: direction   = 'higher_is_better'
    character(len=160) :: description = ''
    character(len=32)  :: family      = 'general'
end type cavg_quality_feature_def

type :: cavg_quality_model_spec
    character(len=64) :: name                         = ''
    character(len=32) :: context                      = 'chunk'
    character(len=64) :: feature_policy               = 'microchunk_plus_score_signal'
    character(len=32) :: model_family                 = CAVG_MODEL_FAMILY_LINEAR
    real              :: weights(CAVG_QUALITY_NFEATS) = 0.0
    real              :: intercept                    = 0.0
    real              :: linear_coefficients(CAVG_QUALITY_NFEATS) = 0.0
    integer           :: n_interactions               = 0
    integer           :: interaction_terms(CAVG_QUALITY_MAX_INTERACTIONS,2) = 0
    real              :: interaction_coefficients(CAVG_QUALITY_MAX_INTERACTIONS) = 0.0
    real              :: prob_threshold               = 0.5
    real              :: regularization_lambda        = 0.0
    real              :: calibration_temperature      = 1.0
    character(len=64) :: relational_feature_schema    = CAVG_RELATIONAL_SCHEMA_NONE
    integer           :: relational_knn               = 0
    real              :: relational_corr_hp           = 0.0
    real              :: relational_corr_lp           = 0.0
    real              :: relational_corr_trs          = 0.0
    real              :: relational_coefficient       = 0.0
    real              :: boundary_margin              = 0.0
    real              :: min_score_separation         = 0.0
    real              :: otsu_min_offset              = 0.0
    real              :: otsu_max_offset              = 0.0
    real              :: cluster_rescue_margin        = 0.0
    real              :: min_accept_frac              = 0.0
    logical           :: use_lowsep_otsu              = .false.
    logical           :: use_otsu_window              = .false.
    logical           :: use_cluster_rescue           = .false.
    logical           :: enforce_min_accept_frac      = .false.
end type cavg_quality_model_spec

type :: cavg_quality_result
    real,    allocatable :: raw(:,:)
    real,    allocatable :: features(:,:)
    real,    allocatable :: scores(:)
    integer, allocatable :: states(:)
    integer, allocatable :: labels(:)
    integer, allocatable :: medoids(:)
    integer, allocatable :: reasons(:)
    logical, allocatable :: hard_reject(:)
    real                 :: threshold        = 0.0
    real                 :: raw_threshold    = 0.0
    real                 :: threshold_offset = 0.0
    real                 :: separation       = 0.0
    integer              :: nclust           = 0
    integer              :: good_label       = 0
    logical              :: used_threshold   = .false.
    character(len=64)    :: model_name       = ''
    character(len=32)    :: soft_decision    = ''
    character(len=64)    :: soft_reason      = ''
contains
    procedure :: kill => reset_cavg_quality_result
end type cavg_quality_result

type :: cavg_quality_training_dataset
    character(len=LONGSTRLEN) :: fname      = ''
    character(len=LONGSTRLEN) :: dataset_id = ''
    character(len=64)         :: relational_feature_schema = CAVG_RELATIONAL_SCHEMA_NONE
    integer                   :: relational_knn = 0
    real                      :: relational_corr_hp = 0.0
    real                      :: relational_corr_lp = 0.0
    real                      :: relational_corr_trs = 0.0
    integer                   :: ncls       = 0
    real,    allocatable      :: features(:,:)
    real,    allocatable      :: relational_feature(:)
    integer, allocatable      :: manual_states(:)
    logical, allocatable      :: hard_reject(:)
end type cavg_quality_training_dataset

type :: cavg_quality_learn_diagnostics
    integer :: n_datasets          = 0
    integer :: n_scored_datasets   = 0
    integer :: n_weight_datasets   = 0
    integer :: n_recall_only       = 0
    integer :: n_specificity_only  = 0
    integer :: n_skipped           = 0
    integer :: total_fp            = 0
    integer :: total_fn            = 0
    integer :: n_lowsep            = 0
    integer :: lowsep_fp           = 0
    integer :: lowsep_fn           = 0
    integer :: n_single_cluster    = 0
    integer :: single_cluster_fp   = 0
    integer :: single_cluster_fn   = 0
    integer :: n_otsu_like         = 0
    integer :: otsu_like_fp        = 0
    integer :: otsu_like_fn        = 0
    integer :: n_rescue_like       = 0
    integer :: rescue_like_fp      = 0
    integer :: rescue_like_fn      = 0
    integer :: n_min_accept_like   = 0
    integer :: min_accept_like_fp  = 0
    integer :: min_accept_like_fn  = 0
    integer :: n_overfit_focus     = 0
    integer :: overfit_focus_bad   = 0
    integer :: overfit_focus_fp    = 0
end type cavg_quality_learn_diagnostics

contains

    subroutine reset_cavg_quality_result( quality )
        class(cavg_quality_result), intent(inout) :: quality
        if( allocated(quality%raw)               ) deallocate(quality%raw)
        if( allocated(quality%features)          ) deallocate(quality%features)
        if( allocated(quality%scores)            ) deallocate(quality%scores)
        if( allocated(quality%states)            ) deallocate(quality%states)
        if( allocated(quality%labels)            ) deallocate(quality%labels)
        if( allocated(quality%medoids)           ) deallocate(quality%medoids)
        if( allocated(quality%hard_reject)       ) deallocate(quality%hard_reject)
        if( allocated(quality%reasons)           ) deallocate(quality%reasons)
        quality%threshold        = 0.0
        quality%raw_threshold    = 0.0
        quality%threshold_offset = 0.0
        quality%separation       = 0.0
        quality%nclust           = 0
        quality%good_label       = 0
        quality%used_threshold   = .false.
        quality%model_name       = ''
        quality%soft_decision    = ''
        quality%soft_reason      = ''
    end subroutine reset_cavg_quality_result

end module simple_cavg_quality_types
