!@descr: shared types for class-average quality analysis
module simple_cavg_quality_types
use simple_defs, only: LONGSTRLEN
implicit none
private

public :: CAVG_QUALITY_NFEATS
public :: EPS
public :: CLIP_Z
public :: cavg_quality_feature_def
public :: cavg_quality_model_spec
public :: cavg_quality_result
public :: cavg_quality_training_dataset
public :: cavg_quality_learn_diagnostics
public :: reset_cavg_quality_result

integer, parameter :: CAVG_QUALITY_NFEATS  = 9
real,    parameter :: EPS                  = 1.0e-6
real,    parameter :: CLIP_Z               = 4.0

type :: cavg_quality_feature_def
    character(len=32)  :: name        = ''
    character(len=32)  :: direction   = 'higher_is_better'
    character(len=160) :: description = ''
    character(len=32)  :: family      = 'general'
end type cavg_quality_feature_def

type :: cavg_quality_model_spec
    character(len=64) :: name                         = ''
    character(len=32) :: context                      = 'chunk'
    character(len=32) :: feature_policy               = 'microchunk_plus_score_signal'
    real              :: weights(CAVG_QUALITY_NFEATS) = 0.0
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
    logical, allocatable :: hard_reject(:)
    real                 :: threshold        = 0.0
    real                 :: raw_threshold    = 0.0
    real                 :: threshold_offset = 0.0
    real                 :: separation       = 0.0
    integer              :: nclust           = 0
    integer              :: good_label       = 0
    logical              :: used_threshold   = .false.
    character(len=64)    :: model_name       = ''
    character(len=32)    :: model_context    = ''
    character(len=32)    :: soft_decision    = ''
    character(len=64)    :: soft_reason      = ''
contains
    procedure :: kill => reset_cavg_quality_result
end type cavg_quality_result

type :: cavg_quality_training_dataset
    character(len=LONGSTRLEN) :: fname      = ''
    character(len=LONGSTRLEN) :: dataset_id = ''
    integer                   :: ncls       = 0
    real,    allocatable      :: features(:,:)
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
end type cavg_quality_learn_diagnostics

contains

    subroutine reset_cavg_quality_result( quality )
        class(cavg_quality_result), intent(inout) :: quality
        if( allocated(quality%raw)         ) deallocate(quality%raw)
        if( allocated(quality%features)    ) deallocate(quality%features)
        if( allocated(quality%scores)      ) deallocate(quality%scores)
        if( allocated(quality%states)      ) deallocate(quality%states)
        if( allocated(quality%labels)      ) deallocate(quality%labels)
        if( allocated(quality%medoids)     ) deallocate(quality%medoids)
        if( allocated(quality%hard_reject) ) deallocate(quality%hard_reject)
        quality%threshold        = 0.0
        quality%raw_threshold    = 0.0
        quality%threshold_offset = 0.0
        quality%separation       = 0.0
        quality%nclust           = 0
        quality%good_label       = 0
        quality%used_threshold   = .false.
        quality%model_name       = ''
        quality%model_context    = ''
        quality%soft_decision    = ''
        quality%soft_reason      = ''
    end subroutine reset_cavg_quality_result

end module simple_cavg_quality_types
