!@descr: shared types for class-average quality analysis
module simple_cavg_quality_types
use simple_defs, only: LONGSTRLEN
implicit none
private

public :: CAVG_REJECTION_CHUNK
public :: CAVG_REJECTION_POOL
public :: CAVG_QUALITY_NFEATS
public :: cavg_quality_feature_def
public :: cavg_quality_model_spec
public :: cavg_quality_result
public :: cavg_quality_training_dataset
public :: reset_cavg_quality_result

integer, parameter :: CAVG_REJECTION_CHUNK = 1
integer, parameter :: CAVG_REJECTION_POOL  = 2
integer, parameter :: CAVG_QUALITY_NFEATS  = 12

type :: cavg_quality_feature_def
    character(len=32)  :: name        = ''
    character(len=32)  :: direction   = 'higher_is_better'
    character(len=160) :: description = ''
end type cavg_quality_feature_def

type :: cavg_quality_model_spec
    character(len=64) :: name                    = ''
    character(len=32) :: family                  = 'linear_boundary'
    character(len=32) :: context                 = 'chunk'
    integer           :: rejection_type          = CAVG_REJECTION_CHUNK
    real              :: weights(CAVG_QUALITY_NFEATS) = 0.0
    real              :: boundary_margin         = 0.0
    real              :: min_score_separation    = 0.0
    real              :: hist_dmat_weight        = 0.0
    real              :: otsu_min_offset         = 0.0
    real              :: otsu_max_offset         = 0.0
    real              :: cluster_rescue_margin   = 0.0
    real              :: min_accept_frac         = 0.0
    logical           :: use_lowsep_otsu         = .false.
    logical           :: use_otsu_window         = .false.
    logical           :: use_cluster_rescue      = .false.
    logical           :: enforce_min_accept_frac = .false.
end type cavg_quality_model_spec

type :: cavg_quality_result
    real,    allocatable :: raw(:,:)
    real,    allocatable :: features(:,:)
    real,    allocatable :: scores(:)
    integer, allocatable :: states(:)
    integer, allocatable :: labels(:)
    integer, allocatable :: medoids(:)
    logical, allocatable :: hard_reject(:)
    real,    allocatable :: hist_dmat(:,:)
    real                 :: threshold        = 0.0
    real                 :: raw_threshold    = 0.0
    real                 :: threshold_margin = 0.0
    real                 :: separation       = 0.0
    integer              :: nclust           = 0
    integer              :: good_label       = 0
    integer              :: rejection_type   = CAVG_REJECTION_CHUNK
    logical              :: used_threshold   = .false.
    character(len=64)    :: model_name       = ''
    character(len=32)    :: model_context    = ''
end type cavg_quality_result

type :: cavg_quality_training_dataset
    character(len=LONGSTRLEN) :: fname      = ''
    character(len=LONGSTRLEN) :: dataset_id = ''
    integer                   :: ncls       = 0
    real,    allocatable      :: features(:,:)
    real,    allocatable      :: hist_dmat(:,:)
    integer, allocatable      :: manual_states(:)
    logical, allocatable      :: hard_reject(:)
end type cavg_quality_training_dataset

contains

    subroutine reset_cavg_quality_result( quality )
        type(cavg_quality_result), intent(inout) :: quality
        if( allocated(quality%raw)         ) deallocate(quality%raw)
        if( allocated(quality%features)    ) deallocate(quality%features)
        if( allocated(quality%scores)      ) deallocate(quality%scores)
        if( allocated(quality%states)      ) deallocate(quality%states)
        if( allocated(quality%labels)      ) deallocate(quality%labels)
        if( allocated(quality%medoids)     ) deallocate(quality%medoids)
        if( allocated(quality%hard_reject) ) deallocate(quality%hard_reject)
        if( allocated(quality%hist_dmat)   ) deallocate(quality%hist_dmat)
        quality%threshold        = 0.0
        quality%raw_threshold    = 0.0
        quality%threshold_margin = 0.0
        quality%separation       = 0.0
        quality%nclust           = 0
        quality%good_label       = 0
        quality%rejection_type   = CAVG_REJECTION_CHUNK
        quality%used_threshold   = .false.
        quality%model_name       = ''
        quality%model_context    = ''
    end subroutine reset_cavg_quality_result

end module simple_cavg_quality_types
