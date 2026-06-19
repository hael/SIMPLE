!@descr: instantiable class-average quality decision model
module simple_cavg_quality_model
use simple_defs,               only: LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string_utils,       only: str_is_true, lowercase, uppercase, &
    fortran_symbol_from_string, fortran_quote, fortran_logical
use simple_clustering_utils,   only: cluster_dmat
use simple_srch_sort_loc,      only: hpsort
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, cavg_quality_model_spec, cavg_quality_result
use simple_cavg_quality_stats, only: normalize_quality_dmat
implicit none
private
#include "simple_local_flags.inc"

public :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT
public :: CAVG_QUALITY_MODEL_POOL_DEFAULT
public :: CAVG_QUALITY_MODEL_MICROCHUNK_P1
public :: CAVG_QUALITY_MODEL_MICROCHUNK_P2
public :: CAVG_QUALITY_MODEL_CHUNK_LP4
public :: cavg_quality_model
public :: cavg_quality_model_spec
public :: cavg_quality_classify_cache
public :: build_classify_cache
public :: apply_cached_decision_to_quality
public :: cached_decision_confusion
public :: kill_classify_cache
public :: write_cavg_quality_model_builtin_code

! Classify-cache decision modes. The cache captures everything in the
! classification pipeline that depends only on (features, hard_reject,
! weights) and not on threshold/Otsu spec parameters, so the grid search
! in learn mode can evaluate many candidate specs against the same cache
! without redoing the k-medoids / Otsu work.
integer, parameter :: CACHE_DECISION_EMPTY     = 1  ! no trainable rows
integer, parameter :: CACHE_DECISION_FORCED    = 2  ! single-cluster accept-all path
integer, parameter :: CACHE_DECISION_CLUSTERED = 3  ! two-cluster artifacts populated

! Why the forced (single-cluster) path was taken. Mirrors the reason
! strings used by accept_fit_as_single_cluster so the cached path can
! reproduce the same soft_reason annotations.
integer, parameter :: CACHE_FORCED_TOO_FEW             = 1
integer, parameter :: CACHE_FORCED_FLAT_DMAT           = 2
integer, parameter :: CACHE_FORCED_INVALID_TWO_CLUSTER = 3

type :: cavg_quality_classify_cache
    integer              :: ncls            = 0
    integer              :: nfit            = 0
    integer              :: decision_mode   = 0
    integer              :: forced_reason   = 0
    integer, allocatable :: inds(:)            ! nfit, indices of non-hard-rejected classes
    real,    allocatable :: scores(:)          ! ncls, with hard-reject rows clipped to -CLIP_Z
    real,    allocatable :: score_fit(:)       ! nfit
    logical, allocatable :: hard_reject(:)     ! ncls
    integer, allocatable :: labels_fit(:)      ! nfit (1/2 in CLUSTERED mode)
    integer, allocatable :: medoids_fit(:)     ! medoid indices into score_fit
    real,    allocatable :: score_fit_sorted(:) ! nfit, ascending score_fit values
    integer              :: good_fit_label  = 0
    integer              :: bad_fit_label   = 0
    real                 :: separation      = 0.0
    real                 :: raw_threshold   = 0.0
    real                 :: otsu_threshold  = 0.0
    real                 :: otsu_separation = 0.0
    logical              :: otsu_ok         = .false.
end type cavg_quality_classify_cache

type :: cavg_quality_cached_decision
    real              :: threshold        = 0.0
    real              :: raw_threshold    = 0.0
    real              :: threshold_offset = 0.0
    real              :: rescue_threshold = 0.0
    real              :: floor_threshold  = 0.0
    logical           :: use_rescue       = .false.
    logical           :: floor_active     = .false.
    logical           :: single_cluster   = .false.
    logical           :: used_threshold   = .false.
    character(len=32) :: soft_decision    = 'hard_only'
    character(len=64) :: soft_reason      = 'initial'
end type cavg_quality_cached_decision

! Built-in presets are complete model specifications. To promote a learned
! model into the code, add a named preset and include it in builtin_names.
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT = 'chunk_default_v2'
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_LP4     = 'chunk_lp4'
character(len=*), parameter :: CAVG_QUALITY_MODEL_POOL_DEFAULT  = 'pool_default_v2'
character(len=*), parameter :: CAVG_QUALITY_MODEL_MICROCHUNK_P1 = 'microchunk_p1'
character(len=*), parameter :: CAVG_QUALITY_MODEL_MICROCHUNK_P2 = 'microchunk_p2'
character(len=*), parameter :: BUILTIN_MODEL_NAMES = CAVG_QUALITY_MODEL_CHUNK_DEFAULT//'|'//&
    CAVG_QUALITY_MODEL_CHUNK_LP4//'|'//&
    CAVG_QUALITY_MODEL_POOL_DEFAULT//'|'//&
    CAVG_QUALITY_MODEL_MICROCHUNK_P1//'|'//&
    CAVG_QUALITY_MODEL_MICROCHUNK_P2

real, parameter :: CHUNK_OTSU_MIN_OFFSET   = 0.25
real, parameter :: CHUNK_OTSU_MAX_OFFSET   = 0.50
real, parameter :: CLUSTER_RESCUE_MARGIN   = 0.20

! Pool default for batch-style class-average rejection.
character(len=*), parameter :: POOL_V2_FEATURE_POLICY = 'microchunk_plus_score_signal'
real, parameter :: CAVG_QUALITY_POOL_V2_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    1.826303E-01, 2.379928E-01, 0.000000E+00, 7.328372E-02, &
    6.126274E-02, 1.257789E-01, 1.572098E-01, 2.139509E-02, &
    1.404466E-01, 0.000000E+00, 0.000000E+00, 0.000000E+00 ]
real, parameter :: POOL_V2_BOUNDARY_MARGIN      = 1.70
real, parameter :: POOL_V2_MIN_SCORE_SEPARATION = 0.20
real, parameter :: POOL_V2_OTSU_MIN_OFFSET      = 0.15
real, parameter :: POOL_V2_OTSU_MAX_OFFSET      = 0.50
real, parameter :: POOL_V2_MIN_ACCEPT_FRAC      = 0.95

! Chunk default for stream-style class-average rejection. Hard validity
! failures reject before fitting; the weights below describe the trainable
! quality boundary for classes that pass those validity checks.
character(len=*), parameter :: CHUNK_V2_FEATURE_POLICY = 'microchunk_plus_score'
real, parameter :: CAVG_QUALITY_CHUNK_V2_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    1.355581E-01, 1.808395E-01, 4.942806E-02, 1.635655E-01, &
    1.726983E-01, 2.108072E-01, 0.000000E+00, 8.710331E-02, &
    0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 ]
real, parameter :: CHUNK_V2_BOUNDARY_MARGIN      =  0.15
real, parameter :: CHUNK_V2_MIN_SCORE_SEPARATION =  0.15
real, parameter :: CHUNK_V2_OTSU_MIN_OFFSET      =  0.35
real, parameter :: CHUNK_V2_OTSU_MAX_OFFSET      =  0.50

! Chunk model learned from lp4-style class-average selection references in
! /Users/elmlundho/cavgs_quality/chunk_new_train2, including the *_2* datasets.
character(len=*), parameter :: CHUNK_LP4_FEATURE_POLICY = 'microchunk_plus_signal'
real, parameter :: CAVG_QUALITY_CHUNK_LP4_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    1.256393E-01, 1.026731E-01, 4.134843E-02, 1.367731E-01, &
    1.996999E-01, 0.000000E+00, 2.763570E-02, 1.404427E-01, &
    2.257877E-01, 0.000000E+00, 0.000000E+00, 0.000000E+00 ]
real, parameter :: CHUNK_LP4_BOUNDARY_MARGIN      = -0.15
real, parameter :: CHUNK_LP4_MIN_SCORE_SEPARATION =  0.15
real, parameter :: CHUNK_LP4_OTSU_MIN_OFFSET      =  0.35
real, parameter :: CHUNK_LP4_OTSU_MAX_OFFSET      =  0.40

type :: cavg_quality_model
    character(len=64) :: name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
    character(len=32) :: context                 = 'chunk'
    character(len=64) :: feature_policy          = CHUNK_V2_FEATURE_POLICY
    real              :: weights(CAVG_QUALITY_NFEATS) = CAVG_QUALITY_CHUNK_V2_WEIGHTS
    real              :: boundary_margin         = CHUNK_V2_BOUNDARY_MARGIN
    real              :: min_score_separation    = CHUNK_V2_MIN_SCORE_SEPARATION
    real              :: otsu_min_offset         = CHUNK_V2_OTSU_MIN_OFFSET
    real              :: otsu_max_offset         = CHUNK_V2_OTSU_MAX_OFFSET
    real              :: cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
    real              :: min_accept_frac         = 0.0
    logical           :: use_lowsep_otsu         = .true.
    logical           :: use_otsu_window         = .true.
    logical           :: use_cluster_rescue      = .false.
    logical           :: enforce_min_accept_frac = .false.
contains
    procedure :: init_preset
    procedure :: init_spec
    procedure :: get_spec
    procedure :: normalize
    procedure :: classify
    procedure :: write => write_model
    procedure :: read  => read_model
    ! Destructor-style clear: after kill, call init_preset or init_spec
    ! before using the model again.
    procedure :: kill  => kill_model
end type cavg_quality_model

contains

    subroutine init_preset( self, preset_name )
        class(cavg_quality_model), intent(inout) :: self
        character(len=*),          intent(in)    :: preset_name
        type(cavg_quality_model_spec) :: spec
        spec = builtin_spec(preset_name)
        call self%init_spec(spec)
    end subroutine init_preset

    function builtin_spec( preset_name ) result( spec )
        character(len=*), intent(in) :: preset_name
        type(cavg_quality_model_spec) :: spec
        character(len=LONGSTRLEN) :: errmsg
        select case(trim(preset_name))
            case(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
                spec = chunk_default_v2_model_spec()
            case(CAVG_QUALITY_MODEL_CHUNK_LP4)
                spec = chunk_lp4_model_spec()
            case(CAVG_QUALITY_MODEL_POOL_DEFAULT)
                spec = pool_default_model_spec()
            case(CAVG_QUALITY_MODEL_MICROCHUNK_P1)
                spec = microchunk_p1_model_spec()
            case(CAVG_QUALITY_MODEL_MICROCHUNK_P2)
                spec = microchunk_p2_model_spec()
            case default
                errmsg = 'unknown class-average quality model preset: '//trim(preset_name)//&
                         '; available presets: '//trim(builtin_names())
                THROW_HARD(trim(errmsg))
        end select
    end function builtin_spec

    function builtin_names() result( names )
        character(len=LONGSTRLEN) :: names
        names = BUILTIN_MODEL_NAMES
    end function builtin_names

    subroutine init_spec( self, spec )
        class(cavg_quality_model),     intent(inout) :: self
        type(cavg_quality_model_spec), intent(in)    :: spec
        self%name                    = trim(spec%name)
        self%context                 = lowercase(trim(spec%context))
        call assert_valid_model_context(self%context)
        self%feature_policy          = trim(spec%feature_policy)
        self%weights                 = spec%weights
        self%boundary_margin         = spec%boundary_margin
        self%min_score_separation    = spec%min_score_separation
        self%otsu_min_offset         = spec%otsu_min_offset
        self%otsu_max_offset         = spec%otsu_max_offset
        self%cluster_rescue_margin   = spec%cluster_rescue_margin
        self%min_accept_frac         = spec%min_accept_frac
        self%use_lowsep_otsu         = spec%use_lowsep_otsu
        self%use_otsu_window         = spec%use_otsu_window
        self%use_cluster_rescue      = spec%use_cluster_rescue
        self%enforce_min_accept_frac = spec%enforce_min_accept_frac
        call self%normalize()
    end subroutine init_spec

    function get_spec( self ) result( spec )
        class(cavg_quality_model), intent(in) :: self
        type(cavg_quality_model_spec) :: spec
        spec%name                    = self%name
        spec%context                 = self%context
        spec%feature_policy          = self%feature_policy
        spec%weights                 = self%weights
        spec%boundary_margin         = self%boundary_margin
        spec%min_score_separation    = self%min_score_separation
        spec%otsu_min_offset         = self%otsu_min_offset
        spec%otsu_max_offset         = self%otsu_max_offset
        spec%cluster_rescue_margin   = self%cluster_rescue_margin
        spec%min_accept_frac         = self%min_accept_frac
        spec%use_lowsep_otsu         = self%use_lowsep_otsu
        spec%use_otsu_window         = self%use_otsu_window
        spec%use_cluster_rescue      = self%use_cluster_rescue
        spec%enforce_min_accept_frac = self%enforce_min_accept_frac
    end function get_spec

    function chunk_default_v2_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
        spec%context                 = 'chunk'
        spec%feature_policy          = CHUNK_V2_FEATURE_POLICY
        spec%weights                 = CAVG_QUALITY_CHUNK_V2_WEIGHTS
        spec%boundary_margin         = CHUNK_V2_BOUNDARY_MARGIN
        spec%min_score_separation    = CHUNK_V2_MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = CHUNK_V2_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK_V2_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = 0.0
        spec%use_lowsep_otsu         = .true.
        spec%use_otsu_window         = .true.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function chunk_default_v2_model_spec

    function chunk_lp4_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_CHUNK_LP4
        spec%context                 = 'chunk'
        spec%feature_policy          = CHUNK_LP4_FEATURE_POLICY
        spec%weights                 = CAVG_QUALITY_CHUNK_LP4_WEIGHTS
        spec%boundary_margin         = CHUNK_LP4_BOUNDARY_MARGIN
        spec%min_score_separation    = CHUNK_LP4_MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = CHUNK_LP4_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK_LP4_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = 0.0
        spec%use_lowsep_otsu         = .true.
        spec%use_otsu_window         = .true.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function chunk_lp4_model_spec

    function pool_default_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_POOL_DEFAULT
        spec%context                 = 'pool'
        spec%feature_policy          = POOL_V2_FEATURE_POLICY
        spec%weights                 = CAVG_QUALITY_POOL_V2_WEIGHTS
        spec%boundary_margin         = POOL_V2_BOUNDARY_MARGIN
        spec%min_score_separation    = POOL_V2_MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = POOL_V2_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = POOL_V2_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = POOL_V2_MIN_ACCEPT_FRAC
        spec%use_lowsep_otsu         = .true.
        spec%use_otsu_window         = .false.
        spec%use_cluster_rescue      = .true.
        spec%enforce_min_accept_frac = .true.
    end function pool_default_model_spec

    function microchunk_p1_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_MICROCHUNK_P1
        spec%context                 = 'chunk'
        spec%feature_policy          = 'microchunk_plus_score_signal'
        spec%weights                 = [ &
              1.152264E-02,   9.737433E-02,   6.403663E-02,   2.116542E-01, &
              1.794099E-01,   1.902685E-01,   1.129871E-01,   6.321478E-03, &
              1.264252E-01,   0.000000E+00,   0.000000E+00,   0.000000E+00 ]
        spec%boundary_margin         =   8.000000E-01
        spec%min_score_separation    =   2.000000E-01
        spec%otsu_min_offset         =   2.500000E-01
        spec%otsu_max_offset         =   5.000000E-01
        spec%cluster_rescue_margin   =   2.000000E-01
        spec%min_accept_frac         =   0.000000E+00
        spec%use_lowsep_otsu         = .true.
        spec%use_otsu_window         = .true.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function microchunk_p1_model_spec

    function microchunk_p2_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_MICROCHUNK_P2
        spec%context                 = 'chunk'
        spec%feature_policy          = 'microchunk_plus_score_signal_texture'
        spec%weights                 = [ &
              9.300453E-02,   9.552415E-02,   1.010380E-01,   9.765499E-02, &
              1.325033E-01,   1.486267E-01,   1.026389E-01,   4.956052E-02, &
              1.352742E-01,   4.417469E-02,   0.000000E+00,   0.000000E+00 ]
        spec%boundary_margin         =   1.500000E-01
        spec%min_score_separation    =   3.000000E-01
        spec%otsu_min_offset         =   2.500000E-01
        spec%otsu_max_offset         =   4.000000E-01
        spec%cluster_rescue_margin   =   2.000000E-01
        spec%min_accept_frac         =   0.000000E+00
        spec%use_lowsep_otsu         = .true.
        spec%use_otsu_window         = .true.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function microchunk_p2_model_spec

    subroutine normalize( self )
        class(cavg_quality_model), intent(inout) :: self
        self%weights = max(0.0, self%weights)
        if( sum(self%weights) > EPS )then
            self%weights = self%weights / sum(self%weights)
        else
            self%weights = 1.0 / real(CAVG_QUALITY_NFEATS)
            self%weights = self%weights / sum(self%weights)
        endif
    end subroutine normalize

    subroutine classify( self, quality )
        class(cavg_quality_model), intent(in)    :: self
        type(cavg_quality_result), intent(inout) :: quality
        if( .not. allocated(quality%features)    ) THROW_HARD('classify: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('classify: missing hard-reject mask')
        quality%model_name     = self%name
        quality%model_context  = self%context
        call apply_linear_boundary(quality, self)
    end subroutine classify

    subroutine write_model( self, fname )
        class(cavg_quality_model), intent(in) :: self
        character(len=*),          intent(in) :: fname
        integer :: funit, i
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# model_cavgs_rejection model'
        write(funit,'(A)') 'model_version=6'
        write(funit,'(A,A)') 'name=', trim(self%name)
        write(funit,'(A,A)') 'context=', trim(self%context)
        write(funit,'(A,A)') 'feature_policy=', trim(self%feature_policy)
        write(funit,'(A)', advance='no') 'feature_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') self%weights(i)
        end do
        write(funit,*)
        write(funit,'(A,ES14.6)') 'boundary_margin=', self%boundary_margin
        write(funit,'(A,ES14.6)') 'min_score_separation=', self%min_score_separation
        write(funit,'(A,ES14.6)') 'otsu_min_offset=', self%otsu_min_offset
        write(funit,'(A,ES14.6)') 'otsu_max_offset=', self%otsu_max_offset
        write(funit,'(A,ES14.6)') 'cluster_rescue_margin=', self%cluster_rescue_margin
        write(funit,'(A,ES14.6)') 'min_accept_frac=', self%min_accept_frac
        write(funit,'(A,L1)') 'use_lowsep_otsu=', self%use_lowsep_otsu
        write(funit,'(A,L1)') 'use_otsu_window=', self%use_otsu_window
        write(funit,'(A,L1)') 'use_cluster_rescue=', self%use_cluster_rescue
        write(funit,'(A,L1)') 'enforce_min_accept_frac=', self%enforce_min_accept_frac
        close(funit)
    end subroutine write_model

    subroutine write_cavg_quality_model_builtin_code( model, fname )
        type(cavg_quality_model), intent(in) :: model
        character(len=*),         intent(in) :: fname
        character(len=64)  :: symbol, func_name
        character(len=128) :: const_name
        integer :: funit
        symbol     = fortran_symbol_from_string(model%name, fallback='quality_model', &
            prefix_if_invalid_start='model_', max_symbol_len=40)
        func_name  = trim(symbol)//'_model_spec'
        const_name = 'CAVG_QUALITY_MODEL_'//trim(uppercase(symbol))
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '! model_cavgs_rejection built-in model promotion snippet'
        write(funit,'(A)') '! Generated from learned model: '//trim(model%name)
        write(funit,'(A)') '! Review the validation report before adding this preset to the library.'
        write(funit,'(A)') ''
        write(funit,'(A)') '! 1. Add this constant near the built-in model names in simple_cavg_quality_model.f90:'
        write(funit,'(A,A,A,A,A)') 'character(len=*), parameter :: ', trim(const_name), ' = ', &
            trim(fortran_quote(model%name)), ''
        write(funit,'(A)') ''
        write(funit,'(A)') '! 2. Append this name to BUILTIN_MODEL_NAMES:'
        write(funit,'(A,A)') '!     //''|''//', trim(const_name)
        write(funit,'(A)') ''
        write(funit,'(A)') '! 3. Add this case in builtin_spec:'
        write(funit,'(A,A,A)') '            case(', trim(const_name), ')'
        write(funit,'(A,A,A)') '                spec = ', trim(func_name), '()'
        write(funit,'(A)') ''
        write(funit,'(A)') '! 4. Add this function next to the other built-in model specs:'
        call write_model_spec_function(funit, model, trim(func_name), trim(const_name))
        write(funit,'(A)') ''
        write(funit,'(A)') '! 5. Add the model name to the quality_model UI/help option lists:'
        write(funit,'(A)') '!     src/main/ui/simple_ui_params_common.f90'
        write(funit,'(A)') '!     src/main/params/simple_parameters.f90'
        write(funit,'(A,A,A)') '!     option token: ', trim(model%name), ''
        close(funit)
    end subroutine write_cavg_quality_model_builtin_code

    subroutine write_model_spec_function( funit, model, func_name, const_name )
        integer,                  intent(in) :: funit
        type(cavg_quality_model), intent(in) :: model
        character(len=*),         intent(in) :: func_name, const_name
        write(funit,'(A,A,A)') '    function ', trim(func_name), '() result( spec )'
        write(funit,'(A)') '        type(cavg_quality_model_spec) :: spec'
        write(funit,'(A,A)') '        spec%name                    = ', trim(const_name)
        write(funit,'(A,A)') '        spec%context                 = ', trim(fortran_quote(model%context))
        write(funit,'(A,A)') '        spec%feature_policy          = ', trim(fortran_quote(model%feature_policy))
        call write_weights_assignment(funit, model%weights)
        write(funit,'(A,ES14.6)') '        spec%boundary_margin         = ', model%boundary_margin
        write(funit,'(A,ES14.6)') '        spec%min_score_separation    = ', model%min_score_separation
        write(funit,'(A,ES14.6)') '        spec%otsu_min_offset         = ', model%otsu_min_offset
        write(funit,'(A,ES14.6)') '        spec%otsu_max_offset         = ', model%otsu_max_offset
        write(funit,'(A,ES14.6)') '        spec%cluster_rescue_margin   = ', model%cluster_rescue_margin
        write(funit,'(A,ES14.6)') '        spec%min_accept_frac         = ', model%min_accept_frac
        write(funit,'(A,A)') '        spec%use_lowsep_otsu         = ', trim(fortran_logical(model%use_lowsep_otsu))
        write(funit,'(A,A)') '        spec%use_otsu_window         = ', trim(fortran_logical(model%use_otsu_window))
        write(funit,'(A,A)') '        spec%use_cluster_rescue      = ', trim(fortran_logical(model%use_cluster_rescue))
        write(funit,'(A,A)') '        spec%enforce_min_accept_frac = ', trim(fortran_logical(model%enforce_min_accept_frac))
        write(funit,'(A,A,A)') '    end function ', trim(func_name), ''
    end subroutine write_model_spec_function

    subroutine write_weights_assignment( funit, weights )
        integer, intent(in) :: funit
        real,    intent(in) :: weights(:)
        integer :: i
        write(funit,'(A)') '        spec%weights                 = [ &'
        write(funit,'(A)', advance='no') '            '
        do i = 1, size(weights)
            write(funit,'(ES14.6)', advance='no') weights(i)
            if( i < size(weights) ) write(funit,'(A)', advance='no') ', '
            if( mod(i, 4) == 0 .and. i < size(weights) )then
                write(funit,'(A)') '&'
                write(funit,'(A)', advance='no') '            '
            endif
        end do
        write(funit,'(A)') ' ]'
    end subroutine write_weights_assignment

    subroutine read_model( self, fname )
        class(cavg_quality_model), intent(inout) :: self
        character(len=*),          intent(in)    :: fname
        character(len=XLONGSTRLEN) :: line
        character(len=LONGSTRLEN)  :: key, val, preset_name
        integer :: funit, ios, parse_ios
        logical :: have_preset, ok_line
        ! Model files are complete model definitions. Start from chunk defaults,
        ! apply any preset found in the file, then apply explicit key overrides.
        call self%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
        open(newunit=funit, file=trim(fname), status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('read_model: failed to open '//trim(fname))
        have_preset = .false.
        preset_name = ''
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            call parse_model_key_value(line, key, val, ok_line)
            if( .not. ok_line ) cycle
            if( trim(key) == 'preset' )then
                preset_name = trim(val)
                have_preset = .true.
            endif
        end do
        if( have_preset ) call self%init_preset(trim(preset_name))
        rewind(funit)
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            call parse_model_key_value(line, key, val, ok_line)
            if( .not. ok_line ) cycle
            select case(trim(key))
                case('model_version')
                    cycle
                case('preset')
                    cycle
                case('name')
                    self%name = trim(val)
                case('context')
                    self%context = lowercase(trim(val))
                case('feature_policy', 'feature_family_set')
                    self%feature_policy = trim(val)
                case('feature_weights')
                    call read_feature_weights(val, self%weights)
                case('boundary_margin')
                    read(val,*,iostat=parse_ios) self%boundary_margin
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse boundary_margin')
                case('min_score_separation')
                    read(val,*,iostat=parse_ios) self%min_score_separation
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse min_score_separation')
                case('otsu_min_offset')
                    read(val,*,iostat=parse_ios) self%otsu_min_offset
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse otsu_min_offset')
                case('otsu_max_offset')
                    read(val,*,iostat=parse_ios) self%otsu_max_offset
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse otsu_max_offset')
                case('cluster_rescue_margin')
                    read(val,*,iostat=parse_ios) self%cluster_rescue_margin
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse cluster_rescue_margin')
                case('min_accept_frac')
                    read(val,*,iostat=parse_ios) self%min_accept_frac
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse min_accept_frac')
                case('use_lowsep_otsu')
                    self%use_lowsep_otsu = str_is_true(val)
                case('use_otsu_window')
                    self%use_otsu_window = str_is_true(val)
                case('use_cluster_rescue')
                    self%use_cluster_rescue = str_is_true(val)
                case('enforce_min_accept_frac')
                    self%enforce_min_accept_frac = str_is_true(val)
                case default
                    THROW_HARD('read_model: unknown key in model file: '//trim(key))
            end select
        end do
        close(funit)
        call assert_valid_model_context(self%context)
        call self%normalize()
    end subroutine read_model

    subroutine read_feature_weights( val, weights )
        character(len=*), intent(in)    :: val
        real,             intent(inout) :: weights(CAVG_QUALITY_NFEATS)
        character(len=LONGSTRLEN) :: errmsg
        integer :: nvals
        real    :: parsed(CAVG_QUALITY_NFEATS)
        call parse_feature_weight_values(val, parsed, nvals)
        if( nvals == 0 ) THROW_HARD('read_model: feature_weights has no values')
        if( nvals > CAVG_QUALITY_NFEATS )then
            write(errmsg,'(A,I0,A,I0)') 'read_model: feature_weights expected at most ', CAVG_QUALITY_NFEATS, &
                ' values, got ', nvals
            THROW_HARD(trim(errmsg))
        endif
        weights = parsed
    end subroutine read_feature_weights

    subroutine parse_feature_weight_values( val, weights, nvals )
        character(len=*), intent(in)  :: val
        real,             intent(out) :: weights(CAVG_QUALITY_NFEATS)
        integer,          intent(out) :: nvals
        character(len=LONGSTRLEN) :: work, token
        integer :: isep, ios
        weights = 0.0
        nvals   = 0
        work    = adjustl(trim(val))
        do while( len_trim(work) > 0 )
            isep = scan(work, ', ')
            if( isep == 1 )then
                work = adjustl(work(2:))
                cycle
            else if( isep > 1 )then
                token = work(1:isep-1)
                work  = adjustl(work(isep+1:))
            else
                token = work
                work  = ''
            endif
            if( len_trim(token) == 0 ) cycle
            nvals = nvals + 1
            if( nvals > CAVG_QUALITY_NFEATS ) cycle
            read(token,*,iostat=ios) weights(nvals)
            if( ios /= 0 ) THROW_HARD('read_model: failed to parse feature_weights')
        end do
    end subroutine parse_feature_weight_values

    subroutine parse_model_key_value( line, key, val, ok )
        character(len=*), intent(in)  :: line
        character(len=*), intent(out) :: key, val
        logical,          intent(out) :: ok
        character(len=XLONGSTRLEN) :: tmp
        integer :: ieq
        key = ''
        val = ''
        ok  = .false.
        tmp = adjustl(line)
        if( len_trim(tmp) == 0 ) return
        if( tmp(1:1) == '#' ) return
        ieq = index(tmp, '=')
        if( ieq <= 1 ) return
        key = adjustl(trim(tmp(1:ieq-1)))
        val = adjustl(trim(tmp(ieq+1:)))
        ok  = .true.
    end subroutine parse_model_key_value

    subroutine kill_model( self )
        class(cavg_quality_model), intent(inout) :: self
        self%name                    = ''
        self%context                 = ''
        self%feature_policy          = ''
        self%weights                 = 0.0
        self%boundary_margin         = 0.0
        self%min_score_separation    = 0.0
        self%otsu_min_offset         = 0.0
        self%otsu_max_offset         = 0.0
        self%cluster_rescue_margin   = 0.0
        self%min_accept_frac         = 0.0
        self%use_lowsep_otsu         = .false.
        self%use_otsu_window         = .false.
        self%use_cluster_rescue      = .false.
        self%enforce_min_accept_frac = .false.
    end subroutine kill_model

    subroutine apply_linear_boundary( quality, model )
        type(cavg_quality_result), intent(inout) :: quality
        class(cavg_quality_model), intent(in) :: model
        type(cavg_quality_classify_cache) :: cache
        if( .not. allocated(quality%features)    ) THROW_HARD('apply_linear_boundary: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('apply_linear_boundary: missing hard-reject mask')
        if( size(quality%features, dim=2) /= CAVG_QUALITY_NFEATS ) THROW_HARD('apply_linear_boundary: invalid feature count')
        if( size(quality%hard_reject) /= size(quality%features, dim=1) ) THROW_HARD('apply_linear_boundary: invalid mask size')
        call build_classify_cache(quality%features, quality%hard_reject, model%weights, cache)
        call apply_cached_decision_to_quality(cache, model, quality)
        call kill_classify_cache(cache)
    end subroutine apply_linear_boundary

    ! Heavy spec-independent portion of the classification pipeline:
    ! scores, feature-masked pairwise distances, k-medoids labels, raw
    ! threshold, and Otsu threshold. Within a fixed feature policy (and
    ! therefore fixed weights) the cache is reusable across an arbitrary
    ! number of candidate model specs.
    subroutine build_classify_cache( features, hard_reject, weights, cache )
        real,                              intent(in)    :: features(:,:)
        logical,                           intent(in)    :: hard_reject(:)
        real,                              intent(in)    :: weights(:)
        type(cavg_quality_classify_cache), intent(inout) :: cache
        real,    allocatable :: dmat(:,:)
        real,    allocatable :: feats_fit(:,:)
        integer              :: ncls, nfit, i, j, nclust
        real                 :: d, score1, score2
        logical              :: dmat_ok
        call kill_classify_cache(cache)
        ncls = size(features, dim=1)
        if( size(features, dim=2) /= CAVG_QUALITY_NFEATS ) THROW_HARD('build_classify_cache: invalid feature count')
        if( size(hard_reject)     /= ncls                ) THROW_HARD('build_classify_cache: invalid mask size')
        if( size(weights)         /= CAVG_QUALITY_NFEATS ) THROW_HARD('build_classify_cache: invalid weight size')
        cache%ncls = ncls
        allocate(cache%hard_reject(ncls), source=hard_reject)
        allocate(cache%scores(ncls))
        cache%scores = matmul(features, weights)
        where( hard_reject ) cache%scores = -CLIP_Z
        nfit = count(.not. hard_reject)
        cache%nfit = nfit
        if( nfit == 0 )then
            cache%decision_mode = CACHE_DECISION_EMPTY
            return
        endif
        allocate(cache%inds(nfit))
        cache%inds = pack([(i, i=1,ncls)], .not. hard_reject)
        allocate(cache%score_fit(nfit))
        do i = 1, nfit
            cache%score_fit(i) = cache%scores(cache%inds(i))
        end do
        allocate(cache%score_fit_sorted(nfit), source=cache%score_fit)
        call hpsort(cache%score_fit_sorted)
        if( nfit < 4 )then
            cache%decision_mode = CACHE_DECISION_FORCED
            cache%forced_reason = CACHE_FORCED_TOO_FEW
            return
        endif
        allocate(feats_fit(nfit, CAVG_QUALITY_NFEATS))
        do i = 1, nfit
            feats_fit(i,:) = features(cache%inds(i),:)
        end do
        allocate(dmat(nfit, nfit), source=0.0)
        !$omp parallel do default(shared) private(i,j,d) schedule(static) proc_bind(close)
        do i = 1, nfit - 1
            do j = i + 1, nfit
                d = feature_vector_distance(feats_fit(i,:), feats_fit(j,:), weights)
                dmat(i,j) = d
                dmat(j,i) = d
            end do
        end do
        !$omp end parallel do
        call normalize_quality_dmat(dmat, dmat_ok)
        if( .not. dmat_ok )then
            cache%decision_mode = CACHE_DECISION_FORCED
            cache%forced_reason = CACHE_FORCED_FLAT_DMAT
            deallocate(dmat, feats_fit)
            return
        endif
        nclust = 2
        call cluster_dmat(dmat, 'kmed', nclust, cache%medoids_fit, cache%labels_fit)
        if( .not. two_cluster_result_is_valid(nclust, cache%labels_fit, cache%medoids_fit, nfit) )then
            cache%decision_mode = CACHE_DECISION_FORCED
            cache%forced_reason = CACHE_FORCED_INVALID_TWO_CLUSTER
            if( allocated(cache%labels_fit)  ) deallocate(cache%labels_fit)
            if( allocated(cache%medoids_fit) ) deallocate(cache%medoids_fit)
            deallocate(dmat, feats_fit)
            return
        endif
        cache%decision_mode = CACHE_DECISION_CLUSTERED
        score1 = mean_score_for_label(cache%score_fit, cache%labels_fit, 1)
        score2 = mean_score_for_label(cache%score_fit, cache%labels_fit, 2)
        if( score1 > score2 + EPS )then
            cache%good_fit_label = 1
            cache%bad_fit_label  = 2
        else if( score2 > score1 + EPS )then
            cache%good_fit_label = 2
            cache%bad_fit_label  = 1
        else
            call choose_tied_good_label(cache%score_fit, cache%labels_fit, cache%medoids_fit, &
                cache%good_fit_label, cache%bad_fit_label)
        endif
        cache%separation    = abs(score1 - score2)
        cache%raw_threshold = 0.5 * (mean_score_for_label(cache%score_fit, cache%labels_fit, cache%good_fit_label) + &
                                     mean_score_for_label(cache%score_fit, cache%labels_fit, cache%bad_fit_label))
        call otsu_score_threshold(cache%score_fit, cache%otsu_threshold, cache%otsu_separation, cache%otsu_ok)
        deallocate(dmat, feats_fit)
    end subroutine build_classify_cache

    subroutine kill_classify_cache( cache )
        type(cavg_quality_classify_cache), intent(inout) :: cache
        if( allocated(cache%inds)        ) deallocate(cache%inds)
        if( allocated(cache%scores)      ) deallocate(cache%scores)
        if( allocated(cache%score_fit)   ) deallocate(cache%score_fit)
        if( allocated(cache%hard_reject) ) deallocate(cache%hard_reject)
        if( allocated(cache%labels_fit)  ) deallocate(cache%labels_fit)
        if( allocated(cache%medoids_fit) ) deallocate(cache%medoids_fit)
        if( allocated(cache%score_fit_sorted) ) deallocate(cache%score_fit_sorted)
        cache%ncls            = 0
        cache%nfit            = 0
        cache%decision_mode   = 0
        cache%forced_reason   = 0
        cache%good_fit_label  = 0
        cache%bad_fit_label   = 0
        cache%separation      = 0.0
        cache%raw_threshold   = 0.0
        cache%otsu_threshold  = 0.0
        cache%otsu_separation = 0.0
        cache%otsu_ok         = .false.
    end subroutine kill_classify_cache

    ! Spec-dependent portion of the classification pipeline. Operates on
    ! a populated cache and a fully configured model and produces the
    ! same quality result that apply_linear_boundary used to produce.
    subroutine apply_cached_decision_to_quality( cache, model, quality )
        type(cavg_quality_classify_cache), intent(in)    :: cache
        class(cavg_quality_model),         intent(in)    :: model
        type(cavg_quality_result),         intent(inout) :: quality
        integer :: ncls, nfit, i, k
        type(cavg_quality_cached_decision) :: decision
        ncls = cache%ncls
        nfit = cache%nfit
        if( .not. allocated(cache%hard_reject) ) THROW_HARD('apply_cached_decision_to_quality: missing cache hard-reject mask')
        if( allocated(quality%hard_reject) )then
            if( size(quality%hard_reject) /= ncls ) &
                THROW_HARD('apply_cached_decision_to_quality: quality/cache mask size mismatch')
            if( any(quality%hard_reject .neqv. cache%hard_reject) ) &
                THROW_HARD('apply_cached_decision_to_quality: quality/cache mask mismatch')
        else
            allocate(quality%hard_reject(ncls), source=cache%hard_reject)
        endif
        if( allocated(quality%states)  ) deallocate(quality%states)
        if( allocated(quality%labels)  ) deallocate(quality%labels)
        if( allocated(quality%medoids) ) deallocate(quality%medoids)
        if( allocated(quality%scores)  ) deallocate(quality%scores)
        allocate(quality%states(ncls), quality%labels(ncls), source=0)
        allocate(quality%scores(ncls), source=0.0)
        quality%scores = cache%scores
        quality%threshold        = 0.0
        quality%raw_threshold    = 0.0
        quality%threshold_offset = 0.0
        quality%separation       = 0.0
        quality%nclust           = 0
        quality%good_label       = 0
        quality%used_threshold   = .false.
        quality%model_name       = model%name
        quality%model_context    = model%context
        quality%soft_decision    = 'hard_only'
        quality%soft_reason      = 'initial'
        call prepare_cached_decision(cache, model, decision)
        select case(cache%decision_mode)
        case(CACHE_DECISION_EMPTY)
            quality%soft_reason = decision%soft_reason
            return
        case(CACHE_DECISION_FORCED)
            call accept_fit_as_single_cluster(quality, cache%inds, decision%soft_reason)
            return
        case(CACHE_DECISION_CLUSTERED)
            continue
        case default
            THROW_HARD('apply_cached_decision_to_quality: invalid cache decision mode')
        end select
        quality%nclust        = 2
        quality%separation    = cache%separation
        quality%raw_threshold = cache%raw_threshold
        quality%labels(cache%inds) = cache%labels_fit
        allocate(quality%medoids(size(cache%medoids_fit)))
        do k = 1, size(cache%medoids_fit)
            quality%medoids(k) = cache%inds(cache%medoids_fit(k))
        end do
        if( decision%single_cluster )then
            call accept_fit_as_single_cluster(quality, cache%inds, decision%soft_reason)
            return
        endif
        quality%threshold        = decision%threshold
        quality%raw_threshold    = decision%raw_threshold
        quality%threshold_offset = decision%threshold_offset
        quality%good_label       = cache%good_fit_label
        quality%used_threshold   = decision%used_threshold
        quality%soft_decision    = decision%soft_decision
        quality%soft_reason      = decision%soft_reason
        do i = 1, nfit
            if( cached_score_fit_is_good(cache, i, decision) ) quality%states(cache%inds(i)) = 1
        end do
    end subroutine apply_cached_decision_to_quality

    subroutine cached_decision_confusion( cache, model, manual_states, tp, fp, tn, fn )
        type(cavg_quality_classify_cache), intent(in) :: cache
        class(cavg_quality_model),         intent(in) :: model
        integer,                           intent(in) :: manual_states(:)
        integer,                           intent(out):: tp, fp, tn, fn
        integer :: i, idx
        logical :: pred
        type(cavg_quality_cached_decision) :: decision
        if( size(manual_states) /= cache%ncls ) &
            THROW_HARD('cached_decision_confusion: manual-state/cache size mismatch')
        tp = 0; fp = 0; tn = 0; fn = 0
        select case(cache%decision_mode)
        case(CACHE_DECISION_EMPTY)
            return
        case(CACHE_DECISION_FORCED)
            call accumulate_accept_all_confusion(cache, manual_states, tp, fp, tn, fn)
            return
        case(CACHE_DECISION_CLUSTERED)
            continue
        case default
            THROW_HARD('cached_decision_confusion: invalid cache decision mode')
        end select
        call prepare_cached_decision(cache, model, decision)
        if( decision%single_cluster )then
            call accumulate_accept_all_confusion(cache, manual_states, tp, fp, tn, fn)
            return
        endif
        do i = 1, cache%nfit
            idx = cache%inds(i)
            pred = cached_score_fit_is_good(cache, i, decision)
            if( pred )then
                if( manual_states(idx) > 0 )then
                    tp = tp + 1
                else
                    fp = fp + 1
                endif
            else
                if( manual_states(idx) > 0 )then
                    fn = fn + 1
                else
                    tn = tn + 1
                endif
            endif
        end do
    end subroutine cached_decision_confusion

    subroutine prepare_cached_decision( cache, model, decision )
        type(cavg_quality_classify_cache), intent(in)  :: cache
        class(cavg_quality_model),         intent(in)  :: model
        type(cavg_quality_cached_decision), intent(out) :: decision
        integer :: i, min_accept, naccepted
        real    :: candidate_threshold
        decision = cavg_quality_cached_decision()
        select case(cache%decision_mode)
        case(CACHE_DECISION_EMPTY)
            decision%soft_reason = 'no_trainable_after_hard'
            return
        case(CACHE_DECISION_FORCED)
            decision%single_cluster = .true.
            decision%soft_reason    = forced_reason_name(cache%forced_reason)
            return
        case(CACHE_DECISION_CLUSTERED)
            continue
        case default
            THROW_HARD('prepare_cached_decision: invalid cache decision mode')
        end select
        if( cache%separation < model%min_score_separation )then
            if( model%use_lowsep_otsu .and. cache%otsu_ok .and. &
                cache%otsu_separation >= model%min_score_separation )then
                decision%threshold        = cache%otsu_threshold
                decision%raw_threshold    = cache%otsu_threshold
                decision%threshold_offset = 0.0
                decision%used_threshold   = .true.
                decision%soft_decision    = 'soft_threshold'
                decision%soft_reason      = 'lowsep_otsu'
                call mark_threshold_accepts_all(cache, decision)
            else
                decision%single_cluster = .true.
                decision%soft_reason    = 'low_score_separation'
            endif
        else
            candidate_threshold = cache%raw_threshold - model%boundary_margin
            decision%threshold  = candidate_threshold
            decision%raw_threshold = cache%raw_threshold
            if( model%use_otsu_window .and. cache%otsu_ok .and. &
                cache%otsu_separation >= model%min_score_separation .and. &
                cache%otsu_threshold >= cache%raw_threshold + model%otsu_min_offset .and. &
                cache%otsu_threshold <= cache%raw_threshold + model%otsu_max_offset ) decision%threshold = cache%otsu_threshold
            decision%use_rescue       = model%use_cluster_rescue
            decision%rescue_threshold = decision%threshold - model%cluster_rescue_margin
            if( model%enforce_min_accept_frac )then
                if( .not. allocated(cache%score_fit_sorted) ) &
                    THROW_HARD('prepare_cached_decision: missing sorted scores')
                min_accept = min(cache%nfit, max(1, ceiling(model%min_accept_frac * real(cache%nfit))))
                naccepted  = 0
                do i = 1, cache%nfit
                    if( cached_score_fit_is_good(cache, i, decision) ) naccepted = naccepted + 1
                end do
                if( naccepted < min_accept )then
                    decision%floor_active    = .true.
                    decision%floor_threshold = cache%score_fit_sorted(cache%nfit - min_accept + 1)
                    decision%threshold       = min(decision%threshold, decision%floor_threshold)
                endif
            endif
            decision%threshold_offset = decision%raw_threshold - decision%threshold
            decision%used_threshold   = .true.
            decision%soft_decision    = 'soft_threshold'
            if( abs(decision%threshold - cache%otsu_threshold) <= EPS .and. &
                abs(decision%threshold - candidate_threshold) > EPS )then
                decision%soft_reason = 'otsu_window'
            else
                decision%soft_reason = 'cluster_boundary'
            endif
            call mark_threshold_accepts_all(cache, decision)
        endif
    end subroutine prepare_cached_decision

    subroutine mark_threshold_accepts_all( cache, decision )
        type(cavg_quality_classify_cache), intent(in)    :: cache
        type(cavg_quality_cached_decision), intent(inout) :: decision
        integer :: i
        if( .not. decision%used_threshold ) return
        do i = 1, cache%nfit
            if( .not. cached_score_fit_is_good(cache, i, decision) ) return
        end do
        decision%soft_decision = 'hard_only'
        decision%soft_reason   = 'soft_threshold_accepts_all'
    end subroutine mark_threshold_accepts_all

    logical function cached_score_fit_is_good( cache, ifit, decision )
        type(cavg_quality_classify_cache), intent(in) :: cache
        integer,                           intent(in) :: ifit
        type(cavg_quality_cached_decision), intent(in) :: decision
        cached_score_fit_is_good = cache%score_fit(ifit) >= decision%threshold
        if( decision%use_rescue .and. cache%labels_fit(ifit) == cache%good_fit_label .and. &
            cache%score_fit(ifit) >= decision%rescue_threshold ) cached_score_fit_is_good = .true.
        if( decision%floor_active .and. cache%score_fit(ifit) >= decision%floor_threshold ) &
            cached_score_fit_is_good = .true.
    end function cached_score_fit_is_good

    subroutine accumulate_accept_all_confusion( cache, manual_states, tp, fp, tn, fn )
        type(cavg_quality_classify_cache), intent(in)    :: cache
        integer,                           intent(in)    :: manual_states(:)
        integer,                           intent(inout) :: tp, fp, tn, fn
        integer :: i, idx
        do i = 1, cache%nfit
            idx = cache%inds(i)
            if( manual_states(idx) > 0 )then
                tp = tp + 1
            else
                fp = fp + 1
            endif
        end do
    end subroutine accumulate_accept_all_confusion

    function forced_reason_name( reason ) result( name )
        integer, intent(in) :: reason
        character(len=32)   :: name
        select case(reason)
        case(CACHE_FORCED_TOO_FEW)
            name = 'too_few_trainable'
        case(CACHE_FORCED_FLAT_DMAT)
            name = 'flat_feature_distances'
        case(CACHE_FORCED_INVALID_TWO_CLUSTER)
            name = 'invalid_two_cluster_result'
        case default
            name = 'single_cluster'
        end select
    end function forced_reason_name

    real function feature_vector_distance( feat1, feat2, weights )
        real, intent(in) :: feat1(:), feat2(:), weights(:)
        integer :: ifeat
        real    :: delta
        feature_vector_distance = 0.0
        do ifeat = 1, CAVG_QUALITY_NFEATS
            ! Nonzero weights define both the linear score and the feature set
            ! allowed to shape k-medoids. Distances are intentionally masked,
            ! not weight-scaled: sqrt(sum_{w_i>EPS} (a_i - b_i)^2).
            ! Hard rejections are applied before clustering.
            if( weights(ifeat) <= EPS ) cycle
            delta = feat1(ifeat) - feat2(ifeat)
            feature_vector_distance = feature_vector_distance + delta**2
        end do
        feature_vector_distance = sqrt(feature_vector_distance)
    end function feature_vector_distance

    subroutine accept_fit_as_single_cluster( quality, inds, reason )
        type(cavg_quality_result), intent(inout) :: quality
        integer,                   intent(in)    :: inds(:)
        character(len=*), optional,intent(in)    :: reason
        if( size(inds) == 0 )then
            if( allocated(quality%medoids) ) deallocate(quality%medoids)
            quality%threshold        = 0.0
            quality%raw_threshold    = 0.0
            quality%threshold_offset = 0.0
            quality%nclust           = 0
            quality%good_label       = 0
            quality%used_threshold   = .false.
            quality%soft_decision    = 'hard_only'
            quality%soft_reason      = 'no_trainable_after_hard'
            return
        endif
        quality%states(inds) = 1
        quality%labels(inds) = 1
        if( allocated(quality%medoids) ) deallocate(quality%medoids)
        allocate(quality%medoids(1), source=inds(1))
        quality%threshold        = minval(quality%scores(inds)) - EPS
        quality%raw_threshold    = quality%threshold
        quality%threshold_offset = 0.0
        quality%nclust           = 1
        quality%good_label       = 1
        quality%used_threshold   = .false.
        quality%soft_decision    = 'hard_only'
        if( present(reason) )then
            quality%soft_reason = reason
        else
            quality%soft_reason = 'single_cluster'
        endif
    end subroutine accept_fit_as_single_cluster

    logical function two_cluster_result_is_valid( nclust, labels, medoids, nfit )
        integer,              intent(in) :: nclust, nfit
        integer, allocatable, intent(in) :: labels(:), medoids(:)
        two_cluster_result_is_valid = .false.
        if( nclust /= 2 ) return
        if( .not. allocated(labels)  ) return
        if( .not. allocated(medoids) ) return
        if( size(labels) /= nfit ) return
        if( size(medoids) < 2 ) return
        if( count(labels == 1) == 0 .or. count(labels == 2) == 0 ) return
        two_cluster_result_is_valid = .true.
    end function two_cluster_result_is_valid

    subroutine choose_tied_good_label( scores, labels, medoids, good_label, bad_label )
        real,    intent(in)  :: scores(:)
        integer, intent(in)  :: labels(:), medoids(:)
        integer, intent(out) :: good_label, bad_label
        integer :: n1, n2
        n1 = count(labels == 1)
        n2 = count(labels == 2)
        if( n1 > n2 )then
            good_label = 1
            bad_label  = 2
        else if( n2 > n1 )then
            good_label = 2
            bad_label  = 1
        else if( scores(medoids(1)) >= scores(medoids(2)) )then
            good_label = 1
            bad_label  = 2
        else
            good_label = 2
            bad_label  = 1
        endif
    end subroutine choose_tied_good_label

    subroutine assert_valid_model_context( context )
        character(len=*), intent(in) :: context
        select case(trim(context))
            case('chunk', 'pool')
                return
            case default
                THROW_HARD('invalid class-average quality model context: '//trim(context))
        end select
    end subroutine assert_valid_model_context

    subroutine otsu_score_threshold( scores, threshold, separation, ok )
        real,    intent(in)  :: scores(:)
        real,    intent(out) :: threshold, separation
        logical, intent(out) :: ok
        real, allocatable :: sorted(:), prefix(:)
        integer :: i, n, nlo, nhi
        real    :: candidate, mean_lo, mean_hi, between, best_between, total
        n = size(scores)
        threshold  = minval(scores) - EPS
        separation = 0.0
        ok         = .false.
        if( n < 4 ) return
        allocate(sorted(n), prefix(n))
        sorted = scores
        call hpsort(sorted)
        if( sorted(n) - sorted(1) <= EPS )then
            deallocate(sorted, prefix)
            return
        endif
        prefix(1) = sorted(1)
        do i = 2, n
            prefix(i) = prefix(i-1) + sorted(i)
        end do
        total = prefix(n)
        best_between = -huge(1.0)
        do i = 1, n - 1
            if( sorted(i+1) - sorted(i) <= EPS ) cycle
            nlo = i
            nhi = n - i
            mean_lo = prefix(i) / real(nlo)
            mean_hi = (total - prefix(i)) / real(nhi)
            between = real(nlo) * real(nhi) * (mean_hi - mean_lo)**2
            if( between > best_between )then
                best_between = between
                candidate    = 0.5 * (sorted(i) + sorted(i+1))
                threshold    = candidate
                separation   = mean_hi - mean_lo
                ok           = .true.
            endif
        end do
        deallocate(sorted, prefix)
    end subroutine otsu_score_threshold

    real function mean_score_for_label( scores, labels, label )
        real,    intent(in) :: scores(:)
        integer, intent(in) :: labels(:), label
        integer             :: n
        n = count(labels == label)
        ! Callers validate the two-cluster result before asking for label
        ! means; the guard keeps this helper safe for future direct use.
        if( n == 0 ) THROW_HARD('mean_score_for_label: empty cluster')
        mean_score_for_label = sum(scores, mask=labels == label) / real(n)
    end function mean_score_for_label

end module simple_cavg_quality_model
