!@descr: instantiable class-average quality decision model
module simple_cavg_quality_model
use simple_defs,               only: LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string_utils,       only: str_is_true, lowercase, uppercase, &
    fortran_symbol_from_string, fortran_quote, fortran_logical
use simple_clustering_utils,   only: cluster_dmat
use simple_srch_sort_loc,      only: hpsort
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, CAVG_QUALITY_MAX_INTERACTIONS, EPS, CLIP_Z, &
    CAVG_MODEL_FAMILY_LINEAR, CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC, CAVG_QUALITY_CONTEXT_CHUNK, &
    CAVG_QUALITY_CONTEXT_POOL, cavg_quality_model_spec, cavg_quality_result
use simple_cavg_quality_feats, only: cavg_quality_feature_name
use simple_cavg_quality_stats, only: normalize_quality_dmat
implicit none
private
#include "simple_local_flags.inc"

public :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT
public :: CAVG_QUALITY_MODEL_CHUNK_LINEAR
public :: CAVG_QUALITY_MODEL_POOL
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
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT = 'chunk100mics'
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_LINEAR = 'chunk100mics_linear'
character(len=*), parameter :: CAVG_QUALITY_MODEL_POOL = 'pool'
character(len=*), parameter :: BUILTIN_MODEL_NAMES = CAVG_QUALITY_MODEL_CHUNK_DEFAULT//'|'//&
    CAVG_QUALITY_MODEL_CHUNK_LINEAR//'|'//CAVG_QUALITY_MODEL_POOL

real, parameter :: CLUSTER_RESCUE_MARGIN = 0.20

! Default chunk class-average quality model, promoted from the pairwise
! logistic artifact learned from
! /Users/elmlundho/cavgs_quality/chunk100mic_training_data_v4.
character(len=*), parameter :: CHUNK100MICS_FEATURE_POLICY = 'microchunk_plus_score_signal'
real, parameter :: CAVG_QUALITY_LOGISTIC_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02 ]
character(len=*), parameter :: CHUNK100MICS_LINEAR_FEATURE_POLICY = 'microchunk_plus_score_signal'
real, parameter :: CAVG_QUALITY_CHUNK100MICS_LINEAR_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    9.978756E-02, 1.167914E-01, 3.642511E-02, 1.329548E-01, &
    1.402481E-01, 1.610645E-01, 6.630784E-02, 6.981037E-02, &
    1.294257E-01, 0.000000E+00, 0.000000E+00, 4.718454E-02, &
    0.000000E+00, 0.000000E+00 ]
character(len=*), parameter :: POOL_FEATURE_POLICY = 'microchunk_plus_score_signal'
real, parameter :: CAVG_QUALITY_POOL_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02 ]
real, parameter :: CHUNK100MICS_BOUNDARY_MARGIN      = 0.00
real, parameter :: CHUNK100MICS_MIN_SCORE_SEPARATION = 0.05
real, parameter :: CHUNK100MICS_OTSU_MIN_OFFSET      = 0.00
real, parameter :: CHUNK100MICS_OTSU_MAX_OFFSET      = 0.00
real, parameter :: CHUNK100MICS_MIN_ACCEPT_FRAC      = 0.00
real, parameter :: CHUNK100MICS_LINEAR_BOUNDARY_MARGIN = 0.30
real, parameter :: CHUNK100MICS_LINEAR_OTSU_MIN_OFFSET = 0.05
real, parameter :: CHUNK100MICS_LINEAR_OTSU_MAX_OFFSET = 0.40
real, parameter :: CHUNK100MICS_LINEAR_MIN_ACCEPT_FRAC = 0.60

type :: cavg_quality_model
    character(len=64) :: name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
    character(len=32) :: context                 = 'chunk'
    character(len=64) :: feature_policy          = CHUNK100MICS_FEATURE_POLICY
    character(len=32) :: model_family            = CAVG_MODEL_FAMILY_LINEAR
    real              :: weights(CAVG_QUALITY_NFEATS) = CAVG_QUALITY_LOGISTIC_WEIGHTS
    real              :: intercept               = 0.0
    real              :: linear_coefficients(CAVG_QUALITY_NFEATS) = 0.0
    integer           :: n_interactions          = 0
    integer           :: interaction_terms(CAVG_QUALITY_MAX_INTERACTIONS,2) = 0
    real              :: interaction_coefficients(CAVG_QUALITY_MAX_INTERACTIONS) = 0.0
    real              :: prob_threshold          = 0.5
    real              :: regularization_lambda   = 0.0
    real              :: calibration_temperature = 1.0
    real              :: boundary_margin         = CHUNK100MICS_BOUNDARY_MARGIN
    real              :: min_score_separation    = CHUNK100MICS_MIN_SCORE_SEPARATION
    real              :: otsu_min_offset         = CHUNK100MICS_OTSU_MIN_OFFSET
    real              :: otsu_max_offset         = CHUNK100MICS_OTSU_MAX_OFFSET
    real              :: cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
    real              :: min_accept_frac         = CHUNK100MICS_MIN_ACCEPT_FRAC
    logical           :: use_lowsep_otsu         = .false.
    logical           :: use_otsu_window         = .false.
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
                spec = chunk100mics_model_spec()
            case(CAVG_QUALITY_MODEL_CHUNK_LINEAR)
                spec = chunk100mics_linear_model_spec()
            case(CAVG_QUALITY_MODEL_POOL)
                spec = pool_model_spec()
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
        self%context                 = trim(spec%context)
        self%feature_policy          = trim(spec%feature_policy)
        self%model_family            = trim(spec%model_family)
        self%weights                 = spec%weights
        self%intercept               = spec%intercept
        self%linear_coefficients     = spec%linear_coefficients
        self%n_interactions          = spec%n_interactions
        self%interaction_terms       = spec%interaction_terms
        self%interaction_coefficients = spec%interaction_coefficients
        self%prob_threshold          = spec%prob_threshold
        self%regularization_lambda   = spec%regularization_lambda
        self%calibration_temperature = spec%calibration_temperature
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
        spec%model_family            = self%model_family
        spec%weights                 = self%weights
        spec%intercept               = self%intercept
        spec%linear_coefficients     = self%linear_coefficients
        spec%n_interactions          = self%n_interactions
        spec%interaction_terms       = self%interaction_terms
        spec%interaction_coefficients = self%interaction_coefficients
        spec%prob_threshold          = self%prob_threshold
        spec%regularization_lambda   = self%regularization_lambda
        spec%calibration_temperature = self%calibration_temperature
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

    function chunk100mics_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
        spec%context                 = 'chunk'
        spec%feature_policy          = CHUNK100MICS_FEATURE_POLICY
        spec%model_family            = CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC
        spec%weights                 = CAVG_QUALITY_LOGISTIC_WEIGHTS
        spec%intercept               = 2.681497E+00
        spec%linear_coefficients     = [ &
           -6.286113E-01,  6.247919E-01, -2.818124E-01,  1.076405E+00, &
           -7.464533E-01,  6.145704E-01, -6.273948E-02,  1.407827E-01, &
            2.027739E-01, -4.221159E-01,  7.574302E-01, -2.793084E-01, &
            2.794536E+00,  7.073164E-01 ]
        call set_pairwise_interactions_for_feature_count(spec, 14, [ &
           -5.405482E-01,  2.232487E-02,  6.679537E-01, -9.844879E-01, &
            4.016007E-01,  8.536045E-01,  6.862369E-01, -3.145630E-01, &
           -6.455917E-01,  8.791583E-01, -4.676543E-01,  4.777114E-01, &
            5.521021E-01,  1.271381E-02, -1.254850E-01,  6.664535E-02, &
           -2.221981E-01, -6.056051E-02, -4.815833E-01,  1.268961E-01, &
           -1.271788E-01,  2.494773E-01, -1.269121E-01,  1.350642E+00, &
           -1.692577E-01, -9.738522E-02,  1.265601E-01, -5.409655E-01, &
            2.652861E-01,  4.474671E-01, -2.019620E-01, -4.014040E-01, &
           -3.299867E-01,  3.392396E-01, -2.475521E-01, -2.321465E-01, &
           -2.825888E-01, -3.237399E-01,  8.590719E-01, -1.205967E-01, &
           -2.884976E-01, -4.536750E-01,  2.778712E-01, -7.639435E-01, &
           -2.127371E-01, -4.576955E-01,  2.047243E-01,  7.392398E-02, &
            6.617538E-01, -2.564861E-01, -2.180120E-02,  1.481527E-02, &
            2.854128E-01, -1.482570E-01,  3.309094E-01, -8.375521E-01, &
           -6.958253E-01, -1.989315E-02, -1.469698E-01,  8.772197E-02, &
            2.987876E-02, -4.436344E-01,  2.778504E-01,  1.992271E-01, &
           -5.480839E-01, -5.014452E-01, -1.514034E-01, -7.361591E-01, &
           -6.602764E-01, -1.696232E-01, -1.325360E-01,  9.232772E-02, &
           -5.104917E-01, -3.444873E-01, -8.808707E-01,  7.508239E-02, &
            3.759776E-02,  3.667796E-01,  6.808946E-01,  7.423107E-01, &
           -3.577121E-01, -6.254713E-03,  9.768074E-02, -5.226645E-01, &
            1.452297E-01, -6.939387E-01,  1.408253E-01, -3.367719E-01, &
           -1.496779E-01,  2.686394E-01,  2.973516E-02 ])
        spec%prob_threshold          = 3.500000E-01
        spec%regularization_lambda   = 1.000000E-03
        spec%calibration_temperature = 1.000000E+00
        spec%boundary_margin         = CHUNK100MICS_BOUNDARY_MARGIN
        spec%min_score_separation    = CHUNK100MICS_MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = CHUNK100MICS_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK100MICS_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = CHUNK100MICS_MIN_ACCEPT_FRAC
        spec%use_lowsep_otsu         = .false.
        spec%use_otsu_window         = .false.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function chunk100mics_model_spec

    function chunk100mics_linear_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_CHUNK_LINEAR
        spec%context                 = 'chunk'
        spec%feature_policy          = CHUNK100MICS_LINEAR_FEATURE_POLICY
        spec%model_family            = CAVG_MODEL_FAMILY_LINEAR
        spec%weights                 = CAVG_QUALITY_CHUNK100MICS_LINEAR_WEIGHTS
        spec%boundary_margin         = CHUNK100MICS_LINEAR_BOUNDARY_MARGIN
        spec%min_score_separation    = CHUNK100MICS_MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = CHUNK100MICS_LINEAR_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK100MICS_LINEAR_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = CHUNK100MICS_LINEAR_MIN_ACCEPT_FRAC
        spec%use_lowsep_otsu         = .false.
        spec%use_otsu_window         = .true.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .true.
    end function chunk100mics_linear_model_spec

    function pool_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_POOL
        spec%context                 = 'pool'
        spec%feature_policy          = POOL_FEATURE_POLICY
        spec%model_family            = CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC
        spec%weights                 = CAVG_QUALITY_POOL_WEIGHTS
        spec%intercept               = 2.893284E+00
        spec%linear_coefficients     = [ &
            1.170739E+00,  2.270497E+00, -1.045277E-01,  3.912764E-01, &
            2.637467E-01, -3.482836E-01,  5.717540E-01,  5.872802E-01, &
           -8.924436E-01, -3.850991E-01, -2.258097E-01, -8.288236E-01, &
            2.961799E-01,  4.474984E-02 ]
        call set_pairwise_interactions_for_feature_count(spec, 14, [ &
            8.842925E-01, -2.361528E-01,  1.482243E-01,  1.150898E-01, &
           -5.629415E-01, -5.466182E-01,  4.866918E-01, -3.237223E-01, &
           -1.537362E-01, -1.599723E-01, -6.046203E-02,  2.468477E-01, &
           -1.763398E-01,  1.439890E-01, -2.247636E-01, -1.460124E-01, &
            4.241256E-01,  1.989624E-01, -5.949336E-01, -2.482053E-01, &
            2.171634E-01,  4.548043E-02, -4.890917E-01,  2.822063E-01, &
           -6.085361E-02, -1.576771E-01,  9.350962E-02, -8.126303E-02, &
            5.275285E-02,  1.138827E-01,  2.839310E-01,  1.540730E-01, &
           -7.069448E-02,  3.623002E-01, -1.834264E-01, -2.364621E-02, &
           -3.576021E-01,  5.002445E-01,  2.744766E-01, -3.280614E-01, &
            5.400299E-01, -3.220865E-01,  2.986487E-01, -2.170379E-01, &
            2.146203E-01, -6.204020E-01,  7.821966E-02, -1.992914E-02, &
            4.002428E-01,  1.279922E-01,  3.450857E-01, -1.679638E-01, &
            6.394580E-02,  5.042716E-02, -2.082259E-01, -3.374619E-01, &
            5.608899E-01,  1.005668E+00, -5.060473E-01, -1.200944E-01, &
            4.492203E-01, -9.182516E-01, -1.505287E-01, -1.230667E-01, &
            7.545972E-02, -2.856418E-01, -1.219734E-01, -4.482515E-01, &
           -3.432786E-01,  4.062396E-01, -8.293809E-02,  3.388880E-01, &
           -1.097642E-01,  2.624200E-01, -7.789035E-01,  2.066645E-01, &
           -5.378622E-01,  1.340223E-02, -6.153102E-02,  1.182417E+00, &
            1.317468E-01, -2.861343E-01,  2.079481E-01, -2.212095E-01, &
            6.200773E-01, -7.408983E-03, -1.736890E-01,  3.201505E-01, &
           -2.867308E-01,  8.472726E-02, -1.584640E+00 ])
        spec%prob_threshold          = 3.500000E-01
        spec%regularization_lambda   = 3.000000E-04
        spec%calibration_temperature = 1.000000E+00
        spec%boundary_margin         = 0.0
        spec%min_score_separation    = CHUNK100MICS_MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = 0.0
        spec%otsu_max_offset         = 0.0
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = 0.0
        spec%use_lowsep_otsu         = .false.
        spec%use_otsu_window         = .false.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function pool_model_spec

    subroutine set_pairwise_interactions_for_feature_count( spec, nfeatures, coefficients )
        type(cavg_quality_model_spec), intent(inout) :: spec
        integer,                       intent(in)    :: nfeatures
        real,                          intent(in)    :: coefficients(:)
        integer :: ifeat, jfeat, iterm, expected_terms
        if( nfeatures < 1 .or. nfeatures > CAVG_QUALITY_NFEATS ) &
            THROW_HARD('set_pairwise_interactions_for_feature_count: invalid feature count')
        expected_terms = (nfeatures * (nfeatures - 1)) / 2
        if( size(coefficients) /= expected_terms ) &
            THROW_HARD('set_pairwise_interactions_for_feature_count: coefficient count mismatch')
        spec%n_interactions          = expected_terms
        spec%interaction_terms       = 0
        spec%interaction_coefficients = 0.0
        iterm = 0
        do ifeat = 1, nfeatures - 1
            do jfeat = ifeat + 1, nfeatures
                iterm = iterm + 1
                spec%interaction_terms(iterm,:) = [ifeat, jfeat]
            end do
        end do
        spec%interaction_coefficients(1:expected_terms) = coefficients
    end subroutine set_pairwise_interactions_for_feature_count

    subroutine set_pairwise_interactions_for_feature_indices( spec, feature_indices, coefficients )
        type(cavg_quality_model_spec), intent(inout) :: spec
        integer,                       intent(in)    :: feature_indices(:)
        real,                          intent(in)    :: coefficients(:)
        integer :: i, j, iterm, nfeatures, expected_terms
        nfeatures = size(feature_indices)
        if( nfeatures < 1 .or. nfeatures > CAVG_QUALITY_NFEATS ) &
            THROW_HARD('set_pairwise_interactions_for_feature_indices: invalid feature count')
        do i = 1, nfeatures
            if( feature_indices(i) < 1 .or. feature_indices(i) > CAVG_QUALITY_NFEATS ) &
                THROW_HARD('set_pairwise_interactions_for_feature_indices: invalid feature index')
        end do
        expected_terms = (nfeatures * (nfeatures - 1)) / 2
        if( size(coefficients) /= expected_terms ) &
            THROW_HARD('set_pairwise_interactions_for_feature_indices: coefficient count mismatch')
        spec%n_interactions          = expected_terms
        spec%interaction_terms       = 0
        spec%interaction_coefficients = 0.0
        iterm = 0
        do i = 1, nfeatures - 1
            do j = i + 1, nfeatures
                iterm = iterm + 1
                spec%interaction_terms(iterm,:) = [feature_indices(i), feature_indices(j)]
            end do
        end do
        spec%interaction_coefficients(1:expected_terms) = coefficients
    end subroutine set_pairwise_interactions_for_feature_indices

    subroutine normalize( self )
        class(cavg_quality_model), intent(inout) :: self
        select case(trim(self%model_family))
        case(CAVG_MODEL_FAMILY_LINEAR)
            self%weights = max(0.0, self%weights)
            if( sum(self%weights) > EPS )then
                self%weights = self%weights / sum(self%weights)
            else
                self%weights = 1.0 / real(CAVG_QUALITY_NFEATS)
                self%weights = self%weights / sum(self%weights)
            endif
        case(CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC)
            self%prob_threshold          = min(1.0, max(0.0, self%prob_threshold))
            self%calibration_temperature = max(EPS, self%calibration_temperature)
            if( self%n_interactions < 0 .or. self%n_interactions > CAVG_QUALITY_MAX_INTERACTIONS ) &
                THROW_HARD('normalize: invalid n_interactions for pairwise logistic model')
        case default
            THROW_HARD('normalize: unknown model_family: '//trim(self%model_family))
        end select
    end subroutine normalize

    subroutine classify( self, quality )
        class(cavg_quality_model), intent(in)    :: self
        type(cavg_quality_result), intent(inout) :: quality
        if( .not. allocated(quality%features)    ) THROW_HARD('classify: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('classify: missing hard-reject mask')
        quality%model_name     = self%name
        select case(trim(self%model_family))
        case(CAVG_MODEL_FAMILY_LINEAR)
            call apply_linear_boundary(quality, self)
        case(CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC)
            call apply_pairwise_logistic(quality, self)
        case default
            THROW_HARD('classify: unknown model_family: '//trim(self%model_family))
        end select
    end subroutine classify

    subroutine write_model( self, fname )
        class(cavg_quality_model), intent(in) :: self
        character(len=*),          intent(in) :: fname
        integer :: funit, i
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# model_cavgs_rejection model'
        if( trim(self%model_family) == CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC )then
            write(funit,'(A)') 'model_version=9'
        else
            write(funit,'(A)') 'model_version=8'
        endif
        write(funit,'(A,A)') 'name=', trim(self%name)
        write(funit,'(A,A)') 'model_family=', trim(self%model_family)
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
        if( trim(self%model_family) == CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC )then
            write(funit,'(A,ES14.6)') 'intercept=', self%intercept
            call write_model_real_list(funit, 'linear_coefficients=', self%linear_coefficients)
            call write_interaction_terms(funit, self)
            call write_model_real_list(funit, 'interaction_coefficients=', self%interaction_coefficients, self%n_interactions)
            write(funit,'(A,ES14.6)') 'prob_threshold=', self%prob_threshold
            write(funit,'(A,ES14.6)') 'regularization_lambda=', self%regularization_lambda
            write(funit,'(A,ES14.6)') 'calibration_temperature=', self%calibration_temperature
        endif
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
        write(funit,'(A,A)') '        spec%model_family            = ', trim(fortran_quote(model%model_family))
        call write_weights_assignment(funit, model%weights)
        if( trim(model%model_family) == CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC ) &
            call write_logistic_spec_assignments(funit, model)
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

    subroutine write_logistic_spec_assignments( funit, model )
        integer,                  intent(in) :: funit
        type(cavg_quality_model), intent(in) :: model
        integer :: iterm
        write(funit,'(A,ES14.6)') '        spec%intercept               = ', model%intercept
        call write_real_array_assignment(funit, '        spec%linear_coefficients', &
            model%linear_coefficients, CAVG_QUALITY_NFEATS)
        write(funit,'(A,I0)') '        spec%n_interactions          = ', model%n_interactions
        write(funit,'(A)') '        spec%interaction_terms        = 0'
        write(funit,'(A)') '        spec%interaction_coefficients = 0.0'
        do iterm = 1, model%n_interactions
            write(funit,'(A,I0,A,I0,A,I0,A)') '        spec%interaction_terms(', iterm, ',:) = [', &
                model%interaction_terms(iterm,1), ', ', model%interaction_terms(iterm,2), ' ]'
            write(funit,'(A,I0,A,ES14.6)') '        spec%interaction_coefficients(', iterm, ') = ', &
                model%interaction_coefficients(iterm)
        end do
        write(funit,'(A,ES14.6)') '        spec%prob_threshold          = ', model%prob_threshold
        write(funit,'(A,ES14.6)') '        spec%regularization_lambda   = ', model%regularization_lambda
        write(funit,'(A,ES14.6)') '        spec%calibration_temperature = ', model%calibration_temperature
    end subroutine write_logistic_spec_assignments

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

    subroutine write_real_array_assignment( funit, lhs, values, nvals )
        integer,          intent(in) :: funit, nvals
        character(len=*), intent(in) :: lhs
        real,             intent(in) :: values(:)
        integer :: i
        if( nvals < 0 .or. nvals > size(values) ) THROW_HARD('write_real_array_assignment: invalid value count')
        write(funit,'(A,A)') trim(lhs), ' = [ &'
        write(funit,'(A)', advance='no') '            '
        do i = 1, nvals
            write(funit,'(ES14.6)', advance='no') values(i)
            if( i < nvals ) write(funit,'(A)', advance='no') ', '
            if( mod(i, 4) == 0 .and. i < nvals )then
                write(funit,'(A)') '&'
                write(funit,'(A)', advance='no') '            '
            endif
        end do
        write(funit,'(A)') ' ]'
    end subroutine write_real_array_assignment

    subroutine write_model_real_list( funit, key, vals, nvals )
        integer,          intent(in) :: funit
        character(len=*), intent(in) :: key
        real,             intent(in) :: vals(:)
        integer, optional,intent(in) :: nvals
        integer :: i, nwrite
        nwrite = size(vals)
        if( present(nvals) ) nwrite = nvals
        if( nwrite < 0 .or. nwrite > size(vals) ) THROW_HARD('write_model_real_list: invalid value count')
        write(funit,'(A)', advance='no') trim(key)
        do i = 1, nwrite
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') vals(i)
        end do
        write(funit,*)
    end subroutine write_model_real_list

    subroutine write_interaction_terms( funit, model )
        integer,                  intent(in) :: funit
        type(cavg_quality_model), intent(in) :: model
        integer :: i, ifeat, jfeat
        if( model%n_interactions < 0 .or. model%n_interactions > CAVG_QUALITY_MAX_INTERACTIONS ) &
            THROW_HARD('write_interaction_terms: invalid n_interactions')
        write(funit,'(A)', advance='no') 'interaction_terms='
        do i = 1, model%n_interactions
            ifeat = model%interaction_terms(i,1)
            jfeat = model%interaction_terms(i,2)
            if( ifeat < 1 .or. ifeat > CAVG_QUALITY_NFEATS .or. &
                jfeat < 1 .or. jfeat > CAVG_QUALITY_NFEATS ) &
                THROW_HARD('write_interaction_terms: invalid interaction feature index')
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(I0,A,I0)', advance='no') ifeat, ':', jfeat
        end do
        write(funit,*)
    end subroutine write_interaction_terms

    subroutine read_model( self, fname )
        class(cavg_quality_model), intent(inout) :: self
        character(len=*),          intent(in)    :: fname
        character(len=XLONGSTRLEN) :: line
        character(len=LONGSTRLEN)  :: key, preset_name
        character(len=XLONGSTRLEN) :: val
        integer :: funit, ios, parse_ios, model_version, n_interaction_coefficients
        logical :: have_model_family, have_preset, ok_line
        ! Model files are complete model definitions. Start from chunk defaults,
        ! apply any preset found in the file, then apply explicit key overrides.
        call self%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
        open(newunit=funit, file=trim(fname), status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('read_model: failed to open '//trim(fname))
        have_model_family = .false.
        have_preset       = .false.
        model_version     = 0
        preset_name       = ''
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            call parse_model_key_value(line, key, val, ok_line)
            if( .not. ok_line ) cycle
            select case(trim(key))
            case('model_version')
                read(val,*,iostat=parse_ios) model_version
            case('model_family')
                have_model_family = .true.
            case('preset')
                preset_name = trim(val)
                have_preset = .true.
            end select
        end do
        if( have_preset ) call self%init_preset(trim(preset_name))
        if( .not. have_model_family .and. model_version < 9 )then
            self%model_family             = CAVG_MODEL_FAMILY_LINEAR
            self%intercept                = 0.0
            self%linear_coefficients      = 0.0
            self%n_interactions           = 0
            self%interaction_terms        = 0
            self%interaction_coefficients = 0.0
            self%prob_threshold           = 0.5
            self%regularization_lambda    = 0.0
            self%calibration_temperature  = 1.0
        endif
        rewind(funit)
        n_interaction_coefficients = 0
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
                case('model_family')
                    self%model_family = trim(val)
                case('context')
                    select case(trim(val))
                        case(CAVG_QUALITY_CONTEXT_CHUNK, CAVG_QUALITY_CONTEXT_POOL)
                            self%context = trim(val)
                        case DEFAULT
                            THROW_HARD('read_model: context must be chunk or pool')
                    end select
                case('feature_policy', 'feature_family_set')
                    self%feature_policy = trim(val)
                case('feature_weights')
                    call read_feature_weights(val, self%weights)
                case('intercept')
                    read(val,*,iostat=parse_ios) self%intercept
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse intercept')
                case('linear_coefficients')
                    call read_real_values_keyed(val, self%linear_coefficients, CAVG_QUALITY_NFEATS, &
                        'linear_coefficients')
                case('interaction_terms')
                    call read_interaction_terms(val, self%interaction_terms, self%n_interactions)
                case('interaction_coefficients')
                    call read_real_values_keyed(val, self%interaction_coefficients, CAVG_QUALITY_MAX_INTERACTIONS, &
                        'interaction_coefficients', n_interaction_coefficients)
                case('prob_threshold')
                    read(val,*,iostat=parse_ios) self%prob_threshold
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse prob_threshold')
                case('regularization_lambda')
                    read(val,*,iostat=parse_ios) self%regularization_lambda
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse regularization_lambda')
                case('calibration_temperature')
                    read(val,*,iostat=parse_ios) self%calibration_temperature
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse calibration_temperature')
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
        if( trim(self%model_family) == CAVG_MODEL_FAMILY_PAIRWISE_LOGISTIC )then
            if( n_interaction_coefficients /= self%n_interactions ) &
                THROW_HARD('read_model: interaction_coefficients count must match interaction_terms count')
        endif
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

    subroutine read_real_values_keyed( val, values, maxvals, key, nvals_out )
        character(len=*), intent(in)    :: val, key
        real,             intent(inout) :: values(:)
        integer,          intent(in)    :: maxvals
        integer, optional,intent(out)   :: nvals_out
        character(len=LONGSTRLEN) :: errmsg
        integer :: nvals
        values = 0.0
        call parse_real_values(val, values, nvals)
        if( nvals > maxvals .or. nvals > size(values) )then
            write(errmsg,'(A,A,A,I0,A,I0)') 'read_model: ', trim(key), ' expected at most ', maxvals, &
                ' values, got ', nvals
            THROW_HARD(trim(errmsg))
        endif
        if( present(nvals_out) ) nvals_out = nvals
    end subroutine read_real_values_keyed

    subroutine parse_real_values( val, values, nvals )
        character(len=*), intent(in)  :: val
        real,             intent(out) :: values(:)
        integer,          intent(out) :: nvals
        character(len=XLONGSTRLEN) :: work, token
        integer :: isep, ios
        values = 0.0
        nvals  = 0
        work   = adjustl(trim(val))
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
            if( nvals > size(values) ) cycle
            read(token,*,iostat=ios) values(nvals)
            if( ios /= 0 ) THROW_HARD('read_model: failed to parse real-valued list')
        end do
    end subroutine parse_real_values

    subroutine read_interaction_terms( val, terms, nterms )
        character(len=*), intent(in)  :: val
        integer,          intent(out) :: terms(CAVG_QUALITY_MAX_INTERACTIONS,2)
        integer,          intent(out) :: nterms
        character(len=XLONGSTRLEN) :: work, token
        character(len=LONGSTRLEN)  :: lhs, rhs, errmsg
        integer :: isep, icolon, ifeat, jfeat
        terms  = 0
        nterms = 0
        work   = adjustl(trim(val))
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
            token = adjustl(trim(token))
            if( len_trim(token) == 0 ) cycle
            icolon = index(token, ':')
            if( icolon <= 1 .or. icolon >= len_trim(token) ) &
                THROW_HARD('read_model: interaction_terms entries must be feature_a:feature_b')
            lhs = token(1:icolon-1)
            rhs = token(icolon+1:)
            ifeat = feature_index_from_token(lhs)
            jfeat = feature_index_from_token(rhs)
            if( ifeat < 1 .or. jfeat < 1 )then
                write(errmsg,'(A,A)') 'read_model: unknown interaction feature in ', trim(token)
                THROW_HARD(trim(errmsg))
            endif
            nterms = nterms + 1
            if( nterms > CAVG_QUALITY_MAX_INTERACTIONS ) &
                THROW_HARD('read_model: too many interaction_terms')
            terms(nterms,1) = ifeat
            terms(nterms,2) = jfeat
        end do
    end subroutine read_interaction_terms

    integer function feature_index_from_token( token )
        character(len=*), intent(in) :: token
        character(len=LONGSTRLEN) :: work
        integer :: ifeat, ios
        feature_index_from_token = 0
        work = adjustl(trim(token))
        read(work,*,iostat=ios) ifeat
        if( ios == 0 )then
            if( ifeat >= 1 .and. ifeat <= CAVG_QUALITY_NFEATS ) feature_index_from_token = ifeat
            return
        endif
        work = lowercase(trim(work))
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( trim(work) == trim(lowercase(cavg_quality_feature_name(ifeat))) )then
                feature_index_from_token = ifeat
                return
            endif
        end do
    end function feature_index_from_token

    subroutine parse_feature_weight_values( val, weights, nvals )
        character(len=*), intent(in)  :: val
        real,             intent(out) :: weights(CAVG_QUALITY_NFEATS)
        integer,          intent(out) :: nvals
        character(len=XLONGSTRLEN) :: work, token
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
        self%model_family            = ''
        self%weights                 = 0.0
        self%intercept               = 0.0
        self%linear_coefficients     = 0.0
        self%n_interactions          = 0
        self%interaction_terms       = 0
        self%interaction_coefficients = 0.0
        self%prob_threshold          = 0.0
        self%regularization_lambda   = 0.0
        self%calibration_temperature = 1.0
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

    subroutine apply_pairwise_logistic( quality, model )
        type(cavg_quality_result), intent(inout) :: quality
        class(cavg_quality_model), intent(in)    :: model
        integer :: ncls, icls
        real    :: prob
        if( .not. allocated(quality%features)    ) THROW_HARD('apply_pairwise_logistic: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('apply_pairwise_logistic: missing hard-reject mask')
        if( size(quality%features, dim=2) /= CAVG_QUALITY_NFEATS ) THROW_HARD('apply_pairwise_logistic: invalid feature count')
        ncls = size(quality%features, dim=1)
        if( size(quality%hard_reject) /= ncls ) THROW_HARD('apply_pairwise_logistic: invalid mask size')
        if( allocated(quality%states)  ) deallocate(quality%states)
        if( allocated(quality%labels)  ) deallocate(quality%labels)
        if( allocated(quality%medoids) ) deallocate(quality%medoids)
        if( allocated(quality%scores)  ) deallocate(quality%scores)
        allocate(quality%states(ncls), quality%labels(ncls), source=0)
        allocate(quality%scores(ncls), source=0.0)

        ! Logistic models are direct probability classifiers:
        !
        !   eta = intercept + sum_i beta_i z_i + sum_(i,j) gamma_ij z_i z_j
        !   P(accept) = sigmoid(eta / calibration_temperature)
        !
        ! Standard hard gates have already populated hard_reject. Hard-rejected
        ! rows remain rejected with probability 0; trainable rows are accepted
        ! exactly when P(accept) >= prob_threshold.
        do icls = 1, ncls
            if( quality%hard_reject(icls) ) cycle
            prob = pairwise_logistic_probability(model, quality%features(icls,:))
            quality%scores(icls) = prob
            if( prob >= model%prob_threshold )then
                quality%states(icls) = 1
                quality%labels(icls) = 1
            else
                quality%labels(icls) = 2
            endif
        end do
        quality%threshold        = model%prob_threshold
        quality%raw_threshold    = model%prob_threshold
        quality%threshold_offset = 0.0
        quality%separation       = 0.0
        quality%nclust           = 0
        quality%good_label       = 1
        quality%used_threshold   = .true.
        quality%model_name       = model%name
        quality%soft_decision    = 'probability_threshold'
        quality%soft_reason      = 'pairwise_logistic'
    end subroutine apply_pairwise_logistic

    real function pairwise_logistic_probability( model, feat )
        class(cavg_quality_model), intent(in) :: model
        real,                      intent(in) :: feat(:)
        integer :: iterm, ifeat, jfeat
        real    :: eta
        if( size(feat) /= CAVG_QUALITY_NFEATS ) THROW_HARD('pairwise_logistic_probability: invalid feature count')
        if( model%n_interactions < 0 .or. model%n_interactions > CAVG_QUALITY_MAX_INTERACTIONS ) &
            THROW_HARD('pairwise_logistic_probability: invalid n_interactions')
        eta = model%intercept + dot_product(model%linear_coefficients, feat)
        do iterm = 1, model%n_interactions
            ifeat = model%interaction_terms(iterm,1)
            jfeat = model%interaction_terms(iterm,2)
            if( ifeat < 1 .or. ifeat > CAVG_QUALITY_NFEATS .or. &
                jfeat < 1 .or. jfeat > CAVG_QUALITY_NFEATS ) &
                THROW_HARD('pairwise_logistic_probability: invalid interaction feature index')
            eta = eta + model%interaction_coefficients(iterm) * feat(ifeat) * feat(jfeat)
        end do
        eta = eta / max(EPS, model%calibration_temperature)
        eta = max(-80.0, min(80.0, eta))
        pairwise_logistic_probability = 1.0 / (1.0 + exp(-eta))
    end function pairwise_logistic_probability

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
