!@descr: instantiable class-average quality decision model
module simple_cavg_quality_model
use simple_defs,               only: LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string_utils,       only: lowercase, uppercase, fortran_symbol_from_string, fortran_quote
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, CAVG_QUALITY_MAX_INTERACTIONS, EPS, &
    CAVG_QUALITY_CONTEXT_CHUNK, CAVG_QUALITY_CONTEXT_POOL, CAVG_QUALITY_CONTEXT_SIEVE, CAVG_RELATIONAL_SCHEMA_NONE, &
    CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1, &
    CAVG_RELATIONAL_DEFAULT_KNN, CAVG_RELATIONAL_DEFAULT_CORR_HP, CAVG_RELATIONAL_DEFAULT_CORR_LP, &
    CAVG_RELATIONAL_DEFAULT_CORR_TRS, &
    cavg_quality_model_spec, cavg_quality_result, CAVG_REJECT_REASON_MODEL
use simple_cavg_quality_feats, only: cavg_quality_feature_name
implicit none
private
#include "simple_local_flags.inc"

public :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT
public :: CAVG_QUALITY_MODEL_SIEVE_DEFAULT
public :: CAVG_QUALITY_MODEL_POOL_DEFAULT
public :: CAVG_QUALITY_BUILTIN_MODELS
public :: cavg_quality_model
public :: cavg_quality_model_spec
public :: write_cavg_quality_model_builtin_code


! Built-in presets are complete model specifications. To promote a learned
! model into the code, add a named preset and include it in both public and
! printable built-in model inventories below.
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT = 'chunk100mics'
character(len=*), parameter :: CAVG_QUALITY_MODEL_SIEVE_DEFAULT = 'sieve'
character(len=*), parameter :: CAVG_QUALITY_MODEL_POOL_DEFAULT  = 'pool'
character(len=64), parameter :: CAVG_QUALITY_BUILTIN_MODELS(3) = [character(len=64) :: &
    CAVG_QUALITY_MODEL_CHUNK_DEFAULT, CAVG_QUALITY_MODEL_SIEVE_DEFAULT, CAVG_QUALITY_MODEL_POOL_DEFAULT]
character(len=*), parameter :: BUILTIN_MODEL_NAMES = CAVG_QUALITY_MODEL_CHUNK_DEFAULT // '|' // &
    CAVG_QUALITY_MODEL_SIEVE_DEFAULT // '|' // CAVG_QUALITY_MODEL_POOL_DEFAULT


! Default chunk class-average quality model, promoted from the pairwise
! logistic artifact learned from
! an out-of-tree chunk training set.
character(len=*), parameter :: CHUNK100MICS_FEATURE_POLICY = 'microchunk_plus_score_signal'
real, parameter :: CAVG_QUALITY_LOGISTIC_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02, 7.142857E-02, 7.142857E-02, &
    7.142857E-02, 7.142857E-02 ]

type :: cavg_quality_model
    character(len=64) :: name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
    character(len=32) :: context                 = 'chunk'
    character(len=64) :: feature_policy          = CHUNK100MICS_FEATURE_POLICY
    real              :: weights(CAVG_QUALITY_NFEATS) = CAVG_QUALITY_LOGISTIC_WEIGHTS
    real              :: intercept               = 0.0
    real              :: linear_coefficients(CAVG_QUALITY_NFEATS) = 0.0
    integer           :: n_interactions          = 0
    integer           :: interaction_terms(CAVG_QUALITY_MAX_INTERACTIONS,2) = 0
    real              :: interaction_coefficients(CAVG_QUALITY_MAX_INTERACTIONS) = 0.0
    real              :: prob_threshold          = 0.5
    real              :: regularization_lambda   = 0.0
    real              :: calibration_temperature = 1.0
    character(len=64) :: relational_feature_schema = CAVG_RELATIONAL_SCHEMA_NONE
    integer           :: relational_knn            = 0
    real              :: relational_corr_hp        = 0.0
    real              :: relational_corr_lp        = 0.0
    real              :: relational_corr_trs       = 0.0
    real              :: relational_coefficient    = 0.0
contains
    procedure :: init_preset
    procedure :: init_spec
    procedure :: get_spec
    procedure :: normalize
    procedure :: classify
    procedure :: supports_relational
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
            case(CAVG_QUALITY_MODEL_SIEVE_DEFAULT)
                spec = sieve_model_spec()
            case(CAVG_QUALITY_MODEL_POOL_DEFAULT)
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
        self%weights                 = spec%weights
        self%intercept               = spec%intercept
        self%linear_coefficients     = spec%linear_coefficients
        self%n_interactions          = spec%n_interactions
        self%interaction_terms       = spec%interaction_terms
        self%interaction_coefficients = spec%interaction_coefficients
        self%prob_threshold          = spec%prob_threshold
        self%regularization_lambda   = spec%regularization_lambda
        self%calibration_temperature = spec%calibration_temperature
        self%relational_feature_schema = trim(spec%relational_feature_schema)
        self%relational_knn          = spec%relational_knn
        self%relational_corr_hp      = spec%relational_corr_hp
        self%relational_corr_lp      = spec%relational_corr_lp
        self%relational_corr_trs     = spec%relational_corr_trs
        self%relational_coefficient  = spec%relational_coefficient
        call self%normalize()
    end subroutine init_spec

    function get_spec( self ) result( spec )
        class(cavg_quality_model), intent(in) :: self
        type(cavg_quality_model_spec) :: spec
        spec%name                    = self%name
        spec%context                 = self%context
        spec%feature_policy          = self%feature_policy
        spec%weights                 = self%weights
        spec%intercept               = self%intercept
        spec%linear_coefficients     = self%linear_coefficients
        spec%n_interactions          = self%n_interactions
        spec%interaction_terms       = self%interaction_terms
        spec%interaction_coefficients = self%interaction_coefficients
        spec%prob_threshold          = self%prob_threshold
        spec%regularization_lambda   = self%regularization_lambda
        spec%calibration_temperature = self%calibration_temperature
        spec%relational_feature_schema = self%relational_feature_schema
        spec%relational_knn          = self%relational_knn
        spec%relational_corr_hp      = self%relational_corr_hp
        spec%relational_corr_lp      = self%relational_corr_lp
        spec%relational_corr_trs     = self%relational_corr_trs
        spec%relational_coefficient  = self%relational_coefficient
    end function get_spec

    function chunk100mics_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
        spec%context                 = 'chunk'
        spec%feature_policy          = CHUNK100MICS_FEATURE_POLICY
        spec%weights                 = CAVG_QUALITY_LOGISTIC_WEIGHTS
        spec%intercept               = 2.802381E+00
        spec%linear_coefficients     = [ &
           -5.841123E-01,  6.280109E-01, -2.764060E-01,  1.036179E+00, &
           -7.276409E-01,  5.576301E-01, -1.771272E-01, -4.475693E-02, &
            2.538559E-01, -4.135939E-01,  7.407248E-01, -3.170226E-01, &
            2.768010E+00,  5.132661E-01 ]
        call set_pairwise_interactions_for_feature_count(spec, 14, [ &
           -5.016133E-01,  2.286104E-03,  6.567349E-01, -9.844457E-01, &
            4.001724E-01,  8.693924E-01,  7.284134E-01, -2.955010E-01, &
           -5.959087E-01,  8.986003E-01, -4.522021E-01,  5.249431E-01, &
            4.998736E-01,  7.446390E-03, -1.947948E-01,  7.215320E-02, &
           -1.647213E-01, -1.089008E-01, -5.241460E-01,  1.061381E-01, &
           -1.016107E-01,  2.078417E-01, -1.008714E-01,  1.364482E+00, &
           -1.701328E-01, -5.651891E-02,  1.481122E-01, -5.923793E-01, &
            2.441176E-01,  3.966064E-01, -1.921423E-01, -3.747033E-01, &
           -3.324885E-01,  3.396430E-01, -2.124085E-01, -2.583523E-01, &
           -2.673632E-01, -3.769121E-01,  9.726956E-01, -8.286975E-02, &
           -2.729041E-01, -4.568700E-01,  2.382788E-01, -6.906177E-01, &
           -1.946201E-01, -4.144796E-01,  2.073843E-01,  4.982973E-02, &
            6.386595E-01, -3.007336E-01, -1.482573E-02,  4.862371E-02, &
            2.569429E-01, -1.767143E-01,  2.898083E-01, -9.003410E-01, &
           -6.977128E-01,  6.913373E-02, -8.214340E-02,  5.055159E-02, &
            1.103404E-01, -4.825652E-01,  3.560530E-01,  2.035851E-01, &
           -5.634665E-01, -5.736839E-01, -1.286691E-01, -8.029981E-01, &
           -7.644567E-01, -4.376983E-02, -1.081729E-01, -1.196917E-02, &
           -4.971317E-01, -3.948403E-01, -1.028552E+00,  1.349212E-01, &
            2.116886E-02,  4.029729E-01,  7.093447E-01,  6.204602E-01, &
           -3.332528E-01,  1.358440E-02,  1.395677E-01, -5.286202E-01, &
            1.052782E-01, -6.519729E-01,  1.541351E-01, -2.955521E-01, &
           -1.872015E-01,  2.358026E-01,  2.090225E-02 ])
        spec%prob_threshold          = 3.500000E-01
        spec%regularization_lambda   = 1.000000E-03
        spec%calibration_temperature = 1.000000E+00
        spec%relational_feature_schema = CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1
        spec%relational_knn          = CAVG_RELATIONAL_DEFAULT_KNN
        spec%relational_corr_hp      = CAVG_RELATIONAL_DEFAULT_CORR_HP
        spec%relational_corr_lp      = CAVG_RELATIONAL_DEFAULT_CORR_LP
        spec%relational_corr_trs     = CAVG_RELATIONAL_DEFAULT_CORR_TRS
        spec%relational_coefficient  = -3.941003E-01
    end function chunk100mics_model_spec

    function pool_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_POOL_DEFAULT
        spec%context                 = 'pool'
        spec%feature_policy          = 'microchunk_plus_signal'
        spec%weights                 = [ &
              7.692308E-02,   7.692308E-02,   7.692308E-02,   7.692308E-02, &
              7.692308E-02,   0.000000E+00,   7.692308E-02,   7.692308E-02, &
              7.692308E-02,   7.692308E-02,   7.692308E-02,   7.692308E-02, &
              7.692308E-02,   7.692308E-02 ]
        spec%intercept               =   4.523899E+00
        spec%linear_coefficients = [ &
              9.266560E-01,   3.771429E+00,   3.566606E-01,   9.810680E-01, &
             -9.356901E-01,   0.000000E+00,   9.143119E-01,   1.431512E+00, &
             -8.987821E-01,  -1.057487E+00,   1.118601E+00,  -7.717300E-01, &
              2.762437E-01,   5.524297E-01 ]
        spec%n_interactions          = 78
        spec%interaction_terms        = 0
        spec%interaction_coefficients = 0.0
        spec%interaction_terms(1,:) = [1, 2 ]
        spec%interaction_coefficients(1) =   9.478272E-01
        spec%interaction_terms(2,:) = [1, 3 ]
        spec%interaction_coefficients(2) =   1.424371E-01
        spec%interaction_terms(3,:) = [1, 4 ]
        spec%interaction_coefficients(3) =   1.333783E+00
        spec%interaction_terms(4,:) = [1, 5 ]
        spec%interaction_coefficients(4) =  -1.301163E+00
        spec%interaction_terms(5,:) = [1, 7 ]
        spec%interaction_coefficients(5) =  -7.984090E-01
        spec%interaction_terms(6,:) = [1, 8 ]
        spec%interaction_coefficients(6) =  -4.896447E-01
        spec%interaction_terms(7,:) = [1, 9 ]
        spec%interaction_coefficients(7) =   1.605485E+00
        spec%interaction_terms(8,:) = [1, 10 ]
        spec%interaction_coefficients(8) =  -1.237449E+00
        spec%interaction_terms(9,:) = [1, 11 ]
        spec%interaction_coefficients(9) =   1.313880E+00
        spec%interaction_terms(10,:) = [1, 12 ]
        spec%interaction_coefficients(10) =  -7.300295E-01
        spec%interaction_terms(11,:) = [1, 13 ]
        spec%interaction_coefficients(11) =  -3.169714E-01
        spec%interaction_terms(12,:) = [1, 14 ]
        spec%interaction_coefficients(12) =  -1.348646E+00
        spec%interaction_terms(13,:) = [2, 3 ]
        spec%interaction_coefficients(13) =  -2.033189E-02
        spec%interaction_terms(14,:) = [2, 4 ]
        spec%interaction_coefficients(14) =   3.749562E-01
        spec%interaction_terms(15,:) = [2, 5 ]
        spec%interaction_coefficients(15) =  -4.521743E-03
        spec%interaction_terms(16,:) = [2, 7 ]
        spec%interaction_coefficients(16) =  -2.518075E-01
        spec%interaction_terms(17,:) = [2, 8 ]
        spec%interaction_coefficients(17) =   7.352542E-01
        spec%interaction_terms(18,:) = [2, 9 ]
        spec%interaction_coefficients(18) =   9.499735E-01
        spec%interaction_terms(19,:) = [2, 10 ]
        spec%interaction_coefficients(19) =  -2.957729E-01
        spec%interaction_terms(20,:) = [2, 11 ]
        spec%interaction_coefficients(20) =   1.594937E-01
        spec%interaction_terms(21,:) = [2, 12 ]
        spec%interaction_coefficients(21) =  -3.624845E-01
        spec%interaction_terms(22,:) = [2, 13 ]
        spec%interaction_coefficients(22) =   1.325966E-01
        spec%interaction_terms(23,:) = [2, 14 ]
        spec%interaction_coefficients(23) =  -8.330233E-01
        spec%interaction_terms(24,:) = [3, 4 ]
        spec%interaction_coefficients(24) =  -5.046661E-01
        spec%interaction_terms(25,:) = [3, 5 ]
        spec%interaction_coefficients(25) =   3.593823E-01
        spec%interaction_terms(26,:) = [3, 7 ]
        spec%interaction_coefficients(26) =   1.356772E-01
        spec%interaction_terms(27,:) = [3, 8 ]
        spec%interaction_coefficients(27) =   1.327510E-01
        spec%interaction_terms(28,:) = [3, 9 ]
        spec%interaction_coefficients(28) =   5.638146E-01
        spec%interaction_terms(29,:) = [3, 10 ]
        spec%interaction_coefficients(29) =   4.765443E-01
        spec%interaction_terms(30,:) = [3, 11 ]
        spec%interaction_coefficients(30) =  -3.199774E-01
        spec%interaction_terms(31,:) = [3, 12 ]
        spec%interaction_coefficients(31) =   6.913511E-01
        spec%interaction_terms(32,:) = [3, 13 ]
        spec%interaction_coefficients(32) =   1.444847E-01
        spec%interaction_terms(33,:) = [3, 14 ]
        spec%interaction_coefficients(33) =  -3.585712E-01
        spec%interaction_terms(34,:) = [4, 5 ]
        spec%interaction_coefficients(34) =  -1.855719E+00
        spec%interaction_terms(35,:) = [4, 7 ]
        spec%interaction_coefficients(35) =   4.347877E-02
        spec%interaction_terms(36,:) = [4, 8 ]
        spec%interaction_coefficients(36) =   1.792246E-01
        spec%interaction_terms(37,:) = [4, 9 ]
        spec%interaction_coefficients(37) =   5.590391E-01
        spec%interaction_terms(38,:) = [4, 10 ]
        spec%interaction_coefficients(38) =  -2.152137E+00
        spec%interaction_terms(39,:) = [4, 11 ]
        spec%interaction_coefficients(39) =   1.479134E+00
        spec%interaction_terms(40,:) = [4, 12 ]
        spec%interaction_coefficients(40) =  -8.113728E-01
        spec%interaction_terms(41,:) = [4, 13 ]
        spec%interaction_coefficients(41) =   5.293784E-01
        spec%interaction_terms(42,:) = [4, 14 ]
        spec%interaction_coefficients(42) =  -4.795363E-01
        spec%interaction_terms(43,:) = [5, 7 ]
        spec%interaction_coefficients(43) =  -1.689706E-01
        spec%interaction_terms(44,:) = [5, 8 ]
        spec%interaction_coefficients(44) =  -2.129593E-01
        spec%interaction_terms(45,:) = [5, 9 ]
        spec%interaction_coefficients(45) =   1.167485E+00
        spec%interaction_terms(46,:) = [5, 10 ]
        spec%interaction_coefficients(46) =   1.933542E+00
        spec%interaction_terms(47,:) = [5, 11 ]
        spec%interaction_coefficients(47) =  -3.856001E+00
        spec%interaction_terms(48,:) = [5, 12 ]
        spec%interaction_coefficients(48) =   1.291968E+00
        spec%interaction_terms(49,:) = [5, 13 ]
        spec%interaction_coefficients(49) =  -1.965936E-01
        spec%interaction_terms(50,:) = [5, 14 ]
        spec%interaction_coefficients(50) =  -9.229491E-01
        spec%interaction_terms(51,:) = [7, 8 ]
        spec%interaction_coefficients(51) =  -1.591857E-01
        spec%interaction_terms(52,:) = [7, 9 ]
        spec%interaction_coefficients(52) =   2.809146E-01
        spec%interaction_terms(53,:) = [7, 10 ]
        spec%interaction_coefficients(53) =  -3.236808E-02
        spec%interaction_terms(54,:) = [7, 11 ]
        spec%interaction_coefficients(54) =  -9.184004E-02
        spec%interaction_terms(55,:) = [7, 12 ]
        spec%interaction_coefficients(55) =  -3.821822E-01
        spec%interaction_terms(56,:) = [7, 13 ]
        spec%interaction_coefficients(56) =   3.135403E-01
        spec%interaction_terms(57,:) = [7, 14 ]
        spec%interaction_coefficients(57) =   4.656276E-01
        spec%interaction_terms(58,:) = [8, 9 ]
        spec%interaction_coefficients(58) =   2.569876E-01
        spec%interaction_terms(59,:) = [8, 10 ]
        spec%interaction_coefficients(59) =  -1.775942E-01
        spec%interaction_terms(60,:) = [8, 11 ]
        spec%interaction_coefficients(60) =   4.600332E-01
        spec%interaction_terms(61,:) = [8, 12 ]
        spec%interaction_coefficients(61) =   3.049327E-01
        spec%interaction_terms(62,:) = [8, 13 ]
        spec%interaction_coefficients(62) =  -1.334592E-01
        spec%interaction_terms(63,:) = [8, 14 ]
        spec%interaction_coefficients(63) =   2.507228E-01
        spec%interaction_terms(64,:) = [9, 10 ]
        spec%interaction_coefficients(64) =  -5.893188E-01
        spec%interaction_terms(65,:) = [9, 11 ]
        spec%interaction_coefficients(65) =  -7.615877E-01
        spec%interaction_terms(66,:) = [9, 12 ]
        spec%interaction_coefficients(66) =  -1.150061E+00
        spec%interaction_terms(67,:) = [9, 13 ]
        spec%interaction_coefficients(67) =   1.613582E+00
        spec%interaction_terms(68,:) = [9, 14 ]
        spec%interaction_coefficients(68) =   8.980396E-02
        spec%interaction_terms(69,:) = [10, 11 ]
        spec%interaction_coefficients(69) =  -1.556036E+00
        spec%interaction_terms(70,:) = [10, 12 ]
        spec%interaction_coefficients(70) =   7.161171E-01
        spec%interaction_terms(71,:) = [10, 13 ]
        spec%interaction_coefficients(71) =  -4.010158E-01
        spec%interaction_terms(72,:) = [10, 14 ]
        spec%interaction_coefficients(72) =   4.667039E-01
        spec%interaction_terms(73,:) = [11, 12 ]
        spec%interaction_coefficients(73) =  -1.506662E+00
        spec%interaction_terms(74,:) = [11, 13 ]
        spec%interaction_coefficients(74) =  -2.137732E-02
        spec%interaction_terms(75,:) = [11, 14 ]
        spec%interaction_coefficients(75) =   1.287851E+00
        spec%interaction_terms(76,:) = [12, 13 ]
        spec%interaction_coefficients(76) =  -1.019148E+00
        spec%interaction_terms(77,:) = [12, 14 ]
        spec%interaction_coefficients(77) =   7.272188E-01
        spec%interaction_terms(78,:) = [13, 14 ]
        spec%interaction_coefficients(78) =  -1.638013E+00
        spec%prob_threshold          =   3.000000E-01
        spec%regularization_lambda   =   1.000000E-04
        spec%calibration_temperature =   1.000000E+00
        spec%relational_feature_schema = 'corr_knn_signal_v1'
        spec%relational_knn          = 5
        spec%relational_corr_hp      =   1.000000E+02
        spec%relational_corr_lp      =   1.500000E+01
        spec%relational_corr_trs     =   1.000000E+01
        spec%relational_coefficient  =   2.449666E-02
    end function pool_model_spec

    function sieve_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_SIEVE_DEFAULT
        spec%context                 = 'sieve'
        spec%feature_policy          = 'microchunk_plus_score_signal'
        spec%weights                 = [ &
              7.142857E-02,   7.142857E-02,   7.142857E-02,   7.142857E-02, &
              7.142857E-02,   7.142857E-02,   7.142857E-02,   7.142857E-02, &
              7.142857E-02,   7.142857E-02,   7.142857E-02,   7.142857E-02, &
              7.142857E-02,   7.142857E-02 ]
        spec%intercept               =   9.832671E-01
        spec%linear_coefficients     = [ &
             -6.918444E-02,   1.174217E+00,   2.082862E-01,   9.721558E-01, &
              2.745598E-01,  -6.372386E-01,  -2.369738E-01,   3.745508E-01, &
             -8.853315E-01,  -8.241642E-01,  -1.025514E-01,  -2.550836E-02, &
              3.438007E-01,   1.061425E+00 ]
        spec%n_interactions           = 91
        spec%interaction_terms        = 0
        spec%interaction_coefficients = 0.0
        spec%interaction_terms(1,:) = [1, 2 ]
        spec%interaction_coefficients(1) =   6.158664E-01
        spec%interaction_terms(2,:) = [1, 3 ]
        spec%interaction_coefficients(2) =  -7.278305E-02
        spec%interaction_terms(3,:) = [1, 4 ]
        spec%interaction_coefficients(3) =  -5.793010E-02
        spec%interaction_terms(4,:) = [1, 5 ]
        spec%interaction_coefficients(4) =  -1.581391E-01
        spec%interaction_terms(5,:) = [1, 6 ]
        spec%interaction_coefficients(5) =   2.160046E-01
        spec%interaction_terms(6,:) = [1, 7 ]
        spec%interaction_coefficients(6) =  -1.268776E-01
        spec%interaction_terms(7,:) = [1, 8 ]
        spec%interaction_coefficients(7) =   3.217513E-01
        spec%interaction_terms(8,:) = [1, 9 ]
        spec%interaction_coefficients(8) =  -1.239059E-01
        spec%interaction_terms(9,:) = [1, 10 ]
        spec%interaction_coefficients(9) =  -2.689630E-01
        spec%interaction_terms(10,:) = [1, 11 ]
        spec%interaction_coefficients(10) =   4.064969E-01
        spec%interaction_terms(11,:) = [1, 12 ]
        spec%interaction_coefficients(11) =   2.079959E-01
        spec%interaction_terms(12,:) = [1, 13 ]
        spec%interaction_coefficients(12) =   8.464044E-02
        spec%interaction_terms(13,:) = [1, 14 ]
        spec%interaction_coefficients(13) =   2.688188E-02
        spec%interaction_terms(14,:) = [2, 3 ]
        spec%interaction_coefficients(14) =   2.528967E-01
        spec%interaction_terms(15,:) = [2, 4 ]
        spec%interaction_coefficients(15) =  -4.502688E-02
        spec%interaction_terms(16,:) = [2, 5 ]
        spec%interaction_coefficients(16) =   3.397096E-01
        spec%interaction_terms(17,:) = [2, 6 ]
        spec%interaction_coefficients(17) =   1.369977E-01
        spec%interaction_terms(18,:) = [2, 7 ]
        spec%interaction_coefficients(18) =  -3.075861E-01
        spec%interaction_terms(19,:) = [2, 8 ]
        spec%interaction_coefficients(19) =  -4.387426E-01
        spec%interaction_terms(20,:) = [2, 9 ]
        spec%interaction_coefficients(20) =   2.256129E-01
        spec%interaction_terms(21,:) = [2, 10 ]
        spec%interaction_coefficients(21) =   8.311773E-02
        spec%interaction_terms(22,:) = [2, 11 ]
        spec%interaction_coefficients(22) =  -1.010890E-01
        spec%interaction_terms(23,:) = [2, 12 ]
        spec%interaction_coefficients(23) =  -6.351334E-01
        spec%interaction_terms(24,:) = [2, 13 ]
        spec%interaction_coefficients(24) =   6.612171E-01
        spec%interaction_terms(25,:) = [2, 14 ]
        spec%interaction_coefficients(25) =  -1.799542E-01
        spec%interaction_terms(26,:) = [3, 4 ]
        spec%interaction_coefficients(26) =  -6.657630E-01
        spec%interaction_terms(27,:) = [3, 5 ]
        spec%interaction_coefficients(27) =   5.755506E-01
        spec%interaction_terms(28,:) = [3, 6 ]
        spec%interaction_coefficients(28) =  -7.396340E-02
        spec%interaction_terms(29,:) = [3, 7 ]
        spec%interaction_coefficients(29) =   3.756891E-02
        spec%interaction_terms(30,:) = [3, 8 ]
        spec%interaction_coefficients(30) =  -7.736346E-02
        spec%interaction_terms(31,:) = [3, 9 ]
        spec%interaction_coefficients(31) =  -1.425177E-02
        spec%interaction_terms(32,:) = [3, 10 ]
        spec%interaction_coefficients(32) =   1.281168E-01
        spec%interaction_terms(33,:) = [3, 11 ]
        spec%interaction_coefficients(33) =   3.084263E-02
        spec%interaction_terms(34,:) = [3, 12 ]
        spec%interaction_coefficients(34) =   1.866026E-01
        spec%interaction_terms(35,:) = [3, 13 ]
        spec%interaction_coefficients(35) =  -1.912241E-01
        spec%interaction_terms(36,:) = [3, 14 ]
        spec%interaction_coefficients(36) =   1.800826E-01
        spec%interaction_terms(37,:) = [4, 5 ]
        spec%interaction_coefficients(37) =  -2.024839E-01
        spec%interaction_terms(38,:) = [4, 6 ]
        spec%interaction_coefficients(38) =  -2.072305E-01
        spec%interaction_terms(39,:) = [4, 7 ]
        spec%interaction_coefficients(39) =   2.061003E-01
        spec%interaction_terms(40,:) = [4, 8 ]
        spec%interaction_coefficients(40) =  -4.541796E-01
        spec%interaction_terms(41,:) = [4, 9 ]
        spec%interaction_coefficients(41) =   5.552803E-01
        spec%interaction_terms(42,:) = [4, 10 ]
        spec%interaction_coefficients(42) =  -1.353896E+00
        spec%interaction_terms(43,:) = [4, 11 ]
        spec%interaction_coefficients(43) =   6.565881E-01
        spec%interaction_terms(44,:) = [4, 12 ]
        spec%interaction_coefficients(44) =  -5.576389E-01
        spec%interaction_terms(45,:) = [4, 13 ]
        spec%interaction_coefficients(45) =  -5.728185E-02
        spec%interaction_terms(46,:) = [4, 14 ]
        spec%interaction_coefficients(46) =  -3.663005E-01
        spec%interaction_terms(47,:) = [5, 6 ]
        spec%interaction_coefficients(47) =   9.450029E-01
        spec%interaction_terms(48,:) = [5, 7 ]
        spec%interaction_coefficients(48) =   3.618473E-01
        spec%interaction_terms(49,:) = [5, 8 ]
        spec%interaction_coefficients(49) =  -5.069597E-02
        spec%interaction_terms(50,:) = [5, 9 ]
        spec%interaction_coefficients(50) =   8.555743E-01
        spec%interaction_terms(51,:) = [5, 10 ]
        spec%interaction_coefficients(51) =   1.143771E-01
        spec%interaction_terms(52,:) = [5, 11 ]
        spec%interaction_coefficients(52) =   4.424291E-01
        spec%interaction_terms(53,:) = [5, 12 ]
        spec%interaction_coefficients(53) =   4.453060E-01
        spec%interaction_terms(54,:) = [5, 13 ]
        spec%interaction_coefficients(54) =   6.675552E-01
        spec%interaction_terms(55,:) = [5, 14 ]
        spec%interaction_coefficients(55) =   1.489348E-01
        spec%interaction_terms(56,:) = [6, 7 ]
        spec%interaction_coefficients(56) =  -2.402572E-01
        spec%interaction_terms(57,:) = [6, 8 ]
        spec%interaction_coefficients(57) =   3.175540E-01
        spec%interaction_terms(58,:) = [6, 9 ]
        spec%interaction_coefficients(58) =   2.686897E-02
        spec%interaction_terms(59,:) = [6, 10 ]
        spec%interaction_coefficients(59) =  -1.394199E-01
        spec%interaction_terms(60,:) = [6, 11 ]
        spec%interaction_coefficients(60) =   2.202978E-01
        spec%interaction_terms(61,:) = [6, 12 ]
        spec%interaction_coefficients(61) =   3.135217E-02
        spec%interaction_terms(62,:) = [6, 13 ]
        spec%interaction_coefficients(62) =  -1.774485E-01
        spec%interaction_terms(63,:) = [6, 14 ]
        spec%interaction_coefficients(63) =  -6.164746E-01
        spec%interaction_terms(64,:) = [7, 8 ]
        spec%interaction_coefficients(64) =  -3.134043E-01
        spec%interaction_terms(65,:) = [7, 9 ]
        spec%interaction_coefficients(65) =   2.526013E-01
        spec%interaction_terms(66,:) = [7, 10 ]
        spec%interaction_coefficients(66) =  -6.858773E-01
        spec%interaction_terms(67,:) = [7, 11 ]
        spec%interaction_coefficients(67) =   4.007939E-01
        spec%interaction_terms(68,:) = [7, 12 ]
        spec%interaction_coefficients(68) =  -7.243691E-01
        spec%interaction_terms(69,:) = [7, 13 ]
        spec%interaction_coefficients(69) =   1.833894E-01
        spec%interaction_terms(70,:) = [7, 14 ]
        spec%interaction_coefficients(70) =  -1.969232E-01
        spec%interaction_terms(71,:) = [8, 9 ]
        spec%interaction_coefficients(71) =  -3.968418E-01
        spec%interaction_terms(72,:) = [8, 10 ]
        spec%interaction_coefficients(72) =   3.011461E-01
        spec%interaction_terms(73,:) = [8, 11 ]
        spec%interaction_coefficients(73) =   1.477741E-01
        spec%interaction_terms(74,:) = [8, 12 ]
        spec%interaction_coefficients(74) =  -1.669138E-01
        spec%interaction_terms(75,:) = [8, 13 ]
        spec%interaction_coefficients(75) =   4.175667E-01
        spec%interaction_terms(76,:) = [8, 14 ]
        spec%interaction_coefficients(76) =   7.915357E-01
        spec%interaction_terms(77,:) = [9, 10 ]
        spec%interaction_coefficients(77) =  -6.477614E-01
        spec%interaction_terms(78,:) = [9, 11 ]
        spec%interaction_coefficients(78) =   1.153760E-01
        spec%interaction_terms(79,:) = [9, 12 ]
        spec%interaction_coefficients(79) =   9.445300E-01
        spec%interaction_terms(80,:) = [9, 13 ]
        spec%interaction_coefficients(80) =  -4.692881E-01
        spec%interaction_terms(81,:) = [9, 14 ]
        spec%interaction_coefficients(81) =  -1.464467E-01
        spec%interaction_terms(82,:) = [10, 11 ]
        spec%interaction_coefficients(82) =  -5.707061E-01
        spec%interaction_terms(83,:) = [10, 12 ]
        spec%interaction_coefficients(83) =   3.766019E-01
        spec%interaction_terms(84,:) = [10, 13 ]
        spec%interaction_coefficients(84) =  -1.581933E-01
        spec%interaction_terms(85,:) = [10, 14 ]
        spec%interaction_coefficients(85) =   1.928267E-01
        spec%interaction_terms(86,:) = [11, 12 ]
        spec%interaction_coefficients(86) =  -4.707946E-01
        spec%interaction_terms(87,:) = [11, 13 ]
        spec%interaction_coefficients(87) =  -8.185964E-01
        spec%interaction_terms(88,:) = [11, 14 ]
        spec%interaction_coefficients(88) =   6.579612E-01
        spec%interaction_terms(89,:) = [12, 13 ]
        spec%interaction_coefficients(89) =   4.306115E-01
        spec%interaction_terms(90,:) = [12, 14 ]
        spec%interaction_coefficients(90) =  -7.650308E-01
        spec%interaction_terms(91,:) = [13, 14 ]
        spec%interaction_coefficients(91) =  -1.115168E-01
        spec%prob_threshold          =   2.500000E-01
        spec%regularization_lambda   =   1.000000E-03
        spec%calibration_temperature =   1.000000E+00
        spec%relational_feature_schema = 'corr_knn_signal_v1'
        spec%relational_knn          = 5
        spec%relational_corr_hp      =   1.000000E+02
        spec%relational_corr_lp      =   1.500000E+01
        ! Promoted snippet (a93af895b) read 1.000000E+0, a truncated ES14.6 literal.
        ! The sieve training tables were produced at the shared default (10.0 px);
        ! the learner propagates that value, so the shift range is the shared default.
        spec%relational_corr_trs     = CAVG_RELATIONAL_DEFAULT_CORR_TRS
        spec%relational_coefficient  =  -3.960855E-01
    end function sieve_model_spec

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

    subroutine normalize( self )
        class(cavg_quality_model), intent(inout) :: self
        select case(trim(self%context))
        case(CAVG_QUALITY_CONTEXT_CHUNK, CAVG_QUALITY_CONTEXT_POOL, CAVG_QUALITY_CONTEXT_SIEVE)
            continue
        case default
            THROW_HARD('normalize: model context must be chunk, pool, or sieve')
        end select
        select case(trim(self%relational_feature_schema))
        case(CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1)
            if( self%relational_knn < 1 ) THROW_HARD('normalize: relational_knn must be positive')
            if( self%relational_corr_hp <= 0.0 .or. self%relational_corr_lp <= 0.0 .or. &
                self%relational_corr_hp < self%relational_corr_lp ) &
                THROW_HARD('normalize: invalid relational correlation limits')
            if( self%relational_corr_trs < 0.0 ) THROW_HARD('normalize: relational shift range must be nonnegative')
        case default
            THROW_HARD('normalize: relational_feature_schema=corr_knn_signal_v1 is required')
        end select
        self%prob_threshold          = min(1.0, max(0.0, self%prob_threshold))
        self%calibration_temperature = max(EPS, self%calibration_temperature)
        if( self%n_interactions < 0 .or. self%n_interactions > CAVG_QUALITY_MAX_INTERACTIONS ) &
            THROW_HARD('normalize: invalid pairwise interaction count')
    end subroutine normalize

    logical function supports_relational( self )
        class(cavg_quality_model), intent(in) :: self
        supports_relational = trim(self%relational_feature_schema) == CAVG_RELATIONAL_SCHEMA_CORR_KNN_SIGNAL_V1
    end function supports_relational

    subroutine classify( self, quality, relational_feature )
        class(cavg_quality_model), intent(in)    :: self
        type(cavg_quality_result), intent(inout) :: quality
        real, optional,            intent(in)    :: relational_feature(:)
        if( .not. allocated(quality%features)    ) THROW_HARD('classify: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('classify: missing hard-reject mask')
        quality%model_name     = self%name
        if( .not. present(relational_feature) ) THROW_HARD('classify: missing relational feature')
        call apply_pairwise_logistic(quality, self, relational_feature)
    end subroutine classify

    subroutine write_model( self, fname )
        class(cavg_quality_model), intent(in) :: self
        character(len=*),          intent(in) :: fname
        integer :: funit, i
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# model_cavgs_rejection model'
        write(funit,'(A)') 'model_version=11'
        write(funit,'(A,A)') 'name=', trim(self%name)
        write(funit,'(A,A)') 'context=', trim(self%context)
        write(funit,'(A,A)') 'feature_policy=', trim(self%feature_policy)
        write(funit,'(A)', advance='no') 'feature_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') self%weights(i)
        end do
        write(funit,*)
        write(funit,'(A,ES14.6)') 'intercept=', self%intercept
        call write_model_real_list(funit, 'linear_coefficients=', self%linear_coefficients)
        call write_interaction_terms(funit, self)
        call write_model_real_list(funit, 'interaction_coefficients=', self%interaction_coefficients, self%n_interactions)
        write(funit,'(A,ES14.6)') 'prob_threshold=', self%prob_threshold
        write(funit,'(A,ES14.6)') 'regularization_lambda=', self%regularization_lambda
        write(funit,'(A,ES14.6)') 'calibration_temperature=', self%calibration_temperature
        write(funit,'(A,A)') 'relational_feature_schema=', trim(self%relational_feature_schema)
        write(funit,'(A,I0)') 'relational_knn=', self%relational_knn
        write(funit,'(A,ES14.6)') 'relational_corr_hp=', self%relational_corr_hp
        write(funit,'(A,ES14.6)') 'relational_corr_lp=', self%relational_corr_lp
        write(funit,'(A,ES14.6)') 'relational_corr_trs=', self%relational_corr_trs
        write(funit,'(A,ES14.6)') 'relational_coefficient=', self%relational_coefficient
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
        write(funit,'(A)') '! 2. Append this name to CAVG_QUALITY_BUILTIN_MODELS and BUILTIN_MODEL_NAMES:'
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
        call write_logistic_spec_assignments(funit, model)
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
        write(funit,'(A,A)') '        spec%relational_feature_schema = ', &
            trim(fortran_quote(model%relational_feature_schema))
        write(funit,'(A,I0)') '        spec%relational_knn          = ', model%relational_knn
        write(funit,'(A,ES14.6)') '        spec%relational_corr_hp      = ', model%relational_corr_hp
        write(funit,'(A,ES14.6)') '        spec%relational_corr_lp      = ', model%relational_corr_lp
        write(funit,'(A,ES14.6)') '        spec%relational_corr_trs     = ', model%relational_corr_trs
        write(funit,'(A,ES14.6)') '        spec%relational_coefficient  = ', model%relational_coefficient
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
        logical :: have_relational_schema, have_preset, ok_line
        ! Model files are complete model definitions. Start from chunk defaults,
        ! apply any preset found in the file, then apply explicit key overrides.
        call self%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
        open(newunit=funit, file=trim(fname), status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('read_model: failed to open '//trim(fname))
        have_relational_schema = .false.
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
            case('relational_feature_schema')
                have_relational_schema = .true.
            case('preset')
                preset_name = trim(val)
                have_preset = .true.
            end select
        end do
        if( have_preset ) call self%init_preset(trim(preset_name))
        if( model_version /= 11 .or. .not. have_relational_schema ) &
            THROW_HARD('read_model: relational logistic model_version=11 is required')
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
                case('context')
                    select case(trim(val))
                        case(CAVG_QUALITY_CONTEXT_CHUNK, CAVG_QUALITY_CONTEXT_POOL, CAVG_QUALITY_CONTEXT_SIEVE)
                            self%context = trim(val)
                        case DEFAULT
                            THROW_HARD('read_model: context must be chunk, pool, or sieve')
                    end select
                case('feature_policy')
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
                case('relational_feature_schema')
                    self%relational_feature_schema = trim(val)
                case('relational_knn')
                    read(val,*,iostat=parse_ios) self%relational_knn
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse relational_knn')
                case('relational_corr_hp')
                    read(val,*,iostat=parse_ios) self%relational_corr_hp
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse relational_corr_hp')
                case('relational_corr_lp')
                    read(val,*,iostat=parse_ios) self%relational_corr_lp
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse relational_corr_lp')
                case('relational_corr_trs')
                    read(val,*,iostat=parse_ios) self%relational_corr_trs
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse relational_corr_trs')
                case('relational_coefficient')
                    read(val,*,iostat=parse_ios) self%relational_coefficient
                    if( parse_ios /= 0 ) THROW_HARD('read_model: failed to parse relational_coefficient')
                case default
                    THROW_HARD('read_model: unknown key in model file: '//trim(key))
            end select
        end do
        close(funit)
        if( n_interaction_coefficients /= self%n_interactions ) &
            THROW_HARD('read_model: interaction_coefficients count must match interaction_terms count')
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
        self%weights                 = 0.0
        self%intercept               = 0.0
        self%linear_coefficients     = 0.0
        self%n_interactions          = 0
        self%interaction_terms       = 0
        self%interaction_coefficients = 0.0
        self%prob_threshold          = 0.0
        self%regularization_lambda   = 0.0
        self%calibration_temperature = 1.0
        self%relational_feature_schema = CAVG_RELATIONAL_SCHEMA_NONE
        self%relational_knn          = 0
        self%relational_corr_hp      = 0.0
        self%relational_corr_lp      = 0.0
        self%relational_corr_trs     = 0.0
        self%relational_coefficient  = 0.0
    end subroutine kill_model

    subroutine apply_pairwise_logistic( quality, model, relational_feature )
        type(cavg_quality_result), intent(inout) :: quality
        class(cavg_quality_model), intent(in)    :: model
        real, optional,            intent(in)    :: relational_feature(:)
        integer :: ncls, icls
        real    :: prob
        if( .not. allocated(quality%features)    ) THROW_HARD('apply_pairwise_logistic: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('apply_pairwise_logistic: missing hard-reject mask')
        if( size(quality%features, dim=2) /= CAVG_QUALITY_NFEATS ) THROW_HARD('apply_pairwise_logistic: invalid feature count')
        ncls = size(quality%features, dim=1)
        if( size(quality%hard_reject) /= ncls ) THROW_HARD('apply_pairwise_logistic: invalid mask size')
        if( model%supports_relational() )then
            if( .not. present(relational_feature) ) &
                THROW_HARD('apply_pairwise_logistic: model requires relational feature')
            if( size(relational_feature) /= ncls ) &
                THROW_HARD('apply_pairwise_logistic: relational feature size mismatch')
        end if
        if( allocated(quality%states)  ) deallocate(quality%states)
        if( allocated(quality%labels)  ) deallocate(quality%labels)
        if( allocated(quality%medoids) ) deallocate(quality%medoids)
        if( allocated(quality%scores)  ) deallocate(quality%scores)
        allocate(quality%states(ncls), quality%labels(ncls), source=0)
        allocate(quality%scores(ncls), source=0.0)
        ! quality%reasons is populated by the upstream hard-gate pass (extract_cavg_quality_features)
        ! with the specific reason (e.g. population, bad pixels, mask geometry) for every hard-rejected
        ! row. Do not wipe it here: hard-rejected rows are skipped below (they never reach the
        ! model-reject branch), so resetting the array to CAVG_REJECT_REASON_NONE would silently
        ! discard their real rejection reason.
        if( .not. allocated(quality%reasons) ) then
            allocate(quality%reasons(ncls), source=0)
        else if( size(quality%reasons) /= ncls ) then
            THROW_HARD('apply_pairwise_logistic: invalid reasons size')
        end if

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
            if( model%supports_relational() )then
                prob = pairwise_logistic_probability(model, quality%features(icls,:), relational_feature(icls))
            else
                prob = pairwise_logistic_probability(model, quality%features(icls,:))
            end if
            quality%scores(icls) = prob
            if( prob >= model%prob_threshold )then
                quality%states(icls) = 1
                quality%labels(icls) = 1
            else
                quality%reasons(icls) = CAVG_REJECT_REASON_MODEL
                quality%labels(icls)  = 2
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
        if( model%supports_relational() ) quality%soft_reason = 'pairwise_logistic_relational'
    end subroutine apply_pairwise_logistic

    real function pairwise_logistic_probability( model, feat, relational_feature )
        class(cavg_quality_model), intent(in) :: model
        real,                      intent(in) :: feat(:)
        real, optional,            intent(in) :: relational_feature
        integer :: iterm, ifeat, jfeat
        real    :: eta
        if( size(feat) /= CAVG_QUALITY_NFEATS ) THROW_HARD('pairwise_logistic_probability: invalid feature count')
        if( model%n_interactions < 0 .or. model%n_interactions > CAVG_QUALITY_MAX_INTERACTIONS ) &
            THROW_HARD('pairwise_logistic_probability: invalid n_interactions')
        eta = model%intercept + dot_product(model%linear_coefficients, feat)
        if( model%supports_relational() )then
            if( .not. present(relational_feature) ) &
                THROW_HARD('pairwise_logistic_probability: missing relational feature')
            eta = eta + model%relational_coefficient * relational_feature
        end if
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

end module simple_cavg_quality_model
