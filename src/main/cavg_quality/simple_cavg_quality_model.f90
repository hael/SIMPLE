!@descr: instantiable class-average quality decision model
module simple_cavg_quality_model
use simple_defs,               only: LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string_utils,       only: str_is_true, csv_field, lowercase, uppercase
use simple_clustering_utils,   only: cluster_dmat
use simple_srch_sort_loc,      only: hpsort
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, cavg_quality_model_spec, cavg_quality_result
use simple_cavg_quality_stats, only: normalize_quality_dmat
implicit none
private
#include "simple_local_flags.inc"

public :: CAVG_QUALITY_DEFAULT_WEIGHTS
public :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT
public :: CAVG_QUALITY_MODEL_POOL_DEFAULT
public :: cavg_quality_model
public :: cavg_quality_model_spec
public :: write_cavg_quality_model_builtin_code

! Built-in presets are complete model specifications. To promote a learned
! model into the code, add a named preset and include it in builtin_names.
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT = 'chunk_default_v2'
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT_V1 = 'chunk_default_v1'
character(len=*), parameter :: CAVG_QUALITY_MODEL_POOL_DEFAULT  = 'pool_default_v1'
character(len=*), parameter :: BUILTIN_MODEL_NAMES = CAVG_QUALITY_MODEL_CHUNK_DEFAULT//'|'//&
    CAVG_QUALITY_MODEL_CHUNK_DEFAULT_V1//'|'//CAVG_QUALITY_MODEL_POOL_DEFAULT

real, parameter :: MIN_SCORE_SEPARATION    = 0.15
real, parameter :: BOUNDARY_MARGIN_DEFAULT = 0.05
real, parameter :: CHUNK_BOUNDARY_OFFSET   = 0.25
real, parameter :: CHUNK_OTSU_MIN_OFFSET   = 0.25
real, parameter :: CHUNK_OTSU_MAX_OFFSET   = 0.50
real, parameter :: POOL_MIN_ACCEPT_FRAC    = 0.65
real, parameter :: CLUSTER_RESCUE_MARGIN   = 0.20

! Reference v1/base linear score coefficients. These are retained for the
! legacy preset and as the learner fallback; feature definitions and extraction
! live in simple_cavg_quality_feats.
real, parameter :: CAVG_QUALITY_DEFAULT_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    3.953488E-01, 2.093023E-01, 0.000000E+00, 0.000000E+00, &
    1.860465E-01, 2.093023E-01, 0.000000E+00, 0.000000E+00, &
    0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, &
    0.000000E+00, 0.000000E+00, 0.000000E+00 ]

! Batch-trained chunk default promoted from the representative batch10chunk
! learning round after hard rejection was separated from model fitting and Otsu
! threshold policy was added to the internal learn grid. The selected policy
! keeps the cheap scalar diagnostics that generalized best on stream chunks,
! while zeroing mask_inside, single_component, and cc_area_frac because those
! soft geometry flags duplicate hard geometry rejection or were inconsistent
! across the representative chunk set. Pairwise distance matrices and unstable
! spectrum, ice, and foreground/background ratio diagnostics have been removed
! from the default scalar model space.
real, parameter :: CAVG_QUALITY_CHUNK_V2_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    8.404870E-02, 1.065039E-01, 0.000000E+00, 3.393937E-02, &
    8.439852E-02, 8.831557E-02, 0.000000E+00, 1.135785E-01, &
    6.035764E-02, 7.611332E-02, 0.000000E+00, 5.907804E-02, &
    9.650873E-02, 1.008953E-01, 9.626245E-02 ]
real, parameter :: CHUNK_V2_BOUNDARY_MARGIN      =  0.10
real, parameter :: CHUNK_V2_MIN_SCORE_SEPARATION =  0.15
real, parameter :: CHUNK_V2_OTSU_MIN_OFFSET      =  0.15
real, parameter :: CHUNK_V2_OTSU_MAX_OFFSET      =  0.50

type :: cavg_quality_model
    character(len=64) :: name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
    character(len=32) :: family                  = 'linear_boundary'
    character(len=32) :: context                 = 'chunk'
    character(len=32) :: feature_policy          = 'all15_no_geom_softs'
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
    procedure :: kill  => reset_model
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
            case(CAVG_QUALITY_MODEL_CHUNK_DEFAULT_V1)
                spec = chunk_default_v1_model_spec()
            case(CAVG_QUALITY_MODEL_POOL_DEFAULT)
                spec = pool_default_model_spec()
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
        if( trim(spec%family) /= 'linear_boundary' ) &
            THROW_HARD('unsupported class-average quality model family: '//trim(spec%family))
        self%name                    = trim(spec%name)
        self%family                  = trim(spec%family)
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
        spec%family                  = self%family
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
        spec%family                  = 'linear_boundary'
        spec%context                 = 'chunk'
        spec%feature_policy          = 'all15_no_geom_softs'
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

    function chunk_default_v1_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT_V1
        spec%family                  = 'linear_boundary'
        spec%context                 = 'chunk'
        spec%feature_policy          = 'legacy_core'
        spec%weights                 = CAVG_QUALITY_DEFAULT_WEIGHTS
        spec%boundary_margin         = -CHUNK_BOUNDARY_OFFSET
        spec%min_score_separation    = MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = CHUNK_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = 0.0
        spec%use_lowsep_otsu         = .true.
        spec%use_otsu_window         = .true.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function chunk_default_v1_model_spec

    function pool_default_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_POOL_DEFAULT
        spec%family                  = 'linear_boundary'
        spec%context                 = 'pool'
        spec%feature_policy          = 'legacy_core'
        spec%weights                 = CAVG_QUALITY_DEFAULT_WEIGHTS
        spec%boundary_margin         = BOUNDARY_MARGIN_DEFAULT
        spec%min_score_separation    = MIN_SCORE_SEPARATION
        spec%otsu_min_offset         = CHUNK_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = POOL_MIN_ACCEPT_FRAC
        spec%use_lowsep_otsu         = .false.
        spec%use_otsu_window         = .false.
        spec%use_cluster_rescue      = .true.
        spec%enforce_min_accept_frac = .true.
    end function pool_default_model_spec

    subroutine normalize( self )
        class(cavg_quality_model), intent(inout) :: self
        self%weights = max(0.0, self%weights)
        if( sum(self%weights) > EPS )then
            self%weights = self%weights / sum(self%weights)
        else
            self%weights = CAVG_QUALITY_DEFAULT_WEIGHTS
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
        write(funit,'(A)') '# cluster_cavgs_quality model'
        write(funit,'(A)') 'model_version=2'
        write(funit,'(A,A)') 'name=', trim(self%name)
        write(funit,'(A,A)') 'family=', trim(self%family)
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
        symbol     = model_symbol_from_name(model%name)
        func_name  = trim(symbol)//'_model_spec'
        const_name = 'CAVG_QUALITY_MODEL_'//trim(uppercase(symbol))
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '! cluster_cavgs_quality built-in model promotion snippet'
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
        write(funit,'(A,A)') '        spec%family                  = ', trim(fortran_quote(model%family))
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

    function model_symbol_from_name( name ) result( symbol )
        character(len=*), intent(in) :: name
        character(len=64) :: symbol
        character(len=LONGSTRLEN) :: lower_name
        character(len=1) :: ch
        integer :: i, n
        lower_name = lowercase(adjustl(trim(name)))
        symbol = ''
        n = 0
        do i = 1, len_trim(lower_name)
            ch = lower_name(i:i)
            if( is_model_symbol_char(ch) )then
                call append_symbol_char(symbol, n, ch)
            else if( n > 0 .and. symbol(n:n) /= '_' )then
                call append_symbol_char(symbol, n, '_')
            endif
        end do
        do while( n > 1 .and. symbol(n:n) == '_' )
            symbol(n:n) = ' '
            n = n - 1
        end do
        if( n == 0 )then
            symbol = 'quality_model'
        else if( symbol(1:1) >= '0' .and. symbol(1:1) <= '9' )then
            symbol = 'model_'//trim(symbol)
        endif
    end function model_symbol_from_name

    subroutine append_symbol_char( symbol, n, ch )
        character(len=*), intent(inout) :: symbol
        integer,          intent(inout) :: n
        character(len=1), intent(in)    :: ch
        integer, parameter :: MAX_SYMBOL_LEN = 40
        if( n >= min(len(symbol), MAX_SYMBOL_LEN) ) return
        n = n + 1
        symbol(n:n) = ch
    end subroutine append_symbol_char

    logical function is_model_symbol_char( ch )
        character(len=1), intent(in) :: ch
        is_model_symbol_char = (ch >= 'a' .and. ch <= 'z') .or. (ch >= '0' .and. ch <= '9') .or. ch == '_'
    end function is_model_symbol_char

    function fortran_quote( str ) result( quoted )
        character(len=*), intent(in) :: str
        character(len=XLONGSTRLEN) :: quoted
        integer :: i, n
        quoted = ''
        n = 1
        quoted(n:n) = ''''
        do i = 1, len_trim(str)
            if( n + 2 > len(quoted) ) exit
            if( str(i:i) == '''' )then
                quoted(n+1:n+2) = ''''''
                n = n + 2
            else
                quoted(n+1:n+1) = str(i:i)
                n = n + 1
            endif
        end do
        if( n < len(quoted) )then
            quoted(n+1:n+1) = ''''
        else
            quoted(n:n) = ''''
        endif
    end function fortran_quote

    function fortran_logical( val ) result( literal )
        logical, intent(in) :: val
        character(len=8) :: literal
        if( val )then
            literal = '.true.'
        else
            literal = '.false.'
        endif
    end function fortran_logical

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
                case('preset')
                    cycle
                case('name')
                    self%name = trim(val)
                case('family')
                    self%family = trim(val)
                case('context')
                    self%context = lowercase(trim(val))
                case('feature_policy')
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
            end select
        end do
        close(funit)
        if( trim(self%family) /= 'linear_boundary' ) &
            THROW_HARD('read_model: unsupported model family: '//trim(self%family))
        call assert_valid_model_context(self%context)
        call self%normalize()
    end subroutine read_model

    subroutine read_feature_weights( val, weights )
        character(len=*), intent(in)    :: val
        real,             intent(inout) :: weights(CAVG_QUALITY_NFEATS)
        character(len=LONGSTRLEN) :: field
        integer :: ios, i, nvals
        real    :: parsed(CAVG_QUALITY_NFEATS)
        nvals = csv_weight_count(val)
        if( nvals == CAVG_QUALITY_NFEATS )then
            parsed = 0.0
            do i = 1, CAVG_QUALITY_NFEATS
                field = csv_field(val, i)
                read(field,*,iostat=ios) parsed(i)
                if( ios /= 0 ) THROW_HARD('read_model: failed to parse feature_weights')
            end do
            weights = parsed
            return
        endif
        if( nvals > 1 )then
            if( nvals < CAVG_QUALITY_NFEATS ) THROW_HARD('read_model: feature_weights has too few values')
            THROW_HARD('read_model: feature_weights has too many values')
        endif
        parsed = weights
        read(val,*,iostat=ios) parsed
        if( ios == 0 )then
            weights = parsed
            return
        endif
        THROW_HARD('read_model: failed to parse feature_weights')
    end subroutine read_feature_weights

    integer function csv_weight_count( val )
        character(len=*), intent(in) :: val
        integer :: i
        csv_weight_count = 0
        do i = 1, CAVG_QUALITY_NFEATS + 1
            if( len_trim(csv_field(val, i)) == 0 ) exit
            csv_weight_count = i
        end do
    end function csv_weight_count

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

    subroutine reset_model( self )
        class(cavg_quality_model), intent(inout) :: self
        call self%init_preset(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
    end subroutine reset_model

    subroutine apply_linear_boundary( quality, model )
        type(cavg_quality_result), intent(inout) :: quality
        class(cavg_quality_model), intent(in) :: model
        real,    allocatable :: dmat(:,:), feats_fit(:,:), score_fit(:)
        integer, allocatable :: inds(:), labels_fit(:), medoids_fit(:)
        integer              :: ncls, nfit, i, j, k, good_fit_label, bad_fit_label
        real                 :: d, score1, score2, rescue_threshold, otsu_threshold, otsu_separation
        real                 :: candidate_threshold
        logical              :: dmat_ok, otsu_ok
        if( .not. allocated(quality%features)    ) THROW_HARD('apply_linear_boundary: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('apply_linear_boundary: missing hard-reject mask')
        ncls = size(quality%features, dim=1)
        if( size(quality%features, dim=2) /= CAVG_QUALITY_NFEATS ) THROW_HARD('apply_linear_boundary: invalid feature count')
        if( size(quality%hard_reject) /= ncls ) THROW_HARD('apply_linear_boundary: invalid mask size')
        if( allocated(quality%states)  ) deallocate(quality%states)
        if( allocated(quality%labels)  ) deallocate(quality%labels)
        if( allocated(quality%medoids) ) deallocate(quality%medoids)
        if( allocated(quality%scores)  ) deallocate(quality%scores)
        allocate(quality%states(ncls), quality%labels(ncls), source=0)
        allocate(quality%scores(ncls), source=0.0)
        quality%scores = matmul(quality%features, model%weights)
        where( quality%hard_reject ) quality%scores = -CLIP_Z
        quality%threshold        = 0.0
        quality%raw_threshold    = 0.0
        quality%threshold_margin = 0.0
        quality%separation       = 0.0
        quality%nclust           = 0
        quality%good_label       = 0
        quality%used_threshold   = .false.
        quality%model_name       = model%name
        quality%model_context    = model%context
        nfit = count(.not. quality%hard_reject)
        if( nfit == 0 ) return
        inds = pack([(i, i=1,ncls)], .not. quality%hard_reject)
        if( nfit < 4 ) then
            call accept_fit_as_single_cluster(quality, inds)
            deallocate(inds)
            return
        end if
        allocate(feats_fit(nfit, CAVG_QUALITY_NFEATS), score_fit(nfit))
        do i = 1, nfit
            feats_fit(i,:) = quality%features(inds(i),:)
            score_fit(i)   = quality%scores(inds(i))
        end do
        allocate(dmat(nfit, nfit), source=0.0)
        !$omp parallel do default(shared) private(i,j,d) schedule(static) proc_bind(close)
        do i = 1, nfit - 1
            do j = i + 1, nfit
                d = feature_vector_distance(feats_fit(i,:), feats_fit(j,:), model%weights)
                dmat(i,j) = d
                dmat(j,i) = d
            end do
        end do
        !$omp end parallel do
        call normalize_quality_dmat(dmat, dmat_ok)
        if( .not. dmat_ok ) then
            call accept_fit_as_single_cluster(quality, inds)
            deallocate(dmat, feats_fit, score_fit, inds)
            return
        end if
        quality%nclust = 2
        call cluster_dmat(dmat, 'kmed', quality%nclust, medoids_fit, labels_fit)
        if( .not. two_cluster_result_is_valid(quality%nclust, labels_fit, medoids_fit, nfit) )then
            call accept_fit_as_single_cluster(quality, inds)
            if( allocated(labels_fit)  ) deallocate(labels_fit)
            if( allocated(medoids_fit) ) deallocate(medoids_fit)
            deallocate(dmat, feats_fit, score_fit, inds)
            return
        endif
        score1 = mean_score_for_label(score_fit, labels_fit, 1)
        score2 = mean_score_for_label(score_fit, labels_fit, 2)
        if( score1 > score2 + EPS ) then
            good_fit_label = 1
            bad_fit_label  = 2
        else if( score2 > score1 + EPS )then
            good_fit_label = 2
            bad_fit_label  = 1
        else
            call choose_tied_good_label(score_fit, labels_fit, medoids_fit, good_fit_label, bad_fit_label)
        endif
        quality%separation = abs(score1 - score2)
        quality%raw_threshold = 0.5 * (mean_score_for_label(score_fit, labels_fit, good_fit_label) + &
                                                mean_score_for_label(score_fit, labels_fit, bad_fit_label))
        quality%labels(inds) = labels_fit
        allocate(quality%medoids(size(medoids_fit)))
        do k = 1, size(medoids_fit)
            quality%medoids(k) = inds(medoids_fit(k))
        end do
        call otsu_score_threshold(score_fit, otsu_threshold, otsu_separation, otsu_ok)
        if( quality%separation < model%min_score_separation ) then
            if( model%use_lowsep_otsu .and. otsu_ok .and. otsu_separation >= model%min_score_separation )then
                quality%raw_threshold     = otsu_threshold
                quality%threshold         = otsu_threshold
                quality%threshold_margin  = 0.0
                do i = 1, nfit
                    if( quality%scores(inds(i)) >= quality%threshold ) quality%states(inds(i)) = 1
                end do
                quality%good_label     = good_fit_label
                quality%used_threshold = .true.
            else
                call accept_fit_as_single_cluster(quality, inds)
            endif
        else
            candidate_threshold = quality%raw_threshold - model%boundary_margin
            quality%threshold = candidate_threshold
            if( model%use_otsu_window .and. otsu_ok .and. otsu_separation >= model%min_score_separation .and. &
                otsu_threshold >= quality%raw_threshold + model%otsu_min_offset .and. &
                otsu_threshold <= quality%raw_threshold + model%otsu_max_offset ) quality%threshold = otsu_threshold
            if( model%use_cluster_rescue )then
                rescue_threshold = quality%threshold - model%cluster_rescue_margin
                do i = 1, nfit
                    if( quality%scores(inds(i)) >= quality%threshold .or. &
                       (labels_fit(i) == good_fit_label .and. quality%scores(inds(i)) >= rescue_threshold) ) &
                       quality%states(inds(i)) = 1
                end do
            else
                do i = 1, nfit
                    if( quality%scores(inds(i)) >= quality%threshold ) quality%states(inds(i)) = 1
                end do
            endif
            if( model%enforce_min_accept_frac ) &
                call enforce_min_accept_fraction(quality%scores, quality%hard_reject, quality%states, quality%threshold, &
                                                 model%min_accept_frac)
            quality%threshold_margin = quality%raw_threshold - quality%threshold
            quality%good_label     = good_fit_label
            quality%used_threshold = .true.
        end if
        deallocate(dmat, feats_fit, score_fit, inds, labels_fit, medoids_fit)
    end subroutine apply_linear_boundary

    real function feature_vector_distance( feat1, feat2, weights )
        real, intent(in) :: feat1(:), feat2(:), weights(:)
        integer :: ifeat
        real    :: delta
        if( size(feat1) /= CAVG_QUALITY_NFEATS .or. size(feat2) /= CAVG_QUALITY_NFEATS .or. &
            size(weights) /= CAVG_QUALITY_NFEATS ) &
            THROW_HARD('feature_vector_distance: invalid feature vector size')
        feature_vector_distance = 0.0
        do ifeat = 1, CAVG_QUALITY_NFEATS
            ! Nonzero weights define both the linear score and the feature set
            ! allowed to shape k-medoids. Hard rejections stay outside this
            ! policy and are applied before clustering.
            if( weights(ifeat) <= EPS ) cycle
            delta = feat1(ifeat) - feat2(ifeat)
            feature_vector_distance = feature_vector_distance + delta**2
        end do
        feature_vector_distance = sqrt(feature_vector_distance)
    end function feature_vector_distance

    subroutine accept_fit_as_single_cluster( quality, inds )
        type(cavg_quality_result), intent(inout) :: quality
        integer,                   intent(in)    :: inds(:)
        if( size(inds) == 0 ) return
        quality%states(inds) = 1
        quality%labels(inds) = 1
        if( allocated(quality%medoids) ) deallocate(quality%medoids)
        allocate(quality%medoids(1), source=inds(1))
        quality%threshold        = minval(quality%scores(inds)) - EPS
        quality%raw_threshold    = quality%threshold
        quality%threshold_margin = 0.0
        quality%nclust           = 1
        quality%good_label       = 1
        quality%used_threshold   = .false.
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

    subroutine enforce_min_accept_fraction( scores, hard_reject, states, threshold, min_accept_frac )
        real,    intent(in)    :: scores(:)
        logical, intent(in)    :: hard_reject(:)
        integer, intent(inout) :: states(:)
        real,    intent(inout) :: threshold
        real,    intent(in)    :: min_accept_frac
        real, allocatable :: vals(:)
        integer :: nfit, min_accept, naccepted, i, loc
        real    :: floor_threshold
        if( size(hard_reject) /= size(scores) ) THROW_HARD('enforce_min_accept_fraction: mask size mismatch')
        if( size(states)      /= size(scores) ) THROW_HARD('enforce_min_accept_fraction: state size mismatch')
        nfit = count(.not. hard_reject)
        if( nfit <= 0 ) return
        min_accept = min(nfit, max(1, ceiling(min_accept_frac * real(nfit))))
        naccepted  = count(states > 0 .and. .not. hard_reject)
        if( naccepted >= min_accept ) return
        vals = pack(scores, .not. hard_reject)
        floor_threshold = minval(vals)
        do i = 1, min_accept
            loc = maxloc(vals, dim=1)
            floor_threshold = vals(loc)
            vals(loc) = -huge(1.0)
        end do
        do i = 1, size(scores)
            if( .not. hard_reject(i) .and. scores(i) >= floor_threshold ) states(i) = 1
        end do
        threshold = min(threshold, floor_threshold)
        deallocate(vals)
    end subroutine enforce_min_accept_fraction

    real function mean_score_for_label( scores, labels, label )
        real,    intent(in) :: scores(:)
        integer, intent(in) :: labels(:), label
        integer             :: n
        n = count(labels == label)
        if( n == 0 ) THROW_HARD('mean_score_for_label: empty cluster')
        mean_score_for_label = sum(scores, mask=labels == label) / real(n)
    end function mean_score_for_label

end module simple_cavg_quality_model
