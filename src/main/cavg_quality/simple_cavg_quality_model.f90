!@descr: instantiable class-average quality decision model
module simple_cavg_quality_model
use simple_defs,               only: LONGSTRLEN, XLONGSTRLEN
use simple_error,              only: simple_exception
use simple_string_utils,       only: str_is_true
use simple_clustering_utils,   only: cluster_dmat
use simple_cavg_quality_feats, only: EPS, CLIP_Z, normalize_quality_dmat
use simple_cavg_quality_types, only: CAVG_REJECTION_CHUNK, CAVG_REJECTION_POOL, CAVG_QUALITY_NFEATS, &
    cavg_quality_model_spec, cavg_quality_result
implicit none
private
#include "simple_local_flags.inc"

public :: CAVG_REJECTION_CHUNK
public :: CAVG_REJECTION_POOL
public :: CAVG_QUALITY_DEFAULT_WEIGHTS
public :: CAVG_QUALITY_DEFAULT_HIST_DMAT_WEIGHT
public :: CAVG_QUALITY_MODEL_DEFAULT
public :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT
public :: CAVG_QUALITY_MODEL_POOL_DEFAULT
public :: cavg_quality_model
public :: cavg_quality_model_spec
public :: cavg_quality_model_builtin_names
public :: cavg_quality_model_default_name
public :: cavg_rejection_type_from_name
public :: cavg_rejection_type_name
public :: cluster_cavg_quality
public :: normalize_cavg_quality_model

type :: cavg_quality_model_catalog_entry
    character(len=64) :: name           = ''
    integer           :: rejection_type = CAVG_REJECTION_CHUNK
end type cavg_quality_model_catalog_entry

! Built-in presets are complete model specifications. To promote a learned
! model into the code, add a named preset and include it in the catalog.
character(len=*), parameter :: CAVG_QUALITY_MODEL_DEFAULT       = 'default'
character(len=*), parameter :: CAVG_QUALITY_MODEL_CHUNK_DEFAULT = 'chunk_default_v1'
character(len=*), parameter :: CAVG_QUALITY_MODEL_POOL_DEFAULT  = 'pool_default_v1'

integer, parameter :: CAVG_QUALITY_NBUILTIN_MODELS = 2
type(cavg_quality_model_catalog_entry), parameter :: CAVG_QUALITY_BUILTIN_MODELS(CAVG_QUALITY_NBUILTIN_MODELS) = [ &
    cavg_quality_model_catalog_entry(CAVG_QUALITY_MODEL_CHUNK_DEFAULT, CAVG_REJECTION_CHUNK), &
    cavg_quality_model_catalog_entry(CAVG_QUALITY_MODEL_POOL_DEFAULT,  CAVG_REJECTION_POOL) ]

real, parameter :: MIN_SCORE_SEPARATION    = 0.15
real, parameter :: BOUNDARY_MARGIN_DEFAULT = 0.05
real, parameter :: CHUNK_BOUNDARY_OFFSET   = 0.25
real, parameter :: CHUNK_OTSU_MIN_OFFSET   = 0.25
real, parameter :: CHUNK_OTSU_MAX_OFFSET   = 0.50
real, parameter :: POOL_MIN_ACCEPT_FRAC    = 0.65
real, parameter :: HIST_DMAT_WEIGHT        = 0.50
real, parameter :: CLUSTER_RESCUE_MARGIN   = 0.20

! Default interpretable linear score coefficients. These are model parameters,
! while the feature definitions and extraction live in simple_cavg_quality_feats.
real, parameter :: CAVG_QUALITY_DEFAULT_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    0.34, 0.18, 0.00, 0.00, 0.16, 0.18, 0.00, 0.00, 0.00, 0.00, 0.00, 0.14 ]
real, parameter :: CAVG_QUALITY_DEFAULT_HIST_DMAT_WEIGHT = HIST_DMAT_WEIGHT

type :: cavg_quality_model
    character(len=64) :: name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
    character(len=32) :: family                  = 'linear_boundary'
    character(len=32) :: context                 = 'chunk'
    integer           :: rejection_type          = CAVG_REJECTION_CHUNK
    real              :: weights(CAVG_QUALITY_NFEATS) = CAVG_QUALITY_DEFAULT_WEIGHTS
    real              :: boundary_margin         = -CHUNK_BOUNDARY_OFFSET
    real              :: min_score_separation    = MIN_SCORE_SEPARATION
    real              :: hist_dmat_weight        = HIST_DMAT_WEIGHT
    real              :: otsu_min_offset         = CHUNK_OTSU_MIN_OFFSET
    real              :: otsu_max_offset         = CHUNK_OTSU_MAX_OFFSET
    real              :: cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
    real              :: min_accept_frac         = 0.0
    logical           :: use_lowsep_otsu         = .true.
    logical           :: use_otsu_window         = .true.
    logical           :: use_cluster_rescue      = .false.
    logical           :: enforce_min_accept_frac = .false.
contains
    procedure :: init_default => cavg_quality_model_init_default
    procedure :: init_preset  => cavg_quality_model_init_preset
    procedure :: init_spec    => cavg_quality_model_init_spec
    procedure :: get_spec     => cavg_quality_model_get_spec
    procedure :: normalize    => cavg_quality_model_normalize
    procedure :: classify     => cavg_quality_model_classify
    procedure :: write        => cavg_quality_model_write
    procedure :: read         => cavg_quality_model_read
    procedure :: kill         => cavg_quality_model_kill
end type cavg_quality_model

contains

    subroutine cavg_quality_model_init_default( self, rejection_type )
        class(cavg_quality_model), intent(inout) :: self
        integer,                   intent(in)    :: rejection_type
        call assert_valid_rejection_type(rejection_type)
        call self%init_preset(cavg_quality_model_default_name(rejection_type))
    end subroutine cavg_quality_model_init_default

    subroutine cavg_quality_model_init_preset( self, preset_name )
        class(cavg_quality_model), intent(inout) :: self
        character(len=*),          intent(in)    :: preset_name
        type(cavg_quality_model_spec) :: spec
        spec = cavg_quality_model_builtin_spec(preset_name)
        call self%init_spec(spec)
    end subroutine cavg_quality_model_init_preset

    function cavg_quality_model_builtin_spec( preset_name ) result( spec )
        character(len=*), intent(in) :: preset_name
        type(cavg_quality_model_spec) :: spec
        character(len=LONGSTRLEN) :: errmsg
        select case(trim(preset_name))
            case(CAVG_QUALITY_MODEL_CHUNK_DEFAULT)
                spec = chunk_default_model_spec()
            case(CAVG_QUALITY_MODEL_POOL_DEFAULT)
                spec = pool_default_model_spec()
            case DEFAULT
                errmsg = 'unknown class-average quality model preset: '//trim(preset_name)//&
                         '; available presets: '//trim(cavg_quality_model_builtin_names())
                THROW_HARD(trim(errmsg))
        end select
    end function cavg_quality_model_builtin_spec

    function cavg_quality_model_builtin_names() result( names )
        character(len=LONGSTRLEN) :: names
        integer :: i
        names = ''
        do i = 1, CAVG_QUALITY_NBUILTIN_MODELS
            if( i > 1 ) names = trim(names)//'|'
            names = trim(names)//trim(CAVG_QUALITY_BUILTIN_MODELS(i)%name)
        end do
    end function cavg_quality_model_builtin_names

    function cavg_quality_model_default_name( rejection_type ) result( name )
        integer, intent(in) :: rejection_type
        character(len=64) :: name
        call assert_valid_rejection_type(rejection_type)
        select case(rejection_type)
            case(CAVG_REJECTION_CHUNK)
                name = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
            case(CAVG_REJECTION_POOL)
                name = CAVG_QUALITY_MODEL_POOL_DEFAULT
        end select
    end function cavg_quality_model_default_name

    subroutine cavg_quality_model_init_spec( self, spec )
        class(cavg_quality_model), intent(inout) :: self
        type(cavg_quality_model_spec), intent(in) :: spec
        call assert_valid_rejection_type(spec%rejection_type)
        if( trim(spec%family) /= 'linear_boundary' ) &
            THROW_HARD('unsupported class-average quality model family: '//trim(spec%family))
        self%name                    = trim(spec%name)
        self%family                  = trim(spec%family)
        self%context                 = trim(spec%context)
        self%rejection_type          = spec%rejection_type
        self%weights                 = spec%weights
        self%boundary_margin         = spec%boundary_margin
        self%min_score_separation    = spec%min_score_separation
        self%hist_dmat_weight        = spec%hist_dmat_weight
        self%otsu_min_offset         = spec%otsu_min_offset
        self%otsu_max_offset         = spec%otsu_max_offset
        self%cluster_rescue_margin   = spec%cluster_rescue_margin
        self%min_accept_frac         = spec%min_accept_frac
        self%use_lowsep_otsu         = spec%use_lowsep_otsu
        self%use_otsu_window         = spec%use_otsu_window
        self%use_cluster_rescue      = spec%use_cluster_rescue
        self%enforce_min_accept_frac = spec%enforce_min_accept_frac
        call self%normalize()
    end subroutine cavg_quality_model_init_spec

    function cavg_quality_model_get_spec( self ) result( spec )
        class(cavg_quality_model), intent(in) :: self
        type(cavg_quality_model_spec) :: spec
        spec%name                    = self%name
        spec%family                  = self%family
        spec%context                 = self%context
        spec%rejection_type          = self%rejection_type
        spec%weights                 = self%weights
        spec%boundary_margin         = self%boundary_margin
        spec%min_score_separation    = self%min_score_separation
        spec%hist_dmat_weight        = self%hist_dmat_weight
        spec%otsu_min_offset         = self%otsu_min_offset
        spec%otsu_max_offset         = self%otsu_max_offset
        spec%cluster_rescue_margin   = self%cluster_rescue_margin
        spec%min_accept_frac         = self%min_accept_frac
        spec%use_lowsep_otsu         = self%use_lowsep_otsu
        spec%use_otsu_window         = self%use_otsu_window
        spec%use_cluster_rescue      = self%use_cluster_rescue
        spec%enforce_min_accept_frac = self%enforce_min_accept_frac
    end function cavg_quality_model_get_spec

    function chunk_default_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_CHUNK_DEFAULT
        spec%family                  = 'linear_boundary'
        spec%context                 = 'chunk'
        spec%rejection_type          = CAVG_REJECTION_CHUNK
        spec%weights                 = CAVG_QUALITY_DEFAULT_WEIGHTS
        spec%boundary_margin         = -CHUNK_BOUNDARY_OFFSET
        spec%min_score_separation    = MIN_SCORE_SEPARATION
        spec%hist_dmat_weight        = HIST_DMAT_WEIGHT
        spec%otsu_min_offset         = CHUNK_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = 0.0
        spec%use_lowsep_otsu         = .true.
        spec%use_otsu_window         = .true.
        spec%use_cluster_rescue      = .false.
        spec%enforce_min_accept_frac = .false.
    end function chunk_default_model_spec

    function pool_default_model_spec() result( spec )
        type(cavg_quality_model_spec) :: spec
        spec%name                    = CAVG_QUALITY_MODEL_POOL_DEFAULT
        spec%family                  = 'linear_boundary'
        spec%context                 = 'pool'
        spec%rejection_type          = CAVG_REJECTION_POOL
        spec%weights                 = CAVG_QUALITY_DEFAULT_WEIGHTS
        spec%boundary_margin         = BOUNDARY_MARGIN_DEFAULT
        spec%min_score_separation    = MIN_SCORE_SEPARATION
        spec%hist_dmat_weight        = HIST_DMAT_WEIGHT
        spec%otsu_min_offset         = CHUNK_OTSU_MIN_OFFSET
        spec%otsu_max_offset         = CHUNK_OTSU_MAX_OFFSET
        spec%cluster_rescue_margin   = CLUSTER_RESCUE_MARGIN
        spec%min_accept_frac         = POOL_MIN_ACCEPT_FRAC
        spec%use_lowsep_otsu         = .false.
        spec%use_otsu_window         = .false.
        spec%use_cluster_rescue      = .true.
        spec%enforce_min_accept_frac = .true.
    end function pool_default_model_spec

    subroutine cavg_quality_model_normalize( self )
        class(cavg_quality_model), intent(inout) :: self
        self%weights = max(0.0, self%weights)
        if( sum(self%weights) > EPS )then
            self%weights = self%weights / sum(self%weights)
        else
            self%weights = CAVG_QUALITY_DEFAULT_WEIGHTS
        endif
        self%hist_dmat_weight = max(0.0, min(1.0, self%hist_dmat_weight))
    end subroutine cavg_quality_model_normalize

    subroutine normalize_cavg_quality_model( model )
        type(cavg_quality_model), intent(inout) :: model
        call model%normalize()
    end subroutine normalize_cavg_quality_model

    subroutine cavg_quality_model_classify( self, quality )
        class(cavg_quality_model), intent(in)    :: self
        type(cavg_quality_result), intent(inout) :: quality
        if( .not. allocated(quality%features)    ) THROW_HARD('cavg_quality_model_classify: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('cavg_quality_model_classify: missing hard-reject mask')
        if( .not. allocated(quality%hist_dmat)   ) THROW_HARD('cavg_quality_model_classify: missing histogram distance matrix')
        quality%rejection_type = self%rejection_type
        quality%model_name     = self%name
        quality%model_context  = self%context
        call cluster_cavg_quality(quality, self)
    end subroutine cavg_quality_model_classify

    subroutine cavg_quality_model_write( self, fname )
        class(cavg_quality_model), intent(in) :: self
        character(len=*),          intent(in) :: fname
        integer :: funit, i
        open(newunit=funit, file=trim(fname), status='replace', action='write')
        write(funit,'(A)') '# cluster_cavgs_quality model'
        write(funit,'(A)') 'model_version=1'
        write(funit,'(A,A)') 'name=', trim(self%name)
        write(funit,'(A,A)') 'family=', trim(self%family)
        write(funit,'(A,A)') 'context=', trim(self%context)
        write(funit,'(A,I0)') 'rejection_type=', self%rejection_type
        write(funit,'(A)', advance='no') 'feature_weights='
        do i = 1, CAVG_QUALITY_NFEATS
            if( i > 1 ) write(funit,'(A)', advance='no') ','
            write(funit,'(ES14.6)', advance='no') self%weights(i)
        end do
        write(funit,*)
        write(funit,'(A,ES14.6)') 'boundary_margin=', self%boundary_margin
        write(funit,'(A,ES14.6)') 'min_score_separation=', self%min_score_separation
        write(funit,'(A,ES14.6)') 'hist_dmat_weight=', self%hist_dmat_weight
        write(funit,'(A,ES14.6)') 'otsu_min_offset=', self%otsu_min_offset
        write(funit,'(A,ES14.6)') 'otsu_max_offset=', self%otsu_max_offset
        write(funit,'(A,ES14.6)') 'cluster_rescue_margin=', self%cluster_rescue_margin
        write(funit,'(A,ES14.6)') 'min_accept_frac=', self%min_accept_frac
        write(funit,'(A,L1)') 'use_lowsep_otsu=', self%use_lowsep_otsu
        write(funit,'(A,L1)') 'use_otsu_window=', self%use_otsu_window
        write(funit,'(A,L1)') 'use_cluster_rescue=', self%use_cluster_rescue
        write(funit,'(A,L1)') 'enforce_min_accept_frac=', self%enforce_min_accept_frac
        close(funit)
    end subroutine cavg_quality_model_write

    subroutine cavg_quality_model_read( self, fname )
        class(cavg_quality_model), intent(inout) :: self
        character(len=*),          intent(in)    :: fname
        character(len=XLONGSTRLEN) :: line
        character(len=LONGSTRLEN)  :: key, val
        integer :: funit, ios, ieq, rej_type
        call self%init_default(CAVG_REJECTION_CHUNK)
        open(newunit=funit, file=trim(fname), status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('cavg_quality_model_read: failed to open '//trim(fname))
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            line = adjustl(line)
            if( len_trim(line) == 0 ) cycle
            if( line(1:1) == '#' ) cycle
            ieq = index(line, '=')
            if( ieq <= 1 ) cycle
            key = adjustl(trim(line(1:ieq-1)))
            val = adjustl(trim(line(ieq+1:)))
            select case(trim(key))
                case('preset')
                    call self%init_preset(trim(val))
                case('name')
                    self%name = trim(val)
                case('family')
                    self%family = trim(val)
                case('context')
                    self%context = trim(val)
                    self%rejection_type = cavg_rejection_type_from_name(trim(val))
                case('rejection_type')
                    read(val,*,iostat=ios) rej_type
                    if( ios == 0 ) self%rejection_type = rej_type
                case('feature_weights')
                    call read_feature_weights(val, self%weights)
                case('boundary_margin')
                    read(val,*,iostat=ios) self%boundary_margin
                case('min_score_separation')
                    read(val,*,iostat=ios) self%min_score_separation
                case('hist_dmat_weight')
                    read(val,*,iostat=ios) self%hist_dmat_weight
                case('otsu_min_offset')
                    read(val,*,iostat=ios) self%otsu_min_offset
                case('otsu_max_offset')
                    read(val,*,iostat=ios) self%otsu_max_offset
                case('cluster_rescue_margin')
                    read(val,*,iostat=ios) self%cluster_rescue_margin
                case('min_accept_frac')
                    read(val,*,iostat=ios) self%min_accept_frac
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
        call assert_valid_rejection_type(self%rejection_type)
        self%context = cavg_rejection_type_name(self%rejection_type)
        call self%normalize()
    end subroutine cavg_quality_model_read

    subroutine read_feature_weights( val, weights )
        character(len=*), intent(in)    :: val
        real,             intent(inout) :: weights(CAVG_QUALITY_NFEATS)
        character(len=LONGSTRLEN) :: field
        integer :: ios, i
        real    :: parsed(CAVG_QUALITY_NFEATS)
        parsed = weights
        read(val,*,iostat=ios) parsed
        if( ios == 0 )then
            weights = parsed
            return
        endif
        parsed = weights
        do i = 1, CAVG_QUALITY_NFEATS
            field = csv_field(val, i)
            if( len_trim(field) == 0 ) exit
            read(field,*,iostat=ios) parsed(i)
            if( ios /= 0 ) THROW_HARD('cavg_quality_model_read: failed to parse feature_weights')
        end do
        weights = parsed
    end subroutine read_feature_weights

    function csv_field( line, ifield ) result( field )
        character(len=*), intent(in) :: line
        integer,          intent(in) :: ifield
        character(len=LONGSTRLEN) :: field
        integer :: i, start, finish, nfield, llen
        field = ''
        start = 1
        nfield = 1
        llen = len_trim(line)
        do i = 1, llen + 1
            if( i == llen + 1 .or. line(i:i) == ',' )then
                finish = i - 1
                if( nfield == ifield )then
                    if( finish >= start ) field = adjustl(trim(line(start:finish)))
                    return
                endif
                nfield = nfield + 1
                start = i + 1
            endif
        end do
    end function csv_field

    subroutine cavg_quality_model_kill( self )
        class(cavg_quality_model), intent(inout) :: self
        call self%init_default(CAVG_REJECTION_CHUNK)
    end subroutine cavg_quality_model_kill

    subroutine cluster_cavg_quality( quality, model )
        type(cavg_quality_result), intent(inout) :: quality
        class(cavg_quality_model), intent(in) :: model
        real,    allocatable :: dmat(:,:), dmat_hist_fit(:,:), feats_fit(:,:), score_fit(:)
        integer, allocatable :: inds(:), labels_fit(:), medoids_fit(:)
        integer              :: ncls, nfit, i, j, k, good_fit_label, bad_fit_label
        real                 :: d, score1, score2, rescue_threshold, otsu_threshold, otsu_separation
        real                 :: candidate_threshold
        logical              :: dmat_ok, hist_ok, otsu_ok
        if( .not. allocated(quality%features)    ) THROW_HARD('cluster_cavg_quality: missing features')
        if( .not. allocated(quality%hard_reject) ) THROW_HARD('cluster_cavg_quality: missing hard-reject mask')
        if( .not. allocated(quality%hist_dmat)   ) THROW_HARD('cluster_cavg_quality: missing histogram distance matrix')
        ncls = size(quality%features, dim=1)
        if( size(quality%features, dim=2) /= CAVG_QUALITY_NFEATS ) THROW_HARD('cluster_cavg_quality: invalid feature count')
        if( size(quality%hard_reject) /= ncls ) THROW_HARD('cluster_cavg_quality: invalid mask size')
        call assert_valid_rejection_type(model%rejection_type)
        if( size(quality%hist_dmat, dim=1) /= ncls .or. size(quality%hist_dmat, dim=2) /= ncls ) &
            THROW_HARD('cluster_cavg_quality: invalid histogram distance matrix size')
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
        quality%rejection_type   = model%rejection_type
        quality%used_threshold   = .false.
        quality%model_name       = model%name
        quality%model_context    = model%context
        nfit = count(.not. quality%hard_reject)
        if( nfit == 0 ) return
        inds = pack([(i, i=1,ncls)], .not. quality%hard_reject)
        if( nfit < 4 ) then
            quality%states(inds) = 1
            quality%labels(inds) = 1
            allocate(quality%medoids(1), source=inds(1))
            quality%threshold        = minval(quality%scores(inds)) - EPS
            quality%raw_threshold    = quality%threshold
            quality%threshold_margin = 0.0
            quality%nclust           = 1
            quality%good_label       = 1
            return
        end if
        allocate(feats_fit(nfit, CAVG_QUALITY_NFEATS), score_fit(nfit))
        do i = 1, nfit
            feats_fit(i,:) = quality%features(inds(i),:)
            score_fit(i)   = quality%scores(inds(i))
        end do
        allocate(dmat(nfit, nfit), source=0.0)
        do i = 1, nfit - 1
            do j = i + 1, nfit
                d = sqrt(sum((feats_fit(i,:) - feats_fit(j,:))**2))
                dmat(i,j) = d
                dmat(j,i) = d
            end do
        end do
        call normalize_quality_dmat(dmat, dmat_ok)
        if( .not. dmat_ok ) then
            quality%states(inds) = 1
            quality%labels(inds) = 1
            allocate(quality%medoids(1), source=inds(1))
            quality%threshold        = minval(quality%scores(inds)) - EPS
            quality%raw_threshold    = quality%threshold
            quality%threshold_margin = 0.0
            quality%nclust           = 1
            quality%good_label       = 1
            deallocate(dmat, feats_fit, score_fit, inds)
            return
        end if
        if( model%hist_dmat_weight > EPS )then
            allocate(dmat_hist_fit(nfit, nfit), source=0.0)
            do i = 1, nfit
                do j = 1, nfit
                    dmat_hist_fit(i,j) = quality%hist_dmat(inds(i), inds(j))
                end do
            end do
            call normalize_quality_dmat(dmat_hist_fit, hist_ok)
            if( hist_ok )then
                dmat = (1.0 - model%hist_dmat_weight) * dmat + model%hist_dmat_weight * dmat_hist_fit
                call normalize_quality_dmat(dmat, dmat_ok)
            endif
            deallocate(dmat_hist_fit)
        endif
        quality%nclust = 2
        call cluster_dmat(dmat, 'kmed', quality%nclust, medoids_fit, labels_fit)
        if( quality%nclust /= 2 ) THROW_HARD('cluster_cavg_quality: expected two k-medoids clusters')
        score1 = mean_score_for_label(score_fit, labels_fit, 1)
        score2 = mean_score_for_label(score_fit, labels_fit, 2)
        if( score1 >= score2 ) then
            good_fit_label = 1
            bad_fit_label  = 2
        else
            good_fit_label = 2
            bad_fit_label  = 1
        end if
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
                quality%states(inds)     = 1
                quality%labels(inds)     = 1
                quality%medoids          = [inds(1)]
                quality%nclust           = 1
                quality%good_label       = 1
                quality%threshold        = minval(quality%scores(inds)) - EPS
                quality%raw_threshold    = quality%threshold
                quality%threshold_margin = 0.0
                quality%used_threshold   = .false.
            endif
        else
            candidate_threshold = quality%raw_threshold - model%boundary_margin
            quality%threshold = candidate_threshold
            if( model%use_otsu_window .and. otsu_ok .and. otsu_separation >= model%min_score_separation .and. &
                otsu_threshold >= quality%raw_threshold + model%otsu_min_offset .and. &
                otsu_threshold <= quality%raw_threshold + model%otsu_max_offset ) quality%threshold = otsu_threshold
            quality%threshold_margin = quality%raw_threshold - quality%threshold
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
    end subroutine cluster_cavg_quality

    integer function cavg_rejection_type_from_name( name )
        character(len=*), intent(in) :: name
        select case(trim(name))
            case('chunk')
                cavg_rejection_type_from_name = CAVG_REJECTION_CHUNK
            case('pool')
                cavg_rejection_type_from_name = CAVG_REJECTION_POOL
            case DEFAULT
                THROW_HARD('invalid class-average rejection type name: '//trim(name))
        end select
    end function cavg_rejection_type_from_name

    function cavg_rejection_type_name( rejection_type ) result( name )
        integer, intent(in) :: rejection_type
        character(len=16)   :: name
        select case(rejection_type)
            case(CAVG_REJECTION_CHUNK)
                name = 'chunk'
            case(CAVG_REJECTION_POOL)
                name = 'pool'
            case DEFAULT
                THROW_HARD('invalid class-average rejection type')
        end select
    end function cavg_rejection_type_name

    subroutine assert_valid_rejection_type( rejection_type )
        integer, intent(in) :: rejection_type
        select case(rejection_type)
            case(CAVG_REJECTION_CHUNK, CAVG_REJECTION_POOL)
                return
            case DEFAULT
                THROW_HARD('invalid class-average rejection type')
        end select
    end subroutine assert_valid_rejection_type

    subroutine otsu_score_threshold( scores, threshold, separation, ok )
        real,    intent(in)  :: scores(:)
        real,    intent(out) :: threshold, separation
        logical, intent(out) :: ok
        integer :: i, n, nlo, nhi
        real    :: candidate, mean_lo, mean_hi, between, best_between
        n = size(scores)
        threshold  = minval(scores) - EPS
        separation = 0.0
        ok         = .false.
        if( n < 4 ) return
        best_between = -huge(1.0)
        do i = 1, n
            candidate = scores(i)
            nlo = count(scores <  candidate)
            nhi = count(scores >= candidate)
            if( nlo == 0 .or. nhi == 0 ) cycle
            mean_lo = sum(scores, mask=scores <  candidate) / real(nlo)
            mean_hi = sum(scores, mask=scores >= candidate) / real(nhi)
            between = real(nlo) * real(nhi) * (mean_hi - mean_lo)**2
            if( between > best_between )then
                best_between = between
                threshold    = candidate
                separation   = mean_hi - mean_lo
                ok           = .true.
            endif
        end do
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
