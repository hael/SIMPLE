!@descr: feature inventory and feature extraction for 2D class-average quality
module simple_cavg_quality_feats
!$ use omp_lib
use simple_defs,               only: nthr_glob
use simple_error,              only: simple_exception
use simple_image,              only: image
use simple_image_bin,          only: image_bin
use simple_oris,               only: oris
use simple_segmentation,       only: otsu_img
use simple_stat,               only: median, mad_gau
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, CAVG_QUALITY_CONTEXT_CHUNK, &
    CAVG_QUALITY_CONTEXT_POOL, CAVG_QUALITY_CONTEXT_SIEVE, cavg_quality_feature_def, &
    cavg_quality_result, CAVG_REJECT_REASON_NONE, CAVG_REJECT_REASON_POP_NONPOS, &
    CAVG_REJECT_REASON_POP_LOWFRAC, CAVG_REJECT_REASON_BAD_PIXELS, CAVG_REJECT_REASON_NO_COMPONENT, &
    CAVG_REJECT_REASON_MASK_GEOMETRY, CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW
implicit none
private

public :: CAVG_QUALITY_NFEATS
public :: LOG_EPS
public :: I_LOG_POP
public :: I_NEG_LOG_RES
public :: I_CENTERED
public :: I_LOCVAR_FG
public :: I_LOCVAR_BG
public :: I_CORR_FRC
public :: I_CENTER_EDGE_SNR
public :: I_CC_AREA_FRAC
public :: I_PRESENCE
public :: I_NEG_LOCVAR_FG
public :: I_NEG_LOCVAR_BG
public :: I_LOG_LOCVAR_FG_BG_RATIO
public :: I_BP40_100_CENTER_EDGE_VAR
public :: I_FUZZY_BALL_SIGNAL
public :: cavg_quality_feature_def
public :: cavg_quality_feature_name
public :: cavg_quality_feature_description
public :: cavg_quality_feature_direction
public :: cavg_quality_feature_family
public :: write_cavg_quality_feature_inventory
public :: extract_cavg_quality_features
public :: normalize_cavg_quality_features
public :: CAVG_RES_HARD_REJECT_A
public :: POP_FRACTION_HARD_REJECT
public :: CHUNK_LOCVAR_FG_HARD_REJECT_MAX
#include "simple_local_flags.inc"

real,    parameter :: LOG_EPS                   = 1.0e-12
! Low-pass scales used by the foreground and signal metrics documented in
! doc/microchunk_and_rejection/model_cavgs_rejection.md.
real,    parameter :: FOREGROUND_SEG_LP         = 30.0
real,    parameter :: SIGNAL_METRIC_LP          = 10.0
real,    parameter :: OVERFIT_SIGNAL_BP_HP      = 100.0
real,    parameter :: OVERFIT_SIGNAL_BP_LP      = 40.0
real,    parameter :: CAVG_RES_HARD_REJECT_A    = 40.0
real,    parameter :: POP_FRACTION_HARD_REJECT  = 0.0035
real,    parameter :: BP_CENTER_EDGE_VAR_HARD_REJECT_MIN = 1.5
real,    parameter :: CHUNK_LOCVAR_FG_HARD_REJECT_MAX = exp(-4.5)
real,    parameter :: POOL_RES_HARD_REJECT_A    = 25.0
real,    parameter :: POOL_POP_FRACTION_HARD_REJECT = 5.0e-4
real,    parameter :: POOL_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN = 10.0
integer, parameter :: LOCVAR_WINDOW             = 10
integer, parameter :: MASK_HARD_OUTSIDE_PIXELS  = 10
integer, parameter :: ANALYSIS_BOXSIZE          = 128

integer, parameter :: I_LOG_POP                  = 1
integer, parameter :: I_NEG_LOG_RES              = 2
integer, parameter :: I_CENTERED                 = 3
integer, parameter :: I_LOCVAR_FG                = 4
integer, parameter :: I_LOCVAR_BG                = 5
integer, parameter :: I_CORR_FRC                 = 6
integer, parameter :: I_CENTER_EDGE_SNR          = 7
integer, parameter :: I_CC_AREA_FRAC             = 8
integer, parameter :: I_PRESENCE                 = 9
integer, parameter :: I_NEG_LOCVAR_FG            = 10
integer, parameter :: I_NEG_LOCVAR_BG            = 11
integer, parameter :: I_LOG_LOCVAR_FG_BG_RATIO   = 12
integer, parameter :: I_BP40_100_CENTER_EDGE_VAR = 13
integer, parameter :: I_FUZZY_BALL_SIGNAL        = 14

type(cavg_quality_feature_def), parameter :: FEATURE_DEFS(CAVG_QUALITY_NFEATS) = [ &
    cavg_quality_feature_def('log_pop', 'higher_is_better', &
        'log class population; larger classes have more particle support', 'microchunk'), &
    cavg_quality_feature_def('neg_log_res', 'higher_is_better', &
        'negative log class resolution estimate; higher values indicate better nominal resolution', 'microchunk'), &
    cavg_quality_feature_def('centered', 'higher_is_better', &
        'negative normalized centroid displacement of segmented class-average signal', 'microchunk'), &
    cavg_quality_feature_def('log_locvar_fg', 'higher_is_better', &
        'log local variance measured in the foreground Otsu mask', 'microchunk'), &
    cavg_quality_feature_def('log_locvar_bg', 'higher_is_better', &
        'log local variance measured in the complementary background mask', 'microchunk'), &
    cavg_quality_feature_def('corr_frc_proxy', 'higher_is_better', &
        'stored class correlation or FRC-like quality score when available in cls2D data', 'score'), &
    cavg_quality_feature_def('log_center_edge_snr', 'higher_is_better', &
        'log central signal variance relative to edge variance in the class average', 'signal'), &
    cavg_quality_feature_def('cc_area_frac', 'higher_is_better', &
        'main connected-component area divided by the expected circular mask area', 'microchunk'), &
    cavg_quality_feature_def('presence', 'higher_is_better', &
        'central signal excess relative to border background noise', 'signal'), &
    cavg_quality_feature_def('neg_log_locvar_fg', 'higher_is_better', &
        'negative log local variance in the foreground Otsu mask; higher values indicate lower foreground variation', 'overfit'), &
    cavg_quality_feature_def('neg_log_locvar_bg', 'higher_is_better', &
        'negative log local variance outside the foreground Otsu mask; higher values indicate quieter background', 'overfit'), &
    cavg_quality_feature_def('log_locvar_fg_bg_ratio', 'higher_is_better', &
        'log foreground/background local variance ratio; higher values indicate support-localized detail', 'overfit'), &
    cavg_quality_feature_def('log_bp40_100_center_edge_var', 'higher_is_better', &
        'log center/edge variance ratio after 100-40 A band-pass; low values flag overfit fuzzy balls', 'overfit'), &
    cavg_quality_feature_def('fuzzy_ball_signal', 'higher_is_better', &
        'combined foreground texture, presence, and 100-40 A localization; low values flag fuzzy balls', 'overfit') ]

contains

    function cavg_quality_feature_name( i ) result( name )
        integer, intent(in) :: i
        character(len=32)   :: name
        if( i < 1 .or. i > CAVG_QUALITY_NFEATS ) THROW_HARD('invalid cavg quality feature index')
        name = FEATURE_DEFS(i)%name
    end function cavg_quality_feature_name

    function cavg_quality_feature_description( i ) result( description )
        integer, intent(in) :: i
        character(len=160)  :: description
        if( i < 1 .or. i > CAVG_QUALITY_NFEATS ) THROW_HARD('invalid cavg quality feature index')
        description = FEATURE_DEFS(i)%description
    end function cavg_quality_feature_description

    function cavg_quality_feature_direction( i ) result( direction )
        integer, intent(in) :: i
        character(len=32)   :: direction
        if( i < 1 .or. i > CAVG_QUALITY_NFEATS ) THROW_HARD('invalid cavg quality feature index')
        direction = FEATURE_DEFS(i)%direction
    end function cavg_quality_feature_direction

    function cavg_quality_feature_family( i ) result( family )
        integer, intent(in) :: i
        character(len=32)   :: family
        if( i < 1 .or. i > CAVG_QUALITY_NFEATS ) THROW_HARD('invalid cavg quality feature index')
        family = FEATURE_DEFS(i)%family
    end function cavg_quality_feature_family

    subroutine write_cavg_quality_feature_inventory( funit )
        integer, intent(in) :: funit
        integer :: i
        write(funit,'(A)') '# feature_inventory_header=index,name,direction,family,description'
        do i = 1, CAVG_QUALITY_NFEATS
            write(funit,'(A,I0,A,A,A,A,A,A,A,A)') '# feature_inventory,', i, ',', trim(FEATURE_DEFS(i)%name), ',', &
                trim(FEATURE_DEFS(i)%direction), ',', trim(FEATURE_DEFS(i)%family), ',', &
                trim(FEATURE_DEFS(i)%description)
        end do
    end subroutine write_cavg_quality_feature_inventory

    subroutine extract_cavg_quality_features( imgs, cls_oris, mskdiam, quality, quality_context )
        class(image),              intent(inout) :: imgs(:)
        type(oris),                intent(in)    :: cls_oris
        type(cavg_quality_result), intent(inout) :: quality
        real,                 intent(in)    :: mskdiam
        ! real,    allocatable, intent(inout) :: raw(:,:)
        ! logical, allocatable, intent(inout) :: hard_reject(:)
        character(len=*), optional, intent(in) :: quality_context
        integer, allocatable :: pop(:), disc_area(:)
        real,    allocatable :: res(:), corr(:), corr_in(:)
        integer              :: ncls, i, ldim(3), pop_hard_threshold, pool_pop_hard_threshold
        real                 :: smpd, rad_px
        character(len=32)    :: context
        type(image_bin), allocatable :: bin_img(:), cc_img(:), disc_img(:)
        real, allocatable    :: rmat_cc(:,:,:,:), rmat_disc(:,:,:,:)
        real                 :: centroid_norm, locvar_fg, locvar_bg
        real                 :: cc_area_frac
        real                 :: center_edge_snr, presence_score, bp_center_edge_var
        integer              :: ithr
        logical              :: no_component, mask_hard_reject, bad_pixels
        ncls = size(imgs)
        if( ncls == 0 ) THROW_HARD('extract_cavg_quality_features: no class averages')
        if( cls_oris%get_noris() /= ncls ) THROW_HARD('extract_cavg_quality_features: # cls oris /= # cavgs')
        if( mskdiam <= 0.0 ) THROW_HARD('extract_cavg_quality_features: mskdiam must be positive')
        context = CAVG_QUALITY_CONTEXT_CHUNK
        if( present(quality_context) ) context = trim(quality_context)
        call validate_quality_context(context)
        if( allocated(quality%raw)         ) deallocate(quality%raw)
        if( allocated(quality%hard_reject) ) deallocate(quality%hard_reject)
        if( allocated(quality%reasons)     ) deallocate(quality%reasons)
        allocate(quality%raw(ncls, CAVG_QUALITY_NFEATS), source=0.0)
        allocate(quality%hard_reject(ncls),              source=.false.)
        allocate(quality%reasons(ncls),                  source=CAVG_REJECT_REASON_NONE)
        pop  = cls_oris%get_all_asint('pop')
        res  = cls_oris%get_all('res')
        if( size(pop) /= ncls ) THROW_HARD('extract_cavg_quality_features: invalid pop size')
        if( size(res) /= ncls ) THROW_HARD('extract_cavg_quality_features: invalid res size')
        pop_hard_threshold      = ceiling(real(sum(pop)) * POP_FRACTION_HARD_REJECT)
        pool_pop_hard_threshold = ceiling(real(sum(pop)) * POOL_POP_FRACTION_HARD_REJECT)
        allocate(corr(ncls), source=0.0)
        if( cls_oris%isthere('corr') )then
            corr_in = cls_oris%get_all('corr')
            if( size(corr_in) /= ncls ) THROW_HARD('extract_cavg_quality_features: invalid corr size')
            corr = corr_in
            deallocate(corr_in)
        endif
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( smpd <= 0.0 ) THROW_HARD('extract_cavg_quality_features: non-positive smpd')
        rad_px = (mskdiam / smpd) / 2.0
        call imgs(1)%memoize_mask_coords()
        allocate(bin_img(nthr_glob), cc_img(nthr_glob), disc_img(nthr_glob))
        allocate(disc_area(nthr_glob), source=0)
        allocate(rmat_cc(ldim(1), ldim(2), ldim(3), nthr_glob), &
                 rmat_disc(ldim(1), ldim(2), ldim(3), nthr_glob))
        do ithr = 1, nthr_glob
            call bin_img(ithr)%new_bimg(ldim,  smpd, wthreads=.false.)
            call cc_img(ithr)%new_bimg(ldim,   smpd, wthreads=.false.)
            call disc_img(ithr)%disc(ldim,     smpd, rad_px)
            call disc_img(ithr)%get_rmat_sub(rmat_disc(:,:,:,ithr))
            disc_area(ithr) = count(rmat_disc(:,:,:,ithr) > 0.0)
        end do
        !$omp parallel do default(shared) private(i,ithr,centroid_norm,no_component,&
        !$omp& mask_hard_reject,bad_pixels,locvar_fg,locvar_bg,center_edge_snr,presence_score,bp_center_edge_var,&
        !$omp& cc_area_frac)&
        !$omp proc_bind(close) schedule(static)
        do i = 1, ncls
            ithr = omp_get_thread_num() + 1
            bad_pixels = imgs(i)%contains_nans()
            if( bad_pixels )then
                centroid_norm       = 2.0
                cc_area_frac        = 0.0
                no_component        = .true.
                mask_hard_reject    = .true.
                locvar_fg           = 0.0
                locvar_bg           = 0.0
                center_edge_snr     = 0.0
                presence_score      = 0.0
                bp_center_edge_var  = 0.0
            else
                call measure_cavg_foreground_geometry(imgs(i), bin_img(ithr), cc_img(ithr), rmat_cc(:,:,:,ithr), &
                                                      rmat_disc(:,:,:,ithr), disc_area(ithr), rad_px, &
                                                      centroid_norm, cc_area_frac, no_component, mask_hard_reject)
                call measure_cavg_image_metrics(imgs(i), rad_px, locvar_fg, locvar_bg, center_edge_snr, &
                                                presence_score, bp_center_edge_var)
            endif
            quality%raw(i, I_LOG_POP)                  = log(real(max(pop(i), 0)) + 1.0)
            quality%raw(i, I_NEG_LOG_RES)              = resolution_feature(res(i))
            quality%raw(i, I_CENTERED)                 = -centroid_norm
            quality%raw(i, I_LOCVAR_FG)                = log(max(locvar_fg, LOG_EPS))
            quality%raw(i, I_LOCVAR_BG)                = log(max(locvar_bg, LOG_EPS))
            quality%raw(i, I_CORR_FRC)                 = corr(i)
            quality%raw(i, I_CENTER_EDGE_SNR)          = log(max(center_edge_snr, LOG_EPS))
            quality%raw(i, I_CC_AREA_FRAC)             = cc_area_frac
            quality%raw(i, I_PRESENCE)                 = presence_score
            quality%raw(i, I_NEG_LOCVAR_FG)            = -log(max(locvar_fg, LOG_EPS))
            quality%raw(i, I_NEG_LOCVAR_BG)            = -log(max(locvar_bg, LOG_EPS))
            quality%raw(i, I_LOG_LOCVAR_FG_BG_RATIO)   = log(max(locvar_fg, LOG_EPS)) - log(max(locvar_bg, LOG_EPS))
            quality%raw(i, I_BP40_100_CENTER_EDGE_VAR) = log(max(bp_center_edge_var, LOG_EPS))
            quality%raw(i, I_FUZZY_BALL_SIGNAL)        = quality%raw(i, I_LOCVAR_FG) + quality%raw(i, I_PRESENCE) + &
                                                         quality%raw(i, I_BP40_100_CENTER_EDGE_VAR)
            quality%hard_reject(i) = cavg_hard_reject_for_context(context, pop(i), pop_hard_threshold, &
                pool_pop_hard_threshold, res(i), bad_pixels, no_component, mask_hard_reject, &
                locvar_fg, locvar_bg, bp_center_edge_var, quality%reasons(i))
        end do
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call bin_img(ithr)%kill_bimg()
            call cc_img(ithr)%kill_bimg()
            call disc_img(ithr)%kill_bimg()
        end do
        deallocate(bin_img, cc_img, disc_img, rmat_cc, rmat_disc, disc_area, pop, res, corr)
    end subroutine extract_cavg_quality_features

    subroutine validate_quality_context( quality_context )
        character(len=*), intent(in) :: quality_context
        select case(trim(quality_context))
            case(CAVG_QUALITY_CONTEXT_CHUNK, CAVG_QUALITY_CONTEXT_POOL, CAVG_QUALITY_CONTEXT_SIEVE)
                ! valid contexts
            case DEFAULT
                THROW_HARD('extract_cavg_quality_features: quality_context must be chunk, pool, or sieve')
        end select
    end subroutine validate_quality_context

    logical function cavg_hard_reject_for_context( quality_context, pop, pop_hard_threshold, pool_pop_hard_threshold, &
                                                   res, bad_pixels, no_component, mask_hard_reject, locvar_fg, locvar_bg, &
                                                   bp_center_edge_var, reason )
        character(len=*), intent(in)  :: quality_context
        integer,          intent(in)  :: pop, pop_hard_threshold, pool_pop_hard_threshold
        real,             intent(in)  :: res, locvar_fg, locvar_bg, bp_center_edge_var
        logical,          intent(in)  :: bad_pixels, no_component, mask_hard_reject
        integer,          intent(out) :: reason
        reason = CAVG_REJECT_REASON_NONE
        cavg_hard_reject_for_context = pop <= 0 .or. &
                                       res > CAVG_RES_HARD_REJECT_A .or. &
                                       bad_pixels .or. no_component .or. mask_hard_reject .or. &
                                       (locvar_fg <= EPS .and. locvar_bg <= EPS)
        select case(trim(quality_context))
            case(CAVG_QUALITY_CONTEXT_CHUNK)
                ! Chunk is the first learned-model stage: sieve-cleaned particles are
                ! reclassified in larger 10-30k-particle chunks, then logistic
                ! rejection handles the trainable cases. These gates only remove
                ! obvious non-trainable failures before the model stage.
                cavg_hard_reject_for_context = cavg_hard_reject_for_context .or. &
                    pop < pop_hard_threshold .or. &
                    bp_center_edge_var < BP_CENTER_EDGE_VAR_HARD_REJECT_MIN .or. &
                    locvar_fg < CHUNK_LOCVAR_FG_HARD_REJECT_MAX
            case(CAVG_QUALITY_CONTEXT_SIEVE)
                ! Sieve is the small-chunk screening route. It intentionally uses
                ! conservative hard gates only and does not inherit the shared
                ! chunk/pool validity gates, because no learned model is meant to
                ! rescue or reinterpret these very small 2D runs.
                cavg_hard_reject_for_context = .false.
                if( pop < pop_hard_threshold ) then
                    cavg_hard_reject_for_context = .true.
                    reason = CAVG_REJECT_REASON_POP_LOWFRAC
                else if( bp_center_edge_var < BP_CENTER_EDGE_VAR_HARD_REJECT_MIN ) then
                    cavg_hard_reject_for_context = .true.
                    reason = CAVG_REJECT_REASON_BP_CENTER_EDGE_LOW
                end if
            case(CAVG_QUALITY_CONTEXT_POOL)
                ! Pool is the final pre-3D 2D selection stage after highly cleaned
                ! chunk outputs are merged. These pool-specific gates remove
                ! obvious low-population, weak-localization, or low-resolution
                ! failures before the pool model; subtler cases remain learned.
                cavg_hard_reject_for_context = cavg_hard_reject_for_context .or. &
                    pop < pool_pop_hard_threshold .or. &
                    res > POOL_RES_HARD_REJECT_A .or. &
                    bp_center_edge_var < POOL_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN
        end select
    end function cavg_hard_reject_for_context

    subroutine normalize_cavg_quality_features( raw, hard_reject, features )
        real,                 intent(in)    :: raw(:,:)
        logical,              intent(in)    :: hard_reject(:)
        real,    allocatable, intent(inout) :: features(:,:)
        real, allocatable :: vals(:)
        integer           :: ncls, nfeats, j
        real              :: med, dev
        ncls   = size(raw, dim=1)
        nfeats = size(raw, dim=2)
        if( nfeats /= CAVG_QUALITY_NFEATS ) THROW_HARD('normalize_cavg_quality_features: invalid feature count')
        if( size(hard_reject) /= ncls ) THROW_HARD('normalize_cavg_quality_features: invalid mask size')
        if( allocated(features) ) deallocate(features)
        allocate(features(ncls, nfeats), source=0.0)
        do j = 1, nfeats
            if( count(.not. hard_reject) > 1 ) then
                vals = pack(raw(:,j), .not. hard_reject)
                med  = median(vals)
                dev  = mad_gau(vals, med)
                if( dev > EPS ) then
                    features(:,j) = max(-CLIP_Z, min(CLIP_Z, (raw(:,j) - med) / dev))
                else
                    features(:,j) = 0.0
                end if
                deallocate(vals)
            else
                features(:,j) = 0.0
            end if
            where( hard_reject ) features(:,j) = -CLIP_Z
        end do
    end subroutine normalize_cavg_quality_features

    subroutine measure_cavg_foreground_geometry( img, bin_img, cc_img, rmat_cc, rmat_disc, disc_area, rad_px, &
                                                 centroid_norm, cc_area_frac, no_component, mask_hard_reject )
        class(image),     intent(inout) :: img
        type(image_bin),  intent(inout) :: bin_img, cc_img
        real,             intent(inout) :: rmat_cc(:,:,:)
        real,             intent(in)    :: rmat_disc(:,:,:)
        integer,          intent(in)    :: disc_area
        real,             intent(in)    :: rad_px
        real,             intent(out)   :: centroid_norm, cc_area_frac
        logical,          intent(out)   :: no_component, mask_hard_reject
        real, allocatable :: ccsizes(:)
        integer           :: j, loc, nccs, nccs_valid, area, outside
        real              :: cc_diam, xy(2)
        centroid_norm = 2.0
        cc_area_frac  = 0.0
        nccs_valid    = 0
        no_component  = .false.
        mask_hard_reject = .false.
        call bin_img%copy(img)
        call bin_img%zero_edgeavg()
        call bin_img%bp(0.0, FOREGROUND_SEG_LP)
        call otsu_img(bin_img)
        call bin_img%set_imat()
        call bin_img%find_ccs(cc_img)
        call cc_img%get_nccs(nccs)
        nccs_valid = nccs
        do j = 1, nccs
            call cc_img%diameter_cc(j, cc_diam)
            if( cc_diam > real(size(rmat_cc, dim=1)) )then
                call cc_img%elim_cc(j, update=.false.)
                nccs_valid = nccs_valid - 1
            endif
        end do
        if( nccs_valid <= 0 ) then
            no_component = .true.
            mask_hard_reject = .true.
            return
        end if
        centroid_norm = 0.0
        call cc_img%order_ccs()
        call cc_img%update_img_rmat()
        call cc_img%get_nccs(nccs)
        do j = 1, nccs
            call cc_img%masscen_cc(j, xy)
            centroid_norm = max(centroid_norm, sqrt(xy(1)**2 + xy(2)**2) / max(rad_px, 1.0))
            if( sqrt(xy(1)**2 + xy(2)**2) > rad_px ) mask_hard_reject = .true.
        end do
        ccsizes = cc_img%size_ccs()
        loc     = maxloc(ccsizes, dim=1)
        deallocate(ccsizes)
        call cc_img%cc2bin(loc)
        call cc_img%get_rmat_sub(rmat_cc)
        area    = count(rmat_cc > 0.0)
        outside = count(rmat_cc - rmat_disc > 0.0)
        cc_area_frac = real(area) / real(max(disc_area, 1))
        if( outside > MASK_HARD_OUTSIDE_PIXELS ) mask_hard_reject = .true.
    end subroutine measure_cavg_foreground_geometry

    subroutine measure_cavg_image_metrics( img_src, rad_px, locvar_fg, locvar_bg, center_edge_snr, presence_score, &
                                           bp_center_edge_var )
        class(image), intent(inout) :: img_src
        real,         intent(in)    :: rad_px
        real,         intent(out)   :: locvar_fg, locvar_bg, center_edge_snr
        real,         intent(out)   :: presence_score, bp_center_edge_var
        type(image)       :: img, img_bin, img_bp
        real, allocatable :: bin_mask(:,:,:)
        real              :: bp_center_edge_std
        integer           :: ldim(3), ldim_target(3)
        img  = img_src
        ldim = img%get_ldim()
        call img%zero_edgeavg()
        call img%bp(0.0, SIGNAL_METRIC_LP)
        presence_score = img%presence()
        center_edge_snr = img%center_edge_snr(rad_px)
        img_bp = img
        ldim_target = [ANALYSIS_BOXSIZE, ANALYSIS_BOXSIZE, ldim(3)]
        call img_bp%fft()
        if( ldim(1) > ldim_target(1) .or. ldim(2) > ldim_target(2) )then
            call img_bp%clip_inplace(ldim_target)
        else if( ldim(1) < ldim_target(1) .or. ldim(2) < ldim_target(2) )then
            call img_bp%pad_inplace(ldim_target)
        end if
        call img_bp%ifft()
        call img_bp%set_smpd(img_bp%get_smpd() * real(ldim(1)) / real(ldim_target(1)))
        call img_bp%bp(OVERFIT_SIGNAL_BP_HP, OVERFIT_SIGNAL_BP_LP)
        bp_center_edge_std = img_bp%center_edge_snr((real(ldim_target(1)) * 0.4))
        bp_center_edge_var = bp_center_edge_std * bp_center_edge_std
        img_bin = img
        call otsu_img(img_bin)
        allocate(bin_mask(ldim(1), ldim(2), ldim(3)))
        call img_bin%get_rmat_sub(bin_mask)
        call img%loc_var_masked(bin_mask(:,:,1), LOCVAR_WINDOW, locvar_fg, locvar_bg)
        deallocate(bin_mask)
        call img%kill()
        call img_bin%kill()
        call img_bp%kill()
    end subroutine measure_cavg_image_metrics

    real function resolution_feature( res )
        real, intent(in) :: res
        if( res > EPS ) then
            resolution_feature = -log(res)
        else
            resolution_feature = -log(1000.0)
        end if
    end function resolution_feature

end module simple_cavg_quality_feats
