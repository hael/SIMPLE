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
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, CAVG_QUALITY_TARGET_OVERFIT, &
    cavg_quality_feature_def
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
public :: I_LOG_INT_BOUNDARY_TEX
public :: I_INT_TEX_COVERAGE
public :: I_TEX_EFFECTIVE_AREA
public :: I_NEG_LOCVAR_FG
public :: I_NEG_LOCVAR_BG
public :: I_LOG_LOCVAR_FG_BG_RATIO
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
public :: CAVG_OVERFIT_SCORE_REJECT_THRESHOLD
public :: CAVG_OVERFIT_SCORE_REJECT_WEIGHTS
public :: cavg_overfit_score

#include "simple_local_flags.inc"

real,    parameter :: LOG_EPS                   = 1.0e-12
! Low-pass scales used by the foreground and signal metrics documented in
! doc/microchunk_and_rejection/model_cavgs_rejection.md.
real,    parameter :: FOREGROUND_SEG_LP         = 30.0
real,    parameter :: SIGNAL_METRIC_LP          = 10.0
real,    parameter :: TEXTURE_METRIC_HP         = 30.0
real,    parameter :: TEXTURE_METRIC_LP         = 10.0
real,    parameter :: TEXTURE_BG_MAD_SIGMA      = 3.0
real,    parameter :: CAVG_RES_HARD_REJECT_A    = 40.0
real,    parameter :: POP_FRACTION_HARD_REJECT  = 0.0035
integer, parameter :: LOCVAR_WINDOW             = 10
integer, parameter :: MASK_HARD_OUTSIDE_PIXELS  = 10
integer, parameter :: TEXTURE_RING_WIDTH        = 2

integer, parameter :: I_LOG_POP                 = 1
integer, parameter :: I_NEG_LOG_RES             = 2
integer, parameter :: I_CENTERED                = 3
integer, parameter :: I_LOCVAR_FG               = 4
integer, parameter :: I_LOCVAR_BG               = 5
integer, parameter :: I_CORR_FRC                = 6
integer, parameter :: I_CENTER_EDGE_SNR         = 7
integer, parameter :: I_CC_AREA_FRAC            = 8
integer, parameter :: I_PRESENCE                = 9
integer, parameter :: I_LOG_INT_BOUNDARY_TEX    = 10
integer, parameter :: I_INT_TEX_COVERAGE        = 11
integer, parameter :: I_TEX_EFFECTIVE_AREA      = 12
integer, parameter :: I_NEG_LOCVAR_FG           = 13
integer, parameter :: I_NEG_LOCVAR_BG           = 14
integer, parameter :: I_LOG_LOCVAR_FG_BG_RATIO  = 15

! Optional hard-gate score for rejecting overfitted fuzzy-ball class averages.
! This is deliberately not a learned model at application time. It is the
! promoted four-feature linear score from overfit_default_v1, evaluated on the
! dataset-normalized z-features produced after the standard quality hard gates.
!
! The score intentionally excludes resolution, population, and centering:
!   score =
!       0.1930057 * z_log_center_edge_snr
!     + 0.2138741 * z_cc_area_frac
!     + 0.4490834 * z_presence
!     + 0.1440368 * z_log_locvar_fg_bg_ratio
!
! Classes with score < CAVG_OVERFIT_SCORE_REJECT_THRESHOLD are rejected when
! overfit_score_reject=yes. Joe/microchunking can reuse this block directly:
! compute the same normalized features, call cavg_overfit_score, and compare
! against the threshold below.
real, parameter :: CAVG_OVERFIT_SCORE_REJECT_THRESHOLD = 0.0
real, parameter :: CAVG_OVERFIT_SCORE_REJECT_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, &
    0.0000000, 0.1930057, 0.2138741, 0.4490834, 0.0000000, &
    0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1440368 ]

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
    cavg_quality_feature_def('log_int_boundary_tex', 'higher_is_better', &
        'log medium-frequency texture energy in the foreground interior relative to its boundary ring', 'texture'), &
    cavg_quality_feature_def('int_tex_coverage', 'higher_is_better', &
        'fraction of eroded foreground interior with medium-frequency texture above robust background', 'texture'), &
    cavg_quality_feature_def('tex_effective_area', 'higher_is_better', &
        'foreground fraction effectively occupied by above-background medium-frequency texture', 'texture'), &
    cavg_quality_feature_def('neg_log_locvar_fg', 'higher_is_better', &
        'negative log local variance in the foreground Otsu mask; higher values indicate less overfit-like texture', 'overfit'), &
    cavg_quality_feature_def('neg_log_locvar_bg', 'higher_is_better', &
        'negative log local variance outside the foreground Otsu mask; higher values indicate quieter background', 'overfit'), &
    cavg_quality_feature_def('log_locvar_fg_bg_ratio', 'higher_is_better', &
        'log foreground/background local variance ratio; higher values indicate support-localized detail', 'overfit') ]

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

    function cavg_overfit_score( z_features ) result( score )
        real, intent(in) :: z_features(:)
        real             :: score
        if( size(z_features) /= CAVG_QUALITY_NFEATS ) THROW_HARD('cavg_overfit_score: invalid feature count')
        score = sum(z_features * CAVG_OVERFIT_SCORE_REJECT_WEIGHTS)
    end function cavg_overfit_score

    subroutine extract_cavg_quality_features( imgs, cls_oris, mskdiam, raw, hard_reject, quality_target )
        class(image),         intent(inout) :: imgs(:)
        type(oris),           intent(in)    :: cls_oris
        real,                 intent(in)    :: mskdiam
        real,    allocatable, intent(inout) :: raw(:,:)
        logical, allocatable, intent(inout) :: hard_reject(:)
        character(len=*), optional, intent(in) :: quality_target
        integer, allocatable :: pop(:), disc_area(:)
        real,    allocatable :: res(:), corr(:), corr_in(:)
        integer              :: ncls, i, ldim(3), pop_hard_threshold
        real                 :: smpd, rad_px
        type(image_bin), allocatable :: bin_img(:), cc_img(:), disc_img(:)
        real, allocatable    :: rmat_cc(:,:,:,:), rmat_disc(:,:,:,:)
        real                 :: centroid_norm, locvar_fg, locvar_bg
        real                 :: cc_area_frac
        real                 :: center_edge_snr, presence_score
        real                 :: int_boundary_tex_ratio, int_tex_coverage, tex_effective_area
        integer              :: ithr
        logical              :: no_component, mask_hard_reject, bad_pixels, overfit_target
        ncls = size(imgs)
        if( ncls == 0 ) THROW_HARD('extract_cavg_quality_features: no class averages')
        if( cls_oris%get_noris() /= ncls ) THROW_HARD('extract_cavg_quality_features: # cls oris /= # cavgs')
        if( mskdiam <= 0.0 ) THROW_HARD('extract_cavg_quality_features: mskdiam must be positive')
        if( allocated(raw)         ) deallocate(raw)
        if( allocated(hard_reject) ) deallocate(hard_reject)
        allocate(raw(ncls, CAVG_QUALITY_NFEATS), source=0.0)
        allocate(hard_reject(ncls),              source=.false.)
        pop  = cls_oris%get_all_asint('pop')
        res  = cls_oris%get_all('res')
        if( size(pop) /= ncls ) THROW_HARD('extract_cavg_quality_features: invalid pop size')
        if( size(res) /= ncls ) THROW_HARD('extract_cavg_quality_features: invalid res size')
        pop_hard_threshold = ceiling(real(sum(pop)) * POP_FRACTION_HARD_REJECT)
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
        overfit_target = .false.
        if( present(quality_target) ) overfit_target = trim(quality_target) == CAVG_QUALITY_TARGET_OVERFIT
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
        !$omp& mask_hard_reject,bad_pixels,locvar_fg,locvar_bg,center_edge_snr,presence_score,&
        !$omp& int_boundary_tex_ratio,int_tex_coverage,tex_effective_area,&
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
                int_boundary_tex_ratio = 0.0
                int_tex_coverage    = 0.0
                tex_effective_area  = 0.0
            else
                call measure_cavg_foreground_geometry(imgs(i), bin_img(ithr), cc_img(ithr), rmat_cc(:,:,:,ithr), &
                                                      rmat_disc(:,:,:,ithr), disc_area(ithr), rad_px, &
                                                      centroid_norm, cc_area_frac, no_component, mask_hard_reject)
                call measure_cavg_image_metrics(imgs(i), rad_px, locvar_fg, locvar_bg, center_edge_snr, &
                                                presence_score, int_boundary_tex_ratio, int_tex_coverage, &
                                                tex_effective_area)
            endif
            raw(i, I_LOG_POP)                 = log(real(max(pop(i), 0)) + 1.0)
            raw(i, I_NEG_LOG_RES)             = resolution_feature(res(i))
            raw(i, I_CENTERED)                = -centroid_norm
            raw(i, I_LOCVAR_FG)               = log(max(locvar_fg, LOG_EPS))
            raw(i, I_LOCVAR_BG)               = log(max(locvar_bg, LOG_EPS))
            raw(i, I_CORR_FRC)                = corr(i)
            raw(i, I_CENTER_EDGE_SNR)         = log(max(center_edge_snr, LOG_EPS))
            raw(i, I_CC_AREA_FRAC)            = cc_area_frac
            raw(i, I_PRESENCE)                = presence_score
            raw(i, I_LOG_INT_BOUNDARY_TEX)    = log(max(int_boundary_tex_ratio, LOG_EPS))
            raw(i, I_INT_TEX_COVERAGE)        = int_tex_coverage
            raw(i, I_TEX_EFFECTIVE_AREA)      = tex_effective_area
            raw(i, I_NEG_LOCVAR_FG)          = -log(max(locvar_fg, LOG_EPS))
            raw(i, I_NEG_LOCVAR_BG)          = -log(max(locvar_bg, LOG_EPS))
            raw(i, I_LOG_LOCVAR_FG_BG_RATIO) = log(max(locvar_fg, LOG_EPS)) - log(max(locvar_bg, LOG_EPS))
            if( overfit_target )then
                ! The overfit target keeps generic-quality failures trainable.
                ! Only true validity failures remain hard rejects.
                hard_reject(i) = pop(i) <= 0 .or. bad_pixels
            else
                ! Catastrophic population, resolution, and foreground-geometry
                ! failures are hard validity rejects. The population fraction and
                ! connected-component pruning mirror the microchunk rejector, while
                ! ordinary variation remains active scalar evidence for the model.
                hard_reject(i) = pop(i) <= 0 .or. pop(i) < pop_hard_threshold .or. &
                                                  res(i) > CAVG_RES_HARD_REJECT_A .or. &
                                                  bad_pixels .or. no_component .or. mask_hard_reject .or. &
                                                  (locvar_fg <= EPS .and. locvar_bg <= EPS)
            endif
        end do
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call bin_img(ithr)%kill_bimg()
            call cc_img(ithr)%kill_bimg()
            call disc_img(ithr)%kill_bimg()
        end do
        deallocate(bin_img, cc_img, disc_img, rmat_cc, rmat_disc, disc_area, pop, res, corr)
    end subroutine extract_cavg_quality_features

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
                                           int_boundary_tex_ratio, int_tex_coverage, tex_effective_area )
        class(image), intent(inout) :: img_src
        real,         intent(in)    :: rad_px
        real,         intent(out)   :: locvar_fg, locvar_bg, center_edge_snr
        real,         intent(out)   :: presence_score
        real,         intent(out)   :: int_boundary_tex_ratio, int_tex_coverage, tex_effective_area
        type(image)       :: img, img_bin, img_tex
        real, allocatable :: bin_mask(:,:,:), tex_rmat(:,:,:)
        integer           :: ldim(3)
        img  = img_src
        call img%zero_edgeavg()
        presence_score = img%presence()
        center_edge_snr = img%center_edge_snr(rad_px)
        call img%bp(0.0, SIGNAL_METRIC_LP)
        img_bin = img
        call otsu_img(img_bin)
        ldim = img%get_ldim()
        allocate(bin_mask(ldim(1), ldim(2), ldim(3)))
        call img_bin%get_rmat_sub(bin_mask)
        call img%loc_var_masked(bin_mask(:,:,1), LOCVAR_WINDOW, locvar_fg, locvar_bg)
        img_tex = img_src
        call img_tex%zero_edgeavg()
        call img_tex%bp(TEXTURE_METRIC_HP, TEXTURE_METRIC_LP)
        allocate(tex_rmat(ldim(1), ldim(2), ldim(3)))
        call img_tex%get_rmat_sub(tex_rmat)
        call measure_texture_metrics(tex_rmat(:,:,1), bin_mask(:,:,1), int_boundary_tex_ratio, &
                                     int_tex_coverage, tex_effective_area)
        deallocate(tex_rmat)
        deallocate(bin_mask)
        call img%kill()
        call img_bin%kill()
        call img_tex%kill()
    end subroutine measure_cavg_image_metrics

    subroutine measure_texture_metrics( tex, bin_mask, int_boundary_tex_ratio, int_tex_coverage, tex_effective_area )
        real, intent(in)  :: tex(:,:), bin_mask(:,:)
        real, intent(out) :: int_boundary_tex_ratio, int_tex_coverage, tex_effective_area
        logical, allocatable :: fg_mask(:,:), interior_mask(:,:), boundary_mask(:,:), bg_mask(:,:)
        real,    allocatable :: abs_tex(:,:), bg_vals(:)
        real :: bg_med, bg_dev, threshold, interior_mean, boundary_mean
        real :: sum_pos, sum_pos2, val
        integer :: nr, nc, i, j, n_fg, n_int, n_bg
        nr = size(tex, dim=1)
        nc = size(tex, dim=2)
        if( size(bin_mask, dim=1) /= nr .or. size(bin_mask, dim=2) /= nc ) &
            THROW_HARD('measure_texture_metrics: mask size mismatch')
        allocate(abs_tex(nr,nc), source=abs(tex))
        allocate(fg_mask(nr,nc), interior_mask(nr,nc), boundary_mask(nr,nc), bg_mask(nr,nc))
        fg_mask = bin_mask /= 0.0
        bg_mask = .not. fg_mask
        n_fg = count(fg_mask)
        if( n_fg <= 0 )then
            int_boundary_tex_ratio = 0.0
            int_tex_coverage       = 0.0
            tex_effective_area     = 0.0
            deallocate(abs_tex, fg_mask, interior_mask, boundary_mask, bg_mask)
            return
        endif
        call erode_logical_mask(fg_mask, TEXTURE_RING_WIDTH, interior_mask)
        if( count(interior_mask) < max(4, n_fg / 20) ) &
            call erode_logical_mask(fg_mask, 1, interior_mask)
        boundary_mask = fg_mask .and. .not. interior_mask
        n_int = count(interior_mask)
        n_bg  = count(bg_mask)
        if( n_bg > 1 )then
            bg_vals = pack(abs_tex, bg_mask)
            bg_med  = median(bg_vals)
            bg_dev  = mad_gau(bg_vals, bg_med)
            threshold = bg_med + TEXTURE_BG_MAD_SIGMA * max(bg_dev, LOG_EPS)
            deallocate(bg_vals)
        else
            bg_vals = reshape(abs_tex, (/nr * nc/))
            threshold = median(bg_vals)
            deallocate(bg_vals)
        endif
        threshold = max(threshold, LOG_EPS)
        interior_mean = masked_mean(abs_tex, interior_mask)
        boundary_mean = masked_mean(abs_tex, boundary_mask)
        int_boundary_tex_ratio = (interior_mean + LOG_EPS) / (boundary_mean + LOG_EPS)
        if( n_int > 0 )then
            int_tex_coverage = real(count((abs_tex > threshold) .and. interior_mask)) / real(n_int)
        else
            int_tex_coverage = 0.0
        endif
        sum_pos  = 0.0
        sum_pos2 = 0.0
        do j = 1, nc
            do i = 1, nr
                if( .not. fg_mask(i,j) ) cycle
                val = max(abs_tex(i,j) - threshold, 0.0)
                sum_pos  = sum_pos  + val
                sum_pos2 = sum_pos2 + val * val
            end do
        end do
        if( sum_pos2 > LOG_EPS )then
            tex_effective_area = min(1.0, max(0.0, (sum_pos * sum_pos) / (real(n_fg) * sum_pos2)))
        else
            tex_effective_area = 0.0
        endif
        deallocate(abs_tex, fg_mask, interior_mask, boundary_mask, bg_mask)
    end subroutine measure_texture_metrics

    subroutine erode_logical_mask( mask, radius, eroded )
        logical, intent(in)  :: mask(:,:)
        integer, intent(in)  :: radius
        logical, intent(out) :: eroded(:,:)
        integer :: nr, nc, i, j, ii, jj, r
        nr = size(mask, dim=1)
        nc = size(mask, dim=2)
        if( size(eroded, dim=1) /= nr .or. size(eroded, dim=2) /= nc ) &
            THROW_HARD('erode_logical_mask: output size mismatch')
        eroded = .false.
        r = max(0, radius)
        do j = 1, nc
            do i = 1, nr
                if( .not. mask(i,j) ) cycle
                eroded(i,j) = .true.
                do jj = max(1, j-r), min(nc, j+r)
                    do ii = max(1, i-r), min(nr, i+r)
                        if( mask(ii,jj) ) cycle
                        eroded(i,j) = .false.
                        exit
                    end do
                    if( .not. eroded(i,j) ) exit
                end do
                if( i-r < 1 .or. i+r > nr .or. j-r < 1 .or. j+r > nc ) eroded(i,j) = .false.
            end do
        end do
    end subroutine erode_logical_mask

    real function masked_mean( vals, mask )
        real,    intent(in) :: vals(:,:)
        logical, intent(in) :: mask(:,:)
        integer :: n
        if( size(vals, dim=1) /= size(mask, dim=1) .or. size(vals, dim=2) /= size(mask, dim=2) ) &
            THROW_HARD('masked_mean: size mismatch')
        n = count(mask)
        if( n > 0 )then
            masked_mean = sum(vals, mask=mask) / real(n)
        else
            masked_mean = 0.0
        endif
    end function masked_mean

    real function resolution_feature( res )
        real, intent(in) :: res
        if( res > EPS ) then
            resolution_feature = -log(res)
        else
            resolution_feature = -log(1000.0)
        end if
    end function resolution_feature

end module simple_cavg_quality_feats
