!@descr: feature inventory and feature extraction for 2D class-average quality
module simple_cavg_quality_feats
!$ use omp_lib
use simple_defs,         only: nthr_glob
use simple_error,        only: simple_exception
use simple_image,        only: image
use simple_image_bin,    only: image_bin
use simple_oris,         only: oris
use simple_segmentation, only: otsu_img
use simple_stat,         only: median, mad_gau
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, cavg_quality_feature_def
implicit none
private

public :: CAVG_QUALITY_NFEATS
public :: LOG_EPS
public :: I_LOG_POP
public :: I_NEG_LOG_RES
public :: I_MASK_INSIDE
public :: I_CENTERED
public :: I_LOCVAR_FG
public :: I_LOCVAR_BG
public :: I_CC_SINGLE
public :: I_CORR_FRC
public :: I_CENTER_EDGE_SNR
public :: I_CC_AREA_FRAC
public :: I_DETAIL_TEXTURE
public :: I_PRESENCE
public :: I_LOG_DETAIL_BG_SNR
public :: I_LOG_DETAIL_SIGNAL_RATIO
public :: I_DETAIL_COVERAGE
public :: I_DETAIL_EDGE_DENSITY
public :: cavg_quality_feature_def
public :: cavg_quality_feature_name
public :: cavg_quality_feature_description
public :: cavg_quality_feature_direction
public :: cavg_quality_feature_family
public :: cavg_quality_feature_model_eligible
public :: write_cavg_quality_feature_inventory
public :: extract_cavg_quality_features
public :: normalize_cavg_quality_features
public :: apply_cavg_quality_model_feature_exclusions
public :: apply_cavg_quality_model_mask_exclusions
public :: CAVG_RES_HARD_REJECT_A

#include "simple_local_flags.inc"

real,    parameter :: LOG_EPS             = 1.0e-12
real,    parameter :: FOREGROUND_SEG_LP   = 30.0
real,    parameter :: SIGNAL_METRIC_LP    = 10.0
real,    parameter :: DETAIL_BAND_HP_A    = 20.0
real,    parameter :: DETAIL_BAND_LP_A    =  6.0
real,    parameter :: CAVG_RES_HARD_REJECT_A = 40.0
integer, parameter :: LOCVAR_WINDOW       = 10
integer, parameter :: MASK_HARD_OUTSIDE_PIXELS = 10

integer, parameter :: I_LOG_POP           = 1
integer, parameter :: I_NEG_LOG_RES       = 2
integer, parameter :: I_MASK_INSIDE       = 3
integer, parameter :: I_CENTERED          = 4
integer, parameter :: I_LOCVAR_FG         = 5
integer, parameter :: I_LOCVAR_BG         = 6
integer, parameter :: I_CC_SINGLE         = 7
integer, parameter :: I_CORR_FRC          = 8
integer, parameter :: I_CENTER_EDGE_SNR   = 9
integer, parameter :: I_CC_AREA_FRAC      = 10
integer, parameter :: I_DETAIL_TEXTURE    = 11
integer, parameter :: I_PRESENCE          = 12
integer, parameter :: I_LOG_DETAIL_BG_SNR = 13
integer, parameter :: I_LOG_DETAIL_SIGNAL_RATIO = 14
integer, parameter :: I_DETAIL_COVERAGE   = 15
integer, parameter :: I_DETAIL_EDGE_DENSITY = 16

type(cavg_quality_feature_def), parameter :: FEATURE_DEFS(CAVG_QUALITY_NFEATS) = [ &
    cavg_quality_feature_def('log_pop', 'higher_is_better', &
        'log class population; larger classes have more particle support', 'population_support'), &
    cavg_quality_feature_def('neg_log_res', 'higher_is_better', &
        'negative log class resolution estimate; higher values indicate better nominal resolution', 'resolution'), &
    cavg_quality_feature_def('mask_inside', 'higher_is_better', &
        'negative outside-mask fraction for the main segmented component', 'foreground_geometry'), &
    cavg_quality_feature_def('centered', 'higher_is_better', &
        'negative normalized centroid displacement of segmented class-average signal', 'foreground_geometry'), &
    cavg_quality_feature_def('log_locvar_fg', 'higher_is_better', &
        'log local variance measured in the foreground Otsu mask', 'local_variance'), &
    cavg_quality_feature_def('log_locvar_bg', 'higher_is_better', &
        'log local variance measured in the complementary background mask', 'local_variance'), &
    cavg_quality_feature_def('single_component', 'higher_is_better', &
        'negative distance from exactly one valid connected component', 'foreground_geometry'), &
    cavg_quality_feature_def('corr_frc_proxy', 'higher_is_better', &
        'stored class correlation or FRC-like score when available in cls2D metadata', 'quality_proxy'), &
    cavg_quality_feature_def('log_center_edge_snr', 'higher_is_better', &
        'log central signal variance relative to edge variance in the class average', 'center_edge_signal'), &
    cavg_quality_feature_def('cc_area_frac', 'diagnostic', &
        'main connected-component area divided by the expected circular mask area', 'foreground_geometry'), &
    cavg_quality_feature_def('detail_texture', 'higher_is_better', &
        '20-6 A local texture heterogeneity and radial-residual detail', 'internal_detail'), &
    cavg_quality_feature_def('presence', 'higher_is_better', &
        'central signal excess relative to border background noise', 'presence'), &
    cavg_quality_feature_def('log_detail_bg_snr', 'diagnostic', &
        'log in-mask band-passed detail RMS relative to outside-mask detail RMS', 'internal_detail'), &
    cavg_quality_feature_def('log_detail_signal_ratio', 'diagnostic', &
        'log in-mask band-passed detail RMS relative to in-mask total signal RMS', 'internal_detail'), &
    cavg_quality_feature_def('detail_coverage', 'diagnostic', &
        'fraction of in-mask pixels with band-passed detail above local background and image-scale thresholds', &
        'internal_detail'), &
    cavg_quality_feature_def('detail_edge_density', 'diagnostic', &
        'fraction of in-mask pixels with band-passed gradient magnitude above background and image-scale thresholds', &
        'internal_detail') ]

contains

    subroutine apply_cavg_quality_model_feature_exclusions( weights )
        real, intent(inout) :: weights(:)
        integer :: ifeat
        if( size(weights) /= CAVG_QUALITY_NFEATS ) &
            THROW_HARD('apply_cavg_quality_model_feature_exclusions: invalid weight vector size')
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( .not. cavg_quality_feature_model_eligible(ifeat) ) weights(ifeat) = 0.0
        end do
    end subroutine apply_cavg_quality_model_feature_exclusions

    subroutine apply_cavg_quality_model_mask_exclusions( mask )
        logical, intent(inout) :: mask(:)
        integer :: ifeat
        if( size(mask) /= CAVG_QUALITY_NFEATS ) &
            THROW_HARD('apply_cavg_quality_model_mask_exclusions: invalid feature mask size')
        do ifeat = 1, CAVG_QUALITY_NFEATS
            if( .not. cavg_quality_feature_model_eligible(ifeat) ) mask(ifeat) = .false.
        end do
    end subroutine apply_cavg_quality_model_mask_exclusions

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

    logical function cavg_quality_feature_model_eligible( i )
        integer, intent(in) :: i
        if( i < 1 .or. i > CAVG_QUALITY_NFEATS ) THROW_HARD('invalid cavg quality feature index')
        ! Model eligibility is family-based so genuinely misleading feature
        ! families can be excluded in one place. The current inventory contains
        ! only hard-reject diagnostics and active scalar model candidates.
        cavg_quality_feature_model_eligible = .true.
    end function cavg_quality_feature_model_eligible

    subroutine write_cavg_quality_feature_inventory( funit )
        integer, intent(in) :: funit
        integer :: i
        write(funit,'(A)') '# feature_inventory_header=index,name,direction,family,description'
        do i = 1, CAVG_QUALITY_NFEATS
            write(funit,'(A,I0,A,A,A,A,A,A,A,A)') '# feature_inventory,', i, ',', trim(FEATURE_DEFS(i)%name), ',', &
                trim(FEATURE_DEFS(i)%direction), ',', trim(FEATURE_DEFS(i)%family), ',', trim(FEATURE_DEFS(i)%description)
        end do
    end subroutine write_cavg_quality_feature_inventory

    subroutine extract_cavg_quality_features( imgs, cls_oris, mskdiam, raw, hard_reject )
        class(image),         intent(inout) :: imgs(:)
        type(oris),           intent(in)    :: cls_oris
        real,                 intent(in)    :: mskdiam
        real,    allocatable, intent(inout) :: raw(:,:)
        logical, allocatable, intent(inout) :: hard_reject(:)
        integer, allocatable :: pop(:), disc_area(:)
        real,    allocatable :: res(:), corr(:), corr_in(:)
        integer              :: ncls, i, ldim(3)
        real                 :: smpd, rad_px
        type(image_bin), allocatable :: bin_img(:), cc_img(:), disc_img(:)
        real, allocatable    :: rmat_cc(:,:,:,:), rmat_disc(:,:,:,:)
        real                 :: outside_frac, centroid_norm, locvar_fg, locvar_bg
        real                 :: cc_area_frac
        real                 :: center_edge_snr, presence_score
        real                 :: detail_texture, detail_bg_snr, detail_signal_ratio, detail_coverage
        real                 :: detail_edge_density
        integer              :: nccs_valid, ithr
        logical              :: no_component, mask_hard_reject, bad_pixels
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
        !$omp parallel do default(shared) private(i,ithr,outside_frac,centroid_norm,nccs_valid,no_component,&
        !$omp& mask_hard_reject,bad_pixels,locvar_fg,locvar_bg,center_edge_snr,presence_score,&
        !$omp& cc_area_frac,detail_texture,detail_bg_snr,detail_signal_ratio,&
        !$omp& detail_coverage,detail_edge_density)&
        !$omp proc_bind(close) schedule(static)
        do i = 1, ncls
            ithr = omp_get_thread_num() + 1
            bad_pixels = imgs(i)%contains_nans()
            if( bad_pixels )then
                outside_frac     = 1.0
                centroid_norm    = 2.0
                cc_area_frac     = 0.0
                nccs_valid       = 0
                no_component     = .true.
                mask_hard_reject = .true.
                locvar_fg        = 0.0
                locvar_bg        = 0.0
                center_edge_snr  = 0.0
                presence_score   = 0.0
                detail_texture   = 0.0
                detail_bg_snr    = 0.0
                detail_signal_ratio = 0.0
                detail_coverage  = 0.0
                detail_edge_density = 0.0
            else
                call measure_cavg_foreground_geometry(imgs(i), bin_img(ithr), cc_img(ithr), rmat_cc(:,:,:,ithr), &
                                                      rmat_disc(:,:,:,ithr), disc_area(ithr), rad_px, outside_frac, &
                                                      centroid_norm, cc_area_frac, nccs_valid, &
                                                      no_component, mask_hard_reject)
                call measure_cavg_image_metrics(imgs(i), rad_px, locvar_fg, locvar_bg, center_edge_snr, &
                                                presence_score, detail_texture, detail_bg_snr, detail_signal_ratio, &
                                                detail_coverage, detail_edge_density)
            endif
            raw(i, I_LOG_POP)       = log(real(max(pop(i), 0)) + 1.0)
            raw(i, I_NEG_LOG_RES)   = resolution_feature(res(i))
            raw(i, I_MASK_INSIDE)   = -outside_frac
            raw(i, I_CENTERED)      = -centroid_norm
            raw(i, I_LOCVAR_FG)     = log(max(locvar_fg, LOG_EPS))
            raw(i, I_LOCVAR_BG)     = log(max(locvar_bg, LOG_EPS))
            raw(i, I_CC_SINGLE)     = -abs(real(nccs_valid - 1))
            raw(i, I_CORR_FRC)      = corr(i)
            raw(i, I_CENTER_EDGE_SNR) = log(max(center_edge_snr, LOG_EPS))
            raw(i, I_CC_AREA_FRAC)  = cc_area_frac
            raw(i, I_DETAIL_TEXTURE) = detail_texture
            raw(i, I_PRESENCE)      = presence_score
            raw(i, I_LOG_DETAIL_BG_SNR) = detail_bg_snr
            raw(i, I_LOG_DETAIL_SIGNAL_RATIO) = detail_signal_ratio
            raw(i, I_DETAIL_COVERAGE) = detail_coverage
            raw(i, I_DETAIL_EDGE_DENSITY) = detail_edge_density
            ! Catastrophic resolution and geometry failures are hard validity
            ! rejects; ordinary resolution variation remains an active scalar
            ! model feature through neg_log_res.
            hard_reject(i) = pop(i) <= 0 .or. res(i) > CAVG_RES_HARD_REJECT_A .or. bad_pixels .or. &
                                              no_component .or. mask_hard_reject .or. &
                                              (locvar_fg <= EPS .and. locvar_bg <= EPS)
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

    subroutine measure_cavg_foreground_geometry( img, bin_img, cc_img, rmat_cc, rmat_disc, disc_area, rad_px, outside_frac, &
                                                 centroid_norm, cc_area_frac, nccs_valid, no_component, &
                                                 mask_hard_reject )
        class(image),     intent(inout) :: img
        type(image_bin),  intent(inout) :: bin_img, cc_img
        real,             intent(inout) :: rmat_cc(:,:,:)
        real,             intent(in)    :: rmat_disc(:,:,:)
        integer,          intent(in)    :: disc_area
        real,             intent(in)    :: rad_px
        real,             intent(out)   :: outside_frac, centroid_norm, cc_area_frac
        integer,          intent(out)   :: nccs_valid
        logical,          intent(out)   :: no_component, mask_hard_reject
        real, allocatable :: ccsizes(:)
        integer           :: j, loc, nccs, area, outside
        real              :: xy(2)
        outside_frac  = 1.0
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
        outside_frac = real(outside) / real(max(area, 1))
        if( outside > MASK_HARD_OUTSIDE_PIXELS ) mask_hard_reject = .true.
    end subroutine measure_cavg_foreground_geometry

    subroutine measure_cavg_image_metrics( img_src, rad_px, locvar_fg, locvar_bg, center_edge_snr, presence_score, &
                                           detail_texture, detail_bg_snr, detail_signal_ratio, detail_coverage, &
                                           detail_edge_density )
        class(image), intent(inout) :: img_src
        real,         intent(in)    :: rad_px
        real,         intent(out)   :: locvar_fg, locvar_bg, center_edge_snr
        real,         intent(out)   :: presence_score
        real,         intent(out)   :: detail_texture
        real,         intent(out)   :: detail_bg_snr, detail_signal_ratio, detail_coverage, detail_edge_density
        type(image)       :: img, img_bin
        real, allocatable :: bin_mask(:,:,:)
        integer           :: ldim(3)
        img  = img_src
        call img%zero_edgeavg()
        presence_score = img%presence()
        center_edge_snr = img%center_edge_snr(rad_px)
        call measure_cavg_detail_metrics(img, rad_px, detail_texture, detail_bg_snr, detail_signal_ratio, &
                                         detail_coverage, detail_edge_density)
        call img%bp(0.0, SIGNAL_METRIC_LP)
        img_bin = img
        call otsu_img(img_bin)
        ldim = img%get_ldim()
        allocate(bin_mask(ldim(1), ldim(2), ldim(3)))
        call img_bin%get_rmat_sub(bin_mask)
        call img%loc_var_masked(bin_mask(:,:,1), LOCVAR_WINDOW, locvar_fg, locvar_bg)
        deallocate(bin_mask)
        call img%kill()
        call img_bin%kill()
    end subroutine measure_cavg_image_metrics

    subroutine measure_cavg_detail_metrics( img_src, rad_px, detail_texture, detail_bg_snr, detail_signal_ratio, &
                                            detail_coverage, detail_edge_density )
        class(image), intent(inout) :: img_src
        real,         intent(in)    :: rad_px
        real,         intent(out)   :: detail_texture, detail_bg_snr, detail_signal_ratio
        real,         intent(out)   :: detail_coverage, detail_edge_density
        type(image) :: detail_img
        real, allocatable :: orig(:,:,:), detail(:,:,:)
        real(8), allocatable :: ring_sum(:), ring_mean(:), block_energy(:,:)
        integer, allocatable :: ring_count(:), block_count(:,:)
        integer :: ldim(3), cx, cy, ix, iy, ir, nrings
        integer :: ib, jb, nxblocks, nyblocks, block_size, nblocks
        integer :: n_in, n_bg, n_grad_in, n_grad_bg, n_cov, n_edge
        real :: dx, dy, r2, rad2, core_rad2, val, grad, detail_rms_in, detail_rms_bg
        real :: signal_rms_in, grad_rms_in, grad_rms_bg, detail_thr, grad_thr
        real :: block_cv, residual_frac
        real(8) :: detail_s2_in, detail_s2_bg, signal_s2_in, grad_s2_in, grad_s2_bg
        real(8) :: residual_s2_in, resid
        real(8) :: block_mean, block_mean_all, block_var_all
        detail_texture      = 0.0
        detail_bg_snr       = 0.0
        detail_signal_ratio = 0.0
        detail_coverage     = 0.0
        detail_edge_density = 0.0
        ldim = img_src%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('measure_cavg_detail_metrics: 2D only')
        cx   = ldim(1) / 2 + 1
        cy   = ldim(2) / 2 + 1
        rad2 = rad_px**2
        if( rad2 <= EPS ) return
        detail_img = img_src
        ! The detail diagnostics are deliberately cheap scalar probes: one
        ! broad 20-6 A band-pass, then real-space in-mask/background ratios,
        ! coverage estimates, and a texture score. The texture score removes
        ! the radial mean first, then combines coarse local residual-energy
        ! heterogeneity in the core mask with total radial-residual signal.
        ! This avoids punishing symmetric top views merely for being radially
        ! organized, while still asking for non-flat local variance that
        ! smooth ice blobs should lack.
        call detail_img%bp(DETAIL_BAND_HP_A, DETAIL_BAND_LP_A)
        orig   = img_src%get_rmat()
        detail = detail_img%get_rmat()
        nrings    = max(1, int(rad_px) + 2)
        block_size = max(4, int(rad_px / 4.0))
        nxblocks  = (ldim(1) + block_size - 1) / block_size
        nyblocks  = (ldim(2) + block_size - 1) / block_size
        core_rad2 = (0.85 * rad_px)**2
        allocate(ring_sum(nrings), ring_mean(nrings), ring_count(nrings))
        allocate(block_energy(nxblocks, nyblocks), block_count(nxblocks, nyblocks))
        ring_sum   = 0.0d0
        ring_mean  = 0.0d0
        ring_count = 0
        block_energy = 0.0d0
        block_count  = 0
        detail_s2_in = 0.0d0
        detail_s2_bg = 0.0d0
        signal_s2_in = 0.0d0
        grad_s2_in   = 0.0d0
        grad_s2_bg   = 0.0d0
        residual_s2_in = 0.0d0
        n_in         = 0
        n_bg         = 0
        n_grad_in    = 0
        n_grad_bg    = 0
        do iy = 1, ldim(2)
            dy = real(iy - cy)
            do ix = 1, ldim(1)
                dx  = real(ix - cx)
                r2  = dx * dx + dy * dy
                val = detail(ix,iy,1)
                if( r2 <= rad2 )then
                    ir = min(nrings, int(sqrt(r2)) + 1)
                    ring_sum(ir)   = ring_sum(ir) + real(val,kind=8)
                    ring_count(ir) = ring_count(ir) + 1
                    detail_s2_in = detail_s2_in + real(val,kind=8)**2
                    signal_s2_in = signal_s2_in + real(orig(ix,iy,1),kind=8)**2
                    n_in = n_in + 1
                else
                    detail_s2_bg = detail_s2_bg + real(val,kind=8)**2
                    n_bg = n_bg + 1
                endif
                if( ix <= 1 .or. ix >= ldim(1) .or. iy <= 1 .or. iy >= ldim(2) ) cycle
                grad = sqrt(0.25 * (detail(ix+1,iy,1) - detail(ix-1,iy,1))**2 + &
                            0.25 * (detail(ix,iy+1,1) - detail(ix,iy-1,1))**2)
                if( r2 <= rad2 )then
                    grad_s2_in = grad_s2_in + real(grad,kind=8)**2
                    n_grad_in = n_grad_in + 1
                else
                    grad_s2_bg = grad_s2_bg + real(grad,kind=8)**2
                    n_grad_bg = n_grad_bg + 1
                endif
            end do
        end do
        if( n_in <= 0 )then
            call detail_img%kill()
            deallocate(orig, detail, ring_sum, ring_mean, ring_count, block_energy, block_count)
            return
        endif
        do ir = 1, nrings
            if( ring_count(ir) > 0 ) ring_mean(ir) = ring_sum(ir) / real(ring_count(ir),kind=8)
        end do
        do iy = 1, ldim(2)
            dy = real(iy - cy)
            do ix = 1, ldim(1)
                dx = real(ix - cx)
                r2 = dx * dx + dy * dy
                if( r2 > rad2 ) cycle
                ir = min(nrings, int(sqrt(r2)) + 1)
                resid = real(detail(ix,iy,1),kind=8) - ring_mean(ir)
                residual_s2_in = residual_s2_in + resid**2
                if( r2 <= core_rad2 )then
                    ib = min(nxblocks, (ix - 1) / block_size + 1)
                    jb = min(nyblocks, (iy - 1) / block_size + 1)
                    block_energy(ib,jb) = block_energy(ib,jb) + resid**2
                    block_count(ib,jb)  = block_count(ib,jb) + 1
                endif
            end do
        end do
        nblocks = 0
        block_mean_all = 0.0d0
        do jb = 1, nyblocks
            do ib = 1, nxblocks
                if( block_count(ib,jb) <= 0 ) cycle
                block_energy(ib,jb) = block_energy(ib,jb) / real(block_count(ib,jb),kind=8)
                block_mean_all = block_mean_all + block_energy(ib,jb)
                nblocks = nblocks + 1
            end do
        end do
        block_cv = 0.0
        if( nblocks > 1 )then
            block_mean_all = block_mean_all / real(nblocks,kind=8)
            block_var_all = 0.0d0
            do jb = 1, nyblocks
                do ib = 1, nxblocks
                    if( block_count(ib,jb) <= 0 ) cycle
                    block_mean = block_energy(ib,jb)
                    block_var_all = block_var_all + (block_mean - block_mean_all)**2
                end do
            end do
            block_cv = sqrt(real(block_var_all / real(nblocks - 1,kind=8))) / &
                       max(real(block_mean_all), LOG_EPS)
        endif
        if( detail_s2_in > 0.0d0 )then
            residual_frac = sqrt(real(residual_s2_in / detail_s2_in))
        else
            residual_frac = 0.0
        endif
        detail_texture = log(max(block_cv + residual_frac, LOG_EPS))
        detail_rms_in = sqrt(real(detail_s2_in / real(max(n_in, 1),kind=8)))
        signal_rms_in = sqrt(real(signal_s2_in / real(max(n_in, 1),kind=8)))
        if( n_bg > 0 )then
            detail_rms_bg = sqrt(real(detail_s2_bg / real(n_bg,kind=8)))
        else
            detail_rms_bg = detail_rms_in
        endif
        if( n_grad_in > 0 )then
            grad_rms_in = sqrt(real(grad_s2_in / real(n_grad_in,kind=8)))
        else
            grad_rms_in = 0.0
        endif
        if( n_grad_bg > 0 )then
            grad_rms_bg = sqrt(real(grad_s2_bg / real(n_grad_bg,kind=8)))
        else
            grad_rms_bg = grad_rms_in
        endif
        detail_bg_snr       = log(max(detail_rms_in / max(detail_rms_bg, LOG_EPS), LOG_EPS))
        detail_signal_ratio = log(max(detail_rms_in / max(signal_rms_in, LOG_EPS), LOG_EPS))
        detail_thr = max(2.0 * detail_rms_bg, 0.5 * detail_rms_in)
        grad_thr   = max(2.0 * grad_rms_bg,   0.5 * grad_rms_in)
        n_cov  = 0
        n_edge = 0
        do iy = 2, ldim(2) - 1
            dy = real(iy - cy)
            do ix = 2, ldim(1) - 1
                dx = real(ix - cx)
                if( dx * dx + dy * dy > rad2 ) cycle
                if( abs(detail(ix,iy,1)) > detail_thr ) n_cov = n_cov + 1
                grad = sqrt(0.25 * (detail(ix+1,iy,1) - detail(ix-1,iy,1))**2 + &
                            0.25 * (detail(ix,iy+1,1) - detail(ix,iy-1,1))**2)
                if( grad > grad_thr ) n_edge = n_edge + 1
            end do
        end do
        detail_coverage     = real(n_cov)  / real(max(n_grad_in, 1))
        detail_edge_density = real(n_edge) / real(max(n_grad_in, 1))
        call detail_img%kill()
        deallocate(orig, detail, ring_sum, ring_mean, ring_count, block_energy, block_count)
    end subroutine measure_cavg_detail_metrics

    real function resolution_feature( res )
        real, intent(in) :: res
        if( res > EPS ) then
            resolution_feature = -log(res)
        else
            resolution_feature = -log(1000.0)
        end if
    end function resolution_feature

end module simple_cavg_quality_feats
