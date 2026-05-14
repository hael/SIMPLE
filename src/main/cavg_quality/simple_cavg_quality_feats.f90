!@descr: feature inventory and feature extraction for 2D class-average quality
module simple_cavg_quality_feats
!$ use omp_lib
use simple_defs,         only: ICE_BAND1, nthr_glob
use simple_error,        only: simple_exception
use simple_histogram,    only: histogram
use simple_image,        only: image
use simple_image_bin,    only: image_bin
use simple_jiffys,       only: swap
use simple_oris,         only: oris
use simple_math_ft,      only: calc_fourier_index
use simple_segmentation, only: otsu_img
use simple_stat,         only: median, mad_gau
use simple_cavg_quality_types, only: CAVG_QUALITY_NFEATS, EPS, CLIP_Z, cavg_quality_feature_def
use simple_cavg_quality_stats, only: normalize_quality_dmat
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
public :: I_SPEC_DYNRANGE
public :: I_CC_SINGLE
public :: I_CORR_FRC
public :: I_CENTER_EDGE_SNR
public :: I_NEG_ICE_SCORE
public :: I_HIST_ENTROPY
public :: I_CC_AREA_FRAC
public :: I_CC_DIAMETER_NORM
public :: I_PRESENCE
public :: I_LOG_CONTRAST
public :: I_LOG_HIST_VARIANCE
public :: I_LOG_FG_BG_LOCVAR_RATIO
public :: cavg_quality_feature_def
public :: cavg_quality_feature_name
public :: cavg_quality_feature_description
public :: cavg_quality_feature_direction
public :: write_cavg_quality_feature_inventory
public :: extract_cavg_quality_features
public :: normalize_cavg_quality_features
public :: calc_histogram_quality_signal
public :: calc_power_spectrum_quality_signal

#include "simple_local_flags.inc"

real,    parameter :: LOG_EPS             = 1.0e-12
real,    parameter :: HP_SPEC             = 20.0
real,    parameter :: LP_SPEC             = 6.0
real,    parameter :: ICE_BG_START        = 15.0
real,    parameter :: ICE_BG_END          = 6.0
real,    parameter :: FOREGROUND_SEG_LP   = 30.0
real,    parameter :: SIGNAL_METRIC_LP    = 10.0
integer, parameter :: LOCVAR_WINDOW       = 10
integer, parameter :: NHISTBINS           = 128
integer, parameter :: MASK_HARD_OUTSIDE_PIXELS = 10

integer, parameter :: I_LOG_POP           = 1
integer, parameter :: I_NEG_LOG_RES       = 2
integer, parameter :: I_MASK_INSIDE       = 3
integer, parameter :: I_CENTERED          = 4
integer, parameter :: I_LOCVAR_FG         = 5
integer, parameter :: I_LOCVAR_BG         = 6
integer, parameter :: I_SPEC_DYNRANGE     = 7
integer, parameter :: I_CC_SINGLE         = 8
integer, parameter :: I_CORR_FRC          = 9
integer, parameter :: I_CENTER_EDGE_SNR   = 10
integer, parameter :: I_NEG_ICE_SCORE     = 11
integer, parameter :: I_HIST_ENTROPY      = 12
integer, parameter :: I_CC_AREA_FRAC      = 13
integer, parameter :: I_CC_DIAMETER_NORM  = 14
integer, parameter :: I_PRESENCE          = 15
integer, parameter :: I_LOG_CONTRAST      = 16
integer, parameter :: I_LOG_HIST_VARIANCE = 17
integer, parameter :: I_LOG_FG_BG_LOCVAR_RATIO = 18

type(cavg_quality_feature_def), parameter :: FEATURE_DEFS(CAVG_QUALITY_NFEATS) = [ &
    cavg_quality_feature_def('log_pop', 'higher_is_better', &
        'log class population; larger classes have more particle support'), &
    cavg_quality_feature_def('neg_log_res', 'higher_is_better', &
        'negative log class resolution estimate; higher values indicate better nominal resolution'), &
    cavg_quality_feature_def('mask_inside', 'higher_is_better', &
        'negative outside-mask fraction for the main segmented component'), &
    cavg_quality_feature_def('centered', 'higher_is_better', &
        'negative normalized centroid displacement of segmented class-average signal'), &
    cavg_quality_feature_def('log_locvar_fg', 'higher_is_better', &
        'log local variance measured in the foreground Otsu mask'), &
    cavg_quality_feature_def('log_locvar_bg', 'higher_is_better', &
        'log local variance measured in the complementary background mask'), &
    cavg_quality_feature_def('spectrum_dynrange', 'diagnostic', &
        'difference between high- and low-frequency spectral samples; retained as a diagnostic'), &
    cavg_quality_feature_def('single_component', 'higher_is_better', &
        'negative distance from exactly one valid connected component'), &
    cavg_quality_feature_def('corr_frc_proxy', 'higher_is_better', &
        'stored class correlation or FRC-like score when available in cls2D metadata'), &
    cavg_quality_feature_def('log_center_edge_snr', 'higher_is_better', &
        'log central signal variance relative to edge variance in the class average'), &
    cavg_quality_feature_def('neg_ice_score', 'higher_is_better', &
        'negative class-average ice-ring excess score around the 3.7 Angstrom band'), &
    cavg_quality_feature_def('hist_entropy', 'diagnostic', &
        'bounded masked intensity-histogram entropy from the same histograms used for Hellinger distances'), &
    cavg_quality_feature_def('cc_area_frac', 'diagnostic', &
        'main connected-component area divided by the expected circular mask area'), &
    cavg_quality_feature_def('cc_diameter_norm', 'diagnostic', &
        'main connected-component diameter divided by the expected mask diameter'), &
    cavg_quality_feature_def('presence', 'higher_is_better', &
        'central signal excess relative to border background noise'), &
    cavg_quality_feature_def('log_contrast', 'diagnostic', &
        'log full-image intensity standard deviation after edge background subtraction'), &
    cavg_quality_feature_def('log_hist_variance', 'diagnostic', &
        'log masked intensity-histogram variance from the same histograms used for Hellinger distances'), &
    cavg_quality_feature_def('log_fg_bg_locvar_ratio', 'higher_is_better', &
        'log ratio of foreground to background local variance') ]

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

    subroutine write_cavg_quality_feature_inventory( funit )
        integer, intent(in) :: funit
        integer :: i
        write(funit,'(A)') '# feature_inventory_header=index,name,direction,description'
        do i = 1, CAVG_QUALITY_NFEATS
            write(funit,'(A,I0,A,A,A,A,A,A)') '# feature_inventory,', i, ',', trim(FEATURE_DEFS(i)%name), ',', &
                trim(FEATURE_DEFS(i)%direction), ',', trim(FEATURE_DEFS(i)%description)
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
        real                 :: cc_area_frac, cc_diameter_norm
        real                 :: spec_dynrange, center_edge_snr, ice_score, presence_score, contrast
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
        !$omp& mask_hard_reject,bad_pixels,locvar_fg,locvar_bg,spec_dynrange,center_edge_snr,ice_score,presence_score,&
        !$omp& contrast,cc_area_frac,cc_diameter_norm)&
        !$omp proc_bind(close) schedule(static)
        do i = 1, ncls
            ithr = omp_get_thread_num() + 1
            bad_pixels = imgs(i)%contains_nans()
            if( bad_pixels )then
                outside_frac     = 1.0
                centroid_norm    = 2.0
                cc_area_frac     = 0.0
                cc_diameter_norm = 0.0
                nccs_valid       = 0
                no_component     = .true.
                mask_hard_reject = .true.
                locvar_fg        = 0.0
                locvar_bg        = 0.0
                spec_dynrange    = 0.0
                center_edge_snr  = 0.0
                ice_score        = 0.0
                presence_score   = 0.0
                contrast         = 0.0
            else
                call measure_cavg_foreground_geometry(imgs(i), bin_img(ithr), cc_img(ithr), rmat_cc(:,:,:,ithr), &
                                                      rmat_disc(:,:,:,ithr), disc_area(ithr), rad_px, outside_frac, &
                                                      centroid_norm, cc_area_frac, cc_diameter_norm, nccs_valid, &
                                                      no_component, mask_hard_reject)
                call measure_cavg_image_metrics(imgs(i), rad_px, locvar_fg, locvar_bg, spec_dynrange, &
                                                center_edge_snr, ice_score, presence_score, contrast)
            endif
            raw(i, I_LOG_POP)       = log(real(max(pop(i), 0)) + 1.0)
            raw(i, I_NEG_LOG_RES)   = resolution_feature(res(i))
            raw(i, I_MASK_INSIDE)   = -outside_frac
            raw(i, I_CENTERED)      = -centroid_norm
            raw(i, I_LOCVAR_FG)     = log(max(locvar_fg, LOG_EPS))
            raw(i, I_LOCVAR_BG)     = log(max(locvar_bg, LOG_EPS))
            raw(i, I_SPEC_DYNRANGE) = spec_dynrange
            raw(i, I_CC_SINGLE)     = -abs(real(nccs_valid - 1))
            raw(i, I_CORR_FRC)      = corr(i)
            raw(i, I_CENTER_EDGE_SNR) = log(max(center_edge_snr, LOG_EPS))
            raw(i, I_NEG_ICE_SCORE) = -ice_score
            raw(i, I_CC_AREA_FRAC)  = cc_area_frac
            raw(i, I_CC_DIAMETER_NORM) = cc_diameter_norm
            raw(i, I_PRESENCE)      = presence_score
            raw(i, I_LOG_CONTRAST)  = log(max(contrast, LOG_EPS))
            raw(i, I_LOG_FG_BG_LOCVAR_RATIO) = log(max(locvar_fg, LOG_EPS)) - log(max(locvar_bg, LOG_EPS))
            ! Geometry failures are hard validity rejects; the model should not
            ! rescue classes whose segmented foreground lies outside the mask.
            hard_reject(i) = pop(i) <= 0 .or. bad_pixels .or. no_component .or. mask_hard_reject .or. &
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

    subroutine calc_histogram_quality_signal( imgs, mskdiam, hard_reject, hist_dmat, hist_entropy, hist_variance )
        class(image),         intent(inout) :: imgs(:)
        real,                 intent(in)    :: mskdiam
        logical,              intent(in)    :: hard_reject(:)
        real,    allocatable, intent(inout) :: hist_dmat(:,:)
        real,    allocatable, intent(inout), optional :: hist_entropy(:)
        real,    allocatable, intent(inout), optional :: hist_variance(:)
        type(histogram), allocatable :: hists(:)
        real, allocatable :: dmat_hd(:,:)
        real :: smpd, rad_px, mm(2), mm_i(2)
        integer :: ncls, nfit, i, j
        logical :: ok
        ncls = size(imgs)
        if( size(hard_reject) /= ncls ) THROW_HARD('calc_histogram_quality_signal: invalid mask size')
        if( allocated(hist_dmat) ) deallocate(hist_dmat)
        allocate(hist_dmat(ncls, ncls), source=0.0)
        if( present(hist_entropy) )then
            if( allocated(hist_entropy) ) deallocate(hist_entropy)
            allocate(hist_entropy(ncls), source=0.0)
        endif
        if( present(hist_variance) )then
            if( allocated(hist_variance) ) deallocate(hist_variance)
            allocate(hist_variance(ncls), source=0.0)
        endif
        nfit = count(.not. hard_reject)
        if( nfit < 4 ) return
        smpd = imgs(1)%get_smpd()
        if( smpd <= 0.0 ) THROW_HARD('calc_histogram_quality_signal: non-positive smpd')
        rad_px = (mskdiam / smpd) / 2.0
        mm = [huge(1.0), -huge(1.0)]
        do i = 1, ncls
            if( hard_reject(i) ) cycle
            mm_i  = imgs(i)%minmax(radius=rad_px)
            mm(1) = min(mm(1), mm_i(1))
            mm(2) = max(mm(2), mm_i(2))
        end do
        if( mm(2) - mm(1) <= EPS ) return
        allocate(hists(ncls))
        do i = 1, ncls
            if( hard_reject(i) ) cycle
            call hists(i)%new(imgs(i), NHISTBINS, minmax=mm, radius=rad_px)
            if( present(hist_entropy) ) hist_entropy(i) = hists(i)%entropy(bounded=.true.)
            if( present(hist_variance) ) hist_variance(i) = log(max(hists(i)%variance(), LOG_EPS))
        end do
        allocate(dmat_hd(ncls,ncls), source=0.0)
        !$omp parallel do default(shared) private(i,j) schedule(dynamic) proc_bind(close)
        do i = 1, ncls - 1
            if( hard_reject(i) ) cycle
            do j = i + 1, ncls
                if( hard_reject(j) ) cycle
                dmat_hd(i,j)  = hists(i)%HD(hists(j))
                dmat_hd(j,i)  = dmat_hd(i,j)
            end do
        end do
        !$omp end parallel do
        ! Hellinger is the histogram signal used by cluster_cavgs_quality.
        ! Earlier averaging with TVD/JSD could dilute a useful, stable
        ! bounded metric with noisier or implementation-sensitive distances.
        call normalize_quality_dmat(dmat_hd, ok)
        if( ok ) hist_dmat = dmat_hd
        call hists(:)%kill
        deallocate(hists, dmat_hd)
    end subroutine calc_histogram_quality_signal

    subroutine calc_power_spectrum_quality_signal( imgs, mskdiam, hard_reject, spec_dmat )
        class(image),         intent(inout) :: imgs(:)
        real,                 intent(in)    :: mskdiam
        logical,              intent(in)    :: hard_reject(:)
        real,    allocatable, intent(inout) :: spec_dmat(:,:)
        type(image), allocatable :: cavg_threads(:)
        real, allocatable :: profiles(:,:), pspec(:)
        integer, allocatable :: inds(:)
        integer :: ncls, nfit, ldim(3), lfny, kfrom, kto, nspec, i, j, ithr, tmpi
        real    :: smpd, rad_px, dist
        logical :: ok
        ncls = size(imgs)
        if( size(hard_reject) /= ncls ) THROW_HARD('calc_power_spectrum_quality_signal: invalid mask size')
        if( allocated(spec_dmat) ) deallocate(spec_dmat)
        allocate(spec_dmat(ncls, ncls), source=0.0)
        nfit = count(.not. hard_reject)
        if( nfit < 4 ) return
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( smpd <= 0.0 ) THROW_HARD('calc_power_spectrum_quality_signal: non-positive smpd')
        lfny = imgs(1)%get_lfny(1)
        kfrom = max(1, min(lfny, calc_fourier_index(HP_SPEC, ldim(1), smpd)))
        kto   = max(1, min(lfny, calc_fourier_index(LP_SPEC, ldim(1), smpd)))
        if( kto < kfrom )then
            tmpi  = kfrom
            kfrom = kto
            kto   = tmpi
        endif
        nspec = kto - kfrom + 1
        if( nspec < 2 ) return
        rad_px = (mskdiam / smpd) / 2.0
        allocate(inds(nfit))
        allocate(profiles(nfit, nspec), source=0.0)
        inds = pack([(i, i=1,ncls)], .not. hard_reject)
        allocate(cavg_threads(nthr_glob))
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%new(ldim, smpd)
        end do
        call cavg_threads(1)%memoize_mask_coords
        !$omp parallel do default(shared) private(i,ithr,pspec) schedule(static) proc_bind(close)
        do i = 1, nfit
            ithr = omp_get_thread_num() + 1
            call cavg_threads(ithr)%copy(imgs(inds(i)))
            call cavg_threads(ithr)%norm
            call cavg_threads(ithr)%mask2D_soft(rad_px)
            call cavg_threads(ithr)%spectrum('sqrt', pspec)
            profiles(i,:) = pspec(kfrom:kto)
            deallocate(pspec)
        end do
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%kill
        end do
        deallocate(cavg_threads)
        !$omp parallel do default(shared) private(i,j,dist) schedule(dynamic) proc_bind(close)
        do i = 1, nfit - 1
            do j = i + 1, nfit
                dist = sqrt(sum((profiles(i,:) - profiles(j,:))**2))
                spec_dmat(inds(i), inds(j)) = dist
                spec_dmat(inds(j), inds(i)) = dist
            end do
        end do
        !$omp end parallel do
        call normalize_quality_dmat(spec_dmat, ok)
        if( .not. ok ) spec_dmat = 0.0
        deallocate(inds, profiles)
    end subroutine calc_power_spectrum_quality_signal

    subroutine measure_cavg_foreground_geometry( img, bin_img, cc_img, rmat_cc, rmat_disc, disc_area, rad_px, outside_frac, &
                                                 centroid_norm, cc_area_frac, cc_diameter_norm, nccs_valid, no_component, &
                                                 mask_hard_reject )
        class(image),     intent(inout) :: img
        type(image_bin),  intent(inout) :: bin_img, cc_img
        real,             intent(inout) :: rmat_cc(:,:,:)
        real,             intent(in)    :: rmat_disc(:,:,:)
        integer,          intent(in)    :: disc_area
        real,             intent(in)    :: rad_px
        real,             intent(out)   :: outside_frac, centroid_norm, cc_area_frac, cc_diameter_norm
        integer,          intent(out)   :: nccs_valid
        logical,          intent(out)   :: no_component, mask_hard_reject
        real, allocatable :: ccsizes(:)
        integer           :: j, loc, nccs, area, outside
        integer           :: ldim(3)
        real              :: cc_diam, xy(2)
        ldim = img%get_ldim()
        outside_frac  = 1.0
        centroid_norm = 2.0
        cc_area_frac  = 0.0
        cc_diameter_norm = 0.0
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
            if( cc_diam > ldim(1) ) then
                call cc_img%elim_cc(j, update=.false.)
                nccs_valid = nccs_valid - 1
            end if
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
        call cc_img%diameter_cc(loc, cc_diam)
        cc_diameter_norm = cc_diam / max(2.0 * rad_px, 1.0)
        call cc_img%cc2bin(loc)
        call cc_img%get_rmat_sub(rmat_cc)
        area    = count(rmat_cc > 0.0)
        outside = count(rmat_cc - rmat_disc > 0.0)
        cc_area_frac = real(area) / real(max(disc_area, 1))
        outside_frac = real(outside) / real(max(area, 1))
        if( outside > MASK_HARD_OUTSIDE_PIXELS ) mask_hard_reject = .true.
    end subroutine measure_cavg_foreground_geometry

    subroutine measure_cavg_image_metrics( img_src, rad_px, locvar_fg, locvar_bg, spec_dynrange, center_edge_snr, &
                                           ice_score, presence_score, contrast )
        class(image), intent(inout) :: img_src
        real,         intent(in)    :: rad_px
        real,         intent(out)   :: locvar_fg, locvar_bg, spec_dynrange, center_edge_snr, ice_score
        real,         intent(out)   :: presence_score, contrast
        type(image)       :: img, img_bin, img_ice
        real, allocatable :: bin_mask(:,:,:), pspec(:), ice_pspec(:)
        integer           :: ldim(3), k_hp, k_lp
        real              :: smpd
        ice_score = 0.0
        if( img_src%get_smpd() <= (ICE_BAND1 / 2.0) )then
            img_ice = img_src
            call img_ice%zero_edgeavg()
            call img_ice%mask2D_soft(rad_px, backgr=0.0)
            call img_ice%spectrum('power', ice_pspec)
            ice_score = class_average_ice_score(ice_pspec, img_ice%get_box(), img_ice%get_smpd())
            call img_ice%kill()
        endif
        img  = img_src
        call img%zero_edgeavg()
        presence_score = img%presence()
        contrast       = img%contrast()
        center_edge_snr = img%center_edge_snr(rad_px)
        call img%bp(0.0, SIGNAL_METRIC_LP)
        img_bin = img
        call otsu_img(img_bin)
        ldim = img%get_ldim()
        allocate(bin_mask(ldim(1), ldim(2), ldim(3)))
        call img_bin%get_rmat_sub(bin_mask)
        call img%loc_var_masked(bin_mask(:,:,1), LOCVAR_WINDOW, locvar_fg, locvar_bg)
        smpd = img%get_smpd()
        call img%mask2D_soft(rad_px, backgr=0.0)
        call img%spectrum('sqrt', pspec)
        k_hp = max(1, min(size(pspec), calc_fourier_index(HP_SPEC, ldim(1), smpd)))
        k_lp = max(1, min(size(pspec), calc_fourier_index(LP_SPEC, ldim(1), smpd)))
        spec_dynrange = pspec(k_hp) - pspec(k_lp)
        deallocate(bin_mask, pspec)
        if( allocated(ice_pspec) ) deallocate(ice_pspec)
        call img%kill()
        call img_bin%kill()
    end subroutine measure_cavg_image_metrics

    real function class_average_ice_score( pspec, box, smpd )
        real,    intent(in) :: pspec(:)
        integer, intent(in) :: box
        real,    intent(in) :: smpd
        integer :: k_ice, k_ice1, k_ice2, k_bg1, k_bg2, k, nbg
        real    :: ice_avg, bg_avg
        class_average_ice_score = 0.0
        if( smpd <= 0.0 ) return
        if( smpd > (ICE_BAND1 / 2.0) ) return
        if( size(pspec) < 5 ) return
        k_ice = max(1, min(size(pspec), calc_fourier_index(ICE_BAND1, box, smpd)))
        if( k_ice <= 1 .or. k_ice >= size(pspec) ) return
        k_ice1 = max(1, k_ice - 1)
        k_ice2 = min(size(pspec), k_ice + 1)
        ice_avg = sum(pspec(k_ice1:k_ice2)) / real(k_ice2 - k_ice1 + 1)
        k_bg1 = max(1, min(size(pspec), calc_fourier_index(ICE_BG_START, box, smpd)))
        k_bg2 = max(1, min(size(pspec), calc_fourier_index(ICE_BG_END,   box, smpd)))
        if( k_bg1 > k_bg2 ) call swap(k_bg1, k_bg2)
        bg_avg = 0.0
        nbg    = 0
        do k = k_bg1, k_bg2
            if( abs(k - k_ice) <= 3 ) cycle
            bg_avg = bg_avg + pspec(k)
            nbg    = nbg + 1
        end do
        if( nbg < 1 ) return
        bg_avg = bg_avg / real(nbg)
        if( bg_avg <= EPS ) return
        class_average_ice_score = max(0.0, ice_avg / bg_avg - 1.0)
    end function class_average_ice_score

    real function resolution_feature( res )
        real, intent(in) :: res
        if( res > EPS ) then
            resolution_feature = -log(res)
        else
            resolution_feature = -log(1000.0)
        end if
    end function resolution_feature

end module simple_cavg_quality_feats
