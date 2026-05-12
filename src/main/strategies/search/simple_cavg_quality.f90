!@descr: feature-space quality analysis for 2D class averages
module simple_cavg_quality
use simple_defs,             only: logfhandle
use simple_error,            only: simple_exception
use simple_image,            only: image
use simple_image_bin,        only: image_bin
use simple_oris,             only: oris
use simple_clustering_utils, only: cluster_dmat
use simple_math_ft,          only: calc_fourier_index
use simple_segmentation,     only: otsu_img
use simple_stat,             only: median, mad_gau
implicit none
private

public :: CAVG_QUALITY_NFEATS
public :: cavg_quality_result
public :: cavg_quality_feature_name
public :: extract_cavg_quality_features
public :: normalize_cavg_quality_features
public :: cluster_cavg_quality
public :: evaluate_cavg_quality

#include "simple_local_flags.inc"

integer, parameter :: CAVG_QUALITY_NFEATS = 8

real, parameter :: EPS                  = 1.0e-6
real, parameter :: LOG_EPS              = 1.0e-12
real, parameter :: CLIP_Z               = 4.0
real, parameter :: MIN_SCORE_SEPARATION = 0.15
real, parameter :: HP_SPEC              = 20.0
real, parameter :: LP_SPEC              = 6.0

integer, parameter :: I_LOG_POP         = 1
integer, parameter :: I_NEG_LOG_RES     = 2
integer, parameter :: I_MASK_INSIDE     = 3
integer, parameter :: I_CENTERED        = 4
integer, parameter :: I_LOCVAR_FG       = 5
integer, parameter :: I_LOCVAR_BG       = 6
integer, parameter :: I_SPEC_DYNRANGE   = 7
integer, parameter :: I_CC_SINGLE       = 8

character(len=32), parameter :: FEATURE_NAMES(CAVG_QUALITY_NFEATS) = [character(len=32) :: &
    'log_pop', 'neg_log_res', 'mask_inside', 'centered', &
    'log_locvar_fg', 'neg_log_locvar_bg', 'spectrum_dynrange', 'single_component']

real, parameter :: FEATURE_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    0.15, 0.15, 0.18, 0.10, 0.14, 0.08, 0.14, 0.06 ]

type :: cavg_quality_result
    real,    allocatable :: raw(:,:)
    real,    allocatable :: features(:,:)
    real,    allocatable :: scores(:)
    integer, allocatable :: states(:)
    integer, allocatable :: labels(:)
    integer, allocatable :: medoids(:)
    logical, allocatable :: hard_reject(:)
    real                 :: threshold  = 0.0
    real                 :: separation = 0.0
    integer              :: nclust     = 0
    integer              :: good_label = 0
    logical              :: used_threshold = .false.
end type cavg_quality_result

contains

    function cavg_quality_feature_name( i ) result( name )
        integer, intent(in) :: i
        character(len=32)   :: name
        if( i < 1 .or. i > CAVG_QUALITY_NFEATS ) THROW_HARD('invalid cavg quality feature index')
        name = FEATURE_NAMES(i)
    end function cavg_quality_feature_name

    subroutine evaluate_cavg_quality( imgs, cls_oris, mskdiam, quality )
        class(image),              intent(inout) :: imgs(:)
        type(oris),                intent(in)    :: cls_oris
        real,                      intent(in)    :: mskdiam
        type(cavg_quality_result), intent(inout) :: quality
        call reset_quality_result(quality)
        call extract_cavg_quality_features(imgs, cls_oris, mskdiam, quality%raw, quality%hard_reject)
        call normalize_cavg_quality_features(quality%raw, quality%hard_reject, quality%features)
        call cluster_cavg_quality(quality%features, quality%hard_reject, quality%states, quality%labels, &
                                                            quality%medoids, quality%scores, quality%threshold, quality%separation, &
                                                            quality%nclust, quality%good_label, quality%used_threshold)
    end subroutine evaluate_cavg_quality

    subroutine extract_cavg_quality_features( imgs, cls_oris, mskdiam, raw, hard_reject )
        class(image),         intent(inout) :: imgs(:)
        type(oris),           intent(in)    :: cls_oris
        real,                 intent(in)    :: mskdiam
        real,    allocatable, intent(inout) :: raw(:,:)
        logical, allocatable, intent(inout) :: hard_reject(:)
        integer, allocatable :: pop(:)
        real,    allocatable :: res(:)
        integer              :: ncls, i, ldim(3)
        real                 :: smpd, rad_px
        type(image_bin)      :: bin_img, cc_img, disc_img
        real, allocatable    :: rmat_cc(:,:,:), rmat_disc(:,:,:)
        real                 :: outside_frac, centroid_norm, locvar_fg, locvar_bg
        real                 :: spec_dynrange
        integer              :: nccs_valid
        logical              :: no_component
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
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( smpd <= 0.0 ) THROW_HARD('extract_cavg_quality_features: non-positive smpd')
        rad_px = (mskdiam / smpd) / 2.0
        call bin_img%new_bimg(ldim,  smpd, wthreads=.false.)
        call cc_img%new_bimg(ldim,   smpd, wthreads=.false.)
        call disc_img%disc(ldim,     smpd, rad_px)
        allocate(rmat_cc(ldim(1), ldim(2), ldim(3)), rmat_disc(ldim(1), ldim(2), ldim(3)))
        call disc_img%get_rmat_sub(rmat_disc)
        do i = 1, ncls
            call calc_mask_features(imgs(i), bin_img, cc_img, rmat_cc, rmat_disc, rad_px, &
                                                            outside_frac, centroid_norm, nccs_valid, no_component)
            call calc_signal_features(imgs(i), rad_px, locvar_fg, locvar_bg, spec_dynrange)
            raw(i, I_LOG_POP)       = log(real(max(pop(i), 0)) + 1.0)
            raw(i, I_NEG_LOG_RES)   = resolution_feature(res(i))
            raw(i, I_MASK_INSIDE)   = -outside_frac
            raw(i, I_CENTERED)      = -centroid_norm
            raw(i, I_LOCVAR_FG)     = log(max(locvar_fg, LOG_EPS))
            raw(i, I_LOCVAR_BG)     = -log(max(locvar_bg, LOG_EPS))
            raw(i, I_SPEC_DYNRANGE) = spec_dynrange
            raw(i, I_CC_SINGLE)     = -abs(real(nccs_valid - 1))
            hard_reject(i) = pop(i) <= 0 .or. no_component .or. &
                                              (locvar_fg <= EPS .and. locvar_bg <= EPS)
        end do
        call bin_img%kill_bimg()
        call cc_img%kill_bimg()
        call disc_img%kill_bimg()
        deallocate(rmat_cc, rmat_disc, pop, res)
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

    subroutine cluster_cavg_quality( features, hard_reject, states, labels, medoids, scores, threshold, &
                                                                      separation, nclust, good_label, used_threshold )
        real,                 intent(in)    :: features(:,:)
        logical,              intent(in)    :: hard_reject(:)
        integer, allocatable, intent(inout) :: states(:), labels(:), medoids(:)
        real,    allocatable, intent(inout) :: scores(:)
        real,                 intent(out)   :: threshold, separation
        integer,              intent(out)   :: nclust, good_label
        logical,              intent(out)   :: used_threshold
        real,    allocatable :: dmat(:,:), feats_fit(:,:), score_fit(:)
        integer, allocatable :: inds(:), labels_fit(:), medoids_fit(:)
        integer              :: ncls, nfit, i, j, k, good_fit_label, bad_fit_label
        real                 :: d, dmin, dmax, score1, score2
        ncls = size(features, dim=1)
        if( size(features, dim=2) /= CAVG_QUALITY_NFEATS ) THROW_HARD('cluster_cavg_quality: invalid feature count')
        if( size(hard_reject) /= ncls ) THROW_HARD('cluster_cavg_quality: invalid mask size')
        if( allocated(states) ) deallocate(states)
        if( allocated(labels) ) deallocate(labels)
        if( allocated(medoids)) deallocate(medoids)
        if( allocated(scores) ) deallocate(scores)
        allocate(states(ncls), labels(ncls), source=0)
        allocate(scores(ncls), source=0.0)
        scores = matmul(features, FEATURE_WEIGHTS)
        where( hard_reject ) scores = -CLIP_Z
        threshold      = 0.0
        separation     = 0.0
        nclust         = 0
        good_label     = 0
        used_threshold = .false.
        nfit = count(.not. hard_reject)
        if( nfit == 0 ) return
        inds = pack([(i, i=1,ncls)], .not. hard_reject)
        if( nfit < 4 ) then
            states(inds) = 1
            labels(inds) = 1
            allocate(medoids(1), source=inds(1))
            threshold  = minval(scores(inds)) - EPS
            nclust     = 1
            good_label = 1
            return
        end if
        allocate(feats_fit(nfit, CAVG_QUALITY_NFEATS), score_fit(nfit))
        do i = 1, nfit
            feats_fit(i,:) = features(inds(i),:)
            score_fit(i)   = scores(inds(i))
        end do
        allocate(dmat(nfit, nfit), source=0.0)
        do i = 1, nfit - 1
            do j = i + 1, nfit
                d = sqrt(sum((feats_fit(i,:) - feats_fit(j,:))**2))
                dmat(i,j) = d
                dmat(j,i) = d
            end do
        end do
        dmin = minval(dmat)
        dmax = maxval(dmat)
        if( dmax - dmin <= EPS ) then
            states(inds) = 1
            labels(inds) = 1
            allocate(medoids(1), source=inds(1))
            threshold  = minval(scores(inds)) - EPS
            nclust     = 1
            good_label = 1
            deallocate(dmat, feats_fit, score_fit, inds)
            return
        end if
        dmat = (dmat - dmin) / (dmax - dmin)
        nclust = 2
        call cluster_dmat(dmat, 'kmed', nclust, medoids_fit, labels_fit)
        if( nclust /= 2 ) THROW_HARD('cluster_cavg_quality: expected two k-medoids clusters')
        score1 = mean_score_for_label(score_fit, labels_fit, 1)
        score2 = mean_score_for_label(score_fit, labels_fit, 2)
        if( score1 >= score2 ) then
            good_fit_label = 1
            bad_fit_label  = 2
        else
            good_fit_label = 2
            bad_fit_label  = 1
        end if
        separation = abs(score1 - score2)
        threshold  = 0.5 * (mean_score_for_label(score_fit, labels_fit, good_fit_label) + &
                                                mean_score_for_label(score_fit, labels_fit, bad_fit_label))
        labels(inds) = labels_fit
        allocate(medoids(size(medoids_fit)))
        do k = 1, size(medoids_fit)
            medoids(k) = inds(medoids_fit(k))
        end do
        if( separation < MIN_SCORE_SEPARATION ) then
            states(inds)    = 1
            labels(inds)    = 1
            medoids         = [inds(1)]
            nclust          = 1
            good_label      = 1
            threshold       = minval(scores(inds)) - EPS
            used_threshold  = .false.
        else
            do i = 1, nfit
                if( scores(inds(i)) >= threshold ) states(inds(i)) = 1
            end do
            good_label     = good_fit_label
            used_threshold = .true.
        end if
        deallocate(dmat, feats_fit, score_fit, inds, labels_fit, medoids_fit)
    end subroutine cluster_cavg_quality

    subroutine calc_mask_features( img, bin_img, cc_img, rmat_cc, rmat_disc, rad_px, outside_frac, &
                                                                  centroid_norm, nccs_valid, no_component )
        class(image),     intent(inout) :: img
        type(image_bin),  intent(inout) :: bin_img, cc_img
        real,             intent(inout) :: rmat_cc(:,:,:)
        real,             intent(in)    :: rmat_disc(:,:,:)
        real,             intent(in)    :: rad_px
        real,             intent(out)   :: outside_frac, centroid_norm
        integer,          intent(out)   :: nccs_valid
        logical,          intent(out)   :: no_component
        real, allocatable :: ccsizes(:)
        integer           :: j, loc, nccs, area, outside
        integer           :: ldim(3)
        real              :: cc_diam, xy(2)
        ldim = img%get_ldim()
        outside_frac  = 1.0
        centroid_norm = 2.0
        nccs_valid    = 0
        no_component  = .false.
        call bin_img%copy(img)
        call bin_img%zero_edgeavg()
        call bin_img%bp(0.0, 30.0)
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
            return
        end if
        centroid_norm = 0.0
        call cc_img%order_ccs()
        call cc_img%update_img_rmat()
        call cc_img%get_nccs(nccs)
        do j = 1, nccs
            call cc_img%masscen_cc(j, xy)
            centroid_norm = max(centroid_norm, sqrt(xy(1)**2 + xy(2)**2) / max(rad_px, 1.0))
        end do
        ccsizes = cc_img%size_ccs()
        loc     = maxloc(ccsizes, dim=1)
        deallocate(ccsizes)
        call cc_img%cc2bin(loc)
        call cc_img%get_rmat_sub(rmat_cc)
        area    = count(rmat_cc > 0.0)
        outside = count(rmat_cc - rmat_disc > 0.0)
        outside_frac = real(outside) / real(max(area, 1))
    end subroutine calc_mask_features

    subroutine calc_signal_features( img_src, rad_px, locvar_fg, locvar_bg, spec_dynrange )
        class(image), intent(inout) :: img_src
        real,         intent(in)    :: rad_px
        real,         intent(out)   :: locvar_fg, locvar_bg, spec_dynrange
        type(image)       :: img, img_bin
        real, allocatable :: bin_mask(:,:,:), pspec(:)
        integer           :: ldim(3), k_hp, k_lp
        real              :: smpd
        img  = img_src
        call img%zero_edgeavg()
        call img%bp(0.0, 10.0)
        img_bin = img
        call otsu_img(img_bin)
        ldim = img%get_ldim()
        allocate(bin_mask(ldim(1), ldim(2), ldim(3)))
        call img_bin%get_rmat_sub(bin_mask)
        call img%loc_var_masked(bin_mask(:,:,1), 10, locvar_fg, locvar_bg)
        smpd = img%get_smpd()
        call img%mask2D_soft(rad_px, backgr=0.0)
        call img%spectrum('sqrt', pspec)
        k_hp = max(1, min(size(pspec), calc_fourier_index(HP_SPEC, ldim(1), smpd)))
        k_lp = max(1, min(size(pspec), calc_fourier_index(LP_SPEC, ldim(1), smpd)))
        spec_dynrange = pspec(k_hp) - pspec(k_lp)
        deallocate(bin_mask, pspec)
        call img%kill()
        call img_bin%kill()
    end subroutine calc_signal_features

    real function resolution_feature( res )
        real, intent(in) :: res
        if( res > EPS ) then
            resolution_feature = -log(res)
        else
            resolution_feature = -log(1000.0)
        end if
    end function resolution_feature

    real function mean_score_for_label( scores, labels, label )
        real,    intent(in) :: scores(:)
        integer, intent(in) :: labels(:), label
        integer             :: n
        n = count(labels == label)
        if( n == 0 ) THROW_HARD('mean_score_for_label: empty cluster')
        mean_score_for_label = sum(scores, mask=labels == label) / real(n)
    end function mean_score_for_label

    subroutine reset_quality_result( quality )
        type(cavg_quality_result), intent(inout) :: quality
        if( allocated(quality%raw)         ) deallocate(quality%raw)
        if( allocated(quality%features)    ) deallocate(quality%features)
        if( allocated(quality%scores)      ) deallocate(quality%scores)
        if( allocated(quality%states)      ) deallocate(quality%states)
        if( allocated(quality%labels)      ) deallocate(quality%labels)
        if( allocated(quality%medoids)     ) deallocate(quality%medoids)
        if( allocated(quality%hard_reject) ) deallocate(quality%hard_reject)
        quality%threshold      = 0.0
        quality%separation     = 0.0
        quality%nclust         = 0
        quality%good_label     = 0
        quality%used_threshold = .false.
    end subroutine reset_quality_result

end module simple_cavg_quality
