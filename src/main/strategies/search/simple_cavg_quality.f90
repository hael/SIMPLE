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
public :: write_cavg_quality_reference_analysis

#include "simple_local_flags.inc"

integer, parameter :: CAVG_QUALITY_NFEATS = 8

real, parameter :: EPS                     = 1.0e-6
real, parameter :: LOG_EPS                 = 1.0e-12
real, parameter :: CLIP_Z                  = 4.0
real, parameter :: MIN_SCORE_SEPARATION    = 0.15
real, parameter :: BOUNDARY_MARGIN_DEFAULT = -0.30
real, parameter :: HP_SPEC                 = 20.0
real, parameter :: LP_SPEC                 = 6.0

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
    'log_locvar_fg', 'log_locvar_bg', 'spectrum_dynrange', 'single_component']

! Tuning guide:
! - FEATURE_WEIGHTS is the main calibration surface. It maps robust-normalized
!   features into the scalar score used to identify the high-quality cluster.
! - BOUNDARY_MARGIN_DEFAULT is an internal threshold offset. Positive values
!   retain more borderline classes; negative values raise the boundary and
!   reject more junk. This is deliberately not exposed as a command-line knob.
! - MIN_SCORE_SEPARATION controls when two clusters are trusted at all; if the
!   cluster mean scores are closer than this, the selector falls back to keeping
!   all non-hard-rejected classes.
! - CLIP_Z limits the influence of feature outliers after median/MAD
!   normalization. The hard-reject rule in extract_cavg_quality_features is the
!   only pre-clustering veto and should remain conservative.
! - HP_SPEC/LP_SPEC only affect the spectrum_dynrange diagnostic. That feature
!   is retained for analysis but not scored by default because high spectral
!   dynamic range can also indicate ring/ice artifacts.
real, parameter :: FEATURE_WEIGHTS(CAVG_QUALITY_NFEATS) = [ &
    0.27, 0.13, 0.00, 0.15, 0.30, 0.15, 0.00, 0.00 ]

type :: cavg_quality_result
    real,    allocatable :: raw(:,:)
    real,    allocatable :: features(:,:)
    real,    allocatable :: scores(:)
    integer, allocatable :: states(:)
    integer, allocatable :: labels(:)
    integer, allocatable :: medoids(:)
    logical, allocatable :: hard_reject(:)
    real                 :: threshold        = 0.0
    real                 :: raw_threshold    = 0.0
    real                 :: threshold_margin = 0.0
    real                 :: separation       = 0.0
    integer              :: nclust           = 0
    integer              :: good_label       = 0
    logical              :: used_threshold = .false.
end type cavg_quality_result

contains

    function cavg_quality_feature_name( i ) result( name )
        integer, intent(in) :: i
        character(len=32)   :: name
        if( i < 1 .or. i > CAVG_QUALITY_NFEATS ) THROW_HARD('invalid cavg quality feature index')
        name = FEATURE_NAMES(i)
    end function cavg_quality_feature_name

    subroutine evaluate_cavg_quality( imgs, cls_oris, mskdiam, quality, boundary_margin )
        class(image),              intent(inout) :: imgs(:)
        type(oris),                intent(in)    :: cls_oris
        real,                      intent(in)    :: mskdiam
        type(cavg_quality_result), intent(inout) :: quality
        real, optional,            intent(in)    :: boundary_margin
        real :: margin
        call reset_quality_result(quality)
        margin = BOUNDARY_MARGIN_DEFAULT
        if( present(boundary_margin) ) margin = boundary_margin
        call extract_cavg_quality_features(imgs, cls_oris, mskdiam, quality%raw, quality%hard_reject)
        call normalize_cavg_quality_features(quality%raw, quality%hard_reject, quality%features)
        call cluster_cavg_quality(quality%features, quality%hard_reject, quality%states, quality%labels, &
                                  quality%medoids, quality%scores, quality%threshold, quality%raw_threshold, &
                                  quality%threshold_margin, quality%separation, quality%nclust, quality%good_label, &
                                  quality%used_threshold, margin)
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
            raw(i, I_LOCVAR_BG)     = log(max(locvar_bg, LOG_EPS))
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

    subroutine cluster_cavg_quality( features, hard_reject, states, labels, medoids, scores, threshold, raw_threshold, &
                                                                      threshold_margin, separation, nclust, good_label, &
                                                                      used_threshold, boundary_margin )
        real,                 intent(in)    :: features(:,:)
        logical,              intent(in)    :: hard_reject(:)
        integer, allocatable, intent(inout) :: states(:), labels(:), medoids(:)
        real,    allocatable, intent(inout) :: scores(:)
        real,                 intent(out)   :: threshold, raw_threshold, threshold_margin, separation
        integer,              intent(out)   :: nclust, good_label
        logical,              intent(out)   :: used_threshold
        real, optional,       intent(in)    :: boundary_margin
        real,    allocatable :: dmat(:,:), feats_fit(:,:), score_fit(:)
        integer, allocatable :: inds(:), labels_fit(:), medoids_fit(:)
        integer              :: ncls, nfit, i, j, k, good_fit_label, bad_fit_label
        real                 :: d, dmin, dmax, score1, score2, margin
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
        threshold        = 0.0
        raw_threshold    = 0.0
        threshold_margin = 0.0
        separation       = 0.0
        nclust           = 0
        good_label       = 0
        used_threshold   = .false.
        margin = BOUNDARY_MARGIN_DEFAULT
        if( present(boundary_margin) ) margin = boundary_margin
        nfit = count(.not. hard_reject)
        if( nfit == 0 ) return
        inds = pack([(i, i=1,ncls)], .not. hard_reject)
        if( nfit < 4 ) then
            states(inds) = 1
            labels(inds) = 1
            allocate(medoids(1), source=inds(1))
            threshold        = minval(scores(inds)) - EPS
            raw_threshold    = threshold
            threshold_margin = 0.0
            nclust           = 1
            good_label       = 1
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
            threshold        = minval(scores(inds)) - EPS
            raw_threshold    = threshold
            threshold_margin = 0.0
            nclust           = 1
            good_label       = 1
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
        ! The weighted feature score is fixed by design; this margin is the single
        ! exposed decision-boundary knob. Positive margins lower the effective
        ! threshold to retain more classes; negative margins raise it to reject
        ! more junk.
        raw_threshold = 0.5 * (mean_score_for_label(score_fit, labels_fit, good_fit_label) + &
                                                mean_score_for_label(score_fit, labels_fit, bad_fit_label))
        threshold = raw_threshold - margin
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
            threshold_margin = 0.0
            used_threshold  = .false.
        else
            threshold_margin = margin
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

    subroutine write_cavg_quality_reference_analysis( quality, reference_states, prefix )
        type(cavg_quality_result), intent(in) :: quality
        integer,                   intent(in) :: reference_states(:)
        character(len=*),          intent(in) :: prefix
        logical, allocatable :: auto(:), ref(:)
        real,    allocatable :: suggested_weights(:)
        integer :: ncls, ifeat, tp, fp, tn, fn, ngood, nbad, funit
        real    :: precision, recall, specificity, f1, balacc, accuracy
        real    :: best_bal_thr, best_balacc, best_f1_thr, best_f1
        real    :: auc, med_good, med_bad, mad_good, mad_bad, sep
        ncls = size(reference_states)
        if( size(quality%states) /= ncls ) THROW_HARD('write_cavg_quality_reference_analysis: state size mismatch')
        if( size(quality%scores) /= ncls ) THROW_HARD('write_cavg_quality_reference_analysis: score size mismatch')
        if( size(quality%features, dim=1) /= ncls ) THROW_HARD('write_cavg_quality_reference_analysis: feature size mismatch')
        ref  = reference_states > 0
        auto = quality%states > 0
        ngood = count(ref)
        nbad  = ncls - ngood
        call calc_confusion(auto, ref, tp, fp, tn, fn)
        call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
        call write_threshold_scan(trim(prefix)//'_threshold_scan.txt', quality%scores, reference_states, &
                                  best_bal_thr, best_balacc, best_f1_thr, best_f1)
        allocate(suggested_weights(CAVG_QUALITY_NFEATS), source=0.0)
        do ifeat = 1, CAVG_QUALITY_NFEATS
            suggested_weights(ifeat) = max(0.0, auc_for_values(quality%features(:,ifeat), reference_states) - 0.5)
        end do
        if( sum(suggested_weights) > EPS ) then
            suggested_weights = suggested_weights / sum(suggested_weights)
        else
            suggested_weights = FEATURE_WEIGHTS
        end if
        open(newunit=funit, file=trim(prefix)//'_summary.txt', status='replace', action='write')
        write(funit,'(A)') '# cluster_cavgs_quality reference analysis'
        write(funit,'(A,I0)') 'n_classes=', ncls
        write(funit,'(A,I0)') 'manual_selected=', ngood
        write(funit,'(A,I0)') 'manual_rejected=', nbad
        write(funit,'(A,I0)') 'auto_selected=', count(auto)
        write(funit,'(A,I0)') 'true_positive=', tp
        write(funit,'(A,I0)') 'false_positive=', fp
        write(funit,'(A,I0)') 'true_negative=', tn
        write(funit,'(A,I0)') 'false_negative=', fn
        write(funit,'(A,F10.5)') 'precision=', precision
        write(funit,'(A,F10.5)') 'recall=', recall
        write(funit,'(A,F10.5)') 'specificity=', specificity
        write(funit,'(A,F10.5)') 'f1=', f1
        write(funit,'(A,F10.5)') 'balanced_accuracy=', balacc
        write(funit,'(A,F10.5)') 'accuracy=', accuracy
        write(funit,'(A,F10.5)') 'score_auc=', auc_for_values(quality%scores, reference_states)
        write(funit,'(A,F10.5)') 'raw_score_threshold=', quality%raw_threshold
        write(funit,'(A,F10.5)') 'threshold_boundary_margin=', quality%threshold_margin
        write(funit,'(A,F10.5)') 'current_score_threshold=', quality%threshold
        write(funit,'(A,F10.5)') 'best_balacc_threshold=', best_bal_thr
        write(funit,'(A,F10.5)') 'best_balacc=', best_balacc
        write(funit,'(A,F10.5)') 'best_f1_threshold=', best_f1_thr
        write(funit,'(A,F10.5)') 'best_f1=', best_f1
        write(funit,'(A)') ''
        write(funit,'(A)') 'feature,auc,median_manual_good,median_manual_bad,robust_separation,current_weight,suggested_weight'
        do ifeat = 1, CAVG_QUALITY_NFEATS
            auc      = auc_for_values(quality%features(:,ifeat), reference_states)
            med_good = median_by_state(quality%features(:,ifeat), ref)
            med_bad  = median_by_state(quality%features(:,ifeat), .not. ref)
            mad_good = mad_by_state(quality%features(:,ifeat), ref, med_good)
            mad_bad  = mad_by_state(quality%features(:,ifeat), .not. ref, med_bad)
            sep      = safe_div(med_good - med_bad, 0.5 * (mad_good + mad_bad) + EPS)
            write(funit,'(A,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5)') &
                trim(cavg_quality_feature_name(ifeat)), ',', auc, ',', med_good, ',', med_bad, ',', sep, ',', &
                FEATURE_WEIGHTS(ifeat), ',', suggested_weights(ifeat)
        end do
        close(funit)
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(prefix)//'_summary.txt'
        write(logfhandle,'(A,A)') '>>> WROTE ', trim(prefix)//'_threshold_scan.txt'
        write(logfhandle,'(A,F8.3,A,F8.3)') '>>> BEST QUALITY THRESHOLD BALACC/F1: ', best_bal_thr, ' / ', best_f1_thr
        deallocate(auto, ref, suggested_weights)
    end subroutine write_cavg_quality_reference_analysis

    subroutine write_threshold_scan( fname, scores, reference_states, best_bal_thr, best_balacc, best_f1_thr, best_f1 )
        character(len=*), intent(in)  :: fname
        real,             intent(in)  :: scores(:)
        integer,          intent(in)  :: reference_states(:)
        real,             intent(out) :: best_bal_thr, best_balacc, best_f1_thr, best_f1
        logical, allocatable :: ref(:), pred(:)
        integer :: funit, i, tp, fp, tn, fn
        real    :: threshold, precision, recall, specificity, f1, balacc, accuracy
        ref = reference_states > 0
        allocate(pred(size(scores)), source=.false.)
        best_bal_thr = 0.0
        best_f1_thr  = 0.0
        best_balacc  = -huge(1.0)
        best_f1      = -huge(1.0)
        open(newunit=funit, file=fname, status='replace', action='write')
        write(funit,'(A)') 'score_threshold,selected,tp,fp,tn,fn,precision,recall,specificity,f1,balanced_accuracy,accuracy'
        call write_one_threshold(minval(scores) - EPS)
        do i = 1, size(scores)
            call write_one_threshold(scores(i))
        end do
        call write_one_threshold(maxval(scores) + EPS)
        close(funit)
        deallocate(ref, pred)

    contains

        subroutine write_one_threshold( thresh )
            real, intent(in) :: thresh
            threshold = thresh
            pred = scores >= threshold
            call calc_confusion(pred, ref, tp, fp, tn, fn)
            call calc_binary_metrics(tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy)
            if( balacc > best_balacc ) then
                best_balacc  = balacc
                best_bal_thr = threshold
            end if
            if( f1 > best_f1 ) then
                best_f1     = f1
                best_f1_thr = threshold
            end if
            write(funit,'(F12.6,A,I0,A,I0,A,I0,A,I0,A,I0,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5,A,F10.5)') &
                threshold, ',', count(pred), ',', tp, ',', fp, ',', tn, ',', fn, ',', precision, ',', recall, ',', &
                specificity, ',', f1, ',', balacc, ',', accuracy
        end subroutine write_one_threshold

    end subroutine write_threshold_scan

    subroutine calc_confusion( pred, ref, tp, fp, tn, fn )
        logical, intent(in)  :: pred(:), ref(:)
        integer, intent(out) :: tp, fp, tn, fn
        if( size(pred) /= size(ref) ) THROW_HARD('calc_confusion: size mismatch')
        tp = count( pred .and.  ref)
        fp = count( pred .and. .not. ref)
        tn = count(.not. pred .and. .not. ref)
        fn = count(.not. pred .and.  ref)
    end subroutine calc_confusion

    subroutine calc_binary_metrics( tp, fp, tn, fn, precision, recall, specificity, f1, balacc, accuracy )
        integer, intent(in)  :: tp, fp, tn, fn
        real,    intent(out) :: precision, recall, specificity, f1, balacc, accuracy
        precision   = safe_div(real(tp), real(tp + fp))
        recall      = safe_div(real(tp), real(tp + fn))
        specificity = safe_div(real(tn), real(tn + fp))
        f1          = safe_div(2.0 * precision * recall, precision + recall)
        balacc      = 0.5 * (recall + specificity)
        accuracy    = safe_div(real(tp + tn), real(tp + fp + tn + fn))
    end subroutine calc_binary_metrics

    real function auc_for_values( vals, reference_states )
        real,    intent(in) :: vals(:)
        integer, intent(in) :: reference_states(:)
        integer :: i, j, ngood, nbad
        real    :: wins
        if( size(vals) /= size(reference_states) ) THROW_HARD('auc_for_values: size mismatch')
        ngood = count(reference_states > 0)
        nbad  = size(reference_states) - ngood
        if( ngood == 0 .or. nbad == 0 ) then
            auc_for_values = 0.5
            return
        end if
        wins = 0.0
        do i = 1, size(vals)
            if( reference_states(i) <= 0 ) cycle
            do j = 1, size(vals)
                if( reference_states(j) > 0 ) cycle
                if( vals(i) > vals(j) ) then
                    wins = wins + 1.0
                else if( abs(vals(i) - vals(j)) <= EPS ) then
                    wins = wins + 0.5
                end if
            end do
        end do
        auc_for_values = wins / real(ngood * nbad)
    end function auc_for_values

    real function median_by_state( vals, mask )
        real,    intent(in) :: vals(:)
        logical, intent(in) :: mask(:)
        real, allocatable :: packed(:)
        if( size(vals) /= size(mask) ) THROW_HARD('median_by_state: size mismatch')
        if( count(mask) == 0 ) then
            median_by_state = 0.0
            return
        end if
        packed = pack(vals, mask)
        median_by_state = median(packed)
        deallocate(packed)
    end function median_by_state

    real function mad_by_state( vals, mask, med )
        real,    intent(in) :: vals(:), med
        logical, intent(in) :: mask(:)
        real, allocatable :: packed(:)
        if( size(vals) /= size(mask) ) THROW_HARD('mad_by_state: size mismatch')
        if( count(mask) < 2 ) then
            mad_by_state = 0.0
            return
        end if
        packed = pack(vals, mask)
        mad_by_state = mad_gau(packed, med)
        deallocate(packed)
    end function mad_by_state

    real function safe_div( num, den )
        real, intent(in) :: num, den
        if( abs(den) <= EPS ) then
            safe_div = 0.0
        else
            safe_div = num / den
        end if
    end function safe_div

    subroutine reset_quality_result( quality )
        type(cavg_quality_result), intent(inout) :: quality
        if( allocated(quality%raw)         ) deallocate(quality%raw)
        if( allocated(quality%features)    ) deallocate(quality%features)
        if( allocated(quality%scores)      ) deallocate(quality%scores)
        if( allocated(quality%states)      ) deallocate(quality%states)
        if( allocated(quality%labels)      ) deallocate(quality%labels)
        if( allocated(quality%medoids)     ) deallocate(quality%medoids)
        if( allocated(quality%hard_reject) ) deallocate(quality%hard_reject)
        quality%threshold        = 0.0
        quality%raw_threshold    = 0.0
        quality%threshold_margin = 0.0
        quality%separation       = 0.0
        quality%nclust           = 0
        quality%good_label       = 0
        quality%used_threshold   = .false.
    end subroutine reset_quality_result

end module simple_cavg_quality
