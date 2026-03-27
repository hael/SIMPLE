!@descr: utilities for microchunk-based 2D clustering in stream
!==============================================================================
! MODULE: simple_stream_microchunk_utils
!
! PURPOSE:
!   Provides scoring and rejection routines for 2D class averages produced
!   during stream microchunk processing. Three rejection strategies are
!   implemented, intended to be applied in sequence:
!
!     1. reject_outliers — removes statistically anomalous classes based on
!                          variance, maximum pixel value, centre/edge SNR, and
!                          image skewness using robust z-scores.
!     2. reject_auto     — ranks surviving classes by a composite quality score
!                          and rejects those below a percentile threshold within
!                          the top half of the distribution.
!     3. reject_basic    — removes classes with low population or poor
!                          resolution directly from a cls2D oris object.
!
! SCORING:
!   calc_rejection_score — computes a weighted sum of SNR, signal presence,
!                          and contrast for a single class-average image.
!
! REJECTION THRESHOLDS (module parameters):
!   VARIANCE_MAX_Z  — z-score upper bound on class variance             (3.0)
!   VARIANCE_MIN_Z  — z-score lower bound on class variance            (-1.0)
!   MAXPIXEL_Z      — z-score upper bound on maximum pixel value        (3.0)
!   CEN_EDGE_SNR_Z  — z-score lower bound on centre/edge SNR          (-2.0)
!   SKEW_Z          — z-score lower bound on image skewness            (-2.0)
!   AUTO_SNR_W      — weight for SNR in composite score                (0.5)
!   AUTO_PRESENCE_W — weight for signal presence in composite score   (0.02)
!   AUTO_CONTRAST_W — weight for contrast in composite score           (5.0)
!   AUTO_THRESHOLD  — percentile cut within top-50% for auto reject   (10.0)
!   BASIC_MIN_POP   — minimum class population for reject_basic        (10.0)
!   BASIC_MAX_RES   — maximum class resolution (Å) for reject_basic   (30.0)
!
! DEPENDENCIES:
!   simple_stream_api, simple_oris, simple_image, simple_stat,
!   simple_srch_sort_loc
!==============================================================================
module simple_stream_microchunk_utils
  use simple_stream_api
  use simple_oris,          only: oris
  use simple_image,         only: image
  use simple_stat,          only: robust_z_scores
  use simple_srch_sort_loc, only: hpsort, reverse

  implicit none
  private
  public :: reject_outliers, reject_auto, reject_basic, calc_rejection_score
#include "simple_local_flags.inc"

  ! Outlier rejection z-score thresholds
  real, parameter :: VARIANCE_MAX_Z  =  3.0
  real, parameter :: VARIANCE_MIN_Z  = -1.0
  real, parameter :: MAXPIXEL_Z      =  3.0
  real, parameter :: CEN_EDGE_SNR_Z  = -2.0
  real, parameter :: SKEW_Z          = -2.0

  ! Auto rejection composite score weights and percentile cut
  real, parameter :: AUTO_SNR_W      =  0.5
  real, parameter :: AUTO_PRESENCE_W =  0.02
  real, parameter :: AUTO_CONTRAST_W =  5.0
  real, parameter :: AUTO_THRESHOLD  =  10.0  ! percentile within top-50%

  ! Basic rejection hard limits
  real, parameter :: BASIC_MIN_POP   =  10.0
  real, parameter :: BASIC_MAX_RES   =  30.0

contains

  ! Rejects statistically anomalous class averages based on four image
  ! statistics, each assessed via robust z-scores:
  !   - Variance:        zero-variance (empty) classes, and extreme high/low
  !                      outliers beyond VARIANCE_MAX_Z / VARIANCE_MIN_Z.
  !   - Maximum pixel:   classes with anomalously bright pixels (MAXPIXEL_Z).
  !   - Centre/edge SNR: classes with poor signal localisation (CEN_EDGE_SNR_Z).
  !   - Skewness:        classes with anomalously negative skew (SKEW_Z).
  ! Sets the corresponding l_rejected elements to .true.; never clears them.
  subroutine reject_outliers( imgs, mskrad, l_rejected )
    type(image), allocatable, intent(inout) :: imgs(:)
    logical,     allocatable, intent(inout) :: l_rejected(:)
    real,                     intent(in)    :: mskrad
    real, allocatable :: stats(:), zscores(:)
    real    :: mm(2)
    integer :: icls, ncls
    if( .not. allocated(l_rejected)    ) THROW_HARD('l_rejected is not allocated')
    if( size(l_rejected) /= size(imgs) ) THROW_HARD('l_rejected does not match the number of images')
    ncls = size(imgs)

    ! Variance: reject empty classes and high/low outliers
    allocate(stats(ncls))
    do icls = 1, ncls
      stats(icls) = imgs(icls)%variance()
    end do
    zscores = robust_z_scores(stats)
    where( stats    == 0.0           ) l_rejected = .true.
    where( zscores  >  VARIANCE_MAX_Z) l_rejected = .true.
    where( zscores  <  VARIANCE_MIN_Z) l_rejected = .true.
    write(logfhandle,'(A,I4)') '>>> CLASSES REJECTED (VARIANCE=0)  :', count(stats   == 0.0)
    write(logfhandle,'(A,I4)') '>>> CLASSES REJECTED (VARIANCE>)   :', count(zscores >  VARIANCE_MAX_Z)
    write(logfhandle,'(A,I4)') '>>> CLASSES REJECTED (VARIANCE<)   :', count(zscores <  VARIANCE_MIN_Z)
    deallocate(stats, zscores)

    ! Maximum pixel value: reject classes with anomalously bright pixels
    allocate(stats(ncls))
    do icls = 1, ncls
      mm          = imgs(icls)%minmax()
      stats(icls) = mm(2)
    end do
    zscores = robust_z_scores(stats)
    where( zscores > MAXPIXEL_Z ) l_rejected = .true.
    write(logfhandle,'(A,I4)') '>>> CLASSES REJECTED (MAXPIX)      :', count(zscores > MAXPIXEL_Z)
    deallocate(stats, zscores)

    ! Centre/edge SNR: reject classes with poor signal localisation
    allocate(stats(ncls))
    do icls = 1, ncls
      stats(icls) = imgs(icls)%center_edge_snr(mskrad)
    end do
    zscores = robust_z_scores(stats)
    where( zscores < CEN_EDGE_SNR_Z ) l_rejected = .true.
    write(logfhandle,'(A,I4)') '>>> CLASSES REJECTED (CEN/EDGE SNR):', count(zscores < CEN_EDGE_SNR_Z)
    deallocate(stats, zscores)

    ! Skewness: reject classes with anomalously negative skew
    allocate(stats(ncls))
    do icls = 1, ncls
      stats(icls) = imgs(icls)%skew()
    end do
    zscores = robust_z_scores(stats)
    where( zscores < SKEW_Z ) l_rejected = .true.
    write(logfhandle,'(A,I4)') '>>> CLASSES REJECTED (SKEW)        :', count(zscores < SKEW_Z)
    deallocate(stats, zscores)
  end subroutine reject_outliers

  ! Ranks surviving (non-rejected) classes by a weighted composite quality
  ! score (SNR, signal presence, contrast) and rejects those scoring below the
  ! AUTO_THRESHOLD percentile within the top half of the distribution.
  ! Sets the corresponding l_rejected elements to .true.; never clears them.
  subroutine reject_auto( imgs, l_rejected )
    type(image), allocatable, intent(inout) :: imgs(:)
    logical,     allocatable, intent(inout) :: l_rejected(:)
    integer, allocatable :: indices(:), indices_all(:)
    real,    allocatable :: stats(:), sorted_stats(:)
    real    :: threshold
    integer :: icls, ncls
    if( .not. allocated(l_rejected)    ) THROW_HARD('l_rejected is not allocated')
    if( size(l_rejected) /= size(imgs) ) THROW_HARD('l_rejected does not match the number of images')
    ncls = size(imgs)

    ! Compute composite quality score for every class
    allocate(indices_all(ncls), stats(ncls))
    do icls = 1, ncls
      indices_all(icls) = icls
      stats(icls)       = calc_rejection_score(imgs(icls), AUTO_SNR_W, AUTO_PRESENCE_W, AUTO_CONTRAST_W)
    end do

    ! Restrict to surviving classes and rank scores descending
    indices = pack(indices_all, .not. l_rejected)
    stats   = pack(stats,       .not. l_rejected)
    call hpsort(stats, indices)
    call reverse(stats)
    call reverse(indices)

    ! Derive percentile threshold from the top half of the surviving classes
    sorted_stats = stats
    call hpsort(sorted_stats)
    call reverse(sorted_stats)
    threshold = percentile_value(sorted_stats(1 : size(sorted_stats) / 2), &
                                 size(sorted_stats) / 2, AUTO_THRESHOLD)

    ! Reject classes below the threshold
    do icls = 1, size(indices)
      if( stats(icls) < threshold ) l_rejected(indices(icls)) = .true.
    end do
    write(logfhandle,'(A,I4)') '>>> CLASSES REJECTED (AUTO)        :', count(stats < threshold)

    deallocate(indices, indices_all, stats, sorted_stats)

  contains

    ! Returns the value at percentile p (0–100) from a descending-sorted array.
    ! p=0 maps to the largest value (index 1); p=100 maps to the smallest.
    pure real function percentile_value( sorted_desc, n, p )
      integer, intent(in) :: n
      real,    intent(in) :: sorted_desc(n)
      real,    intent(in) :: p
      integer :: idx
      idx = max(1, min(n, int(real(n) * (100.0 - p) / 100.0) + 1))
      percentile_value = sorted_desc(idx)
    end function percentile_value

  end subroutine reject_auto

  ! Applies hard-limit rejection directly from a cls2D oris object.
  ! Rejects classes with population below BASIC_MIN_POP or resolution
  ! worse than BASIC_MAX_RES Å. Sets corresponding l_rejected elements
  ! to .true.; never clears them.
  subroutine reject_basic( cls2D, l_rejected )
    type(oris),            intent(in)    :: cls2D
    logical, allocatable,  intent(inout) :: l_rejected(:)
    real, allocatable :: stats(:)
    stats = cls2D%get_all('pop')
    if( size(stats) /= size(l_rejected) ) THROW_HARD('pop and rejected arrays differ in size')
    where( stats < BASIC_MIN_POP ) l_rejected = .true.
    deallocate(stats)
    stats = cls2D%get_all('res')
    if( size(stats) /= size(l_rejected) ) THROW_HARD('res and rejected arrays differ in size')
    where( stats > BASIC_MAX_RES ) l_rejected = .true.
    deallocate(stats)
  end subroutine reject_basic

  ! Computes a weighted composite quality score for a single class-average
  ! image. Higher scores indicate better quality. The score is a weighted sum
  ! of SNR, signal presence (clamped to zero from below), and contrast.
  real function calc_rejection_score( img, snr_w, presence_w, contrast_w )
    type(image), intent(inout) :: img
    real,        intent(in)    :: snr_w, presence_w, contrast_w
    calc_rejection_score = snr_w      * img%snr()                   &
                         + presence_w * max(0.0, img%presence())    &
                         + contrast_w * img%contrast()
  end function calc_rejection_score

end module simple_stream_microchunk_utils