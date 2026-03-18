!@descr: utilities for microchunk-based 2D clustering in stream
module simple_stream_microchunk_utils
use simple_stream_api
use simple_image,         only: image
use simple_stat,          only: robust_z_scores
use simple_srch_sort_loc, only: hpsort, reverse

implicit none

#include "simple_local_flags.inc"

contains

    subroutine reject_outliers( imgs, mskrad, l_rejected )
        type(image), allocatable, intent(inout) :: imgs(:)
        logical,     allocatable, intent(inout) :: l_rejected(:)
        real,                        intent(in) :: mskrad 
        real,          parameter :: VARIANCE_MAX_Z_THRESHOLD = 3.0
        real,          parameter :: VARIANCE_MIN_Z_THRESHOLD = -1.0
        real,          parameter :: MAXPIXEL_Z_THRESHOLD     = 3.0
        real,          parameter :: CEN_EDGE_SNR_Z_THRESHOLD = -2.0
        real,          parameter :: SKEW_Z_THRESHOLD         = -2.0
        integer,     allocatable :: indices(:), indices_all(:)
        real,        allocatable :: stats(:), zscores(:)
        integer :: icls
        real    :: mm(2)
        if( .not. allocated(l_rejected)    ) THROW_HARD('l_rejected is not allocated')
        if( size(l_rejected) /= size(imgs) ) THROW_HARD('l_rejected does not match the number of images')
        allocate(indices_all(size(imgs)))
        do icls=1, size(imgs)
            indices_all(icls) = icls
        end do
        ! reject where variance == 0.0 (empty classes)
        allocate(stats(size(imgs)))
        do icls=1, size(imgs)
            stats(icls) = imgs(icls)%variance()
        end do
        zscores = robust_z_scores(stats)
        where( stats == 0.0                      ) l_rejected = .true.
        where( zscores > VARIANCE_MAX_Z_THRESHOLD) l_rejected = .true.
        where( zscores < VARIANCE_MIN_Z_THRESHOLD) l_rejected = .true.
        write(logfhandle, '(A,I4)') '>>> CLASSES REJECTED (VARIANCE=0)  :', count(stats == 0.0)
        write(logfhandle, '(A,I4)') '>>> CLASSES REJECTED (VARIANCE>)   :', count(zscores > VARIANCE_MAX_Z_THRESHOLD) 
        write(logfhandle, '(A,I4)') '>>> CLASSES REJECTED (VARIANCE<)   :', count(zscores < VARIANCE_MIN_Z_THRESHOLD) 
        deallocate(stats, zscores)
        ! reject where max pixel value z score > MAXPIXEL_Z_THRESHOLD
        allocate(stats(size(imgs)))
        do icls=1, size(imgs)
            mm = imgs(icls)%minmax()
            stats(icls) = mm(2)
        end do
        zscores = robust_z_scores(stats)
        where( zscores > MAXPIXEL_Z_THRESHOLD ) l_rejected = .true.
        write(logfhandle, '(A,I4)') '>>> CLASSES REJECTED (MAXPIX)      :', count(zscores > MAXPIXEL_Z_THRESHOLD)
        deallocate(stats, zscores)
        ! reject where center/edge snr value z score < CEN_EDGE_SNR_Z_THRESHOLD
        allocate(stats(size(imgs)))
        do icls=1, size(imgs)
            stats(icls) = imgs(icls)%center_edge_snr(mskrad)
            write(*,*) icls, stats(icls), imgs(icls)%center_edge_snr(), imgs(icls)%snr(), imgs(icls)%variance(), imgs(icls)%presence(), imgs(icls)%contrast()
        end do
        zscores = robust_z_scores(stats)
        where( zscores < CEN_EDGE_SNR_Z_THRESHOLD ) l_rejected = .true.
        write(logfhandle, '(A,I4)') '>>> CLASSES REJECTED (CEN/EDGE SNR):', count(zscores < CEN_EDGE_SNR_Z_THRESHOLD)
        deallocate(stats, zscores)
        ! reject where skew value z score < CEN_EDGE_SNR_Z_THRESHOLD
        allocate(stats(size(imgs)))
        do icls=1, size(imgs)
            stats(icls) = imgs(icls)%skew()
        end do
        zscores = robust_z_scores(stats)
         write(*,*) zscores
        where( zscores < SKEW_Z_THRESHOLD ) l_rejected = .true.
        write(logfhandle, '(A,I4)') '>>> CLASSES REJECTED (SKEW)        :', count(zscores < SKEW_Z_THRESHOLD)
        deallocate(stats, zscores)
    end subroutine reject_outliers

    subroutine reject_auto( imgs, l_rejected )
        type(image), allocatable, intent(inout) :: imgs(:)
        logical,     allocatable, intent(inout) :: l_rejected(:)
        REAL, PARAMETER :: DEFAULT_SNR_W       = 0.5
        REAL, PARAMETER :: DEFAULT_PRESENCE_W  = 0.02
        REAL, PARAMETER :: DEFAULT_CONTRAST_W  = 5.0
        REAL, PARAMETER :: DEFAULT_THRESHOLD   = 10.0   ! percentile of top-50% retained
        integer,     allocatable :: indices(:), indices_all(:)
        real,        allocatable :: stats(:), sorted_stats(:)
        integer :: icls
        real    :: threshold
        if( .not. allocated(l_rejected)    ) THROW_HARD('l_rejected is not allocated')
        if( size(l_rejected) /= size(imgs) ) THROW_HARD('l_rejected does not match the number of images')
        allocate(indices_all(size(imgs)))
        do icls=1, size(imgs)
            indices_all(icls) = icls
        end do
        allocate(stats(size(imgs)))
        do icls=1, size(imgs)
            stats(icls) = calc_rejection_score(imgs(icls), DEFAULT_SNR_W, DEFAULT_PRESENCE_W, DEFAULT_CONTRAST_W)
        end do
        ! mask stats
        indices = pack(indices_all, .not.l_rejected)
        stats   = pack(stats,       .not.l_rejected)
        ! Rank scores (descending) to identify top-50%
        call hpsort(stats, indices)
        call reverse(stats)
        call reverse(indices)
        ! Compute percentile threshold within the top-50% of classes
        sorted_stats = stats
       ! call sort_descending(sorted_stats, size(stats))
        call hpsort(sorted_stats)
        call reverse(sorted_stats)
        threshold = percentile_value(sorted_stats(1 : size(sorted_stats)/2), size(sorted_stats)/2, DEFAULT_THRESHOLD)
        ! Apply mask
        where( stats < threshold ) l_rejected(indices) = .true.
        write(logfhandle, '(A,I4)') '>>> CLASSES REJECTED (AUTO)        :', count(stats < threshold)
        deallocate(indices, indices_all)
        deallocate(stats, sorted_stats)

    contains

        PURE FUNCTION percentile_value(sorted_desc, n, p) RESULT(val)
            INTEGER, INTENT(IN) :: n
            REAL,    INTENT(IN) :: sorted_desc(n)
            REAL,    INTENT(IN) :: p          ! percentile, 0–100
            REAL                :: val
            INTEGER             :: idx
            ! p=0 → largest (index 1), p=100 → smallest (index n)
            idx = MAX(1, MIN(n, INT(REAL(n) * (100.0 - p) / 100.0) + 1))
            val = sorted_desc(idx)
        END FUNCTION percentile_value

    end subroutine reject_auto

    real function calc_rejection_score( img, snr_w, presence_w, contrast_w )
        type(image), intent(inout) :: img
        real,           intent(in) :: snr_w, presence_w, contrast_w
        real :: snr, presence, contrast
        ! calc image stats
        snr      = img%snr()
        presence = img%presence()
        contrast = img%contrast()
        ! combine into rejection score
        calc_rejection_score = snr_w * snr &
         + presence_w * max(0.0, presence) &
         + contrast_w * contrast
    end function calc_rejection_score
   

end module simple_stream_microchunk_utils
