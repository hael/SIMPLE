!@descr: scoring and rejection utilities for microchunk-based 2D classification
!
! Provides the cluster2D_rejector type, which holds an array of class-average images
! and a per-class rejection mask. Rejection criteria are applied independently and
! accumulate: a class rejected by any criterion remains rejected.
!
! Rejection strategies:
!   reject_pop            — removes under-populated classes (< POP_PERCENT_THRESHOLD of total);
!                           optional thres overrides POP_PERCENT_THRESHOLD for that call
!   reject_res            — removes classes with poor FSC resolution (> RES_THRESHOLD Å)
!   reject_mask           — removes classes with Otsu foreground signal outside the mask disc
!   reject_local_variance — removes structurally flat classes: zero inside/outside variance
!                           scores are rejected unconditionally; remaining classes are rejected
!                           when robust z-scores in both regions fall below the dual threshold;
!                           optional strong_thresh/weak_thresh override the defaults
!
! Constants:
!   POP_PERCENT_THRESHOLD — minimum population fraction to keep a class          (0.005)
!   RES_THRESHOLD         — maximum FSC resolution to keep a class, Å            (40.0)
!   MASK_THRESHOLD        — maximum pixels outside mask disc for the largest CC  (10.0)
!   LOCVAR_STRONG_THRESH  — strong z-score cut for local-variance rejection       (-0.5)
!   LOCVAR_WEAK_THRESH    — weak z-score cut for local-variance rejection         (-0.1)
!
! Typical usage:
!   call rejector%new(cavg_imgs, mskdiam)
!   call rejector%reject_pop(os_cls2D)
!   call rejector%reject_res(os_cls2D)
!   call rejector%reject_mask()
!   call rejector%reject_local_variance()
!   l_rejected = rejector%get_rejected()   ! .true. = rejected, .false. = kept
!   states     = rejector%get_states()     ! 0 = rejected, 1 = kept

module simple_cluster2D_rejector
  use simple_stream_api
  use simple_oris,          only: oris
  use simple_image,         only: image
  use simple_stat,          only: robust_z_scores
  use simple_image_msk,     only: density_inoutside_mask
  use simple_image_bin,     only: image_bin
  use simple_segmentation,  only: otsu_img

  implicit none
  private
  public :: cluster2D_rejector
#include "simple_local_flags.inc"
  logical, parameter :: DEBUG = .true.

  ! Reject classes whose population is below this fraction of the total particle count.
  real, parameter :: POP_PERCENT_THRESHOLD = 0.005
  ! Reject classes whose FSC resolution estimate is worse (higher Angstrom value) than this.
  real, parameter :: RES_THRESHOLD         = 40.0
  ! Reject classes where the number of foreground pixels outside the mask disc exceeds this.
  real, parameter :: MASK_THRESHOLD        = 10.0
  ! Reject classes where local variance z-scores are below these thresholds in both regions.
  real, parameter :: LOCVAR_STRONG_THRESH  = -0.5
  real, parameter :: LOCVAR_WEAK_THRESH    = -0.1

  type :: cluster2D_rejector
    private
    logical,     allocatable :: l_rejected(:)
    type(image), allocatable :: imgs(:)
    real                     :: mskdiam = 0.0
  contains
    procedure :: new
    procedure :: kill
    procedure :: get_rejected
    procedure :: get_states
    procedure :: reject_pop
    procedure :: reject_res
    procedure :: reject_mask
    procedure :: reject_local_variance
  end type cluster2D_rejector

contains

  ! Initialises the rejector from an allocated image array and a mask diameter in
  ! Angstroms. Kills any previous state, copies the images, and clears all
  ! rejection flags so each class starts as kept.
  subroutine new( self, imgs, mskdiam )
    class(cluster2D_rejector),    intent(inout) :: self
    type(image),     allocatable, intent(in)    :: imgs(:)
    real,                         intent(in)    :: mskdiam
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    call self%kill()
    self%mskdiam = mskdiam
    allocate(self%imgs, source=imgs)
    allocate(self%l_rejected(size(imgs)), source=.false.)
    call timer_stop(t0, string('new'))
  end subroutine new

  ! Deallocates the image array and rejection mask and resets mskdiam to 0.
  subroutine kill( self )
    class(cluster2D_rejector), intent(inout) :: self
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    if( allocated(self%imgs)       ) deallocate(self%imgs)
    if( allocated(self%l_rejected) ) deallocate(self%l_rejected)
    self%mskdiam = 0.0
    call timer_stop(t0, string('kill'))
  end subroutine kill

  ! Returns a copy of the per-class rejection mask (.true. = rejected, .false. = kept).
  function get_rejected( self ) result( rejected )
    class(cluster2D_rejector), intent(in) :: self
    logical, allocatable :: rejected(:)
    allocate(rejected, source=self%l_rejected)
  end function get_rejected

  ! Returns per-class states as integers (0 = rejected, 1 = kept), matching the
  ! convention used by sp_project os_cls2D.
  function get_states( self ) result( states )
    class(cluster2D_rejector), intent(in) :: self
    integer, allocatable :: states(:)
    allocate(states(size(self%l_rejected)), source=1)
    where(self%l_rejected) states = 0
  end function get_states

  ! Rejects classes whose particle count falls below POP_PERCENT_THRESHOLD of the
  ! total, eliminating junk classes that attracted almost no particles. Optional
  ! thres overrides POP_PERCENT_THRESHOLD for this call only.
  subroutine reject_pop( self, cls_oris, thres )
    class(cluster2D_rejector), intent(inout) :: self
    type(oris),                intent(in)    :: cls_oris
    real,    optional,         intent(in)    :: thres
    integer, allocatable    :: pop(:)
    integer                 :: i, noris, nrejected, threshold
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    noris = cls_oris%get_noris()
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) then
      call timer_stop(t0, string('reject_pop'))
      return
    end if
    pop = cls_oris%get_all_asint('pop')
    if( present(thres) ) then
      threshold = ceiling(sum(pop) * thres)
    else
      threshold = ceiling(sum(pop) * POP_PERCENT_THRESHOLD)
    end if
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> POPULATION REJECTION THRESHOLD :', threshold
    nrejected = 0
    do i = 1, noris
      if( pop(i) < threshold ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4,I4)') '>>> POPULATION REJECTION OF CLASS :', i, pop(i)
      end if
    end do
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> # CLASSES REJECTED ON POPULATION :', nrejected
    deallocate(pop)
    call timer_stop(t0, string('reject_pop'))
  end subroutine reject_pop

  ! Rejects classes with an FSC resolution estimate worse than RES_THRESHOLD Angstroms,
  ! indicating the class average contains mostly noise.
  subroutine reject_res( self, cls_oris )
    class(cluster2D_rejector), intent(inout) :: self
    type(oris),                intent(in)    :: cls_oris
    real, allocatable       :: res(:)
    integer                 :: i, noris, nrejected
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    noris = cls_oris%get_noris()
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) then
      call timer_stop(t0, string('reject_res'))
      return
    end if
    res = cls_oris%get_all('res')
    if( DEBUG ) write(logfhandle,'(A,F4.1)') '>>> RESOLUTION REJECTION THRESHOLD :', RES_THRESHOLD
    nrejected = 0
    do i = 1, noris
      if( res(i) > RES_THRESHOLD ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4)') '>>> RESOLUTION REJECTION OF CLASS :', i
      end if
    end do
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> # CLASSES REJECTED ON RESOLUTION :', nrejected
    deallocate(res)
    call timer_stop(t0, string('reject_res'))
  end subroutine reject_res

  ! Rejects classes where the largest Otsu-thresholded connected component (CC) either
  ! has its mass centre outside the mask disc, or has more than MASK_THRESHOLD pixels
  ! extending beyond it — indicating artefact signal (carbon edge, ice, detector glow).
  ! CCs spanning the full image dimension are first pruned as spurious background blobs.
  subroutine reject_mask( self )
    class(cluster2D_rejector), intent(inout) :: self
    real,    allocatable                     :: ccsizes(:), rmat_cc(:,:,:), rmat_msk(:,:,:)
    type(image_bin)                          :: img_bin, cc_img, img_mask
    integer                                  :: i, j, noris, nrejected
    integer                                  :: ldim(3), loc, nccs, nccs_updated
    integer(timer_int_kind)                  :: t0
    real                                     :: smpd, cc_diam, rad_px, xy(2)
    logical                                  :: l_rejected
    t0 = timer_start()
    noris = size(self%imgs)
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) then
      call timer_stop(t0, string('reject_mask'))
      return
    end if
    ldim      = self%imgs(1)%get_ldim()
    smpd      = self%imgs(1)%get_smpd()
    rad_px    = (self%mskdiam / smpd) / 2.0
    nrejected = 0
    ! construct image buffers once; img_bin%copy reinitialises content each iteration
    call img_bin%new_bimg(ldim, smpd, wthreads=.false.)
    call cc_img%new_bimg(ldim,  smpd, wthreads=.false.)
    call img_mask%disc(ldim,    smpd, rad_px)
    ! pre-allocate rmat buffers; avoids two heap allocations per iteration from get_rmat()
    allocate(rmat_cc(ldim(1), ldim(2), ldim(3)), rmat_msk(ldim(1), ldim(2), ldim(3)))
    call img_mask%get_rmat_sub(rmat_msk)
    do i = 1, noris
      call img_bin%copy(self%imgs(i))
      call img_bin%zero_edgeavg()
      call img_bin%bp(0., 30.0)
      call otsu_img(img_bin)
      call img_bin%set_imat()
      call img_bin%find_ccs(cc_img)
      call cc_img%get_nccs(nccs)
      ! prune CCs whose diameter spans the full image — these are background blobs,
      ! not particles; labels are stable across calls because update=.false. defers reordering
      nccs_updated = nccs
      do j = 1, nccs
        call cc_img%diameter_cc(j, cc_diam)
        if( cc_diam > ldim(1) ) then
          call cc_img%elim_cc(j, update=.false.)
          nccs_updated = nccs_updated - 1
        end if
      end do
      ! reject if any CC's centroid lies outside the mask radius, or if the
      ! largest CC has more than MASK_THRESHOLD pixels outside the disc
      if( nccs_updated > 0 ) then
        l_rejected = .false.
        call cc_img%order_ccs()
        call cc_img%update_img_rmat()
        call cc_img%get_nccs(nccs)
        do j = 1, nccs
          call cc_img%masscen_cc(j, xy)
          if( sqrt(xy(1)**2 + xy(2)**2) > rad_px ) l_rejected = .true.
        end do
        ccsizes = cc_img%size_ccs()
        loc     = maxloc(ccsizes, dim=1)
        deallocate(ccsizes)   ! done with sizes; deallocate here, not conditionally after the loop
        call cc_img%cc2bin(loc)
        call cc_img%get_rmat_sub(rmat_cc)
        if( count(rmat_cc - rmat_msk > 0.0) > MASK_THRESHOLD ) l_rejected = .true.
      else
        l_rejected = .true.
      end if
      if( l_rejected ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4)') '>>> MASK REJECTION OF CLASS :', i
      end if
    end do
    call img_bin%kill_bimg()
    call cc_img%kill_bimg()
    call img_mask%kill_bimg()
    deallocate(rmat_cc, rmat_msk)
    if( DEBUG ) write(logfhandle,'(A,I8)') '>>> # CLASSES REJECTED ON MASK :', nrejected
    call timer_stop(t0, string('reject_mask'))
  end subroutine reject_mask

  ! Rejects classes where local variance is anomalously low in both the foreground and
  ! background simultaneously, indicating a structurally flat class average with no real
  ! signal (e.g. pure noise or empty classes). Requiring low z-scores in BOTH regions
  ! avoids rejecting valid low-variance particles on a single weak signal. Classes with
  ! zero scores in both regions are rejected unconditionally and excluded from z-score
  ! computation so they do not skew the distribution for the remaining classes.
  subroutine reject_local_variance( self, strong_thresh, weak_thresh )
    class(cluster2D_rejector), intent(inout) :: self
    real,    optional,         intent(in)    :: strong_thresh, weak_thresh
    real,    allocatable    :: scores_inside(:), scores_outside(:)
    real,    allocatable    :: zscores_inside(:), zscores_outside(:)
    real,    allocatable    :: bin_mask(:,:,:)
    logical, allocatable    :: l_zero(:)
    type(image)             :: img, img1
    integer                 :: i, noris, nrejected, ldim(3)
    integer(timer_int_kind) :: t0
    real                    :: strong_threshold, weak_threshold
    t0 = timer_start()
    noris = size(self%imgs)
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) then
      call timer_stop(t0, string('reject_local_variance'))
      return
    end if
    ldim = self%imgs(1)%get_ldim()
    allocate(bin_mask(ldim(1), ldim(2), 1))
    allocate(scores_inside(noris), scores_outside(noris))
    do i = 1, noris
      img  = self%imgs(i)
      call img%zero_edgeavg()
      call img%bp(0., 10.)
      img1 = img
      call otsu_img(img1)
      call img1%get_rmat_sub(bin_mask)
      call img%loc_var_masked(bin_mask(:,:,1), 10, scores_outside(i), scores_inside(i))
      call img%kill()
      call img1%kill()
    end do
    strong_threshold = LOCVAR_STRONG_THRESH
    weak_threshold   = LOCVAR_WEAK_THRESH
    if( present(strong_thresh) ) strong_threshold = strong_thresh
    if( present(weak_thresh)   ) weak_threshold   = weak_thresh
    ! Unconditionally reject classes with zero variance in both regions and exclude
    ! them from z-scoring so blank/empty classes do not skew the distribution.
    allocate(l_zero(noris),         source=.false.)
    allocate(zscores_inside(noris), source=0.0)
    allocate(zscores_outside(noris),source=0.0)
    nrejected = 0
    do i = 1, noris
      if( scores_inside(i) == 0.0 .and. scores_outside(i) == 0.0 ) then
        l_zero(i)          = .true.
        self%l_rejected(i) = .true.
        nrejected          = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4)') '>>> ZERO-VARIANCE REJECTION OF CLASS :', i
      end if
    end do
    if( count(.not. l_zero) > 1 ) then
      zscores_inside( pack([(i,i=1,noris)], .not. l_zero) ) = &
        robust_z_scores(pack(scores_inside,  .not. l_zero))
      zscores_outside(pack([(i,i=1,noris)], .not. l_zero) ) = &
        robust_z_scores(pack(scores_outside, .not. l_zero))
    end if
    do i = 1, noris
      if( l_zero(i) ) cycle
      if( (zscores_inside(i) < strong_threshold .and. zscores_outside(i) < weak_threshold  ) .or. &
          (zscores_inside(i) < weak_threshold   .and. zscores_outside(i) < strong_threshold) ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4,4F8.4)') '>>> LOCAL VARIANCE REJECTION OF CLASS :', &
          i, scores_outside(i), scores_inside(i), zscores_outside(i), zscores_inside(i)
      end if
    end do
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> # CLASSES REJECTED ON LOCAL VARIANCE :', nrejected
    deallocate(bin_mask, scores_inside, scores_outside, zscores_inside, zscores_outside, l_zero)
    call timer_stop(t0, string('reject_local_variance'))
  end subroutine reject_local_variance

  integer(timer_int_kind) function timer_start()
    timer_start = tic()
  end function timer_start

  subroutine timer_stop( t0, routine_name )
    integer(timer_int_kind), intent(in) :: t0
    type(string),            intent(in) :: routine_name
    if( .not. DEBUG ) return
    write(logfhandle,'(A,A,A,F8.1)') 'cluster2D_rejector->', routine_name%to_char(), ' execution time:', toc(t0)
    call flush(logfhandle)
  end subroutine timer_stop

end module simple_cluster2D_rejector
