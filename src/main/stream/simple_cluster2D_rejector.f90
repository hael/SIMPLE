!@descr: scoring and rejection utilities for microchunk-based 2D classification
!
! Provides the cluster2D_rejector type, which holds an array of class-average images
! and a per-class rejection mask. Rejection criteria are applied independently and
! accumulate: a class rejected by any criterion remains rejected.
!
! Rejection strategies:
!   reject_pop            — removes under-populated classes (< POP_PERCENT_THRESHOLD of total)
!   reject_res            — removes classes with poor FSC resolution (> RES_THRESHOLD Å)
!   reject_mask           — removes classes with Otsu foreground signal outside the mask disc
!   reject_brightness     — removes classes with abnormally high mean signal (contaminants)
!   reject_local_variance — removes structurally flat classes with low fg/bg local variance
!
! Typical usage:
!   call rejector%new(cavg_imgs, mskdiam)
!   call rejector%reject_pop(os_cls2D)
!   call rejector%reject_res(os_cls2D)
!   call rejector%reject_mask()
!   states = rejector%get_states()   ! 1 = kept, 0 = rejected

module simple_cluster2D_rejector
  use simple_stream_api
  use simple_oris,          only: oris
  use simple_image,         only: image
  use simple_stat,          only: robust_z_scores
  use simple_segmentation

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
  real, parameter :: MASK_THRESHOLD   = 10.0
  ! Reject classes whose foreground brightness robust z-score exceeds this.
  real, parameter :: BRIGHTNESS_ZSCORE_THRESH = 2.5
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
    procedure :: reject_brightness
    procedure :: reject_local_variance
  end type cluster2D_rejector

contains

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

  subroutine kill( self )
    class(cluster2D_rejector), intent(inout) :: self
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    if( allocated(self%imgs)       ) deallocate(self%imgs)
    if( allocated(self%l_rejected) ) deallocate(self%l_rejected)
    self%mskdiam = 0.0
    call timer_stop(t0, string('kill'))
  end subroutine kill

  function get_rejected( self ) result( rejected )
    class(cluster2D_rejector), intent(in) :: self
    logical, allocatable :: rejected(:)
    allocate(rejected, source=self%l_rejected)
  end function get_rejected

  function get_states( self ) result( states )
    class(cluster2D_rejector), intent(in) :: self
    integer, allocatable :: states(:)
    allocate(states(size(self%l_rejected)), source=1)
    where(self%l_rejected) states = 0
  end function get_states

  ! Rejects classes whose particle count falls below POP_PERCENT_THRESHOLD of the
  ! total, eliminating junk classes that attracted almost no particles.
  subroutine reject_pop( self, cls_oris )
    class(cluster2D_rejector), intent(inout) :: self
    type(oris),                intent(in)    :: cls_oris
    integer, allocatable    :: pop(:)
    integer                 :: i, noris, nrejected, threshold
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    noris = cls_oris%get_noris()
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) return
    pop = cls_oris%get_all_asint('pop')
    threshold = ceiling(sum(pop) * POP_PERCENT_THRESHOLD)
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> POPULATION REJECTION THRESHOLD :', threshold
    nrejected = 0
    do i = 1, noris
      if( pop(i) < threshold ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4)') '>>> POPULATION REJECTION OF CLASS :', i
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
    if( noris == 0 ) return
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

  ! Rejects classes that have significant Otsu-thresholded signal outside the mask disc,
  ! indicating ice contamination, carbon edge, or other non-particle artefacts.
  subroutine reject_mask( self )
    class(cluster2D_rejector), intent(inout) :: self
    real, allocatable       :: rmat(:,:,:)
    logical, allocatable    :: l_mask(:,:,:)
    type(image)             :: mask, img
    integer                 :: i, noris, nrejected, ldim(3)
    integer(timer_int_kind) :: t0
    real                    :: mskrad_px, smpd
    t0 = timer_start()
    noris = size(self%imgs)
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) return
    ldim      = self%imgs(1)%get_ldim()
    smpd      = self%imgs(1)%get_smpd()
    mskrad_px = (self%mskdiam / smpd) / 2.0
    allocate(rmat(ldim(1), ldim(2), ldim(3)))
    call mask%disc(ldim, smpd, mskrad_px)
    l_mask = mask%get_rmat() == 0.0   ! true outside the disc
    call mask%kill()
    nrejected = 0
    do i = 1, noris
      img = self%imgs(i)
      call img%bp(0., 10.)
      call otsu_img(img)
      call img%get_rmat_sub(rmat)
      call img%kill()
      if( sum(rmat, mask=l_mask) > MASK_THRESHOLD ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4)') '>>> MASK REJECTION OF CLASS :', i
      end if
    end do
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> # CLASSES REJECTED ON MASK :', nrejected
    deallocate(rmat, l_mask)
    call timer_stop(t0, string('reject_mask'))
  end subroutine reject_mask

  ! Rejects classes whose mean signal inside the Otsu foreground mask is a positive
  ! outlier (robust z-score > 2.5), flagging abnormally bright classes caused by gold
  ! beads, crystalline ice, or other high-contrast contaminants.
  subroutine reject_brightness( self )
    class(cluster2D_rejector), intent(inout) :: self
    real, allocatable       :: scores(:), zscores(:)
    real, allocatable       :: rmat(:,:,:), rmat_bin(:,:,:)
    type(image)             :: img
    integer                 :: i, noris, nrejected, ldim(3)
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    noris = size(self%imgs)
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) return
    ldim = self%imgs(1)%get_ldim()
    allocate(rmat(ldim(1), ldim(2), ldim(3)), rmat_bin(ldim(1), ldim(2), ldim(3)))
    allocate(scores(noris))
    do i = 1, noris
      img = self%imgs(i)
      call img%get_rmat_sub(rmat)      ! original signal, captured before filtering
      call img%bp(0., 30.)
      call otsu_img(img)
      call img%get_rmat_sub(rmat_bin)  ! binary foreground mask
      call img%kill()
      scores(i) = sum(rmat, mask=(rmat_bin > 0.5))
    end do
    zscores = robust_z_scores(scores)
    nrejected = 0
    do i = 1, noris
      if( zscores(i) > BRIGHTNESS_ZSCORE_THRESH ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4,F8.4)') '>>> BRIGHTNESS REJECTION OF CLASS :', i, zscores(i)
      end if
    end do
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> # CLASSES REJECTED ON BRIGHTNESS :', nrejected
    deallocate(rmat, rmat_bin, scores, zscores)
    call timer_stop(t0, string('reject_brightness'))
  end subroutine reject_brightness

  ! Rejects classes where local variance is anomalously low in both the foreground and
  ! background simultaneously, indicating a structurally flat class average with no real
  ! signal (e.g. pure noise or empty classes). Requiring low z-scores in BOTH regions
  ! avoids rejecting valid low-variance particles on a single weak signal.
  subroutine reject_local_variance( self )
    class(cluster2D_rejector), intent(inout) :: self
    real, allocatable       :: scores_inside(:), scores_outside(:)
    real, allocatable       :: zscores_inside(:), zscores_outside(:)
    real, allocatable       :: bin_mask(:,:,:)
    type(image)             :: img, img1
    integer                 :: i, noris, nrejected, ldim(3)
    integer(timer_int_kind) :: t0
    t0 = timer_start()
    noris = size(self%imgs)
    if( noris /= size(self%l_rejected) ) THROW_HARD("number cls oris does not match rejected")
    if( noris == 0 ) return
    ldim = self%imgs(1)%get_ldim()
    allocate(bin_mask(ldim(1), ldim(2), 1))
    allocate(scores_inside(noris), scores_outside(noris))
    do i = 1, noris
      img  = self%imgs(i)
      call img%bp(0., 10.)
      img1 = img
      call otsu_img(img1)
      call img1%get_rmat_sub(bin_mask)
      call img%loc_var_masked(bin_mask(:,:,1), 10, scores_outside(i), scores_inside(i))
      call img%kill()
      call img1%kill()
    end do
    zscores_inside  = robust_z_scores(scores_inside)
    zscores_outside = robust_z_scores(scores_outside)
    nrejected = 0
    do i = 1, noris
      write(*,*) i, zscores_inside(i), zscores_outside(i)
      if( (zscores_inside(i) < LOCVAR_STRONG_THRESH .and. zscores_outside(i) < LOCVAR_WEAK_THRESH) .or. &
          (zscores_inside(i) < LOCVAR_WEAK_THRESH   .and. zscores_outside(i) < LOCVAR_STRONG_THRESH) ) then
        self%l_rejected(i) = .true.
        nrejected = nrejected + 1
        if( DEBUG ) write(logfhandle,'(A,I4,4F8.4)') '>>> LOCAL VARIANCE REJECTION OF CLASS :', &
          i, scores_outside(i), scores_inside(i), zscores_outside(i), zscores_inside(i)
      end if
    end do
    if( DEBUG ) write(logfhandle,'(A,I4)') '>>> # CLASSES REJECTED ON LOCAL VARIANCE :', nrejected
    deallocate(bin_mask, scores_inside, scores_outside, zscores_inside, zscores_outside)
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
