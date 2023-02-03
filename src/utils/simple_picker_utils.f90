module simple_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: picker_utils
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_TARGET1 = 4.0, SMPD_TARGET2 = 2.0, LAMBDA = 5.
integer, parameter :: OFFSET  = 3
logical, parameter :: L_WRITE = .true.

type picker_utils
    private
    integer                  :: ldim_raw(3) = 0 , ldim_shrink1(3) = 0 , ldim_shrink2(3) = 0
    integer                  :: ldim_box(3) = 0 , ldim_box1(3)    = 0 , ldim_box2(3)    = 0
    real                     :: smpd_raw    = 0., smpd_shrink1    = 0., smpd_shrink2    = 0.
    integer                  :: nboxes1 = 0, nboxes2 = 0, nx1 = 0, ny1 = 0, nx2 = 0, ny2 = 0
    type(image)              :: mic_shrink1, mic_shrink2
    type(image), pointer     :: mic_raw => null()
    type(image), allocatable :: boximgs1(:), boximgs2(:)
    integer,     allocatable :: positions1(:,:)
    logical                  :: pos_set = .false.
  contains
    procedure          :: set_mics
    procedure          :: bin_mic_shrink1
    procedure, private :: set_positions_1
    procedure, private :: set_positions_2
    generic            :: set_positions => set_positions_1, set_positions_2
    procedure, private :: set_pos_priv
    procedure, private :: extract_boximgs1
end type picker_utils

contains

    subroutine set_mics( self, mic, smpd, bp_lp )
        class(picker_utils),  intent(inout) :: self
        class(image), target, intent(in)    :: mic
        real,                 intent(in)    :: smpd     !< sampling distance in A
        real, optional,       intent(in)    :: bp_lp(2) !< high- and low-pass limits in A
        real :: scale1, scale2
        ! set raw micrograph info
        self%ldim_raw = mic%get_ldim()
        if( self%ldim_raw(3) /= 1 ) THROW_HARD('Only for 2D images')
        self%smpd_raw = smpd
        self%mic_raw  => mic
        ! shrink micrograph
        scale1 = self%smpd_raw / SMPD_TARGET1
        scale2 = self%smpd_raw / SMPD_TARGET2
        self%ldim_shrink1(1) = round2even(real(self%ldim_raw(1)) * scale1)
        self%ldim_shrink1(2) = round2even(real(self%ldim_raw(2)) * scale1)
        self%ldim_shrink1(3) = 1
        self%ldim_shrink2(1) = round2even(real(self%ldim_raw(1)) * scale2)
        self%ldim_shrink2(2) = round2even(real(self%ldim_raw(2)) * scale2)
        self%ldim_shrink2(3) = 1
        call self%mic_shrink1%new(self%ldim_shrink1, SMPD_TARGET1)
        call self%mic_shrink2%new(self%ldim_shrink2, SMPD_TARGET2)
        call self%mic_shrink1%set_ft(.true.)
        call self%mic_shrink2%set_ft(.true.)
        call self%mic_raw%fft
        call self%mic_raw%clip(self%mic_shrink1)
        call self%mic_raw%clip(self%mic_shrink2)
        if( present(bp_lp) )then
            call self%mic_shrink1%bp(bp_lp(1),bp_lp(2))
            call self%mic_shrink2%bp(bp_lp(1),bp_lp(2))
        endif
        call self%mic_raw%ifft
        call self%mic_shrink1%ifft
        call self%mic_shrink2%ifft
        if( L_WRITE )then
            call self%mic_shrink1%write('mic_shrink1.mrc')
            call self%mic_shrink2%write('mic_shrink2.mrc')
        endif
    end subroutine set_mics

    ! type stats_struct
    !     real :: avg  = 0.
    !     real :: med  = 0.
    !     real :: sdev = 0.
    !     real :: maxv = 0
    !     real :: minv = 0.
    ! end type stats_struct

    subroutine bin_mic_shrink1( self, bp_lp )
        use simple_segmentation
        use simple_tvfilter, only: tvfilter
        class(picker_utils), intent(in) :: self
        real,                intent(in) :: bp_lp(2) !< high- and low-pass limits in A
        real, allocatable  :: vec(:), vec2(:)
        type(tvfilter)     :: tvf
        type(image)        :: img, img_copy
        type(stats_struct) :: stats_ptcl, stats_bg
        logical            :: lmsk(self%ldim_shrink1(1),self%ldim_shrink1(2),1)
        real :: t
        call img%copy(self%mic_shrink1) ! remplace with copy_fast in optimized code (pre-construct)
        call img%fft
        call img%bp(bp_lp(1),bp_lp(2))
        call tvf%new()
        call tvf%apply_filter(img, LAMBDA)
        call tvf%kill
        call img%ifft
        call img%write('tv_mic.mrc')
        vec = img%serialize()
        call otsu(size(vec), vec, t)
        call img_copy%copy(img)
        call img_copy%binarize(t)
        call img_copy%write('otsu_mic.mrc')
        vec = self%mic_shrink1%serialize()
        call calc_stats(vec, stats_bg,   mask=vec >= t) ! background stats, assuming black particle density
        call calc_stats(vec, stats_ptcl, mask=vec <  t) ! particle   stats, assuming black particle density

        print *, 'background #pix/avg/sdev ', count(vec >= t), stats_bg%avg,   stats_bg%sdev
        print *, 'particle   #pix/avg/sdev ', count(vec <  t), stats_ptcl%avg, stats_ptcl%sdev
        print *, 'std_mean_diff ', std_mean_diff(stats_bg%avg, stats_ptcl%avg, stats_bg%sdev, stats_ptcl%sdev)

        ! results from filtered
        ! background #pix/avg/sdev       435326  0.964909315       7.91966170E-03
        ! particle   #pix/avg/sdev       261674  0.940358341       1.32814981E-02
        ! std_mean_diff    2.24531102

        ! results from raw
        ! background #pix/avg/sdev       368903  0.995647132       3.23507152E-02
        ! particle   #pix/avg/sdev       328097  0.910767972       3.42143402E-02
        ! std_mean_diff    2.54926300
        
        call img%kill
        call img_copy%kill
    end subroutine bin_mic_shrink1

    subroutine set_positions_1( self, maxdiam )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: maxdiam !< maximum diameter in A
        ! set logical dimensions of boxes
        self%ldim_box(1)  = round2even(maxdiam / self%smpd_raw)
        self%ldim_box(2)  = self%ldim_box(1)
        self%ldim_box(3)  = 1
        self%ldim_box1(1) = round2even(maxdiam / SMPD_TARGET1)
        self%ldim_box1(2) = self%ldim_box1(1)
        self%ldim_box1(3) = 1
        self%ldim_box2(1) = round2even(maxdiam / SMPD_TARGET2)
        self%ldim_box2(2) = self%ldim_box2(1)
        self%ldim_box2(3) = 1
        call self%set_pos_priv
    end subroutine set_positions_1

    subroutine set_positions_2( self, box_raw, box12 )
        class(picker_utils), intent(inout) :: self
        real,                intent(in)    :: box_raw, box12(2)
        ! set logical dimensions of boxes
        self%ldim_box(1)  = box_raw
        self%ldim_box(2)  = self%ldim_box(1)
        self%ldim_box(3)  = 1
        self%ldim_box1(1) = box12(1)
        self%ldim_box1(2) = self%ldim_box1(1)
        self%ldim_box1(3) = 1
        self%ldim_box2(1) = box12(2)
        self%ldim_box2(2) = self%ldim_box2(1)
        self%ldim_box2(3) = 1
        call self%set_pos_priv
    end subroutine set_positions_2

    subroutine set_pos_priv( self )
        class(picker_utils), intent(inout) :: self
        integer :: xind, yind
        ! set # pixels in x/y for both box sizes
        self%nx1 = self%ldim_shrink1(1) - self%ldim_box1(1)
        self%ny1 = self%ldim_shrink1(2) - self%ldim_box1(2)
        self%nx2 = self%ldim_shrink2(1) - self%ldim_box2(1)
        self%ny2 = self%ldim_shrink2(2) - self%ldim_box2(2)
        ! count # boxes
        self%nboxes1 = 0
        do xind = 0,self%nx1,OFFSET
            do yind = 0,self%ny1,OFFSET
                self%nboxes1 = self%nboxes1 + 1
            end do
        end do
        ! allocate and set positions1 
        if( allocated(self%positions1) ) deallocate(self%positions1)
        allocate(self%positions1(self%nboxes1,2), source=0)
        self%nboxes1 = 0
        do xind = 0,self%nx1,OFFSET
            do yind = 0,self%ny1,OFFSET
                self%nboxes1 = self%nboxes1 + 1
                self%positions1(self%nboxes1,:) = [xind,yind]
            end do
        end do
        self%pos_set = .true.
    end subroutine set_pos_priv

    subroutine extract_boximgs1( self )
        class(picker_utils), intent(inout) :: self
        integer :: ibox
        logical :: outside
        if( .not. self%pos_set ) THROW_HARD('positions need to be set before constructing boximgs1')
        if( allocated(self%boximgs1) )then
            do ibox = 1,self%nboxes1
                call self%boximgs1(ibox)%kill
            end do
            deallocate(self%boximgs1)
        endif
        allocate(self%boximgs1(self%nboxes1))
        !$omp parallel do schedule(static) default(shared) private(ibox,outside) proc_bind(close)
        do ibox = 1,self%nboxes1
            call self%boximgs1(ibox)%new(self%ldim_box1, SMPD_TARGET1)
            call self%mic_shrink1%window_slim(self%positions1(ibox,:), self%ldim_box1(1), self%boximgs1(ibox), outside)
        end do
        !$omp end parallel do   
    end subroutine extract_boximgs1

end module simple_picker_utils
