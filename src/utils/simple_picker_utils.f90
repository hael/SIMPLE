module simple_picker_utils
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: picker_utils
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_TARGET1 = 4.0, SMPD_TARGET2 = 2.0
integer, parameter :: OFFSET = 3

type picker_utils
    private
    integer                  :: ldim_raw(3) = 0 , ldim_shrink1(3) = 0 , ldim_shrink2(3) = 0
    integer                  :: ldim_box(3) = 0 , ldim_box1(3)    = 0 , ldim_box2(3)    = 0
    real                     :: smpd_raw    = 0., smpd_shrink1    = 0., smpd_shrink2    = 0.
    type(image), pointer     :: mic_raw => null()
    type(image)              :: mic_shrink1, mic_shrink2
    integer                  :: nboxes1 = 0, nboxes2 = 0
    type(image), allocatable :: boximgs1(:), boximgs2(:)
    integer,     allocatable :: positions1(:,:)
  contains
    procedure :: set_mics
    procedure :: set_positions
end type picker_utils

contains

    subroutine set_mics( self, mic, smpd, lp )
        class(picker_utils),  intent(inout) :: self
        class(image), target, intent(in)    :: mic
        real,                 intent(in)    :: smpd    !< sampling distance in A
        real, optional,       intent(in)    :: lp      !< low-pass limit in A
        real :: scale1, scale2
        ! set raw micrograph info
        self%ldim_raw = mic%get_ldim()
        if( self%ldim_raw(3) /= 1 ) THROW_HARD('Only for 2D images')
        self%smpd_raw = smpd
        self%mic_raw => mic
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
        if( present(lp) )then
            call self%mic_shrink1%bp(0.,lp)
            call self%mic_shrink2%bp(0.,lp)
        endif
        call self%mic_raw%ifft
        call self%mic_shrink1%ifft
        call self%mic_shrink2%ifft
    end subroutine set_mics

    subroutine set_positions( self, maxdiam )
        class(picker_utils),  intent(inout) :: self
        real,                 intent(in)    :: maxdiam !< maximum diameter in A
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
    end subroutine set_positions

end module simple_picker_utils
