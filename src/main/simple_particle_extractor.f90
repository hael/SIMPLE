module simple_particle_extractor
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,                         only: image, image_ptr
use simple_parameters,                    only: params_glob
use simple_starfile_wrappers
implicit none
! private
#include "simple_local_flags.inc"

public :: ptcl_extractor

integer, parameter :: POLYDIM = 18
real,    parameter :: NSIGMAS = 6.
logical, parameter :: DEBUG   = .true.

type :: ptcl_extractor
    type(image),      allocatable :: frames(:)
    real,             allocatable :: weights(:), isoshifts(:,:)
    character(len=:), allocatable :: gainrefname, stkname, moviename, docname
    real(dp)                      :: polyx(POLYDIM), polyy(POLYDIM)
    real                          :: doseperframe, scale, smpd, kv
    integer                       :: ldim(3), ldim_sc(2)
    integer                  :: boxsz, nptcls, nframes, ref_frame
    logical                  :: l_doseweighting = .false.
    logical                  :: l_scale = .false.
    logical                  :: l_gain  = .false.
    logical                  :: l_eer   = .false.
    logical                  :: l_error = .false.
    logical                  :: exists
  contains
    procedure          :: read_doc
    procedure          :: prep_frames
    procedure          :: extract_ptcl
    procedure, private :: pix2polycoords
    procedure, private :: get_local_shift
    ! Destructor
    procedure :: kill
end type ptcl_extractor

contains

    subroutine read_doc( self, docname )
        use simple_oris, only: oris
        class(ptcl_extractor), intent(inout) :: self
        character(len=*),      intent(in)    :: docname
        call self%kill
    end subroutine read_doc

    subroutine extract_ptcl( self, ptcl_pos_in, box_in, ptcl )
        class(ptcl_extractor), intent(inout) :: self
        integer,                 intent(in)    :: ptcl_pos_in(2), box_in
        type(image),             intent(inout) :: ptcl
        type(image) :: frame_ptcl(self%nframes)
        real(dp)    :: x,y
        real        :: subpixel_shift(2),total_shift(2), aniso_shift(2), center(2)
        real        :: scale
        integer     :: pos(2), discrete_shift(2), foo, box, t
        ! init dimensions
        box = box_in
        pos = ptcl_pos_in
        if( self%l_scale )then
          box = round2even(real(box_in)/self%scale)
          pos = nint(real(ptcl_pos_in)/self%scale)
        endif
        center = real(pos+box/2+1)
        x      = real(center(1),dp)
        y      = real(center(2),dp)
        ! init images
        call ptcl%new([box,box,1], self%smpd)
        call ptcl%zero_and_flag_ft
        ! work...
        !$omp parallel do private(t,aniso_shift,total_shift,discrete_shift,subpixel_shift,foo)&
        !$omp default(shared) proc_bind(close) schedule(static)
        do t = 1,self%nframes
            call frame_ptcl(t)%new([box,box,1], self%smpd, wthreads=.false.)
            ! manage shifts
            call self%get_local_shift(t, x, y, aniso_shift)
            total_shift    = self%isoshifts(t,:) + aniso_shift
            discrete_shift = nint(total_shift)
            subpixel_shift = total_shift - real(discrete_shift)
            ! extract particle frame
            call self%frames(t)%window(pos-discrete_shift, box, frame_ptcl(t), foo)
            ! shift for subpixel accuracy
            call frame_ptcl(t)%fft
            call frame_ptcl(t)%shift2Dserial(-subpixel_shift)
            ! weight
            call frame_ptcl(t)%mul(self%weights(t))
        enddo
        !$omp end parallel do
        ! sum
        do t = 1,self%nframes
            call ptcl%add_workshare(frame_ptcl(t))
            call frame_ptcl(t)%kill
        enddo
        ! return to real space, account for scale
        if( self%l_scale ) call ptcl%clip_inplace([box_in,box_in,1])
        call ptcl%ifft
    end subroutine extract_ptcl

    !>  read movie, corrects for gain, dose & outliers
    subroutine prep_frames( self )
        class(ptcl_extractor), intent(inout) :: self

        
    end subroutine prep_frames

    !>  pixels to coordinates for polynomial evaluation (scaled in/out)
    elemental subroutine pix2polycoords( self, xin, yin, x, y )
        class(ptcl_extractor), intent(in) :: self
        real(dp),                intent(in)  :: xin, yin
        real(dp),                intent(out) :: x, y
        x = (xin-1.d0) / real(self%ldim(1)-1,dp) - 0.5d0
        y = (yin-1.d0) / real(self%ldim(2)-1,dp) - 0.5d0
    end subroutine pix2polycoords

    pure subroutine get_local_shift( self, iframe, x, y, shift )
        class(ptcl_extractor), intent(in)  :: self
        integer,                 intent(in)  :: iframe
        real(dp),                intent(in)  :: x, y
        real,                    intent(out) :: shift(2)
        real(dp) :: t, xx, yy
        t = real(iframe-self%ref_frame, dp)
        call self%pix2polycoords(x,y, xx,yy)
        shift(1) = polyfun(self%polyx(:), xx,yy,t)
        shift(2) = polyfun(self%polyy(:), xx,yy,t)
    end subroutine get_local_shift

    pure real function polyfun(c, x, y, t)
        real(dp), intent(in) :: c(POLYDIM), x, y, t
        real(dp) :: res, t2, t3
        t2 = t * t
        t3 = t2 * t
        res =       dot_product( c(1:3),         [t,t2,t3])
        res = res + dot_product( c(4:6),     x * [t,t2,t3]) + dot_product( c(7:9),   x*x * [t,t2,t3])
        res = res + dot_product( c(10:12),   y * [t,t2,t3]) + dot_product( c(13:15), y*y * [t,t2,t3])
        res = res + dot_product( c(16:18), x*y * [t,t2,t3])
        polyfun = real(res)
    end function polyfun

    subroutine kill(self)
        class(ptcl_extractor), intent(inout) :: self
        integer :: i
        if( allocated(self%weights) )   deallocate(self%weights)
        if( allocated(self%isoshifts) ) deallocate(self%isoshifts)
        if( allocated(self%frames) )then
            do i=1,self%nframes
                call self%frames(i)%kill
            enddo
            deallocate(self%frames)
        endif
        self%l_doseweighting = .false.
        self%l_gain          = .false.
        self%l_error         = .false.
        self%l_eer           = .false.
        self%nframes         = 0
        self%doseperframe    = 0.
        self%scale           = 1.
        self%moviename   = ''
        self%stkname     = ''
        self%docname     = ''
        self%gainrefname = ''
        self%exists = .false.
    end subroutine kill

end module simple_particle_extractor