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
    real                          :: total_dose, doseperframe, preexposure, scale, smpd, kv
    integer                       :: ldim(3), ldim_sc(2)
    integer                  :: boxsz, nptcls, nframes, ref_frame, start_frame
    integer                  :: eer_fraction, eer_upsampling
    logical                  :: l_doseweighing = .false.
    logical                  :: l_scale        = .false.
    logical                  :: l_gain         = .false.
    logical                  :: l_eer          = .false.
    logical                  :: exists         = .false.
  contains
    procedure          :: init
    procedure          :: display
    procedure          :: prep_frames
    procedure          :: extract_ptcl
    procedure, private :: pix2polycoords
    procedure, private :: get_local_shift
    ! Destructor
    procedure :: kill
end type ptcl_extractor

contains

    !>  Constructor
    subroutine init( self, docname )
        use simple_oris, only: oris
        class(ptcl_extractor), intent(inout) :: self
        character(len=*),      intent(in)    :: docname
        type(str4arr),    allocatable :: names(:)
        type(starfile_table_type) :: table
        integer(C_long)  :: num_objs, object_id
        integer          :: i,iframe,n,nmics,ind, motion_model
        logical          :: err
        call self%kill
        self%docname = trim(docname)
        call starfile_table__new(table)
        call starfile_table__getnames(table, trim(self%docname)//C_NULL_CHAR, names)
        n = size(names)
        do i = 1,n
            call starfile_table__read(table, trim(self%docname)//C_NULL_CHAR, names(i)%str )
            select case(trim(names(i)%str))
            case('general')
                self%ldim(1)        = parse_int(table, EMDL_IMAGE_SIZE_X, err)
                self%ldim(2)        = parse_int(table, EMDL_IMAGE_SIZE_Y, err)
                self%nframes        = parse_int(table, EMDL_IMAGE_SIZE_Z, err)
                call parse_string(table, EMDL_MICROGRAPH_MOVIE_NAME, self%moviename, err)
                self%l_eer          = fname2format(self%moviename) == 'K'
                call parse_string(table, EMDL_MICROGRAPH_GAIN_NAME, self%gainrefname, err)
                self%l_gain         = .not.err
                self%scale          = parse_double(table, EMDL_MICROGRAPH_BINNING, err)
                self%scale          = 1./self%scale
                self%l_scale        = (.not.err) .and. (abs(self%scale - 1.0) > 0.01)
                self%smpd           = parse_double(table, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, err)
                self%doseperframe   = parse_double(table, EMDL_MICROGRAPH_DOSE_RATE, err)
                self%l_doseweighing = (.not.err) .and. (self%doseperframe > 0.0001)
                self%preexposure    = parse_double(table, EMDL_MICROGRAPH_PRE_EXPOSURE, err)
                self%kv             = parse_double(table, EMDL_CTF_VOLTAGE, err)
                self%start_frame    = parse_int(table, EMDL_MICROGRAPH_START_FRAME, err)
                if( self%l_eer )then
                    self%eer_upsampling =  parse_int(table, EMDL_MICROGRAPH_EER_UPSAMPLING, err)
                    self%eer_fraction   = parse_int(table, EMDL_MICROGRAPH_EER_GROUPING, err)
                endif
                motion_model = parse_int(table, EMDL_MICROGRAPH_MOTION_MODEL_VERSION, err)
            case('global_shift')
                object_id  = starfile_table__firstobject(table)
                num_objs   = starfile_table__numberofobjects(table)
                if( int(num_objs - object_id) /= self%nframes ) THROW_HARD('Inconsistent # of shift entries and frames')
                allocate(self%isoshifts(self%nframes,2),source=0.)
                iframe = 0
                do while( (object_id < num_objs) .and. (object_id >= 0) )
                    iframe = iframe + 1
                    self%isoshifts(iframe,1) = parse_double(table, EMDL_MICROGRAPH_SHIFT_X, err)
                    self%isoshifts(iframe,2) = parse_double(table, EMDL_MICROGRAPH_SHIFT_Y, err)
                    object_id = starfile_table__nextobject(table)
                end do
                self%isoshifts = self%isoshifts / self%scale
            case('local_motion_model')
                ! todo
            case('hot_pixels')
                ! todo
            case DEFAULT
                THROW_HARD('Invalid table: '//trim(names(i)%str))
            end select
        enddo
        call starfile_table__delete(table)
        self%total_dose = real(self%nframes) * self%doseperframe
    end subroutine init

    subroutine display( self )
        class(ptcl_extractor), intent(in) :: self
        integer :: i
        print *, 'docname        ', self%docname
        print *, 'nframes        ', self%nframes
        print *, 'dimensions     ', self%ldim
        print *, 'smpd           ', self%smpd
        print *, 'voltage        ', self%kv
        print *, 'doseperframe   ', self%doseperframe
        print *, 'gainrefname    ', trim(self%gainrefname)
        print *, 'moviename      ', trim(self%moviename)
        print *, 'doseweighting  ', self%l_doseweighing
        print *, 'total dose     ', self%total_dose
        print *, 'scale          ', self%l_scale
        print *, 'gain           ', self%l_gain
        print *, 'eer            ', self%l_eer
        print *, 'eer_fraction   ', self%eer_fraction
        print *, 'eer_upsampling ', self%eer_upsampling
        if( allocated(self%isoshifts) )then
            do i = 1,size(self%isoshifts,dim=1)
                print *,'isoshifts ',i,self%isoshifts(i,:)
            enddo
        endif
    end subroutine display

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

    integer function parse_int( table, emdl_id, err )
        class(starfile_table_type) :: table
        integer, intent(in)        :: emdl_id
        logical, intent(out)       :: err
        err = starfile_table__getValue_int(table, emdl_id, parse_int)
        err = .not.err
    end function parse_int

    real function parse_double( table, emdl_id, err )
        class(starfile_table_type) :: table
        integer, intent(in)        :: emdl_id
        logical, intent(out)       :: err
        real(dp) :: v
        err = starfile_table__getValue_double(table, emdl_id, v)
        err = .not.err
        parse_double = real(v)
    end function parse_double

    subroutine parse_string( table, emdl_id, string, err )
        class(starfile_table_type)                 :: table
        integer,                       intent(in)  :: emdl_id
        character(len=:), allocatable, intent(out) :: string
        logical, intent(out)       :: err
        err = starfile_table__getValue_string(table, emdl_id, string)
        err = .not.err
    end subroutine parse_string

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
        self%l_doseweighing = .false.
        self%l_gain         = .false.
        self%l_eer          = .false.
        self%nframes        = 0
        self%doseperframe   = 0.
        self%scale          = 1.
        self%moviename   = ''
        self%stkname     = ''
        self%docname     = ''
        self%gainrefname = ''
        self%exists = .false.
    end subroutine kill

end module simple_particle_extractor