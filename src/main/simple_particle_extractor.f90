module simple_particle_extractor
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,                         only: image, image_ptr
use simple_eer_factory,                   only: eer_decoder
use simple_motion_correct,                only: correct_gain
use simple_starfile_wrappers
implicit none
! private
#include "simple_local_flags.inc"

public :: ptcl_extractor

integer, parameter :: POLYDIM = 18
logical, parameter :: DEBUG   = .true.

type :: ptcl_extractor
    type(image),      allocatable :: frames(:)
    type(image),      allocatable :: particle(:), frame_particle(:)
    type(eer_decoder)             :: eer
    real,             allocatable :: doses(:,:,:),weights(:), isoshifts(:,:)
    integer,          allocatable :: hotpix_coords(:,:)
    logical,          allocatable :: particle_mask(:,:,:)
    character(len=:), allocatable :: gainrefname, moviename, docname
    real(dp)                      :: polyx(POLYDIM), polyy(POLYDIM)
    real                          :: smpd, smpd_out
    real                          :: total_dose, doseperframe, preexposure, scale, kv
    integer                       :: ldim(3), box, box_out
    integer                       :: nframes, ref_frame, start_frame, nhotpix
    integer                       :: eer_fraction, eer_upsampling
    logical                       :: l_doseweighing = .false.
    logical                       :: l_scale        = .false.
    logical                       :: l_gain         = .false.
    logical                       :: l_eer          = .false.
    logical                       :: l_poly         = .false.
    logical                       :: l_neg          = .true.
    logical                       :: exists         = .false.
  contains
    procedure          :: init
    procedure, private :: generate_dose_weighing
    procedure          :: display
    procedure          :: extract_ptcl
    procedure, private :: pix2polycoords
    procedure, private :: get_local_shift
    ! Destructor
    procedure :: kill
end type ptcl_extractor

contains

    !>  Constructor
    subroutine init( self, docname, box, neg )
        class(ptcl_extractor), intent(inout) :: self
        character(len=*),      intent(in)    :: docname
        integer,               intent(in)    :: box
        logical,               intent(in)    :: neg
        type(str4arr), allocatable :: names(:)
        type(image)                :: gain, tmp
        type(starfile_table_type)  :: table
        real            :: radius
        integer(C_long) :: num_objs, object_id
        integer         :: i,j,iframe,n,nmics,ind, motion_model
        logical         :: err
        call self%kill
        self%l_neg   = neg
        self%docname = trim(docname)
        ! parsing individual movie meta-data
        call starfile_table__new(table)
        call starfile_table__getnames(table, trim(self%docname)//C_NULL_CHAR, names)
        n = size(names)
        do i = 1,n
            call starfile_table__read(table, trim(self%docname)//C_NULL_CHAR, names(i)%str )
            select case(trim(names(i)%str))
            case('general')
                ! global variables, movie at original size
                self%ldim(1)        = parse_int(table, EMDL_IMAGE_SIZE_X, err)
                self%ldim(2)        = parse_int(table, EMDL_IMAGE_SIZE_Y, err)
                self%ldim(3)        = 0
                self%nframes        = parse_int(table, EMDL_IMAGE_SIZE_Z, err)
                call parse_string(table, EMDL_MICROGRAPH_MOVIE_NAME, self%moviename, err)
                self%l_eer          = fname2format(self%moviename) == 'K'
                call parse_string(table, EMDL_MICROGRAPH_GAIN_NAME, self%gainrefname, err)
                self%l_gain         = .not.err
                self%scale          = 1./ parse_double(table, EMDL_MICROGRAPH_BINNING, err)
                self%l_scale        = (.not.err) .and. (abs(self%scale - 1.0) > 0.001)
                self%smpd           = parse_double(table, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, err)
                self%doseperframe   = parse_double(table, EMDL_MICROGRAPH_DOSE_RATE, err)
                self%l_doseweighing = (.not.err) .and. (self%doseperframe > 0.0001)
                self%preexposure    = parse_double(table, EMDL_MICROGRAPH_PRE_EXPOSURE, err)
                self%kv             = parse_double(table, EMDL_CTF_VOLTAGE, err)
                self%start_frame    = parse_int(table, EMDL_MICROGRAPH_START_FRAME, err)
                if( self%l_eer )then
                    self%eer_upsampling = parse_int(table, EMDL_MICROGRAPH_EER_UPSAMPLING, err)
                    self%eer_fraction   = parse_int(table, EMDL_MICROGRAPH_EER_GROUPING, err)
                endif
                motion_model = parse_int(table, EMDL_MICROGRAPH_MOTION_MODEL_VERSION, err)
                self%l_poly  = motion_model == 1
            case('global_shift')
                ! parse isotropic shifts
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
            case('local_motion_model')
                ! parse polynomial coefficients
                object_id  = starfile_table__firstobject(table)
                num_objs   = starfile_table__numberofobjects(table)
                if( int(num_objs - object_id) /= 2*POLYDIM ) THROW_HARD('Inconsistent # polynomial coefficient')
                ind = 0
                do while( (object_id < num_objs) .and. (object_id >= 0) )
                    ind = ind+1
                    j = parse_int(table, EMDL_MICROGRAPH_MOTION_COEFFS_IDX, err)
                    if( j < POLYDIM)then
                        self%polyx(j+1) = parse_double(table, EMDL_MICROGRAPH_MOTION_COEFF, err)
                    else
                        self%polyy(j-POLYDIM+1) = parse_double(table, EMDL_MICROGRAPH_MOTION_COEFF, err)
                    endif
                    object_id = starfile_table__nextobject(table)
                end do
            case('hot_pixels')
                object_id  = starfile_table__firstobject(table)
                num_objs   = starfile_table__numberofobjects(table)
                self%nhotpix = int(num_objs)
                allocate(self%hotpix_coords(2,int(self%nhotpix)),source=-1)
                j = 0
                do while( (object_id < num_objs) .and. (object_id >= 0) )
                    j = j+1
                    self%hotpix_coords(1,j) = nint(parse_double(table, EMDL_IMAGE_COORD_X, err))
                    self%hotpix_coords(2,j) = nint(parse_double(table, EMDL_IMAGE_COORD_Y, err))
                    object_id = starfile_table__nextobject(table)
                end do
            case DEFAULT
                THROW_HARD('Invalid table: '//trim(names(i)%str))
            end select
        enddo
        call starfile_table__delete(table)
        ! other variables
        allocate(self%weights(self%nframes),source=1.0/real(self%nframes)) ! needs to be updated
        self%total_dose = real(self%nframes) * self%doseperframe
        ! updates dimensions and pixel size
        if( self%l_eer )then
            select case(self%eer_upsampling)
                case(1)
                    ! 4K
                case(2)
                    ! 8K
                    self%ldim(1:2) = 2 * self%ldim(1:2)
                    self%smpd      = self%smpd / 2.0
                case DEFAULT
                    THROW_HARD('Unsupported up-sampling: '//int2str(self%eer_upsampling))
            end select
        endif
        self%smpd_out = self%smpd / self%scale
        ! allocate & read frames
        allocate(self%frames(self%nframes))
        if( self%l_eer )then
            call self%eer%new(self%moviename, self%smpd, self%eer_upsampling)
            call self%eer%decode(self%frames, self%eer_fraction)
        else
            !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
            do iframe=1,self%nframes
                call self%frames(iframe)%new(self%ldim, self%smpd, wthreads=.false.)
            enddo
            !$omp end parallel do
            do iframe=1,self%nframes
                call self%frames(iframe)%read(self%moviename, iframe)
            end do
        endif
        ! gain correction
        if( self%l_gain )then
            if( .not.file_exists(self%gainrefname) )then
                THROW_HARD('gain reference: '//trim(self%gainrefname)//' not found')
            endif
            if( self%l_eer )then
                call correct_gain(self%frames, self%gainrefname, gain, eerdecoder=self%eer)
            else
                call correct_gain(self%frames, self%gainrefname, gain)
            endif
        endif
        ! outliers curation
        if( self%nhotpix > 0 )then
            ! TODO
        endif
        call gain%kill
        ! dimensions of the particle & frame particle
        self%box_out = box
        self%box     = box
        if( self%l_scale ) self%box = round2even(real(self%box)/self%scale)
        self%box = find_larger_magic_box(self%box+2) ! subpixel shift & fftw friendly
        allocate(self%frame_particle(nthr_glob),self%particle(nthr_glob))
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1,nthr_glob
            call self%particle(i)%new(      [self%box,self%box,1], self%smpd, wthreads=.false.)
            call self%frame_particle(i)%new([self%box,self%box,1], self%smpd, wthreads=.false.)
        enddo
        !$omp end parallel do
        ! dose weighting
        call self%generate_dose_weighing
        ! mask for post-extraction normalizations
        radius = RADFRAC_NORM_EXTRACT * real(self%box_out/2)
        call tmp%disc([self%box_out,self%box_out,1], 1., radius, self%particle_mask)
        call tmp%kill
        ! index to reference frame
        self%ref_frame = 1 ! this remains to be checked? 
        ! all done
        self%exists = .true.
    end subroutine init

    subroutine generate_dose_weighing( self )
        class(ptcl_extractor), intent(inout) :: self
        real,   parameter :: A=0.245, B=-1.665, C=2.81
        real, allocatable :: qs(:), accumulated_doses(:)
        real              :: spaFreqk, dose_per_frame, twoNe, spafreq, limsq
        integer           :: cshape(3), nrflims(3,2), ldim(3), hphys,kphys, iframe, h,k
        if( .not.self%l_doseweighing ) return
        cshape = self%frame_particle(1)%get_array_shape()
        allocate(self%doses(cshape(1),cshape(2),self%nframes),&
            &qs(self%nframes), accumulated_doses(self%nframes), source=0.)
        do iframe=1,self%nframes
            accumulated_doses(iframe) = real(iframe) * self%doseperframe
        end do
        if( is_equal(self%kV,200.) )then
            accumulated_doses = accumulated_doses / 0.8
        else if( is_equal(self%kV,100.) )then
            accumulated_doses = accumulated_doses / 0.64
        endif
        limsq = (real(self%box)*self%smpd)**2.
        !$omp parallel private(h,k,spafreq,spafreqk,twone,kphys,hphys,qs)&
        !$omp default(shared) proc_bind(close)
        !$omp do schedule(static)
        do k = nrflims(2,1),nrflims(2,2)
            kphys    = k + 1 + merge(ldim(2),0,k<0)
            spaFreqk = real(k*k)/limsq
            do h = nrflims(1,1),nrflims(1,2)
                hphys   = h + 1
                spaFreq = sqrt( real(h*h)/limsq + spaFreqk )
                twoNe   = 2.*(A*spaFreq**B + C)
                qs = exp(-accumulated_doses/twoNe)
                qs = qs / sqrt(sum(qs*qs))
                self%doses(hphys,kphys,:) = qs
            enddo
        enddo
        !$omp end do
        !$omp end parallel
    end subroutine generate_dose_weighing

    subroutine display( self )
        class(ptcl_extractor), intent(in) :: self
        integer :: i
        print *, 'docname        ', self%docname
        print *, 'nframes        ', self%nframes
        print *, 'dimensions     ', self%ldim
        print *, 'smpd           ', self%smpd
        print *, 'smpd_out       ', self%smpd_out
        print *, 'voltage        ', self%kv
        print *, 'doseperframe   ', self%doseperframe
        print *, 'gainrefname    ', trim(self%gainrefname)
        print *, 'moviename      ', trim(self%moviename)
        print *, 'doseweighting  ', self%l_doseweighing
        print *, 'total dose     ', self%total_dose
        print *, 'scale          ', self%l_scale
        print *, 'gain           ', self%l_gain
        print *, 'nhotpix        ', self%nhotpix
        print *, 'eer            ', self%l_eer
        print *, 'eer_fraction   ', self%eer_fraction
        print *, 'eer_upsampling ', self%eer_upsampling
        if( allocated(self%isoshifts) )then
            do i = 1,size(self%isoshifts,dim=1)
                print *,'isoshifts ',i,self%isoshifts(i,:)
            enddo
        endif
        do i = 1,POLYDIM
            print *,'polycoeffs    ',i,self%polyx(i),self%polyy(i)
        enddo
        ! if( self%nhotpix > 0 )then
        !     do i = 1,self%nhotpix
        !         print *,'hotpix    ',i,self%hotpix_coords(:,i)
        !     enddo
        ! endif
    end subroutine display

    !>  on a single thread
    subroutine extract_ptcl( self, ptcl_pos_in, ptcl_out )
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: ptcl_pos_in(2) ! top left corner
        class(image),          intent(inout) :: ptcl_out
        real(dp)    :: x,y
        real        :: subpixel_shift(2),total_shift(2), aniso_shift(2), center(2)
        real        :: scale, sdev_noise
        integer     :: pos(2), discrete_shift(2), foo, t, ithr
        ithr = omp_get_thread_num() + 1
        ! sanity check
        if( any(ptcl_out%get_ldim() /= [self%box_out,self%box_out,1]) )then
            THROW_HARD('Inconsistent dimensions!')
        endif
        ! particle coordinates
        pos = ptcl_pos_in
        if( self%l_scale ) pos = nint(real(ptcl_pos_in)/self%scale)
        center = real(pos + self%box/2 + 1)
        x      = real(center(1),dp)
        y      = real(center(2),dp)
        ! work...
        call self%particle(ithr)%zero_and_flag_ft
        do t = 1,self%nframes
            ! manage shifts
            if( self%l_poly )then
                call self%get_local_shift(t, x, y, aniso_shift)
                total_shift = self%isoshifts(t,:) + aniso_shift
            else
                total_shift = self%isoshifts(t,:)
            endif
            discrete_shift = nint(total_shift)
            subpixel_shift = total_shift - real(discrete_shift)
            ! extract particle frame
            call self%frames(t)%window(pos-discrete_shift, self%box, self%frame_particle(ithr), foo)
            ! subpixel shift
            call self%frame_particle(ithr)%fft
            call self%frame_particle(ithr)%shift2Dserial(-subpixel_shift)
            ! dose-weighing
            if( self%l_doseweighing )then
                call self%frame_particle(ithr)%mul_cmat(self%doses(:,:,t:t))
            endif
            ! frame weight
            call self%frame_particle(ithr)%mul(self%weights(t))
            ! sum
            call self%particle(ithr)%add(self%frame_particle(ithr))
        enddo
        ! clipping to correct size
        call self%particle(ithr)%ifft 
        call self%particle(ithr)%clip(ptcl_out)
        ! post-processing
        if( self%l_neg ) call ptcl_out%neg()
        call ptcl_out%subtr_backgr_ramp(self%particle_mask)
        call ptcl_out%norm_noise(self%particle_mask, sdev_noise)
    end subroutine extract_ptcl

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
        if( allocated(self%particle) )then
            do i = 1,nthr_glob
                call self%particle(i)%kill
                call self%frame_particle(i)%kill
            enddo
            deallocate(self%particle,self%frame_particle)
        endif
        if( allocated(self%doses) )         deallocate(self%doses)
        if( allocated(self%hotpix_coords) ) deallocate(self%hotpix_coords)
        if( allocated(self%particle_mask) ) deallocate(self%particle_mask)
        call self%eer%kill
        self%l_doseweighing = .false.
        self%l_gain         = .false.
        self%l_eer          = .false.
        self%l_neg          = .true.
        self%nframes        = 0
        self%nhotpix        = 0
        self%doseperframe   = 0.
        self%scale          = 1.
        self%moviename   = ''
        self%docname     = ''
        self%gainrefname = ''
        self%exists = .false.
    end subroutine kill

end module simple_particle_extractor