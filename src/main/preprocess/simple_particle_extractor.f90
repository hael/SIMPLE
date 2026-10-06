
!@descr: core functionality for extracting particles from micrographs
module simple_particle_extractor
use simple_core_module_api
use simple_image,                only: image
use simple_eer_factory,          only: eer_decoder
use simple_motion_correct_utils, only: correct_gain, pix2polycoords, apply_patch_poly
use simple_starfile_wrappers
use simple_motion_model,         only: motion_model
implicit none
private
#include "simple_local_flags.inc"

public :: ptcl_extractor
public :: dealloc_particles_frames

integer, parameter :: POLYDIM    = 18
real,    parameter :: NSIGMAS    = 6.0
logical, parameter :: DEBUG_HERE = .false.

type(image), allocatable :: particles_frames(:,:)

type :: ptcl_extractor
    type(image),      allocatable :: frames(:)
    type(image),      allocatable :: particle(:), frame_particle(:)
    type(motion_model)            :: mmodel
    real,             allocatable :: doses(:,:,:),weights(:), isoshifts(:,:)
    integer,          allocatable :: hotpix_coords(:,:)
    logical,          allocatable :: particle_mask(:,:,:)
    type(string)                  :: gainrefname, moviename, docname
    type(image)                   :: gain
    type(eer_decoder)             :: eer
    real(dp)                      :: polyx(POLYDIM), polyy(POLYDIM)
    real                          :: smpd, smpd_out
    real                          :: scale
    real                          :: total_dose, doseperframe, preexposure, kv
    integer                       :: ldim(3), ldim_sc(3), box, box_pd
    integer                       :: nframes, start_frame, nhotpix, align_frame
    integer                       :: eer_fraction, eer_upsampling=1
    logical                       :: l_doseweighing = .false.
    logical                       :: l_scale        = .false.
    logical                       :: l_gain         = .false.
    logical                       :: l_eer          = .false.
    logical                       :: l_poly         = .false.
    logical                       :: l_neg          = .true.
    logical                       :: l_mov          = .true.
    logical                       :: l_model        = .false.
    logical                       :: exists         = .false.
  contains
    procedure          :: init_model
    procedure          :: init_mov
    procedure          :: init_mic
    procedure, private :: init_mask
    procedure, private :: parse_movie_metadata
    procedure          :: display
    procedure          :: extract_particles
    procedure          :: extract_particles_from_mic
    procedure          :: extract_all_particles_frames
    procedure          :: generate_cumulative_particles_stack
    procedure, private :: extract_ptcl
    procedure, private :: extract_single_particle_frame
    procedure, private :: post_process
    procedure, private :: cure_outliers
    procedure, private :: add_eer_gain_defects
    procedure, private :: cure_outliers_from_coords, cure_outlier_coords
    procedure, private :: evaluate_local_shift
    procedure          :: get_nframes
    ! Destructor
    procedure :: kill
end type ptcl_extractor

contains

    subroutine init_model( self, omov, params )
        use simple_ori,        only: ori
        use simple_parameters, only: parameters
        class(ptcl_extractor),     intent(inout) :: self
        class(ori),                intent(in)    :: omov
        class(parameters), target, intent(in)    :: params
        integer :: i, iframe
        call self%kill
        if( .not. omov%isthere('mcmodel') ) THROW_HARD('Motion model entry is absent!')
        call omov%getter('mcmodel',self%docname)
        if( .not.file_exists(self%docname) ) THROW_HARD('Motion model doc is absent!')
        call self%mmodel%read(self%docname, params)
        if( .not.file_exists(self%mmodel%movie) ) THROW_HARD('Movie is absent!')
        if( self%mmodel%nframes < 1 ) THROW_HARD('Motion model contains no movie frames!')
        self%l_model = .true.
        self%l_mov   = .false.
        self%l_neg   = (params%pcontrast .eq. 'black')
        self%box     = params%box
        ! Movie and image geometry
        self%moviename = self%mmodel%movie
        self%smpd      = self%mmodel%smpd_movie
        self%smpd_out  = self%mmodel%smpd
        self%scale     = 1. / self%mmodel%binning
        self%l_scale   = abs(self%scale - 1.) > 0.001
        self%ldim      = [self%mmodel%ldim_movie, 1]
        self%ldim_sc   = [self%mmodel%ldim,       1]
        self%nframes   = self%mmodel%nframes
        self%start_frame = 1
        if( params%tof > self%nframes) THROW_HARD('TOF is large than the number of frames!')
        ! Dose and detector metadata
        self%kv              = self%mmodel%voltage
        self%doseperframe    = self%mmodel%dose_per_frame
        self%total_dose      = self%mmodel%accumulated_dose
        self%preexposure     = 0.
        self%l_doseweighing  = self%mmodel%dw
        self%l_eer           = self%mmodel%eer
        self%eer_fraction    = self%mmodel%eer_fraction
        self%eer_upsampling  = self%mmodel%eer_upsampling
        ! Global and local motion
        allocate(self%isoshifts(2,self%nframes), self%weights(self%nframes), source=0.)
        if( allocated(self%mmodel%drift_offsets_x) .and. allocated(self%mmodel%drift_offsets_y) )then
            if( size(self%mmodel%drift_offsets_x) /= self%nframes .or.&
                &size(self%mmodel%drift_offsets_y) /= self%nframes )then
                THROW_HARD('Motion model drift-offset dimensions do not match!')
            endif
            self%isoshifts(1,:) = self%mmodel%drift_offsets_x
            self%isoshifts(2,:) = self%mmodel%drift_offsets_y
        endif
        if( allocated(self%mmodel%frameweights) )then
            if( size(self%mmodel%frameweights) /= self%nframes )then
                THROW_HARD('Motion model frame-weight dimensions do not match!')
            endif
            self%weights = self%mmodel%frameweights
        else
            self%weights = 1. / real(self%nframes)
        endif
        self%l_poly = self%mmodel%npatch > 0 .and. allocated(self%mmodel%local_offsets_x) .and.&
            &allocated(self%mmodel%local_offsets_y) .and.self%mmodel%patch_accepted
        self%align_frame = self%mmodel%fixed_frame
        self%polyx = 0.d0
        self%polyy = 0.d0
        if( self%l_poly )then
            if( self%align_frame < 1 .or. self%align_frame > self%nframes )then
                THROW_HARD('Motion model reference frame is out of range!')
            endif
            call self%mmodel%refit_polynomial( self%align_frame )
            self%polyx = self%mmodel%model_coeffs_x
            self%polyy = self%mmodel%model_coeffs_y
        endif
        write(logfhandle,'(A,A)')'>> PARSED MODEL: ',self%docname%to_char()
        ! Standardize offsets. No scaling need be applied as stage drift and local
        ! offsets are determined and written as scaled, unlike in the star format.
        self%isoshifts(1,:) = self%isoshifts(1,:) - self%isoshifts(1,self%align_frame)
        self%isoshifts(2,:) = self%isoshifts(2,:) - self%isoshifts(2,self%align_frame)
        ! Gain reference and known detector outliers
        self%gainrefname = self%mmodel%gain
        self%l_gain      = .not.self%gainrefname%is_blank()
        if( allocated(self%mmodel%outlier_coords) )then
            self%nhotpix = size(self%mmodel%outlier_coords, dim=2)
            self%hotpix_coords = self%mmodel%outlier_coords
        endif
        ! Allocate and read frames
        allocate(self%frames(self%nframes))
        if( self%l_eer )then
            call self%eer%new(self%moviename, self%smpd, self%eer_upsampling)
            call self%eer%decode(self%frames, self%eer_fraction)
        else
            !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
            do iframe = 1,self%nframes
                call self%frames(iframe)%new(self%ldim, self%smpd, wthreads=.false.)
            enddo
            !$omp end parallel do
            do iframe = 1,self%nframes
                call self%frames(iframe)%read(self%moviename, iframe)
            enddo
        endif
        write(logfhandle,'(A,A)')'>> PARSED MOVIE: ',self%moviename%to_char()
        ! Gain correction and outlier curation
        if( self%l_gain )then
            if( .not.file_exists(self%gainrefname) )then
                THROW_HARD('gain reference: '//self%gainrefname%to_char()//' not found')
            endif
            if( self%l_eer )then
                call correct_gain(self%frames, self%gainrefname, self%gain,&
                    &eerdecoder=self%eer, flipgain=self%mmodel%flipgain)
                call self%add_eer_gain_defects
            else
                call correct_gain(self%frames, self%gainrefname, self%gain, flipgain=self%mmodel%flipgain)
            endif
        endif
        if( self%nhotpix > 0 )then
            call self%cure_outliers_from_coords(self%hotpix_coords)
        endif
        call self%gain%kill
        write(logfhandle,'(A,A)')'>> CORRECTED FOR GAIN AND OUTLIERS: ',self%moviename%to_char()
        ! Downscale frames
        if( any(self%ldim /= self%ldim_sc ) )then
            !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
            do iframe = 1,self%nframes
                call self%frames(iframe)%fft
                call self%frames(iframe)%clip_inplace(self%ldim_sc)
                call self%frames(iframe)%ifft
            enddo
            !$omp end parallel do
        endif
        ! dose weighing, dev only
        ! call self%frames(1)%apply_dose_weighing(self%nframes, self%frames,&
        !         & [1,self%nframes], self%total_dose, self%kv)
        ! Per-thread particle buffers and normalization mask
        self%box_pd = find_larger_magic_box(self%box+1)
        allocate(self%frame_particle(nthr_glob), self%particle(nthr_glob))
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1,nthr_glob
            call self%particle(i)%new(      [self%box_pd,self%box_pd,1], self%smpd_out, wthreads=.false.)
            call self%frame_particle(i)%new([self%box_pd,self%box_pd,1], self%smpd_out, wthreads=.false.)
        enddo
        !$omp end parallel do
        call self%init_mask
        ! all done
        call self%eer%kill
        self%exists = .true.
    end subroutine init_model

    !>  Constructor
    subroutine init_mov( self, omic, box, neg )
        use simple_ori, only: ori
        class(ptcl_extractor), intent(inout) :: self
        class(ori),            intent(in)    :: omic
        integer,               intent(in)    :: box
        logical,               intent(in)    :: neg
        real(dp), allocatable :: poly(:)
        type(string) :: poly_fname
        integer      :: i,iframe
        call self%kill
        if( .not. omic%isthere('mc_starfile') )then
            THROW_HARD('Movie star doc is absent 1, reverting to micrograph extraction')
            self%l_mov = .false.
        endif
        self%l_mov   = .true.
        self%docname = omic%get('mc_starfile')
        self%l_neg   = neg
        self%box     = box
        if( .not.file_exists(self%docname) )then
            THROW_HARD('Movie star doc is absent 2, reverting to micrograph extraction')
            ! revert to mic extraction
            self%l_mov = .false.
        endif
        if( self%l_mov )then
            ! get movie info
            call self%parse_movie_metadata
            if( .not.file_exists(self%moviename) )then
                THROW_HARD('Movie is absent, reverting to micrograph extraction')
                ! revert to mic extraction
                self%l_mov = .false.
            else
                ! frame of reference
                self%isoshifts(1,:) = self%isoshifts(1,:) - self%isoshifts(1,self%align_frame)
                self%isoshifts(2,:) = self%isoshifts(2,:) - self%isoshifts(2,self%align_frame)
                ! downscaling shifts
                self%isoshifts = self%isoshifts * self%scale
                ! polynomial coefficients
                poly_fname = fname_new_ext(self%docname,string('poly'))
                if( file_exists(poly_fname) )then
                    poly = file2drarr(poly_fname)
                    self%polyx = poly(:POLYDIM)
                    self%polyy = poly(POLYDIM+1:)
                endif
                self%polyx = self%polyx * real(self%scale,dp)
                self%polyy = self%polyy * real(self%scale,dp)
                ! dose-weighing
                self%total_dose = real(self%nframes) * self%doseperframe
                ! updates dimensions and pixel size
                if( self%l_eer )then
                    select case(self%eer_upsampling)
                        case(1)
                            ! 4K
                        case(2)
                            ! 8K: no updates to dimensions and pixel size are required
                            ! unlike in motion correction to accomodate relion convention
                        case DEFAULT
                            THROW_HARD('Unsupported up-sampling: '//int2str(self%eer_upsampling))
                    end select
                endif
                self%smpd_out   = self%smpd / self%scale
                self%ldim_sc    = round2even(real(self%ldim)*self%scale)
                self%ldim_sc(3) = 1
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
                        THROW_HARD('gain reference: '//self%gainrefname%to_char()//' not found')
                    endif
                    if( self%l_eer )then
                        call correct_gain(self%frames, self%gainrefname, self%gain, eerdecoder=self%eer)
                    else
                        call correct_gain(self%frames, self%gainrefname, self%gain)
                    endif
                endif
                ! outliers curation
                if( self%nhotpix > 0 ) call self%cure_outliers
                call self%gain%kill
                ! downscale frames & dose-weighing
                !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
                do iframe=1,self%nframes
                    call self%frames(iframe)%fft
                    call self%frames(iframe)%clip_inplace(self%ldim_sc)
                enddo
                !$omp end parallel do
                call self%frames(1)%apply_dose_weighing(self%nframes, self%frames,&
                    &[1,self%nframes], self%total_dose, self%kv)
                !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
                do iframe=1,self%nframes
                    call self%frames(iframe)%ifft
                enddo
                !$omp end parallel do
                ! dimensions of the particle & frame particle
                self%box_pd = find_larger_magic_box(self%box+1) ! subpixel shift & fftw friendly
                allocate(self%frame_particle(nthr_glob),self%particle(nthr_glob))
                !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
                do i = 1,nthr_glob
                    call self%particle(i)%new(      [self%box_pd,self%box_pd,1], self%smpd_out, wthreads=.false.)
                    call self%frame_particle(i)%new([self%box_pd,self%box_pd,1], self%smpd_out, wthreads=.false.)
                enddo
                !$omp end parallel do
                ! mask for post-extraction normalizations
                call self%init_mask
            endif
        endif
        ! micrograph init
        if( .not.self%l_mov ) call self%init_mic( self%box, self%l_neg)
        ! all done
        call self%eer%kill
        self%exists = .true.
    end subroutine init_mov

    subroutine init_mic( self, box, neg )
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: box
        logical,               intent(in)    :: neg
        self%box   = box
        self%l_neg = neg
        call self%init_mask
    end subroutine init_mic

    ! mask for post-extraction normalizations
    subroutine init_mask( self )
        class(ptcl_extractor), intent(inout) :: self
        type(image) :: tmp
        real        :: radius
        if( allocated(self%particle_mask) ) deallocate(self%particle_mask)
        radius = RADFRAC_NORM_EXTRACT * real(self%box/2)
        call tmp%disc([self%box,self%box,1], 1., radius, self%particle_mask)
        call tmp%kill
    end subroutine init_mask

    subroutine parse_movie_metadata( self )
        class(ptcl_extractor), intent(inout) :: self
        type(string), allocatable     :: names(:)
        type(starfile_table_type)     :: table
        character(len=:), allocatable :: buffer
        integer(C_long) :: num_objs, object_id
        integer         :: i,j,iframe,n,ind, motion_model
        logical         :: err
        ! parsing individual movie meta-data
        call starfile_table__new(table)
        call starfile_table__getnames(table, self%docname, names)
        n = size(names)
        do i = 1,n
            call starfile_table__read(table, self%docname, names(i)%to_char() )
            select case(trim(names(i)%to_char()))
            case('general')
                ! global variables, movie at original size
                self%ldim(1)        = parse_int(table, EMDL_IMAGE_SIZE_X, err)
                self%ldim(2)        = parse_int(table, EMDL_IMAGE_SIZE_Y, err)
                self%ldim(3)        = 1
                self%nframes        = parse_int(table, EMDL_IMAGE_SIZE_Z, err)
                call parse_string(table, EMDL_MICROGRAPH_MOVIE_NAME, buffer, err)
                self%moviename      = buffer
                self%l_eer          = fname2format(self%moviename) == 'K'
                call parse_string(table, EMDL_MICROGRAPH_GAIN_NAME, buffer, err)
                self%gainrefname    = buffer
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
                motion_model     = parse_int(table, EMDL_MICROGRAPH_MOTION_MODEL_VERSION, err)
                self%l_poly      = motion_model == 1
                self%align_frame = 1
                if( self%l_poly ) self%align_frame = parse_int(table, SMPL_MOVIE_FRAME_ALIGN, err)
            case('global_shift')
                ! parse isotropic shifts
                object_id  = starfile_table__firstobject(table)
                num_objs   = starfile_table__numberofobjects(table)
                if( int(num_objs - object_id) /= self%nframes ) THROW_HARD('Inconsistent # of shift entries and frames')
                allocate(self%isoshifts(2,self%nframes),self%weights(self%nframes),source=0.)
                iframe = 0
                do while( (object_id < num_objs) .and. (object_id >= 0) )
                    iframe = iframe + 1
                    self%isoshifts(1,iframe) = parse_double(table, EMDL_MICROGRAPH_SHIFT_X, err)
                    self%isoshifts(2,iframe) = parse_double(table, EMDL_MICROGRAPH_SHIFT_Y, err)
                    self%weights(iframe)     = parse_double(table, SMPL_MOVIE_FRAME_WEIGHT, err)
                    if( err ) self%weights(iframe) = 1./real(self%nframes)
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
                THROW_HARD('Invalid table: '//trim(names(i)%to_char()))
            end select
        enddo
        call starfile_table__delete(table)
        if( DEBUG_HERE ) print *,'movie doc parsed'
    end subroutine parse_movie_metadata

    subroutine display( self )
        class(ptcl_extractor), intent(in) :: self
        integer :: i
        print *, 'docname        ', self%docname%to_char()
        print *, 'nframes        ', self%nframes
        print *, 'dimensions     ', self%ldim
        print *, 'smpd           ', self%smpd
        print *, 'smpd_out       ', self%smpd_out
        print *, 'box            ', self%box
        print *, 'box_pd         ', self%box_pd
        print *, 'voltage        ', self%kv
        print *, 'doseperframe   ', self%doseperframe
        print *, 'gainrefname    ', self%gainrefname%to_char()
        if( self%l_model ) print *, 'flipgain       ', trim(self%mmodel%flipgain)
        print *, 'moviename      ', self%moviename%to_char()
        print *, 'doseweighting  ', self%l_doseweighing
        print *, 'total dose     ', self%total_dose
        print *, 'l_scale        ', self%l_scale
        print *, 'scale          ', self%scale
        print *, 'gain           ', self%l_gain
        print *, 'nhotpix        ', self%nhotpix
        print *, 'eer            ', self%l_eer
        print *, 'eer_fraction   ', self%eer_fraction
        print *, 'eer_upsampling ', self%eer_upsampling
        print *, 'align_frame    ', self%align_frame
        if( allocated(self%isoshifts) )then
            do i = 1,size(self%isoshifts,dim=2)
                print *,'isoshifts ',i,self%isoshifts(:,i),self%weights(i)
            enddo
        endif
        do i = 1,POLYDIM
            print *,'polycoeffs    ',i,self%polyx(i),self%polyy(i)
        enddo
    end subroutine display

    subroutine extract_particles( self, pinds, coords, particles, vmin, vmax, vmean, vsdev )
        class(ptcl_extractor),    intent(inout) :: self
        integer,     allocatable, intent(in)    :: pinds(:)
        integer,                  intent(in)    :: coords(:,:)
        type(image), allocatable, intent(inout) :: particles(:)
        real,                     intent(out)   :: vmin, vmax, vmean, vsdev
        real    :: rmin, rmax, rmean, rsdev
        integer :: i,j,n,cnt
        logical :: l_err
        n = size(pinds)
        if( size(coords,dim=2) < n ) THROW_HARD('Inconsistent dimensions 1')
        if( size(particles)    < n ) THROW_HARD('Inconsistent dimensions 2')
        cnt    = 0
        vmin   = huge(vmin)
        vmax   = -vmin
        vmean  = 0.0
        vsdev  = 0.0
        !$omp parallel do schedule(static) default(shared) private(i,j,l_err,rmean,rsdev,rmax,rmin)&
        !$omp reduction(+:vmean,vsdev,cnt) reduction(min:vmin) reduction(max:vmax) proc_bind(close)
        do i = 1,n
            j = pinds(i)
            call self%extract_ptcl(coords(:,j),particles(i))
            call self%post_process(particles(i))
            call particles(i)%stats(rmean, rsdev, rmax, rmin, errout=l_err)
            if( .not.l_err )then
                cnt = cnt + 1
                vmin   = min(vmin,rmin)
                vmax   = max(vmax,rmax)
                vmean  = vmean + rmean
                vsdev  = vsdev + rsdev*rsdev
            endif
        end do
        !$omp end parallel do
        if( cnt > 0 )then
            vmean = vmean / real(cnt)
            vsdev = sqrt(vsdev / real(cnt))
        endif
    end subroutine extract_particles

    subroutine extract_particles_from_mic( self, mic, pinds, coords, particles, vmin, vmax, vmean, vsdev )
        class(ptcl_extractor),    intent(inout) :: self
        class(image),             intent(in)    :: mic
        integer,     allocatable, intent(in)    :: pinds(:)
        integer,                  intent(in)    :: coords(:,:)
        type(image),              intent(inout) :: particles(:)
        real,                     intent(out)   :: vmin, vmax, vmean, vsdev
        real    :: rmin, rmax, rmean, rsdev, sdev_noise
        integer :: i,j,n,cnt,noutside
        logical :: l_err
        n = size(pinds)
        if( size(coords,dim=2) < n ) THROW_HARD('Inconsistent dimensions 1')
        if( size(particles)    < n ) THROW_HARD('Inconsistent dimensions 2')
        cnt   = 0
        vmin  = huge(vmin)
        vmax  = -vmin
        vmean = 0.0
        vsdev = 0.0
        !$omp parallel do schedule(static) default(shared) proc_bind(close)&
        !$omp private(i,j,noutside,sdev_noise,l_err,rmin,rmax,rmean,rsdev)&
        !$omp reduction(+:vmean,vsdev,cnt) reduction(min:vmin) reduction(max:vmax)
        do i = 1,n
            j = pinds(i)
            noutside = 0
            call mic%window(coords(1:2,j), self%box, particles(i), noutside)
            call self%post_process(particles(i))
            call particles(i)%stats(rmean, rsdev, rmax, rmin, errout=l_err)
            if( .not.l_err )then
                cnt   = cnt + 1
                vmin  = min(vmin,rmin)
                vmax  = max(vmax,rmax)
                vmean = vmean + rmean
                vsdev = vsdev + rsdev*rsdev
            endif
        end do
        !$omp end parallel do
        if( cnt > 0 )then
            vmean = vmean / real(cnt)
            vsdev = sqrt(vsdev / real(cnt))
        endif
    end subroutine extract_particles_from_mic

    subroutine extract_all_particles_frames( self, coords, ffromto, fstep)
        class(ptcl_extractor),    intent(inout) :: self
        integer,                  intent(in)    :: coords(:,:)
        integer,                  intent(in)    :: ffromto(2), fstep
        integer :: i, f, nptcls, lb, ub
        logical :: l_allocated
        nptcls = size(coords,dim=2)
        if( ffromto(1) <= 0 .or. ffromto(2) > self%nframes .or.&
            & ffromto(1) > ffromto(2)) THROW_HARD('Invalid frame range!')
        ! Reallocate images
        l_allocated = allocated(particles_frames)
        if( l_allocated )then
            lb = lbound(particles_frames,dim=1)
            ub = ubound(particles_frames,dim=1)
            if( lb /= ffromto(1) .or. ub /= ffromto(2) .or.&
                ubound(particles_frames,dim=2) < nptcls ) then
                call dealloc_particles_frames
                l_allocated = .false.
            endif
        endif
        if( .not.l_allocated )then
            allocate(particles_frames(ffromto(1):ffromto(2), 1:nptcls))
            !$omp parallel do default(shared) private(i,f) schedule(guided) proc_bind(close) collapse(2)
            do i = 1,nptcls
                do f = ffromto(1),ffromto(2)
                    call particles_frames(f,i)%new([self%box,self%box,1],&
                        & self%smpd_out, wthreads=.false.)
                enddo
            enddo
            !$omp end parallel do
        endif
        ! extract
        !$omp parallel do default(shared) private(i,f) schedule(guided) proc_bind(close) collapse(2)
        do i = 1,nptcls
            do f = ffromto(1),ffromto(2)
                call self%extract_single_particle_frame(coords(:,i), particles_frames(f,i), f)
            enddo
        enddo
        !$omp end parallel do
    end subroutine extract_all_particles_frames

    ! Sum the supplied inclusive frame range; output is unweighted and unfiltered, then post-processed.
    subroutine generate_cumulative_particles_stack( self, nptcls, particles, ffromto, vmin, vmax, vmean, vsdev )
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: nptcls
        class(image),          intent(inout) :: particles(:)
        integer,               intent(in)    :: ffromto(2)
        real,                  intent(out)   :: vmin, vmax, vmean, vsdev
        real    :: rmin, rmax, rmean, rsdev
        integer :: i, f, cnt
        logical :: l_err
        if( ffromto(1) <= 0 .or. ffromto(2) > self%nframes .or.&
            & ffromto(1) > ffromto(2) )then
            THROW_HARD('Invalid frame range!')
        endif
        cnt    = 0
        vmin   = huge(vmin)
        vmax   = -vmin
        vmean  = 0.0
        vsdev  = 0.0
        !$omp parallel do default(shared) private(i,f,rmin,rmax,rmean,rsdev,l_err)&
        !$omp schedule(guided) reduction(+:vmean,vsdev,cnt) reduction(min:vmin)&
        !$omp reduction(max:vmax)
        do i = 1,nptcls
            call particles(i)%zero_and_unflag_ft
            do f = ffromto(1), ffromto(2)
                call particles(i)%add(particles_frames(f,i))
            enddo
            call self%post_process(particles(i))
            call particles(i)%stats(rmean, rsdev, rmax, rmin, errout=l_err)
            if( .not.l_err )then
                cnt   = cnt + 1
                vmin  = min(vmin,rmin)
                vmax  = max(vmax,rmax)
                vmean = vmean + rmean
                vsdev = vsdev + rsdev*rsdev
            endif
        enddo
        !$omp end parallel do
        if( cnt > 0 )then
            vmean = vmean / real(cnt)
            vsdev = sqrt(vsdev / real(cnt))
        endif
    end subroutine generate_cumulative_particles_stack

    subroutine extract_single_particle_frame( self, ptcl_pos_in, ptcl_out, t)
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: ptcl_pos_in(2)  ! top left corner, scaled
        class(image),          intent(inout) :: ptcl_out
        integer,               intent(in)    :: t               ! frame number/time
        real(dp) :: cx,cy
        real     :: shift(2), aniso_shift(2), center(2), rpos(2)
        integer  :: pos(2), foo, ithr
        if( any(ptcl_out%get_ldim() /= [self%box,self%box,1]) )then
            THROW_HARD('Inconsistent dimensions!')
        endif
        ithr = omp_get_thread_num() + 1
        ! particle center
        center = real(ptcl_pos_in + self%box/2) ! base 0
        cx     = center(1) + 1.d0               ! base 1
        cy     = center(2) + 1.d0               ! base 1
        ! particle corner
        rpos = center - real(self%box_pd/2)
        ! shift & coordinates
        shift = rpos - self%isoshifts(:,t)
        if( self%l_poly )then
            call self%evaluate_local_shift(t, cx, cy, aniso_shift)
            shift = shift - aniso_shift
        endif
        pos   = nint(shift)         ! extraction coordinate
        shift = shift - real(pos)   ! sub-pixel
        ! extract particle frame
        foo = 0
        call self%particle(ithr)%set_ft(.false.)
        call self%frames(t)%window(pos, self%box_pd, self%particle(ithr), foo)
        ! sub-pixel shift
        call self%particle(ithr)%fft
        call self%particle(ithr)%shift2Dserial(shift)
        ! clipping to correct size
        call self%particle(ithr)%ifft
        call self%particle(ithr)%clip(ptcl_out)
    end subroutine extract_single_particle_frame

    !>  on a single thread
    subroutine extract_ptcl( self, ptcl_pos_in, ptcl_out)
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: ptcl_pos_in(2) ! top left corner, scaled
        class(image),          intent(inout) :: ptcl_out
        real(dp)    :: cx,cy
        real        :: shift(2), aniso_shift(2), center(2), rpos(2)
        integer     :: pos(2), foo, t, ithr
        ithr = omp_get_thread_num() + 1
        ! sanity check
        if( any(ptcl_out%get_ldim() /= [self%box,self%box,1]) )then
            THROW_HARD('Inconsistent dimensions!')
        endif
        ! particle center
        center = real(ptcl_pos_in + self%box/2) ! base 0
        cx     = center(1) + 1.d0               ! base 1
        cy     = center(2) + 1.d0               ! base 1
        ! particle corner
        rpos = center - real(self%box_pd/2)
        ! extraction
        call self%particle(ithr)%zero_and_flag_ft
        do t = 1,self%nframes
            if( self%weights(t) < 1.e-6 ) cycle
            ! shift & coordinates
            shift = rpos - self%isoshifts(:,t)
            if( self%l_poly )then
                call self%evaluate_local_shift(t, cx, cy, aniso_shift)
                shift = shift - aniso_shift
            endif
            pos   = nint(shift)         ! extraction coordinate
            shift = shift - real(pos)   ! sub-pixel
            ! extract particle frame
            foo = 0
            call self%frame_particle(ithr)%set_ft(.false.)
            call self%frames(t)%window(pos, self%box_pd, self%frame_particle(ithr), foo)
            ! sub-pixel shift
            call self%frame_particle(ithr)%fft
            call self%frame_particle(ithr)%shift2Dserial(shift)
            ! weighted sum
            call self%particle(ithr)%add(self%frame_particle(ithr), w=self%weights(t))
        enddo
        ! clipping to correct size
        call self%particle(ithr)%ifft
        call self%particle(ithr)%clip(ptcl_out)
    end subroutine extract_ptcl

    !>  Sign change & normalizations
    subroutine post_process( self, img )
        class(ptcl_extractor), intent(in)    :: self
        class(image),          intent(inout) :: img
        real :: sdev_noise
        if( self%l_neg ) call img%neg()
        call img%subtr_backgr_ramp(self%particle_mask)
        call img%norm_noise(self%particle_mask, sdev_noise)
    end subroutine post_process

    !> Identify statistical outliers and cure them.
    subroutine cure_outliers( self )
        class(ptcl_extractor), intent(inout) :: self
        real,    allocatable :: rsum(:,:)
        integer, allocatable :: outlier_coords(:,:)
        real    :: ave, sdev, var, lthresh, uthresh
        integer :: iframe, noutliers, noutliers_detected, i, j, k
        logical :: outliers(self%ldim(1),self%ldim(2)), err
        allocate(rsum(self%ldim(1),self%ldim(2)), source=0.)
        write(logfhandle,'(a)') '>>> REMOVING DEAD/HOT PIXELS'
        do iframe = 1,self%nframes
            call self%frames(iframe)%add_rmat2mat_workshare(rsum)
        enddo
        call moment(rsum, ave, sdev, var, err)
        if( err .or. sdev < TINY )then
            deallocate(rsum)
            return
        endif
        lthresh = ave - NSIGMAS * sdev
        uthresh = ave + NSIGMAS * sdev
        !$omp workshare
        where( rsum < lthresh .or. rsum > uthresh )
            outliers = .true.
        elsewhere
            outliers = .false.
        end where
        !$omp end workshare
        noutliers_detected = count(outliers)
        if( self%l_eer .and. self%l_gain ) call self%gain%add_zero2mask(outliers)
        noutliers = count(outliers)
        if( noutliers > 0 )then
            write(logfhandle,'(a,1x,i10)') '>>> # DEAD/HOT PIXELS:', noutliers_detected
            if( noutliers /= noutliers_detected )then
                write(logfhandle,'(a,1x,i8)') '>>> # DEAD/HOT PIXELS + EER GAIN DEFECTS:', noutliers
            endif
            write(logfhandle,'(a,1x,2f10.1)') '>>> AVERAGE (STDEV):  ', ave, sdev
            allocate(outlier_coords(2,noutliers))
            k = 0
            do j = 1,self%ldim(2)
                do i = 1,self%ldim(1)
                    if( outliers(i,j) )then
                        k = k + 1
                        outlier_coords(:,k) = [i,j]
                    endif
                enddo
            enddo
            call self%cure_outlier_coords(outlier_coords, ave, sdev)
            deallocate(outlier_coords)
        endif
        deallocate(rsum)
    end subroutine cure_outliers

    !> Add zero-valued EER gain pixels to the known detector outlier coordinates.
    subroutine add_eer_gain_defects( self )
        class(ptcl_extractor), intent(inout) :: self
        logical, allocatable :: outliers(:,:)
        integer :: i, j, k, nmodel_outliers, noutliers
        if( .not.self%gain%exists() ) THROW_HARD('EER gain image is absent!')
        allocate(outliers(self%ldim(1),self%ldim(2)), source=.false.)
        if( allocated(self%hotpix_coords) )then
            if( size(self%hotpix_coords, dim=1) /= 2 )then
                THROW_HARD('Outlier coordinates must be a 2-by-N array!')
            endif
            do k = 1,size(self%hotpix_coords, dim=2)
                i = self%hotpix_coords(1,k)
                j = self%hotpix_coords(2,k)
                if( i < 1 .or. i > self%ldim(1) .or. j < 1 .or. j > self%ldim(2) )then
                    THROW_HARD('Outlier coordinates are outside movie dimensions!')
                endif
                outliers(i,j) = .true.
            enddo
        endif
        nmodel_outliers = count(outliers)
        call self%gain%add_zero2mask(outliers)
        noutliers = count(outliers)
        if( noutliers > nmodel_outliers )then
            write(logfhandle,'(a,1x,i10)') '>>> # EER GAIN DEFECTS ADDED:', noutliers - nmodel_outliers
        endif
        if( allocated(self%hotpix_coords) ) deallocate(self%hotpix_coords)
        self%nhotpix = noutliers
        if( noutliers > 0 )then
            allocate(self%hotpix_coords(2,noutliers))
            k = 0
            do j = 1,self%ldim(2)
                do i = 1,self%ldim(1)
                    if( outliers(i,j) )then
                        k = k + 1
                        self%hotpix_coords(:,k) = [i,j]
                    endif
                enddo
            enddo
        endif
        deallocate(outliers)
    end subroutine add_eer_gain_defects

    !> Cure the supplied one-based detector coordinates without detecting outliers.
    subroutine cure_outliers_from_coords( self, outlier_coords )
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: outlier_coords(:,:)
        real, allocatable :: rsum(:,:)
        real    :: ave, sdev, var
        integer :: iframe, noutliers
        logical :: err
        if( size(outlier_coords, dim=1) /= 2 )then
            THROW_HARD('Outlier coordinates must be a 2-by-N array!')
        endif
        noutliers = size(outlier_coords, dim=2)
        if( noutliers == 0 ) return
        allocate(rsum(self%ldim(1),self%ldim(2)), source=0.)
        do iframe = 1,self%nframes
            call self%frames(iframe)%add_rmat2mat_workshare(rsum)
        enddo
        call moment(rsum, ave, sdev, var, err)
        deallocate(rsum)
        if( err .or. sdev < TINY ) return
        write(logfhandle,'(a)') '>>> REMOVING PROVIDED DEAD/HOT PIXELS'
        write(logfhandle,'(a,1x,i10)') '>>> # DEAD/HOT PIXELS:', noutliers
        write(logfhandle,'(a,1x,2f10.1)') '>>> AVERAGE (STDEV):  ', ave, sdev
        call self%cure_outlier_coords(outlier_coords, ave, sdev)
    end subroutine cure_outliers_from_coords

    !> Replace the supplied outliers using non-outlier values from their neighborhoods.
    subroutine cure_outlier_coords( self, outlier_coords, sum_ave, sum_sdev )
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: outlier_coords(:,:)
        real,                  intent(in)    :: sum_ave, sum_sdev
        integer, parameter :: HWINSZ = 5
        integer, parameter :: NVALS  = (2 * HWINSZ + 1)**2
        real, allocatable :: new_vals(:,:)
        real    :: vals(NVALS), ave, sdev, lthresh, uthresh, l, u, localave
        integer :: iframe, noutliers, i, j, k, ii, jj, n
        logical :: outliers(self%ldim(1),self%ldim(2))
        if( size(outlier_coords, dim=1) /= 2 )then
            THROW_HARD('Outlier coordinates must be a 2-by-N array!')
        endif
        noutliers = size(outlier_coords, dim=2)
        if( noutliers == 0 ) return
        outliers = .false.
        do k = 1,noutliers
            i = outlier_coords(1,k)
            j = outlier_coords(2,k)
            if( i < 1 .or. i > self%ldim(1) .or. j < 1 .or. j > self%ldim(2) )then
                THROW_HARD('Outlier coordinates are outside movie dimensions!')
            endif
            outliers(i,j) = .true.
        enddo
        ave     = sum_ave  / real(self%nframes)
        sdev    = sum_sdev / real(self%nframes)
        lthresh = ave - NSIGMAS * sdev
        uthresh = ave + NSIGMAS * sdev
        allocate(new_vals(noutliers,self%nframes))
        !$omp parallel do default(shared) private(iframe,k,i,j,n,ii,jj,vals,l,u,localave)&
        !$omp proc_bind(close) schedule(static)
        do iframe = 1,self%nframes
            do k = 1,noutliers
                i = outlier_coords(1,k)
                j = outlier_coords(2,k)
                n = 0
                do jj = j-HWINSZ,j+HWINSZ
                    if( jj < 1 .or. jj > self%ldim(2) ) cycle
                    do ii = i-HWINSZ,i+HWINSZ
                        if( ii < 1 .or. ii > self%ldim(1) ) cycle
                        if( outliers(ii,jj) ) cycle
                        n       = n + 1
                        vals(n) = self%frames(iframe)%get([ii,jj,1])
                    enddo
                enddo
                if( n > 1 )then
                    if( real(n)/real(NVALS) < 0.85 )then
                        l = minval(vals(:n))
                        u = maxval(vals(:n))
                        if( abs(u-l) < sdev/1000.0 ) u = uthresh
                        localave = sum(vals(:n)) / real(n)
                        new_vals(k,iframe) = gasdev(localave, sdev, [l,u])
                    else
                        new_vals(k,iframe) = median_nocopy(vals(:n))
                    endif
                else
                    new_vals(k,iframe) = gasdev(ave, sdev, [lthresh,uthresh])
                endif
            enddo
            do k = 1,noutliers
                i = outlier_coords(1,k)
                j = outlier_coords(2,k)
                call self%frames(iframe)%set([i,j,1], new_vals(k,iframe))
            enddo
        enddo
        !$omp end parallel do
        deallocate(new_vals)
    end subroutine cure_outlier_coords

    pure subroutine evaluate_local_shift( self, iframe, x, y, shift )
        class(ptcl_extractor), intent(in)  :: self
        integer,               intent(in)  :: iframe
        real(dp),              intent(in)  :: x, y
        real,                  intent(out) :: shift(2)
        real(dp) :: t, xx, yy
        t = real(iframe-self%align_frame, dp)
        xx = pix2polycoords(x, self%ldim_sc(1))
        yy = pix2polycoords(y, self%ldim_sc(2))
        shift(1) = apply_patch_poly(self%polyx(:), xx,yy,t)
        shift(2) = apply_patch_poly(self%polyy(:), xx,yy,t)
    end subroutine evaluate_local_shift

    pure integer function get_nframes(self)
        class(ptcl_extractor), intent(in) :: self
        get_nframes = self%nframes
    end function get_nframes

    integer function parse_int( table, emdl_id, err )
        class(starfile_table_type) :: table
        integer, intent(in)        :: emdl_id
        logical, intent(out)       :: err
        err = .not.starfile_table__getValue_int(table, emdl_id, parse_int)
    end function parse_int

    real function parse_double( table, emdl_id, err )
        class(starfile_table_type) :: table
        integer, intent(in)        :: emdl_id
        logical, intent(out)       :: err
        real(dp) :: v
        err = .not.starfile_table__getValue_double(table, emdl_id, v)
        parse_double = real(v)
    end function parse_double

    subroutine parse_string( table, emdl_id, string, err )
        class(starfile_table_type)                 :: table
        integer,                       intent(in)  :: emdl_id
        character(len=:), allocatable, intent(out) :: string
        logical,                       intent(out) :: err
        err = .not.starfile_table__getValue_string(table, emdl_id, string)
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
        self%l_doseweighing = .false.
        self%l_gain         = .false.
        self%l_eer          = .false.
        self%l_neg          = .true.
        self%l_mov          = .true.
        self%l_model        = .false.
        self%nframes        = 0
        self%nhotpix        = 0
        self%doseperframe   = 0.
        self%scale          = 1.
        self%align_frame    = 0
        self%eer_upsampling = 1
        call self%mmodel%kill
        call self%moviename%kill
        call self%docname%kill
        call self%gainrefname%kill
        self%exists = .false.
    end subroutine kill

    subroutine dealloc_particles_frames()
        integer :: i, f, lb, ub
        if( allocated(particles_frames) ) then
            lb = lbound(particles_frames,dim=1)
            ub = ubound(particles_frames,dim=1)
            do i = 1, size(particles_frames,dim=2)
                do f = lb, ub
                    call particles_frames(f,i)%kill
                enddo
            enddo
            deallocate(particles_frames)
        endif
    end subroutine dealloc_particles_frames

end module simple_particle_extractor
