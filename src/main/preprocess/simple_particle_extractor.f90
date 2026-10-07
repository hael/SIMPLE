!@descr: particle extraction from micrographs and motion-corrected movie frames
module simple_particle_extractor
use simple_core_module_api
use simple_image,                only: image
use simple_eer_factory,          only: eer_decoder
use simple_motion_correct_utils, only: correct_gain, pix2polycoords, apply_patch_poly
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
    private
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
    real                          :: smpd_physical = 0. ! original input sampling, before any EER upsampling or Fourier cropping
    real                          :: smpd_out
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
    procedure, private :: load_model_metadata, prepare_motion
    procedure, private :: read_movie_frames, apply_gain_correction, prepare_frames
    procedure, private :: init_particle_buffers
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
    procedure          :: get_box
    procedure          :: from_mov, from_model
    procedure          :: get_l_neg
    procedure          :: does_exist
    procedure          :: get_particle_mask
    procedure          :: has_frames, has_isoshifts
    procedure          :: has_weights, has_model_frameweights
    procedure          :: has_hotpix_coords
    procedure          :: set_test_movie_state
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
        call self%kill
        if( .not. omov%isthere('mcmodel') ) THROW_HARD('Motion model entry is absent!')
        call omov%getter('mcmodel',self%docname)
        if( .not.file_exists(self%docname) ) THROW_HARD('Motion model doc is absent!')
        self%l_mov = .false.
        self%l_neg = (params%pcontrast .eq. 'black')
        self%box   = params%box
        call self%load_model_metadata(params)
        if( .not.file_exists(self%moviename) ) THROW_HARD('Movie is absent!')
        if( params%tof > self%nframes ) THROW_HARD('TOF is large than the number of frames!')
        self%total_dose = self%mmodel%accumulated_dose
        call self%prepare_motion(self%mmodel%fixed_frame)
        write(logfhandle,'(A,A)')'>> PARSED MODEL: ',self%docname%to_char()
        call self%read_movie_frames
        write(logfhandle,'(A,A)')'>> PARSED MOVIE: ',self%moviename%to_char()
        call self%apply_gain_correction
        ! Model extraction repairs recorded defects, including zeros in the oriented EER gain.
        if( self%l_gain .and. self%l_eer ) call self%add_eer_gain_defects
        if( self%nhotpix > 0 ) call self%cure_outliers_from_coords(self%hotpix_coords)
        call self%gain%kill
        write(logfhandle,'(A,A)')'>> CORRECTED FOR GAIN AND OUTLIERS: ',self%moviename%to_char()
        call self%prepare_frames(apply_dose_weighting=.false.)
        call self%init_particle_buffers
        call self%eer%kill
        self%exists = .true.
    end subroutine init_model

    !>  Constructor
    subroutine init_mov( self, omic, params )
        use simple_ori,        only: ori
        use simple_parameters, only: parameters
        class(ptcl_extractor),     intent(inout) :: self
        class(ori),                intent(in)    :: omic
        class(parameters), target, intent(in)    :: params
        integer :: reference_frame
        call self%kill
        self%l_mov = .false.
        if( omic%isthere('mcmodel') )then
            call omic%getter('mcmodel', self%docname)
            if( .not.self%docname%is_blank() ) self%l_mov = file_exists(self%docname)
        endif
        if( .not.self%l_mov )then
            THROW_WARN('Motion model absent; extracting from the integrated micrograph')
            call self%init_mic(params%box, (params%pcontrast .eq. 'black'))
            self%exists = .true.
            return
        endif
        self%l_neg = (params%pcontrast .eq. 'black')
        self%box   = params%box
        call self%load_model_metadata(params)
        if( .not.file_exists(self%moviename) )then
            THROW_HARD('Motion model movie is absent: '//self%moviename%to_char())
        endif
        ! Global-only movie extraction retains its first-frame reference.
        reference_frame = 1
        if( self%l_poly ) reference_frame = self%mmodel%fixed_frame
        call self%prepare_motion(reference_frame)
        self%total_dose = real(self%nframes) * self%doseperframe
        call self%read_movie_frames
        call self%apply_gain_correction
        ! Movie extraction re-detects defects when the model recorded outliers.
        if( self%nhotpix > 0 ) call self%cure_outliers
        call self%gain%kill
        call self%prepare_frames(apply_dose_weighting=.true.)
        call self%init_particle_buffers
        call self%eer%kill
        self%exists = .true.
    end subroutine init_mov

    subroutine init_mic( self, box, neg )
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: box
        logical,               intent(in)    :: neg
        self%box   = box
        self%l_neg = neg
        self%l_mov = .false.
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

    subroutine load_model_metadata( self, params )
        use simple_parameters, only: parameters
        class(ptcl_extractor),     intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        call self%mmodel%read(self%docname, params)
        if( self%mmodel%nframes < 1 ) THROW_HARD('Motion model contains no movie frames!')
        if( self%mmodel%binning <= 0. ) THROW_HARD('Invalid motion model binning!')
        self%l_model       = .true.
        self%moviename     = self%mmodel%movie
        self%gainrefname   = self%mmodel%gain
        self%l_gain        = .not.self%gainrefname%is_blank()
        self%ldim          = [self%mmodel%ldim_movie, 1]
        self%ldim_sc       = [self%mmodel%ldim,       1]
        self%smpd_out      = self%mmodel%smpd
        self%scale         = 1. / self%mmodel%binning
        self%l_scale       = abs(self%scale - 1.) > 0.001
        self%nframes       = self%mmodel%nframes
        self%doseperframe  = self%mmodel%dose_per_frame
        self%l_doseweighing = self%mmodel%dw
        self%preexposure   = 0.
        self%kv            = self%mmodel%voltage
        self%start_frame   = 1
        self%l_eer         = self%mmodel%eer
        self%eer_upsampling = self%mmodel%eer_upsampling
        self%eer_fraction   = self%mmodel%eer_fraction
        ! The EER decoder needs physical sampling, not the saved decoded sampling.
        self%smpd_physical = self%mmodel%smpd_movie
        if( self%l_eer )then
            if( self%eer_upsampling /= 1 .and. self%eer_upsampling /= 2 ) THROW_HARD('Unsupported EER up-sampling!')
            if( self%eer_upsampling == 2 ) self%smpd_physical = 2. * self%smpd_physical
        endif
        allocate(self%weights(self%nframes))
        if( allocated(self%mmodel%frameweights) )then
            if( size(self%mmodel%frameweights) /= self%nframes )then
                THROW_HARD('Motion model frame-weight dimensions do not match!')
            endif
            self%weights = self%mmodel%frameweights
        else
            self%weights = 1. / real(self%nframes)
        endif
        self%l_poly = self%mmodel%npatch > 0 .and. allocated(self%mmodel%local_offsets_x) .and.&
            &allocated(self%mmodel%local_offsets_y) .and. self%mmodel%patch_accepted
        self%nhotpix = 0
        if( allocated(self%mmodel%outlier_coords) )then
            self%nhotpix = size(self%mmodel%outlier_coords, dim=2)
            self%hotpix_coords = self%mmodel%outlier_coords
        endif
        if( DEBUG_HERE ) print *,'movie model parsed'
    end subroutine load_model_metadata

    subroutine prepare_motion( self, reference_frame )
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: reference_frame
        if( reference_frame < 1 .or. reference_frame > self%nframes )then
            THROW_HARD('Motion model reference frame is out of range!')
        endif
        self%align_frame = reference_frame
        allocate(self%isoshifts(2,self%nframes), source=0.)
        if( allocated(self%mmodel%drift_offsets_x) .and. allocated(self%mmodel%drift_offsets_y) )then
            if( size(self%mmodel%drift_offsets_x) /= self%nframes .or.&
                &size(self%mmodel%drift_offsets_y) /= self%nframes )then
                THROW_HARD('Motion model drift-offset dimensions do not match!')
            endif
            self%isoshifts(1,:) = self%mmodel%drift_offsets_x
            self%isoshifts(2,:) = self%mmodel%drift_offsets_y
        endif
        self%polyx = 0.d0
        self%polyy = 0.d0
        if( self%l_poly )then
            call self%mmodel%refit_polynomial(self%align_frame)
            self%polyx = self%mmodel%model_coeffs_x
            self%polyy = self%mmodel%model_coeffs_y
        endif
        ! Shifts and refitted coefficients already use scaled pixels.
        self%isoshifts(1,:) = self%isoshifts(1,:) - self%isoshifts(1,self%align_frame)
        self%isoshifts(2,:) = self%isoshifts(2,:) - self%isoshifts(2,self%align_frame)
    end subroutine prepare_motion

    subroutine read_movie_frames( self )
        class(ptcl_extractor), intent(inout) :: self
        integer :: iframe
        allocate(self%frames(self%nframes))
        if( self%l_eer )then
            call self%eer%new(self%moviename, self%smpd_physical, self%eer_upsampling)
            if( any(self%eer%get_ldim() /= self%ldim) ) THROW_HARD('EER dimensions do not match the motion model!')
            if( .not.is_equal(self%eer%get_smpd_out(), self%mmodel%smpd_movie) )then
                THROW_HARD('EER sampling does not match the motion model!')
            endif
            call self%eer%decode(self%frames, self%eer_fraction)
        else
            !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
            do iframe = 1,self%nframes
                call self%frames(iframe)%new(self%ldim, self%smpd_physical, wthreads=.false.)
            enddo
            !$omp end parallel do
            do iframe = 1,self%nframes
                call self%frames(iframe)%read(self%moviename, iframe)
            enddo
        endif
    end subroutine read_movie_frames

    ! Keep the oriented gain alive until the constructor has repaired detector defects.
    subroutine apply_gain_correction( self )
        class(ptcl_extractor), intent(inout) :: self
        if( .not.self%l_gain ) return
        if( .not.file_exists(self%gainrefname) )then
            THROW_HARD('gain reference: '//self%gainrefname%to_char()//' not found')
        endif
        if( self%l_eer )then
            call correct_gain(self%frames, self%gainrefname, self%gain,&
                &eerdecoder=self%eer, flipgain=self%mmodel%flipgain)
        else
            call correct_gain(self%frames, self%gainrefname, self%gain, flipgain=self%mmodel%flipgain)
        endif
    end subroutine apply_gain_correction

    subroutine prepare_frames( self, apply_dose_weighting )
        class(ptcl_extractor), intent(inout) :: self
        logical,               intent(in)    :: apply_dose_weighting
        integer :: iframe
        if( .not.apply_dose_weighting .and. all(self%ldim == self%ldim_sc) ) return
        !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
        do iframe = 1,self%nframes
            call self%frames(iframe)%fft
            call self%frames(iframe)%clip_inplace(self%ldim_sc)
            if( .not.apply_dose_weighting ) call self%frames(iframe)%ifft
        enddo
        !$omp end parallel do
        if( self%l_eer )then
            if( .not.is_equal(self%frames(1)%get_smpd(), self%smpd_out) )then
                THROW_HARD('Scaled EER sampling does not match the motion model!')
            endif
        endif
        if( apply_dose_weighting )then
            call self%frames(1)%apply_dose_weighing(self%nframes, self%frames,&
                &[1,self%nframes], self%total_dose, self%kv)
            !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
            do iframe = 1,self%nframes
                call self%frames(iframe)%ifft
            enddo
            !$omp end parallel do
        endif
    end subroutine prepare_frames

    subroutine init_particle_buffers( self )
        class(ptcl_extractor), intent(inout) :: self
        integer :: i
        self%box_pd = find_larger_magic_box(self%box+1)
        allocate(self%frame_particle(nthr_glob), self%particle(nthr_glob))
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1,nthr_glob
            call self%particle(i)%new(      [self%box_pd,self%box_pd,1], self%smpd_out, wthreads=.false.)
            call self%frame_particle(i)%new([self%box_pd,self%box_pd,1], self%smpd_out, wthreads=.false.)
        enddo
        !$omp end parallel do
        call self%init_mask
    end subroutine init_particle_buffers

    subroutine display( self )
        class(ptcl_extractor), intent(in) :: self
        integer :: i
        print *, 'docname        ', self%docname%to_char()
        print *, 'nframes        ', self%nframes
        print *, 'dimensions     ', self%ldim
        print *, 'smpd_movie     ', self%mmodel%smpd_movie
        print *, 'smpd_physical  ', self%smpd_physical
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

    pure integer function get_box(self)
        class(ptcl_extractor), intent(in) :: self
        get_box = self%box
    end function get_box

    pure logical function from_mov(self)
        class(ptcl_extractor), intent(in) :: self
        from_mov = self%l_mov
    end function from_mov

    pure logical function from_model(self)
        class(ptcl_extractor), intent(in) :: self
        from_model = self%l_model
    end function from_model

    pure logical function get_l_neg(self)
        class(ptcl_extractor), intent(in) :: self
        get_l_neg = self%l_neg
    end function get_l_neg

    pure logical function does_exist(self)
        class(ptcl_extractor), intent(in) :: self
        does_exist = self%exists
    end function does_exist

    pure subroutine get_particle_mask(self, mask)
        class(ptcl_extractor), intent(in) :: self
        logical, allocatable, intent(out) :: mask(:,:,:)
        if( allocated(self%particle_mask) ) mask = self%particle_mask
    end subroutine get_particle_mask

    pure logical function has_frames(self)
        class(ptcl_extractor), intent(in) :: self
        has_frames = allocated(self%frames)
    end function has_frames

    pure logical function has_isoshifts(self)
        class(ptcl_extractor), intent(in) :: self
        has_isoshifts = allocated(self%isoshifts)
    end function has_isoshifts

    pure logical function has_weights(self)
        class(ptcl_extractor), intent(in) :: self
        has_weights = allocated(self%weights)
    end function has_weights

    pure logical function has_hotpix_coords(self)
        class(ptcl_extractor), intent(in) :: self
        has_hotpix_coords = allocated(self%hotpix_coords)
    end function has_hotpix_coords

    pure logical function has_model_frameweights(self)
        class(ptcl_extractor), intent(in) :: self
        has_model_frameweights = allocated(self%mmodel%frameweights)
    end function has_model_frameweights

    ! Allocation-only fixture for lifecycle tests; not a movie extraction initializer.
    subroutine set_test_movie_state(self, nframes)
        class(ptcl_extractor), intent(inout) :: self
        integer,               intent(in)    :: nframes
        if( nframes < 1 ) THROW_HARD('Invalid frame count for particle extractor test fixture')
        call self%kill
        self%nframes = nframes
        self%l_mov   = .true.
        self%l_model = .true.
        self%exists  = .true.
        allocate(self%frames(nframes))
        allocate(self%isoshifts(2,nframes), source=0.)
        allocate(self%weights(nframes), source=1./real(nframes))
        allocate(self%hotpix_coords(2,1), source=1)
        self%nhotpix = 1
        self%mmodel%exists  = .true.
        self%mmodel%nframes = nframes
        allocate(self%mmodel%frameweights(nframes), source=1./real(nframes))
    end subroutine set_test_movie_state

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
        self%smpd_physical  = 0.
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
