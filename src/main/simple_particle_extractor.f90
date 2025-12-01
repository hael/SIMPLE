module simple_particle_extractor
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,                         only: image, image_ptr
use simple_eer_factory,                   only: eer_decoder
use simple_motion_correct_utils,          only: correct_gain
use simple_starfile_wrappers
implicit none
private
#include "simple_local_flags.inc"

public :: ptcl_extractor

integer, parameter :: POLYDIM    = 18
real,    parameter :: NSIGMAS    = 6.0
logical, parameter :: DEBUG_HERE = .false.

type :: ptcl_extractor
    type(image),      allocatable :: frames(:)
    type(image),      allocatable :: particle(:), frame_particle(:)
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
    integer                       :: eer_fraction, eer_upsampling
    logical                       :: l_doseweighing = .false.
    logical                       :: l_scale        = .false.
    logical                       :: l_gain         = .false.
    logical                       :: l_eer          = .false.
    logical                       :: l_poly         = .false.
    logical                       :: l_neg          = .true.
    logical                       :: l_mov          = .true.
    logical                       :: exists         = .false.
  contains
    procedure          :: init_mov
    procedure          :: init_mic
    procedure, private :: init_mask
    procedure, private :: parse_movie_metadata
    procedure, private :: apply_dose_weighing
    procedure          :: display
    procedure          :: extract_particles
    procedure          :: extract_particles_from_mic
    procedure, private :: extract_ptcl
    procedure, private :: post_process
    procedure, private :: cure_outliers
    procedure, private :: pix2polycoords
    procedure, private :: get_local_shift
    ! Destructor
    procedure :: kill
end type ptcl_extractor

contains

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
                call self%apply_dose_weighing
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
                call self%get_local_shift(t, cx, cy, aniso_shift)
                shift = shift - aniso_shift
            endif
            pos   = nint(shift)         ! extraction coordinate
            shift = (shift - real(pos)) ! sub-pixel
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

    subroutine cure_outliers( self )
        class(ptcl_extractor), intent(inout)  :: self
        integer, parameter   :: hwinsz = 5
        type(image_ptr)      :: prmats(self%nframes)
        real,        pointer :: prmat(:,:,:)
        real,    allocatable :: rsum(:,:), new_vals(:,:), vals(:)
        integer, allocatable :: pos_outliers(:,:), pos_outliers_here(:,:)
        real    :: ave, sdev, var, lthresh,uthresh, l,u,localave
        integer :: iframe, noutliers, i,j,k,ii,jj, nvals, winsz, n
        logical :: outliers(self%ldim(1),self%ldim(2)), err
        allocate(rsum(self%ldim(1),self%ldim(2)),source=0.)
        write(logfhandle,'(a)') '>>> REMOVING DEAD/HOT PIXELS'
        ! sum
        do iframe = 1,self%nframes
            call self%frames(iframe)%get_rmat_ptr(prmat)
            !$omp parallel workshare
            rsum(:,:) = rsum(:,:) + prmat(:self%ldim(1),:self%ldim(2),1)
            !$omp end parallel workshare
        enddo
        nullify(prmat)
        ! outliers detection
        call moment( rsum, ave, sdev, var, err )
        if( sdev<TINY )return
        lthresh   = ave - NSIGMAS * sdev
        uthresh   = ave + NSIGMAS * sdev
        !$omp workshare
        where( rsum<lthresh .or. rsum>uthresh )
            outliers = .true.
        elsewhere
            outliers = .false.
        end where
        !$omp end workshare
        deallocate(rsum)
        ! cure
        noutliers = count(outliers)
        if( noutliers > 0 )then
            write(logfhandle,'(a,1x,i10)') '>>> # DEAD/HOT PIXELS:', noutliers
            write(logfhandle,'(a,1x,2f10.1)') '>>> AVERAGE (STDEV):  ', ave, sdev
            winsz = 2*HWINSZ+1
            nvals = winsz*winsz
            allocate(pos_outliers(2,noutliers))
            ! gather positions for star output only
            k = 0
            do j = 1,self%ldim(2)
                do i = 1,self%ldim(1)
                    if( outliers(i,j) )then
                        k = k + 1
                        pos_outliers(:,k) = [i,j]
                    endif
                enddo
            enddo
            ! add eer gain defects for curation
            if( self%l_eer .and. self%l_gain )then
                call self%gain%get_rmat_ptr(prmat)
                where( is_zero(prmat(:self%ldim(1),:self%ldim(2),1)) )
                    outliers = .true.
                end where
                nullify(prmat)
                noutliers = count(outliers)
                if( noutliers > 0 )then
                    write(logfhandle,'(a,1x,i8)') '>>> # DEAD/HOT PIXELS + EER GAIN DEFFECTS:', noutliers
                    ! gather defect positions again for curation
                    allocate(pos_outliers_here(2,noutliers),source=-1)
                    k = 0
                    do j = 1,self%ldim(2)
                        do i = 1,self%ldim(1)
                            if( outliers(i,j) )then
                                k = k + 1
                                pos_outliers_here(:,k) = [i,j]
                            endif
                        enddo
                    enddo
                else
                    ! nothing to do
                    return
                endif
            else
                pos_outliers_here = pos_outliers
            endif
            allocate(new_vals(noutliers,self%nframes),vals(nvals))
            ave  = ave / real(self%nframes)
            sdev = sdev / real(self%nframes)
            uthresh = uthresh / real(self%nframes)
            lthresh = lthresh / real(self%nframes)
            !$omp parallel do default(shared) private(iframe,k,i,j,n,ii,jj,vals,l,u,localave)&
            !$omp proc_bind(close) schedule(static)
            do iframe=1,self%nframes
                call self%frames(iframe)%get_rmat_ptr(prmats(iframe)%rmat)
                ! calulate new values
                do k = 1,noutliers
                    i = pos_outliers_here(1,k)
                    j = pos_outliers_here(2,k)
                    n = 0
                    do jj = j-HWINSZ,j+HWINSZ
                        if( jj < 1 .or. jj > self%ldim(2) ) cycle
                        do ii = i-HWINSZ,i+HWINSZ
                            if( ii < 1 .or. ii > self%ldim(1) ) cycle
                            if( outliers(ii,jj) ) cycle
                            n = n + 1
                            vals(n) = prmats(iframe)%rmat(ii,jj,1)
                        enddo
                    enddo
                    if( n > 1 )then
                        if( real(n)/real(nvals) < 0.85 )then
                            ! high defect area
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
                ! substitute
                do k = 1,noutliers
                    i = pos_outliers_here(1,k)
                    j = pos_outliers_here(2,k)
                    prmats(iframe)%rmat(i,j,1) = new_vals(k,iframe)
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine cure_outliers

    ! Frames assumed in fourier space
    subroutine apply_dose_weighing( self )
        class(ptcl_extractor), intent(inout) :: self
        real, parameter :: A=0.245, B=-1.665, C=2.81
        real            :: qs(self%nframes), acc_doses(self%nframes)
        real            :: spaFreqk, twoNe, spafreq, limhsq,limksq
        integer         :: nrflims(3,2), ldim(3), hphys,kphys, iframe, h,k
        if( .not.self%l_doseweighing ) return
        if( .not.self%frames(1)%is_ft() ) THROW_HARD('Frames should be in in the Fourier domain')
        nrflims = self%frames(1)%loop_lims(2)
        ldim    = self%frames(1)%get_ldim()
        limhsq  = (real(ldim(1))*self%smpd_out)**2.
        limksq  = (real(ldim(2))*self%smpd_out)**2.
        do iframe=1,self%nframes
            acc_doses(iframe) = real(iframe) * self%doseperframe
        end do
        if( is_equal(self%kV,200.) )then
            acc_doses = acc_doses / 0.8
        else if( is_equal(self%kV,100.) )then
            acc_doses = acc_doses / 0.64
        endif
        !$omp parallel do private(h,k,spafreq,spafreqk,twone,kphys,hphys,iframe,qs)&
        !$omp default(shared) schedule(static) proc_bind(close)
        do k = nrflims(2,1),nrflims(2,2)
            kphys    = k + 1 + merge(ldim(2),0,k<0)
            spaFreqk = real(k*k)/limksq
            do h = nrflims(1,1),nrflims(1,2)
                hphys   = h + 1
                spaFreq = sqrt( real(h*h)/limhsq + spaFreqk )
                twoNe   = 2.*(A*spaFreq**B + C)
                qs = exp(-acc_doses/twoNe)
                qs = qs / sqrt(sum(qs*qs))
                do iframe = 1,self%nframes
                    call self%frames(iframe)%mul_cmat_at([hphys,kphys,1], qs(iframe))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_dose_weighing

    !>  pixels to coordinates for polynomial evaluation (scaled in/out)
    elemental subroutine pix2polycoords( self, xin, yin, x, y )
        class(ptcl_extractor), intent(in)  :: self
        real(dp),              intent(in)  :: xin, yin
        real(dp),              intent(out) :: x, y
        x = (xin-1.d0) / real(self%ldim_sc(1)-1,dp) - 0.5d0
        y = (yin-1.d0) / real(self%ldim_sc(2)-1,dp) - 0.5d0
    end subroutine pix2polycoords

    pure subroutine get_local_shift( self, iframe, x, y, shift )
        class(ptcl_extractor), intent(in)  :: self
        integer,               intent(in)  :: iframe
        real(dp),              intent(in)  :: x, y
        real,                  intent(out) :: shift(2)
        real(dp) :: t, xx, yy
        t = real(iframe-self%align_frame, dp)
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
        self%nframes        = 0
        self%nhotpix        = 0
        self%doseperframe   = 0.
        self%scale          = 1.
        self%align_frame    = 0
        call self%moviename%kill
        call self%docname%kill
        call self%gainrefname%kill
        self%exists = .false.
    end subroutine kill

end module simple_particle_extractor