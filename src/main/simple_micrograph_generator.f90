module simple_micrograph_generator
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

public :: mic_generator

integer, parameter :: POLYDIM    = 18
real,    parameter :: NSIGMAS    = 6.0
logical, parameter :: DEBUG_HERE = .false.

type :: mic_generator
    type(image),      allocatable :: frames(:)
    real,             allocatable :: doses(:,:,:),weights(:), isoshifts(:,:)
    integer,          allocatable :: hotpix_coords(:,:)
    character(len=:), allocatable :: gainrefname, moviename, docname
    type(image)                   :: gain
    type(eer_decoder)             :: eer
    real(dp)                      :: polyx(POLYDIM), polyy(POLYDIM)
    real                          :: smpd, smpd_out
    real                          :: scale
    real                          :: total_dose, doseperframe, preexposure, kv
    integer                       :: fromtof(2)
    integer                       :: ldim(3), ldim_sc(3)
    integer                       :: nframes, start_frame, nhotpix, align_frame
    integer                       :: eer_fraction, eer_upsampling
    logical                       :: l_doseweighing = .false.
    logical                       :: l_scale        = .false.
    logical                       :: l_gain         = .false.
    logical                       :: l_eer          = .false.
    logical                       :: l_poly         = .false.
    logical                       :: l_nn_interp    = .false.
    logical                       :: l_frameweights = .false.
    logical                       :: exists         = .false.
  contains
    procedure          :: new
    procedure, private :: parse_movie_metadata
    procedure, private :: apply_dose_weighing
    procedure          :: generate_micrographs
    procedure          :: display
    procedure, private :: cure_outliers
    procedure          :: get_moviename
    ! Destructor
    procedure :: kill
end type mic_generator

contains

    !>  Constructor
    subroutine new( self, omic, convention, frames_range )
        use simple_ori, only: ori
        class(mic_generator), intent(inout) :: self
        class(ori),           intent(in)    :: omic
        character(len=*),     intent(in)    :: convention
        integer,              intent(in)    :: frames_range(2)
        character(len=:), allocatable :: poly_fname
        real(dp),         allocatable :: poly(:)
        integer  :: iframe
        call self%kill
        ! sanity checks
        if( .not. omic%isthere('mc_starfile') )then
            THROW_HARD('Movie star doc is absent 1')
        endif
        self%docname = trim(omic%get_static('mc_starfile'))
        if( .not.file_exists(self%docname) )then
            THROW_HARD('Movie star doc is absent 2')
        endif
        ! weights & inteporlation convention
        select case(trim(convention))
        case('simple')
            self%l_nn_interp    = .false.
            self%l_frameweights = .true.
        case('motioncorr')
            self%l_nn_interp    = .true.
            self%l_frameweights = .false.
        case('relion')
            self%l_nn_interp    = .false.
            self%l_frameweights = .false.
        case('cryosparc')
            self%l_nn_interp    = .true.
            self%l_frameweights = .false.
        case DEFAULT
            THROW_HARD('Unsupported convention!')
        end select
        ! get movie info
        call self%parse_movie_metadata
        if( .not.file_exists(self%moviename) ) THROW_HARD('Movie cannot be found: '//trim(self%moviename))
        ! dose-weighing
        if( self%l_doseweighing )then
            self%total_dose = real(self%nframes) * self%doseperframe
        else
            THROW_WARN('Dose-weighing metadata cannot be found in: '//trim(self%docname))
        endif
        ! frames range
        if( frames_range(1) < 1 ) THROW_HARD('Invalid starting frame')
        if( frames_range(2) > self%nframes )then
            THROW_HARD('Invalid final frame: '//int2str(frames_range(2))//'/'//int2str(self%nframes))
        endif
        self%fromtof = frames_range
        if( self%fromtof(2) == 0 ) self%fromtof(2) = self%nframes
        ! weights
        if( self%l_frameweights )then
            self%weights = self%weights / sum(self%weights(self%fromtof(1):self%fromtof(2)))
        else
            self%weights = 1.0 / real(self%fromtof(2)-self%fromtof(1)+1)
        endif
        ! frame of reference
        self%isoshifts(1,:) = self%isoshifts(1,:) - self%isoshifts(1,self%align_frame)
        self%isoshifts(2,:) = self%isoshifts(2,:) - self%isoshifts(2,self%align_frame)
        ! downscaling shifts
        self%isoshifts = self%isoshifts * self%scale
        ! polynomial coefficients
        poly_fname = fname_new_ext(self%docname,'poly')
        if( file_exists(poly_fname) )then
            poly = file2drarr(poly_fname)
            self%polyx = poly(:POLYDIM)
            self%polyy = poly(POLYDIM+1:)
        endif
        self%polyx = self%polyx * real(self%scale,dp)
        self%polyy = self%polyy * real(self%scale,dp)
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
            call self%eer%decode(self%frames, self%eer_fraction, frames_range=self%fromtof)
        else
            !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
            do iframe = self%fromtof(1),self%fromtof(2)
                call self%frames(iframe)%new(self%ldim, self%smpd, wthreads=.false.)
            enddo
            !$omp end parallel do
            do iframe = 1,self%fromtof(1),self%fromtof(2)
                call self%frames(iframe)%read(self%moviename, iframe)
            end do
        endif
        ! gain correction
        if( self%l_gain )then
            if( .not.file_exists(self%gainrefname) )then
                THROW_HARD('gain reference: '//trim(self%gainrefname)//' not found')
            endif
            if( self%l_eer )then
                call correct_gain(self%frames, self%gainrefname, self%gain, eerdecoder=self%eer, frames_range=self%fromtof)
            else
                call correct_gain(self%frames, self%gainrefname, self%gain, frames_range=self%fromtof)
            endif
        endif
        ! outliers curation
        if( self%nhotpix > 0 )then
            call self%cure_outliers
        endif
        call self%gain%kill
        ! downscale frames & dose-weighing
        !$omp parallel do schedule(guided) default(shared) private(iframe) proc_bind(close)
        do iframe = self%fromtof(1),self%fromtof(2)
            call self%frames(iframe)%fft
            call self%frames(iframe)%clip_inplace(self%ldim_sc)
        enddo
        !$omp end parallel do
        ! all done
        call self%eer%kill
        self%exists = .true.
    end subroutine new

    subroutine parse_movie_metadata( self )
        class(mic_generator), intent(inout) :: self
        type(str4arr), allocatable :: names(:)
        type(starfile_table_type)  :: table
        integer(C_long) :: num_objs, object_id
        integer         :: i,j,iframe,n,ind, motion_model
        logical         :: err
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
                self%ldim(3)        = 1
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
                THROW_HARD('Invalid table: '//trim(names(i)%str))
            end select
        enddo
        call starfile_table__delete(table)
        if( DEBUG_HERE ) print *,'movie doc parsed'
    end subroutine parse_movie_metadata

    !>  Regenerates micrograph from movie
    subroutine generate_micrographs( self, micrograph_dw, micrograph_nodw )
        class(mic_generator), intent(inout) :: self
        type(image),          intent(inout) :: micrograph_dw, micrograph_nodw
        real,        pointer     :: rmatin(:,:,:), rmatout(:,:,:)
        type(image), allocatable :: local_frames(:)
        type(image)              :: gridcorrimg
        integer                  :: ldim(3), iframe

        integer :: lims(3,2), phys(3), h,k
        real    :: pid, dsq
        allocate(local_frames(self%fromtof(1):self%fromtof(2)))
        ldim = self%ldim_sc
        ! Stage drift
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe = self%fromtof(1),self%fromtof(2)
            call self%frames(iframe)%fft
            call self%frames(iframe)%shift2Dserial(-self%isoshifts(:,iframe))
            call local_frames(iframe)%copy(self%frames(iframe))
            call local_frames(iframe)%ifft
        end do
        !$omp end parallel do
        ! Beam-induced motion
        ! non dose-weighted
        call bimc(micrograph_nodw)
        ! dose-weighted
        if( self%l_doseweighing )then
            call self%apply_dose_weighing
            !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
            do iframe = self%fromtof(1),self%fromtof(2)
                call local_frames(iframe)%copy_fast(self%frames(iframe))
                call local_frames(iframe)%ifft
            end do
            !$omp end parallel do
            call bimc(micrograph_dw)
        endif
        ! cleanup
        !$omp parallel do default(shared) private(iframe) proc_bind(close) schedule(static)
        do iframe = self%fromtof(1),self%fromtof(2)
            call local_frames(iframe)%kill
        end do
        !$omp end parallel do
        deallocate(local_frames)
    contains

        subroutine bimc( mic )
            class(image), intent(inout) :: mic
            real(dp)      :: t,ti, dt,dt2,dt3, x,x2,y,y2,xy, A1,A2
            real(dp)      :: B1x,B1x2,B1xy,B2x,B2x2,B2xy
            integer       :: i, j
            real          :: w, pixx,pixy
            call mic%new(ldim, self%smpd_out)
            call mic%zero_and_unflag_ft
            call mic%get_rmat_ptr(rmatout)
            ti = 0.d0
            do iframe = self%fromtof(1),self%fromtof(2)
                call local_frames(iframe)%get_rmat_ptr(rmatin)
                w = self%weights(iframe)
                t = real(iframe-self%align_frame, dp)
                dt  = ti-t
                dt2 = ti*ti - t*t
                dt3 = ti*ti*ti - t*t*t
                B1x  = sum(self%polyx(4:6)   * [dt,dt2,dt3])
                B1x2 = sum(self%polyx(7:9)   * [dt,dt2,dt3])
                B1xy = sum(self%polyx(16:18) * [dt,dt2,dt3])
                B2x  = sum(self%polyy(4:6)   * [dt,dt2,dt3])
                B2x2 = sum(self%polyy(7:9)   * [dt,dt2,dt3])
                B2xy = sum(self%polyy(16:18) * [dt,dt2,dt3])
                !$omp parallel do default(shared) private(i,j,x,x2,y,y2,xy,A1,A2,pixx,pixy)&
                !$omp proc_bind(close) schedule(static)
                do j = 1, ldim(2)
                    y  = real(j-1,dp) / real(ldim(2)-1,dp) - 0.5d0
                    y2 = y*y
                    A1 =           sum(self%polyx(1:3)   * [dt,dt2,dt3])
                    A1 = A1 + y  * sum(self%polyx(10:12) * [dt,dt2,dt3])
                    A1 = A1 + y2 * sum(self%polyx(13:15) * [dt,dt2,dt3])
                    A2 =           sum(self%polyy(1:3)   * [dt,dt2,dt3])
                    A2 = A2 + y  * sum(self%polyy(10:12) * [dt,dt2,dt3])
                    A2 = A2 + y2 * sum(self%polyy(13:15) * [dt,dt2,dt3])
                    do i = 1, ldim(1)
                        x  = real(i-1,dp) / real(ldim(1)-1,dp) - 0.5d0
                        x2 = x*x
                        xy = x*y
                        pixx = real(i) + real(A1 + B1x*x + B1x2*x2 + B1xy*xy)
                        pixy = real(j) + real(A2 + B2x*x + B2x2*x2 + B2xy*xy)
                        if( self%l_nn_interp )then
                            rmatout(i,j,1) = rmatout(i,j,1) + w*interp_nn(pixx,pixy)
                        else
                            rmatout(i,j,1) = rmatout(i,j,1) + w*interp_bilin(pixx,pixy)
                        endif
                    end do
                end do
                !$omp end parallel do
            enddo
        end subroutine bimc

        pure real function interp_bilin( xval, yval )
            real, intent(in) :: xval, yval
            integer  :: x1_h,  x2_h,  y1_h,  y2_h
            real     :: t, u
            logical  :: outside
            outside = .false.
            x1_h = floor(xval)
            x2_h = x1_h + 1
            if( x1_h<1 .or. x2_h<1 )then
                x1_h    = 1
                outside = .true.
            endif
            if( x1_h>ldim(1) .or. x2_h>ldim(1) )then
                x1_h    = ldim(1)
                outside = .true.
            endif
            y1_h = floor(yval)
            y2_h = y1_h + 1
            if( y1_h<1 .or. y2_h<1 )then
                y1_h    = 1
                outside = .true.
            endif
            if( y1_h>ldim(2) .or. y2_h>ldim(2) )then
                y1_h    = ldim(2)
                outside = .true.
            endif
            if( outside )then
                interp_bilin = rmatin(x1_h, y1_h, 1)
                return
            endif
            t  = xval - real(x1_h)
            u  = yval - real(y1_h)
            interp_bilin =  (1. - t) * (1. - u) * rmatin(x1_h, y1_h, 1) + &
                                 &t  * (1. - u) * rmatin(x2_h, y1_h, 1) + &
                                 &t  *       u  * rmatin(x2_h, y2_h, 1) + &
                           &(1. - t) *       u  * rmatin(x1_h, y2_h, 1)
        end function interp_bilin

        pure real function interp_nn( xval, yval )
            real, intent(in) :: xval, yval
            integer  :: x, y
            x = min(max(nint(xval),1),ldim(1))
            y = min(max(nint(yval),1),ldim(2))
            interp_nn = rmatin(x,y,1)
        end function interp_nn

    end subroutine generate_micrographs

    subroutine display( self )
        class(mic_generator), intent(in) :: self
        integer :: i
        print *, 'docname         ', self%docname
        print *, 'nframes         ', self%nframes
        print *, 'dimensions      ', self%ldim
        print *, 'smpd            ', self%smpd
        print *, 'smpd_out        ', self%smpd_out
        print *, 'voltage         ', self%kv
        print *, 'doseperframe    ', self%doseperframe
        print *, 'gainrefname     ', trim(self%gainrefname)
        print *, 'moviename       ', trim(self%moviename)
        print *, 'doseweighting   ', self%l_doseweighing
        print *, 'total dose      ', self%total_dose
        print *, 'l_scale         ', self%l_scale
        print *, 'scale           ', self%scale
        print *, 'gain            ', self%l_gain
        print *, 'nhotpix         ', self%nhotpix
        print *, 'eer             ', self%l_eer
        print *, 'eer_fraction    ', self%eer_fraction
        print *, 'eer_upsampling  ', self%eer_upsampling
        print *, 'align_frame     ', self%align_frame
        if( allocated(self%isoshifts) )then
            do i = 1,size(self%isoshifts,dim=2)
                print *,'isoshifts ',i,self%isoshifts(:,i),self%weights(i)
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

    subroutine cure_outliers( self )
        class(mic_generator), intent(inout)  :: self
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
        do iframe = self%fromtof(1),self%fromtof(2)
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
            do iframe = self%fromtof(1),self%fromtof(2)
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
        class(mic_generator), intent(inout) :: self
        real, parameter :: A=0.245, B=-1.665, C=2.81
        real            :: qs(self%fromtof(1):self%fromtof(2)), acc_doses(self%nframes)
        real            :: spaFreqk, twoNe, spafreq, limhsq,limksq
        integer         :: nrflims(3,2), ldim(3), hphys,kphys, iframe, h,k
        if( .not.self%l_doseweighing ) return
        if( .not.self%frames(self%fromtof(1))%is_ft() ) THROW_HARD('Frames should be in in the Fourier domain')
        nrflims = self%frames(self%fromtof(1))%loop_lims(2)
        ldim    = self%frames(self%fromtof(1))%get_ldim()
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
                qs = exp(-acc_doses(self%fromtof(1):self%fromtof(2))/twoNe)
                qs = qs / sqrt(sum(qs*qs))
                do iframe = self%fromtof(1),self%fromtof(2)
                    call self%frames(iframe)%mul_cmat_at([hphys,kphys,1], qs(iframe))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine apply_dose_weighing

    pure function get_moviename( self )result( fname )
        class(mic_generator), intent(in)  :: self
        character(len=:), allocatable :: fname
        fname = trim(adjustl(self%moviename))
    end function get_moviename

    ! Utilities

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
        logical, intent(out)                       :: err
        err = .not.starfile_table__getValue_string(table, emdl_id, string)
    end subroutine parse_string
    
    subroutine kill(self)
        class(mic_generator), intent(inout) :: self
        integer :: i
        if( allocated(self%weights) )   deallocate(self%weights)
        if( allocated(self%isoshifts) ) deallocate(self%isoshifts)
        if( allocated(self%frames) )then
            do i=1,self%nframes
                call self%frames(i)%kill
            enddo
            deallocate(self%frames)
        endif
        if( allocated(self%doses) )         deallocate(self%doses)
        if( allocated(self%hotpix_coords) ) deallocate(self%hotpix_coords)
        self%l_doseweighing = .false.
        self%l_gain         = .false.
        self%l_eer          = .false.
        self%l_nn_interp    = .false.
        self%l_frameweights = .false.
        self%nframes        = 0
        self%nhotpix        = 0
        self%doseperframe   = 0.
        self%scale          = 1.
        self%align_frame    = 0
        self%moviename   = ''
        self%docname     = ''
        self%gainrefname = ''
        self%exists = .false.
    end subroutine kill

end module simple_micrograph_generator