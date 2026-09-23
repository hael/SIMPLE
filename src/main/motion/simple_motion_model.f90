!@descr: Motion correction polynomial model
module simple_motion_model
use, intrinsic :: iso_fortran_env, only: int8, int32
use simple_core_module_api
use simple_parameters, only: parameters
use simple_image,      only: image

implicit none

public :: motion_model
private
#include "simple_local_flags.inc"

integer, parameter :: MODEL_VERSION = 0    ! Attempt at versioning the model
integer, parameter :: FILE_VERSION  = 0    ! Attempt at versioning file format
integer, parameter :: MODELSZ       = 18   ! number of the models polynomial coefficients

type :: motion_model
    private
    class(parameters), pointer :: p_ptr => null()
    ! Files
    type(string)         :: file_name
    type(string)         :: movie
    type(string)         :: gain
    ! Images dimensions
    real                 :: smpd_movie
    integer              :: ldim_movie(2)
    real                 :: smpd
    integer              :: ldim(2)
    real                 :: binning
    integer              :: nframes
    integer              :: total_nframes
    ! Dose
    real                 :: voltage
    real                 :: dose_per_frame
    real                 :: target_dose_per_frame
    real                 :: total_dose, accumulated_dose
    logical              :: dw
    ! EER
    logical              :: eer
    integer              :: eer_fraction
    ! Stage drift
    real,    allocatable :: drift_offsets_x(:), drift_offsets_y(:)
    real,    allocatable :: frameweights(:)
    ! Local motion & patches
    integer              :: npatch, nx_patch, ny_patch
    integer              :: ldim_patch(2)
    integer, allocatable :: patch_bounds(:,:,:,:)                               ! Patch boundaries
    real,    allocatable :: patch_coords(:,:,:)                                 ! Patch center coordinates
    real,    allocatable :: local_offsets_x(:,:,:), local_offsets_y(:,:,:)      ! Local offsets
    real                 :: rmsd_fit(2)
    ! Polynomial models
    integer              :: fixed_frame                                         ! Index of reference frame
    real(dp)             :: model_coeffs_x(MODELSZ), model_coeffs_y(MODELSZ)    ! Polynomial model coefficients
    ! outliers
    integer, allocatable :: outlier_coords(:,:)
    logical              :: exists = .false.
  contains
    procedure          :: new
    procedure          :: add_patches
    procedure, private :: determine_patch_coordinates
    procedure          :: fit
    procedure          :: set_drift_offsets
    procedure          :: set_local_offsets
    procedure          :: set_frameweights
    procedure          :: set_model_coeffs
    procedure          :: set_outlier_coords
    procedure          :: eval_local_offset
    procedure          :: pix2coords, pix2coordx, pix2coordy
    procedure          :: write
    procedure, private :: write_doc
    procedure, private :: write_bin
    procedure          :: read
    procedure          :: display
    procedure          :: kill
end type motion_model

contains

    subroutine new( self, params, movie, ldim_movie, smpd_movie, scaled_frames,&
            &total_nframes, fixed_frame, voltage, dose_per_frame, eer_fraction, gain )
        class(motion_model), intent(inout) :: self
        class(parameters), target, intent(in) :: params
        class(string),             intent(in) :: movie
        integer,                   intent(in) :: ldim_movie(2)
        real,                      intent(in) :: smpd_movie
        class(image),     pointer, intent(in) :: scaled_frames(:)
        integer,                   intent(in) :: total_nframes
        integer,                   intent(in) :: fixed_frame
        real,                      intent(in) :: voltage
        real,                      intent(in) :: dose_per_frame
        integer,                   intent(in) :: eer_fraction
        class(string),   optional, intent(in) :: gain
        integer :: ldim(3)
        call self%kill
        self%p_ptr => params
        self%movie = simple_abspath(movie)
        if( present(gain) ) self%gain = gain
        self%eer   = fname2format(self%movie) == 'K'
        self%ldim_movie   = ldim_movie
        self%smpd_movie   = smpd_movie
        self%nframes      = size(scaled_frames)
        self%total_nframes = total_nframes
        ldim              = scaled_frames(1)%get_ldim()
        self%ldim         = ldim(1:2)
        self%smpd         = scaled_frames(1)%get_smpd()
        self%binning      = self%smpd / self%smpd_movie
        self%fixed_frame  = fixed_frame
        self%voltage          = voltage
        self%dose_per_frame   = dose_per_frame
        self%target_dose_per_frame = self%p_ptr%fraction_dose_target
        self%dw               = self%dose_per_frame > 0.001
        self%accumulated_dose = real(self%nframes) * self%dose_per_frame
        self%total_dose       = self%p_ptr%total_dose
        self%eer_fraction     = eer_fraction
        if( self%nframes < 0 .or. self%nx_patch < 0 .or. self%ny_patch < 0 ) THROW_HARD('invalid motion model array dimensions')
        allocate(self%drift_offsets_x(self%nframes), self%drift_offsets_y(self%nframes))
        allocate(self%frameweights(self%nframes))
        self%drift_offsets_x = 0.0
        self%drift_offsets_y = 0.0
        self%frameweights    = 1. / real(self%nframes)
        self%exists = .true.
    end subroutine new

    subroutine add_patches( self, nx, ny )
        class(motion_model), intent(inout) :: self
        integer,             intent(in)    :: nx, ny
        self%model_coeffs_x = 0.0d0
        self%model_coeffs_y = 0.0d0
        if( self%nx_patch == nx .and. self%ny_patch == ny ) then
            self%patch_bounds = 0
            self%patch_coords = 0.
            self%local_offsets_x = 0.
            self%local_offsets_y = 0.
        else
            self%nx_patch = nx
            self%ny_patch = ny
            self%npatch   = self%nx_patch * self%ny_patch
            if( allocated(self%patch_bounds) ) deallocate(self%patch_bounds)
            if( allocated(self%patch_coords) ) deallocate(self%patch_coords)
            if( allocated(self%local_offsets_x) ) deallocate(self%local_offsets_x)
            if( allocated(self%local_offsets_y) ) deallocate(self%local_offsets_y)
            allocate(self%patch_bounds(self%nx_patch, self%ny_patch, 2, 2), source=0)
            allocate(self%patch_coords(self%nx_patch, self%ny_patch, 2),&
            &self%local_offsets_x(self%nframes, self%nx_patch, self%ny_patch),&
            &self%local_offsets_y(self%nframes, self%nx_patch, self%ny_patch),source=0.0)
        endif
        call self%determine_patch_coordinates
    end subroutine add_patches

    subroutine determine_patch_coordinates( self )
        class(motion_model), intent(inout) :: self
        integer :: i,j
        real    :: cen, dist
        if( self%nx_patch <= 0 .or. self%ny_patch <= 0 ) then
            self%ldim_patch = self%ldim
            return
        endif
        self%ldim_patch(1) = round2even(real(self%ldim(1)) / real(self%nx_patch))
        self%ldim_patch(2) = round2even(real(self%ldim(2)) / real(self%ny_patch))
        ! fftw friendly size
        self%ldim_patch(1) = find_larger_magic_box(self%ldim_patch(1))
        self%ldim_patch(2) = find_larger_magic_box(self%ldim_patch(2))
        ! along X, limits & center first patches
        self%patch_bounds(1,:,1,1) = 1
        self%patch_bounds(1,:,1,2) = self%ldim_patch(1)
        self%patch_coords(1,:,1)   = real(sum(self%patch_bounds(1,:,1,1:2),dim=2)) / 2.
        ! limits & center last patches
        self%patch_bounds(self%nx_patch,:,1,1) = self%ldim(1)-self%ldim_patch(1)+1
        self%patch_bounds(self%nx_patch,:,1,2) = self%ldim(1)
        self%patch_coords(self%nx_patch,:,1)  = real(sum(self%patch_bounds(self%nx_patch,:,1,1:2),dim=2)) / 2.
        ! adjust other centers for uniform intervals
        dist = real(self%patch_coords(self%nx_patch,1,1)-self%patch_coords(1,1,1)+1) / real(self%nx_patch-1)
        do i=2,self%nx_patch-1
            cen = self%patch_coords(1,1,1) + real(i-1)*dist
            self%patch_bounds(i,:,1,1) = ceiling(cen) - self%ldim_patch(1)/2
            self%patch_bounds(i,:,1,2) = self%patch_bounds(i,:,1,1) + self%ldim_patch(1) - 1
            self%patch_coords(i,:,1)  = real(sum(self%patch_bounds(i,:,1,1:2),dim=2)) / 2.
        enddo
        ! along Y
        self%patch_bounds(:,1,2,1) = 1
        self%patch_bounds(:,1,2,2) = self%ldim_patch(2)
        self%patch_coords(:,1,2)   = real(sum(self%patch_bounds(:,1,2,1:2),dim=2)) / 2.
        self%patch_bounds(:,self%ny_patch,2,1) = self%ldim(2)-self%ldim_patch(2)+1
        self%patch_bounds(:,self%ny_patch,2,2) = self%ldim(2)
        self%patch_coords(:,self%ny_patch,2)  = real(sum(self%patch_bounds(:,self%ny_patch,2,1:2),dim=2)) / 2.
        dist = real(self%patch_coords(1,self%ny_patch,2)-self%patch_coords(1,1,2)+1) / real(self%ny_patch-1)
        do j=2,self%ny_patch-1
            cen = self%patch_coords(1,1,2) + real(j-1)*dist
            self%patch_bounds(:,j,2,1) = ceiling(cen) - self%ldim_patch(2)/2
            self%patch_bounds(:,j,2,2) = self%patch_bounds(:,j,2,1) + self%ldim_patch(2) - 1
            self%patch_coords(:,j,2)  = real(sum(self%patch_bounds(:,j,2,1:2),dim=2)) /2.
        enddo
    end subroutine determine_patch_coordinates

    ! Operations on offsets

    ! Conversion of physical pixel coordinates to normalized coordinates
    pure subroutine pix2coords( self, x, y, cx, cy )
        class(motion_model), intent(in)  :: self
        real(dp),                       intent(in)  :: x, y
        real(dp),                       intent(out) :: cx, cy
        cx = self%pix2coordx( x )
        cy = self%pix2coordy( y )
    end subroutine pix2coords

    elemental real(dp) function pix2coordx( self, x )
        class(motion_model), intent(in)  :: self
        real(dp),                       intent(in)  :: x
        pix2coordx = (x-1.d0) / real(self%ldim(1)-1,dp) - 0.5d0
    end function pix2coordx

    real(dp) elemental function pix2coordy( self, y )
        class(motion_model), intent(in)  :: self
        real(dp),                       intent(in)  :: y
        pix2coordy = (y-1.d0) / real(self%ldim(2)-1,dp) - 0.5d0
    end function pix2coordy

    pure subroutine eval_local_offset( self, iframe, px, py, shift )
        class(motion_model), intent(in)  :: self
        integer,                        intent(in)  :: iframe, px, py
        real,                           intent(out) :: shift(2)
        real(dp) :: t, x, y
        t        = real(iframe-self%fixed_frame, dp)
        x        = self%pix2coordx( real(self%patch_coords(px,py,1), dp) )
        y        = self%pix2coordy( real(self%patch_coords(px,py,2), dp) )
        shift(1) = eval_poly(self%model_coeffs_x(:), x,y,t)
        shift(2) = eval_poly(self%model_coeffs_y(:), x,y,t)
    end subroutine eval_local_offset

    subroutine fit( self )
        class(motion_model), intent(inout) :: self
        real(dp) :: yx(self%nframes*self%npatch)      ! along x
        real(dp) :: yy(self%nframes*self%npatch)      ! along y
        real(dp) :: x(3,self%nframes*self%npatch)     ! x,y,t
        real(dp) :: sig(self%nframes*self%npatch)
        real(dp) :: v(MODELSZ,MODELSZ), w(MODELSZ), chisq
        real     :: fitted_shift(2)
        integer  :: idx, iframe, i, j, dframe
        ! fit
        sig = 1.d0
        idx = 0
        do iframe = 1, self%nframes
            dframe = real(iframe-self%fixed_frame, dp)
            do i = 1, self%nx_patch
                do j = 1, self%ny_patch
                    idx      = idx+1
                    yx(idx)  = real(self%local_offsets_x(iframe,i,j),dp)
                    yy(idx)  = real(self%local_offsets_y(iframe,i,j),dp)
                    x(1,idx) = self%pix2coordx( real(self%patch_coords(i,j,1),dp) )
                    x(2,idx) = self%pix2coordy( real(self%patch_coords(i,j,2),dp) )
                    x(3,idx) = dframe
                end do
            end do
        end do
        call svd_multifit(x, yx, sig, self%model_coeffs_x(:), v, w, chisq, patch_poly)
        call svd_multifit(x, yy, sig, self%model_coeffs_y(:), v, w, chisq, patch_poly)
        ! goodness of fit
        idx = 0
        self%rmsd_fit = 0.
        do iframe = 1,self%nframes
            do i = 1,self%nx_patch
                do j = 1,self%ny_patch
                    idx = idx+1
                    call self%eval_local_offset(iframe, i, j, fitted_shift)
                    self%rmsd_fit(1) = self%rmsd_fit(1) + real((fitted_shift(1)-yx(idx))**2)
                    self%rmsd_fit(2) = self%rmsd_fit(2) + real((fitted_shift(2)-yy(idx))**2)
                end do
            end do
        end do
        self%rmsd_fit = sqrt(self%rmsd_fit/real(self%nframes*self%npatch))
    end subroutine fit

    subroutine set_drift_offsets( self, drift_x, drift_y )
        class(motion_model), intent(inout) :: self
        real,                intent(in)    :: drift_x(:), drift_y(:)
        character(len=:), allocatable :: error_message
        logical :: valid_dimensions
        valid_dimensions = size(drift_x) == self%nframes .and. size(drift_y) == self%nframes
        if( .not.valid_dimensions )then
            error_message = 'motion model drift-offset dimensions do not match: size(drift_x)='//&
                &int2str(size(drift_x))//', size(drift_y)='//int2str(size(drift_y))//&
                &', nframes='//int2str(self%nframes)
            THROW_HARD(error_message)
        endif
        self%drift_offsets_x = drift_x
        self%drift_offsets_y = drift_y
    end subroutine set_drift_offsets

    subroutine set_local_offsets( self, local_x, local_y )
        class(motion_model), intent(inout) :: self
        real,                intent(in)    :: local_x(:,:,:), local_y(:,:,:)
        integer :: expected_shape(3)
        logical :: valid_dimensions
        expected_shape = [self%nframes, self%nx_patch, self%ny_patch]
        valid_dimensions = all(shape(local_x) == expected_shape) .and. all(shape(local_y) == expected_shape)
        if( .not.valid_dimensions ) THROW_HARD('motion model local-offset dimensions do not match')
        self%local_offsets_x = local_x
        self%local_offsets_y = local_y
    end subroutine set_local_offsets

    subroutine set_model_coeffs( self, coeffs_x, coeffs_y )
        class(motion_model), intent(inout) :: self
        real(dp), intent(in) :: coeffs_x(MODELSZ), coeffs_y(MODELSZ)
        self%model_coeffs_x = coeffs_x
        self%model_coeffs_y = coeffs_y
    end subroutine set_model_coeffs

    subroutine set_frameweights( self, weights )
        class(motion_model), intent(inout) :: self
        real, intent(in) :: weights(:)
        if( size(weights) /= self%nframes ) THROW_HARD('motion model frame-weight dimensions do not match')
        if( .not. allocated(self%frameweights) ) allocate(self%frameweights(self%nframes))
        self%frameweights = weights
    end subroutine set_frameweights

    subroutine set_outlier_coords( self, coords )
        class(motion_model),  intent(inout) :: self
        integer, allocatable, intent(in)    :: coords(:,:)
        if( allocated(coords) )then
            if( size(coords,1) /= 2 ) THROW_HARD('motion model outlier coordinates must have two rows')
            if( allocated(self%outlier_coords) ) deallocate(self%outlier_coords)
            allocate(self%outlier_coords, source=coords)
        else
            if( allocated(self%outlier_coords) ) deallocate(self%outlier_coords)
        end if
    end subroutine set_outlier_coords

    ! I/O

    subroutine write( self, star_fname, bin_fname, write_poly )
        class(motion_model), intent(inout) :: self
        class(string),                  intent(in) :: star_fname, bin_fname
        logical,                        intent(in) :: write_poly
        call self%write_doc(star_fname, write_poly)
        call self%write_bin(bin_fname)
    end subroutine write

    subroutine write_doc( self, star_fname, writepoly )
        use simple_starfile_wrappers
        class(motion_model), intent(inout) :: self
        class(string),                  intent(in) :: star_fname
        logical,                        intent(in) :: writepoly
        real(dp),     allocatable :: poly_coeffs(:)
        type(starfile_table_type) :: table
        real(dp)     :: shift(2), scale
        integer      :: i,iframe, ncoeffs, ndeadpixels, motion_model_version
        logical      :: do_scale, l_outliers, l_gain
        l_outliers = allocated(self%outlier_coords)
        l_gain = .not.self%gain%is_blank()
        scale        = 1.d0 / real(self%binning)
        do_scale     = abs(self%binning-1.0) > 1.e-3
        motion_model_version = merge(1, 0, writepoly)   ! 0: stage drift; 1: stage drift + local
        if( motion_model_version == 1 )then
            ncoeffs = 2*MODELSZ
            allocate(poly_coeffs(ncoeffs))
            poly_coeffs(        1:MODELSZ) = self%model_coeffs_x(1:MODELSZ)
            poly_coeffs(MODELSZ+1:ncoeffs) = self%model_coeffs_y(1:MODELSZ)
            if( do_scale ) poly_coeffs = poly_coeffs / scale
        endif
        call starfile_table__new(table)
        call starfile_table__open_ofile(table, star_fname%to_char())
        call starfile_table__addObject(table)
        call starfile_table__setIsList(table, .true.)
        call starfile_table__setname(table, 'general')
        call starfile_table__setValue_int(table,    EMDL_IMAGE_SIZE_X,                   self%ldim_movie(1))
        call starfile_table__setValue_int(table,    EMDL_IMAGE_SIZE_Y,                   self%ldim_movie(2))
        call starfile_table__setValue_int(table,    EMDL_IMAGE_SIZE_Z,                   self%total_nframes)
        call starfile_table__setValue_string(table, EMDL_MICROGRAPH_MOVIE_NAME,          self%movie%to_char())
        if( l_gain )then
            call starfile_table__setValue_string(table, EMDL_MICROGRAPH_GAIN_NAME,       self%gain%to_char())
        endif
        call starfile_table__setValue_double(table, EMDL_MICROGRAPH_BINNING,             real(self%binning,dp))
        call starfile_table__setValue_double(table, EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, real(self%smpd_movie,dp))
        call starfile_table__setValue_double(table, EMDL_MICROGRAPH_DOSE_RATE,           real(self%target_dose_per_frame, dp))
        call starfile_table__setValue_double(table, EMDL_MICROGRAPH_PRE_EXPOSURE,        0.d0)
        call starfile_table__setValue_double(table, EMDL_CTF_VOLTAGE,                    real(self%voltage, dp))
        call starfile_table__setValue_int(table,    EMDL_MICROGRAPH_START_FRAME,         1)
        if( self%eer )then
            call starfile_table__setValue_int(table, EMDL_MICROGRAPH_EER_UPSAMPLING,     self%p_ptr%eer_upsampling)
            call starfile_table__setValue_int(table, EMDL_MICROGRAPH_EER_GROUPING,       self%eer_fraction)
        endif
        call starfile_table__setValue_int(table, EMDL_MICROGRAPH_MOTION_MODEL_VERSION,   motion_model_version)
        ! if( motion_model_version == 1 )then
        !     call starfile_table__setValue_int(table, SMPL_MOVIE_FRAME_ALIGN,             self%fixed_frame)
        ! endif
        call starfile_table__write_ofile(table)
        ! stage drift
        call starfile_table__clear(table)
        call starfile_table__setIsList(table, .false.)
        call starfile_table__setName(table, "global_shift")
        do iframe = 1, self%total_nframes
            call starfile_table__addObject(table)
            call starfile_table__setValue_int(table,    EMDL_MICROGRAPH_FRAME_NUMBER, iframe)
            if( iframe <= self%nframes )then
                shift(1) = real(self%drift_offsets_x(iframe)-self%drift_offsets_x(1), dp)
                shift(2) = real(self%drift_offsets_y(iframe)-self%drift_offsets_y(1), dp)
                shift    = shift / scale
                call starfile_table__setValue_double(table, EMDL_MICROGRAPH_SHIFT_X, shift(1))
                call starfile_table__setValue_double(table, EMDL_MICROGRAPH_SHIFT_Y, shift(2))
                if(allocated(self%frameweights)) then
                    call starfile_table__setValue_double(table, SMPL_MOVIE_FRAME_WEIGHT, real(self%frameweights(iframe),dp))
                endif
            else
                call starfile_table__setValue_double(table, EMDL_MICROGRAPH_SHIFT_X, -9999.0d0)
                call starfile_table__setValue_double(table, EMDL_MICROGRAPH_SHIFT_Y, -9999.0d0)
                call starfile_table__setValue_double(table, SMPL_MOVIE_FRAME_WEIGHT, 0d0)
            endif
        enddo
        call starfile_table__write_ofile(table)
        if( motion_model_version == 1 )then
            ! local motion model
            call starfile_table__clear(table)
            call starfile_table__setIsList(table, .false.)
            call starfile_table__setName(table, "local_motion_model")
            do iframe = 1, ncoeffs
                call starfile_table__addObject(table)
                call starfile_table__setValue_int(table,    EMDL_MICROGRAPH_MOTION_COEFFS_IDX, iframe-1)
                call starfile_table__setValue_double(table, EMDL_MICROGRAPH_MOTION_COEFF,      poly_coeffs(iframe))
            end do
            call starfile_table__write_ofile(table)
        endif
        ! Defects & hot pixels 0-based coordinates
        if( l_outliers )then
            call starfile_table__clear(table)
            call starfile_table__setIsList(table, .false.)
            call starfile_table__setName(table,   'hot_pixels')
            ndeadpixels = size(self%outlier_coords,dim=2)
            do i = 1, ndeadpixels
                call starfile_table__addObject(table)
                call starfile_table__setValue_double(table, EMDL_IMAGE_COORD_X, real(self%outlier_coords(1,i)-1,dp))
                call starfile_table__setValue_double(table, EMDL_IMAGE_COORD_y, real(self%outlier_coords(2,i)-1,dp))
            end do
            call starfile_table__write_ofile(table)
        endif
        !! Patches shifts,  Unused for now
        ! if( writepoly )then
        !     call starfile_table__clear(mc_starfile)
        !     call starfile_table__setIsList(mc_starfile, .false.)
        !     call starfile_table__setName(mc_starfile, "local_shift")
        !     do i = 1, p_ptr%nxpatch
        !         do j = 1, p_ptr%nypatch
        !             do iframe = 1, nframes
        !                 call starfile_table__addObject(mc_starfile)
        !                 call starfile_table__setValue_int(mc_starfile, EMDL_MICROGRAPH_FRAME_NUMBER, iframe)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_IMAGE_COORD_X, patched_centers(i,j,1)/dpscale)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_IMAGE_COORD_Y, patched_centers(i,j,2)/dpscale)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_X, patched_shifts(1,iframe,i,j)/dpscale)
        !                 call starfile_table__setValue_double(mc_starfile, EMDL_MICROGRAPH_SHIFT_Y, patched_shifts(2,iframe,i,j)/dpscale)
        !             enddo
        !         enddo
        !     enddo
        !     call starfile_table__write_ofile(mc_starfile)
        ! endif
        call starfile_table__close_ofile(table)
        call starfile_table__delete(table)
    end subroutine write_doc

    subroutine write_bin( self, bin_fname )
        class(motion_model), intent(inout) :: self
        class(string),                  intent(in) :: bin_fname
        integer(int8) :: flag
        integer       :: funit, ios, stored_file_version, stored_model_version, noutliers
        open(newunit=funit, file=bin_fname%to_char(), access='stream', form='unformatted', &
            &status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('cannot open motion model binary file for writing')
        stored_file_version = FILE_VERSION
        stored_model_version = MODEL_VERSION
        noutliers = 0
        if( allocated(self%outlier_coords) ) noutliers = size(self%outlier_coords, dim=2)
        self%file_name = bin_fname
        write(funit) stored_file_version, stored_model_version
        call write_string(funit, self%movie)
        call write_string(funit, self%gain)
        write(funit) self%smpd_movie, self%ldim_movie, self%smpd, self%ldim, self%binning
        write(funit) self%nframes, self%total_nframes
        write(funit) self%voltage, self%dose_per_frame, self%target_dose_per_frame,&
            &self%total_dose, self%accumulated_dose
        write(funit) self%rmsd_fit, self%fixed_frame
        flag = merge(1_int8, 0_int8, self%dw);  write(funit) flag
        flag = merge(1_int8, 0_int8, self%eer); write(funit) flag
        write(funit) self%eer_fraction
        if( self%nframes > 0 )then
            write(funit) self%drift_offsets_x, self%drift_offsets_y, self%frameweights
            write(funit) self%npatch, self%nx_patch, self%ny_patch, self%ldim_patch
            if( self%npatch > 0 )then
                write(funit) self%patch_bounds, self%patch_coords
                write(funit) self%model_coeffs_x, self%model_coeffs_y
                write(funit) self%local_offsets_x, self%local_offsets_y
            end if
        endif
        write(funit) noutliers
        if( noutliers > 0 ) write(funit) self%outlier_coords
        close(funit)
    end subroutine write_bin

    subroutine read( self, bin_fname, params )
        class(motion_model), intent(inout) :: self
        class(string),                  intent(in)    :: bin_fname
        class(parameters), target,      intent(in)    :: params
        integer(int8) :: flag
        integer       :: funit, ios, stored_file_version, stored_model_version, noutliers
        call self%kill
        open(newunit=funit, file=bin_fname%to_char(), access='stream', form='unformatted', &
            &status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('cannot open motion model binary file for reading')
        read(funit, iostat=ios) stored_file_version, stored_model_version
        if( ios /= 0 ) THROW_HARD('invalid motion model binary header')
        if( stored_file_version /= FILE_VERSION .or. stored_model_version /= MODEL_VERSION ) THROW_HARD('unsupported motion model binary version')
        call read_string(funit, self%movie)
        call read_string(funit, self%gain)
        read(funit) self%smpd_movie, self%ldim_movie, self%smpd, self%ldim, self%binning
        read(funit) self%nframes, self%total_nframes
        read(funit) self%voltage, self%dose_per_frame, self%target_dose_per_frame,&
            &self%total_dose, self%accumulated_dose
        read(funit) self%rmsd_fit, self%fixed_frame
        read(funit) flag; self%dw  = flag /= 0_int8
        read(funit) flag; self%eer = flag /= 0_int8
        read(funit) self%eer_fraction
        if( self%nframes < 0 .or. self%total_nframes < self%nframes ) THROW_HARD('invalid motion model frame dimensions')
        if( allocated(self%drift_offsets_x) ) deallocate(self%drift_offsets_x)
        if( allocated(self%drift_offsets_y) ) deallocate(self%drift_offsets_y)
        if( allocated(self%frameweights) ) deallocate(self%frameweights)
        if( allocated(self%patch_bounds) ) deallocate(self%patch_bounds)
        if( allocated(self%patch_coords) ) deallocate(self%patch_coords)
        if( allocated(self%local_offsets_x) ) deallocate(self%local_offsets_x)
        if( allocated(self%local_offsets_y) ) deallocate(self%local_offsets_y)
        if( allocated(self%outlier_coords) ) deallocate(self%outlier_coords)
        self%npatch         = 0
        self%nx_patch       = 0
        self%ny_patch       = 0
        self%ldim_patch     = 0
        self%model_coeffs_x = 0.0d0
        self%model_coeffs_y = 0.0d0
        if( self%nframes > 0 )then
            allocate(self%drift_offsets_x(self%nframes), self%drift_offsets_y(self%nframes))
            allocate(self%frameweights(self%nframes))
            read(funit) self%drift_offsets_x, self%drift_offsets_y, self%frameweights
            read(funit) self%npatch, self%nx_patch, self%ny_patch, self%ldim_patch
            if( self%npatch < 0 .or. self%nx_patch < 0 .or. self%ny_patch < 0 )then
                THROW_HARD('invalid motion model patch dimensions')
            endif
            if( self%npatch /= self%nx_patch * self%ny_patch ) THROW_HARD('inconsistent motion model patch dimensions')
            if( self%npatch > 0 )then
                allocate(self%patch_bounds(self%nx_patch, self%ny_patch, 2, 2))
                allocate(self%patch_coords(self%nx_patch, self%ny_patch, 2))
                allocate(self%local_offsets_x(self%nframes, self%nx_patch, self%ny_patch))
                allocate(self%local_offsets_y(self%nframes, self%nx_patch, self%ny_patch))
                read(funit) self%patch_bounds, self%patch_coords
                read(funit) self%model_coeffs_x, self%model_coeffs_y
                read(funit) self%local_offsets_x, self%local_offsets_y
            endif
        endif
        read(funit) noutliers
        if( noutliers < 0 ) THROW_HARD('invalid motion model outlier dimensions')
        if( noutliers > 0 ) allocate(self%outlier_coords(2, noutliers))
        if( noutliers > 0 ) read(funit) self%outlier_coords
        close(funit)
        self%file_name = bin_fname
        self%p_ptr => params
        self%exists = .true.
    end subroutine read

    subroutine display( self )
        class(motion_model), intent(in) :: self
        integer :: i, j
        write(logfhandle,'(a)') '>>> MOTION MODEL:'
        write(logfhandle,'(a25,1x,i0)')      'model size:', MODELSZ
        write(logfhandle,'(a25,1x,l1)')      'parameters associated:', associated(self%p_ptr)
        write(logfhandle,'(a25,1x,a)')       'file name:', trim(self%file_name%to_char())
        write(logfhandle,'(a25,1x,a)')       'movie:', trim(self%movie%to_char())
        write(logfhandle,'(a25,1x,a)')       'gain:', trim(self%gain%to_char())
        write(logfhandle,'(a25,1x,es15.7)')  'movie sampling:', self%smpd_movie
        write(logfhandle,'(a25,2(1x,i0))')   'movie dimensions:', self%ldim_movie
        write(logfhandle,'(a25,1x,es15.7)')  'sampling:', self%smpd
        write(logfhandle,'(a25,2(1x,i0))')   'dimensions:', self%ldim
        write(logfhandle,'(a25,1x,es15.7)')  'binning:', self%binning
        write(logfhandle,'(a25,1x,i0)')      'frames:', self%nframes
        write(logfhandle,'(a25,1x,i0)')      'total frames:', self%total_nframes
        write(logfhandle,'(a25,1x,es15.7)')  'voltage:', self%voltage
        write(logfhandle,'(a25,1x,es15.7)')  'dose per frame:', self%dose_per_frame
        write(logfhandle,'(a25,1x,es15.7)')  'total dose:', self%total_dose
        write(logfhandle,'(a25,1x,es15.7)')  'accumulated dose:', self%accumulated_dose
        write(logfhandle,'(a25,1x,l1)')      'dose weighted:', self%dw
        write(logfhandle,'(a25,1x,l1)')      'EER:', self%eer
        write(logfhandle,'(a25,1x,i0)')      'EER fraction:', self%eer_fraction
        write(logfhandle,'(a25,1x,i0)')      'patches:', self%npatch
        write(logfhandle,'(a25,1x,i0)')      'X patches:', self%nx_patch
        write(logfhandle,'(a25,1x,i0)')      'Y patches:', self%ny_patch
        write(logfhandle,'(a25,2(1x,i0))')   'patch dimensions:', self%ldim_patch
        write(logfhandle,'(a25,2(1x,es15.7))') 'fit RMSD:', self%rmsd_fit
        write(logfhandle,'(a25,1x,i0)')      'fixed frame:', self%fixed_frame
        if( allocated(self%drift_offsets_x) .and. allocated(self%drift_offsets_y) .and.&
            &allocated(self%frameweights) )then
            write(logfhandle,'(a)') '  frame  drift_offsets_x  drift_offsets_y  frameweight'
            do i = 1,self%nframes
                write(logfhandle,'(4x,i0,3(2x,es15.7))') i, self%drift_offsets_x(i),&
                    &self%drift_offsets_y(i), self%frameweights(i)
            enddo
        else
            write(logfhandle,'(a,3(1x,l1))') '  frame arrays allocated (X, Y, weight):',&
                &allocated(self%drift_offsets_x), allocated(self%drift_offsets_y), allocated(self%frameweights)
        endif
        if( allocated(self%patch_bounds) .and. allocated(self%patch_coords) )then
            write(logfhandle,'(a,4(1x,i0),a,3(1x,i0),a)') '  patch shapes (bounds:', shape(self%patch_bounds),&
                &'; coords:', shape(self%patch_coords), ')'
            write(logfhandle,'(a)') '  patch X Y  bounds (X1 X2 Y1 Y2)  coords (X Y)'
            do j = 1,size(self%patch_bounds,2)
                do i = 1,size(self%patch_bounds,1)
                    write(logfhandle,'(4x,2(i0,1x),4(i0,1x),2(es15.7,1x))') i, j,&
                        &self%patch_bounds(i,j,1,1), self%patch_bounds(i,j,1,2),&
                        &self%patch_bounds(i,j,2,1), self%patch_bounds(i,j,2,2),&
                        &self%patch_coords(i,j,1), self%patch_coords(i,j,2)
                enddo
            enddo
        else
            write(logfhandle,'(a,2(1x,l1))') '  patch arrays allocated (bounds, coords):',&
                &allocated(self%patch_bounds), allocated(self%patch_coords)
        endif
        write(logfhandle,'(a)') '  model coefficients (X, Y):'
        do i = 1,MODELSZ
            write(logfhandle,'(4x,i0,2(2x,es23.15))') i, self%model_coeffs_x(i), self%model_coeffs_y(i)
        enddo
    end subroutine display

    ! Destructor

    subroutine kill( self )
        class(motion_model), intent(inout) :: self
        if( self%exists )then
            if( allocated(self%drift_offsets_x) ) deallocate(self%drift_offsets_x)
            if( allocated(self%drift_offsets_y) ) deallocate(self%drift_offsets_y)
            if( allocated(self%frameweights) ) deallocate(self%frameweights)
            if( allocated(self%local_offsets_x) ) deallocate(self%local_offsets_x)
            if( allocated(self%local_offsets_y) ) deallocate(self%local_offsets_y)
            if( allocated(self%patch_coords) ) deallocate(self%patch_coords)
            if( allocated(self%patch_bounds) ) deallocate(self%patch_bounds)
            self%model_coeffs_x = 0.0d0; self%model_coeffs_y = 0.0d0
            call self%file_name%kill
            call self%movie%kill
            call self%gain%kill
            self%smpd_movie       = 0.
            self%ldim_movie       = 0
            self%smpd             = 0.
            self%ldim             = 0
            self%binning          = 0.
            self%nframes          = 0
            self%total_nframes    = 0
            self%dose_per_frame   = 0.
            self%target_dose_per_frame = 0.
            self%total_dose       = 0.
            self%accumulated_dose = 0.
            self%npatch           = 0
            self%nx_patch         = 0
            self%ny_patch         = 0
            self%ldim_patch       = 0
            self%fixed_frame      = 0
            self%rmsd_fit         = 0.
            if( allocated(self%outlier_coords) ) deallocate(self%outlier_coords)
            nullify(self%p_ptr)
        endif
        self%exists = .false.
    end subroutine kill

    ! PRIVATE TYPE UNBOUND HELPERS

    pure function patch_poly(p, n) result(res)
        real(dp), intent(in) :: p(:)
        integer,  intent(in) :: n
        real(dp) :: res(n), x, y, t
        x = p(1); y = p(2); t = p(3)
        res(    1) = t
        res(    2) = t**2
        res(    3) = t**3
        res( 4: 6) = x * res( 1: 3)  ! x   * {t,t^2,t^3}
        res( 7: 9) = x * res( 4: 6)  ! x^2 * {t,t^2,t^3}
        res(10:12) = y * res( 1: 3)  ! y   * {t,t^2,t^3}
        res(13:15) = y * res(10:12)  ! y^2 * {t,t^2,t^3}
        res(16:18) = y * res( 4: 6)  ! x*y * {t,t^2,t^3}
    end function patch_poly

    ! evaluate the patch polynomial at (x, y, t)
    pure real(sp) function eval_poly(C, x, y, t)
        real(dp), intent(in) :: C(MODELSZ), x, y, t
        real(dp) :: res, x2, y2, xy, t2, t3
        xy = x * y
        x2 = x * x; y2 = y * y
        t2 = t * t
        t3 = t2 * t
        res =       C( 1) * t      + C( 2) * t2      + C( 3) * t3
        res = res + C( 4) * t * x  + C( 5) * t2 * x  + C( 6) * t3 * x
        res = res + C( 7) * t * x2 + C( 8) * t2 * x2 + C( 9) * t3 * x2
        res = res + C(10) * t * y  + C(11) * t2 * y  + C(12) * t3 * y
        res = res + C(13) * t * y2 + C(14) * t2 * y2 + C(15) * t3 * y2
        res = res + C(16) * t * xy + C(17) * t2 * xy + C(18) * t3 * xy
        eval_poly = real(res, sp)
    end function eval_poly

    subroutine write_string(funit, value)
        integer, intent(in) :: funit
        class(string), intent(in) :: value
        character(len=:), allocatable :: text
        integer(int32) :: n
        text = value%to_char()
        n = len(text)
        write(funit) n
        if( n > 0 ) write(funit) text
    end subroutine write_string

    subroutine read_string(funit, value)
        integer, intent(in) :: funit
        class(string), intent(inout) :: value
        character(len=:), allocatable :: text
        integer(int32) :: n
        read(funit) n
        if( n < 0 ) THROW_HARD('invalid motion model string length')
        if( n > 0 )then
            allocate(character(len=n) :: text)
            read(funit) text
            value = text
        else
            value = ''
        endif
    end subroutine read_string

end module simple_motion_model
