! OpenSlide
! The library can read virtual slides in the following formats:

! Aperio (.svs, .tif)
! Hamamatsu (.vms, .vmu, .ndpi)
! Leica (.scn)
! MIRAX (.mrxs)
! Philips (.tiff)
! Sakura (.svslide)
! Trestle (.tif)
! Ventana (.bif, .tif)
! Generic tiled TIFF (.tif)

! requires -lopenslide -ltiff -lm

module simple_openslide
use, intrinsic :: iso_c_binding
include 'simple_lib.f08'
implicit none
private
#include "simple_local_flags.inc"
public :: read_openslide
public :: write_openslide
! typedef struct _openslide openslide_t;
type(C_PTR) :: openslide_t

interface
    ! openslide_t *openslide_open(const char *filename);
    type(C_PTR) function openslide_open(filename) bind(C,name="openslide_open")
        import
        character(len=1,kind=C_CHAR), dimension(*), intent(in) :: filename
    end function openslide_open
    integer(C_int32_t) function openslide_get_level_count(osr) bind(C,name="openslide_get_level_count")
        import
        type(C_PTR) :: osr
    end function openslide_get_level_count
    !void openslide_get_level_dimensions(openslide_t *osr, int32_t level, int64_t *w, int64_t *h);
    subroutine openslide_get_level0_dimensions(osr, w, h) bind(C,name="openslide_get_level_dimensions")
        import
        type(C_PTR) :: osr
        integer(C_int64_t), intent(inout) :: w,h
    end subroutine openslide_get_level0_dimensions

    !void openslide_get_level_dimensions(openslide_t *osr, int32_t level, int64_t *w, int64_t *h);
    subroutine openslide_get_level_dimensions(osr, level, w, h) bind(C,name="openslide_get_level_dimensions")
        import
        type(C_PTR) :: osr
        integer(C_int32_t), intent(inout) :: level
        integer(C_int64_t), intent(inout) :: w,h
    end subroutine openslide_get_level_dimensions

    ! double openslide_get_level_downsample(openslide_t *osr, int32_t level);
    real (C_DOUBLE) function openslide_get_level_downsample(osr, level) bind(C,name="openslide_get_best_level_downsample")
        import
        type(C_PTR) :: osr
        integer(C_int32_t), intent(inout) :: level
    end function openslide_get_level_downsample

    ! int32_t openslide_get_best_level_for_downsample(openslide_t *osr,
    ! double downsample);
    integer(C_int32_t) function openslide_get_best_level_for_downsample( osr, downsample) &
        bind(C, name="openslide_get_best_level_for_downsample")
        import
        type(C_PTR) :: osr
        real(C_DOUBLE), intent(inout) :: downsample
    end function openslide_get_best_level_for_downsample
    ! void openslide_read_region(openslide_t *osr,
    !                uint32_t *dest,
    !                int64_t x, int64_t y,
    !                int32_t level,
    !       int64_t w, int64_t h);
    subroutine openslide_read_region(osr, dest, x, y, level, w, h)  bind(C, name="openslide_read_region")
        import
        type(C_PTR) :: osr
        type(C_PTR), intent(inout) :: dest
        integer(C_int64_t), intent(in) :: x,y
        integer(C_int32_t), intent(in) :: level
        integer(C_int64_t), intent(in) :: w,h
    end subroutine openslide_read_region

    !>    void openslide_close(openslide_t *osr);
    subroutine openslide_close(osr) bind(C,name="openslide_close")
        import
        type(C_PTR) :: osr
    end subroutine openslide_close
    !  const char *openslide_get_error(openslide_t *osr);
    function openslide_get_error(osr) bind(C,name="openslide_get_error")
        import
        type(C_PTR)  :: openslide_get_error ! const char *
        type(C_PTR) :: osr
    end function openslide_get_error
    ! const char * const *openslide_get_property_names(openslide_t *osr);
    function openslide_get_property_names(osr) bind(C,name="openslide_get_property_names")
        import
        type(C_PTR) :: openslide_get_property_names ! const char *
        type(C_PTR) :: osr
    end function openslide_get_property_names
    ! const char *openslide_get_property_value(openslide_t *osr, const char *name);
    function openslide_get_property_value(osr,name) bind(C,name="openslide_get_property_value")
        import
       type(C_PTR)  :: openslide_get_property_value ! const char *
        type(C_PTR) :: osr
        character(len=1,kind=C_CHAR), dimension(*), intent(in) :: name
    end function openslide_get_property_value
    ! const char * const *openslide_get_associated_image_names(openslide_t *osr);
    function openslide_get_associated_image_names(osr) bind(C,name="openslide_get_associated_image_names")
        import
        type(C_PTR) :: openslide_get_associated_image_names ! const char *
        type(C_PTR) :: osr
    end function openslide_get_associated_image_names
    ! void openslide_get_associated_image_dimensions(openslide_t *osr,
    !                            const char *name,
    !                            int64_t *w, int64_t *h);
    subroutine openslide_get_associated_image_dimensions(osr, name, w, h) bind(C,name="openslide_get_associated_image_dimensions")
        import
        type(C_PTR) :: osr
        character(len=1,kind=C_CHAR), dimension(*), intent(in) :: name
        integer(C_int64_t), intent(out) :: w,h
    end subroutine openslide_get_associated_image_dimensions
    ! void openslide_read_associated_image(openslide_t *osr,
    !                      const char *name,
    !                      uint32_t *dest);
    subroutine openslide_read_associated_image(osr, name, dest) bind(C,name="openslide_read_associated_image")
        import
        type(C_PTR) :: osr
        character(len=1,kind=C_CHAR), dimension(*), intent(in) :: name
        type(C_PTR), intent(out) :: dest
    end subroutine openslide_read_associated_image
    ! const char *openslide_get_version(void);
    function openslide_get_version() bind(C,name="openslide_get_version")
        import
        type(C_PTR)  :: openslide_get_version ! const char *
    end function openslide_get_version
    !! DEPRECIATED
    ! bool openslide_can_open(const char *filename);
    integer(C_BOOL) function openslide_can_open(filename) bind(C,name="openslide_can_open")
        import
        character(len=1,kind=C_CHAR), dimension(*), intent(in) :: filename
    end function openslide_can_open
    ! int32_t openslide_get_layer_count(openslide_t *osr);
    integer(C_int32_t) function openslide_get_layer_count(osr) bind(C,name="openslide_get_level_count")
        import
        type(C_PTR) :: osr
    end function openslide_get_layer_count

end interface

contains

    subroutine read_openslide(filename, img)
        character(len=*),           intent(in) :: filename   !< target pathname
        real, allocatable, target, intent(out) :: img(:,:)
        type(c_ptr)                          :: osr
        character(kind=c_char,len=STDLEN)    :: infilename_c
        integer(C_INT64_T) :: w, h
        type(C_PTR) :: destptr
        integer(C_int32_t), pointer :: dest(:)
        integer(8) :: i,j
        real :: mx
        infilename_c=trim(filename)//C_NULL_CHAR
        osr = openslide_open(trim(infilename_c))
        if( c_associated(osr) )then
            call openslide_get_level0_dimensions(osr, w, h)
            call openslide_read_region(osr, destptr, 0_8, 0_8, 0, w, h)
            ! call c_f_pointer(dest, destptr, shape=(w*h))
            allocate(img(w,h), source=0.)
            do j=1, h
                do i=1,w
                    img(i,j) = REAL(dest(j*w + i))
                end do
            end do
            mx = maxval(img)
            img = img/mx
        else
            write(*,'(a)') 'OpenSlide read_openslide failed to open '//trim(filename)
            stop
        endif

    end subroutine read_openslide

    subroutine write_openslide(filename, img)
        character(len=*),           intent(in)  :: filename
        real, allocatable, intent(in) :: img(:,:)
        character(kind=c_char,len=STDLEN)    :: outfilename_c




    end subroutine write_openslide


end module simple_openslide
