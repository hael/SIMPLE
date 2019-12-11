! Fortran wrapper for libtiff, edited from Unblur
! Copyright 2014 Howard Hughes Medical Institute
! All rights reserved
! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
module simple_tifflib
#include "simple_local_flags.inc"
use, intrinsic :: iso_c_binding
implicit none

#ifdef USING_TIFF
interface

    function TIFFGetField(tiff, flag, val) bind(C,name="TIFFGetField")
        use iso_c_binding, only : c_int16_t, c_ptr
        integer(kind=c_int16_t)     ::  TIFFGetField
        type(c_ptr),        value   ::  tiff
        integer(kind=c_int16_t)     ::  flag, val
    end function  TIFFGetField

    ! TIFF* TIFFOpen(const char *filename, const char *mode)
    function TIFFOpen(cfilename, cmode) bind(c,name='TIFFOpen')
        use iso_c_binding, only : c_ptr, c_char
        type(c_ptr) ::  TIFFOpen
        character(kind=c_char)  ::  cfilename(*)
        character(kind=c_char)  ::  cmode(*)
    end function TIFFOpen

    ! tsize_t TIFFReadEncodedStrip(TIFF *tif, tstrip_t strip, tdata_t buf, tsize_t size)
    ! tstrip_t is uint32
    ! tdata_t is void*
    ! tsize_t is int32
    function TIFFReadEncodedStrip(tiff,strip,buf,my_size) bind(c,name='TIFFReadEncodedStrip')
        use iso_c_binding, only : c_int32_t, c_ptr
        integer(kind=c_int32_t)         ::  TIFFReadEncodedStrip
        type(c_ptr),            value   ::  tiff
        integer(kind=c_int32_t),value   ::  strip
        type(c_ptr),            value   ::  buf
        integer(kind=c_int32_t)         ::  my_size
    end function

    function TIFFLastDirectory(tiff) bind(c,name='TIFFLastDirectory')
        use iso_c_binding, only : c_int, c_ptr
        integer(kind=c_int)         ::  TIFFLastDirectory
        type(c_ptr),        value   ::  tiff
    end function

    function TIFFIsTiled(tiff) bind(c,name='TIFFIsTiled')
        use iso_c_binding, only : c_int, c_ptr
        integer(kind=c_int)         ::  TIFFIsTiled
        type(c_ptr),        value   ::  tiff
    end function

    ! tstrip_t TIFFNumberOfStrips(TIFF *tif)
    ! tstrip_t is uint32
    function TIFFNumberOfStrips(tiff) bind(c,name='TIFFNumberOfStrips')
        use iso_c_binding, only : c_ptr, c_int
        integer(kind=c_int)         ::  TIFFNumberOfStrips
        type(c_ptr),        value   ::  tiff
    end function

    ! tsize_t TIFFStripSize(TIFF *tif)
    ! tsize_t is int32
    function TIFFStripSize(tiff) bind(c,name='TIFFStripSize')
        use iso_c_binding, only : c_ptr, c_int32_t
        integer(kind=c_int32_t)     ::  TIFFStripSize
        type(c_ptr),        value   ::  tiff
    end function

    ! tdir_t (uint16)
    function TIFFCurrentDirectory(tiff) bind(c,name='TIFFCurrentDirectory')
        use iso_c_binding, only : c_int16_t, c_ptr
        integer(kind=c_int16_t)     ::  TIFFCurrentDirectory
        type(c_ptr),        value   ::  tiff
    end function

    ! int TIFFSetDirectory(TIFF *tif, tdir_t dirnum)
    ! tdir_t is uint16
    function TIFFSetDirectory(tiff,dirnum) bind(c,name='TIFFSetDirectory')
        use iso_c_binding, only : c_ptr, c_int, c_int16_t
        integer(kind=c_int)                 ::  TIFFSetDirectory
        type(c_ptr),        value           ::  tiff
        integer(kind=c_int16_t), value      ::  dirnum
    end function

    ! int TIFFReadDirectory(TIFF *tif)
    !
    function TIFFReadDirectory(tiff) bind(c,name='TIFFReadDirectory')
        use iso_c_binding, only : c_ptr, c_int
        integer(kind=c_int)         ::  TIFFReadDirectory
        type(c_ptr),        value   ::  tiff
    end function

    !void TIFFClose(TIFF *tif)
    subroutine TIFFClose(tiff) bind(c,name='TIFFClose')
        use iso_c_binding, only : c_ptr
        type(c_ptr),        value   ::  tiff
    end subroutine

    !
    ! These functions are defined in the C helper functions
    !
    integer function TIFFGetWidth(my_tiff) bind(c,name='TIFFGetWidth')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    integer function TIFFGetLength(my_tiff) bind(c,name='TIFFGetLength')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    ! int TIFFGetSampleFormat(TIFF *my_tiff)
    integer function TIFFGetSampleFormat(my_tiff) bind(c,name='TIFFGetSampleFormat')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    ! int TIFFGetSamplesPerPixel(TIFF *my_tiff)
    integer function TIFFGetSamplesPerPixel(my_tiff) bind(c,name='TIFFGetSamplesPerPixel')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    ! int TIFFGetBitsPerSample(TIFF *my_tiff)
    integer function TIFFGetBitsPerSample(my_tiff) bind(c,name='TIFFGetBitsPerSample')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    ! intTIFFGetRowsPerStrip(TIFF *my_tiff)
    integer function TIFFGetRowsPerStrip(my_tiff) bind(c,name='TIFFGetRowsPerStrip')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    ! intTIFFGetMinVal(TIFF *my_tiff)
    integer function TIFFGetMinVal(my_tiff) bind(c,name='TIFFGetMinVal')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    ! intTIFFGetMaxVal(TIFF *my_tiff)
    integer function TIFFGetMaxVal(my_tiff) bind(c,name='TIFFGetMaxVal')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value :: my_tiff
    end function
    ! tdata_t TIFFAllocateStripBuffer(TIFF *my_tiff)
    function TIFFAllocateStripBuffer(tiff) bind(c,name='TIFFAllocateStripBuffer')
        use iso_c_binding, only : c_ptr
        type(c_ptr)             ::  TIFFAllocateStripBuffer
        type(c_ptr),    value   ::  tiff
    end function

    ! int TIFFDeallocateStripBuffer(tdata_t buf)
    function TIFFDeallocateStripBuffer(buf) bind(c,name='TIFFDeallocateStripBuffer')
        use iso_c_binding, only : c_ptr, c_int
        type(c_ptr),    value   ::  buf
        integer(kind=c_int)     ::  TIFFDeallocateStripBuffer
    end function

    subroutine TIFFPrintInfo(tiff) bind(c,name='TIFFPrintInfo')
        use iso_c_binding, only : c_ptr
        type(c_ptr),    value   ::  tiff
    end subroutine

end interface
#endif
end module simple_tifflib
