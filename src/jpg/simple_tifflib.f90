!! Fortran wrapper for libtiff
module simple_tifflib
include 'simple_lib.f08'
use, intrinsic :: iso_c_binding
implicit none
private
#include "simple_local_flags.inc"
public :: write_tiff
public :: write_bigtiff
public :: write_tiff2
public :: write_tiff3
public :: write_tiff_bigimg

interface
    integer(C_INT) function bigtiff_write (fname, nx, ny, nz, min_val, max_val,img) bind(C, name="bigtiff_write")
        use,intrinsic :: iso_c_binding
        implicit none
        character(len=1, kind=c_char), dimension(*), intent(in) :: fname
        integer(C_INT), intent(in) :: nx,ny,nz
        real(C_FLOAT), intent(in) :: min_val, max_val
        type(C_PTR), intent(in) :: img
    end function bigtiff_write
end interface

    ! interface
    !     subroutine tileContigRoutine(timg, a,b, c, d, e, f, h) bind(C,name="tileContigRoutine")
    !         import
    !         type(TIFFRGBAImage),pointer:: timg
    !         integer(C_INT32_T), dimension(*) :: a   !uint32*
    !         integer(C_INT32_T) :: b, c, d, e, f, g  ! uint32 *4, int32 *2
    !         uint32*, uint32, uint32, uint32, uint32, int32, int32,
    !         character(kind=c_char, len=1), dimension(*) :: h  !unsigned char*
    !     end subroutine tileContigRoutine
    ! end interface

    ! type, bind(C) :: TIFFput
    !     ! put decoded strip/tile
    !     type(C_FUNPTR) :: any   !void (*any)(TIFFRGBAImage*);
    !     type(C_FUNPTR) :: contig  ! tileContigRoutine
    !     type(C_FUNPTR) :: separate ! tileSeparateRoutine
    ! end type TIFFput


    ! type,bind(C) :: TIFFRGBAImage
    !     type(C_PTR) ::  tif                              !!  image handle
    !     integer(c_int)  :: stoponerr                          !!  stop on read error
    !     integer(c_int) :: isContig                           !!  data is packed/separate
    !     integer(c_int) ::  alpha                              !!  type of alpha data present
    !     integer(c_int) :: width                           !!  image width
    !     integer(c_int) :: height                          !!  image height
    !     integer(c_int16_t) :: bitspersample                   !!  image bits/sample
    !     integer(c_int16_t) :: samplesperpixel                 !!  image samples/pixel
    !     integer(c_int16_t) :: orientation                     !!  image orientation
    !     integer(c_int16_t) :: req_orientation                 !!  requested orientation
    !     integer(c_int16_t) :: photometric                     !!  image photometric interp
    !     integer(c_int16_t), allocatable :: redcmap(:)                        !!  colormap pallete
    !     integer(c_int16_t), allocatable :: greencmap(:)
    !     integer(c_int16_t), allocatable :: bluecmap(:)
    !     ! get image data routine
    !     integer(C_FUNPTR) :: get  ! int (*get)(TIFFRGBAImage*, uint32*, uint32, uint32)
    !     type(TIFFput) :: put
    !     !!   TIFFRGBValue* Map                      !!  sample mapping array
    !     integer(C_INT), allocatable :: BWmap(:,:)                         !!  black&white map
    !     integer(C_INT), allocatable :: PALmap(:,:)                        !!  palette image map
    !     !!    TIFFYCbCrToRGB* ycbcr                  !!  YCbCr conversion state
    !     !!    TIFFCIELabToRGB* cielab                !!  CIE L*a*b conversion state

    !     integer(C_INT8_T), allocatable ::  UaToAa(:)                          !!  Unassociated alpha to associated alpha convertion LUT
    !     integer(C_INT8_T), allocatable ::  Bitdepth16To8 (:)                  !!  LUT for conversion from 16bit to 8bit values

    !     integer(C_INT)  :: row_offset
    !     integer(C_INT) ::  col_offset

    ! end type TIFFRGBAImage
    !! tiff.h
    integer(C_INT), parameter :: TIFF_VERSION_CLASSIC = 42
    integer(C_INT), parameter :: TIFF_VERSION_BIG     = 43
    integer(C_INT), parameter :: TIFF_BIGENDIAN    =  INT(Z'4d4d')
    integer(C_INT), parameter :: TIFF_LITTLEENDIAN =  INT(Z'4949')
    enum, bind(C)  !  TIFFDataType
        enumerator :: TIFF_NOTYPE = 0      !!  placeholder
        enumerator :: TIFF_BYTE = 1        !!  8-bit unsigned integer
        enumerator :: TIFF_ASCII = 2       !!  8-bit bytes w/ last byte null
        enumerator :: TIFF_SHORT = 3       !!  16-bit unsigned integer
        enumerator :: TIFF_LONG = 4        !!  32-bit unsigned integer
        enumerator :: TIFF_RATIONAL = 5   !!  64-bit unsigned fraction
        enumerator :: TIFF_LONG8 = 16      !!  BigTIFF 64-bit unsigned integer
        enumerator :: TIFF_SLONG8 = 17     !!  BigTIFF 64-bit signed integer
        enumerator :: TIFF_IFD8 = 18       !!  BigTIFF 64-bit unsigned integer (offset)
    end enum

    !!
    !! TIFF Tag Definitions.

    integer(C_INT), parameter :: TIFFTAG_SUBFILETYPE            = 254   ! subfile data descriptor
    integer(C_INT), parameter ::   FILETYPE_REDUCEDIMAGE        = INT(Z'1')     ! reduced resolution version
    integer(C_INT), parameter ::   FILETYPE_PAGE                = INT(Z'2')     ! one page of many
    integer(C_INT), parameter ::   FILETYPE_MASK                = INT(Z'4')     ! transparency mask
    integer(C_INT), parameter :: TIFFTAG_OSUBFILETYPE           = 255   ! +kind of data in subfile
    integer(C_INT), parameter ::   OFILETYPE_IMAGE              = 1     ! full resolution image data
    integer(C_INT), parameter ::   OFILETYPE_REDUCEDIMAGE       = 2     ! reduced size image data
    integer(C_INT), parameter ::   OFILETYPE_PAGE               = 3     ! one page of many
#define TIFFTAG_IMAGEWIDTH 256
    ! image width in pixels
#define TIFFTAG_IMAGELENGTH 257
    ! image height in pixels
#define TIFFTAG_BITSPERSAMPLE 258
    ! bits per channel (sample) */
#define TIFFTAG_COMPRESSION 259
    ! data compression technique */
   ! integer(C_INT), parameter :: TIFFTAG_IMAGEWIDTH             = 256   ! image width in pixels
   ! integer(C_INT), parameter :: TIFFTAG_IMAGELENGTH            = 257   ! image height in pixels
   ! integer(C_INT), parameter :: TIFFTAG_BITSPERSAMPLE          = 258   ! bits per channel (sample)
   ! integer(C_INT), parameter :: TIFFTAG_COMPRESSION            = 259   ! data compression technique
#define COMPRESSION_NONE 1
! integer(C_INT), parameter ::   COMPRESSION_NONE             = 1     ! dump mode
    integer(C_INT), parameter ::   COMPRESSION_CCITTRLE         = 2     ! CCITT modified Huffman RLE
    integer(C_INT), parameter ::   COMPRESSION_CCITTFAX3        = 3     ! CCITT Group 3 fax encoding
    integer(C_INT), parameter ::   COMPRESSION_CCITT_T4         = 3     ! CCITT T.4 ( TIFF 6 name)
    integer(C_INT), parameter ::   COMPRESSION_CCITTFAX4        = 4     ! CCITT Group 4 fax encoding
    integer(C_INT), parameter ::   COMPRESSION_CCITT_T6         = 4     ! CCITT T.6 ( TIFF 6 name)
    integer(C_INT), parameter ::   COMPRESSION_LZW              = 5     ! Lempel-Ziv  & Welch
    integer(C_INT), parameter ::   COMPRESSION_OJPEG            = 6     ! 6.0 JPEG
    integer(C_INT), parameter ::   COMPRESSION_JPEG             = 7     ! %JPEG DCT compression
    integer(C_INT), parameter ::   COMPRESSION_T85              = 9     ! ! TIFF/FX T.85 JBIG compression
    integer(C_INT), parameter ::   COMPRESSION_T43              = 10    ! ! TIFF/FX T.43 colour by layered JBIG compression
    integer(C_INT), parameter ::   COMPRESSION_NEXT             = 32766 ! NeXT 2-bit RLE
    integer(C_INT), parameter ::   COMPRESSION_CCITTRLEW        = 32771 ! #1 w/ word alignment
    integer(C_INT), parameter ::   COMPRESSION_PACKBITS         = 32773 ! Macintosh RLE
    integer(C_INT), parameter ::   COMPRESSION_THUNDERSCAN      = 32809 ! ThunderScan RLE
    ! codes 32895-32898 are reserved for ANSI IT8                                    TIFF/IT <dkelly@apago.com)
    integer(C_INT), parameter ::   COMPRESSION_IT8CTPAD         = 32895 ! IT8 CT w/padding
    integer(C_INT), parameter ::   COMPRESSION_IT8LW            = 32896 ! IT8 Linework RLE
    integer(C_INT), parameter ::   COMPRESSION_IT8MP            = 32897 ! IT8 Monochrome picture
    integer(C_INT), parameter ::   COMPRESSION_IT8BL            = 32898 ! IT8 Binary line art
    !!  compression codes 32908-32911 are reserved for Pixar
    integer(C_INT), parameter ::   COMPRESSION_PIXARFILM        = 32908 ! Pixar companded 10bit LZW
    integer(C_INT), parameter ::   COMPRESSION_PIXARLOG         = 32909 ! Pixar companded 11bit ZIP
    integer(C_INT), parameter ::   COMPRESSION_DEFLATE          = 32946 ! Deflate compression
    integer(C_INT), parameter ::   COMPRESSION_ADOBE_DEFLATE    = 8     ! Deflate compression as recognized by Adobe
    !!  compression code 32947 is reserved for Oceana Matrix <dev@oceana.com>
    integer(C_INT), parameter ::   COMPRESSION_DCS              = 32947 ! Kodak DCS encoding
    integer(C_INT), parameter ::   COMPRESSION_JBIG             = 34661 ! ISO JBIG
    integer(C_INT), parameter ::   COMPRESSION_SGILOG           = 34676 ! SGI Log Luminance RLE
    integer(C_INT), parameter ::   COMPRESSION_SGILOG24         = 34677 ! SGI Log 24-bit packed
    integer(C_INT), parameter ::   COMPRESSION_JP2000           = 34712 ! Leadtools JPEG2000
    integer(C_INT), parameter ::   COMPRESSION_LZMA             = 34925 ! LZMA2
#define TIFFTAG_PHOTOMETRIC 262
#define PHOTOMETRIC_MINISBLACK 1
    !integer(C_INT), parameter :: TIFFTAG_PHOTOMETRIC            = 262   ! photometric interpretation
    integer(C_INT), parameter ::   PHOTOMETRIC_MINISWHITE       = 0     ! min value is white
    !integer(C_INT), parameter ::   PHOTOMETRIC_MINISBLACK       = 1     ! min value is black
    integer(C_INT), parameter ::   PHOTOMETRIC_RGB              = 2     ! RGB color model
    integer(C_INT), parameter ::   PHOTOMETRIC_PALETTE          = 3     ! color map indexed
    integer(C_INT), parameter ::   PHOTOMETRIC_MASK             = 4     ! $holdout mask
    integer(C_INT), parameter ::   PHOTOMETRIC_SEPARATED        = 5     ! !color separations
    integer(C_INT), parameter ::   PHOTOMETRIC_YCBCR            = 6     ! !CCIR 601
    integer(C_INT), parameter ::   PHOTOMETRIC_CIELAB           = 8     ! !1976 CIE L*a*b*
    integer(C_INT), parameter ::   PHOTOMETRIC_ICCLAB           = 9     ! ICC L*a*b* [Adobe      TIFF Technote 4]
    integer(C_INT), parameter ::   PHOTOMETRIC_ITULAB           = 10    ! ITU L*a*b*
    integer(C_INT), parameter ::   PHOTOMETRIC_CFA              = 32803 ! color filter array
    integer(C_INT), parameter ::   PHOTOMETRIC_LOGL             = 32844 ! CIE Log2(L)
    integer(C_INT), parameter ::   PHOTOMETRIC_LOGLUV           = 32845 ! CIE Log2(L) (u',v')
    integer(C_INT), parameter :: TIFFTAG_THRESHHOLDING          = 263   ! +thresholding used on data
    integer(C_INT), parameter ::   THRESHHOLD_BILEVEL           = 1     ! b&w art scan
    integer(C_INT), parameter ::   THRESHHOLD_HALFTONE          = 2     ! or dithered scan
    integer(C_INT), parameter ::   THRESHHOLD_ERRORDIFFUSE      = 3     ! usually floyd-steinberg
    integer(C_INT), parameter :: TIFFTAG_CELLWIDTH              = 264   ! +dithering matrix width
    integer(C_INT), parameter :: TIFFTAG_CELLLENGTH             = 265   ! +dithering matrix height
    integer(C_INT), parameter :: TIFFTAG_FILLORDER              = 266   ! data order within a byte
    integer(C_INT), parameter ::   FILLORDER_MSB2LSB            = 1     ! most significant -> least
    integer(C_INT), parameter ::   FILLORDER_LSB2MSB            = 2     ! least significant -> most
    integer(C_INT), parameter :: TIFFTAG_DOCUMENTNAME           = 269   ! name of doc. image is from
    integer(C_INT), parameter :: TIFFTAG_IMAGEDESCRIPTION       = 270   ! info about image
    integer(C_INT), parameter :: TIFFTAG_MAKE                   = 271   ! scanner manufacturer name
    integer(C_INT), parameter :: TIFFTAG_MODEL                  = 272   ! scanner model name/number
    integer(C_INT), parameter :: TIFFTAG_STRIPOFFSETS           = 273   ! offsets to data strips
    integer(C_INT), parameter :: TIFFTAG_ORIENTATION            = 274   ! +image orientation
    integer(C_INT), parameter ::   ORIENTATION_TOPLEFT          = 1     ! row 0 top, col 0 lhs
    integer(C_INT), parameter ::   ORIENTATION_TOPRIGHT         = 2     ! row 0 top, col 0 rhs
    integer(C_INT), parameter ::   ORIENTATION_BOTRIGHT         = 3     ! row 0 bottom, col 0 rhs
    integer(C_INT), parameter ::   ORIENTATION_BOTLEFT          = 4     ! row 0 bottom, col 0 lhs
    integer(C_INT), parameter ::   ORIENTATION_LEFTTOP          = 5     ! row 0 lhs, col 0 top
    integer(C_INT), parameter ::   ORIENTATION_RIGHTTOP         = 6     ! row 0 rhs, col 0 top
    integer(C_INT), parameter ::   ORIENTATION_RIGHTBOT         = 7     ! row 0 rhs, col 0 bottom
    integer(C_INT), parameter ::   ORIENTATION_LEFTBOT          = 8     ! row 0 lhs, col 0 bottom
    integer(C_INT), parameter :: TIFFTAG_SAMPLESPERPIXEL        = 277   ! samples per pixel
    integer(C_INT), parameter :: TIFFTAG_ROWSPERSTRIP           = 278   ! rows per strip of data
    integer(C_INT), parameter :: TIFFTAG_STRIPBYTECOUNTS        = 279   ! bytes counts for strips
    integer(C_INT), parameter :: TIFFTAG_MINSAMPLEVALUE         = 280   ! +minimum sample value
    integer(C_INT), parameter :: TIFFTAG_MAXSAMPLEVALUE         = 281   ! +maximum sample value
    integer(C_INT), parameter :: TIFFTAG_XRESOLUTION            = 282   ! pixels/resolution in x
    integer(C_INT), parameter :: TIFFTAG_YRESOLUTION            = 283   ! pixels/resolution in y
    integer(C_INT), parameter :: TIFFTAG_PLANARCONFIG           = 284   ! storage organization
    integer(C_INT), parameter ::   PLANARCONFIG_CONTIG          = 1     ! single image plane
    integer(C_INT), parameter ::   PLANARCONFIG_SEPARATE        = 2     ! separate planes of data
    integer(C_INT), parameter :: TIFFTAG_PAGENAME               = 285   ! page name image is from
    integer(C_INT), parameter :: TIFFTAG_XPOSITION              = 286   ! x page offset of image lhs
    integer(C_INT), parameter :: TIFFTAG_YPOSITION              = 287   ! y page offset of image lhs
    integer(C_INT), parameter :: TIFFTAG_FREEOFFSETS            = 288   ! +byte offset to free block
    integer(C_INT), parameter :: TIFFTAG_FREEBYTECOUNTS         = 289   ! +sizes of free blocks
    integer(C_INT), parameter :: TIFFTAG_GRAYRESPONSEUNIT       = 290   ! $gray scale curve accuracy
    integer(C_INT), parameter ::   GRAYRESPONSEUNIT_10S         = 1     ! tenths of a unit
    integer(C_INT), parameter ::   GRAYRESPONSEUNIT_100S        = 2     ! hundredths of a unit
    integer(C_INT), parameter ::   GRAYRESPONSEUNIT_1000S       = 3     ! thousandths of a unit
    integer(C_INT), parameter ::   GRAYRESPONSEUNIT_10000S      = 4     ! ten-thousandths of a unit
    integer(C_INT), parameter ::   GRAYRESPONSEUNIT_100000S     = 5     ! hundred-thousandths
    integer(C_INT), parameter :: TIFFTAG_GRAYRESPONSECURVE      = 291   ! $gray scale response curve
    integer(C_INT), parameter :: TIFFTAG_GROUP3OPTIONS          = 292   ! 32 flag bits
    integer(C_INT), parameter :: TIFFTAG_T4OPTIONS              = 292   ! TIFF 6.0 proper name alias
    integer(C_INT), parameter ::   GROUP3OPT_2DENCODING         = 1     ! 2-dimensional coding
    integer(C_INT), parameter ::   GROUP3OPT_UNCOMPRESSED       = 2     ! data not compressed
    integer(C_INT), parameter ::   GROUP3OPT_FILLBITS           = 4     ! fill to byte boundary
    integer(C_INT), parameter :: TIFFTAG_GROUP4OPTIONS          = 293   ! 32 flag bits
    integer(C_INT), parameter :: TIFFTAG_T6OPTIONS              = 293   ! TIFF 6.0 proper name
    integer(C_INT), parameter ::   GROUP4OPT_UNCOMPRESSED       = 2     ! data not compressed
    integer(C_INT), parameter :: TIFFTAG_RESOLUTIONUNIT         = 296   ! units of resolutions
    integer(C_INT), parameter ::   RESUNIT_NONE                 = 1     ! no meaningful units
    integer(C_INT), parameter ::   RESUNIT_INCH                 = 2     ! english
    integer(C_INT), parameter ::   RESUNIT_CENTIMETER           = 3     ! metric
    integer(C_INT), parameter :: TIFFTAG_PAGENUMBER             = 297   ! page numbers of multi-page
    integer(C_INT), parameter :: TIFFTAG_COLORRESPONSEUNIT      = 300   ! $color curve accuracy
    integer(C_INT), parameter ::    COLORRESPONSEUNIT_10S       = 1     ! tenths of a unit
    integer(C_INT), parameter ::    COLORRESPONSEUNIT_100S      = 2     ! hundredths of a unit
    integer(C_INT), parameter ::    COLORRESPONSEUNIT_1000S     = 3     ! thousandths of a unit
    integer(C_INT), parameter ::    COLORRESPONSEUNIT_10000S    = 4     ! ten-thousandths of a unit
    integer(C_INT), parameter ::    COLORRESPONSEUNIT_100000S   = 5     ! hundred-thousandths
    integer(C_INT), parameter :: TIFFTAG_TRANSFERFUNCTION       = 301   ! !colorimetry info
    integer(C_INT), parameter :: TIFFTAG_SOFTWARE               = 305   ! name & release
    integer(C_INT), parameter :: TIFFTAG_DATETIME               = 306   ! creation date and time
    integer(C_INT), parameter :: TIFFTAG_ARTIST                 = 315   ! creator of image
    integer(C_INT), parameter :: TIFFTAG_HOSTCOMPUTER           = 316   ! machine where created
    integer(C_INT), parameter :: TIFFTAG_PREDICTOR              = 317   ! prediction scheme w/ LZW
    integer(C_INT), parameter ::   PREDICTOR_NONE               = 1     ! no prediction scheme used
    integer(C_INT), parameter ::   PREDICTOR_HORIZONTAL         = 2     ! horizontal differencing
    integer(C_INT), parameter ::   PREDICTOR_FLOATINGPOINT      = 3     ! floating point predictor
    integer(C_INT), parameter :: TIFFTAG_WHITEPOINT             = 318   ! image white point
    integer(C_INT), parameter :: TIFFTAG_PRIMARYCHROMATICITIES  = 319   ! !primary chromaticities
    integer(C_INT), parameter :: TIFFTAG_COLORMAP               = 320   ! RGB map for pallette image
    integer(C_INT), parameter :: TIFFTAG_HALFTONEHINTS          = 321   ! !highlight+shadow info
    integer(C_INT), parameter :: TIFFTAG_TILEWIDTH              = 322   ! !tile width in pixels
    integer(C_INT), parameter :: TIFFTAG_TILELENGTH             = 323   ! !tile height in pixels
    integer(C_INT), parameter :: TIFFTAG_TILEOFFSETS            = 324   ! !offsets to data tiles
    integer(C_INT), parameter :: TIFFTAG_TILEBYTECOUNTS         = 325   ! !byte counts for tiles
    integer(C_INT), parameter :: TIFFTAG_BADFAXLINES            = 326   ! lines w/ wrong pixel count
    integer(C_INT), parameter :: TIFFTAG_CLEANFAXDATA           = 327   ! regenerated line info
    integer(C_INT), parameter ::   CLEANFAXDATA_CLEAN           = 0     ! no errors detected
    integer(C_INT), parameter ::   CLEANFAXDATA_REGENERATED     = 1     ! receiver regenerated lines
    integer(C_INT), parameter ::   CLEANFAXDATA_UNCLEAN         = 2     ! uncorrected errors exist
    integer(C_INT), parameter :: TIFFTAG_CONSECUTIVEBADFAXLINES = 328   ! max consecutive bad lines
    integer(C_INT), parameter :: TIFFTAG_SUBIFD                 = 330   ! subimage descriptors
    integer(C_INT), parameter :: TIFFTAG_INKSET                 = 332   ! !inks in separated image
    integer(C_INT), parameter ::   INKSET_CMYK                  = 1     ! !cyan-magenta-yellow-black color
    integer(C_INT), parameter ::   INKSET_MULTIINK              = 2     ! !multi-ink or hi-fi color
    integer(C_INT), parameter :: TIFFTAG_INKNAMES               = 333   ! !ascii names of inks
    integer(C_INT), parameter :: TIFFTAG_NUMBEROFINKS           = 334   ! !number of inks
    integer(C_INT), parameter :: TIFFTAG_DOTRANGE               = 336   ! !0% and 100% dot codes
    integer(C_INT), parameter :: TIFFTAG_TARGETPRINTER          = 337   ! !separation target
    integer(C_INT), parameter :: TIFFTAG_EXTRASAMPLES           = 338   ! !info about extra samples
    integer(C_INT), parameter ::   EXTRASAMPLE_UNSPECIFIED      = 0     ! !unspecified data
    integer(C_INT), parameter ::   EXTRASAMPLE_ASSOCALPHA       = 1     ! !associated alpha data
    integer(C_INT), parameter ::   EXTRASAMPLE_UNASSALPHA       = 2     ! !unassociated alpha data
    integer(C_INT), parameter :: TIFFTAG_SAMPLEFORMAT           = 339   ! !data sample format
    integer(C_INT), parameter ::   SAMPLEFORMAT_UINT            = 1     ! !unsigned integer data
    integer(C_INT), parameter ::   SAMPLEFORMAT_INT             = 2     ! !signed integer data
    integer(C_INT), parameter ::   SAMPLEFORMAT_IEEEFP          = 3     ! !IEEE floating point data
    integer(C_INT), parameter ::   SAMPLEFORMAT_VOID            = 4     ! !untyped data
    integer(C_INT), parameter ::   SAMPLEFORMAT_COMPLEXINT      = 5     ! !complex signed int
    integer(C_INT), parameter ::   SAMPLEFORMAT_COMPLEXIEEEFP   = 6     ! !complex ieee floating
    integer(C_INT), parameter :: TIFFTAG_SMINSAMPLEVALUE        = 340   ! !variable MinSampleValue
    integer(C_INT), parameter :: TIFFTAG_SMAXSAMPLEVALUE        = 341   ! !variable MaxSampleValue
    integer(C_INT), parameter :: TIFFTAG_CLIPPATH               = 343   ! %ClipPath         [Adobe TIFF technote 2]
    integer(C_INT), parameter :: TIFFTAG_XCLIPPATHUNITS         = 344   ! %XClipPathUnits         [Adobe TIFF technote 2]
    integer(C_INT), parameter :: TIFFTAG_YCLIPPATHUNITS         = 345   ! %YClipPathUnits         [Adobe TIFF technote 2]
    integer(C_INT), parameter :: TIFFTAG_INDEXED                = 346   ! %Indexed         [Adobe TIFF Technote 3]
    integer(C_INT), parameter :: TIFFTAG_JPEGTABLES             = 347   ! %JPEG table stream
    integer(C_INT), parameter :: TIFFTAG_OPIPROXY               = 351   ! %OPI Proxy [Adobe TIFF technote]
    !!  Tags 400-435 are from the  TIFF/FX spec


    ! TIFFTAG_ARTIST                  1      char*
    ! TIFFTAG_BADFAXLINES             1      uint32
    ! TIFFTAG_BITSPERSAMPLE           1      uint16             †
    ! TIFFTAG_CLEANFAXDATA            1      uint16
    ! TIFFTAG_COLORMAP                3      uint16*            1<<BitsPerSample arrays
    ! TIFFTAG_COMPRESSION             1      uint16             †
    ! TIFFTAG_CONSECUTIVEBADFAXLINES  1      uint32
    ! TIFFTAG_COPYRIGHT               1      char*
    ! TIFFTAG_DATETIME                1      char*
    ! TIFFTAG_DOCUMENTNAME            1      char*
    ! TIFFTAG_DOTRANGE                2      uint16
    ! TIFFTAG_EXTRASAMPLES            2      uint16,uint16*     † count & types array
    ! TIFFTAG_FAXFILLFUNC             1      TIFFFaxFillFunc    G3/G4 compression pseudo-tag
    ! TIFFTAG_FAXMODE                 1      int                † G3/G4 compression pseudo-tag
    ! TIFFTAG_FILLORDER               1      uint16             †
    ! TIFFTAG_GROUP3OPTIONS           1      uint32             †
    ! TIFFTAG_GROUP4OPTIONS           1      uint32             †
    ! TIFFTAG_HALFTONEHINTS           2      uint16
    ! TIFFTAG_HOSTCOMPUTER            1      char*
    ! TIFFTAG_ICCPROFILE              2      uint32,void*       count, profile data
    ! TIFFTAG_IMAGEDEPTH              1      uint32             †
    ! TIFFTAG_IMAGEDESCRIPTION        1      char*
    ! TIFFTAG_IMAGELENGTH             1      uint32
    ! TIFFTAG_IMAGEWIDTH              1      uint32             †
    ! TIFFTAG_INKNAMES                2      uint16, char*
    ! TIFFTAG_INKSET                  1      uint16             †
    ! TIFFTAG_JPEGCOLORMODE           1      int                † JPEG pseudo-tag
    ! TIFFTAG_JPEGQUALITY             1      int                JPEG pseudo-tag
    ! TIFFTAG_JPEGTABLES              2      uint32*,void*      † count & tables
    ! TIFFTAG_JPEGTABLESMODE          1      int                † JPEG pseudo-tag
    ! TIFFTAG_MAKE                    1      char*
    ! TIFFTAG_MATTEING                1      uint16             †
    ! TIFFTAG_MAXSAMPLEVALUE          1      uint16
    ! TIFFTAG_MINSAMPLEVALUE          1      uint16

    ! TIFFTAG_MODEL                   1      char*
    ! TIFFTAG_ORIENTATION             1      uint16
    ! TIFFTAG_PAGENAME                1      char*
    ! TIFFTAG_PAGENUMBER              2      uint16
    ! TIFFTAG_PHOTOMETRIC             1      uint16
    ! TIFFTAG_PHOTOSHOP               ?      uint32,void*       count, data
    ! TIFFTAG_PLANARCONFIG            1      uint16             †
    ! TIFFTAG_PREDICTOR               1      uint16             †
    ! TIFFTAG_PRIMARYCHROMATICITIES   1      float*             6-entry array
    ! TIFFTAG_REFERENCEBLACKWHITE     1      float*             † 6-entry array
    ! TIFFTAG_RESOLUTIONUNIT          1      uint16
    ! TIFFTAG_RICHTIFFIPTC            2      uint32,void*       count, data
    ! TIFFTAG_ROWSPERSTRIP            1      uint32             † must be > 0
    ! TIFFTAG_SAMPLEFORMAT            1      uint16             †
    ! TIFFTAG_SAMPLESPERPIXEL         1      uint16             † value must be <= 4
    ! TIFFTAG_SMAXSAMPLEVALUE         1      double
    ! TIFFTAG_SMINSAMPLEVALUE         1      double
    ! TIFFTAG_SOFTWARE                1      char*
    ! TIFFTAG_STONITS                 1      double             †
    ! TIFFTAG_SUBFILETYPE             1      uint32
    ! TIFFTAG_SUBIFD                  2      uint16,uint32*     count & offsets array

    ! TIFFTAG_THRESHHOLDING           1      uint16
    ! TIFFTAG_TILEDEPTH               1      uint32             †
    ! TIFFTAG_TILELENGTH              1      uint32             † must be a multiple of 8
    ! TIFFTAG_TILEWIDTH               1      uint32             † must be a multiple of 8
    ! TIFFTAG_TRANSFERFUNCTION        1 or 3‡ uint16*           1<<BitsPerSample entry arrays
    !     TIFFTAG_WHITEPOINT              1      float*             2-entry array

    !     TIFFTAG_XPOSITION               1      float
    !     TIFFTAG_XRESOLUTION             1      float

    !     TIFFTAG_YPOSITION               1      float
    !     TIFFTAG_YRESOLUTION             1      float


    !! TIFFIO.H functions
    interface
        !! TIFFOpen returns a TIFF pointer.  Otherwise, NULL is returned.
        type(C_PTR) function TIFFOpen(fname, mode) bind(C, name="TIFFOpen")
            import c_char, c_ptr
            character(len=1,kind=c_char), dimension(*), intent(in) :: fname
            character(len=1,kind=c_char), dimension(*), intent(in) :: mode
        end function TIFFOpen
        !> int TIFFSetField(TIFF *tif, ttag_t tag, ...)
        integer(C_INT) function TIFFSetField(tif, flag, val) bind(C,name="TIFFSetField")
             import
             type(C_PTR), intent(inout) :: tif
             integer(C_INT), intent(in), value :: flag       !! ttag_t is unsigned int
             integer(C_INT), intent(in) :: val
        end function TIFFSetField
        ! subroutine TIFFSetFieldf(tif, flag, val) bind(C,name="TIFFSetField")
        !     import
        !     type(C_PTR), intent(inout) :: tif
        !     integer(C_INT), intent(in) :: flag
        !     real(C_FLOAT), intent(in) :: val
        ! end subroutine TIFFSetFieldf
        ! subroutine TIFFSetFieldc(tif, flag, val) bind(C,name="TIFFSetField")
        !     import
        !     type(C_PTR), intent(inout) :: tif
        !     integer(C_INT), intent(in) :: flag
        !     character(kind=C_CHAR, len=1), dimension(*), intent(in) :: val
        ! end subroutine TIFFSetFieldc
        integer(C_INT) function TIFFGetField(tif, flag, val) bind(C,name="TIFFGetField")
            import
            type(C_PTR), intent(inout) :: tif
            integer(C_INT), intent(in) :: flag
            integer(C_INT), intent(in) :: val
        end function  TIFFGetField
        integer(C_INT) function TIFFWriteScanline(tif, buffer, row, sample) bind(C, name="TIFFWriteScanline")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            type(C_PTR) :: buffer ! (void*)
            integer(C_INT32_T), intent(in) :: row
            integer(C_INT16_T), intent(in) :: sample ! default =0
        end function TIFFWriteScanline
        !! TIFFWriteTile returns -1 if it detects an error; otherwise the number of bytes in the tile is returned.
        !! tsize_t TIFFWriteTile(TIFF *tif, tdata_t buf, uint32 x, uint32 y, uint32 z, tsample_t sample)
        integer(C_INT) function TIFFWriteTile(tif, buffer, x, y, z, sample) bind(C, name="TIFFWriteTile")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            type(C_PTR) :: buffer ! (void*)
            integer(C_INT32_T), intent(in) :: x,y,z  ! uint32
            integer(C_INT16_T), intent(in) :: sample ! default =0
        end function TIFFWriteTile
        integer(C_INT) function TIFFWriteRawTile(tif, tile, buffer, s) bind(C, name="TIFFWriteTile")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            integer(C_INT32_T) :: tile ! ttile_t
            type(C_PTR) :: buffer ! (void*)
            integer(C_INT64_T), intent(in) :: s ! default =0
        end function TIFFWriteRawTile
        !tsize_t TIFFWriteEncodedStrip(TIFF *tif, tstrip_t strip,  tdata_t  buf, tsize_t size)
       integer(C_INT) function TIFFWriteEncodedStrip(tif, strip, buf, sz) bind(C, name="TIFFWriteEncodedStrip")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            integer(C_INT32_T), intent(in) :: strip
            type(C_PTR) :: buf ! (void*)
            integer(C_INT), intent(in) :: sz
        end function TIFFWriteEncodedStrip
        ! tsize_t  TIFFWriteRawStrip(TIFF  *tif,  tstrip_t  strip,  tdata_t  buf,  tsize_t size)
        integer(C_INT) function TIFFWriteRawStrip(tif, strip, buf, sz) bind(C, name="TIFFWriteRawStrip")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            integer(C_INT32_T), intent(in) :: strip
            type(C_PTR) :: buf ! (void*)
            integer(C_INT), intent(in) :: sz
        end function TIFFWriteRawStrip
        integer(C_INT) function TIFFWriteDirectory(tif) bind(C, name="TIFFWriteDirectory")
            import
            type(C_PTR) :: tif  ! (TIFF*)
        end function TIFFWriteDirectory
        !! TIFFReadTile  returns  -1  if  it  detects  an  error; otherwise the number of bytes in the decoded tile is returned.
        integer(C_INT) function TIFFReadTile(tif, buffer, x, y, z, sample) bind(C, name="TIFFReadTile")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            type(C_PTR) :: buffer ! (void*)
            integer(C_INT32_T), intent(in) :: x,y,z  ! uint32
            integer(C_INT16_T), intent(in) :: sample ! default =0
        end function TIFFReadTile
          integer(C_INT) function TIFFComputeTile(tif,  x, y, z, sample) bind(C, name="TIFFComputeTile")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            integer(C_INT32_T), intent(in) :: x,y,z  ! uint32
            integer(C_INT16_T), intent(in) :: sample ! default =0
        end function TIFFComputeTile
        integer(C_INT64_T) function TIFFTileSize(tif) bind(C, name="TIFFTileSize")
            import
            type(C_PTR) :: tif  ! (TIFF*)
        end function TIFFTileSize
        integer(C_INT) function TIFFCheckTile(tif, x, y, z, sample) bind(C, name="TIFFCheckTile")
            import
            type(C_PTR) :: tif  ! (TIFF*)
            integer(C_INT32_T), intent(in) :: x,y,z  ! uint32
            integer(C_INT16_T), intent(in) :: sample ! default =0
        end function TIFFCheckTile
        function TIFFmalloc(sz) bind(C,name="_TIFFmalloc")
            import
            type(C_PTR) :: TIFFmalloc
            integer(C_INT64_T), intent(in) :: sz
        end function TIFFmalloc
        subroutine TIFFmemcpy(pdata, cdata, sz) bind(C, name="_TIFFfree")
            import
            type(C_PTR) :: pdata
            type(C_PTR) :: cdata
            integer(C_INT64_T), intent(in), value :: sz
        end subroutine TIFFmemcpy
        subroutine TIFFfree(pdata) bind(C, name="_TIFFfree")
            import
            type(C_PTR) :: pdata
        end subroutine TIFFfree
        subroutine TIFFClose(tif) bind(C, name="TIFFClose")
            import
            type(C_PTR), intent(inout) :: tif
        end subroutine TIFFClose

    end interface

interface
    ! Ignore the return value of strncpy -> subroutine
    ! "restrict" is always assumed if we do not pass a pointer
    subroutine strncpy(dest, src, n) bind(C)
      import
      character(kind=c_char),  intent(out) :: dest(*)
      character(kind=c_char),  intent(in)  :: src(*)
      integer(c_size_t), value, intent(in) :: n
    end subroutine strncpy
end interface

contains

    subroutine write_tiff( fname , img , ppc)
        character(len=*), intent(in) :: fname
        real, intent(in) :: img(:,:,:)
        real, optional, intent(in):: ppc
        type(C_PTR) :: tif
        type(C_PTR) :: l_buffer, f_ptr
        integer(C_INT64_T) :: l_buffer_size  !!  this is the machine addressing size type
        integer(C_INT32_T) :: width, height, slices, tileLength, tileWidth
        integer(kind=C_INT16_T), pointer :: data(:)
        integer(C_INT) :: cx, cy
        integer ::  x,y,sl,fname_sz, ret, status, row
        real(C_FLOAT) :: pixels_per_cm
        real(8) :: d
        real :: imgrange(2)
        character(len=1, kind=C_CHAR), dimension(:), allocatable :: fname_ptr

        pixels_per_cm = 100.0
        if(present(ppc)) pixels_per_cm = ppc
        width = size(img,1)
        height = size(img,2)
        slices = size(img,3)
        tileWidth = width
        if(slices > 1)then
            tileLength =height*( slices )
        else
            tileLength =height
        endif

        cx = width / 2
        cy = height / 2
        fname_sz=len_trim(fname)
        allocate(fname_ptr(fname_sz+1))
        do x=1, fname_sz
            fname_ptr(x) = fname(x:x)
        enddo
        fname_ptr(x+1)=achar(0)
        allocate(data(width * height))
        imgrange(1) =  minval(img)
        imgrange(2) = maxval(img)
        sl=1
        do y = 1,height
            do x = 1, width
                d =  img(x, y, sl)
                d = (d - imgrange(1))/(imgrange(2)-imgrange(1))
                data((y-1) * width + x) = INT( d * 65534 , C_INT16_T )
            end do
        end do

        tif = TIFFOpen(fname_ptr, "w8")  ! bigtiff

        ret = TIFFSetField(tif, TIFFTAG_IMAGEWIDTH     , width)                  !uint32
        ret = TIFFSetField(tif, TIFFTAG_IMAGELENGTH    , height)                 !uint32
        ret = TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE  , 16)                     !uint16
        ret = TIFFSetField(tif, TIFFTAG_COMPRESSION    , COMPRESSION_NONE)       !uint16
        ret = TIFFSetField(tif, TIFFTAG_PHOTOMETRIC    , 1)! PHOTOMETRIC_MINISBLACK) !uint16
        ! call TIFFSetField(tif, TIFFTAG_FILLORDER      , FILLORDER_MSB2LSB)      !uint16
        ! call TIFFSetField(tif, TIFFTAG_ORIENTATION    , ORIENTATION_TOPLEFT)    !uint16
        ret = TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1)                      !uint16
        ret = TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP   , 8)                      !uint32
        ret = TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT , RESUNIT_CENTIMETER)     !uint16
        ! call TIFFSetField(tif, TIFFTAG_XRESOLUTION    , pixels_per_cm)          !float
        ! call TIFFSetField(tif, TIFFTAG_YRESOLUTION    , pixels_per_cm)          !float
        ret = TIFFSetField(tif, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG)    !uint16
        ret = TIFFSetField(tif, TIFFTAG_TILEWIDTH, tileWidth)
        ret = TIFFSetField(tif, TIFFTAG_TILELENGTH, tileLength)

        l_buffer_size = width*2  !! LINE buffer for 16bit
        f_ptr = c_loc(data(1))
        l_buffer = TIFFmalloc(l_buffer_size)

        do row = 1, height
            do x = 1, width
                d =  img(x, row, sl)
                d = (d - imgrange(1))/(imgrange(2)-imgrange(1))
                data((y-1) * width + x) = INT( d * 65534 , C_INT16_T )
            end do
            f_ptr = c_loc(data(1))
            call TIFFmemcpy(l_buffer, f_ptr, l_buffer_size)

            ret=TIFFWriteScanline(tif, l_buffer, row, 0_2)
            if (ret==-1)then
                call TIFFClose(tif)
                status= -1
            endif
        end do

        call TIFFfree(l_buffer)

        call TIFFClose(tif)

        status=0

    end subroutine write_tiff





    subroutine read_tiff(fname, img, width, height, xtiles, ytiles)
         character(len=*),  intent(in)  :: fname
         real, allocatable, intent(out) :: img(:,:,:)
         integer, intent(out) :: width, height, xtiles, ytiles
         type(c_ptr) :: tif
         character(kind=C_CHAR, len=1), allocatable :: filename(:)
         integer(C_INT32_T) :: imageWidth, imageLength
         integer(C_INT32_T) :: tileWidth, tileLength, ntiles
         integer(C_INT32_T) :: x, y,z
         type(c_ptr) :: buf   !  tdata_t == (void*)
         integer, pointer :: bufint(:)
         integer(C_INT32_T) ::  npixels, ret, tilenumber
         integer(C_INT64_T) :: bufsize
         integer(C_INT16_T) :: sample =0
         integer :: itile
         filename =trim(adjustl(fname))//achar(0)
         tif = TIFFOpen(fname, "r")
         if ( c_associated(tif) ) then
             ret= TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,  imageWidth)
             ret= TIFFGetField(tif, TIFFTAG_IMAGELENGTH, imageLength)
             ret= TIFFGetField(tif, TIFFTAG_TILEWIDTH,   tileWidth)
             ret= TIFFGetField(tif, TIFFTAG_TILELENGTH,  tileLength)
             bufsize=TIFFTileSize(tif)
             buf = TIFFmalloc(bufsize)
             if(.not. c_associated(buf)) then
                 write(*,*) "Error: insufficient memory", trim(fname)
                 call TIFFClose(tif)
                 return
             endif

             xtiles = NINT(REAL(imageWidth)/REAL(tileWidth))
             allocate(img(tileWidth,tileLength,xtiles))
             npixels = tileWidth*tileLength
             itile = 1
             do y = 0, imageLength, tileLength

                 do x = 0, imageWidth, tileWidth
                     tilenumber= TIFFComputeTile(tif, x, y, 0, sample)
                     ret  = TIFFReadTile(tif, buf, x, y, z, sample)
                     !! fixme                    for(register int i=0; i<npixels;i=i+1) realimg[(y*imageWidth)+x+i] = (float) buf[i];
                     call c_f_pointer(buf, bufint, [bufsize])
                     img(1:tileWidth,1:tileLength,itile) = reshape(bufint, shape=(/ tileWidth, tileLength /))
                     itile= itile+1
                 end do
             end do
             call TIFFfree(buf);
             call TIFFClose(tif);
         endif
         width=tileWidth
         height=tileLength
         xtiles=imageWidth/tileWidth
         ytiles=imageLength/tileLength
     end subroutine read_tiff


    !!
    !! The TIFF writer loosly based on  http://paulbourke.net/dataformats/tiff/
     function c_putc(funit, str) result(istat)
        use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, iostat_eor, iostat_end
        integer, intent(in) :: funit
        character(len=1, kind=C_CHAR), dimension(*), intent(in) :: str
        character(len=1) :: c
        integer :: istat,i
        character(len=256) :: imsg
        do i=1, len(str)
            write(c,'(a)') str(i:i)
            write(funit, fmt='(g0)', advance="no", iostat=istat, iomsg=imsg) c
            if ( istat /= 0 ) then
                write (output_unit, fmt='(*(g0))') 'io failed: ', imsg
                stop
            end if
        end do
    end function c_putc

    ! subroutine WriteHexString(funit, str)
    !     integer, intent(in) :: funit
    !     character(len=*), intent(in) :: str
    !     integer(C_INT) :: RBG
    !     character(len=6) :: hexval
    !     integer :: lstr, R,G,B
    !     ! RGB gets RGB integer value somehow...
    !     write (hexval,'(Z6.6)') RGB
    !     READ (hexval,'(3Z2.2)') R,G,B
    !     write(funit,'(B16)') R,G,B
    ! end subroutine WriteHexString
    subroutine swapImgBytes(img, imgBytes)
        integer(1), pointer, intent(inout) :: img(:)
        integer, intent(in) :: imgBytes
        integer(1) :: sw
        integer :: hi = 2
        do while (hi < imgBytes)
            sw = img(hi)
            img(hi) = img(hi-1)
            img(hi-1) = sw
            hi = hi + 2
        end do
    end subroutine swapImgBytes

    subroutine normalize16bit(img, px)
        integer(C_INT32_T),  intent(inout) :: img(:)
        integer, intent(in) :: px
        !! in general, TIFF readers tend to ignore the Maximum sample value tag
        !! they expect 16-bit images to be in the range 0..65535
        !! be warned this could screw up DICOM scale factors....
        integer :: mx, mn,i, h
        real(8) :: scale
        mn = minval(img)
        if(mn < 0 ) img = img - mn
        mx= maxval(img)
        if (mx == 0) return
        scale = 65535./REAL(mx,8)
        do i=1,px
            img(i) = scale * img(i)
        end do
    end subroutine normalize16bit
    subroutine normalize24bit(img, px)
        integer(C_INT32_T),  intent(inout) :: img(:)
        integer, intent(in) :: px
        !! in general, TIFF readers tend to ignore the Maximum sample value tag
        !! they expect 16-bit images to be in the range 0..65535
        !! be warned this could screw up DICOM scale factors....
        integer :: mx, mn,i, h
        real(8) :: scale
        mn = minval(img)
        if(mn < 0 ) img = img - mn
        mx= maxval(img)
        if (mx == 0) return
        scale =16777215./REAL(mx,8)
        do i=1,px
            img(i) = scale * img(i)
        end do
    end subroutine normalize24bit

    function write_tiff_bigimg(fname, img_input, nx,  ny, bits, frames, isPlanar) result(istat)
        character(len=*), intent(in) :: fname
        real, intent(in) :: img_input(:,:)
        integer, intent(in) :: nx,ny
        integer, optional, intent(in) :: bits, frames, isPlanar
        integer :: istat,imgBytes,iostat, fptr,  bits_here, frames_here, i,j,h
        logical :: isPlanar_here
        integer(C_INT8_T),  pointer :: img_buffer(:) => NULL()
        integer(C_INT32_T) :: pixel, offset
        real :: hi,lo
        bits_here=16
        if(present(bits)) bits_here=bits
        istat=1
        if ((bits < 1) .or. (bits > 16)) return
        frames_here = 3
        if(present(frames))frames_here=frames
        isPlanar_here=.false.
        if(present(isPlanar)) isPlanar_here= isPlanar/=0
        imgBytes = nx*ny*frames
        if(nx /= size(img_input,1) .or. ny /= size(img_input,2))then

        endif
        lo = minval(img_input)
        hi = maxval(img_input)
        if (imgBytes < 1) return
        allocate(img_buffer(imgBytes))
        do i=0, nx-1
            do j=0, ny-1
                if  (frames_here == 3) then
                    pixel = NINT( (2**16) * (img_input(i+1,j+1)-lo)/(hi-lo),kind=4)
                    img_buffer((i*ny + j ) * 3 + 1)  = INT( ISHFT( pixel , -16) ,kind=c_char)
                    img_buffer((i*ny + j ) * 3 + 2) =  INT( IAND( ISHFT( pixel , -8) , 255) ,kind=c_char)
                    img_buffer((i*ny + j ) * 3 + 3) =  INT( IAND( pixel , 255) ,kind=c_char)
                else if  (frames_here == 4) then
                    pixel = NINT( (2**24) * (img_input(i+1,j+1)-lo)/(hi-lo),kind=4)
                    img_buffer((i*ny + j ) * 3 + 1) = INT(      ISHFT( pixel , -24)      , kind=c_char)

                    img_buffer((i*ny + j ) * 3 + 2) = INT( IAND(ISHFT( pixel , -16), 255), kind=c_char)
                    img_buffer((i*ny + j ) * 3 + 3) = INT( IAND(ISHFT( pixel , -8) , 255), kind=c_char)
                    img_buffer((i*ny + j ) * 3 + 4) = INT( IAND(       pixel       , 255), kind=c_char)
                else
                    pixel =  INT( REAL( (img_input(i+1,j+1)-lo)/REAL(hi - lo) ) ,kind=c_int)
                    pixel =  IAND( pixel , 255)
                    img_buffer(i*ny + j + 1) = INT(pixel,kind=1)
                end if
            end do
        end do

        !         if (bits_here > 8) call normalize16bit( img_buffer, imgBytes)
        !         if (bits_here > 8) imgBytes = 2 * imgBytes
        !         if (isLittleEndian() .and. (bits > 8)) call swapImgBytes(img, imgBytes)
        !         if (isPlanar_here) call deplanar(img, nx, ny, bits, frames)
        open(unit=fptr, file=trim(fname), iostat=iostat, access="STREAM", status="REPLACE")
        if (iostat/=0)then
            write(*,*) "Unable to open tif file in unformatted stream mode"
            return
        endif

        !!  Write the header
        write(fptr) Z"4d4d002a"    !!  Big endian & TIFF identifier
        offset = imgBytes + 8
        write(fptr)INT( ISHFT(offset, -24)  , 1)                   !     putc((offset & 0xff000000) / 16777216,fptr)
        write(fptr)INT( IAND(ISHFT(offset, -16), 255)  , 1)        !     putc((offset & 0x00ff0000) / 65536,fptr)
        write(fptr)INT( IAND(ISHFT(offset, -8), 255)  , 1)         !     putc((offset & 0x0000ff00) / 256,fptr)
        write(fptr)INT( IAND(offset, 255)  , 1)                    !     putc((offset & 0x000000ff),fptr)


        !!  Write the binary data
        write(fptr) img_buffer          !(img,1,imgBytes, fptr)

        !!  Write the footer
        write(fptr)Z"000e"  !!  The number of directory entries (14)
        !!  Width tag, short int
        write(fptr)Z"0100000300000001"


        write(fptr)INT( IAND(ISHFT(nx, -8), 255)  , 1)!     fputc((nx & 0xff00) / 256,fptr)    !!  Image width
        write(fptr)INT( IAND(nx, 255)  , 1)           !     fputc((nx & 0x00ff),fptr)
        write(fptr)Z"0000"
        !!  Height tag, short int
        write(fptr)Z"0101000300000001"
        write(fptr)INT( IAND(ISHFT(ny, -8), 255)  , 1) !     fputc((ny & 0xff00) / 256,fptr)    !!  Image height
        write(fptr)INT( IAND(ny, 255)  , 1) !     fputc((ny & 0x00ff),fptr)
        write(fptr)Z"0000"!     call WriteHexString(fptr,Z"0000")

        !!  Bits per sample tag, short int
        if (frames_here == 3) then
            write(fptr)Z"0102000300000003" !         call WriteHexString(fptr,Z"0102000300000003")
            offset = imgBytes + 182
            write(fptr)INT( ISHFT(offset, -24)  , 1)            !         putc((offset & 0xff000000) / 16777216,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -16), 255)  , 1) !         putc((offset & 0x00ff0000) / 65536,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -8), 255)  , 1)  !         putc((offset & 0x0000ff00) / 256,fptr)
            write(fptr)INT( IAND(offset, 255)  , 1)             !         putc((offset & 0x000000ff),fptr)
        else
            if (bits_here > 8)then
                write(fptr)Z"010200030000000100100000"!             call   WriteHexString(fptr,Z"010200030000000100100000")
            else
                write(fptr)Z"010200030000000100080000"!             call   WriteHexString(fptr,Z"010200030000000100080000")
            endif
        endif

        !     !!  Compression flag, short int
        write(fptr)Z"010300030000000100010000"

        !!  Photometric interpolation tag, short int
        if (frames_here == 3)then
            write(fptr)Z"010600030000000100020000"   !         call WriteHexString(fptr,Z"010600030000000100020000") //2 = RGB
        else !! http://www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
            write(fptr)Z"010600030000000100010000"   !         call WriteHexString(fptr,Z"010600030000000100010000") //1 = BlackIsZero
        endif
        !!  Strip offset tag, long int
        write(fptr)Z"011100040000000100000008"     !           call WriteHexString(fptr,Z"011100040000000100000008")

        !!  Orientation flag, short int
        write(fptr)Z"011200030000000100010000"     !           call WriteHexString(fptr,Z"011200030000000100010000")

        !!  Sample per pixel tag, short int
        if (frames == 3) then
            write(fptr)Z"011500030000000100030000" !         call WriteHexString(fptr,Z"011500030000000100030000")
        else
            write(fptr)Z"011500030000000100010000" !         call WriteHexString(fptr,Z"011500030000000100010000")
        endif
        !!  Rows per strip tag, short int
        write(fptr)Z"0116000300000001"             !     call WriteHexString(fptr,Z"0116000300000001")
        write(fptr)  !     fputc((ny & 0xff00) / 256,fptr)
        write(fptr) !     fputc((ny & 0x00ff),fptr)
        write(fptr)Z"0000"                         !     call  WriteHexString(fptr,Z"0000")

        !!  Strip byte count flag, long int
        write(fptr)Z"0117000400000001"             !     call WriteHexString(fptr,Z"0117000400000001")
        offset = imgBytes
        write(fptr)INT( ISHFT(offset, -24)  , 1)             !     putc((offset & 0xff000000) / 16777216,fptr)
        write(fptr)INT( IAND(ISHFT(offset, -16), 255)  , 1)  !     putc((offset & 0x00ff0000) / 65536,fptr)
        write(fptr)INT( IAND(ISHFT(offset, -8), 255)  , 1)   !     putc((offset & 0x0000ff00) / 256,fptr)
        write(fptr)INT( IAND(offset, 255)  , 1)              !     putc((offset & 0x000000ff),fptr)

        !!  Minimum sample value flag, short int
        if (frames == 3)  then
            write(fptr)Z"0118000300000003" !         WriteHexString(fptr,Z"0118000300000003")
            offset = imgBytes + 188
            write(fptr)INT( ISHFT(offset, -24)  , 1)             !         putc((offset & 0xff000000) / 16777216,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -16), 255)  , 1)  !         putc((offset & 0x00ff0000) / 65536,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -8), 255)  , 1)   !         putc((offset & 0x0000ff00) / 256,fptr)
            write(fptr)INT( IAND(offset, 255)  , 1)              !         putc((offset & 0x000000ff),fptr)
        else
            write(fptr)Z"011800030000000100000000" !         call  WriteHexString(fptr,Z"011800030000000100000000")
        endif
        !!  Maximum sample value tag, short int
        if (frames == 3)  then
            write(fptr)Z"0119000300000003" !         call   WriteHexString(fptr,Z"0119000300000003")
            offset = imgBytes + 194
            write(fptr)INT( ISHFT(offset, -24)  , 1)             !         putc((offset & 0xff000000) / 16777216,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -16), 255)  , 1)  !         putc((offset & 0x00ff0000) / 65536,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -8), 255)  , 1)   !         putc((offset & 0x0000ff00) / 256,fptr)
            write(fptr)INT( IAND(offset, 255)  , 1)              !         putc((offset & 0x000000ff),fptr)
        else
            if (bits > 8) then
                write(fptr)Z"0119000300000001FFFF0000" !             call  WriteHexString(fptr,Z"0119000300000001FFFF0000")
            else
                write(fptr)Z"011900030000000100FF0000" !             call  WriteHexString(fptr,Z"011900030000000100FF0000")
            endif

        endif

        !!  Planar configuration tag, short int
        write(fptr)Z"011c00030000000100010000" !     call WriteHexString(fptr,Z"011c00030000000100010000")

        !!  Sample format tag, short int
        if (frames == 3)  then
            write(fptr)Z"0153000300000003"  !         call WriteHexString(fptr,Z"0153000300000003")
            offset = imgBytes + 200 ! 200
            write(fptr)INT( ISHFT(offset, -24)  , 1)             !         putc((offset & 0xff000000) / 16777216,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -16), 255)  , 1)  !         putc((offset & 0x00ff0000) / 65536,fptr)
            write(fptr)INT( IAND(ISHFT(offset, -8), 255)  , 1)   !         putc((offset & 0x0000ff00) / 256,fptr)
            write(fptr)INT( IAND(offset, 255)  , 1)              !         putc((offset & 0x000000ff),fptr)
        else
            write(fptr)Z"015300030000000100010000" !         call   WriteHexString(fptr,Z"015300030000000100010000")
        endif
        !!  End of the directory entry
        write(fptr)Z"00000000" !     call WriteHexString(fptr,Z"00000000")
        if ( frames == 3  ) then
            !!  Bits for each colour channel
            write(fptr)Z"000800080008" !         call WriteHexString(fptr,Z"000800080008")
            !!  Minimum value for each component
            write(fptr)Z"000000000000" !         call WriteHexString(fptr,Z"000000000000")
            !!  Maximum value per channel
            write(fptr)Z"00ff00ff00ff" !         call  WriteHexString(fptr,Z"00ff00ff00ff")
            !!  Samples per pixel for each channel
            write(fptr)Z"000100010001" !         call WriteHexString(fptr,Z"000100010001")
        endif
        close(fptr)
        istat = 0 ! success
    end function write_tiff_bigimg

     subroutine write_bigtiff(fname , img)
         character(len=*), intent(in) :: fname
         real, intent(in) :: img(:,:)
         integer :: status, d(3)
         character(len=1, kind=c_char), allocatable :: c_fname(:)
         real :: rn(2)
         real, allocatable, target :: img_ptr(:)
         type(c_ptr) :: cptr
         allocate(c_fname(len(fname)+1), source=trim(fname)//C_NULL_CHAR)
         d(1:2) = size(img)
         d(3) = 1
         rn(1) = minval(img)
         rn(2) = maxval(img)
         allocate(img_ptr(product(d(1:2))))
         img_ptr = pack(img, .true.)
         cptr=c_loc(img_ptr(1))
         status = bigtiff_write(c_fname, d(1), d(2), d(3), rn(1), rn(2), cptr)
         if(status/=0) then
             write(*,'(a)') "simple_tifflib write_bigtiff failed"
         endif
         deallocate(img_ptr)
     end subroutine write_bigtiff

     subroutine write_tiff2(filename, image_data, dims)
         character(len=*), intent(in) :: filename
         integer, intent(in) :: image_data(:), dims(2)
         type(c_ptr) :: out, bufptr
         character(kind=C_char, len=1), allocatable :: filename_c
         integer(1), pointer :: buf(:)
         integer(C_INT32_T) :: imagelength, imagewidth,row, col, n
         integer(C_INT16_T) :: nsamples
         integer :: istat, temp
         filename_c = trim(filename)//C_NULL_CHAR
         out = TIFFOpen(filename_c,"w")
         if (.not. c_associated(out))then  ! out /= C_NULL_PTR
             write(*,'(a)') "Unable to write tif file "
             return
         endif

         imagewidth = dims(1)
         imagelength = dims(2)
         nsamples = 3
         allocate(buf(imagewidth*nsamples))
         istat=TIFFSetField(out, TIFFTAG_IMAGELENGTH, imagelength)
         istat=TIFFSetField(out, TIFFTAG_IMAGEWIDTH, imagewidth)
         istat=TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG)
         istat=TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1) !nsamples)
         istat=TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE) !LZW)
         istat=TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8)
         istat=TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, 8 ) !TIFFDefaultStripSize(out, imagewidth*nsamples));


         bufptr = c_loc(buf(1))

         do row = 1, imagelength
             do col=1, imagewidth
                 temp = image_data(imagewidth*(row-1)+col) ;
                 do n = 1, nsamples
                     buf((col-1)*nsamples + n-1) = INT(IAND(ISHFT( temp,(n-1)*8 -16 ), 255) , 1)
                 end do
             end do

             if (TIFFWriteScanline(out, bufptr, row, 0_2) /= 1 )then
                 write(*,'(a)') "Unable to write a row."
             endif

         end do

         deallocate(buf)
         call TIFFClose(out)
     end subroutine write_tiff2

     subroutine write_tiff3 (filename, arg)
         character(len=*), intent(in) :: filename
         real, intent(in) :: arg(:,:)
         type(c_ptr) :: out, cptr
         character(kind=C_char, len=1), allocatable :: filename_c
         integer(1), pointer :: im_data(:)
         integer :: iret
         integer(C_INT32_T) :: width, height,image_s, imgBytes
         !! Open the TIFF file
         out = TIFFOpen("image.tiff", "w")
         if( .not. c_associated(out) )then
             write(*,'(a)') "Unable to write tif file "
             THROW_HARD("  write_tiff3 failed ")
         endif
         width = size(arg,1)
         height = size(arg,2)
         iret = TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width)
         iret = TIFFSetField(out, TIFFTAG_IMAGELENGTH, height)
         iret = TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1)
         iret = TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 16)
         iret = TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, 1)
         iret = TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
         iret = TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
         iret = TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
         iret = TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);


         !! Write the information to the file
         im_data = float2Dto8bitint(arg,3)
         imgBytes = width*height*3
         cptr = c_loc(im_data(1))
         image_s = TIFFWriteEncodedStrip(out, 0, cptr, imgBytes)
         if( image_s  == -1 ) then
             write(*,'(a)')  "Unable to write tif file "
             stop
         endif

        iret = TIFFWriteDirectory(out)
         call TIFFClose(out)
         deallocate(im_data)
     end subroutine write_tiff3


     function float2Dto8bitint(img, nsamples, scale) result(outimg)
         real, intent(in) :: img(:,:)
         integer, intent(in) :: nsamples
         real, optional , intent(in) :: scale
         integer(1), pointer :: outimg(:)
         integer :: row, col,n, imagelength,imagewidth, temp_i
         real :: temp_r, scale_here
         scale_here = 65535.0
         if(present(scale)) scale_here = scale
         imagewidth = size(img,1)
         imagelength = size(img,2)
         allocate(outimg( imagewidth*imagelength*nsamples ) )
         do row = 1, imagelength
             do col=1, imagewidth
                 temp_r = img(row,col)
                 temp_i = INT(scale_here * temp_r)
                 do n = 1, nsamples
                     outimg((col-1)*nsamples + n-1) = INT(IAND(ISHFT(temp_i,(n-1)*8 -16 ), 255) , 1)
                 end do
             end do
         end do

     end function float2Dto8bitint


 end module simple_tifflib
