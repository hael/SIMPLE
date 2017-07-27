!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Module with type and routine definitions to deal with image file headers
!!
!!  The following file formats are (will be) supported:
!!  - Imagic: http://imagescience.de/formats/index.htm
!!  - Spider: http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
!!  - MRC: http://www2.mrc-lmb.cam.ac.uk/image2000.html
!!
!! This class is based on a class used in CTFFIND4, developed by Alexis Rohou
!! and Nikolaus Grigorieff at Janelia Farm. The below copyright statement therefore
!! needs to be included here:
!! Copyright 2014 Howard Hughes Medical Institute
!! All rights reserved
!! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
!! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
!!
!! Modifications by Cyril Reboul, Michael Eager & Hans Elmlund
!! 
! The SIMPLE code is distributed with the hope that it will be
! useful, but WITHOUT ANY WARRANTY. Redistribution and modification is regulated
! by the GNU General Public License.
! -----------------------------------------------------------------------------!
module simple_imghead
    use simple_imgheadrec, only: int_imgheadrec, real_imgheadrec, char_imgheadrec
    use simple_defs ! singleton
    implicit none

    public :: ImgHead, MrcImgHead, MrcFeiImgHead, SpiImgHead, test_imghead
    private
#include "simple_local_flags.inc"

    integer, parameter         :: NREALS=104
    integer, parameter, public :: dataRbytes   = 1
    integer, parameter, public :: dataRinteger = 2
    integer, parameter, public :: dataRfloat   = 3
    integer, parameter         :: NREMAINS=14
    integer, parameter         :: CLOSE2THEANSWER=43
    integer, parameter         :: MRCHEADSZ=1024
    integer, parameter         :: MRCFEIHEADSZ=128
    integer, parameter         :: NLABL=10


    type ImgHead
        private
        integer                      :: length          !< length (in bytes) of the header
        integer(kind=1), allocatable :: byte_array(:)   !< array of bytes
        logical                      :: exists=.false.  !< to indicate existence
    contains
        procedure :: new
        procedure :: reset2default
        procedure :: setMinimal
        procedure :: print_imghead
        procedure :: read
        procedure :: write
        procedure :: hasLocalEndianess
        procedure, private :: setMachineStamp
        procedure :: bytesPerPix
        procedure :: pixIsSigned
        procedure :: getPixType
        procedure :: pixIsComplex
        procedure :: firstDataByte
        procedure, private :: assign
        generic   :: assignment(=) => assign
        procedure :: getMinPixVal
        procedure :: getMaxPixVal
        procedure :: setMinPixVal
        procedure :: setMaxPixVal
        procedure :: getPixSz
        procedure :: setPixSz
        procedure :: getStackSz
        procedure :: getDims
        procedure :: getDim
        procedure :: setDims
        procedure :: setDim
        procedure :: getLabrec
        procedure :: getLenbyt
        procedure :: getLabbyt
        procedure :: getMaxim
        procedure :: setMaxim
        procedure :: setMode
        procedure :: setRMSD
        procedure :: getIform
        procedure :: getMode
        procedure :: setIform
        procedure :: kill
    end type ImgHead

    type, extends(ImgHead) :: MrcImgHead
    type(int_imgheadrec)  :: nx           !< number of columns (fastest changing in map)
    type(int_imgheadrec)  :: ny           !< number of rows
    type(int_imgheadrec)  :: nz           !< number of sections (slowest changing in map)
    type(int_imgheadrec)  :: mode         !< data type: 0 image: signed 8-bit bytes range -128 to 127
    !!            1 image: 16-bit halfwords
    !!            2 image: 32-bit reals (DEFAULT MODE)
    !!            3 transform: complex 16-bit integers
    !!            4 transform: complex 32-bit reals (THIS WOULD BE THE DEFAULT FT MODE)
    type(int_imgheadrec)  :: nxstart      !< number of first column in map (default = 0)
    type(int_imgheadrec)  :: nystart      !< number of first row in map
    type(int_imgheadrec)  :: nzstart      !< number of first section in map
    type(int_imgheadrec)  :: mx           !< number of intervals along x
    type(int_imgheadrec)  :: my           !< number of intervals along y
    type(int_imgheadrec)  :: mz           !< number of intervals along z
    type(real_imgheadrec) :: cella1       !< cell dims in angstroms
    type(real_imgheadrec) :: cella2
    type(real_imgheadrec) :: cella3
    type(real_imgheadrec) :: cellb1       !< cell angles in degrees
    type(real_imgheadrec) :: cellb2
    type(real_imgheadrec) :: cellb3
    type(int_imgheadrec)  :: mapc         !< axis corresponding to columns  (1,2,3 for x,y,z)
    type(int_imgheadrec)  :: mapr         !< axis corresponding to rows     (1,2,3 for x,y,z)
    type(int_imgheadrec)  :: maps         !< axis corresponding to sections (1,2,3 for x,y,z)
    type(real_imgheadrec) :: dmin         !< minimum density value
    type(real_imgheadrec) :: dmax         !< maximum density value
    type(real_imgheadrec) :: dmean        !< mean density value
    type(int_imgheadrec)  :: ispg         !< space group number 0 or 1 (default=1)
    type(int_imgheadrec)  :: nsymbt       !< number of bytes used for symmetry data (0 or 80)
    type(int_imgheadrec)  :: extra        !< extra space used for anything - 0 by default
    type(int_imgheadrec)  :: originx      !< origin in x used for transforms
    type(int_imgheadrec)  :: originy      !< origin in y used for transforms
    type(int_imgheadrec)  :: originz      !< origin in z used for transforms
    type(char_imgheadrec) :: map          !< character string 'map' to identify file type
    type(int_imgheadrec)  :: machst       !< machine stamp
    type(real_imgheadrec) :: rms          !< rms deviation of map from mean density
    type(int_imgheadrec)  :: nlabl        !< number of labels being used
    type(char_imgheadrec) :: label(NLABL) !< 10 80-character text labels
end type MrcImgHead

type, extends(ImgHead) :: MrcFeiImgHead
    type(real_imgheadrec) :: alpha_tilt        !< Alpha, degrees
    type(real_imgheadrec) :: beta_tilt         !< Beta, degrees
    type(real_imgheadrec) :: xstage            !< Stage X position, meters
    type(real_imgheadrec) :: ystage            !< Stage Y position, meters
    type(real_imgheadrec) :: zstage            !< Stage Z position, meters
    type(real_imgheadrec) :: xshift            !< Image beam shift X, TEM optical units as in general microscope UI
    type(real_imgheadrec) :: yshift            !< Image beam shift Y, TEM optical units as in general microscope UI
    type(real_imgheadrec) :: defocus           !< Defocus as read from microscope in meters
    type(real_imgheadrec) :: exp_time          !< Exposure time in seconds
    type(real_imgheadrec) :: mean_int          !< Image mean pixel intensity
    type(real_imgheadrec) :: na1               !< Not applicable
    type(real_imgheadrec) :: pixsiz            !< Pixel size in the images, meters
    type(real_imgheadrec) :: mag               !< TEM magnification used for recording images
    type(real_imgheadrec) :: ht                !< Value of high tension, volts
    type(real_imgheadrec) :: binning           !< Camera binning used for exposure
    type(real_imgheadrec) :: defocus_app       !< Application’s applied intended defocus, meters
    type(real_imgheadrec) :: na2               !< Unused
    type(real_imgheadrec) :: na3               !< Unused
    type(real_imgheadrec) :: remains(NREMAINS) !< Unused
end type MrcFeiImgHead

type, extends(ImgHead) :: SpiImgHead
    ! All SPIDER image files consist of unformatted, direct access records.
    ! Each record contains NX 4-byte words which are stored as floating point numbers.
    real(kind=4) :: nz        !< (1) Number of slices (planes) in volume (=1 for an image)
    real(kind=4) :: ny        !< (2) Number of rows per slice
    real(kind=4) :: irec      !< (3) Total number of records (including header records) in each image
    !! of a simple image or stacked image file
    real(kind=4) :: iform     !< (5) File type specifier:   1 = 2D image
    !!                             3 = 3D volume
    !!                           -11 = 2D Fourier odd
    !!                           -12 = 2D Fourier even
    !!                           -21 = 3D Fourier odd
    !!                           -22 = 3D Fourier even
    real(kind=4) :: imami     !< (6) Maximum/minimum flag = 0 when the file is created, and = 1 when
    !! the maximum, minimum, average, and standard deviation have been computed
    !! and stored into this header record (see following locations)
    real(kind=4) :: fmax      !< (7) Maximum data value
    real(kind=4) :: fmin      !< (8) Minimum data value
    real(kind=4) :: av        !< (9) Average data value
    real(kind=4) :: sig       !< (10) Standard deviation of data. A value of -1.0 or 0.0 indicates that
    !! SIG has not been computed
    real(kind=4) :: nx        !< (12) Number of pixels (samples) per line
real(kind=4) :: labrec    !< (13) Number of records in file header (label)
real(kind=4) :: iangle      !< (14) Flag that following three tilt angles are present
real(kind=4) :: phi       !< (15) Tilt angle: phi (See note #2 below)
real(kind=4) :: theta     !< (16) Tilt angle: theta
real(kind=4) :: gamma     !< (17) Tilt angle: gamma (also called psi)
real(kind=4) :: xoff      !< (18) X translation
real(kind=4) :: yoff      !< (19) Y translation
real(kind=4) :: zoff      !< (20) Z translation
real(kind=4) :: scale     !< (21) Scale factor
real(kind=4) :: labbyt    !< (22) Total number of bytes in header
real(kind=4) :: lenbyt    !< (23) Record length in bytes
real(kind=4) :: istack    !< (24) istack has a value of 0 in simple 2D or 3D (non-stack) files. In an
!! "image stack" there is one overall stack header followed by a stack of images, in which
!! each image has its own image header. A value of >0 in this position in the overall stack
!! header indicates a stack of images. A value of <0 in this position in the overall stack
!! header indicates an indexed stack of images
real(kind=4) :: maxim     !< (26) only used in the overall header for a stacked image file. There, this position
!! contains the number of the highest image currently used in the stack. This number is
!! updated, if necessary, when an image is added or deleted from the stack
real(kind=4) :: imgnum    !< (27) Position is only used in a stacked image header. There, it contains the number of
!! the current image or zero if this image is unused
real(kind=4) :: lastindx  !< (28) Position is only used in overall header of indexed stacks. There, this position is the
!! highest index location currently in use
real(kind=4) :: kangle    !< (31) Flag that additional rotation angles follow in header.
!! 1 = one additional angle set is present, 2 = two additional angle sets
real(kind=4) :: phi1      !< (32) angle
real(kind=4) :: theta1    !< (33) angle
real(kind=4) :: psi1      !< (34) angle
real(kind=4) :: phi2      !< (35) angle
real(kind=4) :: theta2    !< (36) angle
real(kind=4) :: psi2      !< (37) angle
real(kind=4) :: pixsiz    !< (38) pixel size (Angstroms)
real(kind=4) :: ev        !< (39) electron voltage
real(kind=4) :: proj      !< (40) project number
real(kind=4) :: mic       !< (41) micrograph number
real(kind=4) :: num       !< (42) micrograph window number
real(kind=4) :: glonum    !< (43) global image number
contains
end type SpiImgHead

contains

    ! POLYMORPHIC CONSTRUCTOR

    !>  \brief  create an imghead object
    subroutine new( self, ldim, length )
        class(ImgHead), target, intent(inout) :: self    !< instance
        integer, intent(in), optional         :: ldim(3) !< logical dims of image
        integer, intent(in), optional         :: length  !< length of the header record.
        integer               :: llength, ierr, i, lenbyt, labrec, labbyt
        character(len=STDLEN) :: err
        call self%kill
        select type( self )
        type is( MrcImgHead )
            if( present(length) ) then
                llength = length
            else
                llength = MRCHEADSZ
            endif
            ! allocate byte array
            if( allocated(self%byte_array) )then
                if( size(self%byte_array) .ne. llength ) deallocate(self%byte_array)
            endif
            if( .not. allocated(self%byte_array) )then
                allocate(self%byte_array(llength),stat=ierr,errmsg=err)
                if( ierr .ne. 0 )then
                    write(*,'(a,i0,2a)') '**error(ImgHead::new): memory allocation failed with error ',&
                         ierr, ': ', trim(adjustl(err))
                    stop 'Memory allocation failed; new; simple_imghead'
                endif
            endif
            ! zero the byte array
            self%byte_array = 0
            ! first argument: index_position, i.e. the position of the record within the file header. starting at 1 and incrementing
            ! second argument: byte_position, i.e. the position of the first byte of the record within the header
            call self%nx     %new(1,   1,  self%byte_array)
            call self%ny     %new(2,   5,  self%byte_array)
            call self%nz     %new(3,   9,  self%byte_array)
            call self%mode   %new(4,  13,  self%byte_array)
            call self%nxstart%new(5,  17,  self%byte_array)
            call self%nystart%new(6,  21,  self%byte_array)
            call self%nzstart%new(7,  25,  self%byte_array)
            call self%mx     %new(8,  29,  self%byte_array)
            call self%my     %new(9,  33,  self%byte_array)
            call self%mz     %new(10, 37,  self%byte_array)
            call self%cella1 %new(11, 41,  self%byte_array)
            call self%cella2 %new(12, 45,  self%byte_array)
            call self%cella3 %new(13, 49,  self%byte_array)
            call self%cellb1 %new(14, 53,  self%byte_array)
            call self%cellb2 %new(15, 57,  self%byte_array)
            call self%cellb3 %new(16, 61,  self%byte_array)
            call self%mapc   %new(17, 65,  self%byte_array)
            call self%mapr   %new(18, 69,  self%byte_array)
            call self%maps   %new(19, 73,  self%byte_array)
            call self%dmin   %new(20, 77,  self%byte_array)
            call self%dmax   %new(21, 81,  self%byte_array)
            call self%dmean  %new(22, 85,  self%byte_array)
            call self%ispg   %new(23, 89,  self%byte_array)
            call self%nsymbt %new(24, 93,  self%byte_array)
            call self%extra  %new(25, 97,  self%byte_array)
            call self%originx%new(50, 197, self%byte_array)
            call self%originy%new(51, 201, self%byte_array)
            call self%originz%new(52, 205, self%byte_array)
            call self%map    %new(53, 209, self%byte_array, length=4)
            call self%machst %new(54, 213, self%byte_array)
            call self%rms    %new(55, 217, self%byte_array)
            call self%nlabl  %new(56, 221, self%byte_array)
            do i=1,NLABL
                call self%label(i)%new(57, 225, self%byte_array, length=80)
            end do
            call self%reset2default
            if( present(ldim) )then
                self%nx = INT( ldim(1) ) ! Intel warning
                self%ny = ldim(2)
                self%nz = ldim(3)
                self%mx = ldim(1)
                self%my = ldim(2)
                self%mz = ldim(3)
            endif
        type is( MrcFeiImgHead )
            llength = MRCFEIHEADSZ
            ! allocate byte array
            if( allocated(self%byte_array) )then
                if( size(self%byte_array) .ne. llength ) deallocate(self%byte_array)
            endif
            if( .not. allocated(self%byte_array) )then
                allocate(self%byte_array(llength),stat=ierr,errmsg=err)
                if( ierr .ne. 0 )then
                    write(*,'(a,i0,2a)') '**error(ImgHead::new): memory allocation failed with error ',&
                         ierr, ': ', trim(adjustl(err))
                    stop 'Memory allocation failed; new; simple_imghead'
                endif
            endif
            ! zero the byte array
            self%byte_array = 0
            ! first argument: index_position, i.e. the position of the record within the file header. starting at 1 and incrementing
            ! second argument: byte_position, i.e. the position of the first byte of the record within the header
            call self%alpha_tilt  %new(1,   1,  self%byte_array) !< Alpha, degrees
            call self%beta_tilt   %new(2,   5,  self%byte_array) !< Beta, degrees
            call self%xstage      %new(3,   9,  self%byte_array) !< Stage X position, meters
            call self%ystage      %new(4,  13,  self%byte_array) !< Stage Y position, meters
            call self%zstage      %new(5,  17,  self%byte_array) !< Stage Z position, meters
            call self%xshift      %new(6,  21,  self%byte_array) !< Image beam shift X, TEM optical units as in general microscope UI
            call self%yshift      %new(7,  25,  self%byte_array) !< Image beam shift Y, TEM optical units as in general microscope UI
            call self%defocus     %new(8,  29,  self%byte_array) !< Defocus as read from microscope in meters
            call self%exp_time    %new(9,  33,  self%byte_array) !< Exposure time in seconds
            call self%mean_int    %new(10, 37,  self%byte_array) !< Image mean pixel intensity
            call self%na1         %new(11, 41,  self%byte_array) !< Not applicable
            call self%pixsiz      %new(12, 45,  self%byte_array) !< Pixel size in the images, meters
            call self%mag         %new(13, 49,  self%byte_array) !< TEM magnification used for recording images
            call self%ht          %new(14, 53,  self%byte_array) !< Value of high tension, volts
            call self%binning     %new(15, 57,  self%byte_array) !< Camera binning used for exposure
            call self%defocus_app %new(16, 61,  self%byte_array) !< Application’s applied intended defocus, meters
            call self%na2         %new(17, 65,  self%byte_array) !< Unused
            call self%na3         %new(18, 69,  self%byte_array) !< Unused
            call self%remains(1)  %new(19, 73,  self%byte_array) !< Unused
            call self%remains(2)  %new(20, 77,  self%byte_array) !< Unused
            call self%remains(3)  %new(21, 81,  self%byte_array) !< Unused
            call self%remains(4)  %new(22, 85,  self%byte_array) !< Unused
            call self%remains(5)  %new(23, 89,  self%byte_array) !< Unused
            call self%remains(6)  %new(24, 93,  self%byte_array) !< Unused
            call self%remains(7)  %new(25, 97,  self%byte_array) !< Unused
            call self%remains(8)  %new(26, 101, self%byte_array) !< Unused
            call self%remains(9)  %new(27, 105, self%byte_array) !< Unused
            call self%remains(10) %new(28, 109, self%byte_array) !< Unused
            call self%remains(11) %new(29, 113, self%byte_array) !< Unused
            call self%remains(12) %new(30, 117, self%byte_array) !< Unused
            call self%remains(13) %new(31, 121, self%byte_array) !< Unused
            call self%remains(14) %new(32, 125, self%byte_array) !< Unused
            call self%reset2default
        type is( SpiImgHead )
            if( present(ldim) )then
                call self%reset2default
                call self%setMinimal(ldim, .false., 1.)
            else
                stop 'Need logical dimension (ldim) to create SPIDER header; new; simple_imghead'
            end if
            class DEFAULT
            stop 'Unsupported header type; initImgHead; simple_imghead'
        end select
        self%exists = .true.
        DebugPrint  'created imghead (imghead::new)'
    end subroutine new

!>  \brief  Reset all the values in the header to default values
    subroutine reset2default(self)
        class(ImgHead), intent(inout) :: self
        integer ::  i
        select type(self)
        type is ( MrcImgHead )
            self%nx      = 0
            self%ny      = 0
            self%nz      = 0
            self%mode    = 2
            self%nxstart = 0
            self%nystart = 0
            self%nzstart = 0
            self%mx      = 1
            self%my      = 1
            self%mz      = 1
            self%cella1  = 1.0
            self%cella2  = 1.0
            self%cella3  = 1.0
            self%cellb1  = 90.0
            self%cellb2  = 90.0
            self%cellb3  = 90.0
            self%mapc    = 1
            self%mapr    = 2
            self%maps    = 3
            self%dmin    = 0.0
            self%dmax    = 0.0
            self%dmean   = 0.0
            self%ispg    = 1
            self%nsymbt  = 0
            self%extra   = 0
            self%originx = 0
            self%originy = 0
            self%originz = 0
            self%map     = 'MAP '
            self%machst  = 0
            self%rms     = 1.0
            self%nlabl   = 0
            do i=1,NLABL
                self%label(i) = ' '
            end do
        type is( MrcFeiImgHead )
            self%alpha_tilt  = 0.0    !< Alpha, degrees
            self%beta_tilt   = 0.0    !< Beta, degrees
            self%xstage      = 0.0    !< Stage X position, meters
            self%ystage      = 0.0    !< Stage Y position, meters
            self%zstage      = 0.0    !< Stage Z position, meters
            self%xshift      = 0.0    !< Image beam shift X, TEM optical units as in general microscope UI
            self%yshift      = 0.0    !< Image beam shift Y, TEM optical units as in general microscope UI
            self%defocus     = 0.0    !< Defocus as read from microscope in meters
            self%exp_time    = 0.0    !< Exposure time in seconds
            self%mean_int    = 0.0    !< Image mean pixel intensity
            self%na1         = 0.0    !< Not applicable
            self%pixsiz      = 0.0    !< Pixel size in the images, meters
            self%mag         = 0.0    !< TEM magnification used for recording images
            self%ht          = 0.0    !< Value of high tension, volts
            self%binning     = 0.0    !< Camera binning used for exposure
            self%defocus_app = 0.0    !< Application’s applied intended defocus, meters
            self%na2         = 0.0    !< Unused
            self%na3         = 0.0    !< Unused
            do i=1,NREMAINS
                self%remains(i) = 0.0 !< Unused
            end do
        type is (SpiImgHead)
            self%nz       = 0.
            self%ny       = 0.
            self%irec     = 0.
            self%iform    = 1.
            self%imami    = 0.
            self%fmax     = 0.
            self%fmin     = 0.
            self%av       = 0.
            self%sig      = 0.
            self%nx       = 0.
            self%labrec   = 0.
            self%iangle   = 0.
            self%phi      = 0.
            self%theta    = 0.
            self%gamma    = 0.
            self%xoff     = 0.
            self%yoff     = 0.
            self%zoff     = 0.
            self%scale    = 1.
            self%labbyt   = 0.
            self%lenbyt   = 0.
            self%istack   = 1.
            self%maxim    = 0.
            self%imgnum   = 0.
            self%lastindx = 0.
            self%kangle   = 0.
            self%phi1     = 0.
            self%theta1   = 0.
            self%psi1     = 0.
            self%phi2     = 0.
            self%theta2   = 0.
            self%psi2     = 0.
            self%pixsiz   = 1.
            self%ev       = 300.
            self%proj     = 0.
            self%mic      = 0.
            self%num      = 0.
            self%glonum   = 0.
            class DEFAULT
            stop 'Format not supported; reset2default; simle_imghead'
        end select
    end subroutine reset2default

!>  \brief  if for setting the minimal info in the image header
    subroutine setMinimal( self, ldim, is_ft, smpd )
        class(ImgHead), intent(inout) :: self
        integer, intent(in)           :: ldim(3)
        logical, intent(in)           :: is_ft
        real, intent(in)              :: smpd
        logical                       :: even_dims, is_3d
        integer                       :: lenbyt, labrec, labbyt
        if( ldim(3) == 1 )then
            even_dims = is_even(ldim(:2))
        else
            even_dims = is_even(ldim)
        endif
        is_3d = .false.
        if( ldim(3) > 1 ) is_3d = .true.
        select type(self)
        type is( MrcImgHead )
            call self%setPixSz(smpd)
            call self%setDims(ldim)
        type is( MrcFeiImgHead )
            call self%setPixSz(smpd)
        type is( SpiImgHead )
            call self%setPixSz(smpd)
            lenbyt = ldim(1)*4     ! record length in bytes
            labrec = MRCHEADSZ/lenbyt   ! nr of records in file header (label)
            if( mod(MRCHEADSZ,lenbyt ) /= 0 ) labrec = labrec+1
            labbyt = labrec*lenbyt ! total number of bytes in header
            self%lenbyt = real(lenbyt, kind=4)
            self%labrec = real(labrec, kind=4)
            self%labbyt = real(labbyt, kind=4)
            self%nx     = real(ldim(1), kind=4)
            self%ny     = real(ldim(2), kind=4)
            self%nz     = real(ldim(3), kind=4)
            if( .not. is_3d .and. .not. is_ft )then
                self%iform = 1.
            else if( is_3d .and. .not. is_ft )then
                self%iform = 3.
            else if( .not. is_3d .and. is_ft .and. .not. even_dims )then
                self%iform = -11.
            else if( .not. is_3d .and. is_ft .and. even_dims )then
                self%iform = -12.
            else if( is_3d .and. is_ft .and. .not. even_dims )then
                self%iform = -21.
            else if( is_3d .and. is_ft .and. even_dims )then
                self%iform = -22.
            else
                stop 'undefined file type, setMinimal; simple_imghead'
            endif 
            self%irec = real(ldim(2), kind=4)+self%labrec ! Total number of records (including header records)
            ! in each image of a simple image or stacked image file
            self%sig = -1.         ! Standard deviation of data. A value of -1.0 or 0.0
            ! indicates that SIG has not been computed
            self%scale = 1.        !  Scale factor
            if( ldim(3) > 1 )then  ! istack has a value of 0 in simple 2D or 3D (non-stack) files. In an
                self%istack = 0.   ! "image stack" there is one overall stack header followed by a stack of images, in which
            else                   ! each image has its own image header. A value of >0 in this position in the overall stack
                self%istack = 2.   ! header indicates a stack of images. A value of <0 in this position in the overall stack
            endif                  ! header indicates an indexed stack of images
            if ( ldim(3) > 1 )then ! maxim is only used in the overall header for a stacked image file. There, this position
                self%maxim = 1.    ! contains the number of the highest image currently used in the stack. This number is
            else                   ! updated, if necessary, when an image is added or deleted from the stack
                self%maxim = 0.
            endif
            self%imgnum = 0.   ! only used in a stacked image header, the number of the current image or zero if this image is unused
            class DEFAULT
            stop 'Format not supported; setMinimal; simle_imghead'
        end select

    contains

        !> \brief  to check if all vals in array are even
        pure function is_even( arr ) result( yep )
            integer, intent(in) :: arr(:)
            logical :: yep, test(size(arr))
            integer :: i
            test = .false.
            do i=1,size(arr)
                test(i) = mod(arr(i),2) == 0
            end do
            yep = all(test)
        end function is_even

    end subroutine setMinimal

!>  \brief  print information contained in the header to the terminal
    subroutine print_imghead( self )
        class(ImgHead), intent(in) :: self
        integer :: i
        real    :: smpd(3)
        integer :: bytes_per_pixel
        select type( self )
        type is( MrcImgHead )
            write(*,'(a,3(i0,1x))')     'Number of columns, rows, sections: ', self%nx%get(), self%ny%get(), self%nz%get()
            write(*,'(a,i0)')           'MRC data mode: ',  self%mode%get()
            bytes_per_pixel = self%bytesPerPix()
            write(*,'(a,i0)')           'Bit depth: ',      bytes_per_pixel*8
            smpd = 0.0
            if (self%mx%get() .ne. 0) smpd(1) = self%cella1%get()/self%mx%get()
            if (self%my%get() .ne. 0) smpd(2) = self%cella2%get()/self%my%get()
            if (self%mz%get() .ne. 0) smpd(3) = self%cella3%get()/self%mz%get()
            write(*,'(a,3(f0.3,1x))')   'Pixel size: ', smpd
        type is( MrcFeiImgHead )
            write(*,'(a,1x,f7.3)') 'Pixel size:    ', self%pixsiz%get()
            write(*,'(a,1x,f7.3)') 'Magnification: ', self%mag%get()
            write(*,'(a,1x,f7.3)') 'High tension:  ', self%ht%get()
        type is( SpiImgHead )
            write(*,'(a,3(1x,f7.0))') 'Number of columns, rows, sections: ', self%nx, self%ny, self%nz
            write(*,'(a,1x,f7.0)') 'SPIDER data mode (iform):                                       ', self%iform
            bytes_per_pixel = self%bytesPerPix()
            write(*,'(a,1x,f7.0)') 'Bit depth:                                                       ', real(bytes_per_pixel*8)
            write(*,'(a,1x,f7.3)') 'Pixel size:                                                      ', self%pixsiz
            write(*,'(a,1x,f7.0)') 'Number of slices (planes) in volume (=1 for an image)            ', self%nz
            write(*,'(a,1x,f7.0)') 'Total number of records (including header records) in each image:', self%irec
            write(*,'(a,1x,f7.0)') 'Maximum/minimum flag=1 if stats have been computed and 0 else:   ', self%imami
            write(*,'(a,1x,f7.3)') 'Maximum data value:                                              ', self%fmax
            write(*,'(a,1x,f7.3)') 'Minimum data value:                                              ', self%fmin
            write(*,'(a,1x,f7.3)') 'Average data value:                                              ', self%av
            write(*,'(a,1x,f7.3)') 'Standard deviation. -1. or 0. indicates SIG not computed:        ', self%sig
            write(*,'(a,1x,f7.0)') 'Number of records in file header (label):                        ', self%labrec
            write(*,'(a,1x,f7.0)') 'Flag that following three tilt angles are present:               ', self%iangle
            write(*,'(a,1x,f7.3)') 'Tilt angle: phi (See note #2 below):                             ', self%phi
            write(*,'(a,1x,f7.3)') 'Tilt angle: theta:                                               ', self%theta
            write(*,'(a,1x,f7.3)') 'Tilt angle: gamma (also called psi):                             ', self%gamma
            write(*,'(a,1x,f7.3)') 'X translation:                                                   ', self%xoff
            write(*,'(a,1x,f7.3)') 'Y translation:                                                   ', self%yoff
            write(*,'(a,1x,f7.3)') 'Z translation:                                                   ', self%zoff
            write(*,'(a,1x,f7.3)') 'Scale factor:                                                    ', self%scale
            write(*,'(a,1x,f7.0)') 'Total number of bytes in header:                                 ', self%labbyt
            write(*,'(a,1x,f7.0)') 'Record length in bytes:                                          ', self%lenbyt
            write(*,'(a,1x,f7.0)') 'istack=0 for 3D (non-stack) files and a value > 0 for stacks:    ', self%istack
            write(*,'(a,1x,f7.0)') 'The number of the highest image currently used in the stack:     ', self%maxim
            write(*,'(a,1x,f7.0)') 'Number of current image in stack:                                ', self%imgnum
            write(*,'(a,1x,f7.0)') 'Highest index location currently in use in indexed stacks:       ', self%lastindx
            write(*,'(a,1x,f7.0)') 'Flag that additional rotation angles follow in header:           ', self%kangle
            write(*,'(a,1x,f7.3)') 'Angle:                                                           ', self%phi1
            write(*,'(a,1x,f7.3)') 'Angle:                                                           ', self%theta1
            write(*,'(a,1x,f7.3)') 'Angle:                                                           ', self%psi1
            write(*,'(a,1x,f7.3)') 'Angle:                                                           ', self%phi2
            write(*,'(a,1x,f7.3)') 'Angle:                                                           ', self%theta2
            write(*,'(a,1x,f7.3)') 'Angle:                                                           ', self%psi2
            write(*,'(a,1x,f7.3)') 'Electron voltage:                                                ', self%ev
            write(*,'(a,1x,f7.0)') 'Project number:                                                  ', self%proj
            write(*,'(a,1x,f7.0)') 'Micrograph number:                                               ', self%mic
            write(*,'(a,1x,f7.0)') 'Micrograph window number:                                        ', self%num
            write(*,'(a,1x,f7.0)') 'Global image number:                                             ', self%glonum
            class DEFAULT
            stop 'Format not supported; print; simple_imghead'
        end select
    end subroutine print_imghead

!>  \brief  Read the header data from disk
    subroutine read( self, lun, pos, print_entire )
        class(ImgHead), intent(inout)         :: self
        integer, intent(in)                   :: lun
        integer(kind=8), intent(in), optional :: pos
        logical, intent(in), optional         :: print_entire
        integer                               :: io_status, ppos, i, cnt
        integer                               :: labrec, dim1, alloc_stat
        character(len=512)                    :: io_message
        real(kind=4), allocatable             :: spihed(:)
        ppos = 1
        if( present(pos) ) ppos = pos
        select type( self )
        type is( SpiImgHead )
            allocate(spihed(self%getLabbyt()/4))
            cnt = 0
            do i=ppos,ppos+self%getLabbyt()-1,4
                cnt = cnt+1
                read(unit=lun,pos=i) spihed(cnt)
                if( present(print_entire) )then
                    write(*,*) i, spihed(cnt)
                endif
            end do
            self%nz       = spihed(1)
            DebugPrint ' nz: ', self%nz
            self%ny       = spihed(2)
            DebugPrint ' ny: ', self%ny
            self%irec     = spihed(3)
            self%iform    = spihed(5)
            DebugPrint  ' iform: ', self%iform
            self%imami    = spihed(6)
            self%fmax     = spihed(7)
            self%fmin     = spihed(8)
            self%av       = spihed(9)
            self%sig      = spihed(10)
            self%nx       = spihed(12)
            DebugPrint ' nx: ', self%nx
            self%labrec   = spihed(13)
            DebugPrint ' labrec: ', self%labrec
            self%iangle   = spihed(14)
            self%phi      = spihed(15)
            self%theta    = spihed(16)
            self%gamma    = spihed(17)
            self%xoff     = spihed(18)
            self%yoff     = spihed(19)
            self%zoff     = spihed(20)
            self%scale    = spihed(21)
            self%labbyt   = spihed(22)
            DebugPrint ' labbyt: ', self%labbyt
            self%lenbyt   = spihed(23)
            DebugPrint ' lenbyt: ', self%lenbyt
            self%istack   = spihed(24)
            self%maxim    = spihed(26)
            self%imgnum   = spihed(27)
            self%lastindx = spihed(28)
            self%kangle   = spihed(31)
            self%phi1     = spihed(32)
            self%theta1   = spihed(33)
            self%psi1     = spihed(34)
            self%phi2     = spihed(35)
            self%theta2   = spihed(36)
            self%psi2     = spihed(37)
            self%pixsiz   = spihed(38)
            self%ev       = spihed(39)
            self%proj     = spihed(40)
            self%mic      = spihed(41)
            self%num      = spihed(42)
            self%glonum   = spihed(43)
            deallocate(spihed)
        type is( MrcImgHead )
            read(unit=lun,pos=ppos,iostat=io_status,iomsg=io_message) self%byte_array
            if (io_status .ne. 0) then
                write(*,'(a,i0,2a)') '**error(ImgHead::read): error ', io_status, &
                     ' when reading header bytes from disk: ', trim(io_message)
                stop 'I/O error; read; simple_imghead'
            endif
        type is( MrcFeiImgHead )
            read(unit=lun,pos=ppos,iostat=io_status,iomsg=io_message) self%byte_array
            if (io_status .ne. 0) then
                write(*,'(a,i0,2a)') '**error(ImgHead::read): error ', io_status, &
                     ' when reading header bytes from disk: ', trim(io_message)
                stop 'I/O error; read; simple_imghead'
            endif
            class DEFAULT
            stop 'Format not supported; print; simle_imghead'
        end select
        DebugPrint '(imghead::read) done'
    end subroutine read

!>  \brief write header data to disk
    subroutine write( self, lun, pos )
        class(ImgHead), intent(inout)         :: self
        integer, intent(in)                   :: lun
        integer(kind=8), intent(in), optional :: pos
        integer(kind=8)                       :: ppos, i
        integer                               :: io_status, cnt
        integer                               :: rsz, labrec, dim1, alloc_stat, ldim(3)
        real(kind=4)                          :: spihed(CLOSE2THEANSWER), zero
        zero = 0.
        if( self%exists )then
            ! all good
        else
            stop 'The header that you are trying to write to disk does not exist! write; simple_imghead'
        endif
        ppos = 1
        if( present(pos) ) ppos = pos
        select type( self )
        type is( SpiImgHead )
            spihed = 0.
            spihed(1)  = self%nz
            DebugPrint '::write nz: ', spihed(1)
            spihed(2)  = self%ny
            DebugPrint '::write ny: ', spihed(2)
            spihed(3)  = self%irec
            spihed(5)  = self%iform
            spihed(6)  = self%imami
            spihed(7)  = self%fmax
            spihed(8)  = self%fmin
            spihed(9)  = self%av
            spihed(10) = self%sig
            spihed(12) = self%nx
            DebugPrint '::write nx: ', spihed(12)
            spihed(13) = self%labrec
            DebugPrint '::write labrec: ', spihed(13)
            spihed(14) = self%iangle
            spihed(15) = self%phi
            spihed(16) = self%theta
            spihed(17) = self%gamma
            spihed(18) = self%xoff
            spihed(19) = self%yoff
            spihed(20) = self%zoff
            spihed(21) = self%scale
            spihed(22) = self%labbyt
            DebugPrint '::write labbyt: ', spihed(22)
            spihed(23) = self%lenbyt
            DebugPrint '::write lenbyt: ', spihed(23)
            spihed(24) = self%istack
            spihed(26) = self%maxim
            spihed(27) = self%imgnum
            spihed(28) = self%lastindx
            spihed(31) = self%kangle
            spihed(32) = self%phi1
            spihed(33) = self%theta1
            spihed(34) = self%psi1
            spihed(35) = self%phi2
            spihed(36) = self%theta2
            spihed(37) = self%psi2
            spihed(38) = self%pixsiz
            spihed(39) = self%ev
            spihed(40) = self%proj
            spihed(41) = self%mic
            spihed(42) = self%num
            spihed(43) = self%glonum
            ! write
            cnt = 0
            do i=ppos,ppos+self%getLabbyt()-1,4
                cnt = cnt+1
                if( cnt > CLOSE2THEANSWER )then
                    write(unit=lun,pos=i) zero
                else
                    write(unit=lun,pos=i) spihed(cnt)
                endif
            end do
        type is( MrcImgHead )
            call self%setMachineStamp()
            write(unit=lun,pos=ppos,iostat=io_status) self%byte_array
            if( io_status .ne. 0 )then
                write(*,'(a,i0,a)') '**error(ImgHead::write): error ', io_status, ' when writing header bytes to disk'
                print *, allocated(self%byte_array)
                print *, lun
                stop 'I/O error; write; simple_imghead'
            endif
        type is( MrcFeiImgHead )
            call self%setMachineStamp()
            write(unit=lun,pos=ppos,iostat=io_status) self%byte_array
            if( io_status .ne. 0 )then
                write(*,'(a,i0,a)') '**error(ImgHead::write): error ', io_status, ' when writing header bytes to disk'
                print *, allocated(self%byte_array)
                print *, lun
                stop 'I/O error; write; simple_imghead'
            endif
            class DEFAULT
            stop 'Format not supported; print; simle_imghead'
        end select
    end subroutine write

!>  \brief  Check whether the header was created by a machine with the same endianess as the machine we're on now
!!          If the header doesn't have a machine stamp (or it is 0), we assume it's in the local endianess.
    logical function hasLocalEndianess(self)
        class(ImgHead), intent(in)      ::  self

        !select type(self)
        !    type is (MrcImgHead)
        !        hasLocalEndianess = (getLocalMachineStamp() .eq. self%machst%get()) .or. self%machst%get() .eq. 0
        !    class DEFAULT
        !        call this_program%TerminateWithFatalError('ImgHead::setMachineStamp','Format not supported')
        !end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING THIS TO ALWAYS RETURN TRUE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! It seems most programs ignore the endianess setting in the header, which should basically all be the !
        ! same these days. As such I am setting this function to always return true, the code is being
        ! left in place so that we can revert this back at a later date if required.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING THIS TO ALWAYS RETURN TRUE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        hasLocalEndianess = .true.
    end function hasLocalEndianess

!> \brief   Set the machine stamp
    subroutine setMachineStamp(self)
        class(ImgHead), intent(inout)   ::  self
        select type(self)
        type is (MrcImgHead)
            self%machst = getLocalMachineStamp()
        end select
    end subroutine setMachineStamp

!>  \brief  Return the local machine's "machinestamp", which is endian-specific
    function getLocalMachineStamp()
        integer(kind=4) :: getLocalMachineStamp
        integer(kind=1)            :: machst(4)
        integer(kind=4), parameter :: a0 = 48
        integer(kind=4), parameter :: a1 = 49
        integer(kind=4), parameter :: a2 = 50
        integer(kind=4), parameter :: a3 = 51
        integer(kind=4)            :: i
        character(len=4)           :: ich
        i=a0+a1*256+a2*(256**2)+a3*(256**3) ! = 858927408 (decimal)
        ! = 0011 0011 0011 0010 0011 0001 0011 0000 (binary, little endian)
        ! when this is converted to ASCII characters (1 byte per character,
        ! with the most significant bit always 0) this will give different
        ! results on little- and big-endian machines. For example, '0' in
        ! ASCII has decimal value 48 and bit value 011 0000 '3' in ASCII has
        ! decimal value 51 and bit value 011 0011. Therefore, the value
        ! computed above, when converted to bytes will have the first byte
        ! read off as ASCII character '0' on little-endian and '3' on big-
        ! endian machines

        ! Take the bit pattern over from the 4byte integer to an array of 4 characters (each 1 byte)
        ich = transfer(i,ich)
        if( ich .eq. '0123' )then
            DebugPrint '**debug(getLocalMachineStamp): machine is little-endian (dec/osf, intel, amd ...)'
            !0100 0100
            machst(1)=68
            !0100 0001
            machst(2)=65
            machst(3)=0
            machst(4)=0
        elseif (ich.eq.'3210') then
            DebugPrint '**debug(getLocalMachineStamp): machine is big-endian (sgi, sun, hp, ibm)'
            !0001 0001
            machst(1)=17
            !0001 0001
            machst(2)=17
            machst(3)=0
            machst(4)=0
        else
            DebugPrint '**debug(getLocalMachineStamp): mixed endianity machine (vax)'
            !0010 0010
            machst(1)=34
            !0010 0001
            machst(2)=33
            machst(3)=0
            machst(4)=0
        endif
        ! Convert machst (4 bytes) to a 4-byte integer
        getLocalMachineStamp = transfer(machst,getLocalMachineStamp)
    end function getLocalMachineStamp

!>  \brief  Return the number of bytes per pixel
!!          All SPIDER image files consist of unformatted, direct access records.
!!          Each record contains NX 4-byte words which are stored as floating point numbers.
    integer function bytesPerPix(self)
        class(ImgHead), intent(in) :: self
        bytesPerPix = 4
        select type( self )
        type is( MrcImgHead )
            select case( self%mode%get() )
            case(0)
                bytesPerPix = 1
            case(1,3,6)
                bytesPerPix = 2
            case(2,4)
                bytesPerPix = 4
            end select
        end select
    end function bytesPerPix

!>  \brief  Does the header indicate that pixel density values are stored in unsigned format?
    logical function pixIsSigned(self)
        class(ImgHead), intent(in) :: self
        pixIsSigned = .true.
        select type( self )
        type is( MrcImgHead )
            select case( self%mode%get() )
            case (6,16)
                pixIsSigned = .false.
            end select
        end select
    end function pixIsSigned

!>  \brief  Work out whether pixel data are byte, integer or float
!!          All SPIDER image files consist of unformatted, direct access records.
!!          Each record contains NX 4-byte words which are stored as floating point numbers.
    integer function getPixType(self)
        class(ImgHead), intent(in) :: self
        getPixType = dataRfloat
        select type( self )
        type is( MrcImgHead )
            select case( self%mode%get() )
            case( 0 )
                getPixType = dataRbytes
            case( 1,3,6 )
                getPixType = dataRinteger
            case( 2,4 )
                getPixType = dataRfloat
            end select
        end select
    end function getPixType

!>  \brief  Does the header indicate that pixel density values are complex numbers
    logical function pixIsComplex(self)
        class(ImgHead), intent(in) ::  self
        pixIsComplex = .false.
        select type( self )
        type is( MrcImgHead )
            select case( self%mode%get() )
            case( 3,4 )
                pixIsComplex = .true.
            end select
        end select
    end function pixIsComplex

!>  \brief  Return the index of the first byte containing image data
    integer function firstDataByte( self )
        class(ImgHead), intent(in) :: self
        select type(self)
        type is( MrcImgHead )
            firstDataByte = MRCHEADSZ+1+self%nsymbt%get()
        type is( MrcFeiImgHead )
            firstDataByte = MRCFEIHEADSZ+1
        type is( SpiImgHead )
            firstDataByte = int(self%labbyt)+1
            class DEFAULT
            stop 'Unsupported header type; firstDataByte; simple_imghead'
        end select
    end function firstDataByte

!>  \brief  assignment overloader
    subroutine assign( lhs, rhs )
        class(ImgHead), intent(inout) :: lhs
        class(ImgHead), intent(in)    :: rhs
        integer :: ldim(3)
        select type( lhs )
        type is (SpiImgHead)
            select type( rhs )
            type is( SpiImgHead )
                ldim(1) = int(rhs%nx)
                ldim(2) = int(rhs%ny)
                ldim(3) = int(rhs%nz)
                call lhs%new(ldim=ldim)
                lhs%nz       = rhs%nz
                lhs%ny       = rhs%ny
                lhs%irec     = rhs%irec
                lhs%iform    = rhs%iform
                lhs%imami    = rhs%imami
                lhs%fmax     = rhs%fmax
                lhs%fmin     = rhs%fmin
                lhs%av       = rhs%av
                lhs%sig      = rhs%sig
                lhs%nx       = rhs%nx
                lhs%labrec   = rhs%labrec
                lhs%iangle   = rhs%iangle
                lhs%phi      = rhs%phi
                lhs%theta    = rhs%theta
                lhs%gamma    = rhs%gamma
                lhs%xoff     = rhs%xoff
                lhs%yoff     = rhs%yoff
                lhs%zoff     = rhs%zoff
                lhs%scale    = rhs%scale
                lhs%labbyt   = rhs%labbyt
                lhs%lenbyt   = rhs%lenbyt
                lhs%istack   = rhs%istack
                lhs%maxim    = rhs%maxim
                lhs%imgnum   = rhs%imgnum
                lhs%lastindx = rhs%lastindx
                lhs%kangle   = rhs%kangle
                lhs%phi1     = rhs%phi1
                lhs%theta1   = rhs%theta1
                lhs%psi1     = rhs%psi1
                lhs%phi2     = rhs%phi2
                lhs%theta2   = rhs%theta2
                lhs%psi2     = rhs%psi2
                lhs%pixsiz   = rhs%pixsiz
                lhs%ev       = rhs%ev
                lhs%proj     = rhs%proj
                lhs%mic      = rhs%mic
                lhs%num      = rhs%num
                lhs%glonum   = rhs%glonum
            end select
        type is (MrcImgHead)
            select type( rhs )
            type is( MrcImgHead )
                call lhs%new
                lhs%byte_array = rhs%byte_array
            end select
        type is (MrcFeiImgHead)
            select type( rhs )
            type is( MrcFeiImgHead )
                call lhs%new
                lhs%byte_array = rhs%byte_array
            end select
            class DEFAULT
            stop 'Format not supported (LHS); assign; simle_imghead'
        end select
    end subroutine assign

!>  \brief  Return the minimum pixel value
    real function getMinPixVal(self)
        class(ImgHead), intent(in) :: self
        getMinPixVal = 0.
        select type(self)
        type is (MrcImgHead)
            getMinPixVal = self%dmin%get()
        type is (SpiImgHead)
            getMinPixVal = self%fmin
        end select
    end function getMinPixVal

!>  \brief  Return the maximum pixel value
    real function getMaxPixVal(self)
        class(ImgHead), intent(in) :: self
        getMaxPixVal = 0.
        select type(self)
        type is (MrcImgHead)
            getMaxPixVal = self%dmax%get()
        type is (SpiImgHead)
            getMaxPixVal = self%fmax
        end select
    end function getMaxPixVal

!>  \brief  Set the minimum pixel value
    subroutine setMinPixVal(self,new_value)
        class(ImgHead), intent(inout)  :: self
        real,               intent(in) :: new_value
        select type(self)
        type is (MrcImgHead)
            self%dmin = new_value
        type is (SpiImgHead)
            self%fmin  = new_value
            self%imami = 1
        end select
    end subroutine setMinPixVal

!>  \brief  Set the maximum pixel value
    subroutine setMaxPixVal(self,new_value)
        class(ImgHead), intent(inout) :: self
        real,           intent(in)    :: new_value
        select type(self)
        type is (MrcImgHead)
            self%dmax = new_value
        type is (SpiImgHead)
            self%fmax  = new_value
            self%imami = 1
        end select
    end subroutine setMaxPixVal

!>  \brief  Return the pixel size
    real function getPixSz(self)
        class(ImgHead), intent(in)  ::  self
        select type(self)
        type is (MrcImgHead)
            if (self%mx%get() .ne. 0) then
                getPixSz = self%cella1%get()/self%mx%get()
            endif
        type is (SpiImgHead)
            getPixSz = self%pixsiz
        type is (MrcFeiImgHead)
            getPixSz = self%pixsiz%get()
        end select
    end function getPixSz

!>  \brief  Set the pixel size
    subroutine setPixSz( self, smpd )
        class(ImgHead), intent(inout) :: self
        real, intent(in)              :: smpd
        select type(self)
        type is (MrcImgHead)
            self%cella1 = smpd*self%mx%get()
            self%cella2 = smpd*self%my%get()
            self%cella3 = smpd*self%mz%get()
        type is (SpiImgHead)
            self%pixsiz = smpd
        type is (MrcFeiImgHead)
            self%pixsiz = smpd
        end select
    end subroutine setPixSz

!>  \brief  Return the number of 2D images in the stack
    integer function getStackSz(self) result(stack_size)
        class(ImgHead), intent(in) ::  self
        stack_size = 0
        select type(self)
        type is (MrcImgHead)
            stack_size = self%nz%get()
        type is (SpiImgHead)
            stack_size = int(self%maxim)
        end select
    end function getStackSz

!>  \brief  Return the dims stored in the header
    function getDims(self) result( dims )
        class(ImgHead), intent(in) ::  self
        integer ::  dims(3)
        dims = 0
        select type(self)
        type is (MrcImgHead)
            dims = [self%nx%get(),self%ny%get(),self%nz%get()]
        type is (SpiImgHead)
            dims = int([self%nx,self%ny,self%nz])
        end select
    end function getDims

!>  \brief  Return one of the dims stored in the header
    function getDim( self, which_dim ) result( dim )
        class(ImgHead), intent(in) :: self
        integer, intent(in)        :: which_dim
        integer :: dim
        dim = 0
        select type( self )
        type is( MrcImgHead )
            select case( which_dim )
            case (1)
                dim = self%nx%get()
            case (2)
                dim = self%ny%get()
            case (3)
                dim = self%nz%get()
            case DEFAULT
                stop 'Dimension should be 1, 2 or 3; getDim; simle_imghead'
            end select
        type is( SpiImgHead )
            select case( which_dim )
            case (1)
                dim = int(self%nx)
            case (2)
                dim = int(self%ny)
            case (3)
                dim = int(self%nz)
            case DEFAULT
                stop 'Dimension should be 1, 2 or 3; getDim; simle_imghead'
            end select
        end select
    end function getDim

!>  \brief  set the dims in the header
    subroutine setDims( self, ldim )
        class(ImgHead), intent(inout) :: self
        integer, intent(in)           :: ldim(3)
        select type( self )
        type is( MrcImgHead )
            self%nx = ldim(1)
            self%ny = ldim(2)
            self%nz = ldim(3)
            self%mx = ldim(1)
            self%my = ldim(2)
            self%mz = ldim(3)
        type is( SpiImgHead )
            self%nx = ldim(1)
            self%ny = ldim(2)
            self%nz = ldim(3)
        end select
    end subroutine setDims

!>  \brief  set the dims in the header
    subroutine setDim( self, which_dim, d )
        class(ImgHead), intent(inout) :: self
        integer, intent(in)           :: which_dim, d
        if( d < 1 )then
            write(*,'(a,1x,f7.0)') 'Dimension: ', real(d)
            write(*,*) 'Trying to set image dimension that is nonconforming; setDim; simple_imghead'
            stop
        endif
        select type( self )
        type is( MrcImgHead )
            select case( which_dim )
            case(1)
                self%nx = d
                self%mx = d
            case(2)
                self%ny = d
                self%my = d
            case(3)
                self%nz = d
                self%mz = d
            case DEFAULT
                stop 'not a valid dimension; setDim; simple_imghead'
            end select
        type is( SpiImgHead )
            select case( which_dim )
            case(1)
                self%nx = d
            case(2)
                self%ny = d
            case(3)
                self%nz = d
            case DEFAULT
                stop 'not a valid dimension; setDim; simple_imghead'
            end select
        end select
    end subroutine setDim

!>  \brief  Return the number of records in file header
    function getLabrec( self ) result( labrec )
        class(ImgHead), intent(in) :: self
        integer :: labrec
        labrec = 0
        select type( self )
        type is( SpiImgHead )
            labrec = int(self%labrec)
        end select
    end function getLabrec

!>  \brief  Return the record length in bytes
    function getLenbyt( self ) result( lenbyt )
        class(ImgHead), intent(in) :: self
        integer :: lenbyt
        lenbyt = 0
        select type( self )
        type is( SpiImgHead )
            lenbyt = int(self%lenbyt)
        end select
    end function getLenbyt

!>  \brief  Return the record length in bytes
    function getLabbyt( self ) result( labbyt )
        class(ImgHead), intent(in) :: self
        integer :: labbyt
        labbyt = 0
        select type( self )
        type is( SpiImgHead )
            labbyt = int(self%lenbyt)
        end select
    end function getLabbyt

!>  \brief  Return the maximum nr of images in stack
    function getMaxim( self ) result( maxim )
        class(ImgHead), intent(in) :: self
        integer :: maxim
        maxim = 0
        select type( self )
        type is( SpiImgHead )
            maxim = int(self%maxim)
        type is( MrcImgHead)
            maxim = self%nz%get()
        end select
    end function getMaxim

!>  \brief  Set the maximum nr of images in stack
    subroutine setMaxim( self, maxim )
        class(ImgHead), intent(inout) :: self
        integer, intent(in) :: maxim
        select type( self )
        type is( SpiImgHead )
            self%maxim = maxim
        end select
    end subroutine setMaxim

!>  \brief  Return the image format tag
    function getIform( self ) result( iform )
        class(ImgHead), intent(in) :: self
        integer :: iform
        iform = 1
        select type( self )
        type is( SpiImgHead )
            iform = int(self%iform)
        end select
    end function getIform

!>  \brief  Return the image format tag
    function getMode( self ) result( mode )
        class(ImgHead), intent(in) :: self
        integer :: mode
        mode = 2
        select type( self )
        type is( MrcImgHead )
            mode = self%mode%get()
        end select
    end function getMode

!>  \brief  Set the maximum nr of images in stack
    subroutine setIform( self, iform )
        class(ImgHead), intent(inout) :: self
        integer, intent(in)           :: iform
        select type( self )
        type is( SpiImgHead )
            self%iform = iform
        end select
    end subroutine setIform

!>  \brief  Set the image format tag
    subroutine setMode( self, mode )
        class(ImgHead), intent(inout) :: self
        integer, intent(in)           :: mode
        select type( self )
        type is( MrcImgHead )
            self%mode = mode
        end select
    end subroutine setMode

!>  \brief  Set the root-mean-square deviation
    subroutine setRMSD( self, RMSD )
        class(ImgHead), intent(inout) :: self
        real,           intent(in)    :: RMSD
        select type( self )
        type is( MrcImgHead )
            self%rms = RMSD
        end select
    end subroutine setRMSD

! POLYMORPHIC DESTRUCTOR

    subroutine kill(self)
        class(ImgHead), intent(inout) :: self
        integer :: i
        if( self%exists )then
            select type(self)
            type is (MrcImgHead)
                call self%nx     %kill
                call self%ny     %kill
                call self%nz     %kill
                call self%mode   %kill
                call self%nxstart%kill
                call self%nystart%kill
                call self%nzstart%kill
                call self%mx     %kill
                call self%my     %kill
                call self%mz     %kill
                call self%cella1 %kill
                call self%cella2 %kill
                call self%cella3 %kill
                call self%cellb1 %kill
                call self%cellb2 %kill
                call self%cellb3 %kill
                call self%mapc   %kill
                call self%mapr   %kill
                call self%maps   %kill
                call self%dmin   %kill
                call self%dmax   %kill
                call self%dmean  %kill
                call self%ispg   %kill
                call self%nsymbt %kill
                call self%extra  %kill
                call self%originx%kill
                call self%originy%kill
                call self%originz%kill
                call self%map    %kill
                call self%machst %kill
                call self%rms    %kill
                call self%nlabl  %kill
                do i=1,NLABL
                    call self%label(i)%kill
                enddo
                if(allocated(self%byte_array)) deallocate(self%byte_array)
            type is (MrcFeiImgHead)
                call self%alpha_tilt  %kill
                call self%beta_tilt   %kill
                call self%xstage      %kill
                call self%ystage      %kill
                call self%zstage      %kill
                call self%xshift      %kill
                call self%yshift      %kill
                call self%defocus     %kill
                call self%exp_time    %kill
                call self%mean_int    %kill
                call self%na1         %kill
                call self%pixsiz      %kill
                call self%mag         %kill
                call self%ht          %kill
                call self%binning     %kill
                call self%defocus_app %kill
                call self%na2         %kill
                call self%na3         %kill
                do i=1,NREMAINS
                    call self%remains(i)%kill
                end do
                if(allocated(self%byte_array)) deallocate(self%byte_array)
            type is (SpiImgHead)
                return
                class DEFAULT
                stop 'Unsupported header type; kill; simple_imghead'
            end select
            self%exists = .false.
        endif
    end subroutine kill

    subroutine test_imghead
        use simple_filehandling
        class(ImgHead), allocatable :: hed, hed2
        integer :: recsz, funit, dims(3), dims2(3)
        write(*,'(a)') '**info(simple_imghead_unit_test): testing read/write capabilities'
        allocate(SpiImgHead :: hed, hed2 )
        call hed%new([120,120,1])
        call hed2%new([120,120,1])
        recsz = 120*4
        funit = get_fileunit()
        open(unit=funit,access='STREAM',file='test_imghed.spi',&
             &action='READWRITE',status='UNKNOWN',&
#ifndef INTEL
            &convert='BIG_ENDIAN'&
#endif
            &)
        call hed%write(funit)
        call hed2%read(funit)
        close(funit)
        write(*,*) '>>> PRINTING HEADER THAT WAS WRITTEN TO DISK'
        call hed%print_imghead
        write(*,*) ''
        write(*,*) '*************************************************'
        write(*,*) ''
        write(*,*) '>>> PRINTING HEADER THAT WAS READ FROM DISK'
        call hed2%print_imghead
        if( all(hed%getDims() == hed2%getDims()) )then
            ! all good
        else
            stop 'test_imghedfailed; simple_imghed'
        endif
        write(*,'(a)') 'SIMPLE_IMGHEAD_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_imghead

end module simple_imghead
