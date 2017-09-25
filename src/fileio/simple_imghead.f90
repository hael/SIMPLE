! type and routine definitions to deal with image file headers supporting
!  - Spider: http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
!  - MRC: http://www2.mrc-lmb.cam.ac.uk/image2000.html
! file formats
! This class is based on a class used in CTFFIND4, developed by Alexis Rohou
! and Nikolaus Grigorieff at Janelia Farm. The below copyright statement therefore
! needs to be included here:
! Copyright 2014 Howard Hughes Medical Institute
! All rights reserved
! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
! Modifications by Cyril Reboul, Michael Eager & Hans Elmlund

module simple_imghead
use simple_defs
use simple_syslib, only: alloc_errchk, file_exists
use simple_fileio, only: fopen, fileio_errmsg, fclose, fname2format
use simple_strings, only: int2str
implicit none

public :: ImgHead, MrcImgHead, SpiImgHead, test_imghead, find_ldim_nptcls, has_ldim_nptcls
private
#include "simple_local_flags.inc"

integer, parameter         :: NREALS          = 104
integer, parameter, public :: dataRbytes      = 1
integer, parameter, public :: dataRinteger    = 2
integer, parameter, public :: dataRfloat      = 3
integer, parameter         :: NREMAINS        = 14
integer, parameter         :: CLOSE2THEANSWER = 43
integer, parameter         :: MRCHEADSZ       = 1024
integer, parameter         :: NLABL           = 10

type ImgHead
    private
    integer                      :: length         !< length (in bytes) of the header
    integer(kind=1), allocatable :: byte_array(:)  !< array of bytes
    logical                      :: exists=.false. !< to indicate existence
contains
    ! polymorphic constructor
    procedure          :: new
    procedure          :: reset2default
    procedure          :: setMinimal
    ! I/O
    procedure          :: print_imghead
    procedure          :: read
    procedure          :: write
    ! byte array conversions
    procedure, private :: transfer_obj2byte_array
    procedure, private :: transfer_byte_array2obj
    ! getters/setters
    procedure          :: bytesPerPix
    procedure          :: pixIsSigned
    procedure          :: getPixType
    procedure          :: pixIsComplex
    procedure          :: firstDataByte
    procedure          :: getMinPixVal
    procedure          :: getMaxPixVal
    procedure          :: setMinPixVal
    procedure          :: setMaxPixVal
    procedure          :: getPixSz
    procedure          :: setPixSz
    procedure          :: getStackSz
    procedure          :: getDims
    procedure          :: getDim
    procedure          :: setDims
    procedure          :: setDim
    procedure          :: getLabrec
    procedure          :: getLenbyt
    procedure          :: getLabbyt
    procedure          :: getMaxim
    procedure          :: setMaxim
    procedure          :: setMode
    procedure          :: setRMSD
    procedure          :: getIform
    procedure          :: getMode
    procedure          :: setIform
    ! destructor
    procedure          :: kill
end type ImgHead

type, extends(ImgHead) :: MrcImgHead
    integer   :: nx           !< number of columns (fastest changing in map)
    integer   :: ny           !< number of rows
    integer   :: nz           !< number of sections (slowest changing in map)
    integer   :: mode         !< data type: 0 image: signed 8-bit bytes range -128 to 127
    !!            1 image: 16-bit halfwords
    !!            2 image: 32-bit reals (DEFAULT MODE)
    !!            3 transform: complex 16-bit integers
    !!            4 transform: complex 32-bit reals (THIS WOULD BE THE DEFAULT FT MODE)
    integer   :: nxstart      !< number of first column in map (default = 0)
    integer   :: nystart      !< number of first row in map
    integer   :: nzstart      !< number of first section in map
    integer   :: mx           !< number of intervals along x
    integer   :: my           !< number of intervals along y
    integer   :: mz           !< number of intervals along z
    real      :: cella1       !< cell dims in angstroms
    real      :: cella2
    real      :: cella3
    real      :: cellb1       !< cell angles in degrees
    real      :: cellb2
    real      :: cellb3
    integer   :: mapc         !< axis corresponding to columns  (1,2,3 for x,y,z)
    integer   :: mapr         !< axis corresponding to rows     (1,2,3 for x,y,z)
    integer   :: maps         !< axis corresponding to sections (1,2,3 for x,y,z)
    real      :: dmin         !< minimum density value
    real      :: dmax         !< maximum density value
    real      :: dmean        !< mean density value
    integer   :: ispg         !< space group number 0 or 1 (default=1)
    integer   :: nsymbt       !< number of bytes used for symmetry data (0 or 80)
    integer   :: extra        !< extra space used for anything - 0 by default
    integer   :: originx      !< origin in x used for transforms
    integer   :: originy      !< origin in y used for transforms
    integer   :: originz      !< origin in z used for transforms
    character(len=4)  :: map  !< character string 'map' to identify file type
    integer   :: machst       !< machine stamp
    real      :: rms          !< rms deviation of map from mean density
    integer   :: nlabl        !< number of labels being used
    character(len=80) :: label(NLABL) !< 10 80-character text labels
end type MrcImgHead

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
    real(kind=4) :: iangle    !< (14) Flag that following three tilt angles are present
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
end type SpiImgHead

contains

    ! polymorphic constructor

    !>  \brief  create an imghead object
    subroutine new( self, ldim, length )
        class(ImgHead), target, intent(inout) :: self    !< instance
        integer, optional,      intent(in)    :: ldim(3) !< logical dims of image
        integer, optional,      intent(in)    :: length  !< length of the header record.
        integer               :: llength
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
                    allocate(self%byte_array(llength),stat=alloc_stat)
                    if(alloc_stat .ne. 0) &
                         call alloc_errchk("simple_imghead::new byte_array ", alloc_stat)
                endif
                ! zero the byte array
                self%byte_array = 0
                ! insert default header values
                call self%reset2default
                if( present(ldim) )then
                    self%nx = ldim(1)
                    self%ny = ldim(2)
                    self%nz = ldim(3)
                    self%mx = ldim(1)
                    self%my = ldim(2)
                    self%mz = ldim(3)
                endif
            type is( SpiImgHead )
                if( present(ldim) )then
                    ! insert default header values
                    call self%reset2default
                    ! set minimal header
                    call self%setMinimal(ldim, .false., 1.)
                else
                    stop 'Need logical dimension (ldim) to create SPIDER header; new; simple_imghead'
                end if
            class DEFAULT
                stop 'Unsupported header type; initImgHead; simple_imghead'
        end select
        self%exists = .true.
        DebugPrint  'created imghead ::new '
    end subroutine new

    ! I/O

    !>  \brief  print information contained in the header to the terminal
    subroutine print_imghead( self )
        class(ImgHead), intent(in) :: self
        real    :: smpd(3)
        integer :: bytes_per_pixel
        select type( self )
            type is( MrcImgHead )
                write(*,'(a,3(i0,1x))')     'Number of columns, rows, sections: ', self%nx, self%ny, self%nz
                write(*,'(a,i0)')           'MRC data mode: ',  self%mode
                bytes_per_pixel = self%bytesPerPix()
                write(*,'(a,i0)')           'Bit depth: ',      bytes_per_pixel*8
                smpd = 0.0
                if (self%mx .ne. 0) smpd(1) = self%cella1/self%mx
                if (self%my .ne. 0) smpd(2) = self%cella2/self%my
                if (self%mz .ne. 0) smpd(3) = self%cella3/self%mz
                write(*,'(a,3(f0.3,1x))')   'Pixel size: ', smpd
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
        class(ImgHead),            intent(inout) :: self
        integer,                   intent(in)    :: lun
        integer(kind=8), optional, intent(in)    :: pos
        logical,         optional, intent(in)    :: print_entire
        integer(kind=8)                   ::  ppos, i, cnt
        integer                   :: io_status
        character(len=512)        :: io_message
        real(kind=4), allocatable :: spihed(:)
        ppos = 1
        if( present(pos) ) ppos = pos
        select type( self )
            type is( SpiImgHead )
                allocate(spihed(self%getLabbyt()/4), stat=alloc_stat)
                if(alloc_stat/=0) call alloc_errchk("In simple_imghead::read spihed ", alloc_stat)
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
                call fileio_errmsg(" simple_imghead::read header bytes to disk , message "&
                    //trim(io_message),io_status)
                call self%transfer_byte_array2obj
            class DEFAULT
                stop 'Format not supported; print; simle_imghead'
        end select
        DebugPrint '(imghead::read) done'
    end subroutine read

    !>  \brief write header data to disk
    subroutine write( self, lun, pos )
        class(ImgHead),            intent(inout) :: self
        integer,                   intent(in)    :: lun
        integer(kind=8), optional, intent(in)    :: pos
        integer(kind=8) :: ppos, i
        integer         :: io_status, cnt
        real(kind=4)    :: spihed(CLOSE2THEANSWER), zero
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
            call self%transfer_obj2byte_array
            write(unit=lun,pos=ppos,iostat=io_status) self%byte_array
            call fileio_errmsg(" simple_imghead::write  writing header bytes to disk , unit "//int2str(lun),io_status)
           
        class DEFAULT
            stop 'Format not supported; print; simle_imghead'
        end select
    end subroutine write

    ! byte array conversions

    subroutine transfer_obj2byte_array( self )
        class(ImgHead), intent(inout) :: self !! ICE error #6780: A dummy argument with the INTENT(IN) attribute shall not be defined nor become undefined.
        integer :: i
        select type( self )
            type is( MrcImgHead )
                self%byte_array(1:4)     = transfer(self%nx,      self%byte_array(1:4))
                self%byte_array(5:8)     = transfer(self%ny,      self%byte_array(5:8))
                self%byte_array(9:12)    = transfer(self%nz,      self%byte_array(9:12))
                self%byte_array(13:16)   = transfer(self%mode,    self%byte_array(13:16))
                self%byte_array(17:20)   = transfer(self%nxstart, self%byte_array(17:20))
                self%byte_array(21:24)   = transfer(self%nystart, self%byte_array(21:24))
                self%byte_array(25:28)   = transfer(self%nzstart, self%byte_array(25:28))
                self%byte_array(29:32)   = transfer(self%mx,      self%byte_array(29:32))
                self%byte_array(33:36)   = transfer(self%my,      self%byte_array(33:36))
                self%byte_array(37:40)   = transfer(self%mz,      self%byte_array(37:40))
                self%byte_array(41:44)   = transfer(self%cella1,  self%byte_array(41:44))
                self%byte_array(45:48)   = transfer(self%cella2,  self%byte_array(45:48))
                self%byte_array(49:52)   = transfer(self%cella3,  self%byte_array(49:52))
                self%byte_array(53:56)   = transfer(self%cellb1,  self%byte_array(53:56))
                self%byte_array(57:60)   = transfer(self%cellb2,  self%byte_array(57:60))
                self%byte_array(61:64)   = transfer(self%cellb3,  self%byte_array(61:64))
                self%byte_array(65:68)   = transfer(self%mapc,    self%byte_array(65:68))
                self%byte_array(69:72)   = transfer(self%mapr,    self%byte_array(69:72))
                self%byte_array(73:76)   = transfer(self%maps,    self%byte_array(73:76))
                self%byte_array(77:80)   = transfer(self%dmin,    self%byte_array(77:80))
                self%byte_array(81:84)   = transfer(self%dmax,    self%byte_array(81:84))
                self%byte_array(85:88)   = transfer(self%dmean,   self%byte_array(85:88))
                self%byte_array(89:92)   = transfer(self%ispg,    self%byte_array(89:92))
                self%byte_array(93:96)   = transfer(self%nsymbt,  self%byte_array(93:96))
                self%byte_array(97:100)  = transfer(self%extra,   self%byte_array(97:100))
                self%byte_array(197:200) = transfer(self%originx, self%byte_array(197:200))
                self%byte_array(201:204) = transfer(self%originy, self%byte_array(201:204))
                self%byte_array(205:208) = transfer(self%originz, self%byte_array(205:208))
                self%byte_array(209:212) = transfer(self%map,     self%byte_array(209:212))
                self%byte_array(213:216) = transfer(self%machst,  self%byte_array(213:216))
                self%byte_array(217:220) = transfer(self%rms,     self%byte_array(217:220))
                self%byte_array(221:224) = transfer(self%nlabl,   self%byte_array(221:224))
                do i=1,NLABL
                    self%byte_array(225:304) = transfer(self%label(i), self%byte_array(225:304))
                end do
        end select
    end subroutine transfer_obj2byte_array

    subroutine transfer_byte_array2obj( self )
        class(ImgHead), intent(inout) :: self !! ICE error #6780: A dummy argument with the INTENT(IN) attribute shall not be defined nor become undefined.
        integer :: i
        select type( self )
            type is( MrcImgHead )
                self%nx      = transfer(self%byte_array(1:4),     self%nx)
                self%ny      = transfer(self%byte_array(5:8),     self%ny)
                self%nz      = transfer(self%byte_array(9:12),    self%nz)
                self%mode    = transfer(self%byte_array(13:16),   self%mode)
                self%nxstart = transfer(self%byte_array(17:20),   self%nxstart)
                self%nystart = transfer(self%byte_array(21:24),   self%nystart)
                self%nzstart = transfer(self%byte_array(25:28),   self%nzstart)
                self%mx      = transfer(self%byte_array(29:32),   self%mx)
                self%my      = transfer(self%byte_array(33:36),   self%my)
                self%mz      = transfer(self%byte_array(37:40),   self%mz)
                self%cella1  = transfer(self%byte_array(41:44),   self%cella1)
                self%cella2  = transfer(self%byte_array(45:48),   self%cella2)
                self%cella3  = transfer(self%byte_array(49:52),   self%cella3)
                self%cellb1  = transfer(self%byte_array(53:56),   self%cellb1)
                self%cellb2  = transfer(self%byte_array(57:60),   self%cellb2)
                self%cellb3  = transfer(self%byte_array(61:64),   self%cellb3)
                self%mapc    = transfer(self%byte_array(65:68),   self%mapc)
                self%mapr    = transfer(self%byte_array(69:72),   self%mapr)
                self%maps    = transfer(self%byte_array(73:76),   self%maps)
                self%dmin    = transfer(self%byte_array(77:80),   self%dmin)
                self%dmax    = transfer(self%byte_array(81:84),   self%dmax)
                self%dmean   = transfer(self%byte_array(85:88),   self%dmean)
                self%ispg    = transfer(self%byte_array(89:92),   self%ispg)
                self%nsymbt  = transfer(self%byte_array(93:96),   self%nsymbt)
                self%extra   = transfer(self%byte_array(97:100),  self%extra)
                self%originx = transfer(self%byte_array(197:200), self%originx)
                self%originy = transfer(self%byte_array(201:204), self%originy)
                self%originz = transfer(self%byte_array(205:208), self%originz)
                self%map     = transfer(self%byte_array(209:212), self%map)
                self%machst  = transfer(self%byte_array(213:216), self%machst)
                self%rms     = transfer(self%byte_array(217:220), self%rms)
                self%nlabl   = transfer(self%byte_array(221:224), self%nlabl)
                do i=1,NLABL
                    self%label(i) = transfer(self%byte_array(225:304), self%label(i))
                end do
        end select
    end subroutine transfer_byte_array2obj

    ! getters/setters

    !>  \brief  Reset all the values in the header to default values
    subroutine reset2default( self )
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
        integer,        intent(in)    :: ldim(3)
        logical,        intent(in)    :: is_ft
        real,           intent(in)    :: smpd
        logical :: even_dims, is_3d
        integer :: lenbyt, labrec, labbyt
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
            type is( SpiImgHead )
                call self%setPixSz(smpd)
                lenbyt = ldim(1)*4     ! record length in bytes
                labrec = MRCHEADSZ/lenbyt   ! nr of records in file header (label)
                if( mod(MRCHEADSZ,lenbyt ) /= 0 ) labrec = labrec + 1
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

    !>  \brief  Return the number of bytes per pixel
    !!          All SPIDER image files consist of unformatted, direct access records.
    !!          Each record contains NX 4-byte words which are stored as floating point numbers.
    integer function bytesPerPix(self)
        class(ImgHead), intent(in) :: self
        bytesPerPix = 4
        select type( self )
            type is( MrcImgHead )
                select case( self%mode )
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
            select case( self%mode )
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
            select case( self%mode )
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
            select case( self%mode )
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
            firstDataByte = MRCHEADSZ+1+self%nsymbt
        type is( SpiImgHead )
            firstDataByte = int(self%labbyt)+1
            class DEFAULT
            stop 'Unsupported header type; firstDataByte; simple_imghead'
        end select
    end function firstDataByte

    !>  \brief  Return the minimum pixel value
    real function getMinPixVal(self)
        class(ImgHead), intent(in) :: self
        getMinPixVal = 0.
        select type(self)
            type is (MrcImgHead)
                getMinPixVal = self%dmin
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
                getMaxPixVal = self%dmax
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
        getPixSz = 0.
        select type(self)
            type is (MrcImgHead)
                if( self%mx .ne. 0 ) getPixSz = self%cella1/self%mx
            type is (SpiImgHead)
                getPixSz = self%pixsiz
        end select
    end function getPixSz

    !>  \brief  Set the pixel size
    subroutine setPixSz( self, smpd )
        class(ImgHead), intent(inout) :: self
        real,           intent(in)    :: smpd
        select type(self)
            type is (MrcImgHead)
                self%cella1 = smpd*self%mx
                self%cella2 = smpd*self%my
                self%cella3 = smpd*self%mz
            type is (SpiImgHead)
                self%pixsiz = smpd
        end select
    end subroutine setPixSz

    !>  \brief  Return the number of 2D images in the stack
    integer function getStackSz(self) result(stack_size)
        class(ImgHead), intent(in) ::  self
        stack_size = 0
        select type(self)
            type is (MrcImgHead)
                stack_size = self%nz
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
                dims = [self%nx,self%ny,self%nz]
            type is (SpiImgHead)
                dims = int([self%nx,self%ny,self%nz])
        end select
    end function getDims

    !>  \brief  Return one of the dims stored in the header
    function getDim( self, which_dim ) result( dim )
        class(ImgHead), intent(in) :: self
        integer,        intent(in) :: which_dim
        integer :: dim
        dim = 0
        select type( self )
            type is( MrcImgHead )
                select case( which_dim )
                    case (1)
                        dim = self%nx
                    case (2)
                        dim = self%ny
                    case (3)
                        dim = self%nz
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
        integer,        intent(in)    :: ldim(3)
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
        integer,        intent(inout)    :: which_dim, d
        if( d < 1 )then
            write(*,'(a,1x,f7.0)') 'Dimension: ', real(d)
            write(*,*) 'Trying to set image dimension that is nonconforming; setDim; simple_imghead'
            stop
        endif
        DebugPrint 'imghead::setDim check type of ImgHead'
        
        select type( self )
        class is( MrcImgHead )
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
        class is( SpiImgHead )
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
        DebugPrint 'imghead::setDim check type of ImgHead done'
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
                maxim = self%nz
        end select
    end function getMaxim

    !>  \brief  Set the maximum nr of images in stack
    subroutine setMaxim( self, maxim )
        class(ImgHead), intent(inout) :: self
        integer,        intent(in) :: maxim
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
                mode = self%mode
        end select
    end function getMode

    !>  \brief  Set the maximum nr of images in stack
    subroutine setIform( self, iform )
        class(ImgHead), intent(inout) :: self
        integer,        intent(in)    :: iform
        select type( self )
            type is( SpiImgHead )
                self%iform = iform
        end select
    end subroutine setIform

    !>  \brief  Set the image format tag
    subroutine setMode( self, mode )
        class(ImgHead), intent(inout) :: self
        integer,        intent(in)    :: mode
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


    !>  \brief is for gettign a part of the info in a MRC image header
    subroutine get_mrcfile_info( fname, ldim, form, smpd, doprint )
        character(len=*), intent(in)  :: fname
        character(len=1), intent(in)  :: form
        integer,          intent(out) :: ldim(3)
        real,             intent(out) :: smpd
        logical,          intent(in)  :: doprint
        class(imghead), allocatable   :: hed
        integer :: filnum,ios
        ldim = 0
        smpd = 0.
        if( file_exists(fname) )then
            select case(form)
                case('M')
                    allocate(MrcImgHead :: hed)
                    call hed%new
                    call fopen(filnum, status='OLD', action='READ', file=fname, access='STREAM', iostat=ios)
                    call fileio_errmsg(" get_mrcfile_info fopen error "//trim(fname),ios)
                    call hed%read(filnum)
                    call fclose(filnum,errmsg=" get_mrcfile_info close error "//trim(fname))
                    ldim = hed%getDims()
                    smpd = hed%getPixSz()
                    if( doprint )then
                        call hed%print_imghead
                        write(*,'(a,3(i0,1x))') 'Number of columns, rows, sections: ', ldim(1), ldim(2), ldim(3)
                        write(*,'(a,1x,f15.8)')  'Pixel size: ', smpd
                    endif
                case('F')
                    allocate(MrcImgHead :: hed)
                    call hed%new
                    call fopen(filnum, status='OLD', action='READ', file=fname, access='STREAM', iostat=ios)
                    call fileio_errmsg(" get_mrcfile_info fopen error "//trim(fname),ios)
                    call hed%read(filnum)
                    call fclose(filnum, errmsg=" get_mrcfile_info fclose error "//trim(fname))
                    if( doprint ) call hed%print_imghead
                case DEFAULT
                    write(*,*) 'The inputted file is not an MRC file; get_mrcfile_info; simple_jiffys'
                    write(*,*) fname
                    stop
            end select
        else
            write(*,*) 'The below file does not exists; get_mrcfile_info; simple_jiffys'
            write(*,*) fname
            stop
        endif
    end subroutine get_mrcfile_info

    !>  \brief is for gettign a part of the info in a SPIDER image header
    subroutine get_spifile_info( fname, ldim, iform, maxim, smpd, conv, doprint )
        character(len=*), intent(in)               :: fname
        integer, intent(out)                       :: ldim(3), maxim, iform
        real, intent(out)                          :: smpd
        character(len=:), allocatable, intent(out) :: conv
        logical, intent(in)                        :: doprint
        real    :: spihed(40)
        integer :: filnum, cnt, i, ios
        if( file_exists(fname) )then
            if( fname2format(fname) .eq. 'S' )then
                if( allocated(conv) ) deallocate(conv)
                call fopen(filnum, status='OLD', action='READ', file=fname, access='STREAM',iostat=ios)
                call fileio_errmsg(" get_spifile_info fopen error "//trim(fname),ios)
                call read_spihed
                call fclose(filnum,errmsg=" get_spifile_info fclose error "//trim(fname))
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='NATIVE')
                    call print_spihed
                    return
                endif
                call fopen(filnum, status='OLD', action='READ', file=fname, access='STREAM', iostat=ios)
                call fileio_errmsg(" get_spifile_info fopen error "//trim(fname),ios)
                call read_spihed
                call fclose(filnum,errmsg=" get_spifile_info fclose error "//trim(fname))
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='BIG_ENDIAN')
                    call print_spihed
                    return
                endif
                call fopen(filnum, status='OLD', action='READ', file=fname,&
                &access='STREAM', iostat=ios)
                call fileio_errmsg(" get_spifile_info fopen error "//trim(fname),ios)
                call read_spihed
                call fclose(filnum,errmsg=" get_spifile_info fclose error "//trim(fname))
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='LITTLE_ENDIAN')
                    call print_spihed
                    return
                endif
            else
                write(*,*) 'The inputted file is not a SPIDER file; get_spifile_info; simple_jiffys'
                write(*,*) fname
                stop
            endif
        else
            write(*,*) 'The below file does not exists; get_spifile_info; simple_jiffys'
            write(*,*) fname
            stop
        endif

        contains

            subroutine read_spihed
                cnt = 0
                do i=1,40*4,4
                    cnt = cnt+1
                    read(unit=filnum ,pos=i) spihed(cnt)
                end do
                ldim  = int([spihed(12), spihed(2), spihed(1)])
                iform = int(spihed(5))
                maxim = int(spihed(26))
                smpd  = spihed(38)
            end subroutine

            subroutine print_spihed
                if( doprint )then
                    write(*,'(a,3(i0,1x))') 'Number of columns, rows, sections: ', int(spihed(12)), int(spihed(2)), int(spihed(1))
                    write(*,'(a,1x,i3)')    'Iform descriptor: ', int(spihed(5))
                    write(*,'(a,1x,f7.0)')  'The number of the highest image currently used in the stack: ', spihed(26)
                    write(*,'(a,1x,f7.3)')  'Pixel size: ', spihed(38)
                endif
            end subroutine

    end subroutine get_spifile_info

    !>  \brief  is for finding logical dimension and number of particles in stack
    subroutine find_ldim_nptcls( fname, ldim, nptcls, doprint, formatchar )
        character(len=*),           intent(in)  :: fname      !< filename
        integer,                    intent(out) :: ldim(3)    !< logical dimension
        integer,                    intent(out) :: nptcls     !< number of particles
        logical,          optional, intent(in)  :: doprint    !< do print or not
        character(len=1), optional, intent(in)  :: formatchar !< input format
        integer                       :: iform, maxim
        real                          :: smpd
        character(len=:), allocatable :: conv
        character(len=1)              :: form
        logical                       :: ddoprint
        ddoprint = .false.
        if( present(doprint) ) ddoprint = doprint
        if( present(formatchar) )then
            form = formatchar
        else
            form = fname2format(fname)
        endif
        nptcls = 0
        select case (form)
            case('M','F')
                call get_mrcfile_info(fname, ldim, form, smpd, ddoprint )
                nptcls = ldim(3)
            case('S')
                call get_spifile_info(fname, ldim, iform, maxim, smpd, conv, ddoprint)
                nptcls = maxim
            case DEFAULT
                write(*,*) 'fname: ', fname
                write(*,*) 'format descriptor: ', fname2format(fname)
                stop 'File format not supported; find_ldim_nptcls; simple_procimgfile'
        end select
    end subroutine find_ldim_nptcls

    !>  \brief  is for checking logical dimension and number of particles in stack
    logical function has_ldim_nptcls( fname, ldim, nptcls )
        character(len=*), intent(in) :: fname   !< filename
        integer,          intent(in) :: ldim(3) !< expected logical dimension
        integer,          intent(in) :: nptcls  !< number of expected particles
        integer :: ldim_found(3), nptcls_found
        call find_ldim_nptcls( fname, ldim_found, nptcls_found )
        if( ldim_found(1) /= ldim(1) .or. ldim_found(2) /= ldim(2) )then
            has_ldim_nptcls = .false.
            return
        endif
        if( nptcls_found /= nptcls )then
            has_ldim_nptcls = .false.
            return
        endif
        has_ldim_nptcls = .true.
    end function has_ldim_nptcls
    

    ! polymorphic destructor

    subroutine kill( self )
        class(ImgHead), intent(inout) :: self
        if( self%exists )then
            select type(self)
                type is (MrcImgHead)
                    if( allocated(self%byte_array) )then
                        deallocate(self%byte_array, stat=alloc_stat)
                        if(alloc_stat /= 0) call alloc_errchk('In: simple_imghead; kill ', alloc_stat)
                    end if
                type is (SpiImgHead)
                    return
                class DEFAULT
                    stop 'Unsupported header type; kill; simple_imghead'
            end select
            self%exists = .false.
        endif
    end subroutine kill

    subroutine test_imghead
        use simple_fileio, only: fopen, fclose, fileio_errmsg
        class(ImgHead), allocatable :: hed, hed2
        integer :: recsz, funit, ios
        write(*,'(a)') '**info(simple_imghead_unit_test): testing read/write capabilities'
        allocate(SpiImgHead :: hed, hed2 )
        call hed%new([120,120,1])
        call hed2%new([120,120,1])
        recsz = 120*4
        call fopen(funit,file='test_imghed.spi',status='UNKNOWN',action='READWRITE',&
             access='STREAM',iostat=ios)
        call fileio_errmsg("test_imghead fopen error",ios)
        call hed%write(funit)
        call hed2%read(funit)
        call fclose(funit,ios,errmsg="test_imghead fclose error")
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
