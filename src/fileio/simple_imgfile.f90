! class to deal with image files on disks. Supported formats are:
!  - Spider: http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
!  - MRC: http://www2.mrc-lmb.cam.ac.uk/image2000.html
! This class is based on a class used in CTFFIND4, developed by Alexis Rohou
! and Nikolaus Grigorieff at Janelia Farm. The below copyright statement therefore
! needs to be included here:
! Copyright 2014 Howard Hughes Medical Institute
! All rights reserved
! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
! Modifications by Cyril Reboul, Michael Eager & Hans Elmlund
module simple_imgfile
include 'simple_lib.f08'
use gnufor2
implicit none

public :: imgfile
private
#include "simple_local_flags.inc"

integer,       parameter :: TIFF_SAMPLEFORMAT_UINT          = 1
integer,       parameter :: TIFF_SAMPLEFORMAT_INT           = 2
integer,       parameter :: TIFF_SAMPLEFORMAT_IEEEFP        = 3
integer,       parameter :: TIFF_SAMPLEFORMAT_VOID          = 4
integer,       parameter :: TIFF_SAMPLEFORMAT_COMPLEXINT    = 5
integer,       parameter :: TIFF_SAMPLEFORMAT_COMPLEXIEEEFP = 6

type imgfile
    private
    class(ImgHead), allocatable :: overall_head              !< Overall image head object
    character(len=STDLEN)       :: fname          = ''       !< Filename
    character(len=1)            :: head_format    = ''       !< 'M' (MRC), 'S' (SPIDER) file format
    integer                     :: funit          = 0        !< Unit number
    logical                     :: was_written_to = .false.  !< Indicates whether data was written to the file since it was opened
    logical                     :: isvol          = .false.  !< Indicates if SPIDER file is volume or stack
    logical                     :: existence      = .false.  !< Set to true when the object exists, false when it's closed
contains
    procedure          :: open
    procedure, private :: open_local
    procedure          :: close
    procedure, private :: slice2recpos
    procedure, private :: slice2bytepos
    procedure          :: rSlices
    procedure, private :: rTiffSlices
    procedure          :: wSlices
    procedure          :: wmrcSlices
    procedure          :: update_MRC_stats
    procedure          :: getDims
    procedure          :: getDim
    procedure          :: getIform
    procedure          :: getMode
    procedure          :: setIform
    procedure          :: setPixSz
    procedure          :: setRMSD
    procedure          :: setMean
    procedure          :: setMinMax
    procedure          :: setDims
    procedure          :: setMode
end type imgfile

! class variables (to avoid excessive allocation)
integer(kind=1), allocatable :: tmp_byte_array(:,:,:)
integer(kind=2), allocatable :: tmp_16bit_int_array(:,:,:)

contains

    subroutine open( self, fname, ldim, smpd, del_if_exists, formatchar, readhead, rwaction )
        class(imgfile),             intent(inout) :: self          !< Imagefile object to be created
        character(len=*),           intent(in)    :: fname         !< Filename
        integer,                    intent(in)    :: ldim(3)       !< logical dimension of image/stack
        real,                       intent(in)    :: smpd          !< Pixel size of image data (in Angstroms)
        logical,          optional, intent(in)    :: del_if_exists !< If the file already exists on disk, replace it
        character(len=1), optional, intent(in)    :: formatchar    !< file format character flag
        logical,          optional, intent(in)    :: readhead      !< header read indicator
        character(len=*), optional, intent(in)    :: rwaction      !< action flag
        character(len=1) :: format_descriptor
        logical          :: rreadhead, write_enabled
        rreadhead = .true.
        if( present(readhead) ) rreadhead = readhead
        write_enabled = .true.
        if( present(rwaction) )then
            if( rwaction .eq. 'READ' ) write_enabled = .false.
        endif
        ! set file name
        self%fname = trim(adjustl(fname))
        ! work out which file format to use
        if( present(formatchar) )then
            format_descriptor = formatchar
        else
            format_descriptor = fname2format(self%fname)
        endif
        if( format_descriptor .ne. 'N' )then
            self%head_format = format_descriptor
        else
            self%head_format  = DEFAULT_FILE_FORMAT
        endif
        ! Allocate head object
        select case(self%head_format)
            case ('M')
                allocate(MrcImgHead :: self%overall_head)
                call self%overall_head%new(ldim)
            case ('S')
                allocate(SpiImgHead :: self%overall_head)
                call self%overall_head%new(ldim)
                if( ldim(3) > 1 ) self%isvol = .true.
            case ('J','L')
                if( write_enabled ) THROW_HARD('TIFF format is only for reading; open')
                allocate(TiffImgHead :: self%overall_head)
                call self%overall_head%new(ldim)
            case DEFAULT
                THROW_HARD('unsupported file format')
        end select
        ! open the file
        call self%open_local(del_if_exists, rwaction)
        ! read/write the header
        if( file_exists(self%fname) )then
            ! read header
            if( rreadhead )then
                select case(self%head_format)
                    case ('M','S')
                        call self%overall_head%read(self%funit)
                    case('J','L')
                        call self%overall_head%read_tiff(self%fname)
                end select
            endif
        else
            ! write header
            call self%overall_head%write(self%funit)
        endif
        if( write_enabled )then
            call self%setPixSz(smpd)
        else
            call self%overall_head%setPixSz(smpd)
        endif
        ! The file-handle now exists
        self%existence = .true.
    end subroutine open

    !>  \brief open the file(s) for the imgfile
    subroutine open_local( self, del_if_exists, rwaction )
        class(imgfile),             intent(inout) :: self          !< Imagefile object
        logical, optional,          intent(in)    :: del_if_exists !< overwrite flag
        character(len=*), optional, intent(in)    :: rwaction      !< read/write flag
        character(len=9) :: rw_str
        character(len=7) :: stat_str
        integer          :: ios
        if( self%overall_head%isTiff() )then
            ! nothing to do, we are using the API
        else
            ! We need to prepare a string for the open statement
            if( present(rwaction) )then
                rw_str = trim(rwaction)
            else
                rw_str = 'READWRITE'
            endif
            ! What is the status of the file?
            stat_str = 'UNKNOWN'
            if( present(del_if_exists) )then
                if( del_if_exists )then
                    call del_file(self%fname)
                endif
            endif
            ! Get an IO unit number
            call fopen(self%funit,access='STREAM',file=self%fname,action=rw_str,status=stat_str,iostat=ios)
            call fileiochk("imgfile::open_local fopen error: "//trim(self%fname),ios)
        endif
        self%was_written_to = .false.
    end subroutine open_local

    !>  \brief  close the file(s) and "de-initialise" the imgfile object
    subroutine close( self )
        class(imgfile), target, intent(inout) :: self
        class(ImgHead), pointer :: ptr=>null()
        integer :: ios
        ptr => self%overall_head
        select type(ptr)
            type is (TiffImgHead)
                call self%overall_head%CloseTiff
            class DEFAULT
                if( is_open(self%funit) )then
                    if( self%was_written_to )then
                        call self%overall_head%write(self%funit)
                    endif
                    call fclose(self%funit)
                endif
        end select
        if( allocated(self%overall_head) )then
            call self%overall_head%kill
            deallocate(self%overall_head)
        endif
        self%was_written_to = .false.
        self%existence      = .false.
    end subroutine close

    !>  \brief  for translating an image index to record indices in the stack
    !! \param[out] hedinds,iminds header and image indices in the stack
    subroutine slice2recpos( self, nr, hedinds, iminds )
        class(imgfile), target, intent(in)  :: self
        integer,                intent(in)  :: nr
        integer(kind=8),        intent(out) :: hedinds(2), iminds(2)
        integer :: cnt, j, dims(3)
        class(ImgHead), pointer :: ptr=>null()
        ptr => self%overall_head
        select type(ptr)
        type is (MrcImgHead)
            THROW_HARD('cannot translate an image index to record indices for MRC files')
        type is (SpiImgHead)
            cnt  = self%overall_head%getLabrec()
            dims = self%overall_head%getDims()
            if( nr == 0 )then
                hedinds(1) = 1
                hedinds(2) = self%overall_head%getLabrec()
                iminds     = 0 ! no image associated with the stack header
            else
                do j=1,nr
                    hedinds(1) = cnt + 1 ! hed from
                    cnt        = cnt + self%overall_head%getLabrec()
                    hedinds(2) = cnt     ! hed to
                    iminds(1)  = cnt + 1 ! im from
                    cnt        = cnt + dims(2)
                    iminds(2)  = cnt     ! im to
                end do
            endif
        type is (TiffImgHead)
            THROW_HARD('Should not be here for TIFF format; slice2recpos')
        class DEFAULT
            THROW_HARD('format not supported')
        end select
    end subroutine slice2recpos

    !>  \brief  for translating an image index to record indices in the stack
    !! \param[out] hedinds,iminds  indices in the stack
    subroutine slice2bytepos( self, nr, hedinds, iminds )
        class(imgfile), target, intent(in)    :: self
        integer,                intent(in)    :: nr
        integer(kind=8),        intent(inout) :: hedinds(2), iminds(2)
        class(ImgHead), pointer :: ptr=>null()
        ptr => self%overall_head
        select type(ptr)
        type is (TiffImgHead)
            THROW_HARD('Should not be here for TIFF format; slice2recpos')
        end select
        if( nr < 0 )then
            THROW_HARD('cannot have negative slice indices')
        else if( nr == 0 )then
            hedinds(1) = 1
            hedinds(2) = self%overall_head%getLabbyt()
            iminds     = 0 ! no image in SPIDER stack header
        else
            call self%slice2recpos( nr, hedinds, iminds )
            hedinds(1) = (hedinds(1) - 1) * self%overall_head%getLenbyt() + 1
            hedinds(2) = hedinds(2) * self%overall_head%getLenbyt()
            iminds(1)  = (iminds(1) - 1) * self%overall_head%getLenbyt() + 1
            iminds(2)  = iminds(2) * self%overall_head%getLenbyt()
        endif
    end subroutine slice2bytepos

    !>  \brief  reads a set of contiguous slices of the image file from disk into memory.
    !!          The array of reals should have +2 elements in the first dimension.
    subroutine rSlices( self, first_slice, last_slice, rarr, is_mrc )
        use, intrinsic :: iso_c_binding
        class(imgfile), target, intent(inout) :: self         !< instance  Imagefile object
        integer,                intent(in)    :: first_slice  !< First slice (the first slice in the file is numbered 1)
        integer,                intent(in)    :: last_slice   !< Last slice
        real,                   intent(inout) :: rarr(:,:,:)  !< Array of reals. Will be (re)allocated if needed
        logical,                intent(in)    :: is_mrc
        character(len=100)          :: io_message
        integer                     :: io_stat,dims(3),tmparrdims(3)
        integer(kind=8)             :: first_byte,hedbyteinds(2),imbyteinds(2),first_hedbyte,byteperpix
        logical                     :: arr_is_ready,alloc_tmparr
        class(ImgHead), pointer     :: ptr=>null()
        byteperpix = int(self%overall_head%bytesPerPix(),kind=8)
        dims       = self%overall_head%getDims()
        dims(3)    = last_slice - first_slice + 1
        if( byteperpix == 4 .and. is_mrc )then
            ! fast read
            first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                &*int(product(dims(1:2)),kind=8)*byteperpix
            read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) rarr(:dims(1),:,:)
        else
            ! slower, polymorphic read with more checkpoints
            ! set pointer to overall header
            ptr => self%overall_head
            ! Tiff format
            select type( ptr )
                type is( TiffImgHead )
                    call self%rTiffSlices(first_slice, rarr)
                    return
            end select
            ! Work out the position of the first byte
            select case(self%head_format)
                case('M','F')
                    first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                    &*int(product(dims(1:2)),kind=8)*byteperpix
                case('S')
                    if( self%isvol )then
                        first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                        &*int(product(dims(1:2)),kind=8)*byteperpix
                    else
                        call self%slice2bytepos(first_slice, hedbyteinds, imbyteinds)
                        first_byte    = imbyteinds(1)  ! first image byte
                        first_hedbyte = hedbyteinds(1) ! first header byte
                    endif
                case DEFAULT
                    THROW_HARD('format not supported')
            end select
            rarr = 0. ! initialize to zero
            select case(byteperpix)
                case(1) ! Byte data
                    if( allocated(tmp_byte_array) )then
                        tmparrdims(1) = size(tmp_byte_array,1)
                        tmparrdims(2) = size(tmp_byte_array,2)
                        tmparrdims(3) = size(tmp_byte_array,3)
                        alloc_tmparr  = .false.
                        if( any(tmparrdims .ne. dims) )then
                            deallocate(tmp_byte_array)
                            alloc_tmparr = .true.
                        endif
                    else
                        alloc_tmparr = .true.
                    endif
                    if( alloc_tmparr )then
                        allocate(tmp_byte_array(dims(1),dims(2),dims(3)))
                    endif
                    read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_byte_array(:dims(1),:dims(2),:dims(3))
                    ! Conversion from unsigned byte integer (which MRC appears to be) is tricky because Fortran
                    ! doesn't do unsigned integer natively. The following IAND trick is courtesy of Jim Dempsey
                    ! at http://software.intel.com/en-us/forums/showthread.php?t=64400 Confusingly, the MRC format
                    ! documentation implies that one should expect signed integers, which seems to be incorrect:
                    ! http://www2.mrc-lmb.cam.ac.uk/image2000.html IMOD documentation indicates that prior to IMOD
                    ! 4.2.23, unsigned bytes were used and that one needs to inspect the imodStamp head to check
                    if( self%overall_head%pixIsSigned() )then
                        rarr(1:dims(1),:,:) = tmp_byte_array(:dims(1),:dims(2),:dims(3))
                    else
                        rarr(1:dims(1),:,:) = real(iand(int(tmp_byte_array(:dims(1),:dims(2),:dims(3)),kind=4),int(255,kind=4)))
                    endif
                case(2) ! 16-bit data, assumed to be integer
                    if( allocated(tmp_16bit_int_array) )then
                        tmparrdims(1) = size(tmp_16bit_int_array,1)
                        tmparrdims(2) = size(tmp_16bit_int_array,2)
                        tmparrdims(3) = size(tmp_16bit_int_array,3)
                        alloc_tmparr  = .false.
                        if( any(tmparrdims .ne. dims) )then
                            deallocate(tmp_16bit_int_array)
                            alloc_tmparr = .true.
                        endif
                    else
                        alloc_tmparr = .true.
                    endif
                    if( alloc_tmparr )then
                        allocate(tmp_16bit_int_array(dims(1),dims(2),dims(3)))
                    endif
                    read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_16bit_int_array(:dims(1),:dims(2),:dims(3))
                    if( self%overall_head%pixIsSigned() )then
                        rarr(1:dims(1),:,:) = real(tmp_16bit_int_array(:dims(1),:dims(2),:dims(3)))
                    else
                        rarr(1:dims(1),:,:) = real(iand(int(tmp_16bit_int_array(:dims(1),:dims(2),:dims(3)),kind=4),&
                            &int(huge(int(1,kind=2)), kind=4)))
                    endif
                case(4)
                    read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) rarr(:dims(1),:,:)
                case DEFAULT
                    write(logfhandle,'(2a)') 'fname: ', trim(self%fname)
                    write(logfhandle,'(a,i0,a)') 'bit depth: ', self%overall_head%bytesPerPix(), ' bytes'
                    THROW_HARD('unsupported bit-depth')
            end select
        endif
        ! Check the read was successful
        if( io_stat .ne. 0 )then
            write(logfhandle,'(a,i0,2a)') '**ERROR(rSlices): I/O error ', io_stat, ' when reading from: ', trim(self%fname)
            write(logfhandle,'(2a)') 'IO error message was: ', trim(io_message)
            THROW_HARD('I/O')
        endif
    end subroutine rSlices

    !>  \brief  reads a set of TIFF contiguous slices of the image file from disk into memory.
    !!          The array of reals should have +2 elements in the first dimension.
    !!          2-Dimensional images are assumed
    subroutine rTiffSlices( self, image_index, rarr )
        use simple_strings, only: int2str, tocstring
        use, intrinsic :: iso_c_binding
        class(imgfile), target, intent(inout) :: self         !< instance  Imagefile object
        integer,                intent(in)    :: image_index  !<
        real,                   intent(inout) :: rarr(:,:,:)  !< Array of reals. Will be (re)allocated if needed
        class(ImgHead),  pointer :: ptr=>null()
        real(kind=4),    pointer :: ftemp_byte_array32(:)
        integer(kind=1), pointer :: itemp_byte_array8(:)
        integer(kind=2), pointer :: itemp_byte_array16(:)
        type(c_ptr)              :: strip_buffer_cptr, string_ptr
        real(kind=dp) :: v
        integer :: dims(3), nbytes_in_buffer, first_byte, first_row
        integer :: io_stat, last_byte, last_row, row_counter, strip_counter
        integer :: last_strip_length, strip_length, nbytes, nrows
        character,pointer,dimension(:) :: string
        ! checks have been performed in rSlices
        select case(self%head_format)
        case('J','L')
            ptr  => self%overall_head
        case DEFAULT
            THROW_HARD('format not supported; rTiffSlices')
        end select
        ! Get the dims of the image file
        dims = self%overall_head%getDims()
        rarr = 0. ! initialize to zero
        select type(ptr)
        type is( TiffImgHead )
#ifdef USING_TIFF
        if( self%head_format == 'L' ) call TIFFMuteWarnings
        io_stat = TIFFSetDirectory(ptr%fhandle,int(image_index-1,kind=c_int16_t))
        if( self%head_format == 'L' ) call TIFFUnMuteWarnings
        if (io_stat .ne. 1) THROW_HARD('rTiffSlices; Failed to set directory to: '//int2str(image_index-1))
        ! We don't support tiled organization, only strips
        if (ptr%IsTiled) THROW_HARD('rTiffSlices; Tile-based TIFF files not supported')
        ! Number of Bits/Sample Format
        nbytes = 1
        select case(ptr%BitsPerSample)
            case(8)
                ! 8-bit: unsigned
                if(ptr%SampleFormat .ne. TIFF_SAMPLEFORMAT_UINT)then
                    THROW_HARD('rTiffSlices; 8-bit Unsupported sample format: '//int2str(ptr%SampleFormat))
                endif
            case(16)
                ! 16-bit: signed/unsigned
                if((ptr%SampleFormat.ne.TIFF_SAMPLEFORMAT_UINT) .and. (ptr%SampleFormat.ne.TIFF_SAMPLEFORMAT_INT) )then
                    THROW_HARD('rTiffSlices; 8/16-bit Unsupported sample format: '//int2str(ptr%SampleFormat))
                endif
                nbytes = 2
            case(32)
                ! 32-bit: float
                if(ptr%SampleFormat .ne. TIFF_SAMPLEFORMAT_IEEEFP)then
                    THROW_HARD('rTiffSlices; Unsupported 32-bit sample format: '//int2str(ptr%SampleFormat))
                endif
                nbytes = 4
            case DEFAULT
                THROW_HARD('rTiffSlices; Unsupported combination of # bits & bit depth: '//int2str(ptr%BitsPerSample)//' & '//int2str(ptr%SampleFormat))
        end select
        if (ptr%SamplesPerPixel .ne. 1) then
            THROW_HARD('rTiffSlices; Unsupported number of samples per pixel: '//int2str(ptr%SamplesPerPixel))
        endif
        strip_length      = nbytes*ptr%RowsPerStrip*dims(1)
        strip_buffer_cptr = ptr%StripBuffer_cptr
        ! Allocating buffer prior to read to circumvent overflow in versions < 4.0.9
        strip_buffer_cptr = TIFFmalloc(strip_length)
        do strip_counter  = 1,ptr%NumberOfStrips
            ! Read a strip from the file
            if( (strip_counter==ptr%NumberOfStrips) .and. (ptr%NumberOfStrips>1) )then
                ! making sure buffer is allocated to correct size for the very last strip
                v = real(ptr%StripSize,dp) / real(ptr%RowsPerStrip,dp) ! order is to avoid integer kind=4 overflow!!
                v = v * real(dims(2),dp)
                last_strip_length = nint(v) - ((strip_counter-1)*ptr%StripSize+1) + 1
                if( last_strip_length /= ptr%StripSize )then
                    io_stat           = TIFFfree(strip_buffer_cptr)
                    strip_length      = nbytes*last_strip_length
                    strip_buffer_cptr = TIFFmalloc(strip_length)
                endif
            endif
            nrows = strip_length / dims(1) / nbytes ! rows to be read
            nbytes_in_buffer = TIFFReadEncodedStrip(ptr%fhandle, int(strip_counter-1,kind=c_int32_t),&
                              &strip_buffer_cptr, int(strip_length,kind=c_int32_t))
            if (nbytes_in_buffer .le. 0) THROW_HARD('rTiffSlices; Failed to read from strip: '//int2str(strip_counter))
            if (nbytes_in_buffer .ne. strip_length)then
                string_ptr = TIFFGetVersion()
                call c_f_pointer(string_ptr, string, [23])
                write(logfhandle,*)string(:)
                write(logfhandle,*)'strip_length     ',strip_length
                write(logfhandle,*)'nbytes_in_buffer ',nbytes_in_buffer
                write(logfhandle,*)'StripSize        ',ptr%Stripsize
                write(logfhandle,*)'RowsPerStrip     ',ptr%RowsPerStrip
                write(logfhandle,*)'NumberOfStrips   ',ptr%NumberOfStrips
                write(logfhandle,*)'SampleFormat     ',ptr%SampleFormat
                write(logfhandle,*)'BitsPerSample    ',ptr%BitsPerSample
                THROW_HARD('rTiffSlices; Unexpected number of bytes in buffer 1')
            endif
            select case(ptr%BitsPerSample)
            case(8)
                ! Setup F pointer to the strip data
                call c_f_pointer(strip_buffer_cptr,itemp_byte_array8,[nbytes_in_buffer])
                ! update row counter
                first_row = ptr%RowsPerStrip * (strip_counter-1) + 1
                last_row  = first_row + nrows - 1
                ! cast and copy
                do row_counter = first_row,last_row
                    first_byte = (row_counter-first_row)*dims(1)+1
                    last_byte  = first_byte + dims(1) - 1
                    rarr(1:dims(1),dims(2)-row_counter+1,1) = real(iand(int(itemp_byte_array8(first_byte:last_byte),kind=4),&
                        &int(huge(int(1,kind=1)), kind=4)*2+1))
                enddo
            case(16)
                call c_f_pointer(strip_buffer_cptr,itemp_byte_array16,[nbytes_in_buffer])
                first_row = ptr%RowsPerStrip * (strip_counter-1) + 1
                last_row  = first_row + strip_length/dims(1)/nbytes - 1
                do row_counter = first_row,last_row
                    first_byte = (row_counter-first_row)*dims(1)+1
                    last_byte  = first_byte + dims(1) - 1
                    rarr(1:dims(1),dims(2)-row_counter+1,1) = real(iand(int(itemp_byte_array16(first_byte:last_byte),kind=4),&
                    &int(huge(int(1,kind=2)), kind=4)*2+1))
                enddo
            case(32)
                call c_f_pointer(strip_buffer_cptr,ftemp_byte_array32,[nbytes_in_buffer])
                first_row = ptr%RowsPerStrip * (strip_counter-1) + 1
                last_row  = first_row + nrows - 1
                do row_counter = first_row,last_row
                    first_byte = (row_counter-first_row)*dims(1)+1
                    last_byte  = first_byte + dims(1) - 1
                    rarr(1:dims(1),dims(2)-row_counter+1,1) = real(ftemp_byte_array32(first_byte:last_byte),kind=4)
                enddo
            end select
        enddo
        ! cleanup
        io_stat = TIFFfree(strip_buffer_cptr)
#endif
        end select
    end subroutine rTiffSlices

    !>  \brief  read/write a set of contiguous slices of the image file from disk into memory.
    !!          The array of reals should have +2 elements in the first dimension.
    subroutine wSlices( self, first_slice, last_slice, rarr, ldim, is_ft, smpd )
        use, intrinsic :: iso_c_binding
        class(imgfile), target, intent(inout) :: self         !< instance  Imagefile object
        integer,                intent(in)    :: first_slice  !< First slice (the first slice in the file is numbered 1)
        integer,                intent(in)    :: last_slice   !< Last slice
        real,                   intent(inout) :: rarr(:,:,:)  !< Array of reals. Will be (re)allocated if needed
        integer,                intent(in)    :: ldim(3)      !< Logical size of the array. This will be written to disk: rarr(1:ldim(1),:,:)
        logical,                intent(in)    :: is_ft        !< to indicate FT status of image
        real,                   intent(in)    :: smpd         !< sampling distance
        integer                     :: io_stat,itmp,dims(3)!,dims_stored(3)
        integer(kind=8)             :: first_byte,hedbyteinds(2),imbyteinds(2),first_hedbyte,byteperpix
        logical                     :: arr_is_ready
        real                        :: min_val,max_val
        class(ImgHead), pointer     :: ptr=>null()
        class(ImgHead), allocatable :: imghed
        ! Check that the first and last slice numbers given make sense
        if( first_slice > 0 .and. (first_slice .gt. last_slice) ) THROW_HARD('last < first slice')
        ! Get the dims of the image file
        dims = self%overall_head%getDims()
        ! set pointer to overall header
        ptr  => self%overall_head
        select type(ptr)
        type is (MrcImgHead)
            ! Redefine the file dims (to keep track of the index of the last image of the stack)
            ! dims_stored = dims
            dims(1)     = ldim(1)
            dims(2)     = size(rarr,2)
            dims(3)     = max(last_slice,dims(3))
            call self%overall_head%setDims(dims)
            ! dims = dims_stored ! safety
        type is (TiffImgHead)
            THROW_HARD('Writing is unsupported in TIFF format')
        end select
        if( is_odd(dims(1)) )then
            arr_is_ready = (size(rarr,1) .eq. dims(1) + 1) .and. size(rarr,2) .eq. dims(2)
        else
            arr_is_ready = (size(rarr,1) .eq. dims(1) + 2) .and. size(rarr,2) .eq. dims(2)
        endif
        if( .not. arr_is_ready )then
            write(logfhandle,*) 'Array size: ', size(rarr,1), size(rarr,2), size(rarr,3)
            write(logfhandle,*) 'Dimensions: ', dims(1), dims(2), dims(3)
            THROW_HARD('array not properly allocated')
        endif
        byteperpix = int(self%overall_head%bytesPerPix(),kind=8)
        ! Work out the position of the first byte
        select case(self%head_format)
        case('M','F')
            first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                &*int(product(dims(1:2)),kind=8)*byteperpix
        case('S')
            if( self%isvol )then
                first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                    &*int(product(dims(1:2)),kind=8)*byteperpix
            else
                call self%slice2bytepos(first_slice, hedbyteinds, imbyteinds)
                first_byte    = imbyteinds(1)  ! first image byte
                first_hedbyte = hedbyteinds(1) ! first header byte
            endif
        case DEFAULT
            THROW_HARD('format not supported')
        end select
        ! find minmax
        max_val = maxval(rarr)
        min_val = minval(rarr)
        select case(self%head_format)
        case('M','F')
            max_val = max(max_val, self%overall_head%getMaxPixVal())
            min_val = min(min_val, self%overall_head%getMinPixVal())
            call self%overall_head%setMinPixVal(min_val)
            call self%overall_head%setMaxPixVal(max_val)
        case('S')
            if( .not. self%isvol )then
                ! for SPIDER stack we also need to create and write an image header
                allocate( SpiImgHead :: imghed )
                call imghed%new(ldim=ldim)
                call imghed%setMinimal(ldim, is_ft, smpd)
                call imghed%setMinPixVal(min_val)
                call imghed%setMaxPixVal(max_val)
                call imghed%write(self%funit, pos=first_hedbyte)
                call imghed%kill
                deallocate(imghed)
            endif
        end select
        write(unit=self%funit,pos=first_byte,iostat=io_stat) rarr(1:dims(1),:,:)
        ! Check the write was successful
        if( io_stat .ne. 0 )then
            write(logfhandle,'(a,i0,2a)') '**ERROR(wSlices): I/O error ', io_stat, ' when writing to: ', trim(self%fname)
            THROW_HARD('I/O')
        endif
        ! May need to update file dims
        select case(self%head_format)
        case('M')
            dims(3) = max(dims(3),last_slice)
            call self%overall_head%setDims(dims)
        case('S')
            if( .not. self%isvol )then
                itmp = self%overall_head%getMaxim()
                call self%overall_head%setMaxim(max(itmp,last_slice))
            endif
        end select
        ! May need to update FT status in the overall header head
        select type(ptr)
        type is (MrcImgHead)
            if( is_ft ) call self%overall_head%setMode(4)
        type is (SpiImgHead)
            dims(1) = size(rarr,1)
            dims(2) = size(rarr,2)
            dims(3) = size(rarr,3)
            if( .not. self%isvol .and. .not. is_ft )then
                call self%overall_head%setIform(1)
            else if( self%isvol .and. .not. is_ft )then
                call self%overall_head%setIform(3)
            else if( .not. self%isvol .and. is_ft .and. .not. is_even(dims(1:2)) )then
                call self%overall_head%setIform(-11)
            else if( .not. self%isvol .and. is_ft .and. is_even(dims(1:2)) )then
                call self%overall_head%setIform(-12)
            else if( self%isvol .and. is_ft .and. .not. is_even(dims(1:2)) )then
                call self%overall_head%setIform(-21)
            else if( self%isvol .and. is_ft .and. is_even(dims(1:2)) )then
                call self%overall_head%setIform(-22)
            else
                THROW_HARD('undefined file type')
            endif
        end select
        ! Remember that we wrote to the file
        self%was_written_to = .true.
    end subroutine wSlices

    !>  \brief  read/write a set of contiguous slices of the image file from disk into memory.
    !!          The array of reals should have +2 elements in the first dimension.
    subroutine wmrcSlices( self, first_slice, last_slice, rarr, ldim, is_ft )
        use, intrinsic :: iso_c_binding
        class(imgfile), target, intent(inout) :: self         !< instance  Imagefile object
        integer,                intent(in)    :: first_slice  !< First slice (the first slice in the file is numbered 1)
        integer,                intent(in)    :: last_slice   !< Last slice
        real,                   intent(inout) :: rarr(:,:,:)  !< Array of reals. Will be (re)allocated if needed
        integer,                intent(in)    :: ldim(3)      !< Logical size of the array
        logical,                intent(in)    :: is_ft        !< to indicate FT status of image
        integer         :: io_stat,itmp,dims(3)
        integer(kind=8) :: first_byte,first_hedbyte,byteperpix
        real            :: min_val,max_val
        ! Redefine the file dims (to keep track of the index of the last image of the stack)
        dims = self%overall_head%getDims()
        dims(1)     = ldim(1)
        dims(2)     = size(rarr,2)
        dims(3)     = max(last_slice,dims(3))
        call self%overall_head%setDims(dims)
        ! dims = dims_stored ! safety
        byteperpix = int(self%overall_head%bytesPerPix(),kind=8)
        ! Work out the position of the first byte
        first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
            &*int(product(dims(1:2)),kind=8)*byteperpix
        ! find minmax
        max_val = maxval(rarr)
        min_val = minval(rarr)
        max_val = max(max_val, self%overall_head%getMaxPixVal())
        min_val = min(min_val, self%overall_head%getMinPixVal())
        call self%overall_head%setMinPixVal(min_val)
        call self%overall_head%setMaxPixVal(max_val)
        write(unit=self%funit,pos=first_byte,iostat=io_stat) rarr(1:dims(1),:,:)
        ! Check the write was successful
        if( io_stat .ne. 0 )then
            write(logfhandle,'(a,i0,2a)') '**ERROR(wSlices): I/O error ', io_stat, ' when writing to: ', trim(self%fname)
            THROW_HARD('I/O')
        endif
        ! May need to update file dims
        dims(3) = max(dims(3),last_slice)
        call self%overall_head%setDims(dims)
        ! May need to update FT status in the overall header head
        if( is_ft ) call self%overall_head%setMode(4)
        ! Remember that we wrote to the file
        self%was_written_to = .true.
    end subroutine wmrcSlices

    !>  \brief  read/write a set of contiguous slices of the image file from disk into memory.
    !!          The array of reals should have +2 elements in the first dimension.
    subroutine update_MRC_stats( self, stats )
        use, intrinsic :: iso_c_binding
        class(imgfile), target, intent(inout) :: self         !< instance  Imagefile object
        real,                   intent(in)    :: stats(4)
        select case(self%head_format)
        case('M','F')
            call self%setMinmax(stats(1), stats(2))
            call self%setMean(stats(3))
            call self%setRMSD(stats(4))
        case('S','J','L')
            ! this routine is for MRC only
            return
        case DEFAULT
            THROW_HARD('undefined file type')
        end select
        ! Remember that we wrote to the file
        self%was_written_to = .true.
    end subroutine update_MRC_stats

    !>  \brief  Return the dimension of the image stack
    function getDims( self )
        class(imgfile), intent(in) :: self   !< Imagefile object
        integer :: getDims(3)
        getDims = self%overall_head%getDims()
    end function getDims

    !>  \brief  Return one of the dimensio of the image stack
    function getDim( self, which_dim )
        class(imgfile), intent(in) :: self
        integer, intent(in)        :: which_dim
        integer :: getDim
        getDim = self%overall_head%getDim(which_dim)
    end function getDim

    !>  \brief  Return the format descriptor of the stack
    integer function getIform( self )
        class(imgfile), intent(in) :: self   !< Imagefile object
        getIform = self%overall_head%getIform()
    end function getIform

    !>  \brief  Return the format descriptor of the stack
    integer function getMode( self )
        class(imgfile), intent(in) :: self   !< Imagefile object
        getMode = self%overall_head%getMode()
    end function getMode

    !>  \brief  Set the format descriptor of the stack
    subroutine setIform( self, iform )
        class(imgfile), intent(inout) :: self   !< Imagefile object
        integer, intent(in)           :: iform   !< format descriptor
        call self%overall_head%setIform(iform)
        self%was_written_to = .true.
    end subroutine setIform

    !>  \brief  Set the pixel size of the stack
    subroutine setPixSz( self, smpd )
        class(imgfile), intent(inout) :: self !< Imagefile object
        real, intent(in)              :: smpd !< sample resolution or pixel size (angstroms)
        call self%overall_head%setPixSz(smpd)
        self%was_written_to = .true.
    end subroutine setPixSz

    !>  \brief  Set the RMS of the stack
    subroutine setRMSD( self, dev )
        class(imgfile), intent(inout) :: self !< Imagefile object
        real,           intent(in)    :: dev !< rmsd
        call self%overall_head%setRMSD(dev)
        self%was_written_to = .true.
    end subroutine setRMSD

    !>  \brief  Set the average value of the stack
    subroutine setMean( self, mean )
        class(imgfile), intent(inout) :: self !< Imagefile object
        real,           intent(in)    :: mean !< mean
        call self%overall_head%setMean(mean)
        self%was_written_to = .true.
    end subroutine setMean

    !>  \brief  Set the minimum & maximum values of the stack
    subroutine setMinMax( self, minv, maxv )
        class(imgfile), intent(inout) :: self !< Imagefile object
        real,           intent(in)    :: minv, maxv
        call self%overall_head%setMinPixVal(minv)
        call self%overall_head%setMaxPixVal(maxv)
        self%was_written_to = .true.
    end subroutine setMinmax

    !>  \brief  Set the mode of the MRC file
    subroutine setMode( self, mode )
        class(imgfile), intent(inout) :: self !< Imagefile object
        integer, intent(in)           :: mode !< MRC mode type
        call self%overall_head%setMode(mode)
        self%was_written_to = .true.
    end subroutine setMode

    !>  \brief  for setting the logical dimensions
    subroutine setDims( self, ldim )
        class(imgfile), intent(inout) :: self !< Imagefile object
        integer, intent(in)           :: ldim(3) !< dimensions
        call self%overall_head%setDims(ldim)
        self%was_written_to = .true.
    end subroutine setDims

end module simple_imgfile
