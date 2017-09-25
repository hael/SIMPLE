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
use simple_defs
use simple_syslib,  only: alloc_errchk, file_exists, is_open
use simple_fileio,  only: fileio_errmsg, fopen, fclose, file_size, fname2format, del_file
use simple_imghead  ! only public ImgHead, MrcImgHead, SpiImgHead, test_imghead, find_ldim_nptcls
use gnufor2
implicit none

public :: imgfile
private
#include "simple_local_flags.inc"

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
    ! CORE FUNCTIONALITY
    procedure          :: open
    procedure          :: print_header
    procedure, private :: open_local
    procedure          :: close
    procedure, private :: close_nowrite
    procedure, private :: exists
    procedure, private :: get_format
    procedure, private :: slice2recpos
    procedure, private :: slice2bytepos
    procedure          :: print_slice2bytepos
    ! SLICE & HEADER I/O
    procedure          :: rSlice
    procedure          :: wSlice
    procedure          :: rHead
    procedure          :: wHead
    procedure          :: rwSlices
    ! GETTERS, SETTERS, PRINTERS
    procedure          :: print_imgfile
    procedure          :: getDims
    procedure          :: getDim
    procedure          :: getFname
    procedure          :: getPixSz
    procedure          :: getIform
    procedure          :: getMode
    procedure          :: setIform
    procedure          :: setPixSz
    procedure          :: setRMSD
    procedure          :: setMode
end type imgfile

contains

    ! CONSTRUCTOR
    !>  \brief  constructs an imgfile object (file-handle)
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
        call self%close_nowrite
        ! Remove leading blanks in fname
        self%fname = adjustl(fname)
        DebugPrint   'trying to open: ', self%fname
        ! Work out which file format to use
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
        DebugPrint   'format: ', self%head_format
        ! Allocate head object
        select case(self%head_format)
        case ('M')
            allocate(MrcImgHead :: self%overall_head)
            DebugPrint   ' allocated MrcImgHead'
            call self%overall_head%new(ldim)
        case ('S')
            ! create header
            allocate(SpiImgHead :: self%overall_head)
            DebugPrint   ' allocated SpiImgHead'
            call self%overall_head%new(ldim)
            if( ldim(3) > 1 ) self%isvol = .true.
        case DEFAULT
            stop 'Unsupported file format; new; simple_imgfile'
        end select
        ! open the file
        call self%open_local(del_if_exists, rwaction)
        ! read/write the header
        if( file_exists(self%fname) )then
            ! read header
            if( rreadhead )then
                call self%overall_head%read(self%funit)
                DebugPrint   'did read header'
            endif
        else
            ! write header
            call self%overall_head%write(self%funit)
            DebugPrint  'did write header'
        endif
        if( write_enabled )then
            ! REPLACED WITH THIS ONE TO MAKE THE FLAG WAS_WRITTEN TO TRUE
            ! NEEDED FOR GETTING THE HEADER INFORMATION CORRECT
            call self%setPixSz(smpd)
        else
            call self%overall_head%setPixSz(smpd)
        endif
        if( debug ) call self%overall_head%print_imghead()
        ! The file-handle now exists
        self%existence = .true.
        DebugPrint  '(imgfile::new) constructed an imgfile object (file-handle)'
    end subroutine open

    !>  \brief is forprinting the header
    subroutine print_header( self )
        class(imgfile), intent(in) :: self   !< Imagefile object 
        call self%overall_head%print_imghead()
    end subroutine print_header

    ! CORE FUNCTIONALITY

    !>  \brief open the file(s) for the imgfile
    subroutine open_local( self, del_if_exists, rwaction )
        class(imgfile),             intent(inout) :: self   !< Imagefile object 
        logical, optional,          intent(in) :: del_if_exists !< overwrite flag
        character(len=*), optional, intent(in) :: rwaction      !< read/write flag
        character(len=9) :: rw_str
        character(len=7) :: stat_str
        integer :: ios
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
        ! Get an IO unit number for the head file
        call fopen(self%funit,access='STREAM',file=self%fname,action=rw_str,&
             status=stat_str,iostat=ios)
        call fileio_errmsg("imgfile::open_local fopen error",ios)
        self%was_written_to = .false.
    end subroutine open_local

    !>  \brief  close the file(s) and "de-initialise" the imgfile object
    subroutine close( self )
        class(imgfile), intent(inout) :: self
        integer :: ios
        if( is_open(self%funit) )then
            if( self%was_written_to )then
                call self%overall_head%write(self%funit)
                DebugPrint  '(simple_imgfile::close) wrote overall_head'
            endif
            call fclose(self%funit, ios,errmsg="simple_imgfile::close error")
        endif
        if( allocated(self%overall_head) ) call self%overall_head%kill
        if( allocated(self%overall_head) ) deallocate(self%overall_head)
        self%was_written_to = .false.
        self%existence = .false.
    end subroutine close

    !>  \brief  close the file(s)
    subroutine close_nowrite( self )
        class(imgfile), intent(inout) :: self   !< Imagefile object
        integer :: ios
        call fclose( self%funit, ios,errmsg="simple_imgfile close nowrite error")
        if( allocated(self%overall_head) ) call self%overall_head%kill
        if( allocated(self%overall_head) ) deallocate(self%overall_head)
        self%was_written_to = .false.
        self%existence = .false.
    end subroutine close_nowrite

    !>  \brief  Check whether the file exists on disk
    logical function exists( self )
        class(imgfile), intent(inout) :: self   !< Imagefile object 
        exists = file_exists(self%fname)
        if( exists) exists = exists .and. file_size(self%fname) .gt. 0
    end function exists

    !>  \brief Return a one letter code for the file format designated by the extension in the fname
    !!         if .mrc: M
    !!         if .spi: S
    !!         if .img: I
    !!         if .hed: I
    !!         else: N
    pure function get_format( self )
        class(imgfile), intent(in) :: self   !< Imagefile object 
        character(len=1)           :: get_format
        get_format = self%head_format
    end function get_format

    !>  \brief  for translating an image index to record indices in the stack
    !! \param[out] hedinds,iminds header and image indices in the stack
    subroutine slice2recpos( self, nr, hedinds, iminds )
        class(imgfile), intent(in), target :: self   !< Imagefile object 
        integer, intent(in)                :: nr      !< num images
        integer(kind=8), intent(out)       :: hedinds(2), iminds(2)
        integer                            :: cnt, j, dims(3)
        class(ImgHead), pointer            :: ptr=>null()
        ptr => self%overall_head
        select type(ptr)
        type is (MrcImgHead)
            stop 'Cannot translate an image index to record indices for MRC files; slice2recpos; simple_imgfile'
        type is (SpiImgHead)
            cnt = self%overall_head%getLabrec()
            dims = self%overall_head%getDims()
            if( nr == 0 )then
                hedinds(1) = 1
                hedinds(2) = self%overall_head%getLabrec()
                iminds     = 0 ! no image associated with the stack header
            else
                do j=1,nr
                    hedinds(1) = cnt+1 ! hed from
                    cnt = cnt+self%overall_head%getLabrec()
                    hedinds(2) = cnt   ! hed to
                    iminds(1) = cnt+1  ! im from
                    cnt = cnt+dims(2)
                    iminds(2) = cnt    ! im to
                end do
            endif
            class DEFAULT
            stop 'Format not supported; slice2recpos; simle_imgfile'
        end select
    end subroutine slice2recpos

    !>  \brief  for translating an image index to record indices in the stack
    !! \param[out] hedinds,iminds indices in the stack
    subroutine slice2bytepos( self, nr, hedinds, iminds )
        class(imgfile), intent(in)     :: self   !< Imagefile object 
        integer, intent(in)            :: nr                    !< num in stack
        integer(kind=8), intent(inout) :: hedinds(2), iminds(2)
        if( nr < 0 )then
            stop 'cannot have negative slice indices; slice2bytepos; simple_imgfile'
        else if( nr == 0 )then
            hedinds(1) = 1
            hedinds(2) = self%overall_head%getLabbyt()
            iminds     = 0 ! no image in SPIDER stack header
        else
            call self%slice2recpos( nr, hedinds, iminds )
            hedinds(1) = (hedinds(1)-1)*self%overall_head%getLenbyt()+1
            hedinds(2) = hedinds(2)*self%overall_head%getLenbyt()
            iminds(1)  = (iminds(1)-1)*self%overall_head%getLenbyt()+1
            iminds(2)  = iminds(2)*self%overall_head%getLenbyt()
        endif
    end subroutine slice2bytepos

    subroutine print_slice2bytepos( self, n )
        class(imgfile), intent(in) :: self
        integer, intent(in)        :: n
        integer :: i
        integer(kind=8) :: hedinds(2), iminds(2)
        do i=1,n
            call self%slice2bytepos(i, hedinds, iminds)
            write(*,*) 'slice: ', i, 'first hedpos: ', hedinds(1), 'last hedpos: ', hedinds(2),&
                 'first imgpos: ', iminds(1), 'last imgpos: ', iminds(2), 'hedsz: ',  hedinds(2)- hedinds(1)+1,&
                 'imgsz: ',  iminds(2)- iminds(1)+1
        end do
    end subroutine print_slice2bytepos

    ! SLICE & HEADER I/O

    !>  \brief  read a slice of the image file from disk into memory
    subroutine rSlice( self, slice_nr, rarr )
        class(imgfile), intent(inout)    :: self        !< Imagefile object 
        integer, intent(in)              :: slice_nr    !< Number of the slice to read in (the first slice in the file is numbered 1)
        real, intent(inout), allocatable :: rarr(:,:,:) !< Array of reals. Will be (re)allocated if needed
        call self%rwSlices('r',slice_nr,slice_nr,rarr)
        DebugPrint  '(imgfile::rSlice) read slice: ', slice_nr
    end subroutine rSlice

    !>  \brief  write a slice of the image file from memory to disk
    subroutine wSlice( self, slice_nr, rarr, ldim )
        class(imgfile), intent(inout) :: self        !< Imagefile object 
        integer, intent(in)           :: slice_nr    !<  Number of the slice to read in (the first slice in the file is numbered 1)
        real, intent(inout)           :: rarr(:,:,:) !<  Array of reals. Will be (re)allocated if needed
        integer, intent(in)           :: ldim(3)     !<  Logical size of the array. This will be written to disk: rarr(1:ldim_1,:,:)
        call self%rwSlices('w',slice_nr,slice_nr,rarr,ldim)
        DebugPrint  '(imgfile::wSlice) wrote slice: ', slice_nr
    end subroutine wSlice

    !>  \brief  reads an image or stack header
    subroutine rHead( self, slice, head, ldim )
        class(imgfile), intent(inout), target :: self      !< Imagefile object 
        integer, intent(in)                   :: slice     !< stack slice or zero for individual
        class(ImgHead), intent(inout)         :: head      !< img head object
        integer, intent(inout), optional      :: ldim(3)   !< for reading the overall stack header of SPIDER files
        class(ImgHead), pointer               :: ptr
        integer(kind=8) :: hedbyteinds(2), imbyteinds(2), first_byte
        ptr => self%overall_head
        if( slice == 0 )then
            select type(ptr)
            type is (MrcImgHead)
                call head%new
            type is (SpiImgHead)
                if( present(ldim) )then
                    call head%new(ldim=ldim)
                else
                    stop 'need ldim for reading the overall stack header of SPIDER files; rHead; simple_imgfile'
                endif
                class DEFAULT
                stop 'Format not supported; slice2recpos; simle_imgfile'
            end select
            ! read overall stack header
            if( debug )then
                call head%read(self%funit, print_entire=.true.)
            else
                call head%read(self%funit)
            endif
        else
            select type(ptr)
            type is (MrcImgHead)
                stop 'Individual images in MRC stacks do not have associated headers; rHead; simple_imgfile'
            type is (SpiImgHead)
                call self%slice2bytepos(slice, hedbyteinds, imbyteinds)
                first_byte = hedbyteinds(1)
                DebugPrint  '(simple_imgfile::rHead) position of first byte: ', first_byte
                if( present(ldim) )then
                    call head%new( ldim=ldim )
                else
                    stop 'need ldim for reading the overall stack header of SPIDER files; rHead; simple_imgfile'
                endif
                call head%new( ldim=ldim )
                ! read image header
                if( debug )then
                    call head%read(self%funit, pos=first_byte, print_entire=.true.)
                else
                    call head%read(self%funit, pos=first_byte)
                endif
                class DEFAULT
                stop 'Format not supported; rHead; simle_imgfile'
            end select
        endif
        DebugPrint  '(simple_imgfile::rHead) read header from file'
    end subroutine rHead

    !>  \brief  writes an image or stack header
    subroutine wHead( self, slice, head )
        class(imgfile), intent(inout), target   :: self    !< Imagefile object 
        integer, intent(in)                     :: slice   !< stack slice
        class(ImgHead), intent(inout), optional :: head    !< img head object
        class(ImgHead), pointer                 :: ptr=>null()
        integer(kind=8) :: hedbyteinds(2), imbyteinds(2), first_byte
        if( slice == 0 )then
            if( present(head) )then
                call head%write(self%funit)
            else
                call self%overall_head%write(self%funit)
            endif
        else
            ptr => self%overall_head
            select type(ptr)
            type is (MrcImgHead)
                stop 'Individual images in MRC stacks do not have associated headers; wHead; simple_imgfile'
            type is (SpiImgHead)
                call self%slice2bytepos(slice, hedbyteinds, imbyteinds)
                first_byte = hedbyteinds(1)
                if( present(head) )then
                    if( same_type_as(head,self%overall_head) )then
                        call head%write(self%funit, pos=first_byte)
                    else
                        stop 'Inconsistent header types; wHead; simple_imgfile'
                    endif
                else
                    call self%overall_head%write(self%funit, pos=first_byte)
                endif
                class DEFAULT
                stop 'unsupported type; wHead; simple_imgfile'
            end select
        endif
    end subroutine wHead

    !>  \brief  read/write a set of contiguous slices of the image file from disk into memory.
    !!          The array of reals should have +2 elements in the first dimension.
    subroutine rwSlices( self, mode, first_slice, last_slice, rarr, ldim, is_ft, smpd, read_failure )
        use simple_strings, only: int2str
        use simple_math,    only: is_odd, is_a_number
#ifdef PGI
        use ISO_C_BINDING
        !include 'lib3f.h'
#else
        use, intrinsic :: iso_c_binding
#endif
        use simple_imghead !, only: ImgHead, SpiImgHead, MrcImgHead, dataRbytes, dataRinteger, dataRfloat
        use simple_math, only: is_even
        class(imgfile), target, intent(inout) :: self         !< instance  Imagefile object 
        character(len=1),       intent(in)    :: mode         !< read (r) or write (w)
        integer,                intent(in)    :: first_slice  !< First slice (the first slice in the file is numbered 1)
        integer,                intent(in)    :: last_slice   !< Last slice
        real,                   intent(inout) :: rarr(:,:,:)  !< Array of reals. Will be (re)allocated if needed
        integer, optional,      intent(in)    :: ldim(3)      !< Logical size of the array. This will be written to disk: rarr(1:ldim(1),:,:)
        logical, optional,      intent(in)    :: is_ft        !< to indicate FT status of image
        real,    optional,      intent(in)    :: smpd         !< sampling distance
        logical, optional,      intent(out)   :: read_failure !< status of file io ops
        real(kind=4),    allocatable :: tmp_32bit_float_array(:,:,:)
        integer(kind=1), allocatable :: tmp_byte_array(:,:,:)
        integer(kind=2), allocatable :: tmp_16bit_int_array(:,:,:)
        character(len=100)           :: io_message
        integer                      :: io_stat,i,j,k,itmp,dims(3)
        integer(kind=8)              :: first_byte,hedbyteinds(2),imbyteinds(2),first_hedbyte, byteperpix
        logical                      :: arr_is_ready, ft_indic
        real                         :: min_val,max_val
        class(ImgHead), pointer      :: ptr=>null()
        class(ImgHead), allocatable  :: imghed
!        character(len=20)            :: conv
#ifdef PGI
        include 'lib3f.h'
#endif
        ! Check that the first and last slice numbers given make sense
        if( first_slice > 0 .and. (first_slice .gt. last_slice) ) stop 'Last < first slice; rwSlices; simple_imgfile'
        if( present(ldim) ) dims = ldim
        ptr => self%overall_head
        select type(ptr)
        type is (MrcImgHead)
            ! Get the dims of the image file
            dims = self%overall_head%getDims()
            DebugPrint  '(rwSlices :: simple_imgfile) dims gotten from overall_head: ', dims(1), dims(2), dims(3)
        type is (SpiImgHead)
            if( self%isvol )then
                ! all good
            else
                if( size(rarr,3) /= 1 )then
                    stop 'third dim of inputted array needs to be 1 for SPIDER stacks; rwSlices; simple_imgfile'
                endif
            endif
            ! Get the dims of the image file
            dims = self%overall_head%getDims()
            DebugPrint  '(rwSlices :: simple_imgfile) dims gotten from overall_head: ', dims(1), dims(2), dims(3)
            class DEFAULT
            stop 'Format not supported; rwSlices; simple_imgfile'
        end select
        if( mode .eq. 'r' )then
            dims(3) = last_slice-first_slice+1
        else if( mode .eq. 'w' )then
            if( present(ldim) )then
                DebugPrint  '(rwSlices :: simple_imgfile) ldims inputted: ', ldim(1), ldim(2), ldim(3)
                ! Check that the array dims and the file dims are compatible.
                select type(ptr)
                type is (MrcImgHead)
                    ! Redefine the file dims (to keep track of the index of the last image of the stack)
                    dims=self%overall_head%getDims()
                    DebugPrint  '(rwSlices :: simple_imgfile) getdims in overall_head: '
                    if(debug) write(*,'("Dim1: ",I3," Dim2: ",I3," Dim3: ",I3)') dims(1), dims(2), dims(3)
                    dims(1) = ldim(1)
                    dims(2) = size(rarr,2)
                    dims(3) = max(last_slice,dims(3))
                    DebugPrint  '(rwSlices :: simple_imgfile) calling setdims overall_head: ', dims(1), dims(2), dims(3)
                    call self%overall_head%setDims(dims)
                    DebugPrint  '(rwSlices :: simple_imgfile) dims set in overall_head: ', dims(1), dims(2), dims(3)
                type is (SpiImgHead)
                    ! Since we need to know the image dimensions in order to create a SPIDER
                    ! file-handler there is no point in redefining the file dims. We should check that
                    ! first_slice == last_slice when we are writing to stacks though
                    if( .not. self%isvol )then
                        if( first_slice /= last_slice )then
                            stop 'first_slice = last_slice required for writing spider&
                                 &images to stack; rwSlices; simple_imgfile'
                        endif
                    endif
                    class default
                    DebugPrint  '(rwSlices :: simple_imgfile)  ptr type not one of MrcImgHead or SpiImgHead'
                end select
            else
                stop 'need logical size of array (ldim); rwSlices; simple_imgfile'
            endif
        else
            stop 'unsupported mode; rwSlices; simple_imgfile'
        endif

        DebugPrint 'Check that the array is properly allocated'
        if( is_odd(dims(1)) )then
            arr_is_ready = size(rarr,1) .eq. dims(1)+1
        else
            arr_is_ready = size(rarr,1) .eq. dims(1)+2
        endif
        DebugPrint 'Is array is properly allocated ',  arr_is_ready

        !arr_is_ready = AND(arr_is_ready , (size(rarr,2) .eq. dims(2)) )
        if( .not. arr_is_ready .or. (size(rarr,2) .ne. dims(2) ))then
            write(*,*) 'Array size: ', size(rarr,1), size(rarr,2), size(rarr,3)
            write(*,*) 'Dimensions: ', dims(1), dims(2), dims(3)
            stop 'Array is not properly allocated; rwSlices; simple_imgfile'
        endif
        byteperpix = int(self%overall_head%bytesPerPix(),kind=8)
        DebugPrint 'Work out the position of the first byte'
        DebugPrint 'byte per pix: ', byteperpix, ' first slice: ',first_slice
        ! Work out the position of the first byte
        select case(self%head_format)
        case('M','F')
            first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
            &*int(product(dims(1:2)),kind=8)*byteperpix
            DebugPrint 'rwSlices; MRC format first byte offset ',first_byte
        case('S')
            if( self%isvol )then
                first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                &*int(product(dims(1:2)),kind=8)*byteperpix
            else
                call self%slice2bytepos(first_slice, hedbyteinds, imbyteinds)
                first_byte = imbyteinds(1)     ! first image byte
                first_hedbyte = hedbyteinds(1) ! first header byte
            endif
            if( first_byte < 1)then
                write(*,*) '(imgfile::rwSlices) self%overall_head%firstDataByte: ', self%overall_head%firstDataByte()
                write(*,*) '(imgfile::rwSlices) self%overall_head%bytesPerPix: ', self%overall_head%bytesPerPix()
                write(*,*) '(imgfile::rwSlices) first dim: ',   dims(1)
                write(*,*) '(imgfile::rwSlices) second dim: ',  dims(2)
                write(*,*) '(imgfile::rwSlices) first_slice: ', first_slice
                write(*,*) '(imgfile::rwSlices) first_byte: ',  first_byte
                write(*,*) '(imgfile::rwSlices) hedbyteinds: ',  hedbyteinds
                write(*,*) '(imgfile::rwSlices) imbyteinds: ',  imbyteinds
                stop
            endif
        case DEFAULT
            stop 'Format not supported; rwSlices; simple_imgfile'
        end select
        DebugPrint 'rwSlice:  process data on disk '
        if( mode .eq. 'r' )then
            DebugPrint 'Reading data rwSlices'
            select case(byteperpix)
            case(1) ! Byte data
                DebugPrint 'Allocating 8-bit array in rwSlices'
                allocate(tmp_byte_array(dims(1),dims(2),dims(3)),stat=alloc_stat)
                if(alloc_stat /= 0) call alloc_errchk("In simple_imgfile:: rwSlices ;  Byte data ", alloc_stat)
                read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_byte_array
                if(io_stat /= 0) call fileio_errmsg("In simple_imgfile:: rwSlices ; Byte data "//trim(io_message), io_stat)
                ! Conversion from unsigned byte integer (which MRC appears to be) is tricky because Fortran doesn't do unsigned integer natively.
                ! The following IAND trick is courtesy of Jim Dempsey at http://software.intel.com/en-us/forums/showthread.php?t=64400
                ! Confusingly, the MRC format documentation implies that one should expect signed integers, which seems to be incorrect: http://www2.mrc-lmb.cam.ac.uk/image2000.html
                ! IMOD documentation indicates that prior to IMOD 4.2.23, unsigned bytes were used and that one needs to inspect the imodStamp head to check
                if( self%overall_head%pixIsSigned() )then
                    rarr(1:dims(1),:,:) = tmp_byte_array(:,:,:)
                else
                    rarr(1:dims(1),:,:) = real(iand(int(tmp_byte_array(:,:,:),kind=4),int(255,kind=4)))
                endif
                deallocate(tmp_byte_array,stat=alloc_stat)
                if(alloc_stat /= 0) call alloc_errchk("In simple_imgfile:: rwSlices ; dealloc tmp_byte_array", alloc_stat )
            case(2) ! 16-bit data
                DebugPrint 'Allocating 32-bit array in rwSlices'
                select case (self%overall_head%getPixType())
                case(dataRinteger)
                    allocate(tmp_16bit_int_array(dims(1),dims(2),dims(3)),stat=alloc_stat)
                    if(alloc_stat /= 0) call alloc_errchk("In simple_imgfile:: rwSlices ; 16-bit data ", alloc_stat )
                    read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_16bit_int_array
                    if(io_stat /= 0) call fileio_errmsg("In simple_imgfile:: rwSlices ; "//trim(io_message), io_stat)
                    if( self%overall_head%pixIsSigned() )then
                        rarr(1:dims(1),:,:) = real(tmp_16bit_int_array(:,:,:))
                    else
                        rarr(1:dims(1),:,:) = real(iand(int(tmp_16bit_int_array(:,:,:),kind=4),&
                             int(huge(int(1,kind=2)),kind=4)))
                    endif
                    deallocate(tmp_16bit_int_array,stat=alloc_stat)
                    if(alloc_stat /= 0) call alloc_errchk("In simple_imgfile:: rwSlices ; tmp_16bit_int_array" , alloc_stat)
                case DEFAULT
                    stop 'Non-integer 16-bit data not supported; rwSlices; simple_imgfile'
                end select
            case(4) ! 32-bit data (SPIDER data always has 4 bytes per pixel)
                DebugPrint 'Allocating 32-bit array in rwSlices'
                allocate(tmp_32bit_float_array(dims(1),dims(2),dims(3)),stat=alloc_stat)
                if(alloc_stat /= 0) call alloc_errchk("In simple_imgfile:: rwSlices ; 32-bit data ", alloc_stat )
                read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_32bit_float_array
                if(io_stat /= 0) call fileio_errmsg("In simple_imgfile:: rwSlices ; 32-bit data "&
                     //trim(io_message), io_stat)
                rarr(1:dims(1),:,:) = tmp_32bit_float_array
                deallocate(tmp_32bit_float_array,stat=alloc_stat)
                if(alloc_stat /= 0) call alloc_errchk("In simple_imgfile:: rwSlices ; 32-bit data ", alloc_stat )

            case DEFAULT
                write(*,'(2a)') 'fname: ', self%fname
                write(*,'(a,i0,a)') 'bit depth: ', self%overall_head%bytesPerPix(), ' bytes'
                stop 'Unsupported bit-depth; rSlices; simple_imgfile'
            end select
            ! Make sure we set the non-used part of the array to 0.0
            if( .not. self%overall_head%pixIsComplex() ) rarr(dims(1)+1:,:,:) = 0.
            ! Check the read was successful
            if( io_stat .ne. 0 )then
                if( present(read_failure) )then
                    read_failure = .true.
                endif
                write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', self%fname
                write(*,'(2a)') 'IO error message was: ', io_message
                stop 'I/O error; rwSlices; simple_imgfile'
            endif
            if( present(read_failure) )then
                read_failure = .false.
                return
            endif
        else  ! WRITE MODE
            ! find minmax
            min_val =  huge(1.)
            max_val = -huge(1.)
            do k=1,size(rarr,3)
                do j=1,size(rarr,2)
                    do i=1,ldim(1)
                        if( is_a_number(rarr(i,j,k)) )then
                            if( rarr(i,j,k) .gt. max_val) max_val = rarr(i,j,k)
                            if( rarr(i,j,k) .lt. min_val) min_val = rarr(i,j,k)
                        endif
                    enddo
                enddo
            enddo
            select case(self%overall_head%bytesPerPix())
            case(4)
                select case(self%head_format)
                case('M','F')
                    ! nothing to do
                case('S')
                    if( .not. self%isvol )then
                        ! for SPIDER stack we also need to create and write an image header
                        allocate( SpiImgHead :: imghed )
                        if( present(ldim) )then
                            if( present(is_ft) )then
                                if( present(smpd) )then
                                    call imghed%new(ldim=ldim)
                                    call imghed%setMinimal(ldim, is_ft, smpd)
                                    call imghed%setMinPixVal(min_val)
                                    call imghed%setMaxPixVal(max_val)
                                else
                                    stop 'need optional parameter sampling distance (smpd) to&
                                         &make SPIDER image header; rwSlices; simple_imgfile'
                                endif
                            else
                                stop 'need optional indicator parameter is_ft to make SPIDER&
                                     &image header; rwSlices; simple_imgfile'
                            endif
                        else
                            stop 'need optional logical dimension (ldim) to make SPIDER image&
                                 &header; rwSlices; simple_imgfile'
                        endif
                        call imghed%write(self%funit, pos=first_hedbyte)
                        call imghed%kill
                        deallocate(imghed)
                    endif
                end select

                if( debug )then
                    print *, 'size(rarr, dim1): ', size(rarr,1)
                    print *, 'size(rarr, dim2): ', size(rarr,2)
                    print *, 'size(rarr, dim3): ', size(rarr,3)
                    print *, 'range of 1st dim: ', 1, dims(1)
                    print *, 'first byte:       ', first_byte
                endif
                tmp_32bit_float_array = rarr(1:dims(1),:,:)
                write(unit=self%funit,pos=first_byte,iostat=io_stat) tmp_32bit_float_array
                deallocate(tmp_32bit_float_array)
            case DEFAULT
                print *, 'bit depth: ', int2str(self%overall_head%bytesPerPix())
                stop 'Unsupported bit-depth; rwSlices; simple_imgfile'
            end select
            ! Check the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(wSlices): I/O error ', io_stat, ' when writing to: ', self%fname
                stop 'I/O error; rwSlices; simple_imgfile'
            endif
            ! May need to update file dims
            select case(self%head_format)
            case('M')
                DebugPrint  '(rwSlices :: simple_imgfile) last_slice: ', last_slice
                dims(3) = max(dims(3),last_slice)
                DebugPrint  '(rwSlices :: simple_imgfile) updated dims to: ', dims(1), dims(2), dims(3)
                call self%overall_head%setDims(dims)
            case('S')
                if( .not. self%isvol )then
                    itmp = self%overall_head%getMaxim()
                    call self%overall_head%setMaxim(max(itmp,last_slice))
                endif
            end select
            if( self%head_format .ne. 'F' )then
                ! May need to update min and max in the overall header head
                if( min_val .lt. self%overall_head%getMinPixVal()) call self%overall_head%setMinPixVal(min_val)
                if( max_val .gt. self%overall_head%getMaxPixVal()) call self%overall_head%setMaxPixVal(max_val)
            endif
            ! May need to update FT status in the overall header head
            select type(ptr)
            type is (MrcImgHead)
                if( present(is_ft) )then
                    if( is_ft ) call self%overall_head%setMode(4)
                endif
            type is (SpiImgHead)
                dims(1) = size(rarr,1)
                dims(2) = size(rarr,2)
                dims(3) = size(rarr,3)
                ft_indic = .false.
                if( present(is_ft) ) ft_indic = is_ft
                if( .not. self%isvol .and. .not. ft_indic )then
                    call self%overall_head%setIform(1)
                else if( self%isvol .and. .not. ft_indic )then
                    call self%overall_head%setIform(3)
                else if( .not. self%isvol .and. ft_indic .and. .not. is_even(dims(1:2)) )then
                    call self%overall_head%setIform(-11)
                else if( .not. self%isvol .and. ft_indic .and. is_even(dims(1:2)) )then
                    call self%overall_head%setIform(-12)
                else if( self%isvol .and. ft_indic .and. .not. is_even(dims(1:2)) )then
                    call self%overall_head%setIform(-21)
                else if( self%isvol .and. ft_indic .and. is_even(dims(1:2)) )then
                    call self%overall_head%setIform(-22)
                else
                    stop 'undefined file type, rwSlices; simple_imgfile'
                endif
            end select
            ! Remember that we wrote to the file
            self%was_written_to = .true.
        endif
        DebugPrint  '(imgfile::rwSlices) completed'
    end subroutine rwSlices

    ! GETTERS, SETTERS, PRINTERS

    !>  \brief  Print out basic information about the file
    subroutine print_imgfile( self )
        class(imgfile), intent(in) :: self   !< Imagefile object 
        write(*,'(/2a)') 'Summary information for file ', trim(adjustl(self%fname))
        call self%overall_head%print_imghead()
        write(*,'(a)') ' '
    end subroutine print_imgfile

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

    !>  \brief Return the fname
    function getFname( self )
        class(imgfile), intent(in) :: self   !< Imagefile object
        character(len=STDLEN)      :: getFname
        getFname = self%fname
    end function getFname

    !>  \brief  Return the pixel size of the image data (in Angstroms)
    real function getPixSz( self )
        class(imgfile), intent(in) :: self
        getPixSz = self%overall_head%getPixSz()
    end function getPixSz

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

    !>  \brief  Set the pixel size of the stack
    subroutine setRMSD( self, dev )
        class(imgfile), intent(inout) :: self !< Imagefile object
        real,           intent(in)    :: dev !< rmsd
        call self%overall_head%setRMSD(dev)
        self%was_written_to = .true.
    end subroutine setRMSD

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


    !>  \brief is for gettign a part of the info in a MRC image header
    subroutine get_mrcfile_info( fname, ldim, form, smpd, doprint )
        use simple_imghead, only: ImgHead, MrcImgHead
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
                    call fclose(filnum,ios,errmsg=" get_mrcfile_info close error "//trim(fname))
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
                    call fclose(filnum, ios,errmsg=" get_mrcfile_info fclose error "//trim(fname))
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
                call fclose(filnum,ios,errmsg=" get_spifile_info fclose error "//trim(fname))
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='NATIVE')
                    call print_spihed
                    return
                endif
                call fopen(filnum, status='OLD', action='READ', file=fname, access='STREAM', iostat=ios)
                call fileio_errmsg(" get_spifile_info fopen error "//trim(fname),ios)
                call read_spihed
                call fclose(filnum,ios,errmsg=" get_spifile_info fclose error "//trim(fname))
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='BIG_ENDIAN')
                    call print_spihed
                    return
                endif
                call fopen(filnum, status='OLD', action='READ', file=fname, &
                         access='STREAM', iostat=ios)
                call fileio_errmsg(" get_spifile_info fopen error "//trim(fname),ios)
                call read_spihed
                call fclose(filnum,iostat=ios,errmsg=" get_spifile_info fclose error "//trim(fname))
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
                    write(*,'(a,3(i0,1x))') 'Number of columns, rows, sections: ', &
                        &int(spihed(12)), int(spihed(2)), int(spihed(1))
                    write(*,'(a,1x,i3)')    'Iform descriptor: ', int(spihed(5))
                    write(*,'(a,1x,f7.0)')  'The number of the highest image currently used in the stack: ', spihed(26)
                    write(*,'(a,1x,f7.3)')  'Pixel size: ', spihed(38)
                endif
            end subroutine

    end subroutine get_spifile_info

end module simple_imgfile
