!>  \brief  Class to deal with image files on disks
!!
!!  The following file formats are (will be) supported:!>  \brief  Class to deal with image files on disks
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
module simple_imgfile
use simple_filehandling
use simple_imghead
use simple_defs
use gnufor2
implicit none

public :: imgfile
private

logical :: debug=.false.
logical :: warn=.false.

type imgfile
    private
    class(imghead), allocatable :: overall_head              !< Overall image head object
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
    procedure          :: print
    procedure          :: getDims
    procedure          :: getDim
    procedure          :: getFname
    procedure          :: getPixSz
    procedure          :: getIform
    procedure          :: getMode
    procedure          :: setIform
    procedure          :: setPixSz
    procedure          :: setMode
end type

contains
    
    ! CONSTRUCTOR
    
    !>  \brief  constructs an imgfile object (file-handle)
    subroutine open( self, fname, ldim, smpd, del_if_exists, formatchar, readhead, rwaction )
        use simple_jiffys, only: alloc_err
        class(imgfile),             intent(inout) :: self          !< Imagefile object to be created
        character(len=*),           intent(in)    :: fname         !< Filename
        integer,                    intent(in)    :: ldim(3)       !< logical dimension of image/stack
        real,                       intent(in)    :: smpd          !< Pixel size of image data (in Angstroms)
        logical,          optional, intent(in)    :: del_if_exists !< If the file already exists on disk, replace it
        character(len=1), optional, intent(in)    :: formatchar    !< file format character flag
        logical,          optional, intent(in)    :: readhead      !< header read indicator
        character(len=*), optional, intent(in)    :: rwaction      !< action flag
        character(len=1) :: format_descriptor
        logical          :: rreadhead
        rreadhead = .true.
        if( present(readhead) ) rreadhead = readhead
        call self%close_nowrite
        ! Remove leading blanks in fname
        self%fname = adjustl(fname)
        if( debug ) print *, 'trying to open: ', self%fname
        ! Work out which file format to use
        if( present(formatchar) )then
            format_descriptor = formatchar
        else
            format_descriptor = fname2format(self%fname)
        endif
        if( format_descriptor .ne. 'N' )then
            self%head_format = format_descriptor
        else
            self%head_format  = default_file_format 
        endif
        if( debug ) print *, 'format: ', self%head_format
        ! Allocate head object
        select case(self%head_format)
            case('F')
                allocate(MrcFeiImgHead :: self%overall_head)
                if( debug ) print *, ' allocated MrcFeiImgHead'
                call self%overall_head%new(ldim)
            case ('M')
                allocate(MrcImgHead :: self%overall_head)
                if( debug ) print *, ' allocated MrcImgHead'
                call self%overall_head%new(ldim)
            case ('S')            
                ! create header
                allocate(SpiImgHead :: self%overall_head)
                if( debug ) print *, ' allocated SpiImgHead'
                call self%overall_head%new(ldim)
                if( ldim(3) > 1 ) self%isvol = .true.          
            case DEFAULT
                stop 'Unsupported file format; new; simple_imgfile' 
        end select
        ! check endconv status
        if( .not. allocated(endconv) )then    
            allocate(endconv, source='NATIVE')
        endif
        ! open the file
        call self%open_local(del_if_exists, rwaction)
        ! read/write the header
        if( file_exists(self%fname) )then
            ! read header
            if( rreadhead )then
                call self%overall_head%read(self%funit)
                if( debug ) print *, 'did read header'
            endif
        else
            ! write header
            call self%overall_head%write(self%funit)
            if( debug ) print *, 'did write header'
        endif
        call self%overall_head%setPixSz(smpd)
        if( debug ) call self%overall_head%print
        ! The file-handle now exists
        self%existence = .true.
        if( debug ) write(*,*) '(imgfile::new) constructed an imgfile object (file-handle)'
    end subroutine open
    
    !>  \brief is forprinting the header
    subroutine print_header( self )
        class(imgfile), intent(in) :: self
        call self%overall_head%print
    end subroutine print_header
    
    ! CORE FUNCTIONALITY
    
    !>  \brief open the file(s) for the imgfile
    subroutine open_local( self, del_if_exists, rwaction )
        class(imgfile),             intent(inout) :: self
        logical, optional,          intent(in) :: del_if_exists
        character(len=*), optional, intent(in) :: rwaction
        character(len=9) :: rw_str
        character(len=7) :: stat_str
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
        self%funit = get_fileunit()
        open(unit=self%funit,access='STREAM',file=self%fname,action=rw_str,status=stat_str,convert=endconv)
        self%was_written_to = .false.
    end subroutine open_local

    !>  \brief  close the file(s) and "de-initialise" the imgfile object
    subroutine close( self )
        class(imgfile), intent(inout) :: self
        integer :: ret
        if( is_open(self%funit) )then
            if( self%was_written_to )then
                call self%overall_head%write(self%funit)                    
                if( debug ) write(*,*) '(simple_imgfile::close) wrote overall_head'
            endif
            close(self%funit)
            call flush(self%funit)
            ret = fsync(fnum(self%funit))
        endif
        if( allocated(self%overall_head) ) call self%overall_head%kill
        if( allocated(self%overall_head) ) deallocate(self%overall_head)
        self%was_written_to = .false.   
        self%existence = .false.
    end subroutine close
    
    !>  \brief  close the file(s)
    subroutine close_nowrite( self )
        class(imgfile), intent(inout) :: self
        integer :: ret
        if( is_open(self%funit) )then
            close(self%funit)
            call flush(self%funit)
            ret = fsync(fnum(self%funit))
        endif
        if( allocated(self%overall_head) ) call self%overall_head%kill
        if( allocated(self%overall_head) ) deallocate(self%overall_head)
        self%was_written_to = .false.   
        self%existence = .false.
    end subroutine close_nowrite

    !>  \brief  Check whether the file exists on disk
    logical function exists( self )
        class(imgfile), intent(inout) :: self
        character(len=STDLEN)         :: new_fname
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
        class(imgfile), intent(in) :: self
        character(len=1)           :: get_format
        get_format = self%head_format
    end function get_format

    !>  \brief  for translating an image index to record indices in the stack
    subroutine slice2recpos( self, nr, hedinds, iminds )
        class(imgfile), intent(in), target :: self 
        integer, intent(in)                :: nr
        integer(kind=8), intent(out)       :: hedinds(2), iminds(2)
        integer                            :: cnt, j, dims(3)
        class(imghead), pointer            :: ptr=>null()
        ptr => self%overall_head
        select type(ptr)
            type is (MrcImgHead)
                stop 'Cannot translate an image index to record indices for MRC files; slice2recpos; simple_imgfile'
            type is (MrcFeiImgHead)
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
    subroutine slice2bytepos( self, nr, hedinds, iminds )
        class(imgfile), intent(in)     :: self
        integer, intent(in)            :: nr
        integer(kind=8), intent(inout) :: hedinds(2), iminds(2)
        integer                        :: cnt, j
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
        class(imgfile), intent(inout)    :: self
        integer, intent(in)              :: slice_nr    !< Number of the slice to read in (the first slice in the file is numbered 1)
        real, intent(inout), allocatable :: rarr(:,:,:) !< Array of reals. Will be (re)allocated if needed
        call self%rwSlices('r',slice_nr,slice_nr,rarr)
        if( debug ) write(*,*) '(imgfile::rSlice) read slice: ', slice_nr
    end subroutine rSlice

    !>  \brief  write a slice of the image file from memory to disk
    subroutine wSlice( self, slice_nr, rarr, ldim )
        class(imgfile), intent(inout) :: self
        integer, intent(in)           :: slice_nr    !<  Number of the slice to read in (the first slice in the file is numbered 1)
        real, intent(inout)           :: rarr(:,:,:) !<  Array of reals. Will be (re)allocated if needed
        integer, intent(in)           :: ldim(3)     !<  Logical size of the array. This will be written to disk: rarr(1:ldim_1,:,:)
        call self%rwSlices('w',slice_nr,slice_nr,rarr,ldim)
        if( debug ) write(*,*) '(imgfile::wSlice) wrote slice: ', slice_nr
    end subroutine wSlice

    !>  \brief  reads an image or stack header
    subroutine rHead( self, slice, head, ldim )
        class(imgfile), intent(inout), target :: self
        integer, intent(in)                   :: slice
        class(imghead), intent(inout)         :: head
        integer, intent(inout), optional      :: ldim(3)
        class(imghead), pointer               :: ptr
        integer(kind=8) :: hedbyteinds(2), imbyteinds(2), first_byte
        ptr => self%overall_head
        if( slice == 0 )then
            select type(ptr)
                type is (MrcImgHead)
                    call head%new
                type is (MrcFeiImgHead)
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
                type is (MrcFeiImgHead)
                    stop 'Individual images in MRC stacks do not have associated headers; rHead; simple_imgfile'
                type is (SpiImgHead)
                    call self%slice2bytepos(slice, hedbyteinds, imbyteinds)
                    first_byte = hedbyteinds(1)
                    if( debug ) write(*,*) '(simple_imgfile::rHead) position of first byte: ', first_byte
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
        if( debug ) write(*,*) '(simple_imgfile::rHead) read header from file'
    end subroutine rHead

    !>  \brief  writes an image or stack header
    subroutine wHead( self, slice, head )
        class(imgfile), intent(inout), target   :: self
        integer, intent(in)                     :: slice
        class(imghead), intent(inout), optional :: head
        class(imghead), pointer                 :: ptr=>null()
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
                type is (MrcFeiImgHead)
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
    subroutine rwSlices( self, mode, first_slice, last_slice, rarr, ldim, is_ft, smpd )
        use simple_strings, only: int2str
        use simple_math, only: is_odd
        use, intrinsic :: iso_c_binding
        use simple_imghead, only: ImgHead, SpiImgHead, MrcImgHead, dataRbytes, dataRinteger, dataRfloat
        use simple_math, only: is_even
        class(imgfile), intent(inout), target :: self         !< instance
        character(len=1), intent(in)          :: mode         !< read (r) or write (w)
        integer, intent(in)                   :: first_slice  !< First slice (the first slice in the file is numbered 1)
        integer, intent(in)                   :: last_slice   !< Last slice
        real, intent(inout)                   :: rarr(:,:,:)  !< Array of reals. Will be (re)allocated if needed
        integer, intent(in), optional         :: ldim(3)      !< Logical size of the array. This will be written to disk: rarr(1:ldim(1),:,:)
        logical, intent(in), optional         :: is_ft        !< to indicate FT status of image
        real, intent(in), optional            :: smpd         !< sampling distance  
        real(kind=4), allocatable             :: tmp_32bit_float_array(:,:,:)
        integer(kind=1), allocatable          :: tmp_byte_array(:,:,:)
        integer(kind=2), allocatable          :: tmp_16bit_int_array(:,:,:)
        character(len=100)                    :: io_message
        integer                               :: io_stat, dims(3),i,j,k,itmp,cnt,ldim_here(3),iform,maxim
        integer(kind=8)                       :: first_byte,hedbyteinds(2),imbyteinds(2),first_hedbyte
        logical                               :: arr_is_ready, ft_indic
        real                                  :: min_val,max_val,smpd_here
        class(imghead), pointer               :: ptr=>null()
        class(imghead), allocatable           :: imghed
        character(len=20)                     :: conv
        ! Check that the first and last slice numbers given make sense
        if( first_slice > 0 .and. (first_slice .gt. last_slice) ) stop 'Last < first slice; rwSlices; simple_imgfile'
        if( present(ldim) ) dims = ldim
        ptr => self%overall_head
        select type(ptr)
            type is (MrcImgHead)
                ! Get the dims of the image file
                dims = self%overall_head%getDims()
                if( debug ) write(*,*) '(rwSlices :: simple_imgfile) dims gotten from overall_head: ', dims(1), dims(2), dims(3)
            type is (MrcFeiImgHead)
                ! all good
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
                if( debug ) write(*,*) '(rwSlices :: simple_imgfile) dims gotten from overall_head: ', dims(1), dims(2), dims(3)
            class DEFAULT
                stop 'Format not supported; rwSlices; simle_imgfile'
        end select
        if( mode .eq. 'r' )then
            dims(3) = last_slice-first_slice+1
        else if( mode .eq. 'w' )then
            if( present(ldim) )then
                if( debug ) write(*,*) '(rwSlices :: simple_imgfile) ldims inputted: ', ldim(1), ldim(2), ldim(3)
                ! Check that the array dims and the file dims are compatible.
                select type(ptr)
                    type is (MrcImgHead)
                        ! Redefine the file dims (to keep track of the index of the last image of the stack)
                        dims=self%overall_head%getDims()
                        dims(1) = ldim(1)
                        dims(2) = size(rarr,2)
                        dims(3) = max(last_slice,dims(3))
                        if( debug ) write(*,*) '(rwSlices :: simple_imgfile) dims set in overall_head: ', dims(1), dims(2), dims(3)
                        call self%overall_head%setDims(dims)
                    type is (MrcFeiImgHead)
                        ! nothing to do
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
                end select
            else
                stop 'need logical size of array (ldim); rwSlices; simple_imgfile'
            endif
        else
            stop 'unsupported mode; rwSlices; simple_imgfile'
        endif
        ! Check so that the array is properly allocated
        if( is_odd(dims(1)) )then
            arr_is_ready = size(rarr,1) .eq. dims(1)+1
        else
            arr_is_ready = size(rarr,1) .eq. dims(1)+2
        endif
        arr_is_ready = arr_is_ready .and. (size(rarr,2) .eq. dims(2))
        if( .not. arr_is_ready )then
            write(*,*) 'Array size: ', size(rarr,1), size(rarr,2), size(rarr,3)        
            write(*,*) 'Dimensions: ', dims(1), dims(2), dims(3)
            stop 'Array is not properly allocated; rwSlices; simple_imgfile'
        endif
        ! Work out the position of the first byte
        select case(self%head_format)
            case('M','F')
                first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                *int(product(dims(1:2)),kind=8)*int(self%overall_head%bytesPerPix(),kind=8)
            case('S')
                if( self%isvol )then
                    first_byte = int(self%overall_head%firstDataByte(),kind=8)+int((first_slice-1),kind=8)&
                    *int(product(dims(1:2)),kind=8)*int(self%overall_head%bytesPerPix(),kind=8)
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
        ! process data on disk
        if( mode .eq. 'r' )then
            select case(self%overall_head%bytesPerPix())
                case(1) ! Byte data
                    allocate(tmp_byte_array(dims(1),dims(2),dims(3)))
                    read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_byte_array
                    ! Conversion from unsigned byte integer (which MRC appears to be) is tricky because Fortran doesn't do unsigned integer natively.
                    ! The following IAND trick is courtesy of Jim Dempsey at http://software.intel.com/en-us/forums/showthread.php?t=64400
                    ! Confusingly, the MRC format documentation implies that one should expect signed integers, which seems to be incorrect: http://www2.mrc-lmb.cam.ac.uk/image2000.html
                    ! IMOD documentation indicates that prior to IMOD 4.2.23, unsigned bytes were used and that one needs to inspect the imodStamp head to check
                    if( self%overall_head%pixIsSigned() )then
                        rarr(1:dims(1),:,:) = tmp_byte_array(:,:,:)
                    else
                        rarr(1:dims(1),:,:) = real(iand(int(tmp_byte_array(:,:,:),kind=4),int(255,kind=4)))
                    endif
                    deallocate(tmp_byte_array)
                case(2) ! 16-bit data
                    select case (self%overall_head%getPixType())
                        case(dataRinteger)
                            allocate(tmp_16bit_int_array(dims(1),dims(2),dims(3)))
                            read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_16bit_int_array
                            if( self%overall_head%pixIsSigned() )then
                                rarr(1:dims(1),:,:) = real(tmp_16bit_int_array(:,:,:))
                            else
                                rarr(1:dims(1),:,:) = real(iand(int(tmp_16bit_int_array(:,:,:),kind=4),&
                                int(huge(int(1,kind=2)),kind=4)))
                            endif
                            deallocate(tmp_16bit_int_array)
                        case DEFAULT
                            stop 'Non-integer 16-bit data not supported; rwSlices; simple_imgfile'
                    end select
                case(4) ! 32-bit data (SPIDER data always has 4 bytes per pixel)
                    allocate(tmp_32bit_float_array(dims(1),dims(2),dims(3)))
                    read(unit=self%funit,pos=first_byte,iostat=io_stat,iomsg=io_message) tmp_32bit_float_array
                    rarr(1:dims(1),:,:) = tmp_32bit_float_array
                    deallocate(tmp_32bit_float_array)
                case DEFAULT
                    write(*,'(2a)') 'fname: ', self%fname
                    write(*,'(a,i0,a)') 'bit depth: ', self%overall_head%bytesPerPix(), ' bytes'
                    stop 'Unsupported bit-depth; rSlices; simple_imgfile'
            end select
            ! Make sure we set the non-used part of the array to 0.0
            if( .not. self%overall_head%pixIsComplex() ) rarr(dims(1)+1:,:,:) = 0.
            ! Check the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', self%fname
                write(*,'(2a)') 'IO error message was: ', io_message
                stop 'I/O error; rwSlices; simple_imgfile'
            endif
            if( self%head_format .eq. 'M' .or. self%head_format .eq. 'F' )then
                ! Is this file from a machine with opposite endinaness?
                if( .not. self%overall_head%hasLocalEndianess() )then
                    stop 'Files created by machines with the opposite endianess are not supported; rwSlices; simple_imgfile'
                endif
            endif
        else
            ! find minmax
            min_val =  huge(1.)
            max_val = -huge(1.)
            do k=1,size(rarr,3)
                do j=1,size(rarr,2)
                    do i=1,ldim(1)
                        if( .not. isnan(rarr(i,j,k)) )then
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
                    write(unit=self%funit,pos=first_byte,iostat=io_stat) rarr(1:dims(1),:,:)
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
                    if( debug ) write(*,*) '(rwSlices :: simple_imgfile) last_slice: ', last_slice
                    dims(3) = max(dims(3),last_slice)
                    if( debug ) write(*,*) '(rwSlices :: simple_imgfile) updated dims to: ', dims(1), dims(2), dims(3)
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
                    ldim_here(1) = size(rarr,1)
                    ldim_here(2) = size(rarr,2)
                    ldim_here(3) = size(rarr,3)
                    ft_indic = .false.
                    if( present(is_ft) ) ft_indic = is_ft
                    if( .not. self%isvol .and. .not. ft_indic )then
                        call self%overall_head%setIform(1)
                    else if( self%isvol .and. .not. ft_indic )then
                        call self%overall_head%setIform(3)
                    else if( .not. self%isvol .and. ft_indic .and. .not. is_even(ldim_here(1:2)) )then
                        call self%overall_head%setIform(-11)
                    else if( .not. self%isvol .and. ft_indic .and. is_even(ldim_here(1:2)) )then
                        call self%overall_head%setIform(-12)
                    else if( self%isvol .and. ft_indic .and. .not. is_even(ldim_here(1:2)) )then
                        call self%overall_head%setIform(-21)
                    else if( self%isvol .and. ft_indic .and. is_even(ldim_here(1:2)) )then
                        call self%overall_head%setIform(-22)
                    else
                        stop 'undefined file type, rwSlices; simple_imgfile'
                    endif
            end select
            ! Remember that we wrote to the file
            self%was_written_to = .true.
        endif
        if( debug ) write(*,*) '(imgfile::rwSlices) completed'
    end subroutine rwSlices

    ! GETTERS, SETTERS, PRINTERS
    
    !>  \brief  Print out basic information about the file
    subroutine print( self )
        class(imgfile), intent(in) :: self
        write(*,'(/2a)') 'Summary information for file ', trim(adjustl(self%fname))
        call self%overall_head%print
        write(*,'(a)') ' '
    end subroutine print

    !>  \brief  Return the dimension of the image stack
    function getDims( self )
        class(imgfile), intent(in) :: self
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
        class(imgfile), intent(in) :: self
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
        class(imgfile), intent(in) :: self
        getIform = self%overall_head%getIform()
    end function getIform
     
    !>  \brief  Return the format descriptor of the stack
    integer function getMode( self )
        class(imgfile), intent(in) :: self
        getMode = self%overall_head%getMode()
    end function getMode

    !>  \brief  Set the format descriptor of the stack
    subroutine setIform( self, iform )
        class(imgfile), intent(inout) :: self
        integer, intent(in)           :: iform
        call self%overall_head%setIform(iform)
        self%was_written_to = .true.
    end subroutine setIform
    
    !>  \brief  Set the pixel size of the stack
    subroutine setPixSz( self, smpd )
        class(imgfile), intent(inout) :: self
        real, intent(in)              :: smpd
        call self%overall_head%setPixSz(smpd)
        self%was_written_to = .true.
    end subroutine setPixSz
    
    !>  \brief  Set the mode of the MRC file
    subroutine setMode( self, mode )
        class(imgfile), intent(inout) :: self
        integer, intent(in)           :: mode
        call self%overall_head%setMode(mode)
        self%was_written_to = .true.
    end subroutine setMode
    
    !>  \brief  for setting the logical dimensions
    subroutine setDims( self, ldim )
        class(imgfile), intent(inout) :: self
        integer, intent(in)           :: ldim(3)
        call self%overall_head%setDims(ldim)
        self%was_written_to = .true.
    end subroutine setDims
    
end module simple_imgfile
