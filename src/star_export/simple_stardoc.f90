
module simple_stardoc
include 'simple_lib.f08'
use simple_star_dict, only: star_dict

implicit none
private
public :: stardoc

type starframes
    type(hash)  :: htab               !< hash table for the parameters
    type(chash) :: chtab              !< hash table for the parameter labels
    integer     :: state
    logical     :: existence =.false. !< to indicate existence
contains
    procedure :: new_frame
    procedure :: kill_frame
end type starframes

interface starframes
    module procedure constructor_frames
end interface starframes

type stardoc
    type(star_dict) :: sdict
    type(str4arr), allocatable :: param_labels(:)
    type(str4arr), allocatable         :: data(:)
    logical      , allocatable         :: param_isstr(:)
    real         , allocatable         :: param_scale(:)
    logical      , allocatable         :: param_converted(:)
    type(starframes), allocatable      :: frames(:)
    integer,public   :: num_data_elements    = 0
    integer,public   :: num_data_lines       = 0
    integer          :: num_frames           = 0
    integer          :: num_frame_elements   = 0
    integer          :: funit
    logical          :: l_open               =.false.
    logical          :: doprint              =.false.
    logical          :: existence            =.false. !< to indicate existence
contains
    procedure        :: new
    procedure        :: read
    procedure,public :: open4import
    procedure,public :: close
    procedure,public :: read_header
    procedure,public :: read_data_lines
    procedure,public :: setdoprint
    procedure        :: write
    procedure        :: get_header
    procedure        :: get_data
    procedure        :: get_r4
    procedure        :: get_i4
    procedure        :: get_str
    !generic         :: get => get_r4, get_i4, get_str
    procedure        :: put_r4
    procedure        :: put_i4
    procedure        :: put_str
    !generic         :: put => put_r4, put_i4, put_str
    procedure        :: kill_stored_params
    procedure,public :: kill_doc
end type stardoc

interface stardoc
    module procedure constructor
end interface stardoc
integer, parameter :: MIN_STAR_NBYTES = 10

enum, bind(C) ! STAR_FORMAT
    enumerator :: STAR_MOVIES=1
    enumerator :: STAR_MICROGRAPHS=2
    enumerator :: STAR_CAVGS=3
    enumerator :: STAR_PTCLS=4
end enum

#include "simple_local_flags.inc"
integer, parameter :: MAX_STAR_ARGS_PER_LINE=16
contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    function constructor_frames( ) result( self )
        type(starframes) :: self
        call self%new_frame
    end function constructor_frames
    !>  \brief  is an abstract constructor
    function constructor( n ) result( self )
        integer, intent(in) :: n
        type(stardoc) :: self
        call self%new(n)
    end function constructor
    subroutine new_frame (self)
        class(starframes), intent(inout) :: self
        call self%kill_frame
        self%htab  = hash()
        self%chtab = chash()
        self%existence = .true.
    end subroutine new_frame

    subroutine new( self, n )
        class(stardoc), intent(inout) :: self
        integer,        intent(in)    :: n
        integer :: i
        call self%kill_doc
        if(.not. self%sdict%exists() ) call self%sdict%new()

        ! if(allocated(self%frames))then
        !     if(self%num_frames > 1) then
        !     do i=1,self%num_frames
        !         call self%frames(i)%kill_frame
        !     end do
        !     endif
        !     deallocate(self%frames)
        ! endif
        ! self%num_frame_elements = n
        ! allocate( self%frames(n), stat=alloc_stat )
        ! if(alloc_stat.ne.0)call allocchk('new; simple_stardoc',alloc_stat)
        ! do i=1,n
        !     call self%frames(i)%new_frame
        ! end do
        DebugPrint ' simple_stardoc::new initiated'
    end subroutine new
    subroutine setdoprint(self)
        class(stardoc), intent(inout):: self
        self%doprint=.true.
    end subroutine setdoprint

    subroutine open4import(self, filename)
        class(stardoc), intent(inout) :: self
        character(len=*),intent(inout) :: filename
        integer :: io_stat, tmpunit
        integer(8) :: filesz
        if(self%l_open) call self%close
        if(allocated(self%param_labels)) deallocate(self%param_labels)
        if(.not. file_exists(trim(filename) ))&
            call simple_stop("simple_stardoc::open ERROR file does not exist "//trim(filename) )
        call fopen(tmpunit, file=trim(filename), action='READ', iostat=io_stat)
        if(io_stat/=0)call fileiochk('star_doc ; open '//trim(filename), io_stat)
        self%funit  = tmpunit
        self%l_open = .true.
        self%existence =.true.

        ! check size
        filesz = funit_size(self%funit)
        if( filesz == -1 )then
            stop 'file_size cannot be inquired; stardoc :: open'
        else if (filesz < MIN_STAR_NBYTES) then
            write(*,*) 'file: ', trim(filename)
            stop 'file size too small to contain a header; stardoc :: open'
        endif

        call self%read_header()
        if(self%num_data_elements > 0 )then
            call self%read_data_lines()
        end if
    end subroutine open4import


    subroutine close( self )
        class(stardoc), intent(inout) :: self !< instance
        integer :: io_stat
        if( self%l_open )then
            call fclose(self%funit,io_stat,errmsg='stardoc ; close ')
            self%l_open = .false.
        end if
    end subroutine close


    subroutine readline(funit, line, ier)
        use iso_fortran_env
        implicit none
        integer, intent(in)                      :: funit
        character(len=:),allocatable,intent(out) :: line
        integer,intent(out)                      :: ier

        integer,parameter                     :: buflen=1024
        character(len=buflen)                 :: buffer
        integer                               :: last
        integer                               :: isize
        logical :: isopened
        line=''
        ier=0
        inquire(unit=funit,opened=isopened,iostat=ier)
        if(ier/=0) call fileiochk("readline isopened failed", ier)
        if(.not. isopened )then
            DiePrint("readline isopened failed")
        endif
        ! read characters from line and append to result
        do
            ! read next buffer (an improvement might be to use stream I/O
            ! for files other than stdin so system line limit is not
            ! limiting)
            read(funit,fmt='(a)',advance='no',size=isize,iostat=ier) buffer
            ! append what was read to result
            !isize=len_trim(buffer)
            if(isize.gt.0)line=line//" "//buffer(:isize)
            ! if hit EOR reading is complete unless backslash ends the line
            !print *, ier == iostat_eor, ier == iostat_end, buffer
            !call fileiochk("readline isopened failed", ier)
            if(is_iostat_eor(ier))then
                last=len(line)
                ! if(last.ne.0)then
                !     ! if line ends in backslash it is assumed a continued line
                !     if(line(last:last).eq.'\')then
                !         ! remove backslash
                !         line=line(:last-1)
                !         ! continue on and read next line and append to result
                !         cycle
                !     endif
                ! endif
                ! hitting end of record is not an error for this routine
                ier=0
                ! end of reading line
                exit
                ! end of file or error
            elseif(ier.ne.0)then
                exit
            endif
        enddo

        line=trim(line)

    end subroutine readline


    !> Parse the STAR file header
    subroutine read_header(self)
        class(stardoc), intent(inout) :: self
        integer          :: n,ios,lenstr,isize,cnt,i, pos1,pos2, nargsline
        character(len=:), allocatable :: line,tmp, starlabel,simplelabel
        logical :: inData, inHeader
        self%num_data_elements=0
        self%num_data_lines=0
        inData=.false.;inHeader=.false.
        if(self%l_open)then
            ios=0
            cnt=1
            n=0
            !! Make sure we are at the start of the file
            rewind( self%funit,IOSTAT=ios)
            if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
            do
                call readline(self%funit, line, ios)
                if(ios /= 0) exit
                line = trim(adjustl(line))
                if(self%doprint) print*, "STAR>> line #",cnt,":", trim(line)
                cnt=cnt+1
                !! Parse the start of the STAR file
                lenstr=len_trim(line)
                if (lenstr == 0 )cycle ! empty line

                !! Count number of fields in header
                if(inHeader)then
                    if (.not. (index(trim(line),"_rln") == 0) )then
                        n=n+1
                        self%num_data_elements=self%num_data_elements+1
                        pos1 = firstNonBlank(trim(line))
                        pos2 = firstBlank(trim(line))
                        DebugPrint " Found STAR field line ",self%num_data_elements, trim(line), " VAR= ",pos2, line(pos1+4:pos2)
                        starlabel = line(pos1+4:pos2)
                        if( self%sdict%isthere(starlabel) )then
                            DebugPrint " STAR field found in dictionary ", n
                            simplelabel = self%sdict%star2simple(starlabel)
                            DebugPrint " STAR field: ", trim(starlabel), " equivalent to ", trim(adjustl(simplelabel))
                        else
                            print *," STAR field: ", trim(starlabel), " star field not in dictionary"

                        end if

                        cycle
                    else
                        inHeader=.false.
                        DebugPrint " End of STAR data field lines "
                        inData = .true.
                    endif
                endif

                !! Count number of data lines
                if(inData)then   !! continue through to end of file
                    self%num_data_lines=self%num_data_lines+1
 !                   DebugPrint " Found STAR data line ", line
                    !call parsestr(line,' ',argline,nargsline)
                    nargsline = cntRecsPerLine(trim(line))
                    if(nargsline /=  self%num_data_elements)then
                        print *, " Records on line mismatch ", nargsline,  self%num_data_elements
                        print *, "Line number ",self%num_data_lines, ":: ", line
                        stop " line has insufficient elements "
                    endif
                    cycle
                end if

                !! does line contain ''data_''
                if ( .not. (index(trim(line),"data_") == 0))then
                    DebugPrint " Found STAR 'data_*' in header ", line
                    !! Quick string length comparison
                    if(lenstr == len_trim("data_pipeline_general"))then
                        print *," Found STAR 'data_pipeline_general' header -- Not supported "
                        exit
                    else if(lenstr == len_trim("data_pipeline_processes"))then
                        print *," Found STAR 'data_pipeline_processes' header -- Not supported "
                        exit
                    else if(lenstr == len_trim("data_sampling_general"))then
                        print *," Found STAR 'data_sampling_general' header -- Not supported "
                        exit
                    else if(lenstr == len_trim("data_optimiser_general"))then
                        print *," Found STAR 'data_optimiser_general' header -- Not supported "
                        exit
                    else if(lenstr == len_trim("data_model_general"))then
                        print *," Found STAR 'data_model_general' header -- Not supported "
                        exit
                    end if
                    !!otherwise
                    cycle
                endif
                if (.not. (index(trim(line), "loop_") == 0) )then
                    inHeader=.true.
                    DebugPrint "Begin STAR field lines ", line
                endif

            end do
            !! End of first pass
            rewind( self%funit,IOSTAT=ios)
            if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
            DebugPrint " STAR field lines: ", self%num_data_elements
            DebugPrint " STAR data lines : ", self%num_data_lines
            if(inHeader) stop "stardoc:: read_header parsing failed"
            inData=.false.; inHeader=.false.
            n=0;cnt=1

            if(allocated(self%param_labels)) deallocate(self%param_labels)
            if(allocated(self%param_isstr)) deallocate(self%param_isstr)
            if(allocated(self%param_scale)) deallocate(self%param_scale)
            if(allocated(self%param_converted)) deallocate(self%param_converted)
            allocate(self%param_labels(self%num_data_elements))
            allocate(self%param_isstr(self%num_data_elements))
            allocate(self%param_scale(self%num_data_elements))
            allocate(self%param_converted(self%num_data_elements))
            self%param_isstr=.false.
            self%param_scale = 1.0
            self%param_converted = .false.
            DebugPrint " Rereading STAR file -- populating params "

            do
                call readline(self%funit, line, ios)
                if(ios /= 0) exit
                line = trim(adjustl(line))
                if(self%doprint) print*, "STAR>> line #",cnt,":", trim(line)
                cnt=cnt+1
                !! Parse the start of the STAR file
                lenstr=len_trim(line)
                if (lenstr == 0 )cycle ! empty line

                !! Process fields in header
                if(inHeader)then
                    if (.not. (index(trim(line),"_rln") == 0) )then
                        n=n+1
                        pos1 = firstNonBlank(trim(line))
                        pos2 = firstBlank(trim(line))
                        DebugPrint " Found STAR field line ", trim(line), " VAR= ",pos2, line(pos1+4:pos2)
                        starlabel = line(pos1+4:pos2)
                        if( self%sdict%isthere(starlabel) )then
                            DebugPrint " STAR field found in dictionary ", n
                            simplelabel = self%sdict%star2simple(starlabel)
                            DebugPrint " STAR field: ", trim(starlabel), " equivalent to ", trim(adjustl(simplelabel))

                            allocate(self%param_labels(n)%str,source=trim(adjustl(simplelabel)))
                            self%param_converted(n) = .true.
                            if (.not.(index(trim(simplelabel),'file' ) == 0 ))then
                                self%param_isstr(n) = .true.
                            else if (.not.(index(trim(simplelabel),'dfx' ) == 0 )  .or.  &
                                &.not.(index(trim(simplelabel),'dfy' ) == 0 ))then
                                self%param_scale(n) = 1000.0
                            endif
                            DebugPrint " STAR param ",n, trim(starlabel), trim(self%param_labels(n)%str),&
                                self%param_scale(n), self%param_isstr(n), self%param_converted(n)
                        else
!                            self%param_labels(n)=" "
                        end if

                        cycle
                    else
                        inHeader=.false.
                        DebugPrint " End of STAR data field lines "
                        exit
                    endif
                endif
            end do
            rewind( self%funit,IOSTAT=ios)
            if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)

        endif
        DebugPrint " End of read_header "
        do i=1, self%num_data_elements
            if(allocated(self%param_labels(i)%str)) then
                print*, i, self%param_labels(i)%str,  self%param_scale(i), self%param_isstr(i), self%param_converted(i)
            else
                print*, i,  self%param_scale(i), self%param_isstr(i), self%param_converted(i)
            endif
        end do

    end subroutine read_header

    subroutine read_data_lines(self)
        class(stardoc), intent(inout) :: self
        integer          :: n, cnt, ios, lenstr, pos1, pos2, i, nargsOnDataline, nDataline, nargsParsed
        character(len=:),allocatable :: line,fname
        character(len=LONGSTRLEN) :: sline
        character(len=STDLEN),allocatable :: lineparts(:)
        logical :: inData, inHeader
        real :: tmpval
        inHeader=.false.;inData=.false.
        n=0;cnt=1;ios=0
        nDataline=0; nargsOnDataline=0

        if(allocated(self%data)) deallocate(self%data)
        allocate(self%data(self%num_data_lines))

        if(self%l_open)then
            do
                call readline(self%funit, line, ios)
                if(ios /= 0) exit
                line=trim(adjustl(line))
                lenstr=len_trim(line)
                !! Skip empty lines
                if (lenstr == 0 )cycle ! empty line

                !! Ignore header lines
                if(inHeader)then
                    if ( .not. (index(trim(line),"_rln") ==0) )then
                        n = n+1
                        cycle
                    else
                        inHeader=.false.
                        DebugPrint " Beginning data lines "
                        inData = .true.
                    endif
                endif

                !! Process Data line
                if(inData)then
                    nDataline=nDataline+1
                    DebugPrint " Found STAR data line ", line
                    nargsOnDataline = cntRecsPerLine(trim(line))
                    if(nargsOnDataline /=  self%num_data_elements) then
                        print *, " Records on line mismatch ", nargsOnDataline,  self%num_data_elements
                        print *, line
                        stop " line has insufficient elements "
                    endif
                    self%data(nDataline)%str = trim(line)
                    if(allocated(lineparts))deallocate(lineparts)
                    allocate(lineparts(nargsOnDataline))
                    call parsestr(line,' ',lineparts,nargsParsed)

                    if(nargsParsed == nargsOnDataline )then
                        DebugPrint " Line Parsing : "
                        !write(*,'(A)', advance='no') " data values ["
                        !do i=1, nargsParsed
                        !    write(*,'(A)', advance='no') trim(adjustl(lineparts(i)))//", "
                        !end do
                        !write(*,*) "]"
                        sline=""
                        do i=1, nargsParsed
                            if ( self%param_isstr(i) ) then
                                if(allocated(fname))deallocate(fname)
                                allocate(fname,source=trim(adjustl(lineparts(i))))
                                !! modify string for :mrc extension
                                if(.not.(index(trim(fname),".ctf:mrc")==0) )then
                                    pos1 = index(trim(fname),":")
                                    fname = trim(fname(1:pos1))
                                endif
                                DebugPrint " Parsing file param : ", trim(fname)
                                if( file_exists (trim(fname)) )then
                                    print *," Found ",trim(self%param_labels(i)%str)," file: ", trim(fname)
                                else
                                    print *," Unable to find file ", trim(lineparts(i))
                                    stop
                                endif

                            else
                                if( self%param_converted (i) )then
                                    !! Data to be inserted
                                    tmpval = str2real(trim(adjustl(lineparts(i))))
                                    tmpval = tmpval * self%param_scale(i)
                                    sline  = sline//trim(self%param_labels(i)%str)//"="//real2str(tmpval)//" "
                                else
                                    DebugPrint " ignoring ",i, trim(self%param_labels(i)%str), trim(adjustl(lineparts(i)))
                                endif
                            end if
                        end do
                        DebugPrint " Processed Data Line : ", trim(sline)

                    else
                        print *, " Parsing records on line failed ", nargsOnDataline,  nargsParsed
                        print *, line
                        stop " line has insufficient elements "
                    endif

                    cycle
                endif

                !! Parse the start of the STAR file
                if ( .not. (index(trim(line), "data_") == 0)) cycle
                if ( .not. (index(trim(line) , "loop_") == 0))then
                    inHeader=.true.
                    DebugPrint "Begin STAR field lines ", line
                endif
            end do
            if(nDataline /= self%num_data_lines)then
                print *," Num data lines mismatch in read_data_lines and read_header"
                stop
            endif

            rewind( self%funit,IOSTAT=ios)
            if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
        endif
    end subroutine read_data_lines

    subroutine  write(self, filename, sp, vars)
        use simple_sp_project
        class(stardoc), intent(inout) :: self
        class(sp_project), intent(inout) :: sp
        character(len=KEYLEN) :: vars(:)          ! fixed-length strings
        character(len=*),intent(inout) :: filename
        character(len=LONGSTRLEN*8) :: starline
        character(len=LONGSTRLEN) :: imagename
        character(len=:), allocatable ::  val,stackfile,state
        integer  :: i, io_stat, iframe, ielem
        real     :: statef, defocus

        if(self%l_open)call self%close
        call fopen(self%funit, trim(filename) ,iostat=io_stat)
        if(io_stat/=0) call fileiochk("In stardoc; write; unable to open "//trim(filename))
        !! standard mode only -- no relion data_...
        write(self%funit,'(A)') ""
        write(self%funit,'(A)') "data_"
        write(self%funit,'(A)') ""
        write(self%funit,'(A)') "loop_"

        do i=1, self%num_data_elements
            write(self%funit,'("_rln",A,3x,"#",I0)') trim(self%param_labels(i)%str), i
        end do
        do iframe=1, self%num_frames
            statef = self%frames(iframe)%state
            if (statef .ne. 0)then

                !! create zero-padded frame number and stackfile
                imagename=int2str_pad(self%get_i4(i,"frameid"),5)//'@'//&
                    self%get_str(iframe,"stackfile")//'s'

                stackfile = self%get_str(i,"stackfile")
                call syslib_symlink(trim(stackfile), trim(stackfile)//'s',status=io_stat)
                if(io_stat/=0)call simple_stop("simple_stardoc::write symlink failed")

                ! starline =imagename//&
                !     &" "//real2str(self%get_r4(i,"kv"))//&
                !     &" "//real2str(self%get_r4(i,"dfx")*1000)//&
                !     &" "//real2str(self%get_r4(i,"dfy")*1000)//&
                !     &" "//real2str(self%get_r4(i,"angast"))//&
                !     &" "//real2str(self%get_r4(i,"cs"))//&
                !     &" "//real2str(self%get_r4(i,"fraca"))//&
                !     &" 10000"//&
                !     &" "//real2str(self%get_r4(i,"smpd"))

                ! do ielem = 2,  self%num_data_elements
                !     write(starline,'(a," ",F13.6)')  starline, self%get(iframe,
                ! end do

                write(self%funit,'(A)') starline
            end if
        end do
    end subroutine write

    subroutine read(self, filename)
        class(stardoc), intent(inout) :: self
        character(len=*),intent(inout) :: filename
        character(len=LONGSTRLEN*8) :: line
        character(len=LONGSTRLEN) :: imagename
        integer :: starfd, io_stat, i, n, endPos, num_labels, elems, iframe, idata
        character(len=LONGSTRLEN)     :: args(MAX_STAR_ARGS_PER_LINE)
        logical :: inFrame, inData
        inData=.false.
        call fopen(starfd, trim(filename),iostat=io_stat)
        if(io_stat/=0) call fileiochk("In stardoc; write; unable to open "//trim(filename))
        self%num_data_elements = 0
        self%num_frames =0; self%num_frame_elements =0
        n = 0; num_labels=0
        do
            read(starfd, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            if(inData)then
                if (line == "" ) inData=.false.
                self%num_data_elements=self%num_data_elements+1
                DebugPrint " Found STAR data line ", line
            else
            if (line == "" )cycle
            if (line == "data_")cycle
            if (line == "loop_" )then
                inData=.true.
                self%num_frames =  self%num_frames +1
            else if (line(1:4) == "_rln")then
                DebugPrint " Found STAR label ", line
                num_labels = num_labels + 1
            endif
            endif
        end do
        rewind(starfd)
        if(self%doprint) then
            print*, "STAR>> Number of labels: ", num_labels
            print*, "STAR>> Number of frames: ",  self%num_frames
            print*, "STAR>> Number of data elements: ", self%num_data_elements
        end if

        if(self%num_frames < 1 .or.  num_labels < 1 .or. self%num_data_elements < 1)&
            call simple_stop(" Star file error -- no frames ")
        call self%new(self%num_frames)
        self%num_frame_elements = num_labels
        inData=.false.
        n=1; elems=0;iframe=1; idata=1
        do
            read(starfd, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            if (line == "" )cycle
            if (line == "data_")cycle
            if (line == "loop_" )then
                inData=.true.
            else if (line(1:4) == "_rln")then
                endPos = firstBlank(line)
                if(endPos <= 5) call simple_stop(" Star file error "//trim(line))
                !call self%datatab%push(line(5:endPos), "#"//int2str(n))
                n = n + 1
            else
                ! Data line
                elems = cntRecsPerLine(line)
                if(elems /= num_labels ) &
                    call simple_stop(" Star file error: num elements in line "//int2str(elems)//&
                    " does not match data req "//int2str(num_labels)//" line:"//trim(line))
                call parsestr(line,' ',args,elems)
                if(self%doprint) then
                    print*, "STAR>> Parsing line: num of elems ", elems
                    print*, "STAR>> Parsing line:  ", args
                end if
                do i=1, elems
                    call self%frames(iframe)%chtab%push(trim(args(i)),int2str(idata))
                end do
                idata = idata + 1

            endif
        end do
        call fclose(starfd)
    end subroutine read

    subroutine putdata(self, strarr)
        class(stardoc), intent(inout) :: self
        character(len=*), intent(in) :: strarr(:)
        integer :: i
        self%num_data_elements = size(strarr)
        do i=1, self%num_data_elements
            !call self%datatab%push(strarr(i), int2str(i))
        end do
        !allocate(self%data(self%num_data_elements))
        !do i=1, self%num_data_elements
        !    allocate(self%data(i)%str, source=strarr(i))
        ! end do
    end subroutine putdata


    subroutine put_r4 (self, iframe, key, val)
        class(stardoc), intent(inout) :: self
        integer, intent(in) :: iframe
        character(len=*), intent(in) :: key
        real, intent(in) :: val
        call self%frames(iframe)%htab%set(key, val)
    end subroutine put_r4


    subroutine put_i4 (self, iframe, key, val)
        class(stardoc), intent(inout) :: self
        integer, intent(in) :: iframe
        character(len=*), intent(in) :: key
        integer, intent(in) :: val
        call self%frames(iframe)%htab%set(key, REAL(val))
    end subroutine put_i4

    subroutine put_str (self, iframe, key, val)
        class(stardoc), intent(inout) :: self
        integer, intent(in) :: iframe
        character(len=*), intent(in) :: key
        character(len=LONGSTRLEN) , intent(in) :: val
        call self%frames(iframe)%chtab%set(key, val)
    end subroutine put_str

    function get_header(self,i) result(val)
        class(stardoc), intent(inout) :: self
        integer, intent(in) :: i
        character(len=:),allocatable :: val
        allocate(val,source=trim(self%param_labels(i)%str))
    end function get_header

    function get_data(self,i) result(val)
        class(stardoc), intent(inout) :: self
        integer, intent(in) :: i
        character(len=:),allocatable :: val
        allocate(val,source=trim(self%data(i)%str))
    end function get_data

    function get_r4 (self, iframe, key) result(val)
        class(stardoc), intent(inout) :: self
        character(len=*), intent(in) :: key
        integer, intent(in) :: iframe
        real :: val
        val = self%frames(iframe)%htab%get(key)
    end function get_r4

    function get_i4 (self, iframe, key) result(val)
        class(stardoc), intent(inout) :: self
        character(len=*), intent(in) :: key
        integer, intent(in) :: iframe
        integer :: val
        val = INT(self%frames(iframe)%htab%get(key))
    end function get_i4

    function get_str (self, iframe, key) result(val)
        class(stardoc), intent(inout) :: self
        character(len=*), intent(in) :: key
        integer, intent(in) :: iframe
        character(len=LONGSTRLEN) :: val
        val = self%frames(iframe)%chtab%get(key)
    end function get_str

    subroutine kill_frame( self )
        class(starframes), intent(inout) :: self
        if( self%existence )then
            call self%htab%kill
            call self%chtab%kill
            self%existence = .false.
        endif
    end subroutine kill_frame


    subroutine kill_stored_params(self)
        class(stardoc), intent(inout) :: self
        integer :: i
        if(allocated(self%param_labels))then
            do i=1,self%num_data_elements
                if(allocated(self%param_labels(i)%str))  deallocate(self%param_labels(i)%str)
            end do
            deallocate(self%param_labels)
        end if
        if(allocated(self%data))then
            do i=1,size(self%data)
                if(allocated(self%data(i)%str))  deallocate(self%data(i)%str)
            end do
            deallocate(self%data)
        end if
        if(allocated(self%param_isstr)) deallocate(self%param_isstr)
        if(allocated(self%param_converted)) deallocate(self%param_converted)
        if(allocated(self%param_scale)) deallocate(self%param_scale)
    end subroutine kill_stored_params


    !>  \brief  is a destructor
    subroutine kill_doc( self )
        class(stardoc), intent(inout) :: self
        integer :: i
        if( self%sdict%exists() ) call self%sdict%kill
        call self%kill_stored_params()
        if( allocated(self%frames) )then
            do i=1,self%num_frames
                call self%frames(i)%kill_frame
            end do
            deallocate( self%frames , stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: kill, module: simple_stardoc ')
            self%num_frame_elements = 0
        endif
        if(self%l_open) call self%close
        self%existence =.false.
    end subroutine kill_doc

end module simple_stardoc
