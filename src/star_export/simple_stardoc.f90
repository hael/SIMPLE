
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
    type(str4arr), allocatable         :: param_starlabels(:)
    type(str4arr), allocatable         :: param_labels(:)
    type(str4arr), allocatable         :: data(:)
    logical      , allocatable         :: param_isstr(:)
    real         , allocatable         :: param_scale(:)
    integer      , allocatable         :: param_converted(:)
    type(starframes), allocatable      :: frames(:)
    character(len=:), allocatable      :: current_file
    integer,public   :: num_valid_elements   = 0
    integer,public   :: num_data_elements    = 0
    integer,public   :: num_data_lines       = 0
    integer          :: num_frames           = 0
    integer          :: num_frame_elements   = 0
    integer          :: funit
    logical          :: l_open               =.false.
    logical          :: doprint              =.false.
    logical          :: existence            =.false. !< to indicate existence
    !! optional params for class and frame numbers
    integer         , allocatable         :: data_framenum(:)
    integer         , allocatable         :: data_classnum(:)
contains
    procedure        :: new
    procedure        :: read
    procedure,public :: open4import
    procedure,public :: close
    procedure,public :: read_header
    procedure,public :: read_data_lines
    procedure,public :: setdoprint
    procedure        :: print
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
    procedure,public :: kill
end type stardoc

interface stardoc
    module procedure constructor
end interface stardoc

enum, bind(C) ! STAR_FORMAT
enumerator :: STAR_MOVIES=1
enumerator :: STAR_MICROGRAPHS=2
enumerator :: STAR_CAVGS=3
enumerator :: STAR_PTCLS=4
end enum

#include "simple_local_flags.inc"
integer, parameter :: MAX_STAR_ARGS_PER_LINE=16
integer, parameter :: MIN_STAR_NBYTES = 10

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
        call self%kill()
        if(.not. self%sdict%exists() ) call self%sdict%new()
     
        if(allocated(self%param_starlabels)) call self%kill_stored_params
        if(allocated(self%param_labels)) call self%kill_stored_params
        if(allocated(self%param_isstr)) call self%kill_stored_params
        if(allocated(self%param_scale)) call self%kill_stored_params
        if(allocated(self%param_converted)) call self%kill_stored_params
     
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
     
        if(.not. file_exists(trim(filename) ))&
            call simple_stop("simple_stardoc::open ERROR file does not exist "//trim(filename) )
        call fopen(tmpunit, file=trim(filename), action='READ', iostat=io_stat)
        if(io_stat/=0)call fileiochk('star_doc ; open '//trim(filename), io_stat)
        self%funit  = tmpunit
        self%l_open = .true.
        self%existence =.true.
        if(.not. self%sdict%exists() ) call self%sdict%new()
        if(allocated(self%current_file)) deallocate(self%current_file)
        allocate(self%current_file,source=trim(adjustl(filename)))
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
            if(allocated(self%current_file)) deallocate(self%current_file)
        end if
    end subroutine close
     
     
    subroutine readline(funit, line, ier)
        use iso_fortran_env
        implicit none
        integer, intent(in)                      :: funit
        character(len=:),allocatable,intent(out) :: line
        integer,intent(out)                      :: ier
        integer,parameter     :: buflen = 1024
        character(len=buflen) :: buffer
        integer               :: last
        integer               :: isize
        logical               :: isopened
        line=''
        ier=0
        inquire(unit=funit,opened=isopened,iostat=ier)
        if(ier/=0) call fileiochk("readline isopened failed", ier)
        if(.not. isopened )then
            HALT_NOW("readline isopened failed")
        endif
        ! read characters from line and append to result
        do
            ! read next buffer (an improvement might be to use stream I/O
            ! for files other than stdin so system line limit is not
            ! limiting)
            read(funit,fmt='(a)',advance='no',size=isize,iostat=ier) buffer
            ! append what was read to result
            if(isize.gt.0)line=line//" "//buffer(:isize)
            ! if hit EOR reading is complete unless backslash ends the line
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
        character(len=:), allocatable :: line,tmp, starlabel,simplelabel
        integer :: ielem,ivalid,idata,ios,lenstr,isize,cnt,i
        integer :: pos1,pos2, nargsline
        logical :: inData, inHeader
        self%num_data_elements=0
        self%num_valid_elements=0
        self%num_data_lines=0
        inData=.false.;inHeader=.false.
        if(.not.self%l_open)then
            print *, " stardoc module read_header file not opened"
            stop " simple_stardoc::read_header "
        endif
        ios=0
        cnt=1
        !! Make sure str4arr params are not allocated
        call self%kill_stored_params
     
        !! First Pass
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
                    self%num_data_elements=self%num_data_elements+1
                    pos1 = firstNonBlank(trim(line))
                    pos2 = firstBlank(trim(line))
                    DebugPrint " Found STAR field line ",self%num_data_elements, trim(line), " VAR= ",pos2, line(pos1+4:pos2)
                    starlabel = trim(adjustl(line(pos1+4:pos2)))
                    if( self%sdict%isthere(starlabel) )then
                        self%num_valid_elements=self%num_valid_elements+1
                        simplelabel = self%sdict%star2simple(starlabel)
                        DebugPrint " STAR field :", self%num_data_elements, " found in dictionary ", trim(starlabel)
                        DebugPrint " New Converted Parameter #",self%num_valid_elements, " string:", trim(adjustl(simplelabel))
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
            !! Count number of  data lines
            if(inData)then   !! continue through to end of file
                self%num_data_lines = self%num_data_lines + 1
                cycle
            end if
     
            ! !! does line contain ''data_''
            ! if ( .not. (index(trim(line),"data_") == 0))then
            !     DebugPrint " Found STAR 'data_*' in header ", line
            !     !! Quick string length comparison
            !     if(lenstr == len_trim("data_pipeline_general"))then
            !         print *," Found STAR 'data_pipeline_general' header -- Not supported "
            !         exit
            !     else if(lenstr == len_trim("data_pipeline_processes"))then
            !         print *," Found STAR 'data_pipeline_processes' header -- Not supported "
            !         exit
            !     else if(lenstr == len_trim("data_sampling_general"))then
            !         print *," Found STAR 'data_sampling_general' header -- Not supported "
            !         exit
            !     else if(lenstr == len_trim("data_optimiser_general"))then
            !         print *," Found STAR 'data_optimiser_general' header -- Not supported "
            !         exit
            !     else if(lenstr == len_trim("data_model_general"))then
            !         print *," Found STAR 'data_model_general' header -- Not supported "
            !         exit
            !     end if
            !     !!otherwise
            !     cycle
            ! endif
            if (.not. (index(trim(line), "loop_") == 0) )then
                inHeader=.true.
                DebugPrint "Begin STAR field lines ", line
            endif
     
        end do
     
        !! End of first pass
        DebugPrint " STAR fields: ", self%num_data_elements, " detected with ", self%num_valid_elements, " in dictionary "
        DebugPrint " STAR data lines : ", self%num_data_lines
        if(inHeader) stop "stardoc:: read_header parsing failed"
        inData=.false.; inHeader=.false.
        cnt=1;
        allocate(self%param_starlabels(self%num_data_elements))
        allocate(self%param_labels(self%num_data_elements))
        allocate(self%param_isstr(self%num_data_elements))
        allocate(self%param_scale(self%num_data_elements))
        allocate(self%param_converted(self%num_data_elements))
        self%param_isstr=.false.
        self%param_scale = 1.0
        self%param_converted = 0
        DebugPrint " Rereading STAR file -- populating params "
        !! Second pass
        ielem=0;ivalid=0;idata=0
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
     
            !! Process fields in header
            if(inHeader)then
                if (.not. (index(trim(line),"_rln") == 0) )then
                    ielem=ielem+1
                    pos1 = firstNonBlank(trim(line))
                    pos2 = firstBlank(trim(line))
                    DebugPrint " Found STAR field line ", trim(line), " VAR= ",pos2, line(pos1+4:pos2)
                    starlabel = trim(adjustl(line(pos1+4:pos2)))
                    allocate(self%param_starlabels(ielem)%str,source=trim(adjustl(starlabel)))
                    if( self%sdict%isthere(starlabel) )then
                        ivalid=ivalid+1
                        DebugPrint " STAR field found in dictionary ", ivalid, " of ", self%num_valid_elements
     
                        simplelabel = self%sdict%star2simple(starlabel)
                        allocate(self%param_labels(ielem)%str,source=trim(adjustl(simplelabel)))
                        self%param_converted(ielem) = 1
                        if (.not.(index(trim(simplelabel),'file' ) == 0 ))then
                            self%param_isstr(ielem) = .true.
                        else if (.not.(index(trim(simplelabel),'dfx' ) == 0 )  .or.  &
                            &.not.(index(trim(simplelabel),'dfy' ) == 0 ))then
                            self%param_scale(ielem) = 0.0001  !! factor microns / angstroms == 1/10^4
                        endif
                        DebugPrint " STAR param #",ielem, trim(starlabel), " ==> #",ivalid, trim(self%param_labels(ielem)%str)
                    else
                        DebugPrint " STAR field: ", trim(starlabel), " star field not in dictionary"
                        allocate(self%param_labels(ielem)%str,source=trim(adjustl(starlabel)))
                    end if
                    cycle
                else
                    inHeader=.false.
                    DebugPrint " End of STAR data field lines "
                    inData=.true.
                endif
            endif
     
            !! Count number of elements in data lines
            if(inData)then   !! continue through to end of file
                idata=idata+1
                !                   DebugPrint " Found STAR data line ", line
                !call parsestr(line,' ',argline,nargsline)
                nargsline = cntRecsPerLine(trim(line))
                if(nargsline /=  self%num_data_elements)then
                    print *, " Records on line mismatch ", nargsline,  self%num_data_elements
                    print *, "Line number ",idata, ":: ", line
                    stop " line has insufficient elements "
                endif
                cycle
            end if
     
            if (.not. (index(trim(line), "loop_") == 0) )then
                inHeader=.true.
                DebugPrint "Begin STAR field lines ", line
            endif
     
        end do
        DebugPrint " End of second pass "
        rewind( self%funit,IOSTAT=ios)
        if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
     
        if(ielem /= self%num_data_elements .or. ivalid /= self%num_valid_elements)then
            print *, " Second pass incorrectly counted  detected parameters "
            HALT_NOW( "simple_stardoc:: read_header invalid element mismatch 1")
        endif
        if(idata /= self%num_data_lines)then
            print *, " Second pass incorrectly counted data lines "
            HALT_NOW( "simple_stardoc:: read_header invalid data line mismatch ")
        endif
        if(sum(int(self%param_converted)) /= self%num_valid_elements) then
            print *, " Second pass incorrectly allocated detected parameters "
            HALT_NOW( "simple_stardoc:: read_header invalid element mismatch 2")
        endif
        !! Check allocation of type variables
        if(.not. allocated(self%param_starlabels) .or. &
            .not. allocated(self%param_labels) .or. &
            .not. allocated(self%param_isstr) .or. &
            .not. allocated(self%param_scale) .or. &
            .not. allocated(self%param_converted)) then
            print *, " Second pass do not allocate correctly "
            HALT_NOW( "simple_stardoc:: read_header allocation unsuccessful")
        endif
     
        print *, ">>>  STAR IMPORT HEADER INFO "
        do i=1, self%num_data_elements
            if(allocated(self%param_labels(i)%str)) then
                if(self%param_converted(i) == 1)then
                    print*, i, trim(self%param_starlabels(i)%str),  trim(self%param_labels(i)%str), &
                        &self%param_scale(i), self%param_isstr(i), self%param_converted(i)
                else
                    print*, i, trim(self%param_starlabels(i)%str), " not converted"
                end if
            else
                print *, " Second pass incorrectly allocated detected parameter label ",i
                HALT_NOW( "simple_stardoc:: read_header invalid label ")
            endif
        end do
     
        DebugPrint " End of read_header "
    end subroutine read_header
     
    subroutine read_data_lines(self)
        class(stardoc), intent(inout) :: self
        integer          :: n, cnt, ios, lenstr, pos1, pos2, i, nargsOnDataline, nDataline, nargsParsed
        character(len=:),allocatable :: line,fname,tmp, projrootdir, experrootdir
        character(len=:),allocatable :: sline
        character(len=STDLEN),allocatable :: lineparts(:)
        logical, allocatable :: fnameselected(:)
        integer :: filetabunit, oritabunit, filetabunit2, imagenamefunit,ctfimagefunit
        logical :: inData, inHeader
        real :: tmpval
        inHeader=.false.;inData=.false.
        n=0;cnt=1;ios=0
        nDataline=0; nargsOnDataline=0
     
        if(.not.self%l_open)then
            print *, " stardoc module read_header file not opened"
            stop " simple_stardoc::read_header "
        endif
        if(.not.allocated(experrootdir)) experrootdir= get_fpath(trim(self%current_file))
        if(allocated(self%data)) deallocate(self%data)
        allocate(self%data(self%num_data_lines))
        if(file_exists('filetab-stardoc.txt')) call del_file("filetab-stardoc.txt")
        call fopen(filetabunit,file="filetab-stardoc.txt")
        if(file_exists('oritab-stardoc.txt')) call del_file("oritab-stardoc.txt")
        call fopen(oritabunit,file="oritab-stardoc.txt") 
        if(file_exists('filetab-stardoc2.txt')) call del_file("filetab-stardoc2.txt")
        call fopen(filetabunit2,file="filetab-stardoc2.txt")
        if(file_exists('imagename-stardoc.txt')) call del_file("imagename-stardoc.txt")
        call fopen(imagenamefunit,file="imagename-stardoc.txt")
        if(file_exists('ctfimage-stardoc.txt')) call del_file("ctfimage-stardoc.txt")
        call fopen(ctfimagefunit,file="ctfimage-stardoc.txt")
        allocate(fnameselected(self%num_data_lines),source=.false.)
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
                    DebugPrint " Line Parsing : ", trim(self%data(nDataline)%str)
                    !write(*,'(A)', advance='no') " data values ["
                    !do i=1, nargsParsed
                    !    write(*,'(A)', advance='no') trim(adjustl(lineparts(i)))//", "
                    !end do
                    !write(*,*) "]"
                    if(allocated(sline))deallocate(sline)
                    do i=1, nargsParsed
                        if( self%param_converted(i) == 0 )then
                            DebugPrint " ignoring ",i, trim(self%param_labels(i)%str)," ", trim(adjustl(lineparts(i)))
                            cycle
                        endif
                        !! Check filename params point to actual files
                        if ( self%param_isstr(i) ) then
                            if(allocated(fname))deallocate(fname)
                            allocate(fname,source=trim(adjustl(lineparts(i))))
     
                            !! check filename for irregular STAR behaviour
                            pos1 = index(trim(fname),"@")
                            if( pos1 /= 0 )then
                                if(pos1 > len_trim(fname) - 1)then
                                    HALT_NOW(" File part "//trim(fname)//" irregular"  )
                                endif
                                !! get the frame number
                                if(.not.allocated(self%data_framenum)) allocate(self%data_framenum(self%num_data_lines))
                                call str2int( trim(fname(1:pos1-1)) , ios, self%data_framenum(nDataline) )
                                if(ios/=0) then
                                    HALT_NOW( 'stardoc; read_datalines frame num str2int error :'//trim(fname(1:pos1-1)) )
                                endif
                                !! shift the filename string
                                fname = trim(fname(pos1+1:))
                                write(imagenamefunit,'(A,1x, I0)') trim(fname), self%data_framenum(nDataline)
                            endif
                            !! modify string for :mrc extension
                            !! Special note: CtfImage files are recorded as "*.ctf:mrc". These are
                            !! not imported into SIMPLE, but we still check to see if they exist
                            pos1 = index(trim(fname),":mrc",back=.true.)
                            if( pos1 /= 0 )then
                                if(pos1 < 4 .or. pos1 > len_trim(fname) - 1)then
                                    HALT_NOW(" File part "//trim(fname)//" irregular"  )
                                endif
                                fname = trim(fname(1:pos1-1))
                            endif
                            if(.not. allocated(projrootdir))then
                                !! Estimate Project root dir
                                call simple_abspath(filepath(get_fpath(self%current_file),PATH_PARENT,PATH_PARENT), projrootdir)
                                DebugPrint " Project root dir (estimate from ../../ of starfile) : ", trim(projrootdir)
                            end if
     
                            DebugPrint " Parsing file param : ", trim(fname)
                            if( file_exists (trim(fname)) )then
                                print *," Found ",trim(self%param_labels(i)%str)," file: ", trim(fname)
                            else if( file_exists (trim(fname)//"s") )then
                                !! Account for MRCS extension
                                print *," Found ",trim(self%param_labels(i)%str)," file: ", trim(fname)//"s"
                                fname = trim(fname)//"s"
                            else
                                tmp = get_fpath(trim(self%current_file))
                                if( allocated(tmp) )then
                                    tmp = filepath(trim(tmp), trim(fname))
                                    if( file_exists (trim(tmp)) )then
                                        print *," Found ",trim(self%param_labels(i)%str)," file: ", trim(tmp)
                                        fname = trim(tmp)
                                    else if( file_exists (trim(tmp)//"s") )then
                                        !! Account for MRCS extension
                                        print *," Found ",trim(self%param_labels(i)%str)," file: ", trim(tmp)//"s"
                                        fname = trim(tmp)//"s"
                                    else if ( allocated(projrootdir)) then
                                        tmp = filepath(trim(projrootdir), trim(fname))
                                        if( file_exists (trim(tmp)) )then
                                            print *," Found ",trim(self%param_labels(i)%str)," file: ", trim(tmp)
                                            fname = trim(tmp)
                                        else if( file_exists (trim(tmp)//"s") )then
                                            !! Account for MRCS extension
                                            print *," Found ",trim(self%param_labels(i)%str)," file: ", trim(tmp)//"s"
                                            fname = trim(tmp)//"s"
                                        else
                                            print *," Unable to find file in path of project "//trim(tmp)
                                            HALT_NOW(" Unable to find file "//trim(tmp) )
                                        endif
                                    else
                                        print *," Unable to find file in path of current star file "//trim(tmp)
                                        HALT_NOW(" Unable to find file "//trim(tmp) )
                                    endif
                                else
                                    print *," Unable to find file or path of current star file "//trim(fname)
                                    HALT_NOW("simple_stardoc failed to get file in starfile")
                                endif
                            endif
                            if(self%param_labels(i)%str == "CtfImage")then
                                write(ctfimagefunit,'(A)') trim(fname)
                            else if (.not.fnameselected(nDataline))then
                                write(filetabunit,'(A)') trim(fname)
                                fnameselected(nDataline) = .true.
                            else
                                write(filetabunit2,'(A)') trim(fname)
                            endif
                        else
                            !! Data to be inserted
                            tmpval = str2real(trim(adjustl(lineparts(i))))
                            DebugPrint " String converted ", trim(adjustl(lineparts(i))), " to ", tmpval
                            tmpval = tmpval * self%param_scale(i)
                            if(.not.allocated(sline))then
                                sline  = trim(self%param_labels(i)%str)//&
                                    &"="//trim(real2str(tmpval))
                            else
                                sline  = trim(adjustl(sline))//" "//trim(self%param_labels(i)%str)//&
                                    &"="//trim(real2str(tmpval))
                            endif
                        endif
                    end do
                    DebugPrint " Processed Data Line #",nDataline,":", trim(sline)
                    write(oritabunit,'(A)') trim(sline)
     
                else
                    print *, " Parsing records on line failed ", nargsOnDataline,  nargsParsed
                    print *, line
                    stop " line has insufficient elements "
                endif
     
                cycle
            endif  !! inData
     
            !! Parse the start of the STAR file
            if ( .not. (index(trim(line), "data_") == 0)) cycle
            if ( .not. (index(trim(line), "loop_") == 0))then
                inHeader=.true.
                DebugPrint "Begin STAR field lines ", line
            endif
        end do !! file loop
        if(nDataline /= self%num_data_lines)then
            HALT_NOW(" Num data lines mismatch in read_data_lines and read_header")
        endif
        !! close all the opened files
        if(is_open(filetabunit)) call fclose(filetabunit, errmsg="star_doc ; read_header filetab")
        if(is_open(filetabunit2)) call fclose(filetabunit2, errmsg="star_doc ; read_header filetab2")
        if(is_open(oritabunit))  call fclose(oritabunit,  errmsg="star_doc ; read_header oritab")
        if(is_open(imagenamefunit)) call fclose(imagenamefunit, errmsg="star_doc ; read_header imagename")
        if(is_open(ctfimagefunit)) call fclose(ctfimagefunit, errmsg="star_doc ; read_header ctfimage")
     
        !! Delete empty files
        if(nlines("filetab-stardoc2.txt")==0) call del_file("filetab-stardoc2.txt")
        if(nlines("ctfimage-stardoc.txt")==0) call del_file("ctfimage-stardoc.txt")
        if(nlines("imagename-stardoc.txt")==0) call del_file("imagename-stardoc.txt")
     
        !! rewind file
        rewind( self%funit,IOSTAT=ios)
        if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
     
        print *, ">>>  STAR IMPORT DATA INFO "
        print *, ">>>  number of data lines imported: ",self%num_data_lines
     
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
     
    subroutine print (self)
        class(stardoc), intent(inout) :: self
        integer :: i
        write(*,*) " Star formatted project information "
        write(*,*) " Star file:               ", trim(adjustl(self%current_file))
        write(*,*) " # header parameters:     ", self%num_data_elements
        write(*,*) " # compatible parameters: ", self%num_valid_elements
        write(*,*) " # data lines:            ", self%num_data_lines
     
        if(self%num_valid_elements>0 .and. allocated(self%param_labels))then
            do i=1, self%num_data_elements
                if(self%param_converted(i)==1)then
                    if(self%param_isstr(i))then
                        write(*,'(A30,1x,A30)') self%param_starlabels(i)%str,self%param_labels(i)%str 
                    else
                        write(*,'(A30,1x,A30,2x,F15.8)') self%param_starlabels(i)%str,self%param_labels(i)%str, self%param_scale(i)
                    end if
                else
                    write(*,'(A30,1x,A30)') self%param_starlabels(i)%str, "N/A" 
                end if
            end do
        end if
    end subroutine print
     
     
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
        if(allocated(self%param_starlabels))then
            do i=1,self%num_data_elements
                if(allocated(self%param_starlabels(i)%str))  deallocate(self%param_starlabels(i)%str)
            end do
            deallocate(self%param_starlabels)
        end if
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
    subroutine kill( self,keeptabs)
        class(stardoc), intent(inout) :: self
        logical, optional :: keeptabs
        logical :: keephere
        integer :: i
        keephere = .false.
        if(present(keeptabs)) keephere=keeptabs
        if( self%l_open ) call self%close
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
        self%existence =.false.
        if(.not.keephere)then
            if(file_exists('filetab-stardoc.txt')) call del_file("filetab-stardoc.txt")
            if(file_exists('oritab-stardoc.txt')) call del_file("oritab-stardoc.txt")
        endif
     
    end subroutine kill

end module simple_stardoc
