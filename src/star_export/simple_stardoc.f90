
module simple_stardoc
    include 'simple_lib.f08'
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
        character(len=KEYLEN), allocatable :: param_labels(:)
        type(str4arr), allocatable         :: data(:)
        type(starframes), allocatable      :: frames(:)
        integer,public   :: num_data_elements    = 0
        integer,public   :: num_data_lines       = 0
        integer          :: num_frames           = 0  
        integer          :: num_frame_elements   = 0
        integer          :: funit
        logical          :: l_open               =.false.
        logical          :: existence            =.false. !< to indicate existence
    contains
        procedure        :: new
        procedure,public :: open
        procedure,public :: close
        procedure,public :: read_header
procedure,public         :: read_data_labels
        procedure        :: write
        procedure        :: get_r4
        procedure        :: get_i4
        procedure        :: get_str
        !generic         :: get => get_r4, get_i4, get_str
        procedure        :: put_r4
        procedure        :: put_i4
        procedure        :: put_str
        !generic         :: put => put_r4, put_i4, put_str
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
        integer,     intent(in)    :: n
        integer :: i
        call self%kill_doc
        self%num_frame_elements = n
        allocate( self%frames(n), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('new; simple_stardoc',alloc_stat)
        do i=1,n
            call self%frames(i)%new_frame
        end do
        DebugPrint ' simple_oris::new initiated'
    end subroutine new

    subroutine open(self, filename) 
        class(stardoc), intent(inout) :: self
        character(len=*),intent(inout) :: filename
        integer :: io_stat, tmpunit,filesz
        if(self%l_open) call self%close
        if(allocated(self%param_labels)) deallocate(self%param_labels)
        if(.not. file_exists(trim(filename) ))&
            call simple_stop("simple_stardoc::open ERROR file does not exist "//trim(filename) )
        call fopen(tmpunit, trim(filename), access='STREAM', action='READWRITE',&
            &status='UNKNOWN', form='UNFORMATTED', iostat=io_stat)
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
        call self%read_header
        if(self%num_data_elements > 0 )then
            allocate(self%param_labels(self%num_data_elements))
            call self%read_data_labels
        end if
    end subroutine open


    subroutine close( self )
        class(stardoc), intent(inout) :: self !< instance
        integer :: io_stat
        if( self%l_open )then
            call fclose(self%funit,io_stat,errmsg='stardoc ; close ')
            self%l_open = .false.
        end if
    end subroutine close
    subroutine read_header(self) 
        class(stardoc), intent(inout) :: self
        integer          :: n,ios,lenstr
        character(len=LINE_MAX_LEN) :: line ! 8192
        logical :: inData, inField
        self%num_data_elements=0
        self%num_data_lines=0
        inData=.false.;inField=.false.
        if(self%l_open)then
            do
                read(self%funit,*,IOSTAT=ios) line
                if(ios /= 0) exit

                !! Count number of fields in header
                if(inField)then
                    if (line(1:4) == "_rln" )then
                        self%num_data_elements=self%num_data_elements+1
                        DebugPrint " Found STAR field line ", line
                        cycle
                    else
                        inField=.false.
                        DebugPrint " End of STAR data field lines "
                        inData = .true.
                    endif
                endif
                !! Count number of data lines
                if(inData)then
                    if (line == "" )then
                        inData=.false.
                        DebugPrint " End of STAR data lines ", line
                        exit
                    endif
                    self%num_data_lines=self%num_data_lines+1
                    DebugPrint " Found STAR data line ", line
                    cycle
                end if
                !! Parse the start of the STAR file
                lenstr=len_trim(line)
                if (lenstr == 0 )cycle ! empty line
                if ( line(1:5) == "data_")then !! e.g. data_
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
                if (line == "loop_" )then
                    inField=.true.
                    DebugPrint "Begin STAR field lines ", line
                endif

            end do
            rewind( self%funit,IOSTAT=ios)
            if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
            DebugPrint " STAR field lines: ", self%num_data_elements
            DebugPrint " STAR data lines : ", self%num_data_lines
        endif

    end subroutine read_header
    subroutine read_data_labels(self)
        class(stardoc), intent(inout) :: self
        integer          :: n,ios,lenstr, pos
        character(len=LINE_MAX_LEN) :: line ! 8192
        logical :: inData, inField
        inField=.false.
        n=1
        if(self%l_open)then
            do
                read(self%funit,*,IOSTAT=ios) line
                if(ios /= 0) exit

                !! Count number of fields in header
                if(inField)then
                    if (lenstr == 0 )cycle ! in case there is an empty line after 'loop_'
                    if (trim(line(1:4)) == "_rln" )then
                        pos = firstNonBlank(line)
                        if(pos <= 5)then
                            write(*,*) 'line: ', line
                            stop 'field too small to contain a header; stardoc :: read_data_labels'
                        end if
                        !! shorten label
                        if(pos >= KEYLEN+5)pos = KEYLEN+4
                        self%param_labels(n) = trim(line(5:pos))
                        DebugPrint " STAR param field : ",n, " : ", self%param_labels(n)
                        n = n+1
                    else
                        inField=.false.
                        DebugPrint " End of STAR data field lines "
                        exit
                    endif
                endif

                !! Parse the start of the STAR file
                lenstr=len_trim(line)
                if (lenstr == 0 )cycle ! empty line
                if ( line(1:4) == "data") cycle
                if (trim(line) == "loop_" )then
                    inField=.true.
                    DebugPrint "Begin STAR field lines ", line
                endif
            end do
            rewind( self%funit,IOSTAT=ios)
            if(ios/=0)call fileiochk('star_doc ; read_header - rewind failed ', ios)
        endif
    end subroutine read_data_labels

    subroutine  write(self, sp, vars, filename)
        use simple_sp_project
        class(stardoc), intent(inout) :: self
        class(sp_project), intent(inout) :: sp
        character(len=KEYLEN) :: vars(:)
        character(len=*),intent(inout) :: filename
        character(len=LONGSTRLEN*8) :: starline
        character(len=LONGSTRLEN) :: imagename
        character(len=:), allocatable ::  val,stackfile,state
        integer  :: i, relionstar, io_stat
        real     :: statef, defocus
        if(self%l_open)call self%close
        call fopen(relionstar, trim(filename)) !,stat=io_stat)
        !if(io_stat/=0) call fileiochk("In stardoc; write; unable to open "//trim(filename))
        write(relionstar,'(A)') ""
        write(relionstar,'(A)') "data_"
        write(relionstar,'(A)') ""
        write(relionstar,'(A)') "loop_"
        do i=1, self%num_data_elements
            write(relionstar,'("_rln",A,3x,"#",I0)') self%data(i)%str, i
        end do

        do i=1, self%num_frames
            statef = self%frames(i)%state
            if (statef .ne. 0)then
                !! create zero-padded frame number and stackfile
                imagename=int2str_pad(self%get_i4(i,"frameid"),5)//'@'//&
                    self%get_str(i,"stackfile")//'s'

                stackfile = self%get_str(i,"stackfile")
                call syslib_symlink(trim(stackfile), trim(stackfile)//'s',status=io_stat)
                if(io_stat/=0)call simple_stop("in write; StarDoc")

                starline =imagename//&
                    &" "//real2str(self%get_r4(i,"kv"))//&
                    &" "//real2str(self%get_r4(i,"dfx")*1000)//&
                    &" "//real2str(self%get_r4(i,"dfy")*1000)//&
                    &" "//real2str(self%get_r4(i,"angast"))//&
                    &" "//real2str(self%get_r4(i,"cs"))//&
                    &" "//real2str(self%get_r4(i,"fraca"))//&
                    &" 10000"//&
                    &" "//real2str(self%get_r4(i,"smpd"))
                write(relionstar,'(A)') starline
            end if
        end do
        call fclose(relionstar)
    end subroutine write

    subroutine putdata(self, strarr)
        class(stardoc), intent(inout) :: self
        character(len=*), intent(in) :: strarr(:)
        integer :: i
        self%num_data_elements = size(strarr)
        allocate(self%data(self%num_data_elements))
        do i=1, self%num_data_elements
            allocate(self%data(i)%str, source=strarr(i))
        end do
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
    !>  \brief  is a destructor
    subroutine kill_doc( self )
        class(stardoc), intent(inout) :: self
        integer :: i
        if( allocated(self%frames) )then
            do i=1,self%num_frames
                call self%frames(i)%kill_frame
            end do
            deallocate( self%frames , stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: kill, module: simple_stardoc ')
            self%num_frame_elements = 0
        endif
        if(self%l_open) call self%close
        if(allocated(self%param_labels)) deallocate(self%param_labels)
        self%existence =.false.
    end subroutine kill_doc

end module simple_stardoc
