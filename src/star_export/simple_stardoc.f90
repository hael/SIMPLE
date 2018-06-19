
module simple_stardoc
include 'simple_lib.f08'
implicit none

type starframes
    type(hash)  :: htab              !< hash table for the parameters
    type(chash) :: chtab
    integer     :: state
    logical     :: existence=.false. !< to indicate existence
contains
    procedure :: new_frame
    procedure :: kill_frame
end type starframes

interface starframes
    module procedure constructor_frames
end interface starframes

type stardoc
    type(str4arr), allocatable :: data(:)
    type(starframes), allocatable :: frames(:)
    integer     :: num_data_elements
    integer     :: num_frames
    integer     :: num_frame_elements
    logical     :: existence=.false. !< to indicate existence
contains
    procedure :: new
    procedure :: write
    procedure :: get_r4
    procedure :: get_i4
    procedure :: get_str
    !generic :: get => get_r4, get_i4, get_str
    procedure :: put_r4
    procedure :: put_i4
    procedure :: put_str
    !generic   :: put => put_r4, put_i4, put_str
    procedure  :: kill
end type stardoc

interface stardoc
    module procedure constructor
end interface
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
        call self%kill
        self%num_frame_elements = n
        allocate( self%frames(n), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('new; simple_stardoc',alloc_stat)
        do i=1,n
            call self%frames(i)%new_frame
        end do
        DebugPrint ' simple_oris::new initiated'
    end subroutine new

    subroutine  write(self,  filename)
        class(stardoc), intent(inout) :: self
        character(len=*),intent(inout) :: filename
        character(len=LONGSTRLEN*8) :: starline
        character(len=LONGSTRLEN) :: imagename
        character(len=:), allocatable ::  val,stackfile,state
        integer  :: i, relionstar, io_stat
        real     :: statef, defocus

        call fopen(relionstar, trim(filename)) !,stat=io_stat)
        !if(io_stat/=0) call fileiochk("In stardoc; write; unable to open "//trim(filename))
        write(relionstar,'(A)') ""
        write(relionstar,'(A)') "data_"
        write(relionstar,'(A)') ""
        write(relionstar,'(A)') "loop_"
        do i=1, self%num_data_elements
            write(relionstar,'("_rln",A,3x,"#",I0)') self%data(i)%str, i
        end do

        ! write(relionstar,'(A)') "_rlnImageName"
        ! write(relionstar,'(A)') "_rlnVoltage"
        ! write(relionstar,'(A)') "_rlnDefocusU"
        ! write(relionstar,'(A)') "_rlnDefocusV"
        ! write(relionstar,'(A)') "_rlnDefocusAngle"
        ! write(relionstar,'(A)') "_rlnSphericalAberration"
        ! write(relionstar,'(A)') "_rlnAmplitudeContrast"
        ! write(relionstar,'(A)') "_rlnMagnification"
        ! write(relionstar,'(A)') "_rlnDetectorPixelSize"

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
    subroutine kill( self )
        class(stardoc), intent(inout) :: self
        integer :: i
        if( allocated(self%frames) )then
            do i=1,self%num_frames
                call self%frames(i)%kill_frame
            end do
            deallocate( self%frames , stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: kill, module: simple_oris o')
            self%num_frame_elements = 0
        endif
    end subroutine kill

  end module simple_stardoc
