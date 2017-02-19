!>  \brief  image header records are contained in image headers. Each contains a piece of information regarding the image file.
module simple_imgheadrec
use simple_defs
implicit none

public :: int_imgheadrec, real_imgheadrec, char_imgheadrec
private

!>  a header record is a labelled value which describes the characteristics of an imagefile
type :: imgheadrec
    private
    integer                  :: index_position          !<  the position of the record within the file header. starting at 1 and incrementing
    integer                  :: byte_position           !<  the position of the first byte of the record within the header.
    integer(kind=1), pointer :: byte_array(:) => null() !<  pointer to the array of bytes containing the actual header values
contains
    procedure ::  new
    procedure ::  kill
end type

type, extends(imgheadrec) :: int_imgheadrec
contains
    procedure, private ::  getIntg
    procedure          ::  get => getIntg
    procedure, private ::  setIntg
    generic            ::  assignment(=) => setIntg
end type

type, extends(imgheadrec) :: real_imgheadrec
contains
    procedure, private :: getReal
    procedure          :: get => getReal
    procedure, private :: SetReal
    generic            :: assignment(=) => SetReal
end type

type, extends(imgheadrec) :: char_imgheadrec
integer :: length=0 !<  length of the string of characters
contains
    procedure, private :: getChar
    procedure          :: get => getChar
    procedure, private :: SetChar
    generic            :: assignment(=) => SetChar
end type

interface imgrec
    module procedure constructor
end interface

interface assignment(=)
    module procedure  getIntgAssign
    module procedure  getRealAssign
    module procedure  getCharAssign
end interface

contains
    
    ! CONSTRUCTORS
    
    !>  \brief constructor
    function constructor( index_position, byte_position, byte_array, length ) result( self )
        integer, intent(in)                 :: index_position !<  the position of the record within the file header. starting at 1 and incrementing
        integer, intent(in)                 :: byte_position  !<  the position of the first byte of the record within the header.
        integer(kind=1), target, intent(in) :: byte_array(:)  !<  byte array to point to
        integer, optional, intent(in)       :: length         !<  length of character string
        type(imgheadrec)                    :: self
        call self%new(index_position, byte_position, byte_array, length)
    end function

    !>  \brief constructor
    subroutine new( self, index_position, byte_position, byte_array, length )
        class(imgheadrec), intent(inout)    :: self            !<  header record
        integer, intent(in)                 ::  index_position !<  the position of the record within the file header. starting at 1 and incrementing
        integer, intent(in)                 ::  byte_position  !<  the position of the first byte of the record within the header.
        integer(kind=1), target, intent(in) ::  byte_array(:)  !<  byte array to point to
        integer, optional, intent(in)       ::  length         !<  length of character string
        self%index_position =  index_position
        self%byte_position  =  byte_position
        self%byte_array     => byte_array
        select type(self)
            type is (char_imgheadrec)
                self%length = 1
                if (present(length)) self%length = length
        end select
    end subroutine
    
    ! SETTERS
    
    subroutine setIntg(self,value)
        use, intrinsic :: iso_c_binding
        class(int_imgheadrec), intent(inout) :: self    !<  header record
        integer(kind=4), target, intent(in)  :: value   !<  value to be written to header
        integer(kind=1), pointer             :: byte_ptr(:)
        type(c_ptr)                          :: cptr
        integer                              :: ptr_shape(1)
        cptr = c_loc(value)
        ptr_shape = [4]
        call c_f_pointer(cptr,byte_ptr,ptr_shape)
        self%byte_array(self%byte_position:self%byte_position+3) = byte_ptr(1:4)
        nullify(byte_ptr)
    end subroutine

    subroutine SetReal(self,value)
        use, intrinsic :: iso_c_binding
        class(real_imgheadrec), intent(inout) :: self    !<  header record
        real(kind=4), target, intent(in)      :: value   !<  value to be written to header
        integer(kind=1), pointer              :: byte_ptr(:)
        type(c_ptr)                           :: cptr
        integer                               :: ptr_shape(1)
        cptr = c_loc(value)
        ptr_shape = [4]
        call c_f_pointer(cptr,byte_ptr,ptr_shape)
        self%byte_array(self%byte_position:self%byte_position+3) = byte_ptr(1:4)
        nullify(byte_ptr)
    end subroutine

    subroutine SetChar(self,value)
        use simple_strings, only: int2str
        class(char_imgheadrec), intent(inout) :: self !<  header record
        character(len=*), intent(in)          :: value   !<  value to be written to header
        character(len=:), allocatable         :: tmp_string
        allocate(character(len=self%length) :: tmp_string)
        write(tmp_string,'(a'//int2str(self%length)//')') value
        self%byte_array(self%byte_position:self%byte_position+self%length-1) = transfer(tmp_string,self%byte_array)            
    end subroutine

    !>  \brief  read value from header record (to be used on RHS of assignments)
    subroutine getIntgAssign(value,self)
        class(int_imgheadrec), intent(in) :: self
        integer(kind=4), intent(out)      :: value
        value =  self%getIntg()
    end subroutine

    !>  \brief  read value from header record (to be used on RHS of assignments)
    subroutine getRealAssign(value,self)
        class(real_imgheadrec), intent(in) :: self
        real(kind=4), intent(out)          :: value
        value =  self%getReal()
    end subroutine

    !>  \brief  Read value from header record (to be used on RHS of assignments)
    subroutine getCharAssign(value,self)
        class(char_imgheadrec), intent(in) :: self
        character(len=:), allocatable,intent(inout) :: value
        if (allocated(value)) deallocate(value)
        value = self%getChar()
    end subroutine

    !>  \brief  Read value from header record
    integer function getIntg(self)
        use, intrinsic :: iso_c_binding
        class(int_imgheadrec), intent(in) :: self
        integer(kind=1), pointer          :: byte_ptr(:)
        type(c_ptr)                       :: cptr
        integer(kind=4),    target        :: tmpval
        integer                           :: ptr_shape(1)
        cptr = c_loc(tmpval)
        ptr_shape = [4]
        call c_f_pointer(cptr,byte_ptr,ptr_shape)
        byte_ptr(1:4) = self%byte_array(self%byte_position:self%byte_position+3)
        nullify(byte_ptr)
        getIntg = tmpval
    end function

    !>  \brief  Read value from header record
    real function getReal(self)
        use, intrinsic :: iso_c_binding
        class(real_imgheadrec), intent(in) :: self
        integer(kind=1), pointer           :: byte_ptr(:)
        type(c_ptr)                        :: cptr
        real(kind=4), target               :: tmpval
        integer                            :: ptr_shape(1)
        cptr = c_loc(tmpval)
        ptr_shape = [4]
        call c_f_pointer(cptr,byte_ptr,ptr_shape)
        byte_ptr(1:4) = self%byte_array(self%byte_position:self%byte_position+3)
        nullify(byte_ptr)
        getReal = tmpval
    end function

    !>  \brief  Read value from header record
    function getChar(self)
        use, intrinsic :: iso_c_binding
        class(char_imgheadrec), intent(in)  :: self
        character(len=:), allocatable       :: getChar
        allocate(character(len=self%length) :: getChar)
        getChar = transfer(self%byte_array(self%byte_position:self%byte_position+self%length-1),getChar)
    end function
    
    !>  \brief destructor
    subroutine kill(self)
        class(imgheadrec), intent(inout) :: self
        if(associated(self%byte_array)) nullify(self%byte_array)
    end subroutine
    
end module simple_imgheadrec
