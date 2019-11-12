module simple_starfile_wrappers
    use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr, C_NULL_CHAR, C_char, C_double, C_float, C_bool, C_long, C_F_pointer
    implicit none

#include "starfile/starfile_enum.inc"

    type :: starfile_table_type
        private
        type(C_ptr), public :: object = C_NULL_ptr
    end type starfile_table_type

    interface

        subroutine C_print_pointer(this) bind(C, name="print_pointer")
            import
            type(C_ptr), value :: this
        end subroutine C_print_pointer

        function C_starfile_table__new() result(this) bind(C,name="StarFileTable__new")
            import
            type(C_ptr) :: this
        end function C_starfile_table__new

        subroutine C_starfile_table__delete(this) bind(C,name="StarFileTable__delete")
            import
            type(C_ptr), value :: this
        end subroutine C_starfile_table__delete

        subroutine C_starfile_table__addObject(this) bind(C,name="StarFileTable__addObject")
            import
            type(C_ptr), value :: this
        end subroutine C_starfile_table__addObject

        subroutine C_starfile_table__setIsList(this, is_list) bind(C,name="StarFileTable__setIsList")
            import
            type(C_ptr),                 value :: this
            logical(C_bool), intent(in), value :: is_list
        end subroutine C_starfile_table__setIsList

        subroutine C_starfile_table__setValue_double(this, EMDL_id, avalue, object_id) bind(C,name="StarFileTable__setValue_double")
            import
            type(C_ptr),                 value :: this
            integer(C_int),  intent(in), value :: EMDL_id
            real(C_double),  intent(in), value :: avalue
            integer(C_long), intent(in), value :: object_id
        end subroutine C_starfile_table__setValue_double

        subroutine C_starfile_table__setValue_float(this, EMDL_id, avalue, object_id) bind(C,name="StarFileTable__setValue_float")
            import
            type(C_ptr),                 value :: this
            integer(C_int),  intent(in), value :: EMDL_id
            real(C_float),   intent(in), value :: avalue
            integer(C_long), intent(in), value :: object_id
        end subroutine C_starfile_table__setValue_float

        subroutine C_starfile_table__setValue_int(this, EMDL_id, avalue, object_id) bind(C,name="StarFileTable__setValue_int")
            import
            type(C_ptr),                 value :: this
            integer(C_int),  intent(in), value :: EMDL_id
            integer(C_int),  intent(in), value :: avalue
            integer(C_long), intent(in), value :: object_id
        end subroutine C_starfile_table__setValue_int

        subroutine C_starfile_table__setValue_bool(this, EMDL_id, avalue, object_id) bind(C,name="StarFileTable__setValue_bool")
            import
            type(C_ptr), value                 :: this
            integer(C_int),  intent(in), value :: EMDL_id
            logical(C_bool), intent(in), value :: avalue
            integer(C_long), intent(in), value :: object_id
        end subroutine C_starfile_table__setValue_bool

        subroutine C_starfile_table__setValue_string(this, EMDL_id, avalue, object_id) bind(C,name="StarFileTable__setValue_string")
            import
            type(C_ptr),                   value :: this
            integer(C_int),    intent(in), value :: EMDL_id
            character(C_char), intent(in)        :: avalue
            integer(C_long),   intent(in), value :: object_id
        end subroutine C_starfile_table__setValue_string

        function   C_starfile_table__getValue_double(this, EMDL_id, avalue, object_id) result(aresult) bind(C,name="StarFileTable__getValue_double")
            import
            type(C_ptr),                   value :: this
            integer(C_int),    intent(in), value :: EMDL_id
            real(C_double),    intent(out)       :: avalue
            integer(C_long),   intent(in), value :: object_id
            logical(C_bool)                      :: aresult
        end function C_starfile_table__getValue_double

        function   C_starfile_table__getValue_float(this, EMDL_id, avalue, object_id) result(aresult) bind(C,name="StarFileTable__getValue_float")
            import
            type(C_ptr),                   value :: this
            integer(C_int),    intent(in), value :: EMDL_id
            real(C_float),     intent(out)       :: avalue
            integer(C_long),   intent(in), value :: object_id
            logical(C_bool)                      :: aresult
        end function C_starfile_table__getValue_float

        function   C_starfile_table__getValue_int(this, EMDL_id, avalue, object_id) result(aresult) bind(C,name="StarFileTable__getValue_int")
            import
            type(C_ptr),                   value :: this
            integer(C_int),    intent(in), value :: EMDL_id
            integer(C_int),    intent(out)       :: avalue
            integer(C_long),   intent(in), value :: object_id
            logical(C_bool)                      :: aresult
        end function C_starfile_table__getValue_int

        function   C_starfile_table__getValue_bool(this, EMDL_id, avalue, object_id) result(aresult) bind(C,name="StarFileTable__getValue_bool")
            import
            type(C_ptr),                   value :: this
            integer(C_int),    intent(in), value :: EMDL_id
            logical(C_bool),   intent(out)       :: avalue
            integer(C_long),   intent(in), value :: object_id
            logical(C_bool)                      :: aresult
        end function C_starfile_table__getValue_bool

        subroutine C_starfile_table__getValue_string(this, EMDL_id, str, object_id, alen, aresult) bind(C,name="StarFileTable__getValue_string")
            import
            type(C_ptr),                 value :: this
            integer(C_int),  intent(in), value :: EMDL_id
            type(C_ptr),     intent(out)       :: str
            integer(C_long), intent(in), value :: object_id
            integer(C_int),  intent(out)       :: alen
            logical(C_bool), intent(out)       :: aresult
        end subroutine C_starfile_table__getValue_string

        subroutine C_dealloc_str(c_str) bind(C,name="dealloc_str")
            import
            type(C_ptr), value :: c_str
        end subroutine C_dealloc_str

        subroutine C_starfile_table__clear(this) bind(C,name="StarFileTable__clear")
            import
            type(C_ptr), value :: this
        end subroutine C_starfile_table__clear

        subroutine C_starfile_table__open_ofile(this, fname, mode) bind(C,name="StarFileTable__open_ofile")
            import
            type(C_ptr),       value             :: this
            character(C_char),        intent(in) :: fname
            integer(C_int),    value, intent(in) :: mode
        end subroutine C_starfile_table__open_ofile

        subroutine C_starfile_table__write_ofile(this) bind(C,name="StarFileTable__write_ofile")
            import
            type(C_ptr), value :: this
        end subroutine C_starfile_table__write_ofile

        subroutine C_starfile_table__close_ofile(this) bind(C,name="StarFileTable__close_ofile")
            import
            type(C_ptr), value :: this
        end subroutine C_starfile_table__close_ofile

        subroutine C_starfile_table__read(this, fname, name) bind(C,name="StarFileTable__read")
            import
            type(C_ptr),       value      :: this
            character(C_char), intent(in) :: fname
            character(C_char), intent(in) :: name
        end subroutine C_starfile_table__read

        subroutine C_starfile_table__setname(this, aname) bind(C, name="StarFileTable__setName")
            import
            type(C_ptr), value :: this
            character(C_char), intent(in) :: aname
        end subroutine C_starfile_table__setname

        subroutine C_starfile_table__setcomment(this, acomment) bind(C, name="StarFileTable__setComment")
            import
            type(C_ptr), value :: this
            character(C_char), intent(in) :: acomment
        end subroutine C_starfile_table__setcomment

        function C_starfile_table__firstobject(this) result(ares) bind(C, name="StarFileTable__firstObject")
            import
            type(C_ptr), value :: this
            integer(C_long)    :: ares
        end function C_starfile_table__firstobject

        function C_starfile_table__numberofobjects(this) result(ares) bind(C, name="StarFileTable__numberOfObjects")
            import
            type(C_ptr), value :: this
            integer(C_long)    :: ares
        end function C_starfile_table__numberofobjects

        function C_starfile_table__nextobject(this) result(ares) bind(C, name="StarFileTable__nextObject")
            import
            type(C_ptr), value :: this
            integer(C_long)    :: ares
        end function C_starfile_table__nextobject

    end interface

contains

    subroutine starfile_table__new(this)
        type(starfile_table_type), intent(out) :: this
        this%object = C_starfile_table__new()
    end subroutine starfile_table__new

    subroutine starfile_table__delete(this)
        type(starfile_table_type), intent(inout) :: this
        call C_starfile_table__delete(this%object)
        this%object = C_NULL_PTR
    end subroutine starfile_table__delete

    subroutine starfile_table__addObject(this)
        type(starfile_table_type), intent(inout) :: this
        call C_starfile_table__addObject(this%object)
    end subroutine starfile_table__addObject

    subroutine starfile_table__setIsList(this, is_list)
        type(starfile_table_type), intent(inout) :: this
        logical,                   intent(in)    :: is_list
        logical(C_bool)                          :: is_list2
        is_list2 = is_list
        call C_starfile_table__setIsList(this%object, is_list2)
    end subroutine starfile_table__setIsList

    subroutine starfile_table__setValue_double(this, EMDL_id, avalue, object_id)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        real(C_double),            intent(in)    :: avalue
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        call C_starfile_table__setValue_double(this%object, EMDL_id, avalue, object_id_here)
    end subroutine starfile_table__setValue_double

    subroutine starfile_table__setValue_float(this, EMDL_id, avalue, object_id)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        real(C_float),             intent(in)    :: avalue
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        call C_starfile_table__setValue_float(this%object, EMDL_id, avalue, object_id_here)
    end subroutine starfile_table__setValue_float

    subroutine starfile_table__setValue_int(this, EMDL_id, avalue, object_id)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        integer(C_int),            intent(in)    :: avalue
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        call C_starfile_table__setValue_int(this%object, EMDL_id, avalue, object_id_here)
    end subroutine starfile_table__setValue_int

    subroutine starfile_table__setValue_bool(this, EMDL_id, avalue, object_id)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        logical,                   intent(in)    :: avalue
        logical(C_bool)                          :: avalue2
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        avalue2 = avalue
        call C_starfile_table__setValue_bool(this%object, EMDL_id, avalue2, object_id_here)
    end subroutine starfile_table__setValue_bool

    subroutine starfile_table__setValue_string(this, EMDL_id, avalue, object_id)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        character(len=*),          intent(in)    :: avalue
        character(len=:), allocatable :: str
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        str = avalue // C_NULL_CHAR
        call C_starfile_table__setValue_string(this%object, EMDL_id, str, object_id_here)
    end subroutine starfile_table__setValue_string

    function starfile_table__getValue_double(this, EMDL_id, avalue, object_id) result(aresult)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        real(C_double),            intent(out)   :: avalue
        logical(C_bool) :: aresult
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        aresult = C_starfile_table__getValue_double(this%object, EMDL_id, avalue, object_id_here)
    end function starfile_table__getValue_double

    function starfile_table__getValue_float(this, EMDL_id, avalue, object_id) result(aresult)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        real(C_float),             intent(out)   :: avalue
        logical(C_bool) :: aresult
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        aresult = C_starfile_table__getValue_float(this%object, EMDL_id, avalue, object_id_here)
    end function starfile_table__getValue_float

    function starfile_table__getValue_int(this, EMDL_id, avalue, object_id) result(aresult)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        integer(C_int),            intent(out)   :: avalue
        logical(C_bool) :: aresult
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        aresult = C_starfile_table__getValue_int(this%object, EMDL_id, avalue, object_id_here)
    end function starfile_table__getValue_int

    function starfile_table__getValue_bool(this, EMDL_id, avalue, object_id) result(aresult)
        type(starfile_table_type), intent(inout) :: this
        integer(C_int),            intent(in)    :: EMDL_id
        logical(C_bool),           intent(out)   :: avalue
        logical(C_bool) :: aresult
        integer(C_long), optional, intent(in)    :: object_id
        integer(C_long) :: object_id_here
        if( present(object_id) ) then
            object_id_here = object_id
        else
            object_id_here = -1
        end if
        aresult = C_starfile_table__getValue_bool(this%object, EMDL_id, avalue, object_id_here)
    end function starfile_table__getValue_bool

    function starfile_table__getValue_string(this, EMDL_id, str, object_id) result(aresult)
      type(starfile_table_type),     intent(inout) :: this
      integer(C_int),                intent(in)    :: EMDL_id
      character(len=:), allocatable, intent(out)   :: str
      logical                                      :: aresult
      type(C_ptr)                                         :: c_str
      character(len=1,kind=C_char), dimension(:), pointer :: f_str
      integer(c_int)                                      :: alen
      integer                                             :: i
      logical(C_bool)                                     :: aresult2
      integer(C_long), optional, intent(in)    :: object_id
      integer(C_long) :: object_id_here
      if( present(object_id) ) then
          object_id_here = object_id
      else
          object_id_here = -1
      end if
      call C_starfile_table__getValue_string(this%object, EMDL_id, c_str, &
                                             object_id_here, alen, aresult2)              ! first call C routine
      if (aresult2) then
         allocate( character(len=alen) :: str )                                           ! allocate memory for the fortran string
         call c_f_pointer(c_str, f_str, [alen])                                           ! cast c-pointer to fortran pointer
         do i = 1, alen
            str(i:i) = f_str(i)                                                           ! copy over the data
         end do
         call C_dealloc_str(c_str)                                                        ! free c-pointer to avoid memory leak
      end if
      aresult = aresult2
    end function starfile_table__getValue_string

    subroutine starfile_table__clear(this)
        type(starfile_table_type), intent(inout) :: this
        call C_starfile_table__clear(this%object)
    end subroutine starfile_table__clear

    subroutine starfile_table__open_ofile(this, fname, mode_)
        type(starfile_table_type), intent(inout) :: this
        character(len=*),          intent(in)    :: fname
        character(len=:),          allocatable   :: fname2
        integer,         optional, intent(in)    :: mode_    ! 0: NEW, 1: APPEND
        integer :: mode
        if (present(mode_)) then
            mode = mode_
        else
            mode = 0
        end if
        fname2 = fname // C_NULL_CHAR
        call C_starfile_table__open_ofile(this%object, fname2, mode)
    end subroutine starfile_table__open_ofile

    subroutine starfile_table__write_ofile(this)
        type(starfile_table_type), intent(inout) :: this
        call C_starfile_table__write_ofile(this%object)
    end subroutine starfile_table__write_ofile

    subroutine starfile_table__close_ofile(this)
        type(starfile_table_type), intent(inout) :: this
        call C_starfile_table__close_ofile(this%object)
    end subroutine starfile_table__close_ofile

    subroutine starfile_table__read(this, fname, name)
        type(starfile_table_type),  intent(inout) :: this
        character(len=*),           intent(in)    :: fname
        character(len=*), optional, intent(in)    :: name
        character(len=:), allocatable :: str_fname
        character(len=:), allocatable :: str_name
        str_fname = fname // C_NULL_CHAR
        if( present(name) ) then
            str_name = name // C_NULL_CHAR
        else
            str_name = C_NULL_CHAR
        end if
        call C_starfile_table__read(this%object, str_fname, str_name)
    end subroutine starfile_table__read

    subroutine starfile_table__setname(this, aname)
        type(starfile_table_type), intent(inout) :: this
        character(len=*),          intent(in)    :: aname
        character(len=:),          allocatable   :: aname2
        aname2 = aname // C_NULL_CHAR
        call C_starfile_table__setname(this%object, aname2)
    end subroutine starfile_table__setname

    subroutine starfile_table__setcomment(this, acomment)
        type(starfile_table_type), intent(inout) :: this
        character(len=*),          intent(in)    :: acomment
        character(len=:),          allocatable   :: acomment2
        acomment2 = acomment // C_NULL_CHAR
        call C_starfile_table__setcomment(this%object, acomment2)
    end subroutine starfile_table__setcomment

    function starfile_table__firstobject(this) result(ares)
        type(starfile_table_type), intent(inout) :: this
        integer(C_long) :: ares
        ares = C_starfile_table__firstobject(this%object)
    end function starfile_table__firstobject

    function starfile_table__numberofobjects(this) result(ares)
        type(starfile_table_type), intent(inout) :: this
        integer(C_long) :: ares
        ares = C_starfile_table__numberofobjects(this%object)
    end function starfile_table__numberofobjects

    function starfile_table__nextobject(this) result(ares)
        type(starfile_table_type), intent(inout) :: this
        integer(C_long) :: ares
        ares = C_starfile_table__nextobject(this%object)
    end function starfile_table__nextobject

end module simple_starfile_wrappers
