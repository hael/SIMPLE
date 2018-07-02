module static_intrinsic_type_m
   implicit none
   type :: static_intrinsic_type_t
      integer :: x
   contains
      generic :: assignment(=) => assign_static_intrinsic_type
      generic :: operator(+) => add_static_intrinsic_type
      procedure, pass(lhs) :: assign_static_intrinsic_type
      procedure, pass(lhs) :: add_static_intrinsic_type
   endtype static_intrinsic_type_t
contains
   subroutine assign_static_intrinsic_type(lhs, rhs)
      ! Operator `=`.
      class(static_intrinsic_type_t), intent(inout) :: lhs
      class(static_intrinsic_type_t), intent(in)    :: rhs
      lhs%x = rhs%x
   endsubroutine assign_static_intrinsic_type

   function add_static_intrinsic_type(lhs, rhs) result(res)
      ! Operator `+`.
      class(static_intrinsic_type_t), intent(in)  :: lhs
      class(static_intrinsic_type_t), intent(in)  :: rhs
      class(static_intrinsic_type_t), allocatable :: res
      allocate (static_intrinsic_type_t :: res)
      res%x = lhs%x + rhs%x
   endfunction add_static_intrinsic_type
endmodule static_intrinsic_type_m


module static_dynamic_intrinsic_type_m
   implicit none
   type :: static_dynamic_intrinsic_type_t
      integer              :: x
      integer, allocatable :: y
   contains
      generic :: assignment(=) => assign_static_dynamic_intrinsic_type
      generic :: operator(+) => add_static_dynamic_intrinsic_type
      procedure, pass(lhs) :: assign_static_dynamic_intrinsic_type
      procedure, pass(lhs) :: add_static_dynamic_intrinsic_type
   endtype static_dynamic_intrinsic_type_t
contains
   subroutine assign_static_dynamic_intrinsic_type(lhs, rhs)
      ! Operator `=`.
      class(static_dynamic_intrinsic_type_t), intent(inout) :: lhs
      class(static_dynamic_intrinsic_type_t), intent(in)    :: rhs
      lhs%x = rhs%x
      if (allocated(rhs%y)) then
         if (.not.allocated(lhs%y)) allocate(lhs%y)
         lhs%y = rhs%y
      endif
   endsubroutine assign_static_dynamic_intrinsic_type

   function add_static_dynamic_intrinsic_type(lhs, rhs) result(res)
      ! Operator `+`.
      class(static_dynamic_intrinsic_type_t), intent(in)  :: lhs
      class(static_dynamic_intrinsic_type_t), intent(in)  :: rhs
      class(static_dynamic_intrinsic_type_t), allocatable :: res
      allocate (static_dynamic_intrinsic_type_t :: res)
      res%x = lhs%x + rhs%x
      if (allocated(lhs%y).and.allocated(rhs%y)) then
         allocate(res%y)
         res%y = lhs%y + rhs%y
      endif
   endfunction add_static_dynamic_intrinsic_type
endmodule static_dynamic_intrinsic_type_m

module dynamic_intrinsic_type_m
   implicit none
   type :: dynamic_intrinsic_type_t
      integer, allocatable :: x
   contains
      generic :: assignment(=) => assign_dynamic_intrinsic_type
      generic :: operator(+) => add_dynamic_intrinsic_type
      procedure, pass(lhs) :: assign_dynamic_intrinsic_type
      procedure, pass(lhs) :: add_dynamic_intrinsic_type
   endtype dynamic_intrinsic_type_t
contains
   subroutine assign_dynamic_intrinsic_type(lhs, rhs)
      ! Operator `=`.
      class(dynamic_intrinsic_type_t), intent(inout) :: lhs
      class(dynamic_intrinsic_type_t), intent(in)    :: rhs
      if (allocated(rhs%x)) then
        if (.not.allocated(lhs%x)) allocate(lhs%x)
        lhs%x = rhs%x
      endif
   endsubroutine assign_dynamic_intrinsic_type

   function add_dynamic_intrinsic_type(lhs, rhs) result(res)
      ! Operator `+`.
      class(dynamic_intrinsic_type_t), intent(in)  :: lhs
      class(dynamic_intrinsic_type_t), intent(in)  :: rhs
      class(dynamic_intrinsic_type_t), allocatable :: res
      allocate (dynamic_intrinsic_type_t :: res)
      if (allocated(lhs%x).and.allocated(rhs%x)) then
        allocate(res%x)
        res%x = lhs%x + rhs%x
      endif
   endfunction add_dynamic_intrinsic_type
endmodule dynamic_intrinsic_type_m


module polymorphic_type_m
   implicit none
   type :: derived_t
      integer :: x
   endtype derived_t
   type :: polymorphic_type_t
      class(derived_t), allocatable :: x
   contains
      generic :: assignment(=) => assign_polymorphic_type
      generic :: operator(+) => add_polymorphic_type
      procedure, pass(lhs) :: assign_polymorphic_type
      procedure, pass(lhs) :: add_polymorphic_type
   endtype polymorphic_type_t
contains
   subroutine assign_polymorphic_type(lhs, rhs)
      ! Operator `=`.
      class(polymorphic_type_t), intent(inout) :: lhs
      class(polymorphic_type_t), intent(in)    :: rhs
      if (allocated(rhs%x)) then
         if (.not.allocated(lhs%x)) allocate(derived_t :: lhs%x)
         select type(lhsx => lhs%x)
         type is(derived_t)
            select type(rhsx => rhs%x)
            type is(derived_t)
               lhsx%x = rhsx%x
            endselect
         endselect
      endif
   endsubroutine assign_polymorphic_type

   function add_polymorphic_type(lhs, rhs) result(res)
      ! Operator `+`.
      class(polymorphic_type_t), intent(in)  :: lhs
      class(polymorphic_type_t), intent(in)  :: rhs
      class(polymorphic_type_t), allocatable :: res
      allocate (polymorphic_type_t :: res)
      if (allocated(lhs%x).and.allocated(rhs%x)) then
         allocate(derived_t :: res%x)
         select type(lhsx => lhs%x)
         type is(derived_t)
            select type(rhsx => rhs%x)
            type is(derived_t)
               select type(resx => res%x)
               type is(derived_t)
                  resx%x = lhsx%x + rhsx%x
               endselect
            endselect
         endselect
      endif
   endfunction add_polymorphic_type
endmodule polymorphic_type_m

module inherit_type_m
   implicit none
   type :: derived_t
      integer :: x
   contains
      generic :: assignment(=) => assign_derived_type
      generic :: operator(+) => add_derived_type
      procedure, pass(lhs) :: assign_derived_type
      procedure, pass(lhs) :: add_derived_type
   endtype derived_t

   type :: inherit_type_t
      type(derived_t), allocatable :: x
   contains
      generic :: assignment(=) => assign_inherit_type
      generic :: operator(+) => add_inherit_type
      procedure, pass(lhs) :: assign_inherit_type
      procedure, pass(lhs) :: add_inherit_type
   endtype inherit_type_t
contains
   subroutine assign_derived_type(lhs, rhs)
      ! Operator `=`.
      class(derived_t), intent(inout) :: lhs
      class(derived_t), intent(in)    :: rhs
      lhs%x = rhs%x
   endsubroutine assign_derived_type

   function add_derived_type(lhs, rhs) result(res)
      ! Operator `+`.
      class(derived_t), intent(in)  :: lhs
      class(derived_t), intent(in)  :: rhs
      class(derived_t), allocatable :: res
      allocate (derived_t :: res)
      res%x = lhs%x + rhs%x
   endfunction add_derived_type

   subroutine assign_inherit_type(lhs, rhs)
      ! Operator `=`.
      class(inherit_type_t), intent(inout) :: lhs
      class(inherit_type_t), intent(in)    :: rhs
      if (allocated(rhs%x)) then
        if (.not.allocated(lhs%x)) allocate(lhs%x)
        lhs%x = rhs%x
      endif
   endsubroutine assign_inherit_type

   function add_inherit_type(lhs, rhs) result(res)
      ! Operator `+`.
      class(inherit_type_t), intent(in)  :: lhs
      class(inherit_type_t), intent(in)  :: rhs
      class(inherit_type_t), allocatable :: res
      allocate (inherit_type_t :: res)
      if (allocated(lhs%x).and.allocated(rhs%x)) then
        allocate(res%x)
        res%x = lhs%x + rhs%x
      endif
   endfunction add_inherit_type
endmodule inherit_type_m

module ancestor_m
implicit none

type, abstract :: ancestor_t
   integer, allocatable :: workaround
   contains
      generic :: assignment(=) => assign_
      generic :: operator(+) => add_
      procedure(assign_interface), pass(lhs), deferred :: assign_
      procedure(add_interface),    pass(lhs), deferred :: add_
endtype ancestor_t

abstract interface
   subroutine assign_interface(lhs, rhs)
   ! Operator `=`.
   import :: ancestor_t
   class(ancestor_t), intent(inout) :: lhs
   class(ancestor_t), intent(in)    :: rhs
   endsubroutine assign_interface

   function add_interface(lhs, rhs) result(opr)
   ! Operator `+`.
   import :: ancestor_t
   class(ancestor_t), intent(in)  :: lhs
   class(ancestor_t), intent(in)  :: rhs
   class(ancestor_t), allocatable :: opr
   endfunction add_interface
endinterface
endmodule ancestor_m

module parent_m
use ancestor_m
implicit none

type, extends(ancestor_t), abstract :: parent_t
endtype parent_t
endmodule parent_m

module child_m
use ancestor_m
use parent_m
implicit none

type, extends(parent_t) :: child_t
   integer              :: x
   integer, allocatable :: y
   contains
      procedure, pass(lhs) :: assign_
      procedure, pass(lhs) :: add_
endtype child_t

contains
   subroutine assign_(lhs, rhs)
   ! Operator `=`.
   class(child_t),    intent(inout) :: lhs
   class(ancestor_t), intent(in)    :: rhs

   select type(rhs)
   class is(child_t)
      lhs%x = rhs%x
      if (allocated(rhs%y)) then
        if (.not.allocated(lhs%y)) allocate(lhs%y)
        lhs%y = rhs%y
      endif
   endselect
   endsubroutine assign_

   function add_(lhs, rhs) result(opr)
   ! Operator `+`.
   class(child_t),    intent(in)  :: lhs
   class(ancestor_t), intent(in)  :: rhs
   class(ancestor_t), allocatable :: opr

   allocate(child_t :: opr)
   select type(opr)
   class is(child_t)
      select type(rhs)
      class is(child_t)
         opr%x = lhs%x + rhs%x
         if (allocated(lhs%y).and.allocated(rhs%y)) then
           allocate(opr%y)
           opr%y = lhs%y + rhs%y
         endif
      endselect
   endselect
   endfunction add_
endmodule child_m

program simple_test_memleaks
implicit none

call static_intrinsic_leaks_raiser
 call static_dynamic_intrinsic_leaks_raiser
    call dynamic_intrinsic_leaks_raiser
    ! call polymorphic_leaks_raiser
    ! call inherit_leaks_raiser
    ! call multiple_inheritance_leaks_raiser
    contains

    subroutine  static_intrinsic_leaks_raiser
        use static_intrinsic_type_m
        implicit none
        type(static_intrinsic_type_t) :: a
        type(static_intrinsic_type_t) :: b

        a%x = 1
        b = a
        b = a + b ! here the `+` operator could generate memory leaks
    end  subroutine static_intrinsic_leaks_raiser

subroutine static_dynamic_intrinsic_leaks_raiser
    use static_dynamic_intrinsic_type_m
    implicit none
    type(static_dynamic_intrinsic_type_t) :: a
    type(static_dynamic_intrinsic_type_t) :: b

    a%x = 1
    allocate(a%y)
    a%y = 1
    b = a
    b = a + b ! here the `+` operator could generate memory leaks
end subroutine static_dynamic_intrinsic_leaks_raiser

subroutine dynamic_intrinsic_leaks_raiser
    use dynamic_intrinsic_type_m
    implicit none
    type(dynamic_intrinsic_type_t) :: a
    type(dynamic_intrinsic_type_t) :: b

    allocate(a%x)
    a%x = 1
    b = a
    b = a + b ! here the `+` operator could generate memory leaks
end subroutine dynamic_intrinsic_leaks_raiser

subroutine polymorphic_leaks_raiser
   use polymorphic_type_m
   implicit none
   type(polymorphic_type_t) :: a
   type(polymorphic_type_t) :: b

   allocate(derived_t :: a%x)
   select type(ax => a%x)
   type is(derived_t)
      ax%x = 1
   endselect
   b = a
   b = a + b ! here the `+` operator could generate memory leaks
end subroutine polymorphic_leaks_raiser

subroutine inherit_leaks_raiser
    use inherit_type_m
    implicit none
    type(inherit_type_t) :: a
    type(inherit_type_t) :: b

    allocate(a%x)
    a%x = derived_t(1)
    b = a
    b = a + b ! here the `+` operator could generate memory leaks
end subroutine inherit_leaks_raiser

subroutine multiple_inheritance_leaks_raiser
   use child_m
   implicit none
   type(child_t) :: a
   type(child_t) :: b

   allocate(a%y)
   a%x = 1
   a%y = 2
   b = a
   b = a + b ! here the `+` operator could generate memory leaks
end subroutine multiple_inheritance_leaks_raiser


end program simple_test_memleaks
