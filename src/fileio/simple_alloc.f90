
#ifdef VARIADIC
#define SIMPLE_ALLOCATE(...) allocate( VA_ARGS , stat=n, errmsg=e1); \
call simple_error_check(n,e1, f=__FILENAME__, l=__LINE)
#endif

module simple_alloc
use simple_defs
use simple_syslib
implicit none
    ! interface simple_allocate
    !     module procedure simple_alloc_r4
    !     module procedure simple_alloc_r8
    !     module procedure simple_alloc_i4
    !     module procedure simple_alloc_i8
    !     module procedure simple_alloc_l
    !     module procedure simple_alloc_c4
    !     module procedure simple_alloc_c8
    !     module procedure simple_alloc_ch
    ! end interface simple_allocate


contains
          ! DEALLOC
      subroutine simple_dealloc_r4 (array)
       real, intent(in) :: array
       character(len=70) e1
       integer n,e
       e1 = 'No error'
       n = 99 ! dummy value
       if (allocated(array)) then
           deallocate(array, stat=n, errmsg=e1)
           if (n==0) return
           ! if (trim(e1) == 'No error') return
           if (trim(e1) == 'Attempt to deallocate an unallocated object') then
               write(*,*) " simple_dealloc; Attempt to deallocate an unallocated object"
               stop
           end if
           if (n < 0)then
               write(*,*) " simple_dealloc; negative err no: ", n
           end if
           if (n /= 0 )then
               e = get_sys_error()
               if( e /= 0) stop
           end if
       end if
   end subroutine simple_dealloc_r4



   ! Allocators
   subroutine simple_alloc_r4 (array,f,l)
       real, intent(in) :: array(:)
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4
   subroutine simple_alloc_r8 (array,f,l)
       real, intent(in) :: array(:)
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4
      subroutine simple_alloc_r4 (array,f,l)
       real, intent(in) :: array
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4
      subroutine simple_alloc_r4 (array,f,l)
       real, intent(in) :: array
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4
      subroutine simple_alloc_r4 (array,f,l)
       real, intent(in) :: array
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4
      subroutine simple_alloc_r4 (array,f,l)
       real, intent(in) :: array
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4
      subroutine simple_alloc_r4 (array,f,l)
       real, intent(in) :: array
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4
      subroutine simple_alloc_r4 (array,f,l)
       real, intent(in) :: array
       character(len=STDLEN), intent(in), optional :: f !< filename of caller
       integer,               intent(in), optional :: l !< line number from calling file
       character(len=70) emsg
       integer n
       e1 = 'No error'
       n = 99 ! dummy value
       if(.not.allocated(array)) allocate(array,  STAT=n, ERRMSG=emsg)
       if (n /= 0)then
           if (present(f).and.present(l))then
               call simple_stop("simple_alloc_r4 "//trim(emsg), f, l)
           else
               call simple_stop("simple_alloc_r4 "//trim(emsg))
           end if
       end if
   end subroutine simple_alloc_r4

end module simple_alloc
