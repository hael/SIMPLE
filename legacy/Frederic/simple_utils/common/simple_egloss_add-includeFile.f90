  !local variables
  integer :: length,isize

  ! Consistency checks.
  isize=egloss_size(egloss)
  length=egloss_len(egloss)
  if (file_err(isize > 0,'Add not allowed for this node',&
       err_id=EGLOSS_INVALID_LIST)) return
  if (file_err(length == -1,'Add not allowed for this node',&
       err_id=EGLOSS_INVALID)) return

  if (present(last_item_ptr)) then
     ! Use sibling as last item to add to.
     if (associated(last_item_ptr)) then
        call init_next(last_item_ptr)
        call set_item(last_item_ptr, length)
     else
        last_item_ptr => egloss // 0
     end if
     call set(last_item_ptr, val)
  else
     ! Add new item at the end.
     call set(egloss//length,val)
  end if
