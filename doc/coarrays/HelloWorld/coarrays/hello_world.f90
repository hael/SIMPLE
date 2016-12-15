program Hello_World
  implicit none
  integer :: i  ! Local variable
  character(len=20) :: name[*] ! scalar coarray
  ! Note: "name" is the local variable while "name[<index>]"
  ! accesses the variable on a remote image
 
  ! Interact with the user on Image 1
  if (this_image() == 1) then
    write(*,'(a)',advance='no') 'Enter your name: '
    read(*,'(a)') name
 
    ! Distribute information to other images
    do i = 2, num_images()
      name[i] = name
    end do
  end if
 
  sync all ! Barrier to make sure the data has arrived
 
  ! I/O from all nodes
  write(*,'(3a,i0)') 'Hello ',trim(name),' from image ', this_image()
end program Hello_world