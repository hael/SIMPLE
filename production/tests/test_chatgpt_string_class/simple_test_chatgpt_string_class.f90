program simple_test_chatgpt_string_class
    use simple_chatgpt_string_class
    implicit none

    ! type(String) :: s, num
    ! type(String), allocatable :: parts(:)
    ! type(String) :: text
    ! type(String), allocatable :: lines(:)
    ! integer :: i, iv
    ! real(8) :: rv

    ! s = string("   Hello, Fortran World!   ")
    ! print *, "Trimmed:", '"'//s%strip()%data//'"'
    ! print *, "Upper:", s%to_upper()%data
    ! print *, "Contains 'tran'?", s%contains("tran")

    ! s = string("1,2,3,4,5")
    ! parts = s%split(",")
    ! print *, "Split result:"
    ! do i = 1, size(parts)
    !     print *, i, ":", trim(parts(i)%data)
    ! end do

    ! num = string("42")
    ! iv = num%to_int()
    ! print *, "Integer value:", iv

    ! num = string("3.14159")
    ! rv = num%to_real()
    ! print *, "Real value:", rv

    ! print *, "Joined:", string("")%join(parts, ";")%data

    ! ! Write simple text
    ! text = string("Hello from Fortran!") // string(new_line('a')) // string("Second line.")
    ! call write_file("example.txt", text)

    ! ! Read entire file
    ! text = read_file("example.txt")
    ! print *, "Full file contents:"
    ! print *, trim(text%data)

    ! ! Read as lines
    ! lines = read_lines("example.txt")
    ! print *, "File lines:"
    ! print *, size(lines), "lines read."
    ! print *, trim(lines(1)%data)
    ! print *, trim(lines(2)%data)

    ! ! Append line
    ! call write_file("example.txt", string("Appended!"), append_mode=.true.)

end program simple_test_chatgpt_string_class
