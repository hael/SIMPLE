program simple_test_find_boundaries
    include 'simple_lib.f08'

    implicit none

    character(len=:), allocatable   :: fn_in, fn_out, python_cmd
    character(*), parameter         :: pyth_fn = "$SIMPLE_PATH/../production/tests/test_find_boundaries/my_python.py"
    character(100)                  :: min_size_char, pen_char
    integer                         :: i, min_size, fnlen
    real                            :: pen

    ! Read in command line
    if( command_argument_count() .NE. 4 )then
        write(logfhandle,'(a)') 'Usage: simple_test_find_boundaries reproj_vs_corrs_file.csv min_size pen'
        write(logfhandle,'(a)') 'Example: simple_test_find_boundaries /Users/.../NP87_corrs.csv 2500 5'
        write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
        fn_in = '/Users/wietfeldthc/Documents/singleSimulations/testBoundaryAlgAug12/time-series-corrs/NP1_71-17_corrs.csv'
        min_size = 2500
        pen = 5
    else
        call get_command_argument(1, length=fnlen)
        allocate(character(fnlen) :: fn_in)
        call get_command_argument(1, fn_in)
        call get_command_argument(2, length=fnlen)
        allocate(character(fnlen) :: fn_out)
        call get_command_argument(2, fn_out)
        call get_command_argument(3, min_size_char)
        read(min_size_char, *)min_size
        call get_command_argument(4, pen_char)
        read(pen_char, *)pen
    endif


    ! Call python program
    allocate(character(2 * fnlen + 100) :: python_cmd)
    python_cmd = "python "//trim(pyth_fn)//" "//trim(fn_in)//" "//trim(fn_out)//" "//trim(min_size_char)//" "//trim(pen_char)
    print *, trim(python_cmd)

    print *, "Calling python from Fortran"
    call execute_command_line(trim(python_cmd), exitstat=i)
    print *, "Exit status of python was ", i

end program simple_test_find_boundaries