!@descr: assertions
module simple_is_check_assert
use simple_defs
use simple_error, only: simple_exception
implicit none

interface is_a_number
    module procedure is_a_number_real, is_a_number_cmplx, is_a_number_dp, is_a_number_cmplx_dp 
end interface

interface is_zero
    module procedure is_zero_int, is_zero_real, is_zero_dp
end interface is_zero

interface is_gt_zero
    module procedure is_gt_zero_int, is_gt_zero_real, is_gt_zero_dp
end interface is_gt_zero

interface is_equal
    module procedure is_equal_real, is_equal_dp, is_equal_cmplx
end interface is_equal

interface is_even
    module procedure is_even_scalar, is_even_array
end interface is_even

interface check4nans3D
    module procedure check4nans3D_real, check4nans3D_cmplx
end interface

interface check4nans2D
    module procedure check4nans2D_real, check4nans2D_cmplx 
end interface

interface check4nans
    module procedure check4nans_real, check4nans_cmplx
end interface

interface assert_eq
    module procedure assert_eq2, assert_eq3, assert_eq4, assert_eqn
end interface assert_eq

contains

    elemental logical function is_even_scalar( val )
        integer, intent(in) :: val  !< query val
        is_even_scalar = mod(val,2) == 0
    end function is_even_scalar

    pure function is_even_array( arr ) result( yep )
        integer, intent(in) :: arr(:)     !< query vector
        logical :: yep
        logical :: test(size(arr))
        integer :: i
        test = .false.
        do i=1,size(arr)
            test(i) = is_even_scalar(arr(i))
        end do
        yep = all(test)
    end function is_even_array

    !>  returns true if the argument is odd
    elemental logical function is_odd(i)
        integer, intent(in) :: i    !< query value
        is_odd = btest(i,0)
    end function is_odd

    !>   checking for is_a_number
    elemental logical function is_a_number_real( number )
        real, intent(in) :: number  !< input variable for checking
        is_a_number_real = .true.
        if( number > 0. )then
        else if( number <= 0. )then
        else
            is_a_number_real = .false.
        endif
    end function is_a_number_real

    !>   validity check of complex number (so that it is not nan)
    elemental logical function is_a_number_cmplx( complex_number )
        complex, intent(in) :: complex_number !< input variable for checking
        is_a_number_cmplx = is_a_number_real(real(complex_number)) .and. is_a_number_real(aimag(complex_number))
    end function is_a_number_cmplx

     !>   checking for is_a_number
    elemental logical function is_a_number_dp( number )
        real(dp), intent(in) :: number  !< input variable for checking
        is_a_number_dp = .true.
        if( number > 0. )then
        else if( number <= 0. )then
        else
            is_a_number_dp = .false.
        endif
    end function is_a_number_dp

    !>   validity check of complex number (so that it is not nan)
    elemental logical function is_a_number_cmplx_dp( complex_number )
        complex(dp), intent(in) :: complex_number !< input variable for checking

        is_a_number_cmplx_dp = is_a_number_dp(real(complex_number)) .and. is_a_number_dp(aimag(complex_number))
    end function is_a_number_cmplx_dp

    !>   to check if val is zero
    elemental logical function is_zero_int( val )
        integer, intent(in) :: val  !< query val
        is_zero_int = abs(val) == 0
    end function is_zero_int

    !>   to check if val is zero
    elemental logical function is_zero_real( val )
        real, intent(in) :: val  !< query val
        is_zero_real = abs(val) < TINY
    end function is_zero_real

        !>   to check if val is zero
    elemental logical function is_zero_dp( val )
        real(8), intent(in) :: val  !< query val
        is_zero_dp = abs(val) < DTINY
    end function is_zero_dp

    !>   to check if val is zero
    elemental logical function is_gt_zero_int( val )
        integer, intent(in) :: val  !< query val
        is_gt_zero_int = val > 0
    end function

    !>   to check if val is zero
    elemental logical function is_gt_zero_real( val )
        real, intent(in) :: val  !< query val
        is_gt_zero_real = val > TINY
    end function is_gt_zero_real

    !>   to check if val is zero
    elemental logical function is_gt_zero_dp( val )
        real(8), intent(in) :: val  !< query val
        is_gt_zero_dp = val > DTINY
    end function is_gt_zero_dp

    !>   to check if val is zero
    elemental logical function is_equal_real( val1 , val2)
        real, intent(in) :: val1, val2  !< query val
        is_equal_real = abs(val1-val2) < TINY
    end function

    !>   to check if val is zero
    elemental logical function is_equal_dp( val1, val2 )
        real(8), intent(in) :: val1, val2  !< query val
        is_equal_dp = abs(val1-val2) < DTINY
    end function is_equal_dp

    !>   to check if magnitude of difference is zero
    elemental logical function is_equal_cmplx( val1, val2 )
        complex(4), intent(in) :: val1, val2  !< query val
        is_equal_cmplx = abs(val1-val2) < TINY
    end function is_equal_cmplx

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans3D_real( arr )
        real, intent(in)  :: arr(:,:,:)   !< query vector
        real, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_real(arr1d)
        deallocate(arr1d)
    end subroutine check4nans3D_real

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans3D_cmplx( arr )
        complex, intent(in)  :: arr(:,:,:)    !< query vector
        complex, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_cmplx(arr1d)
        deallocate(arr1d)
    end subroutine check4nans3D_cmplx

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans2D_real( arr )
        real, intent(in)  :: arr(:,:)   !< query vector
        real, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_real(arr1d)
        deallocate(arr1d)
    end subroutine check4nans2D_real

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans2D_cmplx( arr )
        complex, intent(in)  :: arr(:,:)   !< query vector
        complex, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_cmplx(arr1d)
        deallocate(arr1d)
    end subroutine check4nans2D_cmplx

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans_real( arr )
        real, intent(in) :: arr(:)   !< query vector
        integer :: i, n_nans
        n_nans = 0
        do i=1,size(arr)
            if( is_a_number(arr(i)) )then
                ! alles gut
            else
                n_nans = n_nans+1
            endif
        end do
        if( n_nans > 0 )then
            write(logfhandle,*) 'found NaNs in inputted vector; simple_math::check4nans_real', n_nans
        endif
    end subroutine check4nans_real

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans_cmplx( arr )
        complex, intent(in) :: arr(:)   !< query vector
        integer :: i, n_nans
        n_nans = 0
        do i=1,size(arr)
            if( is_a_number_cmplx(arr(i)) )then
                ! alles gut
            else
                n_nans = n_nans+1
            endif
        end do
        if( n_nans > 0 )then
            write(logfhandle,*) 'found NaNs in inputted vector; simple_math::check4nans_cmplx', n_nans
        endif
    end subroutine check4nans_cmplx

    function assert_eq2(n1,n2,string)
        implicit none
        character(len=*), intent(in) :: string
        integer,          intent(in) :: n1,n2
        integer :: assert_eq2
        if (n1 == n2) then
            assert_eq2=n1
        else
            write(logfhandle,*)'program terminated by assert_eq2: ',trim(string)
            stop 'program terminated by simple_math :: assert_eq2'
        end if
    end function assert_eq2

    function assert_eq3(n1,n2,n3,string)
        implicit none
        character(len=*), intent(in) :: string
        integer,          intent(in) :: n1,n2,n3
        integer :: assert_eq3
        if (n1 == n2 .and. n2 == n3) then
            assert_eq3=n1
        else
            write(logfhandle,*)'program terminated by assert_eq3: ',trim(string)
            stop 'program terminated by simple_math :: assert_eq3'
        end if
    end function assert_eq3

    function assert_eq4(n1,n2,n3,n4,string)
        implicit none
        character(len=*), intent(in) :: string
        integer,          intent(in) :: n1,n2,n3,n4
        integer :: assert_eq4
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
            assert_eq4=n1
        else
           write(logfhandle,*)'program terminated by assert_eq4: ',trim(string)
           stop 'program terminated by assert_eq4'
        end if
    end function assert_eq4

    function assert_eqn(nn,string)
        implicit none
        character(len=*),      intent(in) :: string
        integer, dimension(:), intent(in) :: nn
        integer :: assert_eqn
        if (all(nn(2:) == nn(1))) then
            assert_eqn=nn(1)
        else
            write(logfhandle,*)'program terminated by assert_eqn:', trim(string)
            stop 'program terminated by assert_eqn'
        end if
    end function assert_eqn

end module simple_is_check_assert