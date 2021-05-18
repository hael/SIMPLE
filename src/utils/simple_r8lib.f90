module simple_r8lib
implicit none

public :: r8mat_cholesky_factor, r8mat_cholesky_solve
private

contains

    ! R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix
    ! The matrix must be symmetric and positive semidefinite
    ! For a positive semidefinite symmetric matrix A, the Cholesky factorization is a lower triangular matrix L such that:
    ! A = L * L'
    ! An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M]
    ! Author: John Burkardt
    subroutine r8mat_cholesky_factor( n, a, c, flag )
        integer,         intent(in)    :: n      ! the number of rows and columns of the matrix a
        real(kind=8),    intent(in)    :: a(n,n) ! the n by n matrix
        real(kind=8),    intent(inout) :: c(n,n) ! c(n,n) the N by N lower triangular Cholesky factor
        integer,         intent(inout) :: flag   ! 0, no error occurred, 1, the matrix is not positive definite
        integer      :: i, j
        real(kind=8) :: sum2
        flag = 0
        c    = a
        do j = 1, n
            c(1:j-1,j) = 0.0D+00
            do i = j, n
                sum2 = c(j,i) - dot_product ( c(j,1:j-1), c(i,1:j-1) )
                if ( i == j ) then
                    if ( sum2 <= 0.0D+00 ) then
                        flag = 1
                        return
                    else
                        c(i,j) = sqrt ( sum2 )
                    end if
                else
                    if ( c(j,j) /= 0.0D+00 ) then
                        c(i,j) = sum2 / c(j,j)
                    else
                        c(i,j) = 0.0D+00
                    end if
                end if

            end do

        end do
    end subroutine r8mat_cholesky_factor

    ! R8MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b
    ! An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M]
    ! This routine works with the lower triangular Cholesky factor A = L * L'
    ! Author: John Burkardt
    subroutine r8mat_cholesky_solve( n, l, b, x )
        integer,      intent(in)    :: n      ! the number of rows and columns of the matrix A
        real(kind=8), intent(in)    :: l(n,n) ! the N by N lower Cholesky factor of the system matrix A
        real(kind=8), intent(in)    :: b(n)   ! the right hand side of the linear system
        real(kind=8), intent(inout) :: x(n)   ! the solution of the linear system
        ! Solve L * y = b.
        call r8mat_l_solve( n, l, b, x )
        ! Solve L' * x = y.
        call r8mat_lt_solve( n, l, x, x )
    end subroutine r8mat_cholesky_solve

    ! R8MAT_L_SOLVE solves a lower triangular linear system
    ! An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M]
    ! Author: John Burkardt
    subroutine r8mat_l_solve ( n, a, b, x )
        integer,      intent(in)    :: n      ! the number of rows and columns of the matrix A
        real(kind=8), intent(in)    :: a(n,n) ! the N by N lower triangular matrix
        real(kind=8), intent(in)    :: b(n)   ! the right hand side of the linear system
        real(kind=8), intent(inout) :: x(n)   ! the solution of the linear system
        integer :: i
        ! Solve L * x = b.
        do i = 1, n
            x(i) = ( b(i) - dot_product ( a(i,1:i-1), x(1:i-1) ) ) / a(i,i)
        end do
    end subroutine r8mat_l_solve

    ! R8MAT_LT_SOLVE solves a transposed lower triangular linear system
    ! An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M]
    ! Given the lower triangular matrix A, the linear system to be solved is:
    ! A' * x = b
    ! Author: John Burkardt
    subroutine r8mat_lt_solve( n, a, b, x )
        integer,      intent(in)    :: n      ! the number of rows and columns of the matrix
        real(kind=8), intent(in)    :: a(n,n) ! the N by N lower triangular matrix
        real(kind=8), intent(in)    :: b(n)   ! the right hand side of the linear system
        real(kind=8), intent(inout) :: x(n)   ! the solution of the linear system
        integer :: i
        ! Solve L'*x = b.
        do i = n, 1, -1
            x(i) = ( b(i) - dot_product ( x(i+1:n), a(i+1:n,i) ) ) / a(i,i)
        end do
    end subroutine r8mat_lt_solve

end module simple_r8lib
