!! example 2.10.1  Vector Addition on the GPU
module vecaddmod
use openacc
implicit none
contains
    subroutine vecaddgpu( r, a, b, n )
        real, dimension(:) :: r, a, b
        integer :: n
        integer :: i

#ifdef PGI
        !$acc kernels loop  copyin(a(1:n),b(1:n)) copyout(r(1:n))
        do i = 1, n
            r(i) = a(i) + b(i)
        enddo
#elif GNU
#if __GNUC__ >= 6
        if(openacc_version >= 201306) then
            !$acc kernels loop  copyin(a(1:n),b(1:n)) copyout(r(1:n))
            do i = 1, n
                r(i) = a(i) + b(i)
            enddo
        endif
#else
        ! unsupported OpenACC in older GCC
#endif
#endif

    end subroutine vecaddgpu
    subroutine vecaddgpu_nocopy( r, a, b, n )
        real, dimension(:) :: r, a, b
        integer :: n
        integer :: i
#ifdef PGI
        !$acc kernels loop present(r,a,b)
         do i = 1, n
            r(i) = a(i) + b(i)
        enddo
#elif GNU
#if __GNUC__ >= 6
        if(openacc_version >= 201306) then
            !$acc kernels loop present(r,a,b)
            do i = 1, n
                r(i) = a(i) + b(i)
            enddo
        endif
#else
        ! unsupported OpenACC in older GCC
#endif
#endif
        do i = 1, n
            r(i) = a(i) + b(i)
        enddo
    end subroutine vecaddgpu_nocopy
end module vecaddmod

module simple_oacc_vecadd
    use openacc
    use vecaddmod
    implicit none

    public :: test_oacc_vecadd, test_oacc_vecadd_nocopy
contains
    subroutine test_oacc_vecadd (n_in)
        integer, intent(in), optional :: n_in
        integer :: i, errs, n
        real, dimension(:), allocatable :: a, b, r, e

#ifdef PGI
        print *,'test_oacc_vecadd: PGI implementation'
#elif GNU
if(openacc_version >= 201306) then
        print *,'test_oacc_vecadd: GNU Fortran implementation with OpenACC '
else
        print *,'test_oacc_vecadd: GNU Fortran implementation with older OpenACC '
        ! unsupported OpenACC in older GCC
endif
#endif


        n = 1000000 ! default value
        if(present(n_in)) n = n_in
        if( n <= 0 ) n = 100000

        allocate( a(n), b(n), r(n), e(n) )
        do i = 1, n
            a(i) = i
            b(i) = 1000*i
        enddo
        ! compute on the GPU
        call vecaddgpu( r, a, b, n )
        ! compute on the host to compare
        do i = 1, n
            e(i) = a(i) + b(i)
        enddo
        ! compare results
        errs = 0
        do i = 1, n
            if( r(i) /= e(i) )then
                errs = errs + 1
            endif
        enddo
        print *,'test_oacc_vecadd ',errs, ' errors found'
        if( errs /= 0) call exit(errs)


    end subroutine test_oacc_vecadd

    subroutine test_oacc_vecadd_nocopy (n_in)
        integer, intent(in), optional :: n_in
        integer :: n, i, errs
        real, dimension(:), allocatable :: a, b, r, e
#ifdef PGI
        print *,'test_oacc_vecadd_nocopy: PGI implementation'
#elif GNU

        if(openacc_version >= 201306) then
            print *,'test_oacc_vecadd_nocopy: GNU Fortran implementation with OpenACC '
        else
            print *,'test_oacc_vecadd_nocopy: GNU Fortran implementation with older OpenACC '
            ! unsupported OpenACC in older GCC
        endif
#endif
        n = 1000000 ! default value
        if(present(n_in)) n = n_in
        if( n <= 0 ) n = 100000
        allocate( a(n), b(n), r(n), e(n) )
        do i = 1, n
            a(i) = i
            b(i) = 1000*i
        enddo
        ! compute on the GPU
        !$acc data copyin(a,b) copyout(r)
        call vecaddgpu_nocopy( r, a, b, n )
        !$acc end data
        ! compute on the host to compare
        do i = 1, n
            e(i) = a(i) + b(i)
        enddo
        ! compare results
        errs = 0
        do i = 1, n
            if( r(i) /= e(i) )then
                errs = errs + 1
            endif
        enddo
        print *, 'test_oacc_vecadd_nocopy ', errs, ' errors found'
        if( errs /= 0) call exit(errs)
    end subroutine test_oacc_vecadd_nocopy


    subroutine test_vecadd_ornl
        ! Size of vectors
        integer :: n = 100000

        ! Input vectors
        real(8),dimension(:),allocatable :: a
        real(8),dimension(:),allocatable :: b
        ! Output vector
        real(8),dimension(:),allocatable :: c

        integer :: i
        real(8) :: sum

        ! Allocate memory for each vector
        allocate(a(n))
        allocate(b(n))
        allocate(c(n))

        ! Initialize content of input vectors, vector a[i] = sin(i)^2 vector b[i] = cos(i)^2
        do i=1,n
            a(i) = sin(i*1D0)*sin(i*1D0)
            b(i) = cos(i*1D0)*cos(i*1D0)
        enddo

        ! Sum component wise and save result into vector c

        !$acc kernels copyin(a(1:n),b(1:n)), copyout(c(1:n))
        do i=1,n
            c(i) = a(i) + b(i)
        enddo
        !$acc end kernels

        ! Sum up vector c and print result divided by n, this should equal 1 within error
        do i=1,n
            sum = sum +  c(i)
        enddo
        sum = sum/n
        print *, 'final result: ', sum

        ! Release memory
        deallocate(a)
        deallocate(b)
        deallocate(c)
    end subroutine test_vecadd_ornl


    subroutine test_reduction_ornl
        ! Size of vectors
        integer :: n = 100000

        ! Input vectors
        real(8),dimension(:),allocatable :: a
        real(8),dimension(:),allocatable :: b
        ! Output vector
        real(8),dimension(:),allocatable :: c

        integer :: i
        real(8) :: sum

        ! Allocate memory for each vector
        allocate(a(n))
        allocate(b(n))
        allocate(c(n))

        ! Initialize content of input vectors, vector a[i] = sin(i)^2 vector b[i] = cos(i)^2
        do i=1,n
            a(i) = sin(i*1D0)*sin(i*1D0)
            b(i) = cos(i*1D0)*cos(i*1D0)
        enddo

        ! Sum component wise and save result into vector c, perform reduction sum on c

        !$acc kernels !copyin(a(1:n),b(1:n)), copyout(c(1:n))
        do i=1,n
            c(i) = a(i) + b(i)
            sum = sum +  c(i)
        enddo
        !$acc end kernels

        ! Divide the sum by n, the final result should be 1.0 within error
        sum = sum/n
        print *, 'final result: ', sum

        ! Release memory
        deallocate(a)
        deallocate(b)
        deallocate(c)
    end subroutine test_reduction_ornl


end module simple_oacc_vecadd
