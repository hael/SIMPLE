
module simple_tifflib_test
use simple_tifflib
implicit none
private
public :: test_tiff_write1, test_tiff_write2, test_tiff_write3
public :: test_bigtiff_write, test_bigtiff_write2

contains
    subroutine test_tiff_write1()
        real, allocatable :: img(:,:,:)

        allocate(img(10,10,10))
        call random_number(img)
        call write_tiff('tifftest-1.tif', img )
        deallocate(img)
    end subroutine test_tiff_write1


    subroutine test_tiff_write2()
        integer, allocatable :: img(:)
        real :: temp
        integer :: i
        allocate(img(100))
        do i=1,10
            call random_number(temp)
            img(i) = INT( (2.**24) * temp)
        end do

        call write_tiff2('tifftest-2.tif', img, [10,10] )
        deallocate(img)
    end subroutine test_tiff_write2

    subroutine test_tiff_write3()
        integer, allocatable :: img(:,:)
        real :: temp
        integer :: i,j
        allocate(img(10,10))
        do i=1,10
            do j=1,10
                call random_number(temp)
                img(i,j) = INT( (2.**24) * temp)
            end do
        end do

        call write_tiff3('tifftest-3.tif', img )
        deallocate(img)
    end subroutine test_tiff_write3

    subroutine test_bigtiff_write()
        real, allocatable :: img(:,:)
        real :: temp
        integer :: i,j
        allocate(img(10,10))
        do i=1,10
            do j=1,10
                call random_number(temp)
                img(i,j) =  (2.**26) * temp
            end do
        end do

        call write_bigtiff('tifftest-BIG1.tif', img )
        deallocate(img)
    end subroutine test_bigtiff_write

    subroutine test_bigtiff_write2()
        real, allocatable :: img(:,:)
        integer :: status

        allocate(img(10,10))
        call random_number(img)
        status =  write_tiff_bigimg('tifftest-BIG2.tif', img, 10,  10, 16, 3, 0)
        deallocate(img)
    end subroutine test_bigtiff_write2


end module simple_tifflib_test
