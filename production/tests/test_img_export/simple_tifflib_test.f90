
module simple_tifflib_test
#ifdef USING_TIFF
include 'simple_lib.f08'
use simple_tifflib
use gnufor2
implicit none
private
public :: test_tiff_write1, test_tiff_write2, test_tiff_write3, test_tiff_write4
public :: test_bigtiff_write, test_bigtiff_write1, test_bigtiff_write2, test_bigtiff_write3



contains
    subroutine test_tiff_write1()
        real, allocatable :: img(:,:,:)
        character(len=:), allocatable :: fname
        allocate(img(10,10,10))
        call random_number(img)
        allocate(fname,source='tifftest-1.tif')
        write(logfhandle,*) " Testing write_tiff(fname, img )"
        call write_tiff(fname, img )
        call gnufor_image(img(:,:,1), palette='gray')
        call exec_cmdline('display tifftest-1.tif')
        deallocate(img)
    end subroutine test_tiff_write1

    subroutine test_tiff_write2()
        integer, allocatable :: img(:)
        real :: temp
        integer :: i
        character(len=:), allocatable :: fname
        allocate(img(100))
        do i=1,10
            call random_number(temp)
            img(i) = INT( (2.**24) * temp)
        end do
        allocate(fname, source='tifftest-2.tif')
        write(logfhandle,*) " Testing write_tiff2(fname, img )"
        call write_tiff2(fname, img, [10,10] )

        call gnufor_image(reshape(real(img), shape=(/ 10,10 /)), palette='gray')
        call exec_cmdline('display tifftest-2.tif')

        deallocate(img)
    end subroutine test_tiff_write2

    subroutine test_tiff_write3()
        integer, allocatable :: img(:,:)
        real :: temp
        integer :: i,j
        character(len=:), allocatable :: fname
        allocate(img(10,10))
        do i=1,10
            do j=1,10
                call random_number(temp)
                img(i,j) = INT( (2.**24) * temp)
            end do
        end do
        allocate(fname, source='tifftest-3.tif')
         write(logfhandle,*) " Testing write_tiff3(fname, img )"
      !  call write_tiff3(fname, img )
        call gnufor_image(reshape(real(img), shape=(/ 10,10 /)), palette='gray')
        call exec_cmdline('display tifftest-3.tif')

        deallocate(img)
    end subroutine test_tiff_write3

    subroutine test_tiff_write4()
        real, allocatable :: img(:,:,:)
        character(len=:), allocatable :: fname

        allocate(img(10,10,10))
        call random_number(img)
        allocate(fname, source='tifftest-4.tif')
        write(logfhandle,*) " Testing write_tiff4(fname, img )"
        call write_tiff4(fname, img )

        call gnufor_image(img(:,:,1), palette='gray')
        call exec_cmdline('display tifftest-4.tif')
        deallocate(img)
    end subroutine test_tiff_write4

    subroutine test_bigtiff_write()
        real, allocatable :: img(:,:)
        real :: temp
        integer :: i,j
        character(len=:), allocatable :: fname
        allocate(img(100,100))
        do i=1,100
            do j=1,100
                call random_number(temp)
                img(i,j) =  (2.**26) * temp
            end do
        end do
        allocate(fname, source= 'tifftest-BIG.tif')
        write(logfhandle,*) " Testing write_bigtiff(fname, img )"
        call write_bigtiff(fname, img )
        call gnufor_image(img, palette='gray')
        call exec_cmdline('display tifftest-BIG.tif')

        deallocate(img)
    end subroutine test_bigtiff_write

    subroutine test_bigtiff_write1()
        real, allocatable :: img(:,:)
        integer :: status
        character(len=:), allocatable :: fname

        allocate(img(100,100))
        call random_number(img)
        allocate(fname, source='tifftest-BIG1.tif')
        write(logfhandle,*) " Testing write_tiff_bigimg(fname, img )"
      !  status =  write_tiff_bigimg('tifftest-BIG2.tif', img , 10,  10)
        call gnufor_image(img, palette='gray')
        call exec_cmdline('display tifftest-BIG1.tif')

        deallocate(img)
    end subroutine test_bigtiff_write1

    subroutine test_bigtiff_write2()
        real, allocatable :: img(:,:)
        integer :: status
        character(len=:), allocatable :: fname

        allocate(img(100,100))
        call random_number(img)
        allocate(fname, source='tifftest-BIG2.tif')
        write(logfhandle,*) " Testing write_bigtiff2(fname, img )"
        call write_bigtiff2(fname, img)
        call gnufor_image(img, palette='gray')
        call exec_cmdline('display tifftest-BIG2.tif')

        deallocate(img)
    end subroutine test_bigtiff_write2
    subroutine test_bigtiff_write3()
        real, allocatable :: img(:,:)
        integer :: status
        character(len=:), allocatable :: fname

        allocate(img(100,100))
        call random_number(img)
        allocate(fname, source='tifftest-BIG3.tif')
        write(logfhandle,*) " Testing write_bigtiff3(fname, img )"
        call write_bigtiff3(fname, img)
        call gnufor_image(img, palette='gray')
        call exec_cmdline('display tifftest-BIG3.tif')

        deallocate(img)
    end subroutine test_bigtiff_write3
#endif
end module simple_tifflib_test
