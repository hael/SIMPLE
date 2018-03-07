module simple_test_export_jpg
    include 'simple_lib.f08'
    use simple_jpg,   only: jpg_img
    use simple_image, only: image
implicit none

#include "simple_local_flags.inc"
    public ::  test_jpg_image
contains

    subroutine test_jpg_image (doplot)
        logical, intent(in) :: doplot
        logical passed
        global_verbose=.true.
        write(*,'(a)') '**info(simple_jpg_unit_test): testing simple_jpg '
        call test_jpg_image_local( 91, 91, 100, doplot )
        write(*,'(a)') 'SIMPLE_JPG_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'

    contains

        subroutine test_jpg_image_local( ld1, ld2, ld3, doplot )
            integer, intent(in)  :: ld1, ld2,ld3
            logical, intent(in)  :: doplot
            type(image)          :: img, img3
            type(jpg_img)        :: jpg
            character(len=STDLEN) :: str
            real, allocatable    :: rbuf(:,:,:)
            real, allocatable :: rptr(:,:)
            integer, allocatable    :: ibuf(:,:,:)
            integer status

            write(*,'(a)') '**info(simple_jpg_unit_test, part 1)'
            call img%new([ld1,ld2,1], 1.)
            write(*,'(a)') '**info(simple_jpg_unit_test, part 2): testing real buffer write to JPG'
            passed = .true.
            call img%ran
            if( doplot ) call img%vis
            str= 'test_jpg_ran.jpg'
            rbuf = img%get_rmat()
            print *,rbuf
            allocate(rptr(ld1,ld2),source=rbuf(:,:,1))
            status =  jpg%writeJpgToFile(str,rptr, quality=100, colorspec=1)
            if(status /= 0)  call simple_stop('test_jpg_image_local write_jpeg returned error' )
            call exec_cmdline('display test_jpg_ran.jpg')

            write(*,'(a)') '**info(simple_jpg_unit_test, part 3): testing int buffer write to JPG'
            allocate(ibuf(size(rbuf,1),size(rbuf,2),size(rbuf,3)))
            ibuf = INT(rbuf * (2**12))
            str= 'test_jpg_ran_int.jpg'
            status = jpg%writeJpgToFile(str,ibuf(:,:,1))
            if(status /= 0)  call simple_stop('test_jpg_image_local write_jpeg int buffer failed ' )
            call exec_cmdline('display test_jpg_ran_int.jpg')



            write(*,'(a)') '**info(simple_jpg_unit_test, part 4): testing 3D real buffer write from get_rmat'
            call img3%new([21,20,10], 1.)
            call img3%gauran( 5., 15. )
            if( doplot ) call img%vis
            str= 'test_jpg_gauss.jpg'
            status =  jpg%writeJpgToFile(str,img3%get_rmat())
            if(status /= 0)  call simple_stop('test_jpg_image_local write_jpeg 3D buffer failed' )
            call exec_cmdline('montage -geometry +1+1 test_jpg_gauss*.jpg gauss.jpg && display gauss.jpg')


            if(allocated(rbuf)) deallocate(rbuf)
            if(allocated(ibuf)) deallocate(ibuf)
            if(allocated(rptr)) deallocate(rptr)
            call img%kill
            call img3%kill

        end subroutine test_jpg_image_local

    end subroutine test_jpg_image

end module simple_test_export_jpg
