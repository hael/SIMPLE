
module simple_test_jpg
    include 'simple_lib.f08'
    use simple_jpg
    use simple_image

    public ::  test_jpg_image
contains

    subroutine test_jpg_image (doplot)
        logical, intent(in) :: doplot
        write(*,'(a)') '**info(simple_jpg_unit_test): testing simple_jpg '
        call test_jpg_image_local( 100, 100, 100, doplot )
        write(*,'(a)') 'SIMPLE_JPG_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'

    contains

        subroutine test_jpg_image_local( ld1, ld2,  doplot )
            integer, intent(in)  :: ld1, ld2
            logical, intent(in)  :: doplot
            type(image)          :: img
            type(jpg_img)        :: jpg
            character(len=STDLEN) :: str
            real, allocatable    :: rbuf(:,:,;)
            integer, allocatable    :: ibuf(:,:,:)
            integer status
            write(*,'(a)') '**info(simple_jpg_unit_test, part 1)'
            call img%new([ld1,ld2,1], 1.)
            write(*,'(a)') '**info(simple_jpg_unit_test, part 2): testing real buffer write'
            passed = .true.
            call img%ran
            if( doplot ) call img%vis
            str= 'test_jpg_ran.jpg'//c_null_char
            rbuf = img%get_rmat()
            status =  jpg%save_jpeg(rbuf(:,:,1), str)
            if(status /= 0)  call simple_stop('test_jpg_image_local write_jpeg returned status==1' )

            write(*,'(a)') '**info(simple_jpg_unit_test, part 3): testing int buffer write'
            allocate(ibuf(size(rbuf,1),size(rbuf,2),size(rbuf,3)))
            ibuf = rbuf
            str= 'test_jpg_ran.jpg'//c_null_char
            status = jpg%save_jpeg(ibuf, str)
            if(status /= 0)  call simple_stop('test_jpg_image_local write_jpeg ibuf returned status==1' )

            write(*,'(a)') '**info(simple_jpg_unit_test, part 4): testing real buffer write from get_rmat'
            call img%gauran( 5., 15. )
            if( doplot ) call img%vis
            str= 'test_jpg_gauss.jpg'
            status =  jpg%save_jpeg(img%get_rmat(), str)

            if(allocated(rbuf)) deallocate(rbuf)
            if(allocated(ibuf)) deallocate(ibuf)


        end subroutine test_jpg_image_local

    end subroutine test_jpg_image

end module simple_test_jpg
