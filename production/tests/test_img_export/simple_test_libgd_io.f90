!*******************************************************************************
! Small test program for the fortran module interfacing the GD library.
!*******************************************************************************

!! Modified by Michael Eager, Feb 2018
module simple_test_libgd_io
    include 'simple_lib.f08'
    use simple_img
    implicit none
!    public :: test_png_io
contains

    subroutine test_png_io
        use simple_img
        type(base_img)            :: image
        integer                   :: c1,c2,c3,c4,c5
        integer                   :: a1(4),a2(4)
        real(dp)                  :: ptsize = 30.0_dp
        real(dp)                  :: angle = 0.0_dp
        integer                   :: w,h
        integer                   :: r,x,y,status
        real, allocatable         :: buffer(:,:)
        character(len=256)        :: fout = 'gray-buffer.png'
        real                  :: grange(2)
        a1(:) = (/ 5,47,12,55 /)
        a2(:) = (/ 36,60,2,25 /)
        call create_raw_png_tmp

        call create_img(64,64,image)
        call allocate_color(image,0,0,0,c1)
        call allocate_color(image,255,255,255,c2)
        call allocate_color(image,255,0,0,c3)
        call allocate_color(image,255,0,255,c4)
        call allocate_color(image,0,255,0,c5)

        call draw_line(image,1,1,64,64,c2)
        call draw_rectangle(image,10,10,53,53,c2)
        call draw_filled_rectangle(image,30,30,43,53,c3)
        call draw_filled_polygon(image,4,a1,a2,c4)
        call fill_img(image,1,5,c5)
        call write_img_as_png(image,"xxx.png")
        call destroy_img(image)
        write(*,*) 'test png'

        call create_img_from_png("bbb",image)
        call write_img_as_png(image,"bbb-from-png.png")
        call destroy_img(image)
        write(6,*) 'test png to png'
        call create_img_from_png("bbb",image)
        status = greyscale(image)
        call write_img_as_png(image,"gray.png")
        call destroy_img(image)
        write(6,*) 'test png to grayscale'


        status=0
        call read_png("gray.png", buffer, status)
        if (status /= 0) write(*,*) "read_png failed"
        if (.not. allocated(buffer)) then
            write(*,*) "read_png failed, buffer not allocated"
        else
            print*, "read_png width ", size(buffer,1), " height ",  size(buffer,2)
            print*, "read_png max val ", MAXVAL(buffer)
            print *, buffer
        end if
        grange(:) = (/ 0., 255. /)
        call write_png(buffer, fout)! ,grange , status)
        if (status /= 0) write(*,*) "write_png failed"
        deallocate(buffer)
        write(*,*) 'test read/write png from buffers'
        if( .not. compare_imgs("gray.png", fout) ) write(*,*) "png read/wrbite failed - images differ"


        call create_img_from_png("bbb",image)
        call write_img_as_jpeg(image,"bbb.jpg",90)
        call destroy_img(image)
        call create_img_from_jpeg("bbb.jpg",image)
        call write_img_as_png(image,"bbb-from-jpeg.png")
        call write_img_as_png(image,"bbb-from-jpeg.jpg")
        call destroy_img(image)
        write(*,*) 'test jpg done'


#ifdef WITH_XPM
        call create_img_from_png("bbb",image)
        call write_img_as_xpm(image,"bbb.xpm")
        call destroy_img(image)
        call create_img_from_xpm("bbb.xpm",image)
        call write_img_as_png(image,"bbb-from-xpm.png")
        call destroy_img(image)
        write(*,*) 'test xpm done'
#endif

    end subroutine test_png_io

    subroutine test_jpeg_io
        use simple_img
        type(base_img)            :: image
        integer                   :: c1,c2,c3,c4,c5
        integer                   :: a1(4),a2(4)
        real(dp)                  :: ptsize = 30.0_dp
        real(dp)                  :: angle = 0.0_dp
        integer                   :: w,h
        integer                   :: r,x,y,status
        real, allocatable         :: buffer(:,:)
        character(len=256)        :: fout = 'gray-buffer.jpeg'
        real                  :: grange(2)
        a1(:) = (/ 5,47,12,55 /)
        a2(:) = (/ 36,60,2,25 /)
        call create_raw_png_tmp

        call create_img(64,64,image)
        call allocate_color(image,0,0,0,c1)
        call allocate_color(image,255,255,255,c2)
        call allocate_color(image,255,0,0,c3)
        call allocate_color(image,255,0,255,c4)
        call allocate_color(image,0,255,0,c5)

        call draw_line(image,1,1,64,64,c2)
        call draw_rectangle(image,10,10,53,53,c2)
        call draw_filled_rectangle(image,30,30,43,53,c3)
        call draw_filled_polygon(image,4,a1,a2,c4)
        call fill_img(image,1,5,c5)
        call write_img_as_jpeg(image,"xxx.jpeg")
        call destroy_img(image)
        write(*,*) 'test jpeg'


        call create_img_from_png("bbb",image)
        call write_img_as_jpeg(image,"bbb.jpg",70)
        call destroy_img(image)
        call create_img_from_jpeg("bbb.jpg",image)
        call write_img_as_png(image,"bbb-from-jpeg.png")
        call write_img_as_png(image,"bbb-from-jpeg.jpg")
        call destroy_img(image)
        write(*,*) 'test jpg done'

        call create_img_from_png("bbb",image)
        status = greyscale(image)
        call write_img_as_jpeg(image,"gray.jpg",100)
        call destroy_img(image)
        write(6,*) 'test jpeg to grayscale'

        status=0
        call read_jpeg("gray.jpg", buffer, status)
        if (status /= 0) write(*,*) "read_jpeg failed"
        if (.not. allocated(buffer)) then
            write(*,*) "read_jpeg failed, buffer not allocated"
        else
            print*, "read_jpeg width ", size(buffer,1), " height ",  size(buffer,2)
            print*, "read_jpeg max val ", MAXVAL(buffer)
            print *, buffer
        end if
        grange(:) = (/ 0., 255. /)
        call write_jpeg(buffer, fout)! ,grange , status)
        if (status /= 0) write(*,*) "write_jpeg failed"
        deallocate(buffer)
        write(*,*) 'test read/write jpeg from buffers'
        if( .not. compare_imgs("gray.jpg", fout) ) write(*,*) "jpeg read/write failed - images differ"

    end subroutine test_jpeg_io




    subroutine create_raw_png_tmp
        integer :: fid,ios
        integer(2) :: b
        integer(4) ::reclen

        fid=20
        inquire(iolength=reclen)b
        open(unit=20,file='bbb',form='unformatted',access='stream',status='replace',iostat=ios)
        !        open(NEWUNIT=fid,FILE="bbb.raw",ACCESS='STREAM',STATUS="REPLACE",IOSTAT=ios)!,FORM='UNFORMATTED')
        if (ios/=0) return
        write(fid) INT(Z'474e5089',4), INT(Z'0a1a0a0d',4), INT(Z'0d000000',4), INT(Z'52444849',4)
        write(fid) INT(Z'4b000000',4), INT(Z'4b000000',4), INT(Z'00000208',4), INT(Z'0ed2cb700',4)
        write(fid) INT(Z'000000bd',4), INT(Z'4d416704',4), INT(Z'0b1000041',4), INT(Z'61fc0b8f',4)
        write(fid) INT(Z'00000005',4), INT(Z'474b6206',4), INT(Z'00ff0044',4), INT(Z'a0ff00ff',4)
        write(fid) INT(Z'0093a7bd',4), INT(Z'70090000',4), INT(Z'00735948',4), INT(Z'00120b00',4)
        write(fid) INT(Z'01120b00',4), INT(Z'fc7eddd2',4), INT(Z'07000000',4), INT(Z'454d4974',4)
        write(fid) INT(Z'040cd207',4), INT(Z'952d3813',4), INT(Z'00759090',4), INT(Z'49520f00',4)
        write(fid) INT(Z'78544144',4), INT(Z'5d5bad9c',4), INT(Z'8db92393',4), INT(Z'752904cc',4)
        write(fid) INT(Z'763d78cf',4), INT(Z'e2f70fd8',4), INT(Z'ecebffff',4), INT(Z'67aedc70',4)
        write(fid) INT(Z'022aa5a7',4), INT(Z'048240f7',4), INT(Z'daef5459',4), INT(Z'42850ebe',4)
        write(fid) INT(Z'4aaaa51f',4), INT(Z'a0489122',4), INT(Z'e7f1f9aa',4), INT(Z'1fd5e33f',4)
        write(fid) INT(Z'cf6b14a9',4), INT(Z'a651a4a2',4), INT(Z'3414444f',4), INT(Z'7f9abf59',4)
        write(fid) INT(Z'bda1ed35',4), INT(Z'f5876995',4), INT(Z'97cd7f43',4), INT(Z'a20d22c7',4)
        write(fid) INT(Z'f6e1dbb2',4), INT(Z'c9062258',4), INT(Z'1c6f05bf',4), INT(Z'181dbeb2',4)
        write(fid) INT(Z'a72ec036',4), INT(Z'cd9660b5',4), INT(Z'c45d3f0e',4), INT(Z'fd758e4f',4)
        write(fid) INT(Z'02c73715',4), INT(Z'175c205f',4), INT(Z'6445ec29',4), INT(Z'77ec2af9',4)
        write(fid) INT(Z'5a72c08e',4), INT(Z'e8e5d722',4), INT(Z'30ec373f',4), INT(Z'fa436c8b',4)
        write(fid) INT(Z'41db1162',4), INT(Z'ee52d30d',4), INT(Z'd5e15e10',4), INT(Z'b721f99f',4)
        write(fid) INT(Z'2698deab',4), INT(Z'b58e43bb',4), INT(Z'0a346c68',4), INT(Z'79f8c202',4)
        write(fid) INT(Z'ac788cc3',4), INT(Z'5bb238cf',4), INT(Z'172389c4',4), INT(Z'cb725784',4)
        write(fid) INT(Z'a5a6e576',4), INT(Z'438cdbad',4), INT(Z'08ba9f8a',4), INT(Z'29dd8602',4)
        write(fid) INT(Z'fea21d86',4), INT(Z'b486b684',4), INT(Z'dc72eaff',4), INT(Z'6dd6dc4d',4)
        write(fid) INT(Z'f7168426',4), INT(Z'bf5a5322',4), INT(Z'88d92ecc',4), INT(Z'86e48c47',4)
        write(fid) INT(Z'004b08e1',4), INT(Z'b8b7dba9',4), INT(Z'276df585',4), INT(Z'63238983',4)
        write(fid) INT(Z'a8ea2287',4), INT(Z'36326976',4), INT(Z'5de02fce',4), INT(Z'5f55c581',4)
        write(fid) INT(Z'4f52a7e9',4), INT(Z'156f123a',4), INT(Z'dc6df8cc',4), INT(Z'7f514ce1',4)
        write(fid) INT(Z'b8ba9981',4), INT(Z'6cdb1bc0',4), INT(Z'3e4e3da3',4), INT(Z'698d7877',4)
        write(fid) INT(Z'9ac428b8',4), INT(Z'8c9e016f',4), INT(Z'b43ace43',4), INT(Z'4caed092',4)
        write(fid) INT(Z'783ca593',4), INT(Z'2ea05faa',4), INT(Z'eb145d21',4), INT(Z'952b7525',4)
        write(fid) INT(Z'49d20cc0',4), INT(Z'6deebafc',4), INT(Z'a22a94b4',4), INT(Z'93e30432',4)
        write(fid) INT(Z'74523532',4), INT(Z'bc2b88b3',4), INT(Z'828d6645',4), INT(Z'6d046c0e',4)
        write(fid) INT(Z'351788bf',4), INT(Z'75f241a2',4), INT(Z'00d75e25',4), INT(Z'86e29931',4)
        write(fid) INT(Z'd52fe9a1',4), INT(Z'516aec68',4), INT(Z'45f1b8a3',4), INT(Z'f550ba57',4)
        write(fid) INT(Z'858a1c54',4), INT(Z'71eac797',4), INT(Z'17c91c8b',4), INT(Z'9e25db4d',4)
        write(fid) INT(Z'53c1a188',4), INT(Z'd33d8617',4), INT(Z'b1ccfde6',4), INT(Z'24e7eabd',4)
        write(fid) INT(Z'2aa85130',4), INT(Z'7a7e26ad',4), INT(Z'a890e9e2',4), INT(Z'986d064b',4)
        write(fid) INT(Z'de02fcc8',4), INT(Z'3f3ad345',4), INT(Z'edf53d65',4), INT(Z'6de05ca2',4)
        write(fid) INT(Z'4fd17238',4), INT(Z'a2ea2a59',4), INT(Z'84514a46',4), INT(Z'50ef80bc',4)
        write(fid) INT(Z'3ae6eae4',4), INT(Z'37857c9c',4), INT(Z'51511fd4',4), INT(Z'14c9e038',4)
        write(fid) INT(Z'b8ae11b7',4), INT(Z'996c8bcd',4), INT(Z'58e11ece',4), INT(Z'13e4556c',4)
        write(fid) INT(Z'cd38e8f1',4), INT(Z'db03387f',4), INT(Z'532b2b1f',4), INT(Z'521e22dc',4)
        write(fid) INT(Z'6475f45d',4), INT(Z'4c06b6c5',4), INT(Z'eb28cd1e',4), INT(Z'4dbcfd0e',4)
        write(fid) INT(Z'1380a282',4), INT(Z'75eb3733',4), INT(Z'34a4ff25',4), INT(Z'de56e71d',4)
        write(fid) INT(Z'897cebfa',4), INT(Z'481cd5db',4), INT(Z'c49719f3',4), INT(Z'755c319b',4)
        write(fid) INT(Z'9933c6a5',4), INT(Z'd1534c10',4), INT(Z'6e335883',4), INT(Z'c548656e',4)
        write(fid) INT(Z'38cf8a1a',4), INT(Z'9f576982',4), INT(Z'22acc494',4), INT(Z'fe9c188a',4)
        write(fid) INT(Z'a04894dd',4), INT(Z'f19d5ca8',4), INT(Z'45450055',4), INT(Z'a9ea1254',4)
        write(fid) INT(Z'9ea2cf48',4), INT(Z'bb154d7b',4), INT(Z'15ac0236',4), INT(Z'6add6acb',4)
        write(fid) INT(Z'170df9c1',4), INT(Z'e3c51c6e',4), INT(Z'c062455a',4), INT(Z'8bd912f2',4)
        write(fid) INT(Z'e924f591',4), INT(Z'bc2e7660',4), INT(Z'26493363',4), INT(Z'd8be4902',4)
        write(fid) INT(Z'939f8d09',4), INT(Z'20b576b1',4), INT(Z'0f06803c',4), INT(Z'1a71d026',4)
        write(fid) INT(Z'e177b269',4), INT(Z'8a348021',4), INT(Z'b5cf0762',4), INT(Z'93a985f4',4)
        write(fid) INT(Z'de09cb41',4), INT(Z'0642a6ce',4), INT(Z'c99206bc',4), INT(Z'8c51ee7a',4)
        write(fid) INT(Z'0ccc0902',4), INT(Z'080357e7',4), INT(Z'04f11812',4), INT(Z'35a7c0a8',4)
        write(fid) INT(Z'ee0b63b6',4), INT(Z'4c5cddc7',4), INT(Z'de0afcb4',4), INT(Z'3f239ca4',4)
        write(fid) INT(Z'114bc15f',4), INT(Z'85bf4b5e',4), INT(Z'2392eea8',4), INT(Z'35e1cd3f',4)
        write(fid) INT(Z'40fe81a6',4), INT(Z'15e05b33',4), INT(Z'aea3545e',4), INT(Z'16e17dae',4)
        write(fid) INT(Z'ce115e3e',4), INT(Z'af3e4f75',4), INT(Z'4822627b',4), INT(Z'64839e15',4)
        write(fid) INT(Z'4464c01a',4), INT(Z'20d05fcf',4), INT(Z'a0b660f0',4), INT(Z'b47303b1',4)
        write(fid) INT(Z'ae4011e1',4), INT(Z'7c29a845',4), INT(Z'f9341b15',4), INT(Z'35960258',4)
        write(fid) INT(Z'b2a0caac',4), INT(Z'e17cfd04',4), INT(Z'49629457',4), INT(Z'1c94b606',4)
        write(fid) INT(Z'fd78725b',4), INT(Z'7354c105',4), INT(Z'365ec30f',4), INT(Z'0230b666',4)
        write(fid) INT(Z'd1468120',4), INT(Z'60048f06',4), INT(Z'5862a746',4), INT(Z'f27ba9df',4)
        write(fid) INT(Z'df02627d',4), INT(Z'1863a1ea',4), INT(Z'5d37b1d6',4), INT(Z'a35c3cb1',4)
        write(fid) INT(Z'0d99916f',4), INT(Z'426a245e',4), INT(Z'30477f6d',4), INT(Z'997592db',4)
        write(fid) INT(Z'710ccce9',4), INT(Z'00260383',4), INT(Z'1863009b',4), INT(Z'd64d4ba0',4)
        write(fid) INT(Z'd24478b9',4), INT(Z'1d6e29e5',4), INT(Z'1d8c1086',4), INT(Z'2bd725c0',4)
        write(fid) INT(Z'332cfeaa',4), INT(Z'395cfae0',4), INT(Z'930a022c',4), INT(Z'2a7a4440',4)
        write(fid) INT(Z'4ab12d56',4), INT(Z'4c011acb',4), INT(Z'8b7c43d1',4), INT(Z'06d6f9fe',4)
        write(fid) INT(Z'53630868',4), INT(Z'c1a55e9a',4), INT(Z'a3f2f4bc',4), INT(Z'77c64f76',4)
        write(fid) INT(Z'948fe781',4), INT(Z'42336066',4), INT(Z'2de12611',4), INT(Z'9319ba59',4)
        write(fid) INT(Z'006d7a41',4), INT(Z'1044c285',4), INT(Z'00350286',4), INT(Z'19b0a291',4)
        write(fid) INT(Z'5e367b6a',4), INT(Z'd148fa0b',4), INT(Z'603ce568',4), INT(Z'54ba7bb5',4)
        write(fid) INT(Z'a98b5dd5',4), INT(Z'24a9aca1',4), INT(Z'd092a4a4',4), INT(Z'd9823819',4)
        write(fid) INT(Z'a5290664',4), INT(Z'15b9367b',4), INT(Z'cd185313',4), INT(Z'89702aa8',4)
        write(fid) INT(Z'13313357',4), INT(Z'64bb8c93',4), INT(Z'd15a79ea',4), INT(Z'84005bca',4)
        write(fid) INT(Z'0c5af04d',4), INT(Z'cf7c160d',4), INT(Z'a20c1a55',4), INT(Z'bf5b1832',4)
        write(fid) INT(Z'b5d34a1a',4), INT(Z'ba7d5c14',4), INT(Z'4b354f77',4), INT(Z'21327283',4)
        write(fid) INT(Z'24d8adbe',4), INT(Z'7ee12dc9',4), INT(Z'66aa4d57',4), INT(Z'8da895a2',4)
        write(fid) INT(Z'b8d6cdae',4), INT(Z'99a99a8a',4), INT(Z'd92ee321',4), INT(Z'f47a71ca',4)
        write(fid) INT(Z'4ed21e0c',4), INT(Z'6cd59851',4), INT(Z'848e24bc',4), INT(Z'ac6978aa',4)
        write(fid) INT(Z'45e87f46',4), INT(Z'0c9c3662',4), INT(Z'616dd57a',4), INT(Z'1cc1e540',4)
        write(fid) INT(Z'6cb79514',4), INT(Z'a18b4de9',4), INT(Z'c8a9eb62',4), INT(Z'4424879a',4)
        write(fid) INT(Z'090b5d58',4), INT(Z'15403492',4), INT(Z'5bd2de96',4), INT(Z'8d99a729',4)
        write(fid) INT(Z'd68a4588',4), INT(Z'648b5d5c',4), INT(Z'9c778598',4), INT(Z'b4357e61',4)
        write(fid) INT(Z'5b4652f4',4), INT(Z'e9f8b4f7',4), INT(Z'29c1660a',4), INT(Z'53c13248',4)
        write(fid) INT(Z'19956243',4), INT(Z'2908824a',4), INT(Z'46802326',4), INT(Z'bb096b81',4)
        write(fid) INT(Z'6f1a9a72',4), INT(Z'b0e56ee9',4), INT(Z'1b7c245a',4), INT(Z'ba1bb198',4)
        write(fid) INT(Z'417b7f8a',4), INT(Z'10785dc1',4), INT(Z'0e1fda9b',4), INT(Z'1a0b7696',4)
        write(fid) INT(Z'025a7152',4), INT(Z'fab2a1e7',4), INT(Z'90c2f14f',4), INT(Z'c96f124c',4)
        write(fid) INT(Z'b62ac1e1',4), INT(Z'11449024',4), INT(Z'504e2f34',4), INT(Z'3e97aa08',4)
        write(fid) INT(Z'd6030aad',4), INT(Z'c33353fc',4), INT(Z'79cbdd2d',4), INT(Z'e1d0c6a2',4)
        write(fid) INT(Z'177af60d',4), INT(Z'6e757132',4), INT(Z'6c54a78b',4), INT(Z'c0e3e34a',4)
        write(fid) INT(Z'59829a75',4), INT(Z'563b2687',4), INT(Z'aae577cb',4), INT(Z'fb84a424',4)
        write(fid) INT(Z'd120169b',4), INT(Z'29899092',4), INT(Z'b08a2521',4), INT(Z'dd0e1481',4)
        write(fid) INT(Z'4e0482a2',4), INT(Z'8513b983',4), INT(Z'52dc8483',4), INT(Z'4c41672e',4)
        write(fid) INT(Z'65d5ab3f',4), INT(Z'efac6a22',4), INT(Z'a6d1872e',4), INT(Z'b429eec6',4)
        write(fid) INT(Z'4e9ea3a1',4), INT(Z'2aca28f0',4), INT(Z'2788bd2d',4), INT(Z'22092ef3',4)
        write(fid) INT(Z'a9aacf72',4), INT(Z'96638ae9',4), INT(Z'93124494',4), INT(Z'56e49250',4)
        write(fid) INT(Z'462c1515',4), INT(Z'a06ac255',4), INT(Z'33750660',4), INT(Z'3231430a',4)
        write(fid) INT(Z'0f599a67',4), INT(Z'a03335a0',4), INT(Z'59f7b00d',4), INT(Z'a6c8de01',4)
        write(fid) INT(Z'eeb57cd0',4), INT(Z'6bbdd469',4), INT(Z'8f3a160f',4), INT(Z'b71ac310',4)
        write(fid) INT(Z'88a15a4e',4), INT(Z'65169ba5',4), INT(Z'f0e8f06a',4), INT(Z'4b74a524',4)
        write(fid) INT(Z'fee92fb7',4), INT(Z'adc92226',4), INT(Z'7e271f0c',4), INT(Z'a613803c',4)
        write(fid) INT(Z'680d51d0',4), INT(Z'c9649fa8',4), INT(Z'c6c60505',4), INT(Z'08d8eb4c',4)
        write(fid) INT(Z'08fad505',4), INT(Z'20d552fb',4), INT(Z'ccc9a326',4), INT(Z'a0a65743',4)
        write(fid) INT(Z'5ee7e9af',4), INT(Z'eb6b7a5c',4), INT(Z'c2809f8d',4), INT(Z'8813525b',4)
        write(fid) INT(Z'ce22851a',4), INT(Z'e52924cc',4), INT(Z'a7f5fcb7',4), INT(Z'13292924',4)
        write(fid) INT(Z'a18d48b3',4), INT(Z'be437ee2',4), INT(Z'a61f5fe3',4), INT(Z'3885140f',4)
        write(fid) INT(Z'50806117',4), INT(Z'530ce401',4), INT(Z'66656a33',4), INT(Z'aad469c5',4)
        write(fid) INT(Z'9e5dad37',4), INT(Z'118c88dd',4), INT(Z'3978b74e',4), INT(Z'45e19aa2',4)
        write(fid) INT(Z'ba9bd813',4), INT(Z'8f6a2709',4), INT(Z'4df2c92b',4), INT(Z'd436bd15',4)
        write(fid) INT(Z'd29493c4',4), INT(Z'29fd7d3d',4), INT(Z'e5b96531',4), INT(Z'e52dfbf7',4)
        write(fid) INT(Z'ac2b251b',4), INT(Z'f2cf3ea8',4), INT(Z'2fa7f438',4), INT(Z'69e56138',4)
        write(fid) INT(Z'37514745',4), INT(Z'04c119a3',4), INT(Z'e8a0d981',4), INT(Z'0dab5450',4)
        write(fid) INT(Z'7a866864',4), INT(Z'a18bd84f',4), INT(Z'52ccdea7',4), INT(Z'32b631b6',4)
        write(fid) INT(Z'e85e84f6',4), INT(Z'e7eb4d29',4), INT(Z'21888a18',4), INT(Z'2515d52b',4)
        write(fid) INT(Z'cc6aa451',4), INT(Z'94924bed',4), INT(Z'dfdf96e4',4), INT(Z'5fef6fdf',4)
        write(fid) INT(Z'e617fbf3',4), INT(Z'f3e7666c',4), INT(Z'c987c790',4), INT(Z'5f2cf1f3',4)
        write(fid) INT(Z'462cf1df',4), INT(Z'66748cf4',4), INT(Z'0acc9842',4), INT(Z'02a3548c',4)
        write(fid) INT(Z'6c922f4a',4), INT(Z'3b33769e',4), INT(Z'4a1bd1ef',4), INT(Z'db42b1e3',4)
        write(fid) INT(Z'cccf1a7a',4), INT(Z'332e963f',4), INT(Z'4b490760',4), INT(Z'6d6ba735',4)
        write(fid) INT(Z'22915f6b',4), INT(Z'49efdd39',4), INT(Z'b7efde99',4), INT(Z'b7d7fbdb',4)
        write(fid) INT(Z'7feedf6f',4), INT(Z'f0fe96fd',4), INT(Z'fcf9c7ae',4), INT(Z'6ffc77e7',4)
        write(fid) INT(Z'1ca3862a',4), INT(Z'd66e9a99',4), INT(Z'008307f8',4), INT(Z'90819dd5',4)
        write(fid) INT(Z'97cd6a64',4), INT(Z'c169a1bb',4), INT(Z'6bb8d9ae',4), INT(Z'6270860c',4)
        write(fid) INT(Z'a6677c54',4), INT(Z'875b2c6c',4), INT(Z'4b4c31e7',4), INT(Z'a119d612',4)
        write(fid) INT(Z'71779b09',4), INT(Z'142cc68f',4), INT(Z'f94a5b91',4), INT(Z'dcbfdf9e',4)
        write(fid) INT(Z'cbedfaff',4), INT(Z'7f27fedf',4), INT(Z'1e3ea3fb',4), INT(Z'1cbdfdbc',4)
        write(fid) INT(Z'748e3cfa',4), INT(Z'6522796a',4), INT(Z'831061c3',4), INT(Z'984614d1',4)
        write(fid) INT(Z'550840d2',4), INT(Z'3690bb4c',4), INT(Z'439aa1c9',4), INT(Z'3ef0badd',4)
        write(fid) INT(Z'e96547f6',4), INT(Z'6026ea62',4), INT(Z'170dafb3',4), INT(Z'cfe9400f',4)
        write(fid) INT(Z'74a59194',4), INT(Z'bebd62e3',4), INT(Z'39892053',4), INT(Z'e7bc3fa7',4)
        write(fid) INT(Z'7fcc7f6f',4), INT(Z'fafe93f9',4), INT(Z'9dfcff23',4), INT(Z'96531bf7',4)
        write(fid) INT(Z'abbdea24',4), INT(Z'462a19a9',4), INT(Z'a9a4a66f',4), INT(Z'8c09a2a8',4)
        write(fid) INT(Z'0bc70d4e',4), INT(Z'05aa3e17',4), INT(Z'a9bb9215',4), INT(Z'b9eaa703',4)
        write(fid) INT(Z'77029b92',4), INT(Z'a80ef070',4), INT(Z'206ee9f6',4), INT(Z'6d402494',4)
        write(fid) INT(Z'f49e9af6',4), INT(Z'8f4c1b55',4), INT(Z'fd0f1f53',4), INT(Z'78fe87f5',4)
        write(fid) INT(Z'aa9a71e8',4), INT(Z'c082dcc1',4), INT(Z'4cd52500',4), INT(Z'008b8dd9',4)
        write(fid) INT(Z'f8979a50',4), INT(Z'96ebeb5e',4), INT(Z'2d13671b',4), INT(Z'bcff3f5b',4)
        write(fid) INT(Z'63f627e8',4), INT(Z'30c41a0b',4), INT(Z'e661b866',4), INT(Z'd2210029',4)
        write(fid) INT(Z'5ba481b4',4), INT(Z'8cfad516',4), INT(Z'3cf9cb53',4), INT(Z'cefcff9e',4)
        write(fid) INT(Z'cff3dfdb',4), INT(Z'3f38f5df',4), INT(Z'f1cbf1fe',4), INT(Z'394b3cf1',4)
        write(fid) INT(Z'0a3c918b',4), INT(Z'35c36814',4), INT(Z'8c494798',4), INT(Z'c76bf156',4)
        write(fid) INT(Z'1db714d6',4), INT(Z'0c860d62',4), INT(Z'734a1981',4), INT(Z'b92b04f5',4)
        write(fid) INT(Z'c81a3a7d',4), INT(Z'db646f49',4), INT(Z'42a48200',4), INT(Z'40528248',4)
        write(fid) INT(Z'c3140580',4), INT(Z'782ee0b3',4), INT(Z'0f1e9a9e',4), INT(Z'1cbdfcfc',4)
        write(fid) INT(Z'd4cdfbca',4), INT(Z'c7e78f8e',4), INT(Z'2cf3c78f',4), INT(Z'1c314147',4)
        write(fid) INT(Z'442c0505',4), INT(Z'f4a8d401',4), INT(Z'22da0916',4), INT(Z'8ba1b1d6',4)
        write(fid) INT(Z'983d7b2a',4), INT(Z'a477785c',4), INT(Z'cd4e1ebd',4), INT(Z'2f353bbb',4)
        write(fid) INT(Z'68b18132',4), INT(Z'9b38ab11',4), INT(Z'29404266',4), INT(Z'4a802490',4)
        write(fid) INT(Z'4b010195',4), INT(Z'bf8cebed',4), INT(Z'786f903e',4), INT(Z'14790a1c',4)
        write(fid) INT(Z'3279f3c3',4), INT(Z'52ce0325',4), INT(Z'e3f9e79e',4), INT(Z'797e20a1',4)
        write(fid) INT(Z'409c54f0',4), INT(Z'494c0a01',4), INT(Z'22101435',4), INT(Z'43a358a5',4)
        write(fid) INT(Z'faa0de43',4), INT(Z'ee1665eb',4), INT(Z'd450d46d',4), INT(Z'ac9a7bb9',4)
        write(fid) INT(Z'8ff51f4c',4), INT(Z'f4c20b1f',4), INT(Z'5b0b7ad6',4), INT(Z'6533203b',4)
        write(fid) INT(Z'b6b6ac06',4), INT(Z'38e21f6c',4), INT(Z'7e07cbf0',4), INT(Z'53c782fa',4)
        write(fid) INT(Z'674cf28f',4), INT(Z'c7150869',4), INT(Z'75fe2069',4), INT(Z'1389e3e0',4)
        write(fid) INT(Z'5027603c',4), INT(Z'9083f4e0',4), INT(Z'3cda2038',4), INT(Z'57826d1a',4)
        write(fid) INT(Z'859d67b5',4), INT(Z'a874ccbb',4), INT(Z'2da70f2d',4), INT(Z'a5fb691a',4)
        write(fid) INT(Z'5b2d5eb9',4), INT(Z'bbbbd206',4), INT(Z'84025a4e',4), INT(Z'66293782',4)
        write(fid) INT(Z'a4b50d04',4), INT(Z'943500c1',4), INT(Z'afb09d87',4), INT(Z'3b33705f',4)
        write(fid) INT(Z'8cc3b50f',4), INT(Z'5fc05924',4), INT(Z'c4f1f03e',4), INT(Z'4e04f0d3',4)
        write(fid) INT(Z'402c22b0',4), INT(Z'98a44d01',4), INT(Z'8f06c035',4), INT(Z'2e48a724',4)
        write(fid) INT(Z'9c1f414a',4), INT(Z'7dd152b6',4), INT(Z'09c1b8a9',4), INT(Z'18b9bf63',4)
        write(fid) INT(Z'51c47ef6',4), INT(Z'd7d5e885',4), INT(Z'43344caa',4), INT(Z'415f75aa',4)
        write(fid) INT(Z'81419a31',4), INT(Z'7b4f2862',4), INT(Z'e59ba5aa',4), INT(Z'70281464',4)
        write(fid) INT(Z'a7138e28',4), INT(Z'e00f00e2',4), INT(Z'4b1a3c69',4), INT(Z'50850059',4)
        write(fid) INT(Z'50d69c78',4), INT(Z'c9d182f5',4), INT(Z'c4664589',4), INT(Z'483d3f22',4)
        write(fid) INT(Z'ddca6db5',4), INT(Z'695ef4e0',4), INT(Z'118b182e',4), INT(Z'0679b027',4)
        write(fid) INT(Z'89b86289',4), INT(Z'29259b6a',4), INT(Z'62f46750',4), INT(Z'06290628',4)
        write(fid) INT(Z'19b2d6ab',4), INT(Z'869c42ec',4), INT(Z'9ed13803',4), INT(Z'901f580b',4)
        write(fid) INT(Z'1acf4c54',4), INT(Z'2b089f04',4), INT(Z'63eb393a',4), INT(Z'683efb52',4)
        write(fid) INT(Z'2f38ce54',4), INT(Z'0bc074be',4), INT(Z'352a663e',4), INT(Z'48e71523',4)
        write(fid) INT(Z'45326516',4), INT(Z'7c1817ba',4), INT(Z'e74cb13c',4), INT(Z'8425b521',4)
        write(fid) INT(Z'04ecc592',4), INT(Z'f0c6ba0b',4), INT(Z'90d01634',4), INT(Z'0d9c7820',4)
        write(fid) INT(Z'b2d022ad',4), INT(Z'f3f59636',4), INT(Z'c9eabc0b',4), INT(Z'53e89d62',4)
        write(fid) INT(Z'b0abb58c',4), INT(Z'97193c2b',4), INT(Z'73cceb6c',4), INT(Z'84e070c4',4)
        write(fid) INT(Z'86a4a416',4), INT(Z'b5a29a9b',4), INT(Z'1a50d503',4), INT(Z'd5a766b4',4)
        write(fid) INT(Z'4e61568e',4), INT(Z'3a275ac8',4), INT(Z'f3c6411e',4), INT(Z'b391b3b4',4)
        write(fid) INT(Z'55e1eaa1',4), INT(Z'4f3820a2',4), INT(Z'af107b7a',4), INT(Z'b307a1ae',4)
        write(fid) INT(Z'dd7a314d',4), INT(Z'b1cc3e93',4), INT(Z'665d6bd6',4), INT(Z'5c8f0de1',4)
        write(fid) INT(Z'1840d1e0',4), INT(Z'56f4b4f5',4), INT(Z'53bf5a41',4), INT(Z'30b5f6ae',4)
        write(fid) INT(Z'0f66d683',4), INT(Z'7ab4152e',4), INT(Z'197a3dea',4), INT(Z'0f0ed3cf',4)
        write(fid) INT(Z'2e188d94',4), INT(Z'7cc60f0f',4), INT(Z'f5a96d21',4), INT(Z'4b9e03a3',4)
        write(fid) INT(Z'c437a0a2',4), INT(Z'06893b9b',4), INT(Z'628c3ce9',4), INT(Z'4ea18d7e',4)
        write(fid) INT(Z'e961ca8d',4), INT(Z'62b7205e',4), INT(Z'8db3466c',4), INT(Z'add43d46',4)
        write(fid) INT(Z'c79a8585',4), INT(Z'c4b1b510',4), INT(Z'5261c4f3',4), INT(Z'fd97bebc',4)
        write(fid) INT(Z'978bc3a3',4), INT(Z'59aa03a5',4), INT(Z'874b3a72',4), INT(Z'dfabc128',4)
        write(fid) INT(Z'e553a8a2',4), INT(Z'c82a2f58',4), INT(Z'cf1bc998',4), INT(Z'6f0b9203',4)
        write(fid) INT(Z'43524337',4), INT(Z'6f06f098',4), INT(Z'10f44759',4), INT(Z'3544c493',4)
        write(fid) INT(Z'3d3c9c7b',4), INT(Z'4b03354a',4), INT(Z'a21524e3',4), INT(Z'75bc4280',4)
        write(fid) INT(Z'a6dce592',4), INT(Z'c917ce30',4), INT(Z'7a9d6873',4), INT(Z'0de0d897',4)
        write(fid) INT(Z'26f5337f',4), INT(Z'3879cc54',4), INT(Z'2f69485d',4), INT(Z'950f254e',4)
        write(fid) INT(Z'240ce60c',4), INT(Z'0ed68c83',4), INT(Z'0a02946d',4), INT(Z'3cf66a05',4)
        write(fid) INT(Z'639ab4ed',4), INT(Z'070b2d43',4), INT(Z'f38791da',4), INT(Z'58835b3c',4)
        write(fid) INT(Z'618e0084',4), INT(Z'86d981b1',4), INT(Z'd879eba9',4), INT(Z'c2ef71b3',4)
        write(fid) INT(Z'43004b95',4), INT(Z'19974e2a',4), INT(Z'22a84288',4), INT(Z'42486281',4)
        write(fid) INT(Z'9cfcd649',4), INT(Z'9eaa7029',4), INT(Z'0f192065',4), INT(Z'a3e30c7a',4)
        write(fid) INT(Z'ea30cd5a',4), INT(Z'4170ce5d',4), INT(Z'71e167d8',4), INT(Z'8cd6bc22',4)
        write(fid) INT(Z'23ef0b8d',4), INT(Z'86c04a74',4), INT(Z'189fecf0',4), INT(Z'd4b4688c',4)
        write(fid) INT(Z'b6534931',4), INT(Z'b9c756e7',4), INT(Z'422af1d6',4), INT(Z'f8f64213',4)
        write(fid) INT(Z'866fc687',4), INT(Z'964a2857',4), INT(Z'602134ae',4), INT(Z'2bc06bd7',4)
        write(fid) INT(Z'eb893ac2',4), INT(Z'8c6c7026',4), INT(Z'cd816f76',4), INT(Z'da290d5c',4)
        write(fid) INT(Z'c5ed4767',4), INT(Z'b5fa7de8',4), INT(Z'a8983087',4), INT(Z'3ca656f9',4)
        write(fid) INT(Z'86e88164',4), INT(Z'e7d17d9f',4), INT(Z'd094a40c',4), INT(Z'40041b10',4)
        write(fid) INT(Z'51b092be',4), INT(Z'32fe37b4',4), INT(Z'b7efa311',4), INT(Z'aa40b36d',4)
        write(fid) INT(Z'6056ab77',4), INT(Z'6983d9ed',4), INT(Z'8c2588f3',4), INT(Z'9ce7a3b1',4)
        write(fid) INT(Z'bfa2f610',4), INT(Z'bf56d787',4), INT(Z'e119facc',4), INT(Z'0febc801',4)
        write(fid) INT(Z'4b96f47a',4), INT(Z'f9b2f34e',4), INT(Z'60b4c9db',4), INT(Z'f67d0dc4',4)
        write(fid) INT(Z'fa192319',4), INT(Z'6c086a74',4), INT(Z'eaf31457',4), INT(Z'16da6b5f',4)
        write(fid) INT(Z'7d09ef5e',4), INT(Z'afd2d6d6',4), INT(Z'9769ed44',4), INT(Z'758466e5',4)
        write(fid) INT(Z'aa0c2c3d',4), INT(Z'21671f16',4), INT(Z'aa516192',4), INT(Z'31842382',4)
        write(fid) INT(Z'14b073f0',4), INT(Z'f5365635',4), INT(Z'44368eb2',4), INT(Z'31d273af',4)
        write(fid) INT(Z'e039b24f',4), INT(Z'78365a45',4), INT(Z'f008a28d',4), INT(Z'b6834b2a',4)
        write(fid) INT(Z'50e4a25e',4), INT(Z'b87c2e37',4), INT(Z'5307b739',4), INT(Z'4e8f4de2',4)
        write(fid) INT(Z'a756c748',4), INT(Z'32b31b9d',4), INT(Z'd395da28',4), INT(Z'9b3acb7a',4)
        write(fid) INT(Z'12cb4e90',4), INT(Z'03c0ba28',4), INT(Z'e6d3dde0',4), INT(Z'612e04a6',4)
        write(fid) INT(Z'7cf85d9c',4), INT(Z'908306cb',4), INT(Z'bf011e6c',4), INT(Z'c1b60977',4)
        write(fid) INT(Z'844a8936',4), INT(Z'9f27b977',4), INT(Z'5e1c5443',4), INT(Z'4a1bc1d9',4)
        write(fid) INT(Z'a17babb3',4), INT(Z'fba84478',4), INT(Z'2988037f',4), INT(Z'12b3343e',4)
        write(fid) INT(Z'd330dc15',4), INT(Z'0db1e2c0',4), INT(Z'26d1f4c1',4), INT(Z'296284c6',4)
        write(fid) INT(Z'8617360c',4), INT(Z'237a3af2',4), INT(Z'9378170f',4), INT(Z'0eab854d',4)
        write(fid) INT(Z'dc5c025a',4), INT(Z'298a45dc',4), INT(Z'ed13f13a',4), INT(Z'de198a1e',4)
        write(fid) INT(Z'7ee76468',4), INT(Z'2e4e19af',4), INT(Z'c74ccd37',4), INT(Z'fd55e19a',4)
        write(fid) INT(Z'28bbc3df',4), INT(Z'49dd38cd',4), INT(Z'ea821457',4), INT(Z'e426fd71',4)
        write(fid) INT(Z'5db1e0cf',4), INT(Z'51a8badc',4), INT(Z'4c972844',4), INT(Z'd5c60d5b',4)
        write(fid) INT(Z'0f75ad51',4), INT(Z'25ef0dbf',4), INT(Z'575366c2',4), INT(Z'e4eb1737',4)
        write(fid) INT(Z'fad5cdd6',4), INT(Z'd5c53d79',4), INT(Z'2f67ae40',4), INT(Z'4d77966d',4)
        write(fid) INT(Z'fc21dd13',4), INT(Z'2d66f0fd',4), INT(Z'776e858d',4), INT(Z'b6859437',4)
        write(fid) INT(Z'c49ce3a5',4), INT(Z'45bc2b00',4), INT(Z'3c15ed7e',4), INT(Z'7db6e60a',4)
        write(fid) INT(Z'e0253b9b',4), INT(Z'9fdeae2c',4), INT(Z'e879ebc3',4), INT(Z'dc5390c5',4)
        write(fid) INT(Z'aabccc96',4), INT(Z'3b78086b',4), INT(Z'c3cc5c1b',4), INT(Z'18e91c67',4)
        write(fid) INT(Z'6fb87630',4), INT(Z'af0bfa88',4), INT(Z'dfb106ec',4), INT(Z'8610d704',4)
        write(fid) INT(Z'175784eb',4), INT(Z'2fa43d7e',4), INT(Z'f84852f6',4), INT(Z'88b1bc39',4)
        write(fid) INT(Z'14fbd73e',4), INT(Z'31cca8de',4), INT(Z'3baa0368',4), INT(Z'b085dfb2',4)
        write(fid) INT(Z'ea902dc5',4), INT(Z'f71f5cd3',4), INT(Z'5ae30e73',4), INT(Z'b6299afb',4)
        write(fid) INT(Z'1b13877d',4), INT(Z'9df4e097',4), INT(Z'c01eac2c',4), INT(Z'e886a644',4)
        write(fid) INT(Z'b98787df',4), INT(Z'8f14bc5a',4), INT(Z'70ddbdeb',4), INT(Z'6c611d8d',4)
        write(fid) INT(Z'bc49790a',4), INT(Z'f546665d',4), INT(Z'e3146a8a',4), INT(Z'996f0b56',4)
        write(fid) INT(Z'd3f0f52c',4), INT(Z'd5c216ff',4), INT(Z'006aecca',4), INT(Z'3414f763',4)
        write(fid) INT(Z'06761bc6',4), INT(Z'bea80eea',4), INT(Z'1e7b8b0a',4), INT(Z'd6f9ac9f',4)
        write(fid) INT(Z'fe96a785',4), INT(Z'd0eff315',4), INT(Z'adc20d95',4), INT(Z'78f1372e',4)
        write(fid) INT(Z'97961ff9',4), INT(Z'3f1319a6',4), INT(Z'0f8ea2a3',4), INT(Z'edd9659b',4)
        write(fid) INT(Z'40a64247',4), INT(Z'fe7a5a56',4), INT(Z'f0f007ff',4), INT(Z'5c5a7baa',4)
        write(fid) INT(Z'a6f0e7f8',4), INT(Z'04bbf2a4',4), INT(Z'3dd7b63c',4), INT(Z'82d6022c',4)
        write(fid) INT(Z'f9d1e0b3',4), INT(Z'3c31ffe8',4), INT(Z'e97d2cac',4), INT(Z'd96ca7c2',4)
        write(fid) INT(Z'3b74eba4',4), INT(Z'9a44684e',4), INT(Z'cbaf5971',4), INT(Z'16dd65c6',4)
        write(fid) INT(Z'ff0ab886',4), INT(Z'd9003c0d',4), INT(Z'7c863524',4), INT(Z'7b76affa',4)
        write(fid) INT(Z'b2378a35',4), INT(Z'62e3965c',4), INT(Z'cede5b56',4), INT(Z'e4b1fdff',4)
        write(fid) INT(Z'ff88bae1',4), INT(Z'e18f3c25',4), INT(Z'5f3a9ff6',4), INT(Z'd35ed4d8',4)
        write(fid) INT(Z'17e2c3ae',4), INT(Z'7ac9fae9',4), INT(Z'fff6655d',4), INT(Z'ff801e19',4)
        write(fid) INT(Z'8a5bea03',4), INT(Z'9e9c3b83',4), INT(Z'0000009e',4), INT(Z'4e454900',4)
        write(fid) INT(Z'6042ae44',4), INT(Z'0082',2)
        call fclose(fid)
    end subroutine create_raw_png_tmp

end module simple_test_libgd_io
