program simple_test_fileio
    include 'simple_lib.f08'
    use iso_fortran_env
    use simple_testfuns          ! use all in there
    use simple_image,            only: image

    implicit none
    type( image )         :: cube, img
    real                  :: smpd
    integer               :: box, nspace, msk, i,j,istat,ios
    integer               :: un, u(1200)
    character             :: c
    character(len=8)      :: datestr
    character(len=STDLEN) :: folder, cmd
    character(len=300)    :: command
    CHARACTER(len=100)    :: msg
    character(len=30)     :: fname
    real(dp) :: tmp(5)
    integer(int32) :: ii
    real, dimension(10) :: a, b
    call seed_rnd
    call date_and_time(date=datestr)
    folder = trim('./SIMPLE_TEST_FILEIO_'//datestr)
    command = 'mkdir ' // trim( folder )//'|| true'
    call exec_cmdline( trim(command) )
    call simple_chdir( trim(folder) )

#if defined(PGI)
    print *,">>> Testing PGI STREAMING  "
    open(1,file='myfile.txt',access='stream',form='unformatted',iostat=ios,iomsg=msg)
    write(6,*) 'ios1=',ios
    write(6,*) 'mess1=',trim(msg)
    read(1,pos=20,iostat=ios,iomsg=msg) c
    write(6,*) 'ios1=',ios
    write(6,*) 'mess1=',trim(msg)
    read(1,pos=1,iostat=ios) c
    write(6,*) 'ios2=',ios
    write(6,*) 'mess2=',trim(msg)
    close(1)

    open(newunit=un,file='test.bin',status='unknown',action='write',access='stream',form='unformatted',iostat=istat,iomsg=msg)
    if(istat /= 0) then
        call fileiochk(" PGI STREAMING failed ", istat)
        call simple_error_check()
        call gerror(msg)
        write(stderr,'("SIMPLE_SYSLIB::SYSERROR  ",I0)') istat
        write(stderr,*) trim(adjustl(msg))

    endif

    write(un) 1._dp,2._dp,3._dp
    call flush(un)
    close(un,iostat=istat)
    open(un,file='test.bin',status='old',action='write',access='stream',form='unformatted',position='append',iostat=istat)
    write(un) 4._dp, 5._dp
    call flush(un)
    close(un,iostat=istat)
    open(un,file='test.bin',status='old',action='read',access='stream',form='unformatted',iostat=istat)
    read(un) tmp
    print *, tmp
    call flush(un)
    close(un,iostat=istat)

    print *,">>> Testing PGI STREAMING INT32 record test "
    call random_number(a)
    open (10,file='test.dat',form='unformatted',access='stream')
    inquire (iolength=ii) a
    write (10) ii, a, ii
    close (10)
    open (10,file='test.dat',form='unformatted')
    read (10) b
    if (all (a == b)) print *,'">>> Testing PGI STREAMING INT32  success!!'

    print *,">>> Testing PGI STREAMING  completed"
#endif



    print *,">>> Testing FOPEN  (Expecting Error 24)"
    do i=1,1200
        write(fname,'("tmp_",i0,".txt")') i
        call fopen(un,file=trim(adjustl(fname)),iostat=istat)
        if(istat/=0) then
            print *, "   Maximum number of open file objects: ", i, " newunit: ",un
            !!        call simple_error_check()
            exit
        end if
        u(i)=un
    end do
    print *,">>> Number of files opened: ", i
    print *,">>> Testing FCLOSE"
    do j=i,1,-1
        !!   print *," Closing ", u(j)
        call fclose(u(j))
    end do
    command = '/bin/rm -f tmp_*.txt || true'
    call exec_cmdline( trim(command) )
    print *,">>> Testing FOPEN/FCLOSE completed"







    ! dummy data
    box    = 96
    smpd   = 2.
    nspace = 64
    msk    = nint(real(box)/3.)
    ! volume
    call img%new( [box,box,box], smpd )
    call img%square( nint(real(box)/12.) )
    call cube%new( [box,box,box], smpd )
    call cube%square( nint(real(box)/16.) )
    call cube%shift([16.,16.,16.])
    call img%add( cube )
    call cube%new( [box,box,box], smpd )
    call cube%square( nint(real(box)/10.) )
    call cube%shift([4.,-16.,0.])
    call img%add( cube )
    call cube%kill
    call img%write( 'cubes.mrc' )
    call img%kill
    write(*,*)'>>> WROTE TEST VOLUME cubes.mrc'
    ! test units
    !command = 'simple_test_units'
    !call exec_cmdline( trim(command) )
    ! test search
    !command = 'simple_test_srch vol1=cubes.mrc msk='//int2str(msk)//&
    !    & ' smpd='//real2str(smpd)//' verbose=no'
    !call exec_cmdline( trim(command) )
    ! end
    call simple_chdir('../')
    call simple_end('**** SIMPLE_TEST_FILEIO NORMAL STOP ****')
end program simple_test_fileio
