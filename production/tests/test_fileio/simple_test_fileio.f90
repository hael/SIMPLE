program simple_test_fileio
    include 'simple_lib.f08'
    use iso_fortran_env
    use simple_testfuns          ! use all in there
    use simple_image,            only: image

    implicit none
    type( image )         :: cube, img
    real                  :: smpd
    integer               :: box, nspace, msk, i,j,istat
    integer               :: un, u(1200)
    character(len=8)      :: datestr
    character(len=STDLEN) :: folder,  oldfolder, curDir, cmd
    character(len=30)     :: fname

    global_verbose=.true.
    call seed_rnd
    call date_and_time(date=datestr)
    folder = trim('./SIMPLE_TEST_FILEIO_'//datestr)

    call simple_mkdir( trim(folder) )
    call simple_mkdir( trim(folder) ) !! double checking creation to print ignore statement
    print *," Listing of ", folder
    cmd='ls -al '//trim(folder)
    call execute_command_line(trim(cmd))
    print *," Changing directory to ", folder
    call simple_chdir( trim(folder),  oldfolder)
    call simple_getcwd(curDir)
    print *," Current working directory ", curDir
    print *," Previous working directory ", oldfolder


#if defined(PGI)
    print *,">>> Testing PGI STREAMING  "
    call fopen(fid,file='myfile.txt',access='stream',form='unformatted',iostat=ios,iomsg=msg)
    write(*,*) 'ios1=',ios
    write(*,*) 'mess1=',trim(msg)
    read(fid,pos=20,iostat=ios,iomsg=msg) c
    write(*,*) 'ios1=',ios
    write(*,*) 'mess1=',trim(msg)
    read(fid,pos=1,iostat=ios) c
    write(*,*) 'ios2=',ios
    write(*,*) 'mess2=',trim(msg)
    call fclose(fid)

    call fopen(un,file='test.bin',status='unknown',action='write',access='stream',form='unformatted',iostat=istat,iomsg=msg)
    if(istat /= 0) then
        call fileiochk(" PGI STREAMING failed ", istat)
        call simple_error_check()
        call gerror(msg)
        write(stderr,'("SIMPLE_SYSLIB::SYSERROR  ",I0)') istat
        write(stderr,*) trim(adjustl(msg))

    endif

    write(un) 1._dp,2._dp,3._dp
    call flush(un)
    call fclose(un,iostat=istat)
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

    write(*,*)'>>> Testing Read/Write accuracy (floats)'
    call test_readwrite_accuracy
    call simple_chdir('../')
    call simple_end('**** SIMPLE_TEST_FILEIO NORMAL STOP ****')


contains

    subroutine  test_readwrite_accuracy

        character(len=128) :: teststring
        integer, parameter :: nformats=4
        character(len=20), parameter :: formats(nformats) =   &
            [ '( E11.4)', '( E13.6)', '( E15.8)', '(E17.10)' ]
        real, dimension(nformats) :: errors

        real :: output, back
        real, parameter :: delta=epsilon(output)
        integer :: i,j
        real :: spac(100000)
        call random_number(spac)
        spac = spac / sum(spac)
        errors = 0.
        j=0
        output =  (1.0/exp(1.0))
        do j=1, size(spac,1)
            if (output >= 1) exit
            do i=1,nformats
                write(teststring,FMT=formats(i)) output
                read(teststring,*) back
                if (abs(back-output) > errors(i)) errors(i) = abs(back-output)
            enddo
            output = output + spac(j) ! delta
        end do

        print *, 'Maximum errors: '
        print *, formats
        print *, errors

        print *, 'Trying with default format: '
        j=0
        errors = 0.
        output = (1.0/exp(1.0))
        do j=1, size(spac,1)
            if (output >= 1) exit
            write(teststring,*) output
            read(teststring,*) back
            if (abs(back-output) > errors(1)) errors(1) = abs(back-output)
            output = output + spac(j) ! delta

        end do

        print *, 'Error = ', errors(1)
    end subroutine test_readwrite_accuracy

end program simple_test_fileio
