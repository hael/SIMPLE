program simple_test_fileio
    include 'simple_lib.f08'
    use iso_fortran_env
    use simple_testfuns
    use simple_image, only: image
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

    print *,">>> Testing FOPEN  (Expecting Error 24)"
    do i=1,1200
        write(fname,'("tmp_",i0,".txt")') i
        call fopen(un,file=trim(adjustl(fname)),iostat=istat)
        if(istat/=0) then
            print *, "   Maximum number of open file objects: ", i, " newunit: ",un
            exit
        end if
        u(i)=un
    end do
    print *,">>> Number of files opened: ", i
    print *,">>> Testing FCLOSE"
    do j=i,1,-1
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
    write(logfhandle,*)'>>> WROTE TEST VOLUME cubes.mrc'
    write(logfhandle,*)'>>> Testing Read/Write accuracy (floats)'
    call test_readwrite_accuracy
    call simple_chdir(PATH_PARENT)
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
