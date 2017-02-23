!==Class simple_jiffys
!
! simple_jiffys provides jiffys. The code is distributed with the hope that it will be useful, 
! but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification is regulated by the GNU General 
! Public License. *Author:* Hans Elmlund, 2011-08-18.
! 
!==Changes are documented below
!
module simple_jiffys
use simple_defs         ! singleton
use simple_filehandling ! singleton
implicit none

interface assert_eq
    module procedure assert_eq_2,assert_eq_3,assert_eq_4,assert_eq_n
end interface assert_eq

interface swap
    module procedure swap_i,swap_r,swap_rv,swap_c, swap_cv,swap_cm,&
    masked_swap_rs,masked_swap_rv,masked_swap_rm
end interface swap

contains

    ! FILE-HANDLING JIFFYS

    !>  \brief  is for finding logical dimension and numbe rof particles in stack
    subroutine find_ldim_nptcls( fname, ldim, nptcls, doprint, formatchar, endconv )
        character(len=*), intent(in)                         :: fname      !< filename
        integer, intent(out)                                 :: ldim(3)    !< logical dimension
        integer, intent(out)                                 :: nptcls     !< number of particles
        logical, intent(in), optional                        :: doprint    !< do print or not
        character(len=1), intent(in), optional               :: formatchar !< input format
        character(len=:), intent(out), allocatable, optional :: endconv    !< endian conversion
        integer                       :: mode, iform, maxim
        real                          :: smpd
        character(len=:), allocatable :: conv
        character(len=1)              :: form
        logical                       :: ddoprint
        ddoprint = .false. 
        if( present(doprint) ) ddoprint = doprint
        if( present(formatchar) )then
            form = formatchar
        else
            form = fname2format(fname)
        endif
        nptcls = 0
        select case (form)
            case('M','F')
                call get_mrcfile_info(fname, ldim, form, smpd, ddoprint )
                nptcls = ldim(3)
            case('S')
                call get_spifile_info(fname, ldim, iform, maxim, smpd, conv, ddoprint)
                nptcls = maxim
            case DEFAULT
                write(*,*) 'fname: ', fname
                write(*,*) 'format descriptor: ', fname2format(fname)
                stop 'File format not supported; find_ldim_nptcls; simple_procimgfile'
        end select
        if( present(endconv) )then
            if( allocated(endconv) ) deallocate(endconv)
            select case (form)
                case('M','F')
                    allocate(endconv, source='NATIVE')
                case('S')
                    allocate(endconv, source=conv)
            end select
        endif
    end subroutine find_ldim_nptcls
    
    !>  \brief is for gettign a part of the info in a MRC image header
    subroutine get_mrcfile_info( fname, ldim, form, smpd, doprint )
        use simple_imghead, only: ImgHead, MrcImgHead, MrcFeiImgHead
        character(len=*), intent(in)  :: fname
        character(len=1), intent(in)  :: form
        integer,          intent(out) :: ldim(3)
        real,             intent(out) :: smpd
        logical,          intent(in)  :: doprint
        class(imghead), allocatable   :: hed
        integer :: filnum
        ldim = 0
        smpd = 0.
        if( file_exists(fname) )then
            select case(form)
                case('M')
                    allocate(MrcImgHead :: hed)
                    call hed%new
                    filnum = get_fileunit()
                    open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='NATIVE')
                    call hed%read(filnum)
                    close(filnum)
                    ldim = hed%getDims()
                    smpd = hed%getPixSz()
                    if( doprint )then
                        call hed%print
                        write(*,'(a,3(i0,1x))') 'Number of columns, rows, sections: ', ldim(1), ldim(2), ldim(3)
                        write(*,'(a,1x,f15.8)')  'Pixel size: ', smpd
                    endif
                case('F')
                    allocate(MrcImgHead :: hed)
                    call hed%new
                    filnum = get_fileunit()
                    open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='BIG_ENDIAN')
                    call hed%read(filnum)
                    close(filnum)
                    if( doprint ) call hed%print
                case DEFAULT
                    write(*,*) 'The inputted file is not an MRC file; get_mrcfile_info; simple_jiffys'
                    write(*,*) fname
                    stop
            end select
        else
            write(*,*) 'The below file does not exists; get_mrcfile_info; simple_jiffys'
            write(*,*) fname
            stop
        endif
    end subroutine get_mrcfile_info
    
    !>  \brief is for gettign a part of the info in a SPIDER image header
    subroutine get_spifile_info( fname, ldim, iform, maxim, smpd, conv, doprint )
        character(len=*), intent(in)               :: fname
        integer, intent(out)                       :: ldim(3), maxim, iform
        real, intent(out)                          :: smpd
        character(len=:), allocatable, intent(out) :: conv
        logical, intent(in)                        :: doprint
        real    :: spihed(40)
        integer :: filnum, cnt, i
        if( file_exists(fname) )then
            if( fname2format(fname) .eq. 'S' )then
                if( allocated(conv) ) deallocate(conv)
                filnum = get_fileunit()
                open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='NATIVE')
                call read_spihed
                close(filnum)
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='NATIVE')
                    call print_spihed
                    return
                endif
                open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='BIG_ENDIAN') ! 
                call read_spihed
                close(filnum)
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='BIG_ENDIAN')
                    call print_spihed
                    return
                endif
                open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='LITTLE_ENDIAN') ! 
                call read_spihed
                close(filnum)
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='LITTLE_ENDIAN')
                    call print_spihed
                    return
                endif
            else
                write(*,*) 'The inputted file is not a SPIDER file; get_spifile_info; simple_jiffys'
                write(*,*) fname
                stop
            endif
        else
            write(*,*) 'The below file does not exists; get_spifile_info; simple_jiffys'
            write(*,*) fname
            stop
        endif
        
        contains
            
            subroutine read_spihed
                cnt = 0
                do i=1,40*4,4
                    cnt = cnt+1
                    read(unit=filnum ,pos=i) spihed(cnt)
                end do
                ldim  = int([spihed(12), spihed(2), spihed(1)])
                iform = int(spihed(5))
                maxim = int(spihed(26))
                smpd  = spihed(38)
            end subroutine
            
            subroutine print_spihed
                if( doprint )then
                    write(*,'(a,3(i0,1x))') 'Number of columns, rows, sections: ', int(spihed(12)), int(spihed(2)), int(spihed(1))
                    write(*,'(a,1x,i3)')    'Iform descriptor: ', int(spihed(5))
                    write(*,'(a,1x,f7.0)')  'The number of the highest image currently used in the stack: ', spihed(26)
                    write(*,'(a,1x,f7.3)')  'Pixel size: ', spihed(38)
                endif
            end subroutine
            
    end subroutine get_spifile_info

    !> \brief  for reading raw images using stream access
    subroutine read_raw_image( fname, mat, first_byte )
        character(len=*), intent(in)  :: fname
        double precision, intent(out) :: mat(:,:,:)
        integer, intent(in)           :: first_byte
        integer :: filnum, io_stat
        character(len=100) :: io_message
        filnum = get_fileunit()
        open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='NATIVE')
        read(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
        ! Check the read was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', fname
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; read_raw_image; simple_jiffys'
        endif
        close(filnum)
    end subroutine read_raw_image
    
    !> \brief  for writing raw images using stream access
    subroutine write_raw_image( fname, mat, first_byte )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: mat(:,:,:)
        integer,          intent(in) :: first_byte
        integer :: filnum, io_stat
        character(len=100) :: io_message
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file=fname, access='STREAM')
        write(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
        ! Check the write was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', fname
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; read_raw_image; simple_jiffys'
        endif
        close(filnum)
    end subroutine write_raw_image
    
    ! EXCEPTION-HANDLING JIFFYS
    
    !> \brief  is for raising command line exception
    subroutine cmdline_err( cmdstat, cmdlen, arg, pos )
        integer, intent(in)          :: cmdstat, cmdlen, pos
        character(len=*), intent(in) :: arg
        if( cmdstat == -1 )then
            write(*,*) 'ERROR! while parsing the command line; simple_exec'
            write(*,*) 'The string length of argument: ', arg, 'is: ', cmdlen
            write(*,*) 'which likely exeeds the lenght limit STDLEN'
            write(*,*) 'Create a symbolic link with shorter name in the cwd'
            stop
        endif
        if( arg(:pos-1) .ne. 'prg' )then
            write(*,'(a)') 'ERROR!'
            write(*,'(a)') 'prg=simple_program required to be first on command line'
            write(*,'(a)') 'Please, refer to the manual for a comprehensive '
            write(*,'(a)') 'list of all programs and their specific documentation'
            stop
        endif
    end subroutine cmdline_err

    !> \brief  is for checking allocation 
    subroutine alloc_err( message, alloc_stat )
        character(len=*), intent(in) :: message
        integer, intent(in)          :: alloc_stat
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') message
            stop
        endif
    end subroutine alloc_err
        
    ! PRETTY PROGRESS & ENDING JIFFYS

     subroutine progress( i, n )
        integer, intent(in) :: i, n
        if( .not. l_distr_exec_glob )then
            write(6,'(a1,a,t21,i3,a)',advance="no") achar(13),&
            &"Percent Complete: ", nint((real(i)/real(n))*100.0), "%"
            flush 6
            if( i >= n ) write(*,*) ''
        endif
    end subroutine progress
    
    !> \brief is for printing a progress bar
    subroutine progress_gfortran(i, imax)
        integer, intent(in) :: i, imax
        integer :: k, curr, prev
        character(len=1) :: back
        if( .not. l_distr_exec_glob )then
            back = char(8)
            if( i >= imax )then
                ! delete the bar and the percentage
                write(6,'(256a1)', advance='no') (back, k =1,59)
                ! write new bar and percentage
                write(6,'(2x,1a4,2x,1a1,256a1)', advance='no') '100%','|', ('=', k=1,50)
                write(6,'(a)') '| done.'
                flush 6            
            else
                prev = nint( 100.*real(i-1)/real(imax) )
                curr = nint( 100.*real(i)/real(imax) )
                if( curr>prev )then
                    ! delete the bar and the percentage
                    write(6,'(256a1)', advance='no') (back, k =1,(50*i/imax)+9)
                    ! write new bar and percentage
                    write(6,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') curr,'%','|', ('=', k=1,50*i/imax)
                    flush 6
                endif
            endif
        endif
    end subroutine progress_gfortran

    !> \brief  is for pretty ending
    subroutine simple_end( str, print_simple )
        character(len=*),  intent(in) :: str
        logical, optional, intent(in) :: print_simple
        logical :: pprint_simple
        pprint_simple = .true.
        if( present(print_simple) ) pprint_simple = print_simple
        if( pprint_simple )then
            write(*,'(A)') "       _______ _____ _______  _____         _______"
            write(*,'(A)') "       |______   |   |  |  | |_____] |      |______"
            write(*,'(A)') "       ______| __|__ |  |  | |       |_____ |______"      
            write(*,'(A)') " "                                                            
            write(*,'(A)') " _)_ ( _   _     ) o  _             _   _   _   o  _   _"  
            write(*,'(A)') " (_   ) ) )_)   (  ( ) ) (_( \)    )_) ) ) (_(  ( ) ) )_)"   
            write(*,'(A)') "         (_                  (\   (_         _)      (_"
            write(*,'(A)') ""
        endif
        write(*,'(A)') str
    end subroutine simple_end
  
    !> \brief  for pretty haloween ending
    subroutine haloween_end( str )
        character(len=*), intent(in) :: str
        write(*,'(A)') " #"     
        write(*,'(A)') " ##"
        write(*,'(A)') " ###"   
        write(*,'(A)') "  ####"   
        write(*,'(A)') "   #####              _______ _____ _______  _____         _______"
        write(*,'(A)') "   #######            |______   |   |  |  | |_____] |      |______"
        write(*,'(A)') "    #######           ______| __|__ |  |  | |       |_____ |______"
        write(*,'(A)') "    ########"
        write(*,'(A)') "    #########    _)_ ( _   _     ) o  _             _   _   _   o  _   _"
        write(*,'(A)') "    ##########   (_   ) ) )_)   (  ( ) ) (_( \)    )_) ) ) (_(  ( ) ) )_)"
        write(*,'(A)') "    ##########            (_                  (\   (_         _)      (_"
        write(*,'(A)') "   ###########"
        write(*,'(A)') " ##############                                                    #"          
        write(*,'(A)') "###############                                                    #"
        write(*,'(A)') " ##############                                                   ##" 
        write(*,'(A)') "   #############                                                 ##"
        write(*,'(A)') "    #############                                              ###"
        write(*,'(A)') "    ###############                                         #####"
        write(*,'(A)') "     ###############                                    ########"
        write(*,'(A)') "     ################                                ###########"
        write(*,'(A)') "     ################                              ############"
        write(*,'(A)') "     ################                           #############"
        write(*,'(A)') "    #################       #                 ###############"
        write(*,'(A)') "   ###################     ##    #           ################"
        write(*,'(A)') "  #####################   ###   ##          #################"
        write(*,'(A)') " #######################  ########         ###################"
        write(*,'(A)') "   ######################  #######        #####################"
        write(*,'(A)') "     #############################       #######################"
        write(*,'(A)') "        #########################       #####################"
        write(*,'(A)') "           ################################################"
        write(*,'(A)') "            ##############################################"
        write(*,'(A)') "              ###########################################"
        write(*,'(A)') "               #########################################"
        write(*,'(A)') "                #######################################"
        write(*,'(A)') "                 ######################################"
        write(*,'(A)') "                 ######################################"
        write(*,'(A)') "                 ###########################      ######"
        write(*,'(A)') "                 ####  ###################           ###"
        write(*,'(A)') "                 ###    ################              ##"
        write(*,'(A)') "                 ##     ##  ##########                 #"
        write(*,'(A)') "                 #      #   #### ####"
        write(*,'(A)') "                            ##    ###"
        write(*,'(A)') "                            #     ##"
        write(*,'(A)') "                                  #"
        write(*,'(A)') str
    end subroutine haloween_end
    
    ! ASSERTIONS AND SWAPS FROM NR
    
    function assert_eq_2(n1,n2,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2
        integer :: assert_eq_2
        if (n1 == n2) then
            assert_eq_2=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_2; simple_jiffys'
        end if
    end function assert_eq_2

    function assert_eq_3(n1,n2,n3,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2,n3
        integer :: assert_eq_3
        if (n1 == n2 .and. n2 == n3) then
            assert_eq_3=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_3; simple_jiffys'
        end if
    end function assert_eq_3
          
    function assert_eq_4(n1,n2,n3,n4,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2,n3,n4
        integer :: assert_eq_4
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
            assert_eq_4=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_4; simple_jiffys'
        end if
    end function assert_eq_4
    
    function assert_eq_n(nn,string)
        character(len=*), intent(in) :: string
        integer, dimension(:), intent(in) :: nn
        integer :: assert_eq_n
        if (all(nn(2:) == nn(1))) then
            assert_eq_n=nn(1)
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_n; simple_jiffys'
        end if
    end function assert_eq_n
    
    subroutine swap_i(a,b)
          integer, intent(inout) :: a,b
          integer :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_i

    subroutine swap_r(a,b)
          real, intent(inout) :: a,b
          real :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_r

    subroutine swap_rv(a,b)
          real, dimension(:), intent(inout) :: a,b
          real, dimension(size(a)) :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_rv

    subroutine swap_c(a,b)
          complex, intent(inout) :: a,b
          complex :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_c

    subroutine swap_cv(a,b)
          complex, dimension(:), intent(inout) :: a,b
          complex, dimension(size(a)) :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_cv

    subroutine swap_cm(a,b)
          complex, dimension(:,:), intent(inout)  :: a,b
          complex, dimension(size(a,1),size(a,2)) :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_cm

    subroutine masked_swap_rs(a,b,mask)
          real,         intent(inout) :: a,b
          logical(lgt), intent(in)    :: mask
          real :: swp
          if (mask) then
              swp=a
              a=b
              b=swp
          end if
    end subroutine masked_swap_rs

    subroutine masked_swap_rv(a,b,mask)
          real,         dimension(:), intent(inout) :: a,b
          logical(lgt), dimension(:), intent(in)    :: mask
          real, dimension(size(a)) :: swp
          where (mask)
              swp=a
              a=b
              b=swp
          end where
    end subroutine masked_swap_rv

    subroutine masked_swap_rm(a,b,mask)
          real,         dimension(:,:), intent(inout) :: a,b
          logical(lgt), dimension(:,:), intent(in)    :: mask
          real, dimension(size(a,1),size(a,2)) :: swp
          where (mask)
              swp=a
              a=b
              b=swp
          end where
    end subroutine masked_swap_rm

end module simple_jiffys
