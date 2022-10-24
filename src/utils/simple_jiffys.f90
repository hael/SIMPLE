! jiffy =  the time it takes light to travel one centimeter in vacuum
module simple_jiffys
use simple_defs
use simple_timer
implicit none

interface swap
    module procedure swap_i,swap_r,swap_rv,swap_c, swap_cv,swap_cm,&
    &masked_swap_rs,masked_swap_rv,masked_swap_rm
end interface swap

contains

    ! PRETTY PROGRESS & ENDING JIFFYS

    subroutine progress( i, n )
        integer, intent(in) :: i, n
        if( .not. l_distr_exec_glob .and. L_VERBOSE_GLOB )then
            write(logfhandle,'(a1,a,t21,i3,a)',advance="no") achar(13),&
                &"Percent Complete: ", nint((real(i)/real(n))*100.0), "%"
            flush (OUTPUT_UNIT)
            if( i >= n ) write(logfhandle,*) ''
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

    !> \brief  is for pretty elapsed time printing
    subroutine simple_print_timer( elapsed )
        real(timer_int_kind), intent(in) :: elapsed
        write(logfhandle,'(A,F9.2,A)') "**** Execution time : ", elapsed, " seconds"
    end subroutine simple_print_timer

    subroutine simple_print_git_version( git_hash )
        character(len=*), intent(in) :: git_hash
        write(logfhandle,'(A,A)') "**** SIMPLE Git Commit: ", git_hash
    end subroutine simple_print_git_version

    !> \brief  is for pretty ending
    subroutine simple_end( str, print_simple )
        use simple_syslib,  only: get_process_id, del_file
        use simple_strings, only: int2str
        character(len=*),  intent(in) :: str
        logical, optional, intent(in) :: print_simple
        character(len=:), allocatable :: pid_file
        logical :: pprint_simple
        integer :: pid
        ! delete file indicating active process
        pid = get_process_id()
        allocate(pid_file, source='.'//int2str(pid)//'.simple.pid')
        call del_file(pid_file)
        ! pretty ending
        pprint_simple = .true.
        if( present(print_simple) ) pprint_simple = print_simple
        if( pprint_simple .and. L_VERBOSE_GLOB )then
            write(logfhandle,'(A)') "       _______ _____ _______  _____         _______"
            write(logfhandle,'(A)') "       |______   |   |  |  | |_____] |      |______"
            write(logfhandle,'(A)') "       ______| __|__ |  |  | |       |_____ |______"
            write(logfhandle,'(A)') " "
            write(logfhandle,'(A)') " _)_ ( _   _     ) o  _             _   _   _   o  _   _"
            write(logfhandle,'(A)') " (_   ) ) )_)   (  ( ) ) (_( \)    )_) ) ) (_(  ( ) ) )_)"
            write(logfhandle,'(A)') "         (_                  (\   (_         _)      (_"
            write(logfhandle,'(A)') ""
        endif
        write(logfhandle,'(A)') str
    end subroutine simple_end

    !> \brief  for pretty haloween ending
    subroutine haloween_end( str )
        character(len=*), intent(in) :: str
        write(logfhandle,'(A)') " #"
        write(logfhandle,'(A)') " ##"
        write(logfhandle,'(A)') " ###"
        write(logfhandle,'(A)') "  ####"
        write(logfhandle,'(A)') "   #####              _______ _____ _______  _____         _______"
        write(logfhandle,'(A)') "   #######            |______   |   |  |  | |_____] |      |______"
        write(logfhandle,'(A)') "    #######           ______| __|__ |  |  | |       |_____ |______  v3.0"
        write(logfhandle,'(A)') "    ########"
        write(logfhandle,'(A)') "    #########    _)_ ( _   _     ) o  _             _   _   _   o  _   _"
        write(logfhandle,'(A)') "    ##########   (_   ) ) )_)   (  ( ) ) (_( \)    )_) ) ) (_(  ( ) ) )_)"
        write(logfhandle,'(A)') "    ##########            (_                  (\   (_         _)      (_"
        write(logfhandle,'(A)') "   ###########"
        write(logfhandle,'(A)') " ##############                                                    #"
        write(logfhandle,'(A)') "###############                                                    #"
        write(logfhandle,'(A)') " ##############                                                   ##"
        write(logfhandle,'(A)') "   #############                                                 ##"
        write(logfhandle,'(A)') "    #############                                              ###"
        write(logfhandle,'(A)') "    ###############                                         #####"
        write(logfhandle,'(A)') "     ###############                                    ########"
        write(logfhandle,'(A)') "     ################                                ###########"
        write(logfhandle,'(A)') "     ################                              ############"
        write(logfhandle,'(A)') "     ################                           #############"
        write(logfhandle,'(A)') "    #################       #                 ###############"
        write(logfhandle,'(A)') "   ###################     ##    #           ################"
        write(logfhandle,'(A)') "  #####################   ###   ##          #################"
        write(logfhandle,'(A)') " #######################  ########         ###################"
        write(logfhandle,'(A)') "   ######################  #######        #####################"
        write(logfhandle,'(A)') "     #############################       #######################"
        write(logfhandle,'(A)') "        #########################       #####################"
        write(logfhandle,'(A)') "           ################################################"
        write(logfhandle,'(A)') "            ##############################################"
        write(logfhandle,'(A)') "              ###########################################"
        write(logfhandle,'(A)') "               #########################################"
        write(logfhandle,'(A)') "                #######################################"
        write(logfhandle,'(A)') "                 ######################################"
        write(logfhandle,'(A)') "                 ######################################"
        write(logfhandle,'(A)') "                 ###########################      ######"
        write(logfhandle,'(A)') "                 ####  ###################           ###"
        write(logfhandle,'(A)') "                 ###    ################              ##"
        write(logfhandle,'(A)') "                 ##     ##  ##########                 #"
        write(logfhandle,'(A)') "                 #      #   #### ####"
        write(logfhandle,'(A)') "                            ##    ###"
        write(logfhandle,'(A)') "                            #     ##"
        write(logfhandle,'(A)') "                                  #"
        write(logfhandle,'(A)') str
    end subroutine haloween_end

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
        if(mask)then
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
