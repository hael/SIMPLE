! jiffy =  the time it takes light to travel one centimeter in vacuum
module simple_jiffys
use simple_defs ! singleton
implicit none

interface assert_eq
    module procedure assert_eq_2,assert_eq_3,assert_eq_4,assert_eq_n
end interface assert_eq

interface swap
    module procedure swap_i,swap_r,swap_rv,swap_c, swap_cv,swap_cm,&
    &masked_swap_rs,masked_swap_rv,masked_swap_rm
end interface swap

contains

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
            write(*,'(A)') "       ______| __|__ |  |  | |       |_____ |______  v3.0"
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
        write(*,'(A)') "    #######           ______| __|__ |  |  | |       |_____ |______  v3.0"
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
