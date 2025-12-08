program simple_test_units
include 'simple_lib.f08'
! core library tester modules generated with help from chatgpt
use simple_test_utils
use simple_string_tester
use simple_syslib_tester
use simple_fileio_tester
use simple_chash_tester
use simple_hash_tester
use simple_linked_list_tester
use simple_cmdline_tester
use simple_ori_tester
use simple_oris_tester
use simple_rec_list_tester
! hand-written unit tests
use simple_aff_prop,       only: test_aff_prop
use simple_ftexp_shsrch,   only: test_ftexp_shsrch
use simple_ftiter,         only: test_ftiter
use simple_image,          only: test_image
use simple_online_var,     only: test_online_var
use simple_user_interface, only: validate_ui_json
implicit none
#include "simple_local_flags.inc"
character(8)          :: datestr
character(len=STDLEN) :: folder
call seed_rnd
call date_and_time(date=datestr)
folder = './SIMPLE_TEST_UNITS_'//datestr
call simple_mkdir(folder)
call simple_chdir(folder)
! core library tests generated with help from chatgpt
call run_all_string_tests !#
call run_all_syslib_tests !#
call run_all_fileio_tests !#
call run_all_chash_tests
call run_all_hash_tests
call run_all_list_tests
call run_all_cmdline_tests
call run_all_ori_tests
call run_all_oris_tests
call run_all_rec_list_tests
call report_summary()

stop

! hand-written unit tests
write(*,*)'VALIDATING UI JSON FILE:'
call validate_ui_json
write(*,*)'PASSED UI JSON FILE TEST'
call test_online_var
call test_imghead
call test_oris(.false.)
call test_image(.false.)
call test_ftexp_shsrch
call test_ftiter
! local test functions
call test_multinomal
call test_euler_shift
call simple_test_fit_line
call test_aff_prop
call simple_chdir( "../" )
call simple_end('**** SIMPLE_UNIT_TEST NORMAL STOP ****')

contains

    subroutine test_multinomal
        integer :: i, irnd
        real :: pvec(10), prob
        call seed_rnd
        pvec(1) = 0.8
        do i=2,10
            pvec(i) = 0.2/9.
        end do
        write(logfhandle,*) 'this should be one:', sum(pvec)
        prob=0.
        do i=1,1000
            if( multinomal(pvec) == 1 ) prob = prob+1.
        end do
        prob = prob/1000.
        write(logfhandle,*) 'this should be 0.8:', prob
        pvec = 0.1
        write(logfhandle,*) 'this should be one:', sum(pvec)
        prob=0.
        do i=1,1000
            irnd = multinomal(pvec)
            if( irnd == 1 ) prob = prob+1.
        end do
        prob = prob/1000.
        write(logfhandle,*) 'this should be 0.1:', prob
        write(logfhandle,'(a)') 'SIMPLE_RND: MULTINOMAL TEST COMPLETED WITHOUT TERMINAL BUGS ;-)'
    end subroutine

    subroutine test_euler_shift
        type(ori) :: o
        integer   :: i
        real      :: euls(3), euls_shifted(3)
        logical   :: doshift
        call o%new(is_ptcl=.false.)
        do i=1,100000
            euls(1) = ran3()*800.-400.
            euls(2) = ran3()*500-250.
            euls(3) = ran3()*800.-400.
            call o%set_euler(euls)
            euls_shifted = o%get_euler()
            doshift = .false.
            if( euls_shifted(1) < 0. .or. euls_shifted(1) > 360. ) doshift = .true.
            if( euls_shifted(2) < 0. .or. euls_shifted(2) > 180. ) doshift = .true.
            if( euls_shifted(3) < 0. .or. euls_shifted(3) > 360. ) doshift = .true.
            if( doshift ) THROW_HARD('euler shifting does not work!')
        end do
    end subroutine

    subroutine simple_test_fit_line
        real    :: slope, intercept, datavec(100,2), corr, x
        integer :: i, j
        do i=1,10000
            ! generate the line
            slope = 5.*ran3()
            if( ran3() < 0.5 ) slope = -slope
            intercept = 10.*ran3()
            if( ran3() < 0.5 ) intercept = -intercept
!            write(logfhandle,*) '***********************************'
!            write(logfhandle,*) 'Slope/Intercept:', slope, intercept
            ! generate the data
            x = -1.
            do j=1,100
                datavec(j,1) = x
                datavec(j,2) = slope*datavec(j,1)+intercept
                x = x+0.02
            end do
            ! fit the data
            call fit_straight_line(100, datavec, slope, intercept, corr)
!            write(logfhandle,*) 'Fitted Slope/Intercept:', slope, intercept
!            write(logfhandle,*) 'Corr:', corr
            if( corr < 0.9999 )then
                THROW_HARD('fit_straight_line failed!')
            endif
        end do
        write(logfhandle,'(a)') 'FIT_STRAIGHT_LINE UNIT TEST COMPLETED ;-)'
    end subroutine

end program simple_test_units
