program simple_test_units
include 'simple_lib.f08'
use simple_testfuns      ! use all in there
use simple_ftiter,       only: test_ftiter
use simple_ori,          only: test_ori
use simple_oris,         only: test_oris
use simple_image,        only: test_image
use simple_args,         only: test_args
use simple_online_var,   only: test_online_var
use simple_imghead,      only: test_imghead
use simple_ftexp_shsrch, only: test_ftexp_shsrch
use simple_aff_prop,     only: test_aff_prop
implicit none
#include "simple_local_flags.inc"
character(8)          :: datestr
character(len=STDLEN) :: folder
character(len=300)    :: command
call seed_rnd
call date_and_time(date=datestr)
folder = './SIMPLE_TEST_UNITS_'//datestr
call simple_mkdir(folder)
call simple_chdir(folder)
call test_args
call test_online_var
call test_imghead
call test_ori
call test_oris(.false.)
call test_image(.false.)
!call test_ftexp_shsrch
call test_ftiter
! LOCAL TESTFUNCTIONS
call test_multinomal
call test_testfuns
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

    subroutine test_testfuns
        procedure(testfun), pointer :: ptr
        integer                     :: i
        real                        :: gmin, range(2)
        logical                     :: success
        class(*), pointer           :: fun_self => null()
        success = .false.
        do i=1, 20
            call get_testfun(i, 2, gmin, range, ptr)
            select case(i)
                case(1)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(2)
                    if( abs(ptr(fun_self,[1.,1.],2)-gmin) < 1e-5 ) success = .true.
                case(3)
                    if( abs(ptr(fun_self,[-2.903534,-2.903534],2)-gmin) < 1e-5 ) success = .true.
                case(4)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(5)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(6)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(7)
                    if( abs(ptr(fun_self,[420.9687,420.9687],2)-gmin) < 1e-5 ) success = .true.
                case(8)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(9)
                    if( abs(ptr(fun_self,[1.,1.],2)-gmin) < 1e-5 ) success = .true.
                case(10)
                    if( abs(ptr(fun_self,[3.,0.5],2)-gmin) < 1e-5 ) success = .true.
                case(11)
                    if( abs(ptr(fun_self,[0.,-1.],2)-gmin) < 1e-5 ) success = .true.
                case(12)
                    if( abs(ptr(fun_self,[1.,3.],2)-gmin) < 1e-5 ) success = .true.
                case(13)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(14)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(15)
                    if( abs(ptr(fun_self,[1.,1.],2)-gmin) < 1e-5 ) success = .true.
                case(16)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(17)
                    if( abs(ptr(fun_self,[pi,pi],2)-gmin) < 1e-5 ) success = .true.
                case(18)
                    if( abs(ptr(fun_self,[512.,404.2319],2)-gmin) < 1e-5 ) success = .true.
                case(19)
                    if( abs(ptr(fun_self,[0.,0.],2)-gmin) < 1e-5 ) success = .true.
                case(20)
                    if( abs(ptr(fun_self,[0.,1.25313],2)-gmin) < 1e-5 ) success = .true.
                case DEFAULT
                    THROW_HARD('Unknown function index; test_testfuns')
            end select
            if( success )then
                cycle
            else
                write(logfhandle,*) 'testing of testfun:', i, 'failed!'
                write(logfhandle,*) 'minimum:', gmin
            endif
        end do
        write(logfhandle,'(a)') 'SIMPLE_TESTFUNS: TEST OF TEST FUNCTIONS COMPLETED ;-)'
    end subroutine

    subroutine test_euler_shift
        use simple_ori, only: ori
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
