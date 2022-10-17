module simple_nlopt_tester
include 'simple_lib.f08'
implicit none

public :: exec_nlopt_test
private
#include "simple_local_flags.inc"

! module global constants
integer, parameter        :: NOPTS=6         ! nr of optimizers
integer, parameter        :: NRFUNS=19       ! nr of general test functions
integer, parameter        :: LAST2DFUN=32    ! index of last 2D function
integer, parameter        :: NTESTS=50       ! nr of independent tests
integer, parameter        :: NRDIMS=2        ! number of dimensions

! module global variables
character(len=20)         :: str_opts(NOPTS)      ! string descriptors for the NOPTS optimizers
integer, allocatable      :: success(:,:)         ! sucess counter
integer                   :: totfails(NOPTS)      ! number of total failures

contains

    subroutine exec_nlopt_test()
        write(logfhandle,*) '****optimiser_test, init'
        call test_all_optimizers
        write(logfhandle,*) '****optimiser_test, completed'
    end subroutine exec_nlopt_test

    !>  \brief  master test routine
    subroutine test_all_optimizers
        integer :: i, j, k, cnt, nfuns
        str_opts(1) = 'LN_COBYLA'
        str_opts(2) = 'LN_BOBYQA'
        str_opts(3) = 'LN_NEWUOA'
        str_opts(4) = 'LN_PRAXIS'
        str_opts(5) = 'LN_NELDERMEAD'
        str_opts(6) = 'LN_SBPLX'
        if( NRDIMS == 2 )then
            nfuns = LAST2DFUN
        else
            nfuns = NRFUNS
        endif
        allocate(success(NOPTS,nfuns))
        success  = 0
        totfails = 0
        do i=1,NOPTS
            write(logfhandle,*) ''
            write(logfhandle,'(a,1x,a)') '>>> TESTING OPTIMIZER:', str_opts(i)
            cnt = 0
            do j=1,NTESTS    ! restarts per optimizer
                do k=1,nfuns ! nr of test functions
                    cnt = cnt+1
                    call progress(cnt,nfuns*NTESTS)
                    call test_optimizer(str_opts(i), i, k, NRDIMS)
                end do
            end do
            do k=1,nfuns ! nr of test functions
                if( success(i,k) == NTESTS )then
                    write(logfhandle,'(a,1x,a,1x,i2)') str_opts(i), '    always succeeds with function:', k
                else if( success(i,k) == 0         )then
                    totfails(i) = totfails(i)+1
                    write(logfhandle,'(a,1x,a,1x,i2)') str_opts(i), '     never succeeds with function:', k
                else if( success(i,k) <  NTESTS/2 .and. success(i,k) >= 1  )then
                    write(logfhandle,'(a,1x,a,1x,i2)') str_opts(i), ' sometimes succeeds with function:', k
                else if( success(i,k) >= NTESTS/2  )then
                    write(logfhandle,'(a,1x,a,1x,i2)') str_opts(i), '     often succeeds with function:', k
                else
                    print *, 'weird success:', success(i,k)
                endif
            end do
        end do
        ! print results
        write(logfhandle,*) ''
        write(logfhandle,'(a)') '>>> TEST RESULTS'
        do i=1,NOPTS
            write(logfhandle,'(a,f4.0,a,f4.0,a)') str_opts(i)//' successes: ',&
            100.*(sum(success(i,:))/real(NTESTS*nfuns)), ' totfailures: ', 100.*totfails(i)/real(nfuns)
        end do
        deallocate(success)
    end subroutine test_all_optimizers

    !>  \brief  optimizer test routine, parameterized with respect to optimizer & testfunction
    subroutine test_optimizer( wopt, wopt_ind, wfun, ndim )
        use simple_testfuns
        use nlopt_wrap, only: nlopt_opt, nlopt_func, create, destroy
        use nlopt_enum, only: NLOPT_SUCCESS, algorithm_from_string
        character(len=20),  intent(in) :: wopt
        integer,            intent(in) :: wopt_ind, wfun, ndim !< which optimizer, which test function, dimension of problem
        procedure(testfun), pointer    :: costfun_ptr      !< pointer 2 test function
        class(*),           pointer    :: fun_self
        integer,            parameter  :: wp  = kind(0.0d0)
        real(wp),           parameter  :: TOL = 0.001_wp     ! tolerance for success
        type(nlopt_opt)                :: opt
        real                           :: gmin, range(2)
        real(wp)                       :: lims(ndim,2), x(2), lowest_cost
        integer                        :: stat
        call get_testfun(wfun, ndim, gmin, range, costfun_ptr) ! get testfun, gmin is the global min, range is limits       
        lims(:,1) = range(1)
        lims(:,2) = range(2)
        call create(opt, algorithm_from_string(trim(wopt)), 2)
        call opt%set_lower_bounds(lims(:,1))
        call opt%set_upper_bounds(lims(:,2))
        associate(f => nlopt_func(nloptf_myfunc))
            call opt%set_min_objective(f)
            call opt%set_ftol_rel(TOL)
            x(1) = lims(1,1) + ran3()*(lims(1,2)-lims(1,1))
            x(2) = lims(2,1) + ran3()*(lims(2,2)-lims(2,1))
            call opt%optimize(x, lowest_cost, stat)
        end associate
        if( abs(lowest_cost - gmin) <= TOL )then
            success(wopt_ind, wfun) = success(wopt_ind, wfun) + 1
        end if
    contains
        function nloptf_myfunc(x_in, gradient, func_data) result(f)
            real(wp), intent(in)              :: x_in(:)
            real(wp), intent(inout), optional :: gradient(:)
            class(*), intent(in),    optional :: func_data
            real(wp)                          :: f
            f = costfun_ptr(fun_self, real(x_in), 2)
        end function nloptf_myfunc
    end subroutine test_optimizer
end module simple_nlopt_tester
