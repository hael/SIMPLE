module simple_optimiser_tester
use simple_optimizer,   only: optimizer
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_defs         ! singleton
implicit none

public :: exec_optimiser_test
private

! module global constants
integer, parameter        :: NOPTS=5         ! nr of optimizers
integer, parameter        :: NRFUNS=19       ! nr of general test functions
integer, parameter        :: LAST2DFUN=32    ! index of last 2D function
integer, parameter        :: NTESTS=50       ! nr of independent tests
integer, parameter        :: NRESTARTS=1     ! nr of restarts per optimizer
integer, parameter        :: NRDIMS=2        ! number of dimensions
real, parameter           :: TOL=0.001       ! tolerance for success

! module global variables
character(len=8)          :: str_opts(NOPTS) ! string descriptors for the NOPTS optimizers
type(opt_factory)         :: ofac            ! the optimization factory object
type(opt_spec)            :: spec            ! the optimizer specification object
class(optimizer), pointer :: opt_ptr=>null() ! the generic optimizer object
integer, allocatable      :: success(:,:)    ! sucess counter
integer, allocatable      :: nevals(:,:)     ! number of function evaluations required
integer                   :: totfails(NOPTS) ! number of total failures
logical                   :: verbose=.false. ! verbose output or not

contains

    subroutine exec_optimiser_test( be_verbose )
        logical, optional, intent(in)    :: be_verbose
        verbose = .false.
        if( present(be_verbose) ) verbose = be_verbose
        write(*,*) '****optimiser_test, init'
        call test_bound_routine
        call test_all_optimizers
        write(*,*) '****optimiser_test, completed'
    end subroutine exec_optimiser_test

    subroutine test_bound_routine
        real           :: lims(2,2), vec(2)
        integer        :: i
        if( verbose ) write(*,*) 'testing the bound routine'
        lims = 0.
        lims(1,2) = 360.
        lims(2,2) = 180.
        call spec%specify('simplex', 2, limits=lims, cyclic=[.true.,.true.])
        vec(1) = 361.
        vec(2) = 181.
        call check_and_correct_vec(spec, vec)
        if( (abs(vec(1)-1.)+abs(vec(2)-1.))/2. < 1e-6 )then
        else
            stop '****optimiser_tester bound test1 failed!'
        endif
        vec(1) = -1.
        vec(2) = -1.
        call check_and_correct_vec(spec, vec)
        if( (abs(vec(1)-359.)+abs(vec(2)-179.))/2. < 1e-6 )then
        else
            stop '****optimiser_tester bound test2 failed!'
        endif
        vec(1) = 361.
        vec(2) = -1.
        call check_and_correct_vec(spec, vec)
        if( (abs(vec(1)-1.)+abs(vec(2)-179.))/2. < 1e-6 )then
        else
            stop '****optimiser_tester bound test3 failed!'
        endif
        vec(1) = -1
        vec(2) = 181.
        call check_and_correct_vec(spec, vec)
        if( (abs(vec(1)-359.)+abs(vec(2)-1.))/2. < 1e-6 )then
        else
            stop '****optimiser_tester bound test4 failed!'
        endif
        call spec%specify('simplex', 2, limits=lims)
        do i=1,1000
            vec(1) = 361.
            vec(2) = 181.
            call check_and_correct_vec(spec, vec)
            if( vec(1) >= 180. .and. vec(1) <= 360. .and. vec(2) >= 90. .and. vec(2) <= 180. )then
            else
                stop '****optimiser_tester bound test5 failed!'
            endif 
            vec(1) = -1.
            vec(2) = -1.
            call check_and_correct_vec(spec, vec)
            if( vec(1) >= 0. .and. vec(1) <= 180. .and. vec(2) >= 0. .and. vec(2) <= 90. )then
            else
                stop '****optimiser_testerbound test6 failed!'
            endif
        end do
        
        contains
        
            !> \brief  check the vector with respect to the limits
            subroutine check_and_correct_vec( spec, vec )
                use simple_opt_spec, only: opt_spec
                use simple_rnd,      only: ran3
                class(opt_spec), intent(in) :: spec           !< specification
                real, intent(inout)         :: vec(spec%ndim) !< solution vector
                integer                     :: j
                real                        :: cen, halfwidth, ulim
                do j=1,spec%ndim
                    if( spec%cyclic(j) ) then
                        do while(vec(j) < spec%limits(j,1)) 
                            vec(j) = vec(j)+spec%limits(j,2)
                        end do
                        do while(vec(j) > spec%limits(j,2))
                            vec(j) = vec(j)-spec%limits(j,2)
                        end do
                    else
                        halfwidth = (spec%limits(j,2)-spec%limits(j,1))/2.
                        cen       = spec%limits(j,1)+halfwidth
                        ! generate a point by controlled randomization
                        if( vec(j) < spec%limits(j,1) ) then
                            ! random point in the lower half of the interval
                             vec(j) = spec%limits(j,1)+ran3()*(cen-spec%limits(j,1))
                        else if( vec(j) > spec%limits(j,2) ) then
                            ! random point in the upper half of the interval
                            ulim = spec%limits(j,1)+cen
                            vec(j) = ulim+ran3()*(spec%limits(j,2)-ulim)
                        endif
                    endif
                end do
            end subroutine
        
    end subroutine test_bound_routine
    
    !>  \brief  master test routine
    subroutine test_all_optimizers
        use simple_jiffys, only: progress
        use simple_stat,   only: moment
        integer :: i, j, k, mval, cnt, nfuns
        real    :: maxeval, neval
        logical :: err
        str_opts(1) = 'powell'
        str_opts(2) = 'simplex'
        str_opts(3) = 'oasis'
        str_opts(4) = 'pso'
        str_opts(5) = 'de'
        if( NRDIMS == 2 )then
            nfuns = LAST2DFUN
        else
            nfuns = NRFUNS
        endif
        allocate(success(NOPTS,nfuns),nevals(NOPTS,nfuns))
        success  = 0
        nevals   = 0
        totfails = 0
        do i=1,NOPTS
            write(*,*) ''
            write(*,'(a,1x,a)') '>>> TESTING OPTIMIZER:', str_opts(i)
            cnt = 0
            do j=1,NTESTS    ! restarts per optimizer
                do k=1,nfuns ! nr of test functions
                    cnt = cnt+1
                    call progress(cnt,nfuns*NTESTS)
                    call test_optimizer(i,k,NRDIMS)
                end do
            end do
            do k=1,nfuns ! nr of test functions
                if( success(i,k) == NTESTS )then
                    write(*,'(a,1x,a,1x,i2)') str_opts(i), '    always succeeds with function:', k
                else if( success(i,k) == 0         )then
                    totfails(i) = totfails(i)+1
                    write(*,'(a,1x,a,1x,i2)') str_opts(i), '     never succeeds with function:', k
                else if( success(i,k) <  NTESTS/2 .and. success(i,k) >= 1  )then
                    write(*,'(a,1x,a,1x,i2)') str_opts(i), ' sometimes succeeds with function:', k
                else if( success(i,k) >= NTESTS/2  )then
                    write(*,'(a,1x,a,1x,i2)') str_opts(i), '     often succeeds with function:', k
                else 
                    print *, 'weird success:', success(i,k) 
                endif 
            end do
        end do
        ! normalize the nevals matrix
        maxeval = 0.
        do i=1,NOPTS
            neval = real(sum(nevals(i,:)))/real(NTESTS*nfuns)
            if( neval > maxeval )then
                maxeval = neval
            endif
        end do
        ! print results    
        write(*,*) ''
        write(*,'(a)') '>>> TEST RESULTS'
        do i=1,NOPTS
            write(*,'(a,f4.0,a,f4.0,a,f4.0)') str_opts(i)//' successes: ',&
            100.*(sum(success(i,:))/real(NTESTS*nfuns)), ' totfailures: ', 100.*totfails(i)/real(nfuns),&
            ' nevals: ', 100.*real(sum(nevals(i,:)))/(real(NTESTS*nfuns)*maxeval)
        end do
        deallocate(success,nevals)
        ! Finally, we thest the brute-force one
        call simple_test_bforce_opt
    end subroutine test_all_optimizers

    !>  \brief  optimizer test routine, parameterized with respect to optimizer & testfunction
    subroutine test_optimizer( wopt, wfun, ndim )
        use simple_rnd,  only: ran3
        use simple_testfuns
        integer, intent(in)         :: wopt, wfun, ndim !< which optimizer, which test function, dimension of problem
        procedure(testfun), pointer :: costfun_ptr      !< pointer 2 test function
        real :: h(ndim), gmin, range(2), lowest_cost, lims(ndim,2), rtol
        h = 1.      
        call get_testfun(wfun, ndim, gmin, range, costfun_ptr) ! get testfun, gmin is the global min, range is limits
        lims(:,1) = range(1)
        lims(:,2) = range(2)
        call spec%specify(str_opts(wopt),ndim,limits=lims,nrestarts=NRESTARTS) ! make optimizer spec
        call spec%set_costfun(costfun_ptr)                                     ! set pointer to costfun
        call ofac%new(spec, opt_ptr)                                           ! generate optimizer object with the factory
        ! generate starting point
        spec%x(1) = spec%limits(1,1)+ran3()*(spec%limits(1,2)-spec%limits(1,1))
        spec%x(2) = spec%limits(2,1)+ran3()*(spec%limits(2,2)-spec%limits(2,1))
        call opt_ptr%minimize(spec, lowest_cost)                         ! minimize the test function
        rtol=2.0*abs(gmin-lowest_cost)/(abs(gmin)+abs(lowest_cost)+TINY) ! relative tolerance
        if( sqrt((lowest_cost-gmin)**2.) <= TOL )then
            success(wopt,wfun) = success(wopt,wfun)+1
            nevals(wopt,wfun)  = nevals(wopt,wfun)+spec%nevals
        end if    
    end subroutine test_optimizer

    subroutine simple_test_bforce_opt
        use simple_bforce_opt, only: bforce_opt
        use simple_opt_spec,   only: opt_spec
        use simple_testfuns
        implicit none
        type(opt_spec)              :: spec
        type(bforce_opt)            :: bforce
        real                        :: lowest_cost, limits(2,2), range(2), gmin, dist
        procedure(testfun), pointer :: pfun
        call get_testfun(7, 2, gmin, range, pfun)
        limits(1,1) = range(1)
        limits(1,2) = range(2)
        limits(2,1) = range(1)
        limits(2,2) = range(2)
        call spec%specify('bforce', 2, limits=limits, stepsz=[0.01,0.01])
        call spec%set_costfun(pfun)
        call bforce%new(spec)
        call bforce%minimize(spec, lowest_cost)
        dist = sqrt(sum(spec%x**2.0))
        if( dist < 1e-4 .and. abs(lowest_cost) < 1e-9 )then
            ! test passed
        else
            write(*,*) 'dist from global opt (0,0): ', dist
            write(*,*) 'cost obtained (lowest=0):   ', lowest_cost 
            stop '****optimiser_test FAILURE bforce_opt'
        endif
    end subroutine simple_test_bforce_opt

end module simple_optimiser_tester
