program simple_test_nlopt_lp_spectrum
include 'simple_lib.f08'
use simple_cmdline,     only: cmdline
use simple_image,       only: image
use simple_parameters,  only: parameters
use nlopt_wrap,         only: nlopt_opt, nlopt_func, create, destroy
use nlopt_enum,         only: NLOPT_SUCCESS, algorithm_from_string
implicit none

type(parameters)              :: p
type(cmdline)                 :: cline
character(len=:), allocatable :: cmd
real,             allocatable :: spec(:), y(:)
integer,          parameter   :: wp  = kind(0.0d0)
real(wp),         parameter   :: TOL = 0.001_wp     ! tolerance for success
logical                       :: mrc_exists
integer                       :: ifoo, rc, stat, i, lp_ind
real(wp)                      :: x(2), lowest_cost
type(image)                   :: img
type(nlopt_opt)               :: opt
real                          :: sum_spec
if( command_argument_count() < 4 )then
    write(logfhandle,'(a)') 'Usage: simple_test_nlopt_lp_spectrum smpd=xx nthr=yy vol1=volume.mrc mskdiam=zz'
    write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
    inquire(file="1JYX.mrc", exist=mrc_exists)
    if( .not. mrc_exists )then
        write(*, *) 'Downloading the example dataset...'
        cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
        write(*, *) 'Converting .pdb to .mrc...'
        cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
        call execute_command_line(cmd, exitstat=rc)
        cmd = 'rm 1JYX.pdb'
        call execute_command_line(cmd, exitstat=rc)
    endif
    call cline%set('smpd'   , 1.)
    call cline%set('nthr'   , 16.)
    call cline%set('vol1'   , '1JYX.mrc')
    call cline%set('mskdiam', 180.)
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd',    1)
call cline%checkvar('nthr',    2)
call cline%checkvar('vol1',    3)
call cline%checkvar('mskdiam', 4)
call cline%check
call p%new(cline)
call find_ldim_nptcls(p%vols(1), p%ldim, ifoo)
call img%new( p%ldim, p%smpd)
call img%read(p%vols(1))
call img%spectrum('power', spec)
! fitting a decaying function to spec y = A*exp(-alpha*spec)
call create(opt, algorithm_from_string('LN_COBYLA'), 2)
associate(f => nlopt_func(nloptf_myfunc))
    call opt%set_min_objective(f)
    call opt%set_ftol_rel(TOL)
    x(1) = maxval(spec)
    x(2) = 1.
    call opt%optimize(x, lowest_cost, stat)
end associate        
! best lp is at the fourier index that contributes less than 0.1% of the energy
allocate(y(size(spec)), source=0.)
sum_spec = 0.
lp_ind   = -1
do i = 1, size(spec)
    y(i)     = x(1)*exp(-x(2)*i)
    sum_spec = sum_spec + y(i)
    if( y(i)/sum_spec < 0.001 .and. lp_ind == -1 ) lp_ind = i
enddo
print *, 'lp_ind = ', lp_ind, '; cut-off lp = ', calc_lowpass_lim(lp_ind, p%ldim(1), p%smpd)

contains
    function nloptf_myfunc(x_in, gradient, func_data) result(f)
        real(wp), intent(in)              :: x_in(:)
        real(wp), intent(inout), optional :: gradient(:)
        class(*), intent(in),    optional :: func_data
        real(wp)                          :: f
        integer                           :: N, i
        N = size(spec)
        f = 0.
        do i = 1, N
            f = f + (spec(i) - x_in(1)*exp(-x_in(2)*i))**2
        enddo
    end function nloptf_myfunc
end program simple_test_nlopt_lp_spectrum