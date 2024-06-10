program simple_test_steep_grad
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
type(parameters)   :: p
type(cmdline)      :: cline

if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_steep_grad smpd=xx nthr=yy'
    write(logfhandle,'(a)') 'Example: simple_test_steep_grad smpd=1. nthr=4'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd', 1)
call cline%checkvar('nthr', 2)
call cline%check
call p%new(cline)


contains

    subroutine func( x, f, grad )
        real, intent(in)  :: x(2)
        real, intent(out) :: f
        real, intent(out) :: grad(2)
        f       = (x(1) - 2.)**2 - (x(1) - 2.) * (x(2) - 2.) + 0.5*(x(2) - 2.)**2
        grad(1) = 2. * (x(1) - 2.) - (x(2) - 2.)
        grad(2) =     -(x(1) - 2.) + (x(2) - 2.)
    end subroutine func

    function polymod( q0, qp0, lamc, qc, blow, bhigh, lamm, qm ) result(lplus)
        real,           intent(in) :: q0, qp0, lamc, qc, blow, bhigh
        real, optional, intent(in) :: lamm, qm
        real :: lleft, lright, lplus, a(2,2), b(2), c(2), denom
        lleft  = lamc * blow
        lright = lamc * bhigh
        if( present(lamm) .and. present(qm) )then
            lplus = -(qp0 * lamc * lamc) / (2. * (qc - q0 - qp0*lamc))
        else
            a(1,:) = [lamc**2, lamc**3]
            a(2,:) = [lamm**2, lamm**3]
            b      = [qc, qm] - [q0 + qp0*lamc, q0 + qp0*lamm]
            denom  = (a(1,1)*a(2,2) - a(2,1)*a(1,2))
            c(1)   = (b(1) * a(2,2) - b(2) *a(1,2))/denom
            c(2)   = (b(2) * a(1,1) - b(1) *a(2,1))/denom
            lplus  = (-c(1) + sqrt(c(1) * c(1) - 3. * c(2) * qp0)) / (3. * c(2))
        endif
        if( lplus < lleft  ) lplus = lleft
        if( lplus > lright ) lplus = lright
    end function polymod

end program simple_test_steep_grad
