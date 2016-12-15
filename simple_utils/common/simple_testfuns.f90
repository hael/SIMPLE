!==Class simple_testfuns
!
! the simple_testfuns singleton provides 20 mathematical test functions for evaluating unconstrained optimization
! procedures. 1-9 are generalized test functions (of arbitrary dimension) whereas functions 10-20 are two-dimensional
! test functions. The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2014-01-07.
! 
!==Changes are documented below
!
! RELEVANT INFO:
! http://en.wikipedia.org/wiki/Test_functions_for_optimization
! http://www.sfu.ca/~ssurjano/optimization.html
!
module simple_testfuns
use simple_defs
implicit none
    
!>  \brief  defines the test function interface
abstract interface
    function testfun( vec, D ) result( cost )
        integer, intent(in) :: D
        real, intent(in)    :: vec(D)
        real                :: cost
    end function 
end interface

contains

    !>  \brief  for getting global minimum value, search domain (range), and pointer to testfun i
    function get_testfun( i, d, gmin, range ) result( ptr )
        integer, intent(in)         :: i, d
        real, intent(out)           :: gmin, range(2)
        procedure(testfun), pointer :: ptr
        if( d < 2 ) stop 'none of these functions are defined for d < 2; get_testfun; simple_testfuns'
        if( i < 20 ) then
        else
            if( d /= 2 ) stop 'testfuns 20-32 are intended for d == 2; get_testfun; simple_testfuns'
        endif
        gmin     = 0.
        range(2) = 10.
        range(1) = -range(2)
        select case(i)
        
            ! MANY LOCAL MINIMA (MULTIMODAL)
        
            case(1)
                !>  \brief  Ackley's function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-32,32]
                ptr => testfun1
                range(2) = 32.
                range(1) = -range(2)
            case(2)
                !>  \brief  Griewank's function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-600,600]
                ptr => testfun2
                range(2) = 600.
                range(1) = -range(2)
            case(3)
                !>  \brief  Levy function
                !!          global minimum: f(1,...,1) = 0
                !!          search domain: [-10,10]
                ptr => testfun3
            case(4)
                !>  \brief  Rastrigin's function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-5.12,5.12]
                ptr => testfun4
                range(2) = 5.12
                range(1) = -range(2)
            case(5)
                !>  \brief  Schwefel's function
                !!          global minimum: f(420.9687,...,420.9687) = 0
                !!          search domain: [-500,500]
                !!          characteristics: highly complex, with many local minima
                ptr => testfun5
                range(2) = 500.
                range(1) = -range(2)
            case(6)
                !>  \brief  Salomon's function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-100,100]
                !!          characteristics: highly multimodal
                ptr => testfun6
                range(2) = 100.
                range(1) = -range(2)
                
            ! BOWL-SHAPED (UNIMODAL)
                
            case(7)
                !>  \brief  Sphere function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-5.12,5.12]
                ptr => testfun7
                range(2) = 5.12
                range(1) = -range(2)
            case(8)
                !>  \brief  SumSquares function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-10,10]
                ptr => testfun8
            case(9)
                !>  \brief  SumPowers
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-1,1]
                ptr => testfun9
                range(2) = 1.
                range(1) = -range(2)
            case(10)
                !>  \brief  1st Perm function
                !!          global minimum: f(1,1/2,1/3,...,1/d) = 0
                !!          search domain: [-d,d]
                ptr => testfun10
                range(2) = real(d)
                range(1) = -range(2)
                
            ! WAVY BOWLS
            
            case(11)
                !>  \brief  2nd Perm function
                !!          global minimum: f(1,2,3,...,d) = 0
                !!          search domain: [-d,d]
                !!          characteristics: wavy bowl
                ptr => testfun11
                range(2) = real(d)
                range(1) = -range(2)
            case(12)
                !>  \brief  Styblinski-Tang function
                !!          global minimum: f(-2.903534,...,-2.903534) = -39.16599*d
                !!          search domain: [-5,5]
                !!          characteristics: wavy bowl
                ptr => testfun12
                range(2) = 5
                range(1) = -range(2)
                gmin     = -39.16599*real(d)
            case(13)
                !>  \brief  Whitley's function
                !!          global minimum: f(1,...,1) = 0
                !!          search domain: [-10.24,10.24]
                !!          characteristics: wavy bowl
                ptr => testfun13
                range(2) = 10.24
                range(1) = -range(2)
                
            ! PLATE-SHAPED
                
            case(14)
                !>  \brief  Zakharov function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-5,10]
                ptr => testfun14
                range(2) = 10.
                range(1) = -5.
                
            ! VALLEY-SHAPED
                
            case(15)
                !>  \brief  Dixon-Price function
                !!          global minimum: f(x(1),...,x(d)) = 0 at x(i) = 2**(-(2**i-2)/(2**i))
                !!          search domain: [-10,10]
                ptr => testfun15
                range(2) = 10.
                range(1) = -range(2)
            case(16)
                !>  \brief  Rosenbrock's (banana) function
                !!          global minimum: f(1,...,1) = 0
                !!          search domain: [-5,10]
                !!          can be restricted to: [-2.048,2.048]
                ptr => testfun16
                range(2) = 10.
                range(1) = -5.
                
            ! OTHER
            
            case(17)
                !>  \brief  Pinter's function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-100,100]
                ptr => testfun17
                range(2) = 100.
                range(1) = -range(2)
            case(18)
                !>  \brief  MultiMod function
                !!          global minimum: f(0,...,0) = 0
                !!          search domain: [-10,10]
                ptr => testfun18
                range(2) = 10.
                range(1) = -range(2)
            case(19)
                !>  \brief  Step function
                !!          global minimum: f(0.5,...,0.5) = 0
                !!          search domain: [-100,100]
                !!          The presence of many flat plateus and steep ridges presents 
                !!          difficulties for algorithms based on gradient information
                ptr => testfun19
                range(2) = 100.
                range(1) = -range(2)
                
            ! TWO-DIMENSIONAL TEST FUNCTIONS
            
            ! MANY LOCAL MINIMA (MULTIMODAL)
            
             case(20)
                !>  \brief  Drop-wave function
                !!          global minimum: f(0,0) = -1
                !!          search domain: [-5.12,5.12]
                ptr => testfun20
                range(2) = 5.12
                range(1) = -range(2)
                gmin = -1.
            case(21)
                !>  \brief Eggholder function
                !!         global minimum: f(512,404.2319) = -959.6407 
                !!         search domain: [-512,512]
                ptr => testfun21
                range(2) = 512.
                range(1) = -range(2)
                gmin     = -959.6407
            case(22)
                !>  \brief  Holder table function
                !!          global minimum: f(8.05502,9.66459)   = -19.2085
                !!          global minimum: f(8.05502,-9.66459)  = -19.2085
                !!          global minimum: f(-8.05502,9.66459)  = -19.2085
                !!          global minimum: f(-8.05502,-9.66459) = -19.2085
                !!          search domain: [-10,10]
                ptr => testfun22
                gmin = -19.2085
            case(23)
                !>  \brief  Levy function N.13
                !!          global minimum: f(1,1) = 0
                !!          search domain: [-10,10]
                ptr => testfun23
            case(24)
                !>  \brief Schaffer function N.2
                !!         global minimum: f(0,0) = 0 
                !!         search domain: [-100,100]
                ptr => testfun24
                range(2) = 100.
                range(1) = -range(2)
            case(25)
                !>  \brief Schaffer function N.4
                !!         global minimum: f(0,1.25313) = 0.292579 
                !!         search domain: [-100,100]
                ptr => testfun25
                range(2) = 100.
                range(1) = -range(2)
                gmin     = 0.292579
            case(26)
                !>  \brief Shubert function
                !!         global minimum: f(unknown) = -186.7309
                !!         search domain: [-5.12,5.12]
                ptr => testfun26
                range(2) = 5.12
                range(1) = -range(2)
                gmin = -186.7309

            ! PLATE-SHAPED
            
            case(27)
                !>  \brief  Booth's function
                !!          global minimum: f(1,3) = 0
                !!          search domain: [-10,10]
                ptr => testfun27
                range(2) = 10.
                range(1) = -range(2)
            case(28)
                !>  \brief  Matya's function
                !!          global minimum: f(0,0) = 0
                !!          search domain: [-10,10]
                ptr => testfun28

            ! VALLEY-SHAPED
            
            case(29)
                !>  \brief  three-hump camel function
                !!          global minimum: f(0,0) = 0
                !!          search domain: [-5,5]
                ptr => testfun29
                range(2) = 5.
                range(1) = -range(2)
            case(30)
                !>  \brief  six-hump camel function
                !!          global minimum: f(0.0898,-0.7126) = -1.0316
                !!          global minimum: f(-0.0898,0.7126) = -1.0316
                !!          search domain: [-3,-3]
                ptr => testfun30
                range(2) = 3.
                range(1) = -range(2)
                gmin     = -1.0316

            ! OTHER
                
            case(31)
                !>  \brief  Beale's function
                !!          global minimum: f(3,0.5) = 0
                !!          search domain: [-4.5,4.5]
                !!          characteristics: plateued bowl
                ptr => testfun31
                range(2) = 4.5
                range(1) = -range(2)
            case(32)
                !>  \brief Goldstein-Price function
                !!         global minimum: f(0,-1) = 3
                !!         search domain: [-2,2]
                !!         characteristics: plateued humps
                ptr => testfun32
                range(2) = 2.
                range(1) = -range(2)
                gmin     = 3.
            case DEFAULT
                stop 'Unknown function index; get_testfun; simple_testfuns'
        end select
    end function
    
    ! GENERALIZED TEST FUNCTIONS
    
    ! MANY LOCAL MINIMA (MULTIMODAL)
    
    !>  \brief  Ackley's function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-32,32]
    function testfun1( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real :: r, s
        r = 0.
        s = 0.
        do i=1,d
            r = r+x(i)**2.
            s = s+cos(twopi*x(i))
        end do
        r = 20.+exp(1.)-20.*exp(-0.2*sqrt(r/real(d)))-exp(s/real(d))
    end function
    
    !>  \brief  Griewank's function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-600,600]
    function testfun2( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real :: r, s
        r = 0.
        s = 1.
        do i=1,d
            r = r+(x(i)**2.)/4000.
            s = s*cos(x(i)/sqrt(real(i)))
        end do
        r = r-s+1
    end function
    
    !>  \brief  Levy function
    !!          global minimum: f(1,...,1) = 0
    !!          search domain: [-10,10]
    function testfun3( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        real                :: w(d)
        integer :: i
        real    :: r, var, var2, last
        w     = (x-1.)/4.+1
        var   = sin(pi*w(1))
        r     = var*var
        var   = sin(twopi*w(d))
        last  = ((w(d)-1.)**2.)*(1.+var*var)
        do i=1,d-1
            var = sin(pi*w(i)+1.)
            var2 = 1.+10.*var*var
            r = r+(w(i)-1.)**2.*var2
        end do
        r = r+last
    end function
    
    !>  \brief  Rastrigin's function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-5.12,5.12]
    function testfun4( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real :: r
        r = 10.*real(d)
        do i=1,d
            r = r+x(i)**2.-10.*cos(twopi*x(i))
        end do
    end function
    
    !>  \brief  Schwefel's function
    !!          global minimum: f(420.9687,...,420.9687) = 0
    !!          search domain: [-500,500]
    !!          characteristics: highly complex, with many local minima
    function testfun5( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r
        r = 0.
        do i=1,d
            r = r+x(i)*sin(sqrt(abs(x(i))))
        end do
        r = 418.9829*real(d)-r
    end function
    
    !>  \brief  Salomon's function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-100,100]
    !!          characteristics: highly multimodal
    function testfun6( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r
        r = 0.
        do i=1,d
            r = r+x(i)*x(i)
        end do
        r = 1.+0.1*sqrt(r)-cos(twopi*sqrt(r))
    end function
    
    ! BOWL-SHAPED (UNIMODAL)
    
    !>  \brief  Sphere function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-5.12,5.12]
    function testfun7( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real :: r
        r = 0.
        do i=1,d
            r = r+x(i)**2.
        end do
    end function
    
    !>  \brief  SumSquares function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-10,10]
    function testfun8( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r
        r = 0.
        do i=1,d
            r = r+real(i)*x(i)*x(i)
        end do
    end function

    !>  \brief  SumPowers function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-1,1]
    function testfun9( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r
        r = 0.
        do i=1,d
            r = r+abs(x(i))**real(i+1)
        end do
    end function
    
    !>  \brief  1st Perm function
    !!          global minimum: f(1,1/2,1/3,...,1/d) = 0
    !!          search domain: [-d,d]
    function testfun10( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i, j
        real    :: r, ii, jj, r_inner
        r = 0.
        do i=1,d
            ii = real(i)
            r_inner = 0.
            do j=1,d 
                jj = real(j)
                r_inner = r_inner+(jj+10.)*(x(j)**ii-(1./jj)**ii)   
            end do
            r = r+r_inner**2.
        end do
    end function
    
    ! WAVY BOWLS
    
    !>  \brief  2nd Perm function
    !!          global minimum: f(1,2,3,...,d) = 0
    !!          search domain: [-d,d]
    !!          characteristics: wavy bowl
    function testfun11( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i, j
        real    :: r, ii, jj, r_inner
        r = 0.
        do i=1,d
            ii = real(i)
            r_inner = 0.
            do j=1,d 
                jj = real(j)
                r_inner = r_inner+(jj**ii+0.5)*((x(j)/jj)**ii-1.)
            end do
            r = r+r_inner**2.
        end do
    end function
    
    !>  \brief  Styblinski-Tang function
    !!          global minimum: f(-2.903534,...,-2.903534) = -39.16599*d
    !!          search domain: [-5,5]
    !!          characteristics: wavy bowl
    function testfun12( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real :: r
        r = 0.
        do i=1,d
            r = r+x(i)**4.-16.*x(i)**2.+5.*x(i)
        end do
        r = r/2.
    end function
    
    !>  \brief  Whitley's function
    !!          global minimum: f(1,...,1) = 0
    !!          search domain: [-10.24,10.24]
    !!          characteristics: wavy bowl
    function testfun13( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i,j
        real    :: r, yij
        r = 0.
        do j=1,d
            do i=1,d
                yij = 100.*(x(i)**2.-x(j))**2.+(1.-x(i))**2.
                r = r+(yij*yij)/4000.-cos(yij)+1
            end do
        end do
    end function
    
    ! PLATE-SHAPED (UNIMODAL)
    
    !>  \brief  Zakharov function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-5,10]
    function testfun14( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r, ii, s1, s2
        s1 = 0.
        s2 = 0.
        do i=1,d
            ii = real(i)
            s1 = s1+x(i)*x(i)
            s2 = s2+0.5*ii*x(i)
        end do
        r = s1+s2**2.+s2**4
    end function
    
    ! VALLEY-SHAPED 
    
    !>  \brief  Dixon-Price function
    !!          global minimum: f(x(1),...,x(d)) = 0 at x(i) = 2**(-(2**i-2)/(2**i))
    !!          search domain: [-10,10]
    function testfun15( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r
        r = (x(1)-1.)**2.
        do i=2,d
            r = r+real(i)*(2.*x(i)**2.-x(i-1))**2.
        end do
    end function
    
    !>  \brief  Rosenbrock's (banana) function
    !!          global minimum: f(1,...,1) = 0
    !!          search domain: [-5,10]
    !!          can be restricted to: [-2.048,2.048]
    function testfun16( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real :: r
        r = 0.
        do i=1,d-1
            r = r+100.*(x(i+1)-x(i)**2.)**2.+(x(i)-1.)**2.
        end do
    end function
    
    ! OTHER
    
    !>  \brief  Pinter's function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-100,100]
    function testfun17( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r, ir, sum1, sum2, sum3, term1, term2, x_i_min_one, x_i_plus_one
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        do i=1,d
            ir = real(i)
            if( i == 1 )then
                x_i_min_one = x(d)
            else
                x_i_min_one = x(i-1)
            endif
            if( i == d )then
                x_i_plus_one = x(1)
            else
                x_i_plus_one = x(i+1)
            endif
            term1 = sin(x_i_min_one*sin(x(i))-x(i)+sin(x_i_plus_one))
            term2 = x_i_min_one**2.-2.*x(i)+3.*x_i_plus_one-cos(x(i))+1.
            sum1 = sum1+ir*x(i)**2.
            sum2 = sum2+20.*ir*term1*term1
            sum3 = sum3+ir*log10(1.+ir*term2*term2)
        end do
        r = sum1+sum2+sum3
    end function
    
    !>  \brief  MultiMod function
    !!          global minimum: f(0,...,0) = 0
    !!          search domain: [-10,10]
    function testfun18( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r, abx, mabx
        r = 0.
        do i=1,d
            abx = abs(x(i))
            if( i == 1 )then
                mabx = abx
            else
                mabx = mabx*abx
            endif
            r = r+abx*mabx
        end do
    end function
    
    !>  \brief  Step function
    !!          global minimum: f(0.5,...,0.5) = 0
    !!          search domain: [-100,100]
    !!          The presence of many flat plateus and steep ridges presents 
    !!          difficulties for algorithms based on gradient information
    function testfun19( x, d ) result( r )
        integer, intent(in) :: d
        real, intent(in)    :: x(d)
        integer :: i
        real    :: r, var
        r = 0.
        do i=1,d
            var = real(floor(x(i)))+0.5
            r = r+var*var
        end do
    end function

    ! TWO-DIMENSIONAL TEST FUNCTIONS
    
    ! MANY LOCAL MINIMA (MULTIMODAL)
    
    !>  \brief  Drop-wave function
    !!          global minimum: f(0,0) = -1
    !!          search domain: [-5.12,5.12]
    function testfun20( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r, frac1, frac2
        frac1 = 1.+ cos(12.*sqrt(x(1)**2.+x(2)**2.))
        frac2 = 0.5*(x(1)**2.+x(2)**2.)+2.
        r = -frac1/frac2
    end function
    
    !>  \brief Eggholder function
    !!         global minimum: f(512,404.2319) = -959.6407 
    !!         search domain: [-512,512]
    function testfun21( x, d )result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r
        r = -(x(2)+47.)*sin(sqrt(abs(x(2)+x(1)/2.+47.)))-x(1)*sin(sqrt(abs(x(1)-(x(2)+47.))))
    end function
    
    !>  \brief  Holder table function
    !!          global minimum: f(8.05502,9.66459)   = -19.2085
    !!          global minimum: f(8.05502,-9.66459)  = -19.2085
    !!          global minimum: f(-8.05502,9.66459)  = -19.2085
    !!          global minimum: f(-8.05502,-9.66459) = -19.2085
    !!          search domain: [-10,10]
    function testfun22( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r, fact1, fact2
        fact1 = sin(x(1))*cos(x(2))
        fact2 = exp(abs(1.-sqrt(x(1)**2.+x(2)**2.)/pi))
        r = -abs(fact1*fact2)
    end function
        
    !>  \brief  Levy function N.13
    !!          global minimum: f(1,1) = 0
    !!          search domain: [-10,10]
    function testfun23( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r, term1, term2, term3
        term1 = (sin(3.*pi*x(1)))**2.
        term2 = (x(1)-1.)**2*(1.+(sin(3.*pi*x(2)))**2.)
        term3 = (x(2)-1.)**2*(1.+(sin(2.*pi*x(2)))**2.)
        r = term1+term2+term3
    end function
    
    !>  \brief Schaffer function N.2
    !!         global minimum: f(0,0) = 0 
    !!         search domain: [-100,100]
    function testfun24( x, d )result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r
        r = 0.5+(sin(x(1)**2.-x(2)**2.)**2.-0.5)/(1.+0.001*(x(1)**2.+x(2)**2.))**2.
    end function
    
    !>  \brief Schaffer function N.4
    !!         global minimum: f(0,1.25313) = 0.292579 
    !!         search domain: [-100,100]
    function testfun25( x, d )result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r
        r = 0.5+(cos(sin(abs(x(1)**2.-x(2)**2.)))-0.5)/(1.+0.001*(x(1)**2.+x(2)**2.))**2.
    end function
    
    !>  \brief Shubert function
    !!         global minimum: f(unknonw) = -186.7309
    !!         search domain: [-5.12,5.12]
    function testfun26( x, d )result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real    :: r, sum1, sum2, ii
        integer :: i
        sum1 = 0.
        sum2 = 0.
        do i=1,5
            ii = real(i)
            sum1 = sum1+ii*cos((ii+1.)*x(1)+ii)
            sum2 = sum2+ii*cos((ii+1.)*x(2)+ii)
        end do
        r = sum1*sum2
    end function
    
    ! PLATE-SHAPED
    
    !>  \brief  Booth's function
    !!          global minimum: f(1,3) = 0
    !!          search domain: [-10,10]
    function testfun27( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r
        r = (x(1)+2.*x(2)-7.)**2.+(2.*x(1)+x(2)-5.)**2.
    end function
    
    !>  \brief  Matya's function
    !!          global minimum: f(0,0) = 0
    !!          search domain: [-10,-10]
    function testfun28( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r
        r = 0.26*(x(1)**2.+x(2)**2)-0.48*x(1)*x(2)
    end function
    
    ! VALLEY-SHAPED
    
    !>  \brief  three-hump camel function
    !!          global minimum: f(0,0) = 0
    !!          search domain: [-5,-5]
    function testfun29( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r
        r = 2.*x(1)**2.-1.05*x(1)**4.+x(1)**6./6.+x(1)*x(2)+x(2)**2.
    end function
    
    !>  \brief  six-hump camel function
    !!          global minimum: f(0.0898,-0.7126) = -1.0316
    !!          global minimum: f(-0.0898,0.7126) = -1.0316
    !!          search domain: [-3,-3]
    function testfun30( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r, term1, term2, term3
        term1 = (4.-2.1*x(1)**2.+(x(1)**4.)/3.)*x(1)**2.
        term2 = x(1)*x(2)
        term3 = (-4.+4.*x(2)**2.)*x(2)**2.
        r = term1+term2+term3
    end function
    
    ! OTHER
    
    !>  \brief  Beale's function
    !!          global minimum: f(3,0.5) = 0
    !!          search domain: [-4.5,4.5]
    !!          characteristics: plateued bowl
    function testfun31( x, d ) result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r
        r = (1.5-x(1)+x(1)*x(2))**2.+(2.25-x(1)+x(1)*x(2)**2.)**2.+(2.625-x(1)+x(1)*x(2)**3.)**2.
    end function
    
    !>  \brief Goldstein-Price function
    !!         global minimum: f(0,-1) = 3
    !!         search domain: [-2,2]
    !!         characteristics: plateued humps
    function testfun32( x, d )result( r )
        integer, intent(in) :: d ! d=2
        real, intent(in)    :: x(d)
        real :: r, term1, term2, term3, term4
        term1 = x(1)+x(2)+1
        term2 = 19.-14.*x(1)+3.*x(1)**2.-14.*x(2)+6.*x(1)*x(2)+3.*x(2)**2.
        term3 = 2.*x(1)-3.*x(2)
        term4 = 18.-32.*x(1)+12.*x(1)**2.+48*x(2)-36.*x(1)*x(2)+27.*x(2)**2.
        r = (1.+term1*term1*term2)*(30.+term3*term3*term4)        
!        r = (1.+(x(1)+x(2)+1)**2.*(19.-14.*x(1)+&
!        3.*x(1)**2.-14.*x(2)+6.*x(1)*x(2)+3.*x(2)**2.))&
!        *(30.+(2.*x(1)-3.*x(2))**2.*(18.-32.*x(1)+12.*x(1)**2.+&
!        48*x(2)-36.*x(1)*x(2)+27.*x(2)**2.))
    end function
    
end module
