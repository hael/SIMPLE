!==Class simple_mgaus
!
! The simple_mgaus is a class for defining an managing multivariate Gaussians. The multivariate Gaussian
! is governed by a D-dimensional vector avgs and a DXD covariance matrix C that must be symmetric (like the unit matrix)
! and positive definite (transpose(z)Cz > 0 for all z not 0) The inverse of the covariance matrix is the precision matrix.
! 
!==Changes are documented below
!
module simple_mgaus
use simple_math
use simple_defs
use simple_jiffys
implicit none
save

! Implement one dimensional case

type mgaus
    private
    real, allocatable :: C(:,:), Cinv(:,:)
    real, allocatable :: avgs(:)
    real :: norm, Cdet
    integer :: D
    logical :: exists=.false.
end type

interface set_mgaus
    module procedure set_mgaus_1
    module procedure set_mgaus_2
end interface

    contains
        
        function new_mgaus( D, avgs, vars ) result( num )
            integer, intent(in)        :: D
            real, intent(in), optional :: avgs(D), vars(D)
            type(mgaus)                :: num
            integer                    :: alloc_stat 
            num%D = D
            allocate( num%C(D,D), num%Cinv(D,D), num%avgs(D), stat=alloc_stat )
            call alloc_err('In: new_mgaus, module: simple_mgaus', alloc_stat)
            if( present(avgs) )then
                if( present(vars) )then
                    call set_mgaus( num, avgs, vars )
                endif
            endif
            num%exists = .true.
        end function new_mgaus
        
        subroutine kill_mgaus( num )
            type(mgaus) :: num
            if( num%exists )then
                deallocate( num%C, num%Cinv, num%avgs )
                num%exists = .false.
            endif
        end subroutine kill_mgaus
        
        function get_mgaus_cov( num ) result( cov )
            type(mgaus) :: num
            real :: cov(num%D,num%D)
            cov = num%C
        end function get_mgaus_cov
        
        function get_mgaus_avgs( num ) result( avgs )
            type(mgaus) :: num
            real :: avgs(num%D)
            avgs = num%avgs
        end function get_mgaus_avgs

        subroutine set_mgaus_1( num, avgs, vars )
            type(mgaus)      :: num
            real, intent(in) :: avgs(num%D), vars(num%D)
            integer          :: errflg, i
            ! make covariance matrix
            num%C = 0.
            do i=1,num%D ; num%C(i,i) = vars(i) ; end do
            ! calculate its determinant
            num%Cdet = det(num%C, num%D)
            if( num%Cdet > 0. )then
                ! alles ok!
            else
                write(*,*) 'ERROR, covariance matrix must be positive definite!'
                write(*,*) 'In: set_mgaus_1, module: simple_mgaus'
                stop
            endif
            ! calculate normalization constant
            num%norm = 1./((TWOPI**(real(num%D)/2.))*sqrt(num%Cdet))
            ! invert the covariance matrix
            call matinv(num%C, num%Cinv, num%D, errflg)
            if( errflg == -1 )then
                write(*,*) 'ERROR, covariance matrix not invertible!'
                write(*,*) 'In: set_mgaus_1, module: simple_mgaus'
                stop
            endif
            ! set avgs
            num%avgs = avgs
        end subroutine set_mgaus_1
        
        subroutine set_mgaus_2( num, avgs, cov, eps )
            type(mgaus)      :: num
            real, intent(in) :: avgs(num%D), cov(num%D,num%D), eps
            integer          :: errflg, i
            ! set covariance matrix
            num%C = cov(num%D,num%D)
            ! calculate its determinant
            num%Cdet = det(num%C, num%D)
            if( num%Cdet > 0. )then
                ! alles ok!
            else
                ! make covariance matrix
                num%C = 0.
                do i=1,num%D ; num%C(i,i) = eps ; end do
            endif
            ! calculate normalization constant
            num%norm = 1./((TWOPI**(real(num%D)/2.))*sqrt(num%Cdet))
            ! invert the covariance matrix
            call matinv(num%C, num%Cinv, num%D, errflg)
            if( errflg == -1 )then
                write(*,*) 'ERROR, covariance matrix not invertible!'
                write(*,*) 'In: set_mgaus_2, module: simple_mgaus'
                stop
            endif
            ! set avgs
            num%avgs = avgs
        end subroutine set_mgaus_2
        
        function sample_mgaus( num, x ) result( y )
            type(mgaus)      :: num
            real, intent(in) :: x(num%D)
            real             :: y, diffvec(num%D,1), tmp(1,num%D)
            diffvec(:,1) = x-num%avgs
            tmp = matmul(transpose(diffvec),num%Cinv)
            y = num%norm*exp(-0.5*dot_product(tmp(1,:),diffvec(:,1)))
        end function sample_mgaus

end module simple_mgaus