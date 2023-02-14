module simple_linalg
use simple_defs
use simple_is_check_assert
implicit none

! linear algebra

interface eigsrt
    module procedure eigsrt_sp, eigsrt_dp
end interface eigsrt

interface jacobi
    module procedure jacobi_sp, jacobi_dp
end interface jacobi

interface matinv
    module procedure matinv_sp, matinv_dp
end interface matinv

interface norm_2
    module procedure norm_2_sp, norm_2_dp
end interface norm_2

interface outerprod
    module procedure outerprod_r, outerprod_d
end interface outerprod

interface svbksb
    module procedure svbksb_sp, svbksb_dp
end interface svbksb

interface svdcmp
    module procedure svdcmp_sp, svdcmp_dp
end interface svdcmp

interface svdfit
    module procedure svdfit_sp, svdfit_dp
end interface svdfit

interface svd_multifit
    module procedure svd_multifit_sp, svd_multifit_dp
end interface svd_multifit

interface vabs
    module procedure vabs_sp, vabs_dp
end interface vabs

! trigonometry

interface myacos
    module procedure myacos_sp, myacos_dp
end interface myacos

interface deg2rad
    module procedure deg2rad_sp, deg2rad_dp
end interface deg2rad

interface rad2deg
    module procedure rad2deg_1, rad2deg_2
end interface

interface pythag
    module procedure pythag_sp, pythag_dp
end interface pythag

contains

    !>   calculates the argument of a vector
    pure function arg( vec ) result( length )
        real, intent(in) :: vec(:)
        real :: length
        length = sqrt(sum(vec*vec))
    end function arg

    !>  \brief sorts eigenvalues and eigenvectors from jacobi routine in descending order
    ! Given the eigenvalues d and eigenvectors v as output from jacobi (§11.1) or tqli (§11.3),
    ! this routine sorts the eigenvalues into descending order,
    ! and rearranges the columns of v correspondingly. The method is straight insertion.
    subroutine eigsrt_sp(d,v,n,np)
        integer, intent(in)    :: n,np
        real,    intent(inout) :: d(np),v(np,np)
        integer             :: i,j,k
        real                :: p
        do i=1,n-1
            k = i
            p = d(i)
            do j=i+1,n
                if(d(j)>p)then
                    k = j
                    p = d(j)
                endif
            enddo
            if(k.ne.i)then
                d(k) = d(i)
                d(i) = p
                do j=1,n
                    p      = v(j,i)
                    v(j,i) = v(j,k)
                    v(j,k) = p
                enddo
            endif
        enddo
    end subroutine eigsrt_sp

    subroutine eigsrt_dp(d,v,n,np)
        integer,  intent(in)    :: n,np
        real(dp), intent(inout) :: d(np),v(np,np)
        integer  :: i,j,k
        real(dp) :: p
        do i=1,n-1
            k = i
            p = d(i)
            do j=i+1,n
                if(d(j)>p)then
                    k = j
                    p = d(j)
                endif
            enddo
            if(k.ne.i)then
                d(k) = d(i)
                d(i) = p
                do j=1,n
                    p      = v(j,i)
                    v(j,i) = v(j,k)
                    v(j,k) = p
                enddo
            endif
        enddo
    end subroutine eigsrt_dp

    subroutine fit_lsq_plane(n, xyz, A,B,C, err)
        integer, intent(in)  :: n
        real,    intent(in)  :: xyz(n,3)
        real,    intent(out) :: A,B,C ! plane equation: z = Ax + By + C
        logical, intent(out) :: err
        real(dp) :: sx, sy, sz, sxx, syy, sxy, sxz, syz, denom,rn
        err = .false.
        A   = 0.
        B   = 0.
        C   = 0.
        rn  = real(n,dp)
        sx  = sum(real(xyz(:,1),dp))
        sy  = sum(real(xyz(:,2),dp))
        sz  = sum(real(xyz(:,3),dp))
        sxx = sum(real(xyz(:,1),dp)**2)
        syy = sum(real(xyz(:,2),dp)**2)
        sxy = sum(real(xyz(:,1),dp)*real(xyz(:,2),dp))
        sxz = sum(real(xyz(:,1),dp)*real(xyz(:,3),dp))
        syz = sum(real(xyz(:,2),dp)*real(xyz(:,3),dp))
        denom = sx*sx*syy - 2.d0*sxy*sx*sy + sxx*sy*sy + rn*(sxy*sxy - sxx*syy)
        if( abs(denom) < 1.d-10 )then
            err = .true.
            return
        endif
        A = (sy*sy*sxz  - syy*rn*sxz + sxy*rn*syz + sx*syy*sz  - sy*(sx*syz  + sxy*sz) ) / denom
        B = (sxy*rn*sxz + sx*sx*syz  - sxx*rn*syz + sxx*sy*sz  - sx*(sy*sxz  + sxy*sz) ) / denom
        C = (sx*syy*sxz - sxy*sy*sxz - sxy*sx*syz + sxx*sy*syz + sz*(sxy*sxy - sxx*syy)) / denom
    end subroutine fit_lsq_plane

    !>   least squares straight line fit, from Sjors, he took it from
    !!         http://mathworld.wolfram.com/LeastSquaresFitting.html
    !!         ss_xx = sum_i x_i^2 - n * ave_x^2
    !!         ss_yy = sum_i y_i^2 - n * ave_y^2
    !!         ss_xy = sum_i x_i * y_i - n * ave_x * n ave_y
    !! \param  slope  Linear gradient: slope = xx_xy / ss_xx
    !!  \param intercept y-intercept:  intercept = ave_y - slope * ave_x
    !!  \param corr Correlation coefficient: corr_coeff = ss_xy^2 / (ss_xx * ss_yy)
    subroutine fit_straight_line( n, datavec, slope, intercept, corr )
        integer, intent(in) :: n                                            !< size of vector
        real, intent(in)    :: datavec(n,2)                                 !< input vector
        real, intent(out)   :: slope, intercept, corr
        double precision    :: ave_x, ave_y, ss_xx, ss_yy, ss_xy, x, y, dn
        integer             :: i
        ave_x = 0.d0
        ave_y = 0.d0
        ss_xx = 0.d0
        ss_yy = 0.d0
        ss_xy = 0.d0
        do i=1,n
            x = datavec(i,1)
            y = datavec(i,2)
            ave_x = ave_x+x
            ave_y = ave_y+y
            ss_xx = ss_xx+x*x
            ss_yy = ss_yy+y*y
            ss_xy = ss_xy+x*y
        end do
        dn = dble(n)
        ave_x = ave_x/dn
        ave_y = ave_y/dn
        ss_xx = ss_xx-dn*ave_x*ave_x
        ss_yy = ss_yy-dn*ave_y*ave_y
        ss_xy = ss_xy-dn*ave_x*ave_y
        slope     = real(ss_xy/ss_xx)
        intercept = real(ave_y-slope*ave_x)
        corr      = real((ss_xy*ss_xy)/(ss_xx*ss_yy))
    end subroutine fit_straight_line

    !> \brief jacobi SVD, NR
    subroutine jacobi_sp( a, n, np, d, v, nrot)
        integer, intent(in)    :: n,np
        real,    intent(inout) :: a(np,np), v(np,np), d(np)
        integer, intent(inout) :: nrot
        real                   :: c,g,h,s,sm,t,tau,theta,tresh,b(n), z(n)
        integer                :: i,j,ip,iq
        v = 0.
        do i=1,n
            v(i,i) = 1.
            b(i)   = a(i,i)
        enddo
        d    = b
        z    = 0.
        nrot = 0
        do i=1,50
            sm=0.
            do ip=1,n-1
                do iq=ip+1,n
                    sm=sm+abs(a(ip,iq))
                enddo
            enddo
            if(sm.eq.0.)return
            if(i.lt.4)then
                tresh=0.2*sm/n**2
            else
                tresh=0.
            endif
            do ip=1,n-1
                do iq=ip+1,n
                    g=100.*abs(a(ip,iq))
                    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) &
                        .and.(abs(d(iq))+g.eq.abs(d(iq))))then
                        a(ip,iq)=0.
                    else if(abs(a(ip,iq)).gt.tresh)then
                        h=d(iq)-d(ip)
                        if(abs(h)+g.eq.abs(h))then
                            t=a(ip,iq)/h
                        else
                            theta=0.5*h/a(ip,iq)
                            t=1./(abs(theta)+sqrt(1.+theta**2))
                            if(theta.lt.0.)t=-t
                        endif
                        c=1./sqrt(1+t**2)
                        s=t*c
                        tau=s/(1.+c)
                        h=t*a(ip,iq)
                        z(ip)=z(ip)-h
                        z(iq)=z(iq)+h
                        d(ip)=d(ip)-h
                        d(iq)=d(iq)+h
                        a(ip,iq)=0.
                        do j=1,ip-1
                            g=a(j,ip)
                            h=a(j,iq)
                            a(j,ip)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        enddo
                        do j=ip+1,iq-1
                            g=a(ip,j)
                            h=a(j,iq)
                            a(ip,j)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        enddo
                        do j=iq+1,n
                            g=a(ip,j)
                            h=a(iq,j)
                            a(ip,j)=g-s*(h+g*tau)
                            a(iq,j)=h+s*(g-h*tau)
                        enddo
                        do j=1,n
                            g=v(j,ip)
                            h=v(j,iq)
                            v(j,ip)=g-s*(h+g*tau)
                            v(j,iq)=h+s*(g-h*tau)
                        enddo
                        nrot=nrot+1
                    endif
                enddo
            enddo
            do ip=1,n
                b(ip)=b(ip)+z(ip)
                d(ip)=b(ip)
                z(ip)=0.
            enddo
        enddo
        write(*,*)' Too many iterations in Jacobi'
    end subroutine jacobi_sp

    !> \brief jacobi SVD, NR
    subroutine jacobi_dp( a, n, np, d, v, nrot )
        integer,  intent(in)    :: n,np
        real(dp), intent(inout) :: a(np,np), v(np,np), d(np)
        integer,  intent(inout) :: nrot
        real(dp) :: c,g,h,s,sm,t,tau,theta,tresh,b(n), z(n)
        integer  :: i,j,ip,iq
        v = 0.
        do i=1,n
            v(i,i) = 1.
            b(i)   = a(i,i)
        enddo
        d    = b
        z    = 0.
        nrot = 0
        do i=1,50
            sm=0.
            do ip=1,n-1
                do iq=ip+1,n
                    sm=sm+abs(a(ip,iq))
                enddo
            enddo
            if(sm.eq.0.)return
            if(i.lt.4)then
                tresh=0.2*sm/n**2
            else
                tresh=0.
            endif
            do ip=1,n-1
                do iq=ip+1,n
                    g=100.*abs(a(ip,iq))
                    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) &
                        .and.(abs(d(iq))+g.eq.abs(d(iq))))then
                        a(ip,iq)=0.
                    else if(abs(a(ip,iq)).gt.tresh)then
                        h=d(iq)-d(ip)
                        if(abs(h)+g.eq.abs(h))then
                            t=a(ip,iq)/h
                        else
                            theta=0.5*h/a(ip,iq)
                            t=1./(abs(theta)+sqrt(1.+theta**2))
                            if(theta.lt.0.)t=-t
                        endif
                        c=1./sqrt(1+t**2)
                        s=t*c
                        tau=s/(1.+c)
                        h=t*a(ip,iq)
                        z(ip)=z(ip)-h
                        z(iq)=z(iq)+h
                        d(ip)=d(ip)-h
                        d(iq)=d(iq)+h
                        a(ip,iq)=0.
                        do j=1,ip-1
                            g=a(j,ip)
                            h=a(j,iq)
                            a(j,ip)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        enddo
                        do j=ip+1,iq-1
                            g=a(ip,j)
                            h=a(j,iq)
                            a(ip,j)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        enddo
                        do j=iq+1,n
                            g=a(ip,j)
                            h=a(iq,j)
                            a(ip,j)=g-s*(h+g*tau)
                            a(iq,j)=h+s*(g-h*tau)
                        enddo
                        do j=1,n
                            g=v(j,ip)
                            h=v(j,iq)
                            v(j,ip)=g-s*(h+g*tau)
                            v(j,iq)=h+s*(g-h*tau)
                        enddo
                        nrot=nrot+1
                    endif
                enddo
            enddo
            do ip=1,n
                b(ip)=b(ip)+z(ip)
                d(ip)=b(ip)
                z(ip)=0.
            enddo
        enddo
        write(*,*)' Too many iterations in Jacobi'
    end subroutine jacobi_dp

    subroutine LUDCMP(A,N,NP,INDX,D,CODE)
        integer, intent(in)    :: N,NP
        real,    intent(inout) :: A(NP,NP)
        integer, intent(inout) :: INDX(N), D, CODE
        integer, PARAMETER :: NMAX=100
        real,    PARAMETER :: TINY=1E-12
        real    :: AMAX,DUM, SUM, VV(NMAX)
        INTEGER :: I, IMAX, J, K
        D=1; CODE=0
        DO I=1,N
        AMAX=0.
        DO J=1,N
            IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
        END DO ! j loop
        IF(AMAX.LT.TINY) THEN
            CODE = 1
            RETURN
        END IF
        VV(I) = 1. / AMAX
        END DO ! i loop
        DO J=1,N
            DO I=1,J-1
                SUM = A(I,J)
                DO K=1,I-1
                    SUM = SUM - A(I,K)*A(K,J)
                END DO ! k loop
                A(I,J) = SUM
            END DO ! i loop
            AMAX = 0.
            DO I=J,N
                SUM = A(I,J)
                DO K=1,J-1
                    SUM = SUM - A(I,K)*A(K,J)
                END DO ! k loop
                A(I,J) = SUM
                DUM = VV(I)*ABS(SUM)
                IF(DUM.GE.AMAX) THEN
                    IMAX = I
                    AMAX = DUM
                END IF
            END DO ! i loop
            IF(J.NE.IMAX) THEN
                DO K=1,N
                    DUM = A(IMAX,K)
                    A(IMAX,K) = A(J,K)
                    A(J,K) = DUM
                END DO ! k loop
                D = -D
                VV(IMAX) = VV(J)
            END IF
            INDX(J) = IMAX
            IF(ABS(A(J,J)) < TINY) A(J,J) = TINY
            IF(J.NE.N) THEN
                DUM = 1. / A(J,J)
                DO I=J+1,N
                A(I,J) = A(I,J)*DUM
                END DO ! i loop
            END IF
        END DO ! j loop
    end subroutine LUDCMP

    subroutine LUBKSB(A,N,NP,INDX,B)
        ! From Numerical Recipes
        integer, intent(in) :: N,NP
        real,    intent(inout) :: A(NP,NP),B(N)
        integer, intent(inout) :: INDX(N)
        real    :: SUM
        integer :: I, II, J, LL
        II = 0
        DO I=1,N
            LL = INDX(I)
            SUM = B(LL)
            B(LL) = B(I)
            IF(II.NE.0) THEN
                DO J=II,I-1
                    SUM = SUM - A(I,J)*B(J)
                END DO ! j loop
                ELSE IF(SUM.NE.0.) THEN
                    II = I
                END IF
            B(I) = SUM
        END DO ! i loop
        DO I=N,1,-1
            SUM = B(I)
            IF(I < N) THEN
                DO J=I+1,N
                   SUM = SUM - A(I,J)*B(J)
                END DO ! j loop
            END IF
            B(I) = SUM / A(I,I)
        END DO ! i loop
    end subroutine LUBKSB

    !>   subroutine to find the inverse of a square matrix
    !!         author : louisda16th a.k.a ashwith j. rego
    !!         reference : algorithm explained at:
    !!         http://math.uww.edu/~mcfarlat/inverse.htm
    !!         http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
    subroutine matinv_sp(matrix, inverse, n, errflg)
        integer, intent(in)  :: n       !< size of square matrix
        integer, intent(out) :: errflg  !< return error status. -1 for error, 0 for normal
        real, intent(in), dimension(n,n)  :: matrix  !< input matrix
        real, intent(out), dimension(n,n) :: inverse !< inverted matrix
        logical :: flag = .true.
        integer :: i, j, k
        real :: m
        real, dimension(n,2*n) :: augmatrix !< augmented matrix
        ! augment input matrix with an identity matrix
        do i=1,n
            do j=1,2*n
                if(j <= n )then
                    augmatrix(i,j) = matrix(i,j)
                else if((i+n) == j)then
                    augmatrix(i,j) = 1.
                else
                    augmatrix(i,j) = 0.
                endif
            end do
        end do
        ! reduce augmented matrix to upper traingular form
        do k=1,n-1
              if( is_zero(augmatrix(k,k)) )then
                 flag = .false.
                do i=k+1,n
                    if( abs(augmatrix(i,k)) > 0. )then
                        do j=1,2*n
                            augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                        end do
                        flag = .true.
                        exit
                    endif
                end do
                if(flag .eqv. .false.)then
                    inverse = 0.
                    errflg = -1
                    return
                endif
            endif
            do j=k+1, n
                m = augmatrix(j,k)/augmatrix(k,k)
                do i=k,2*n
                    augmatrix(j,i) = augmatrix(j,i)-m*augmatrix(k,i)
                end do
            end do
        end do
        ! test for invertibility
        do i=1,n
            if( is_zero(augmatrix(i,i)) )then
                inverse = 0
                errflg = -1
                return
            endif
        end do
        ! make diagonal elements as 1
        do i=1,n
            m = augmatrix(i,i)
            do j=i,2*n
                augmatrix(i,j) = augmatrix(i,j)/m
            end do
        end do
        ! reduced right side half of augmented matrix to identity matrix
        do k=n-1,1,-1
            do i=1,k
                m = augmatrix(i,k+1)
                do j = k,2*n
                    augmatrix(i,j) = augmatrix(i,j)-augmatrix(k+1,j)*m
                end do
            end do
        end do
        ! store answer
        do i=1,n
            do j=1,n
                inverse(i,j) = augmatrix(i,j+n)
            end do
        end do
        errflg = 0
    end subroutine matinv_sp

    subroutine matinv_dp(matrix, inverse, n, errflg)
        integer, intent(in)  :: n       !< size of square matrix
        integer, intent(out) :: errflg  !< return error status. -1 for error, 0 for normal
        real(dp), intent(in), dimension(n,n)  :: matrix  !< input matrix
        real(dp), intent(out), dimension(n,n) :: inverse !< inverted matrix
        logical :: flag = .true.
        integer :: i, j, k
        real(dp) :: m
        real(dp), dimension(n,2*n) :: augmatrix !< augmented matrix
        ! augment input matrix with an identity matrix
        do i=1,n
            do j=1,2*n
                if(j <= n )then
                    augmatrix(i,j) = matrix(i,j)
                else if((i+n) == j)then
                    augmatrix(i,j) = 1.d0
                else
                    augmatrix(i,j) = 0.d0
                endif
            end do
        end do
        ! reduce augmented matrix to upper traingular form
        do k=1,n-1
            if( is_zero(augmatrix(k,k)) )then
                 flag = .false.
                do i=k+1,n
                    if( abs(augmatrix(i,k)) > 0.d0 )then
                        do j=1,2*n
                            augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                        end do
                        flag = .true.
                        exit
                    endif
                end do
                if(flag .eqv. .false.)then
                    inverse = 0.d0
                    errflg = -1
                    return
                endif
            endif
            do j=k+1, n
                m = augmatrix(j,k)/augmatrix(k,k)
                do i=k,2*n
                    augmatrix(j,i) = augmatrix(j,i)-m*augmatrix(k,i)
                end do
            end do
        end do
        ! test for invertibility
        do i=1,n
            if( is_zero(augmatrix(i,i)) )then
                inverse = 0.d0
                errflg = -1
                return
            endif
        end do
        ! make diagonal elements as 1
        do i=1,n
            m = augmatrix(i,i)
            do j=i,2*n
                augmatrix(i,j) = augmatrix(i,j)/m
            end do
        end do
        ! reduced right side half of augmented matrix to identity matrix
        do k=n-1,1,-1
            do i=1,k
                m = augmatrix(i,k+1)
                do j = k,2*n
                    augmatrix(i,j) = augmatrix(i,j)-augmatrix(k+1,j)*m
                end do
            end do
        end do
        ! store answer
        do i=1,n
            do j=1,n
                inverse(i,j) = augmatrix(i,j+n)
            end do
        end do
        errflg = 0
    end subroutine matinv_dp

    pure function norm_2_sp(v) result(r)
        real, intent(in) :: v(:)
        real             :: r
        r = sqrt(dot_product(v,v))
    end function norm_2_sp

    pure function norm_2_dp(v) result(r)
        real(dp), intent(in) :: v(:)
        real(dp)             :: r
        r = sqrt(dot_product(v,v))
    end function norm_2_dp

    ! imported from numerical recipes
    function outerprod_r(a,b)
        implicit none
        real(sp), dimension(:), intent(in) :: a,b
        real(sp), dimension(size(a),size(b)) :: outerprod_r
        outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
            spread(b,dim=1,ncopies=size(a))
    end function outerprod_r

    ! imported from numerical recipes
    function outerprod_d(a,b)
        implicit none
        real(dp), dimension(:), intent(in) :: a,b
        real(dp), dimension(size(a),size(b)) :: outerprod_d
        outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
            spread(b,dim=1,ncopies=size(a))
    end function outerprod_d

    ! Find the plane that minimises the distance between
    ! a given set of points.
    ! It consists in a solution of a overdetermined system with
    ! the left pseudo inverse.
    ! SOURCE :
    ! https://stackoverflow.com/questions/1400213/3d-least-squares-plane
    ! The output plane will have cartesian equation
    ! vec(1)x + vec(2)y - z = -vec(3).
    ! FORMULA
    ! sol = inv(transpose(M)*M)*transpose(M)*b
    function plane_from_points(points) result(sol)
        real, intent(inout) :: points(:,:) !input
        real    :: sol(3)  !vec(1)x + vec(2)y - z = -vec(3).
        real    :: M(size(points, dim = 2),3), b(size(points, dim = 2)), invM(3,size(points, dim = 2))
        real    :: prod(3,3), prod_inv(3,3), prod1(3,size(points, dim = 2))
        integer :: errflg ! if manages to find inverse matrix
        integer :: p
        integer :: N ! number of points
        if(size(points, dim=1) /=3) then
            write(logfhandle,*) 'Need to input points in 3D!; plane_from_points'
            return
        endif
        if(size(points, dim=2) < 3) then
            write(logfhandle,*) 'Not enough input points for fitting!; plane_from_points'
            return
        endif
        N = size(points, dim=2)
        do p = 1, N
            M(p,1) =  points(1,p)
            M(p,2) =  points(2,p)
            M(p,3) =  1.
            b(p)   =  points(3,p)
        enddo
        prod  = matmul(transpose(M),M)
        call matinv(prod,prod_inv,3,errflg)
        if( errflg /= 0 ) then
            write(logfhandle,*) 'Couldn t find inverse matrix! ;plane_from_points'
            stop
        endif
        prod1 = matmul(prod_inv,transpose(M))
        sol   = matmul(prod1,b)
    end function plane_from_points

    !>   projects a 3d vector in the _z_-direction
    subroutine projz( vec3, vec2 )
        real, intent(in)  :: vec3(3)
        real, intent(out) :: vec2(2)
        vec2(1) = vec3(1)
        vec2(2) = vec3(2)
    end subroutine projz

    ! imported from numerical recipes
    subroutine svbksb_sp(u,w,v,b,x)
        implicit none
        real(sp), dimension(:,:), intent(in)  :: u,v
        real(sp), dimension(:),   intent(in)  :: w,b
        real(sp), dimension(:),   intent(out) :: x
        ! Solves A X=B for a vector X, where A is specified by the arrays u,v,w as returned
        ! by svdcmp. Here u is MxN, v is NxN, and w is of length N. b is the M-dimensional
        ! input right-hand side. x is the N-dimensional output solution vector. No input quantities
        ! are destroyed, so the routine may be called sequentially with different b's.
        integer                      :: mdum,ndum
        real(sp), dimension(size(x)) :: tmp
        mdum=assert_eq(size(u,1),size(b),'svbksb_sp: mdum')
        ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),'svbksb_sp: ndum')
        where (w /= 0.0)
            tmp=matmul(b,u)/w     ! Calculate diag(1/w_j)U^T B,
        elsewhere
            tmp=0.0               ! but replace 1/w_j by zero if w_j=0.
        end where
        x=matmul(v,tmp)           ! Matrix multiply by V to get answer.
    end subroutine svbksb_sp

    ! imported from numerical recipes
    subroutine svbksb_dp(u,w,v,b,x)
        implicit none
        real(dp), dimension(:,:), intent(in ) :: u,v
        real(dp), dimension(:),   intent(in)  :: w,b
        real(dp), dimension(:),   intent(out) :: x
        integer                      :: mdum,ndum
        real(dp), dimension(size(x)) :: tmp
        mdum=assert_eq(size(u,1),size(b),'svbksb_dp: mdum')
        ndum=assert_eq((/size(u,2),size(v,1),size(v,2),size(w),size(x)/),'svbksb_dp: ndum')
        where (w /= 0.0)
            tmp=matmul(b,u)/w
        elsewhere
            tmp=0.0
        end where
        x=matmul(v,tmp)
    end subroutine svbksb_dp

    ! imported from numerical recipes
    ! Given an MxN matrix a, this routine computes its singular value decomposition, A=
    ! U W V^T. The matrix U replaces a on output. The diagonal matrix of singular values
    ! W is output as the N-dimensional vector w. The NxN matrix V (not the transpose V^T)
    ! is output as v.
    subroutine svdcmp_sp(a,w,v)
        implicit none
        real(sp), dimension(:,:), intent(inout) :: a
        real(sp), dimension(:),   intent(out)   :: w
        real(sp), dimension(:,:), intent(out)   :: v
        integer                        :: i,its,j,k,l,m,n,nm
        real(sp)                       :: anorm,c,f,g,h,s,scale,x,y,z
        real(sp), dimension(size(a,1)) :: tempm
        real(sp), dimension(size(a,2)) :: rv1,tempn
        m=size(a,1)
        n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_sp')
        g=0.0
        scale=0.0
        do i=1,n              ! Householder reduction to bidiagonal form.
            l=i+1
            rv1(i)=scale*g
            g=0.0
            scale=0.0
            if (i <= m) then
                scale=sum(abs(a(i:m,i)))
                if (scale /= 0.0) then
                    a(i:m,i)=a(i:m,i)/scale
                    s=dot_product(a(i:m,i),a(i:m,i))
                    f=a(i,i)
                    g=-sign(sqrt(s),f)
                    h=f*g-s
                    a(i,i)=f-g
                    tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                    a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                    a(i:m,i)=scale*a(i:m,i)
                end if
            end if
            w(i)=scale*g
            g=0.0
            scale=0.0
            if ((i <= m) .and. (i /= n)) then
                scale=sum(abs(a(i,l:n)))
                if (scale /= 0.0) then
                    a(i,l:n)=a(i,l:n)/scale
                    s=dot_product(a(i,l:n),a(i,l:n))
                    f=a(i,l)
                    g=-sign(sqrt(s),f)
                    h=f*g-s
                    a(i,l)=f-g
                    rv1(l:n)=a(i,l:n)/h
                    tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                    a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                    a(i,l:n)=scale*a(i,l:n)
                end if
            end if
        end do
        anorm=maxval(abs(w)+abs(rv1))
        do i=n,1,-1            ! Accumulation of right-hand transformations.
            if (i < n) then
                if (g /= 0.0) then
                    v(l:n,i)=(a(i,l:n)/a(i,l))/g     ! Double division to avoid possible underflow
                    tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                    v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
                end if
                v(i,l:n)=0.0
                v(l:n,i)=0.0
            end if
            v(i,i)=1.0
            g=rv1(i)
            l=i
        end do
        do i=min(m,n),1,-1        ! Accumulation of left-hand transformations.
            l=i+1
            g=w(i)
            a(i,l:n)=0.0
            if (g /= 0.0) then
                g=1.0_sp/g
                tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=a(i:m,i)*g
            else
                a(i:m,i)=0.0
            end if
            a(i,i)=a(i,i)+1.0_sp
        end do
        do k=n,1,-1                 ! Diagonalization of the bidiagonal form: Loop over
            do its=1,30             ! singular values, and over allowed iterations.
                do l=k,1,-1         ! Test for splitting.
                    nm=l-1
                    if ((abs(rv1(l))+anorm) == anorm) exit
                    ! Note that rv1(1) is always zero, so can never fall through bottom of loop.
                    if ((abs(w(nm))+anorm) == anorm) then
                        c=0.0       ! Cancellation of rv1(l), if l>1.
                        s=1.0
                        do i=l,k
                            f=s*rv1(i)
                            rv1(i)=c*rv1(i)
                            if ((abs(f)+anorm) == anorm) exit
                            g=w(i)
                            h=pythag(f,g)
                            w(i)=h
                            h=1.0_sp/h
                            c= (g*h)
                            s=-(f*h)
                            tempm(1:m)=a(1:m,nm)
                            a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                            a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                        end do
                        exit
                    end if
                end do
                z=w(k)
                if (l == k) then        ! Convergence.
                    if (z < 0.0) then   ! Singular value is made nonnegative.
                        w(k)=-z
                        v(1:n,k)=-v(1:n,k)
                    end if
                    exit
                end if
                if (its == 30) stop 'svdcmp_sp: no convergence in svdcmp'
                x=w(l)                   ! Shift from bottom 2-by-2 minor.
                nm=k-1
                y=w(nm)
                g=rv1(nm)
                h=rv1(k)
                f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_sp*h*y)
                g=pythag(f,1.0_sp)
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
                c=1.0                    ! Next QR transformation:
                s=1.0
                do j=l,nm
                    i=j+1
                    g=rv1(i)
                    y=w(i)
                    h=s*g
                    g=c*g
                    z=pythag(f,h)
                    rv1(j)=z
                    c=f/z
                    s=h/z
                    f= (x*c)+(g*s)
                    g=-(x*s)+(g*c)
                    h=y*s
                    y=y*c
                    tempn(1:n)=v(1:n,j)
                    v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                    v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                    z=pythag(f,h)
                    w(j)=z            ! Rotation can be arbitrary if z=0.
                    if (z /= 0.0) then
                        z=1.0_sp/z
                        c=f*z
                        s=h*z
                    end if
                    f= (c*g)+(s*y)
                    x=-(s*g)+(c*y)
                    tempm(1:m)=a(1:m,j)
                    a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                    a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                end do
                rv1(l)=0.0
                rv1(k)=f
                w(k)=x
            end do
        end do
    end subroutine svdcmp_sp

    subroutine svdcmp_dp(a,w,v)
        implicit none
        real(dp), dimension(:,:), intent(inout) :: a
        real(dp), dimension(:),   intent(out)   :: w
        real(dp), dimension(:,:), intent(out)   :: v
        integer                        :: i,its,j,k,l,m,n,nm
        real(dp)                       :: anorm,c,f,g,h,s,scale,x,y,z
        real(dp), dimension(size(a,1)) :: tempm
        real(dp), dimension(size(a,2)) :: rv1,tempn
        m=size(a,1)
        n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_dp')
        g=0.0
        scale=0.0
        do i=1,n
            l=i+1
            rv1(i)=scale*g
            g=0.0
            scale=0.0
            if (i <= m) then
                scale=sum(abs(a(i:m,i)))
                if (scale /= 0.0) then
                    a(i:m,i)=a(i:m,i)/scale
                    s=dot_product(a(i:m,i),a(i:m,i))
                    f=a(i,i)
                    g=-sign(sqrt(s),f)
                    h=f*g-s
                    a(i,i)=f-g
                    tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
                    a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                    a(i:m,i)=scale*a(i:m,i)
                end if
            end if
            w(i)=scale*g
            g=0.0
            scale=0.0
            if ((i <= m) .and. (i /= n)) then
                scale=sum(abs(a(i,l:n)))
                if (scale /= 0.0) then
                    a(i,l:n)=a(i,l:n)/scale
                    s=dot_product(a(i,l:n),a(i,l:n))
                    f=a(i,l)
                    g=-sign(sqrt(s),f)
                    h=f*g-s
                    a(i,l)=f-g
                    rv1(l:n)=a(i,l:n)/h
                    tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
                    a(l:m,l:n)=a(l:m,l:n)+outerprod(tempm(l:m),rv1(l:n))
                    a(i,l:n)=scale*a(i,l:n)
                end if
            end if
        end do
        anorm=maxval(abs(w)+abs(rv1))
        do i=n,1,-1
            if (i < n) then
                if (g /= 0.0) then
                    v(l:n,i)=(a(i,l:n)/a(i,l))/g
                    tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
                    v(l:n,l:n)=v(l:n,l:n)+outerprod(v(l:n,i),tempn(l:n))
                end if
                v(i,l:n)=0.0
                v(l:n,i)=0.0
            end if
            v(i,i)=1.0
            g=rv1(i)
            l=i
        end do
        do i=min(m,n),1,-1
            l=i+1
            g=w(i)
            a(i,l:n)=0.0
            if (g /= 0.0) then
                g=1.0_dp/g
                tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
                a(i:m,l:n)=a(i:m,l:n)+outerprod(a(i:m,i),tempn(l:n))
                a(i:m,i)=a(i:m,i)*g
            else
                a(i:m,i)=0.0
            end if
            a(i,i)=a(i,i)+1.0_dp
        end do
        do k=n,1,-1
            do its=1,30
                do l=k,1,-1
                    nm=l-1
                    if ((abs(rv1(l))+anorm) == anorm) exit
                    if ((abs(w(nm))+anorm) == anorm) then
                        c=0.0
                        s=1.0
                        do i=l,k
                            f=s*rv1(i)
                            rv1(i)=c*rv1(i)
                            if ((abs(f)+anorm) == anorm) exit
                            g=w(i)
                            h=pythag(f,g)
                            w(i)=h
                            h=1.0_dp/h
                            c= (g*h)
                            s=-(f*h)
                            tempm(1:m)=a(1:m,nm)
                            a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                            a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                        end do
                        exit
                    end if
                end do
                z=w(k)
                if (l == k) then
                    if (z < 0.0) then
                        w(k)=-z
                        v(1:n,k)=-v(1:n,k)
                    end if
                    exit
                end if
                if (its == 30) stop 'svdcmp_dp: no convergence in svdcmp'
                x=w(l)
                nm=k-1
                y=w(nm)
                g=rv1(nm)
                h=rv1(k)
                f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0_dp*h*y)
                g=pythag(f,1.0_dp)
                f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
                c=1.0
                s=1.0
                do j=l,nm
                    i=j+1
                    g=rv1(i)
                    y=w(i)
                    h=s*g
                    g=c*g
                    z=pythag(f,h)
                    rv1(j)=z
                    c=f/z
                    s=h/z
                    f= (x*c)+(g*s)
                    g=-(x*s)+(g*c)
                    h=y*s
                    y=y*c
                    tempn(1:n)=v(1:n,j)
                    v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
                    v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
                    z=pythag(f,h)
                    w(j)=z
                    if (z /= 0.0) then
                        z=1.0_dp/z
                        c=f*z
                        s=h*z
                    end if
                    f= (c*g)+(s*y)
                    x=-(s*g)+(c*y)
                    tempm(1:m)=a(1:m,j)
                    a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
                    a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
                end do
                rv1(l)=0.0
                rv1(k)=f
                w(k)=x
            end do
        end do
    end subroutine svdcmp_dp

    ! imported from numerical recipes
    ! SVD-based least-squares fit for linear univariate model
    ! Given a set of N data points x,y with individual standard deviations sig, all arrays of length
    ! N, use chi^2 minimization to determine the M coefficients a of a function that depends linearly
    ! on a, y=sum_{i=1}^M a_i * afunc_i(x). Here we solve the fitting equations using singular value
    ! decomposition of the NxM matrix, as in § 2.6. On output, the MxM array v and the
    ! vector w of length M define part of the singular value decomposition, and can be used to
    ! obtain the covariance matrix. The program returns values for the M fit parameters a, and
    ! chi^2, chisq. The user supplies a subroutine funcs(x,afunc) that returns the M basis
    ! functions evaluated at x=X in the array afunc.
    subroutine svdfit_sp(x,y,sig,a,v,w,chisq,funcs)
        implicit none
        real(sp), dimension(:),   intent(in)  :: x,y,sig
        real(sp), dimension(:),   intent(out) :: a,w
        real(sp), dimension(:,:), intent(out) :: v
        real(sp),                 intent(out) :: chisq
        interface
            function funcs(x,n)
                implicit none
                real,    intent(in) :: x
                integer, intent(in) :: n
                real, dimension(n) :: funcs
            end function funcs
        end interface
        real(sp), parameter :: TOL=1.0e-5_sp
        integer                              :: i,ma,n
        real(sp), dimension(size(x))         :: b,sigi
        real(sp), dimension(size(x),size(a)) :: u,usav
        n=assert_eq(size(x),size(y),size(sig),'svdfit: n')
        ma=assert_eq(size(a),size(v,1),size(v,2),size(w),'svdfit: ma')
        sigi=1.0_sp/sig                       ! Accumulate coefficients of the fitting matrix in u.
        b=y*sigi
        do i=1,n
            usav(i,:)=funcs(x(i),ma)
        end do
        u=usav*spread(sigi,dim=2,ncopies=ma)
        usav=u
        call svdcmp(u,w,v)                   ! Singular value decomposition.
        where (w < TOL*maxval(w)) w=0.0      ! Edit the singular values, given TOL from the parameter statement.
        call svbksb(u,w,v,b,a)
        chisq=vabs(matmul(usav,a)-b)**2      ! Evaluate chi-square.
    end subroutine svdfit_sp

    ! imported from numerical recipes
    ! SVD-based least-squares fit for linear univariate model in double precision
    ! Given a set of N data points x,y with individual standard deviations sig, all arrays of length
    ! N, use chi^2 minimization to determine the M coefficients a of a function that depends linearly
    ! on a, y=sum_{i=1}^M a_i * afunc_i(x). Here we solve the fitting equations using singular value
    ! decomposition of the NxM matrix, as in § 2.6. On output, the MxM array v and the
    ! vector w of length M define part of the singular value decomposition, and can be used to
    ! obtain the covariance matrix. The program returns values for the M fit parameters a, and
    ! chi^2, chisq. The user supplies a subroutine funcs(x,afunc) that returns the M basis
    ! functions evaluated at x=X in the array afunc.
    subroutine svdfit_dp(x,y,sig,a,v,w,chisq,funcs)
        implicit none
        real(dp), dimension(:),   intent(in)  :: x,y,sig
        real(dp), dimension(:),   intent(out) :: a,w
        real(dp), dimension(:,:), intent(out) :: v
        real(dp),                 intent(out) :: chisq
        interface
            function funcs(x,n)
                implicit none
                real(kind=8),   intent(in) :: x
                integer,        intent(in) :: n
                real(kind=8), dimension(n) :: funcs
            end function funcs
        end interface
        real(dp), parameter :: TOL=1.0e-14_dp
        integer                              :: i,ma,n
        real(dp), dimension(size(x))         :: b,sigi
        real(dp), dimension(size(x),size(a)) :: u,usav
        n=assert_eq(size(x),size(y),size(sig),'svdfit: n')
        ma=assert_eq(size(a),size(v,1),size(v,2),size(w),'svdfit: ma')
        sigi=1.0_sp/sig                       ! Accumulate coefficients of the fitting matrix in u.
        b=y*sigi
        do i=1,n
            usav(i,:)=funcs(x(i),ma)
        end do
        u=usav*spread(sigi,dim=2,ncopies=ma)
        usav=u
        call svdcmp(u,w,v)                   ! Singular value decomposition.
        where (w < TOL*maxval(w)) w=0.0      ! Edit the singular values, given TOL from the parameter statement.
        call svbksb(u,w,v,b,a)
        chisq=vabs(matmul(usav,a)-b)**2      ! Evaluate chi-square.
    end subroutine svdfit_dp

    ! imported from numerical recipes
    ! modification of svdfit to support multivariate linear model, single precision
    ! Given a set of N data points x,y with individual standard deviations sig, all arrays of length
    ! N, use chi^2 minimization to determine the M coefficients a of a function that depends linearly
    ! on a, y=sum_{i=1}^M a_i * afunc_i(x). Here we solve the fitting equations using singular value
    ! decomposition of the NxM matrix, as in § 2.6. On output, the MxM array v and the
    ! vector w of length M define part of the singular value decomposition, and can be used to
    ! obtain the covariance matrix. The program returns values for the M fit parameters a, and
    ! chi^2, chisq. The user supplies a subroutine funcs(x,afunc) that returns the M basis
    ! functions evaluated at x=X in the array afunc.
    subroutine svd_multifit_sp(x,y,sig,a,v,w,chisq,funcs)
        implicit none
        real(sp), dimension(:,:), intent(in)  :: x
        real(sp), dimension(:),   intent(in)  :: y,sig
        real(sp), dimension(:),   intent(out) :: a,w
        real(sp), dimension(:,:), intent(out) :: v
        real(sp),                 intent(out) :: chisq
        interface
            function funcs(x,n)
                implicit none
                real,     intent(in) :: x(:)
                integer,  intent(in) :: n
                real, dimension(n) :: funcs
            end function funcs
        end interface
        real(sp), parameter :: TOL=1.0e-5_sp
        integer                              :: i,ma,n
        real(sp), dimension(size(sig))       :: b,sigi
        real(sp), dimension(size(x,dim=2),size(a)) :: u,usav
        n=assert_eq(size(x,dim=2),size(y),size(sig),'svd_multifit_sp: n')
        ma=assert_eq(size(a),size(v,1),size(v,2),size(w),'svd_multifit_sp: ma')
        sigi=1.0_sp/sig                       ! Accumulate coefficients of the fitting matrix in u.
        b=y*sigi
        do i=1,n
            usav(i,:)=funcs(x(:,i),ma)
        end do
        u=usav*spread(sigi,dim=2,ncopies=ma)
        usav=u
        call svdcmp(u,w,v)                   ! Singular value decomposition.
        where (w < TOL*maxval(w)) w=0.0      ! Edit the singular values, given TOL from the parameter statement.
        call svbksb(u,w,v,b,a)
        chisq=vabs(matmul(usav,a)-b)**2      ! Evaluate chi-square.
    end subroutine svd_multifit_sp

    ! imported from numerical recipes
    ! modification of svdfit to support multivariate linear model, double precision
    ! Given a set of N data points x,y with individual standard deviations sig, all arrays of length
    ! N, use chi^2 minimization to determine the M coefficients a of a function that depends linearly
    ! on a, y=sum_{i=1}^M a_i * afunc_i(x). Here we solve the fitting equations using singular value
    ! decomposition of the NxM matrix, as in § 2.6. On output, the MxM array v and the
    ! vector w of length M define part of the singular value decomposition, and can be used to
    ! obtain the covariance matrix. The program returns values for the M fit parameters a, and
    ! chi^2, chisq. The user supplies a subroutine funcs(x,afunc) that returns the M basis
    ! functions evaluated at x=X in the array afunc.
    subroutine svd_multifit_dp(x,y,sig,a,v,w,chisq,funcs)
        implicit none
        real(dp), dimension(:,:), intent(in)  :: x
        real(dp), dimension(:),   intent(in)  :: y,sig
        real(dp), dimension(:),   intent(out) :: a,w
        real(dp), dimension(:,:), intent(out) :: v
        real(dp),                 intent(out) :: chisq
        interface
            function funcs(x,n)
                implicit none
                real(kind=8), intent(in)   :: x(:)
                integer,      intent(in )  :: n
                real(kind=8), dimension(n) :: funcs
            end function funcs
        end interface
        real(dp), parameter :: TOL=1.0e-14_dp
        integer                              :: i,ma,n
        real(dp), dimension(size(sig))       :: b,sigi
        real(dp), dimension(size(x,dim=2),size(a)) :: u,usav
        n=assert_eq(size(x,dim=2),size(y),size(sig),'svd_multifit_dp: n')
        ma=assert_eq(size(a),size(v,1),size(v,2),size(w),'svd_multifit_dp: ma')
        sigi=1.0_dp/sig                       ! Accumulate coefficients of the fitting matrix in u.
        b=y*sigi
        do i=1,n
            usav(i,:)=funcs(x(:,i),ma)
        end do
        u=usav*spread(sigi,dim=2,ncopies=ma)
        usav=u
        call svdcmp(u,w,v)                   ! Singular value decomposition.
        where (w < TOL*maxval(w)) w=0.0_dp   ! Edit the singular values, given TOL from the parameter statement.
        call svbksb(u,w,v,b,a)
        chisq=vabs(matmul(usav,a)-b)**2      ! Evaluate chi-square.
    end subroutine svd_multifit_dp

    subroutine svdvar(v,w,cvm)
        implicit none
        real(sp), dimension(:,:), intent(in)  :: v
        real(sp), dimension(:),   intent(in)  :: w
        real(sp), dimension(:,:), intent(out) :: cvm
        ! To evaluate the covariance matrix cvm of the fit for M parameters obtained by svdfit,
        ! call this routine with matrices v,w as returned from svdfit. The dimensions are M for
        ! w and MxM for v and cvm.
        integer                      :: ma
        real(sp), dimension(size(w)) :: wti
        ma=assert_eq((/size(v,1),size(v,2),size(w),size(cvm,1),size(cvm,2)/),'svdvar')
        where (w /= 0.0)
            wti=1.0_sp/(w*w)
        elsewhere
            wti=0.0
        end where
        cvm=v*spread(wti,dim=1,ncopies=ma)
        cvm=matmul(cvm,transpose(v))          ! Covariance matrix is given by (15.4.20).
    end subroutine svdvar

    ! TRACE calculates the trace of a real 2D matrix
    pure function trace(mat) result (tr)
        real, intent(in) :: mat(:,:)
        real             :: tr
        integer          :: i
        tr = 0.
        do i = 1, size(mat, dim = 1)
            tr = tr + mat(i,i)
        enddo
    end function trace

    ! Return the length (ordinary L2 norm) of a vector.
    function vabs_sp(v)
        implicit none
        real(sp), dimension(:), intent(in) :: v
        real(sp) :: vabs_sp
        vabs_sp=sqrt(dot_product(v,v))
    end function vabs_sp

    ! Return the length (ordinary L2 norm) of a vector.
    function vabs_dp(v)
        implicit none
        real(dp), dimension(:), intent(in) :: v
        real(dp) :: vabs_dp
        vabs_dp=sqrt(dot_product(v,v))
    end function vabs_dp

    !> Vector angle between two non-zero vectors in R^n
    !! magnitude of v and w must equal 1
    !! \theta = \arccos \frac {\mathbf v \cdot \mathbf w}{\left\Vert{\mathbf v}\right\Vert \left\Vert{\mathbf w}\right\Vert}
    !! by definition of the cosine formula for dot product.
    !! This function assumes the input vectors are normals!!
    pure function vector_angle_norm(v, w) result(eulerdist)
        real, dimension(3), intent(in) :: v,w
        real :: eulerdist
        eulerdist = acos( dot_product(v, w) )
    end function vector_angle_norm

    ! trigonometry

    !>   returns acos with the argument's absolute value limited to 1.
    !!         this seems to be necessary due to small numerical inaccuracies.
    pure function myacos_sp( arg ) result( r )
        real, intent(in) :: arg     !< input (radians)
        real             :: r, x, y
        x = min(1.,abs(arg))
        y = sign(x,arg)
        r = acos(y)
    end function myacos_sp

    !>   returns acos with the argument's absolute value limited to 1.
    !!         this seems to be necessary due to small numerical inaccuracies.
    pure function myacos_dp( arg ) result( r )
        real(dp), intent(in) :: arg     !< input (radians)
        real(dp)             :: r, x, y
        x = min(1._dp,abs(arg))
        y = sign(x,arg)
        r = acos(y)
    end function myacos_dp

    !>   converts between radians and degrees
    elemental function deg2rad_sp( deg ) result( rad )
        real, intent(in) :: deg  !< angle (degrees)
        real             :: rad  !< angle (radians)
        rad = (deg/180.)*pi
    end function deg2rad_sp

    !>   converts between radians and degrees
    elemental function deg2rad_dp( deg ) result( rad )
        real(dp), intent(in) :: deg  !< angle (degrees)
        real(dp)             :: rad  !< angle (radians)
        rad = (deg/180._dp)*dpi
    end function deg2rad_dp

    !>   converts from radians to degrees
    elemental function rad2deg_1( rad ) result( deg )
        real(sp), intent(in) :: rad  !< angle (radians)
        real(sp)             :: deg  !< angle (degrees)
        deg = (rad/PI)*180.
    end function rad2deg_1

    !>   converts from radians to degrees
    elemental function rad2deg_2( rad ) result( deg )
        real(dp), intent(in) :: rad  !< angle (radians)
        real(dp)             :: deg  !< angle (degrees)
        deg = (rad/DPI)*180.d0
    end function rad2deg_2

    !>    is for calculating the radius
    pure function hyp( x1, x2, x3 ) result( h )
        real, intent(in) :: x1, x2
        real, intent(in), optional :: x3
        real :: h
        if( present(x3) )then
            h = sqrt(x1**2.+x2**2.+x3**2.)
        else
            h = sqrt(x1**2.+x2**2.)
        endif
    end function hyp

    !>   calculates the euclidean distance between two vectors of dimension _n_
    pure function euclid( vec1, vec2 ) result( dist )
        real, intent(in)    :: vec1(:), vec2(:)
        real                :: dist
        dist = sqrt(sum((vec1-vec2)**2))
    end function euclid

    ! imported from numerical recipes
    function pythag_sp(a,b)
        implicit none
        real(sp), intent(in) :: a,b
        real(sp) :: pythag_sp
        ! Computes (a^2+b^2)^(1/2) without destructive underflow or overflow.
        real(sp) :: absa,absb
        absa=abs(a)
        absb=abs(b)
        if (absa > absb) then
            pythag_sp=absa*sqrt(1.0_sp+(absb/absa)**2)
        else
            if (absb == 0.0) then
                pythag_sp=0.0
            else
                pythag_sp=absb*sqrt(1.0_sp+(absa/absb)**2)
            end if
        end if
    end function pythag_sp

    ! imported from numerical recipes
    function pythag_dp(a,b)
        implicit none
        real(dp), intent(in) :: a,b
        real(dp) :: pythag_dp
        real(dp) :: absa,absb
        absa=abs(a)
        absb=abs(b)
        if (absa > absb) then
            pythag_dp=absa*sqrt(1.0_dp+(absb/absa)**2)
        else
            if (absb == 0.0) then
                pythag_dp=0.0
            else
                pythag_dp=absb*sqrt(1.0_dp+(absa/absb)**2)
            end if
        end if
    end function pythag_dp

end module simple_linalg
