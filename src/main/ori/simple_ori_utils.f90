!@descr: Euler angle and rotation matrix utilities for the ori and oris classes
module simple_ori_utils
use simple_defs
use simple_linalg, only: rad2deg, deg2rad, myacos
use simple_math,   only: rotmat2d
use simple_rnd,    only: ran3
implicit none

public :: euler2m, m2euler, dm2euler, euler_dist, euler_inplrotdist, euler_compose, euler_mirror
public :: make_transfmat, rnd_romat, geodesic_frobdev, rotmat2d
public :: genfire2euler, euler2genfire
private

real, parameter :: zvec(3) = [0.,0.,1.]

contains

    !>  \brief  makes a rotation matrix from a Spider format Euler triplet
    pure function euler2m( euls ) result( r )
        real, intent(in)     :: euls(3)
        real, dimension(3,3) :: r1, r2, r3, r, tmp
        r1 = rotmat(euls(1),3) ! rotation around z
        r2 = rotmat(euls(2),2) ! tilt
        r3 = rotmat(euls(3),3) ! rotation around z
        ! order of multiplication is r3r2r1
        tmp = matmul(r3,r2)
        r   = matmul(tmp,r1)
    end function euler2m

    ! Intrinsic zy'z" convention
    ! [r11 r21 r31]   [...         ...         sΘcΦ]
    ! [r12 r22 r32] = [cΘcΨsΦ+cΦsΨ cΦcΨ-cΘsΦsΨ sΘsΦ]
    ! [r13 r23 r33]   [-sΘcΨ       sΘsΨ        cΘ  ]
    pure function m2euler( R ) result( e )
        real, intent(in)  :: R(3,3)
        real :: c, e(3)
        ! clamp to protect acos from tiny numerical drift
        c = max(-1.0, min(1.0, R(3,3)))
        if( c < 0.9999999 )then
            if( c > -0.9999999 )then
                e = rad2deg([atan2(R(3,2),R(3,1)), acos(c), atan2(R(2,3),-R(1,3))])
            else
                e = [rad2deg(-atan2(R(1,2),R(2,2))), 180., 0.]
            endif
        else
            e = [rad2deg(atan2(R(1,2),R(2,2))), 0., 0.]
        endif
        do while( e(1) < 0. )
            e(1) = e(1) + 360.
        end do
        do while( e(3) < 0. )
            e(3) = e(3) + 360.
        end do
    end function m2euler

    ! Double precision version of m2euler above
    pure function dm2euler( R ) result( e )
        real(dp), intent(in)  :: R(3,3)
        real(dp) :: c, e(3)
        c = max(-1.d0, min(1.d0, R(3,3)))
        if( c < 0.9999999d0 )then
            if( c > -0.9999999d0 )then
                e = rad2deg([atan2(R(3,2),R(3,1)), acos(c), atan2(R(2,3),-R(1,3))])
            else
                e = [rad2deg(-atan2(R(1,2),R(2,2))), 180.d0, 0.d0]
            endif
        else
            e = [rad2deg(atan2(R(1,2),R(2,2))), 0.d0, 0.d0]
        endif
        do while( e(1) < 0.d0 )
            e(1) = e(1) + 360.d0
        end do
        do while( e(3) < 0.d0 )
            e(3) = e(3) + 360.d0
        end do
    end function dm2euler

    !>  in-plane parameters to 3x3 transformation matrix
    function make_transfmat( psi, tx, ty ) result( R )
        real,intent(in) :: psi,tx,ty
        real            :: R(3,3),radpsi,cospsi,sinpsi
        radpsi = deg2rad( psi )
        cospsi = cos( radpsi )
        sinpsi = sin( radpsi )
        R      = 0.
        R(1,1) = cospsi
        R(2,2) = cospsi
        R(1,2) = sinpsi
        R(2,1) = -sinpsi
        R(1,3) = tx
        R(2,3) = ty
        R(3,3) = 1.
    end function make_transfmat

    !>  \brief  returns the rotation matrix for _ang_ degrees of rotation
    !! around x,y or z for _choice_ = _1_,_2_ or _3_
    pure function rotmat( ang, choice ) result( r )
        real, intent(in)           :: ang
        integer, intent(in)        :: choice
        real :: r(3,3)
        real :: ang_in_rad, cosang, sinang
        ang_in_rad = deg2rad(ang)
        cosang = cos( ang_in_rad )
        sinang = sin( ang_in_rad )
        if ( choice == 1 ) then
            r( 1,1 ) = 1.
            r( 1,2 ) = 0.
            r( 1,3 ) = 0.
            r( 2,1 ) = 0.
            r( 2,2 ) = cosang
            r( 2,3 ) =-sinang
            r( 3,1 ) = 0.
            r( 3,2 ) = sinang
            r( 3,3 ) = cosang
        elseif ( choice == 2 ) then
            r( 1,1 ) = cosang
            r( 1,2 ) = 0.
            r( 1,3 ) = -sinang
            r( 2,1 ) = 0.
            r( 2,2 ) = 1.
            r( 2,3 ) = 0.
            r( 3,1 ) = sinang
            r( 3,2 ) = 0.
            r( 3,3 ) = cosang
        elseif ( choice == 3 ) then
            r( 1,1 ) = cosang
            r( 1,2 ) = sinang
            r( 1,3 ) = 0.
            r( 2,1 ) = -sinang
            r( 2,2 ) = cosang
            r( 2,3 ) = 0.
            r( 3,1 ) = 0.
            r( 3,2 ) = 0.
            r( 3,3 ) = 1.
            ! beware of the signs:z-rot is really negative
        endif
    end function rotmat

    !> GENFIRE/RESIRE projection angles (phi, theta, psi in degrees; rotMAT = Rz(psi) Ry(theta) Rz(phi),
    !! active rotations of column vectors) to SIMPLE Euler angles. With image x the first array dimension
    !! of a projection, image x is rotMAT(:,3) and image y rotMAT(:,1); the volume frames coincide
    pure function genfire2euler( ang ) result( e )
        real, intent(in) :: ang(3)
        real     :: e(3)
        real(dp) :: m(3,3), r(3,3)
        m      = genfire_rotmat(ang)
        r(1,:) = m(:,3)
        r(2,:) = m(:,1)
        r(3,:) = cross3(r(1,:), r(2,:))
        e      = real(dm2euler(r))
    end function genfire2euler

    !> SIMPLE Euler angles to GENFIRE/RESIRE (phi, theta, psi in degrees), the inverse of genfire2euler;
    !! phi and psi in (-180,180], theta in [0,180], psi = 0 where theta is 0 or 180
    pure function euler2genfire( e ) result( ang )
        real, intent(in) :: e(3)
        real     :: ang(3)
        real(dp) :: r(3,3), m(3,3), c, th
        r      = real(euler2m(e), dp)
        m(:,3) = r(1,:)
        m(:,1) = r(2,:)
        m(:,2) = cross3(m(:,3), m(:,1))
        c  = max(-1._dp, min(1._dp, m(3,3)))
        th = acos(c)
        if( sin(th) > 1.e-6_dp )then
            ang = real([atan2(m(3,2), -m(3,1)), th, atan2(m(2,3), m(1,3))] * 180._dp / DPI)
        else if( c > 0._dp )then
            ang = real([atan2(m(2,1), m(1,1)), 0._dp, 0._dp] * 180._dp / DPI)
        else
            ang = real([atan2(m(1,2), -m(1,1)), DPI, 0._dp] * 180._dp / DPI)
        endif
    end function euler2genfire

    !> rotMAT of GENFIRE/RESIRE: Rz(psi) Ry(theta) Rz(phi), standard active rotations
    pure function genfire_rotmat( ang ) result( m )
        real, intent(in) :: ang(3)
        real(dp) :: m(3,3), a(3), c(3), s(3), rz1(3,3), ry(3,3), rz2(3,3)
        a   = real(ang, dp) * DPI / 180._dp
        c   = cos(a)
        s   = sin(a)
        rz1 = 0._dp
        rz1(1,1) = c(1); rz1(1,2) = -s(1); rz1(2,1) = s(1); rz1(2,2) = c(1); rz1(3,3) = 1._dp
        ry  = 0._dp
        ry(1,1)  = c(2); ry(1,3)  =  s(2); ry(2,2)  = 1._dp; ry(3,1) = -s(2); ry(3,3)  = c(2)
        rz2 = 0._dp
        rz2(1,1) = c(3); rz2(1,2) = -s(3); rz2(2,1) = s(3); rz2(2,2) = c(3); rz2(3,3) = 1._dp
        m = matmul(rz2, matmul(ry, rz1))
    end function genfire_rotmat

    pure function cross3( a, b ) result( c )
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross3

    !>  \brief  for generating a random rotation matrix
    subroutine rnd_romat( rmat )
        ! Fast Random Rotation Matrices, Arvo, Graphics Gems III, 1992
        real, intent(out) :: rmat(3,3)
        real :: theta, phi, z, vx, vy, vz
        real :: r, st, ct, sx, sy
        ! init
        theta = ran3()*TWOPI
        phi   = ran3()*TWOPI
        z     = ran3()*2.
        ! V
        r  = sqrt( z )
        vx = r*sin( phi )
        vy = r*cos( phi )
        vz = sqrt( 2.-z )
        ! S=Vt*R; sz=vz
        st = sin( theta )
        ct = cos( theta )
        sx = vx*ct - vy*st
        sy = vx*st + vy*ct
        ! M
        rmat(1,1) = vx * sx - ct
        rmat(1,2) = vx * sy - st
        rmat(1,3) = vx * vz
        rmat(2,1) = vy * sx + st
        rmat(2,2) = vy * sy - ct
        rmat(2,3) = vy * vz
        rmat(3,1) = vz * sx
        rmat(3,2) = vz * sy
        rmat(3,3) = 1. - z
    end subroutine rnd_romat

    !>  \brief  calculates the distance (in radians) btw two Euler triplets
    pure real function euler_dist( euls1, euls2 )
        real, intent(in) :: euls1(3), euls2(3)
        real             :: normal1(3), normal2(3)
        normal1 = euler_normal(euls1)
        normal2 = euler_normal(euls2)
        euler_dist = myacos(dot_product(normal1, normal2))
    end function euler_dist

    !>  \brief  calculates the normal vector of a Euler triplet
    pure function euler_normal( euls ) result( normal )
        real, intent(in) :: euls(3)
        real             :: normal(3), rmat(3,3)
        rmat   = euler2m(euls)
        normal = matmul(zvec, rmat)
    end function euler_normal

    !>  \brief  calculates the in-plane distance of a euler triplet (so just psi) in radians
    pure real function euler_inplrotdist( euls1, euls2 )
        real, intent(in) :: euls1(3), euls2(3)
        real, parameter  :: u(2) = [0.0, 1.0]
        real             :: mat(2,2), x1(2), x2(2)
        call rotmat2d(euls1(3), mat)
        x1   = matmul(u,mat)
        call rotmat2d(euls2(3), mat)
        x2   = matmul(u,mat)
        euler_inplrotdist = myacos(dot_product(x1,x2))
    end function euler_inplrotdist

    !>  \brief  is for composing Euler triplets
    pure subroutine euler_compose( euls1, euls2, euls_out )
        real, intent(in)    :: euls1(3), euls2(3)
        real, intent(inout) :: euls_out(3)
        real                :: rmat(3,3), rmat1(3,3), rmat2(3,3)
        rmat1 = euler2m(euls1)
        rmat2 = euler2m(euls2)
        ! the composed rotation matrix is constructed
        rmat = matmul(rmat2,rmat1)  ! multiplication of two rotation matrices commute (n>2 not)
        ! convert to euler
        euls_out = m2euler(rmat)
    end subroutine euler_compose

    !>  \brief  Generates mirror euler angles
    pure subroutine euler_mirror( euls, euls_out )
        real, intent(in)    :: euls(3)
        real, intent(inout) :: euls_out(3)
        real :: rmat(3,3)
        euls_out(1) = euls(1)
        euls_out(2) = 180. + euls(2)
        euls_out(3) = 180. - euls(3)
        ! the mirrored rotation matrix is constructed
        rmat = euler2m(euls_out)
        ! convert to euler
        euls_out = m2euler(rmat)
    end subroutine euler_mirror

    !>  \brief  this metric is measuring the frobenius deviation from the identity matrix .in.[0,2*sqrt(2)]
    !!          Larochelle, P.M., Murray, A.P., Angeles, J., A distance metric for finite
    !!          sets of rigid-body displacement in the polar decomposition. ASME J. Mech. Des.
    !!          129, 883--886 (2007)
    pure real function geodesic_frobdev( euls1, euls2 )
        real, intent(in) :: euls1(3), euls2(3)
        real :: Imat(3,3), sumsq, diffmat(3,3)
        real :: rmat1(3,3), rmat2(3,3)
        Imat      = 0.
        Imat(1,1) = 1.
        Imat(2,2) = 1.
        Imat(3,3) = 1.
        rmat1 = euler2m(euls1)
        rmat2 = euler2m(euls2)
        diffmat = Imat - matmul(rmat1, transpose(rmat2))
        sumsq   = sum(diffmat * diffmat)
        if( sumsq > 1e-6 )then
            geodesic_frobdev = sqrt(sumsq)
        else
            geodesic_frobdev = 0.
        endif
    end function geodesic_frobdev

end module simple_ori_utils
