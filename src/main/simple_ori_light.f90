! an orientation
module simple_ori_light
use simple_defs
implicit none

private
public :: ori_light

type :: ori_light
contains
    procedure, nopass, private :: euler2dm_sp
    procedure, nopass, private :: euler2dm_dp
    generic :: euler2dm => euler2dm_sp, euler2dm_dp
    procedure, nopass, private :: euler2m_sp
    procedure, nopass, private :: euler2m_dp
    generic :: euler2m => euler2m_sp, euler2m_dp
end type ori_light

contains

    !> \brief  calculates the derivatives of the rotation matrices w.r.t. one Euler angle
    !! \param e1,e2,e3 Euler triplet
    subroutine euler2dm_sp( euls, drmat )
        real, intent(in)     :: euls(3)
        real, intent(out)    :: drmat(3,3,3)
        real, dimension(3,3) :: r1, r2, r3, dr1, dr2, dr3
        real                 :: e1, e2, e3
        e1 = euls(1)
        e2 = euls(2)
        e3 = euls(3)
        r1  =  rotmat_sp(e1,3) ! rotation around z
        r2  =  rotmat_sp(e2,2) ! tilt
        r3  =  rotmat_sp(e3,3) ! rotation around z
        dr1 = drotmat_sp(e1,3) ! derivative of r1 w.r.t. e1
        dr2 = drotmat_sp(e2,2) ! derivative of r2 w.r.t. e2
        dr3 = drotmat_sp(e3,3) ! derivative of r3 w.r.t. e3
        drmat(:,:,1) = matmul(matmul(r3,r2),dr1)
        drmat(:,:,2) = matmul(matmul(r3,dr2),r1)
        drmat(:,:,3) = matmul(matmul(dr3,r2),r1)
    contains
        !>  \brief  returns the derivative of the rotation matrix for _ang_ degrees of rotation
        !! around x,y or z for _choice_ = _1_,_2_ or _3_
        pure function drotmat_sp( ang, choice ) result( r )
            real,    intent(in) :: ang
            integer, intent(in) :: choice
            real :: r(3,3)
            real :: ang_in_rad
            ang_in_rad = ang*pi/180.
            if ( choice == 1 ) then
                r( 1,1 ) = 0.
                r( 1,2 ) = 0.
                r( 1,3 ) = 0.
                r( 2,1 ) = 0.
                r( 2,2 ) = -sin( ang_in_rad )
                r( 2,3 ) = -cos( ang_in_rad )
                r( 3,1 ) = 0.
                r( 3,2 ) =  cos( ang_in_rad )
                r( 3,3 ) = -sin( ang_in_rad )
            elseif ( choice == 2 ) then
                r( 1,1 ) = -sin( ang_in_rad )
                r( 1,2 ) = 0.
                r( 1,3 ) = -cos( ang_in_rad )
                r( 2,1 ) = 0.
                r( 2,2 ) = 0.
                r( 2,3 ) = 0.
                r( 3,1 ) =  cos( ang_in_rad )
                r( 3,2 ) = 0.
                r( 3,3 ) = -sin( ang_in_rad )
            elseif ( choice == 3 ) then
                r( 1,1 ) = -sin( ang_in_rad )
                r( 1,2 ) =  cos( ang_in_rad )
                r( 1,3 ) = 0.
                r( 2,1 ) = -cos( ang_in_rad )
                r( 2,2 ) = -sin( ang_in_rad )
                r( 2,3 ) = 0.
                r( 3,1 ) = 0.
                r( 3,2 ) = 0.
                r( 3,3 ) = 0.
                ! beware of the signs:z-rot is really negative
            endif
            r = r * pi / 180.
        end function drotmat_sp
    end subroutine euler2dm_sp

    !> \brief  calculates the derivatives of the rotation matrices w.r.t. one Euler angle
    !! \param e1,e2,e3 Euler triplet
    subroutine euler2dm_dp( euls, drmat )
        real(dp), intent(in)     :: euls(3)
        real(dp), intent(out)    :: drmat(3,3,3)
        real(dp), dimension(3,3) :: r1, r2, r3, dr1, dr2, dr3
        real(dp)                 :: e1, e2, e3
        e1 = euls(1)
        e2 = euls(2)
        e3 = euls(3)
        r1  =  rotmat_dp(e1,3) ! rotation around z
        r2  =  rotmat_dp(e2,2) ! tilt
        r3  =  rotmat_dp(e3,3) ! rotation around z
        dr1 = drotmat_dp(e1,3) ! derivative of r1 w.r.t. e1
        dr2 = drotmat_dp(e2,2) ! derivative of r2 w.r.t. e2
        dr3 = drotmat_dp(e3,3) ! derivative of r3 w.r.t. e3
        drmat(:,:,1) = matmul(matmul(r3,r2),dr1)
        drmat(:,:,2) = matmul(matmul(r3,dr2),r1)
        drmat(:,:,3) = matmul(matmul(dr3,r2),r1)
    contains
        !>  \brief  returns the derivative of the rotation matrix for _ang_ degrees of rotation
        !! around x,y or z for _choice_ = _1_,_2_ or _3_
        pure function drotmat_dp( ang, choice ) result( r )
            real(dp), intent(in)           :: ang
            integer, intent(in)        :: choice
            real(dp) :: r(3,3)
            real(dp) :: ang_in_rad
            ang_in_rad = ang*dpi/180._dp
            if ( choice == 1 ) then
                r( 1,1 ) = 0._dp
                r( 1,2 ) = 0._dp
                r( 1,3 ) = 0._dp
                r( 2,1 ) = 0._dp
                r( 2,2 ) = -sin( ang_in_rad )
                r( 2,3 ) = -cos( ang_in_rad )
                r( 3,1 ) = 0._dp
                r( 3,2 ) =  cos( ang_in_rad )
                r( 3,3 ) = -sin( ang_in_rad )
            elseif ( choice == 2 ) then
                r( 1,1 ) = -sin( ang_in_rad )
                r( 1,2 ) = 0._dp
                r( 1,3 ) = -cos( ang_in_rad )
                r( 2,1 ) = 0._dp
                r( 2,2 ) = 0._dp
                r( 2,3 ) = 0._dp
                r( 3,1 ) =  cos( ang_in_rad )
                r( 3,2 ) = 0._dp
                r( 3,3 ) = -sin( ang_in_rad )
            elseif ( choice == 3 ) then
                r( 1,1 ) = -sin( ang_in_rad )
                r( 1,2 ) =  cos( ang_in_rad )
                r( 1,3 ) = 0._dp
                r( 2,1 ) = -cos( ang_in_rad )
                r( 2,2 ) = -sin( ang_in_rad )
                r( 2,3 ) = 0._dp
                r( 3,1 ) = 0._dp
                r( 3,2 ) = 0._dp
                r( 3,3 ) = 0._dp
                ! beware of the signs:z-rot is really negative
            endif
            r = r * dpi / 180._dp
        end function drotmat_dp
    end subroutine euler2dm_dp

    !>  \brief  makes a rotation matrix from a Spider format Euler triplet
    !! \param e1,e2,e3 Euler triplet
    pure function euler2m_sp( euls ) result( r )
        real, intent(in)     :: euls(3)
        real, dimension(3,3) :: r1, r2, r3, r, tmp
        real                 :: e1, e2, e3
        e1 = euls(1)
        e2 = euls(2)
        e3 = euls(3)
        r1 = rotmat_sp(e1,3) ! rotation around z
        r2 = rotmat_sp(e2,2) ! tilt
        r3 = rotmat_sp(e3,3) ! rotation around z
        ! order of multiplication is r3r2r1
        tmp = matmul(r3,r2)
        r = matmul(tmp,r1)
    end function euler2m_sp

    !>  \brief  makes a rotation matrix from a Spider format Euler triplet
    !! \param e1,e2,e3 Euler triplet
    pure function euler2m_dp( euls ) result( r )
        real(dp), intent(in)     :: euls(3)
        real(dp), dimension(3,3) :: r1, r2, r3, r, tmp
        real(dp)                 :: e1, e2, e3
        e1 = euls(1)
        e2 = euls(2)
        e3 = euls(3)
        r1 = rotmat_dp(e1,3) ! rotation around z
        r2 = rotmat_dp(e2,2) ! tilt
        r3 = rotmat_dp(e3,3) ! rotation around z
        ! order of multiplication is r3r2r1
        tmp = matmul(r3,r2)
        r = matmul(tmp,r1)
    end function euler2m_dp

    !>  \brief  returns the rotation matrix for _ang_ degrees of rotation
    !! around x,y or z for _choice_ = _1_,_2_ or _3_
    pure function rotmat_sp( ang, choice ) result( r )
        real, intent(in)           :: ang
        integer, intent(in)        :: choice
        real :: r(3,3)
        real :: ang_in_rad
        ang_in_rad = ang*pi/180.
        if ( choice == 1 ) then
            r( 1,1 ) = 1.
            r( 1,2 ) = 0.
            r( 1,3 ) = 0.
            r( 2,1 ) = 0.
            r( 2,2 ) = cos( ang_in_rad )
            r( 2,3 ) =-sin( ang_in_rad )
            r( 3,1 ) = 0.
            r( 3,2 ) = sin( ang_in_rad )
            r( 3,3 ) = cos( ang_in_rad )
        elseif ( choice == 2 ) then
            r( 1,1 ) = cos( ang_in_rad )
            r( 1,2 ) = 0.
            r( 1,3 ) = -sin( ang_in_rad )
            r( 2,1 ) = 0.
            r( 2,2 ) = 1.
            r( 2,3 ) = 0.
            r( 3,1 ) = sin( ang_in_rad )
            r( 3,2 ) = 0.
            r( 3,3 ) = cos( ang_in_rad )
        elseif ( choice == 3 ) then
            r( 1,1 ) = cos( ang_in_rad )
            r( 1,2 ) = sin( ang_in_rad )
            r( 1,3 ) = 0.
            r( 2,1 ) = -sin( ang_in_rad )
            r( 2,2 ) = cos( ang_in_rad )
            r( 2,3 ) = 0.
            r( 3,1 ) = 0.
            r( 3,2 ) = 0.
            r( 3,3 ) = 1.
            ! beware of the signs:z-rot is really negative
        endif
    end function rotmat_sp

    !>  \brief  returns the rotation matrix for _ang_ degrees of rotation
    !! around x,y or z for _choice_ = _1_,_2_ or _3_
    pure function rotmat_dp( ang, choice ) result( r )
        real(dp), intent(in) :: ang
        integer,  intent(in) :: choice
        real(dp) :: r(3,3)
        real(dp) :: ang_in_rad
        ang_in_rad = ang*dpi/180._dp
        if ( choice == 1 ) then
            r( 1,1 ) = 1._dp
            r( 1,2 ) = 0._dp
            r( 1,3 ) = 0._dp
            r( 2,1 ) = 0._dp
            r( 2,2 ) = cos( ang_in_rad )
            r( 2,3 ) =-sin( ang_in_rad )
            r( 3,1 ) = 0._dp
            r( 3,2 ) = sin( ang_in_rad )
            r( 3,3 ) = cos( ang_in_rad )
        elseif ( choice == 2 ) then
            r( 1,1 ) = cos( ang_in_rad )
            r( 1,2 ) = 0._dp
            r( 1,3 ) = -sin( ang_in_rad )
            r( 2,1 ) = 0._dp
            r( 2,2 ) = 1._dp
            r( 2,3 ) = 0._dp
            r( 3,1 ) = sin( ang_in_rad )
            r( 3,2 ) = 0._dp
            r( 3,3 ) = cos( ang_in_rad )
        elseif ( choice == 3 ) then
            r( 1,1 ) = cos( ang_in_rad )
            r( 1,2 ) = sin( ang_in_rad )
            r( 1,3 ) = 0._dp
            r( 2,1 ) = -sin( ang_in_rad )
            r( 2,2 ) = cos( ang_in_rad )
            r( 2,3 ) = 0._dp
            r( 3,1 ) = 0._dp
            r( 3,2 ) = 0._dp
            r( 3,3 ) = 1._dp
            ! beware of the signs:z-rot is really negative
        endif
    end function rotmat_dp

end module simple_ori_light
