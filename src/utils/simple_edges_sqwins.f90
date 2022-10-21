module simple_edges_sqwins
use simple_defs
use simple_error, only: simple_exception
implicit none

interface cosedge
    module procedure cosedge_1, cosedge_2
end interface

interface cosedge_inner
    module procedure cosedge_inner_1, cosedge_inner_2
end interface

interface hardedge
    module procedure hardedge_1, hardedge_2, hardedge_3, hardedge_4 
end interface

interface hardedge_inner
    module procedure hardedge_inner_1, hardedge_inner_2, hardedge_inner_3, hardedge_inner_4 
end interface

interface sqwin_1d
    module procedure sqwin_1d_1, sqwin_1d_2
end interface

interface sqwin_2d
    module procedure sqwin_2d_1, sqwin_2d_2
end interface

interface sqwin_3d
    module procedure sqwin_3d_1, sqwin_3d_2
end interface

contains

    ! edge functions

    !>   two-dimensional hard edge
    !! \f$r^2 < x^2+y^2\f$.
    !! \param x x position
    !! \param y y position
    !! \param mskrad masking radius
    !! \return w on or off
    !!
    pure function hardedge_1( x, y, mskrad ) result( w )
        real,intent(in) :: x, y, mskrad
        real :: w
        w = 1.
        if( x * x + y * y > mskrad * mskrad ) w = 0.
    end function hardedge_1

    !>   three-dimensional hard edge
    !! \f$r^2 < x^2+y^2+z^2\f$.
    !! \param x x position
    !! \param y y position
    !! \param mskrad masking radius
    !! \return w on or off
    !!
    pure function hardedge_2( x, y, z, mskrad ) result( w )
        real,intent(in) :: x, y, z, mskrad
        real :: w
        w = 1.
        if( x * x + y * y + z * z > mskrad * mskrad ) w = 0.
    end function hardedge_2

   pure function hardedge_3( x, y, mskrad ) result( w )
        integer,intent(in) :: x, y
        real, intent(in)   :: mskrad
        real :: w
        w = 1.
        if( real(x * x + y * y) > mskrad * mskrad ) w = 0.
    end function hardedge_3
    !!
    pure function hardedge_4( x, y, z, mskrad ) result( w )
        integer,intent(in) :: x, y, z
        real, intent(in)   :: mskrad
        real :: w
        w = 1.
        if( real(x * x + y * y + z * z) > mskrad * mskrad ) w = 0.
    end function hardedge_4

    !>   two-dimensional hard edge
    !! \f$r < \sqrt{x^2+y^2}\f$.
    !! \return w on or off
    !!
    pure function hardedge_inner_1( x, y, mskrad ) result( w )
        real,intent(in) :: x, y, mskrad
        real :: w
        w = 0.
        if( x * x + y * y > mskrad * mskrad ) w = 1.
    end function hardedge_inner_1

    pure function hardedge_inner_2( x, y, mskrad ) result( w )
        integer,intent(in) :: x, y
        real, intent(in) ::mskrad
        real :: w
        w = 0.
        if( real(x * x + y * y) > mskrad * mskrad ) w = 1.
    end function hardedge_inner_2

    !>   three-dimensional hard edge
    pure function hardedge_inner_3( x, y, z, mskrad ) result( w )
        real,intent(in) :: x, y, z, mskrad
        real :: w
        w = 0.
        if( x * x + y * y + z * z > mskrad * mskrad ) w = 1.
    end function hardedge_inner_3

    pure function hardedge_inner_4( x, y, z, mskrad ) result( w )
        integer,intent(in) :: x, y,z
        real, intent(in) ::mskrad
        real :: w
        w = 0.
        if( real(x * x + y * y + z * z) > mskrad * mskrad ) w = 1.
    end function hardedge_inner_4

    !>   two-dimensional gaussian edge
    !! \param x x position
    !! \param y y position
   pure function cosedge_1( x, y, box, mskrad ) result( w )
       real, intent(in)    :: x, y     !< input points
       integer, intent(in) :: box      !< window size
       real, intent(in)    :: mskrad   !< mask radius
        real                :: w, rad, width, maxrad
        maxrad = real(box/2)
        rad    = sqrt(x**2.+y**2.)
        width  = 2.*(maxrad-mskrad)
        w      = 1.
        if( rad .ge. maxrad )then
            w = 0.
        else if( rad .ge. (maxrad-width) )then
            w = (cos(((rad-(maxrad-width))/width)*pi)+1.)/2.
        endif
    end function cosedge_1

    !>   three-dimensional gaussian edge
    !! \f$r = \cos{(1+{(\pi{r - (d+2m)/(d-2m)})})}\f$.
    !! \param x x position
    !! \param y y position
    pure function cosedge_2( x, y, z, box, mskrad ) result( w )
        real, intent(in)    :: x, y, z   !< input points
        integer, intent(in) :: box       !< window size
        real, intent(in)    :: mskrad    !< mask radius
        real                :: w, rad, maxrad, width
        maxrad = real(box/2)
        rad    = sqrt(x**2.+y**2.+z**2.)
        width  = 2.*(maxrad-mskrad)
        w      = 1.
        if( rad .ge. maxrad )then
            w = 0.
        else if( rad .ge. (maxrad-width) )then
            w = (cos(((rad-(maxrad-width))/width)*pi)+1.)/2.
        endif
    end function cosedge_2

    !> \brief  two-dimensional gaussian edge
    !! \param x x position
    !! \param y y position
    pure function cosedge_inner_1( x, y, width, mskrad ) result( w )
        real, intent(in) :: x, y, width  !< input points and width
        real, intent(in) :: mskrad       !< mask radius
        real             :: w, rad
        rad = sqrt(x**2.+y**2.)
        if( rad .lt. mskrad-width )then
            w = 0.
        else if( rad .gt. mskrad )then
            w = 1.
        else
            w = (1.+cos(pi*(mskrad-rad)/width))/2.0
        endif
    end function cosedge_inner_1

    !>   two-dimensional gaussian edge
    !! \param x x position
    !! \param y y position
    !! \param z z position
    pure function cosedge_inner_2( x, y, z, width, mskrad ) result( w )
        real, intent(in) :: x, y, z, width !< inner mask radius
        real, intent(in) :: mskrad !< mask radius
        real             :: w, rad
        rad = sqrt(x**2.+y**2.+z**2.)
        if( rad .lt. mskrad-width )then
            w = 0.
        else if( rad .gt. mskrad )then
            w = 1.
        else
            w = (1.+cos(pi*(mskrad-rad)/width))/2.0
        endif
    end function cosedge_inner_2

    !> one-dimensional symmetric hard window
    pure subroutine sqwin_1d_1( x, winsz, lowerlim, upperlim )
        real,    intent(in)  :: x                   !< input point
        real,    intent(in)  :: winsz               !< window size
        integer, intent(out) :: lowerlim, upperlim  !< window bounds
        integer :: iwinsz
        iwinsz   = ceiling(winsz - 0.5)
        lowerlim = nint(x)
        upperlim = lowerlim + iwinsz
        lowerlim = lowerlim - iwinsz
    end subroutine sqwin_1d_1

    !> one-dimensional symmetric hard window with limits
    pure subroutine sqwin_1d_2( x, winsz, lims, win )
        real,    intent(in)  :: x       !< input point
        real,    intent(in)  :: winsz   !< window size
        integer, intent(in)  :: lims(2) !< bounds
        integer, intent(out) :: win(2)  !< window
        integer :: iwinsz
        win(:) = nint(x)
        iwinsz = ceiling(winsz - 0.5)
        win(1) = max(lims(1), win(1) - iwinsz)
        win(2) = min(lims(2), win(2) + iwinsz)
    end subroutine sqwin_1d_2

    !> two-dimensional symmetric hard window
    pure subroutine sqwin_2d_1( x, y, winsz, win )
        real,    intent(in)  :: x,y      !< input point
        real,    intent(in)  :: winsz    !< window size
        integer, intent(out) :: win(2,2) !< window
        call sqwin_1d_1(x, winsz, win(1,1), win(1,2))
        call sqwin_1d_1(y, winsz, win(2,1), win(2,2))
    end subroutine sqwin_2d_1

    !> two-dimensional symmetric hard window with limits
    pure subroutine sqwin_2d_2( x, y, winsz, lims, win )
        real,    intent(in)  :: x,y       !< input point
        real,    intent(in)  :: winsz     !< window size
        integer, intent(in)  :: lims(2,2) !< bounds
        integer, intent(out) :: win(2,2)  !< window
        call sqwin_1d_2(x,winsz,lims(1,:), win(1,:))
        call sqwin_1d_2(y,winsz,lims(2,:), win(2,:))
    end subroutine sqwin_2d_2

    !> three-dimensional symmetric hard window
    pure subroutine sqwin_3d_1( x, y, z, winsz, win )
        real,    intent(in)  :: x,y,z    !< input point
        real,    intent(in)  :: winsz    !< window size
        integer, intent(out) :: win(3,2) !< window
        call sqwin_1d_1(x, winsz, win(1,1), win(1,2))
        call sqwin_1d_1(y, winsz, win(2,1), win(2,2))
        call sqwin_1d_1(z, winsz, win(3,1), win(3,2))
    end subroutine sqwin_3d_1

    !> three-dimensional symmetric hard window with limits
    pure subroutine sqwin_3d_2( x, y, z, winsz, lims, win )
        real,    intent(in)  :: x,y,z     !< input point
        real,    intent(in)  :: winsz     !< window size
        integer, intent(in)  :: lims(3,2) !< bounds
        integer, intent(out) :: win(3,2)  !< window
        call sqwin_1d_2(x,winsz,[lims(1,1), lims(1,2)], win(1,:))
        call sqwin_1d_2(y,winsz,[lims(2,1), lims(2,2)], win(2,:))
        call sqwin_1d_2(z,winsz,[lims(3,1), lims(3,2)], win(3,:))
    end subroutine sqwin_3d_2

end module simple_edges_sqwins