!@descr: stateless utilities over an oris collection: the orientation distribution of a state's particles
! oridist_from_oris bins the projection directions of one state's particles over azimuth -180..180
! and elevation -90..90, with as many bins along each axis as the histogram it fills. The GUI's
! volume metadata holds such a histogram (gui_metadata_vol3D%set_oridist, on a grid of
! ORIDIST_NBINS_X x ORIDIST_NBINS_Y); the bins come from the caller, so this module needs nothing
! of the GUI. (simple_ori_utils holds the Euler-angle utilities of a single orientation; it cannot
! import simple_oris.)
module simple_oris_utils
use simple_oris,   only: oris
use simple_linalg, only: rad2deg
implicit none

public :: oridist_from_oris
private

contains

    !> Counts the particles of @p state in @p os by projection direction into @p hist: azimuth
    !! -180..180 along its first dimension, elevation -90..90 along its second.
    subroutine oridist_from_oris( os, state, hist )
        class(oris), intent(in)  :: os
        integer,     intent(in)  :: state
        integer,     intent(out) :: hist(:,:)
        real    :: normal(3), azimuth, elevation, xwidth, ywidth
        integer :: nx, ny, iptcl, ix, iy
        nx     = size(hist, 1)
        ny     = size(hist, 2)
        xwidth = 360. / real(nx)
        ywidth = 180. / real(ny)
        hist   = 0
        do iptcl = 1,os%get_noris()
            if( os%get_state(iptcl) /= state ) cycle
            normal    = os%get_normal(iptcl)
            azimuth   = rad2deg(atan2(normal(2), normal(1)))
            elevation = rad2deg(asin(max(-1.0, min(1.0, normal(3)))))
            ix = min(nx, max(1, floor((azimuth   + 180.) / xwidth) + 1))
            iy = min(ny, max(1, floor((elevation +  90.) / ywidth) + 1))
            hist(ix,iy) = hist(ix,iy) + 1
        enddo
    end subroutine oridist_from_oris

end module simple_oris_utils
