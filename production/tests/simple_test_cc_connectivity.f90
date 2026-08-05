program simple_test_cc_connectivity
! Pins the connectivity contract of image_bin%find_ccs.
!
! In 3D the flood fill must span the full 26-neighbourhood. It previously omitted
! the two axial neighbours (i,j,k-1) and (i,j,k+1), which split any structure one
! voxel wide in x and y that ran along z into one component per slice. In 2D the
! k-1 and k+1 planes fall outside the box, so the same code must still behave as
! plain 8-connectivity.
use simple_core_module_api
use simple_image_bin, only: image_bin
implicit none
#include "simple_local_flags.inc"
integer, parameter   :: N = 7
integer, parameter   :: C = 4          ! box centre, leaves a margin of 3 on each side
type(image_bin)      :: bimg, ccimg, bimg2D
integer, allocatable :: sz(:)
integer              :: imat(N,N,N), imat2D(N,N,1)
integer              :: nccs, di, dj, dk, nfail

call bimg%new_bimg([N,N,N], 1.)

! ---- every one of the 26 offsets must connect ----------------------------
nfail = 0
do dk = -1,1
    do dj = -1,1
        do di = -1,1
            if( di == 0 .and. dj == 0 .and. dk == 0 ) cycle
            imat = 0
            imat(C,C,C)             = 1
            imat(C+di,C+dj,C+dk)    = 1
            call bimg%set_imat(imat)
            call bimg%find_ccs(ccimg)
            call ccimg%get_nccs(nccs)
            if( nccs /= 1 )then
                write(logfhandle,'(A,3(1X,I2),A,I4)') 'offset', di, dj, dk, ' gave components: ', nccs
                nfail = nfail + 1
            endif
        end do
    end do
end do
if( nfail > 0 ) THROW_HARD('find_ccs is not 26-connected in 3D')

! ---- Chebyshev distance 2 must stay separate -----------------------------
! This is the sharp boundary: widening the neighbourhood must not start merging
! voxels that are genuinely apart.
nfail = 0
do dk = -2,2,2
    do dj = -2,2,2
        do di = -2,2,2
            if( di == 0 .and. dj == 0 .and. dk == 0 ) cycle
            imat = 0
            imat(C,C,C)          = 1
            imat(C+di,C+dj,C+dk) = 1
            call bimg%set_imat(imat)
            call bimg%find_ccs(ccimg)
            call ccimg%get_nccs(nccs)
            if( nccs /= 2 )then
                write(logfhandle,'(A,3(1X,I2),A,I4)') 'offset', di, dj, dk, ' gave components: ', nccs
                nfail = nfail + 1
            endif
        end do
    end do
end do
if( nfail > 0 ) THROW_HARD('find_ccs merged voxels that are more than one step apart')

! ---- the regression case: a 1x1xN column along z -------------------------
imat        = 0
imat(C,C,:) = 1
call bimg%set_imat(imat)
call bimg%find_ccs(ccimg)
call ccimg%get_nccs(nccs)
if( nccs /= 1 ) THROW_HARD('a 1x1xN column along z must be a single connected component')
sz = ccimg%size_ccs()
if( size(sz) /= 1 ) THROW_HARD('column produced more than one component size')
if( sz(1)    /= N ) THROW_HARD('column component has the wrong voxel count')

! ---- controls: the same column along x and along y -----------------------
imat        = 0
imat(:,C,C) = 1
call bimg%set_imat(imat)
call bimg%find_ccs(ccimg)
call ccimg%get_nccs(nccs)
if( nccs /= 1 ) THROW_HARD('a 1x1xN column along x must be a single connected component')
imat        = 0
imat(C,:,C) = 1
call bimg%set_imat(imat)
call bimg%find_ccs(ccimg)
call ccimg%get_nccs(nccs)
if( nccs /= 1 ) THROW_HARD('a 1x1xN column along y must be a single connected component')

! ---- degenerate volumes --------------------------------------------------
imat = 0
call bimg%set_imat(imat)
call bimg%find_ccs(ccimg)
call ccimg%get_nccs(nccs)
if( nccs /= 0 ) THROW_HARD('an empty volume must yield no connected components')
imat = 1
call bimg%set_imat(imat)
call bimg%find_ccs(ccimg)
call ccimg%get_nccs(nccs)
if( nccs /= 1 ) THROW_HARD('a fully occupied volume must yield exactly one component')
sz = ccimg%size_ccs()
if( sz(1) /= N*N*N ) THROW_HARD('full volume component has the wrong voxel count')

! ---- 2D must remain 8-connected -----------------------------------------
call bimg2D%new_bimg([N,N,1], 1.)
nfail = 0
do dj = -1,1
    do di = -1,1
        if( di == 0 .and. dj == 0 ) cycle
        imat2D = 0
        imat2D(C,C,1)         = 1
        imat2D(C+di,C+dj,1)   = 1
        call bimg2D%set_imat(imat2D)
        call bimg2D%find_ccs(ccimg)
        call ccimg%get_nccs(nccs)
        if( nccs /= 1 )then
            write(logfhandle,'(A,2(1X,I2),A,I4)') '2D offset', di, dj, ' gave components: ', nccs
            nfail = nfail + 1
        endif
    end do
end do
if( nfail > 0 ) THROW_HARD('find_ccs is not 8-connected in 2D')
imat2D          = 0
imat2D(C,C,1)   = 1
imat2D(C+2,C,1) = 1
call bimg2D%set_imat(imat2D)
call bimg2D%find_ccs(ccimg)
call ccimg%get_nccs(nccs)
if( nccs /= 2 ) THROW_HARD('find_ccs merged separated pixels in 2D')

call bimg%kill_bimg
call bimg2D%kill_bimg
call ccimg%kill_bimg
if( allocated(sz) ) deallocate(sz)
write(logfhandle,'(A)') 'Connected-component connectivity test passed'

end program simple_test_cc_connectivity
