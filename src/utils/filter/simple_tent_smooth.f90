! Module: simple_tent_smooth
!
! Applies a separable tent (triangle-shaped) smoothing filter to 3-D volumes.
!
! A tent filter is approximated by convolving two box (uniform averaging) filters
! of the same width along each axis.  Because the box filter is a sliding-window
! mean it can be computed in O(N) per row regardless of the filter radius.
!
! Boundary treatment: replicate padding — out-of-bounds positions are clamped to
! the nearest edge value.
module simple_tent_smooth
implicit none

contains

    ! Apply a separable tent-filter smooth to a 3-D real volume.
    !
    ! Arguments:
    !   vol              - volume to smooth (updated in-place on return)
    !   tmp              - scratch array with the same shape as vol
    !   nx, ny, nz       - dimensions of vol and tmp
    !   radius           - half-width of the tent kernel; the base of the tent
    !                      spans 2*radius + 1 samples
    subroutine tent_smooth_3d(vol, tmp, nx, ny, nz, radius)
        integer, intent(in)    :: nx, ny, nz, radius
        real,    intent(inout) :: vol(nx,ny,nz)
        real,    intent(inout) :: tmp(nx,ny,nz)
        integer :: box_width, lext1, rext1, lext2, rext2
        ! Width of each box-filter pass (two passes back-to-back approximate a tent)
        box_width = radius + 1
        ! Split box_width asymmetrically into left/right extents for the first pass,
        ! then mirror those extents for the second pass. The mirror ensures that the
        ! combined response is symmetric around each output sample.
        lext1 = (box_width - 1) / 2     ! floor half-width on the left
        rext1 = (box_width - 1) - lext1 ! ceiling half-width on the right
        lext2 = rext1                   ! swap left/right for the second pass
        rext2 = lext1
        ! Smooth along X: two successive box passes produce a tent-shaped kernel
        call box_filter_x(vol, tmp, nx, ny, nz, lext1, rext1)
        call box_filter_x(tmp, vol, nx, ny, nz, lext2, rext2)
        ! Smooth along Y
        call box_filter_y(vol, tmp, nx, ny, nz, lext1, rext1)
        call box_filter_y(tmp, vol, nx, ny, nz, lext2, rext2)
        ! Smooth along Z
        call box_filter_z(vol, tmp, nx, ny, nz, lext1, rext1)
        call box_filter_z(tmp, vol, nx, ny, nz, lext2, rext2)
    end subroutine tent_smooth_3d

    ! Apply a 1-D box filter to every row along the X axis of a 3-D array.
    ! X rows are contiguous in memory, so no temporary buffer is needed.
    subroutine box_filter_x(src, dst, nx, ny, nz, lext, rext)
        integer, intent(in)  :: nx, ny, nz, lext, rext
        real,    intent(in)  :: src(nx,ny,nz)
        real,    intent(out) :: dst(nx,ny,nz)
        integer :: m, l
        !$omp parallel do collapse(2) schedule(static) default(shared) private(m,l)
        do l = 1, nz
            do m = 1, ny
                call box_filter_1d(src(:,m,l), dst(:,m,l), nx, lext, rext)
            end do
        end do
        !$omp end parallel do
    end subroutine box_filter_x

    ! Apply a 1-D box filter to every column along the Y axis of a 3-D array.
    ! Y columns are non-contiguous, so each column is copied into a contiguous
    ! buffer before filtering to avoid strided memory access inside box_filter_1d.
    subroutine box_filter_y(src, dst, nx, ny, nz, lext, rext)
        integer, intent(in)  :: nx, ny, nz, lext, rext
        real,    intent(in)  :: src(nx,ny,nz)
        real,    intent(out) :: dst(nx,ny,nz)
        integer :: n, l
        real    :: buf_in(ny), buf_out(ny)
        !$omp parallel do collapse(2) schedule(static) default(shared) private(n,l,buf_in,buf_out)
        do l = 1, nz
            do n = 1, nx
                buf_in       = src(n, :, l)
                call box_filter_1d(buf_in, buf_out, ny, lext, rext)
                dst(n, :, l) = buf_out
            end do
        end do
        !$omp end parallel do
    end subroutine box_filter_y

    ! Apply a 1-D box filter to every pillar along the Z axis of a 3-D array.
    ! Z pillars are non-contiguous, so a contiguous buffer is used (see box_filter_y).
    subroutine box_filter_z(src, dst, nx, ny, nz, lext, rext)
        integer, intent(in)  :: nx, ny, nz, lext, rext
        real,    intent(in)  :: src(nx,ny,nz)
        real,    intent(out) :: dst(nx,ny,nz)
        integer :: n, m
        real    :: buf_in(nz), buf_out(nz)
        !$omp parallel do collapse(2) schedule(static) default(shared) private(n,m,buf_in,buf_out)
        do m = 1, ny
            do n = 1, nx
                buf_in       = src(n, m, :)
                call box_filter_1d(buf_in, buf_out, nz, lext, rext)
                dst(n, m, :) = buf_out
            end do
        end do
        !$omp end parallel do
    end subroutine box_filter_z

    ! Apply a 1-D box (sliding-window mean) filter to a 1-D array.
    !
    ! The kernel window around output position i spans src(i - lext : i + rext),
    ! giving a box of width (lext + rext + 1).  Out-of-bounds indices are handled
    ! with replicate (clamp-to-edge) padding.
    !
    ! Complexity: O(N) — a running sum is updated by adding the new leading sample
    ! and subtracting the outgoing trailing sample, so cost is independent of the
    ! kernel width.
    subroutine box_filter_1d(src, dst, n, lext, rext)
        integer, intent(in)  :: n, lext, rext
        real,    intent(in)  :: src(n)
        real,    intent(out) :: dst(n)
        integer :: i, j, box_width
        real    :: running_sum
        box_width = lext + rext + 1
        ! Bootstrap: accumulate the full box sum centred on i = 1
        running_sum = 0.0
        do j = -lext, rext
            running_sum = running_sum + src(clamp_index(1 + j, 1, n))
        end do
        dst(1) = running_sum / real(box_width)
        ! Slide the window one step at a time: add the incoming right-edge sample
        ! and subtract the outgoing left-edge sample
        do i = 2, n
            running_sum = running_sum &
                        + src(clamp_index(i + rext,     1, n)) &
                        - src(clamp_index(i - 1 - lext, 1, n))
            dst(i) = running_sum / real(box_width)
        end do
    end subroutine box_filter_1d

    ! Clamp integer i to the closed interval [lo, hi].
    pure integer function clamp_index(i, lo, hi)
        integer, intent(in) :: i, lo, hi
        clamp_index = max(lo, min(hi, i))
    end function clamp_index

end module simple_tent_smooth
