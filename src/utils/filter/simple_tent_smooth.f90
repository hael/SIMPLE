!@descr: separable tent (triangle) smoothing of 3D volumes by two passes of a sliding box mean
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
    ! X rows are contiguous in memory, so no temporary buffer is needed. A row of zeros (the outside of a mask) stays
    ! zero, so it is written without filtering.
    subroutine box_filter_x(src, dst, nx, ny, nz, lext, rext)
        integer, intent(in)  :: nx, ny, nz, lext, rext
        real,    intent(in)  :: src(nx,ny,nz)
        real,    intent(out) :: dst(nx,ny,nz)
        integer :: m, l
        !$omp parallel do collapse(2) schedule(static) default(shared) private(m,l)
        do l = 1, nz
            do m = 1, ny
                if( is_zero(src(:,m,l), nx) )then
                    dst(:,m,l) = 0.
                else
                    call box_filter_1d(src(:,m,l), dst(:,m,l), nx, lext, rext)
                endif
            end do
        end do
        !$omp end parallel do
    end subroutine box_filter_x

    ! Apply a 1-D box filter to every column along the Y axis of a 3-D array.
    ! The columns of one XY plane are filtered together, a running sum per column, so the loop reads whole X rows;
    ! each column keeps the operation sequence of box_filter_1d. A plane of zeros stays zero.
    subroutine box_filter_y(src, dst, nx, ny, nz, lext, rext)
        integer, intent(in)  :: nx, ny, nz, lext, rext
        real,    intent(in)  :: src(nx,ny,nz)
        real,    intent(out) :: dst(nx,ny,nz)
        integer :: i, j, l, box_width
        real    :: running_sum(nx)
        box_width = lext + rext + 1
        !$omp parallel do schedule(static) default(shared) private(i,j,l,running_sum)
        do l = 1, nz
            if( is_zero(src(:,:,l), nx*ny) )then
                dst(:,:,l) = 0.
                cycle
            endif
            running_sum = 0.0
            do j = -lext, rext
                running_sum = running_sum + src(:,clamp_index(1 + j, 1, ny),l)
            end do
            dst(:,1,l) = running_sum / real(box_width)
            do i = 2, ny
                running_sum = running_sum &
                            + src(:,clamp_index(i + rext,     1, ny),l) &
                            - src(:,clamp_index(i - 1 - lext, 1, ny),l)
                dst(:,i,l) = running_sum / real(box_width)
            end do
        end do
        !$omp end parallel do
    end subroutine box_filter_y

    ! Apply a 1-D box filter to every pillar along the Z axis of a 3-D array.
    ! The pillars of one XZ slab are filtered together, as in box_filter_y. A slab of zeros stays zero.
    subroutine box_filter_z(src, dst, nx, ny, nz, lext, rext)
        integer, intent(in)  :: nx, ny, nz, lext, rext
        real,    intent(in)  :: src(nx,ny,nz)
        real,    intent(out) :: dst(nx,ny,nz)
        integer :: i, j, m, box_width
        real    :: running_sum(nx)
        box_width = lext + rext + 1
        !$omp parallel do schedule(static) default(shared) private(i,j,m,running_sum)
        do m = 1, ny
            if( is_zero_slab(src, nx, ny, nz, m) )then
                dst(:,m,:) = 0.
                cycle
            endif
            running_sum = 0.0
            do j = -lext, rext
                running_sum = running_sum + src(:,m,clamp_index(1 + j, 1, nz))
            end do
            dst(:,m,1) = running_sum / real(box_width)
            do i = 2, nz
                running_sum = running_sum &
                            + src(:,m,clamp_index(i + rext,     1, nz)) &
                            - src(:,m,clamp_index(i - 1 - lext, 1, nz))
                dst(:,m,i) = running_sum / real(box_width)
            end do
        end do
        !$omp end parallel do
    end subroutine box_filter_z

    ! all n values are zero (either sign)
    pure logical function is_zero(v, n)
        integer, intent(in) :: n
        real,    intent(in) :: v(n)
        integer :: i
        is_zero = .false.
        do i = 1, n
            if( v(i) /= 0. ) return
        end do
        is_zero = .true.
    end function is_zero

    ! the XZ slab y = m of a volume is zero
    pure logical function is_zero_slab(v, nx, ny, nz, m)
        integer, intent(in) :: nx, ny, nz, m
        real,    intent(in) :: v(nx,ny,nz)
        integer :: k
        is_zero_slab = .false.
        do k = 1, nz
            if( .not. is_zero(v(:,m,k), nx) ) return
        end do
        is_zero_slab = .true.
    end function is_zero_slab

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
