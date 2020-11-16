! Just a little program to demonstrate how phase-correlations work
! and read off the convention using the image class (direction and offset of phase-correlations)
! Phase-correlations: FT^-1{ FT{A} * FT{B} }

module mod_phasecorr
    use simple_image

    real,    parameter :: smpd = 1.0
    integer, parameter :: Nx   = 1000      ! size of image
    integer, parameter :: Ny   = 1000      ! size of image
    real, parameter :: center_A_x = 100.   ! center of half-circle A
    real, parameter :: center_A_y = 100.   ! center of half-circle A
    real, parameter :: center_B_x = 200.   ! center of half-circle B
    real, parameter :: center_B_y = 200.   ! center of half-circle B
    real, parameter :: radius     =  80.   ! radius of half-circles

    type :: t_phasecorr
        type(image) :: imgA, imgB          ! images to determine phase-correlations from
        type(image) :: imgC                ! used for calculation
    contains
        procedure :: new
        procedure :: run
        procedure :: kill

        procedure :: alloc
        procedure :: gen_images            ! generate images (two half-circles)
        procedure :: vis_imgs
    end type t_phasecorr

contains

    subroutine new( self )
        class(t_phasecorr), intent(inout) :: self
        call self%alloc
        call self%gen_images
    end subroutine new

    subroutine run( self )
        class(t_phasecorr), intent(inout) :: self
        complex, pointer :: cmat_A(:,:,:), cmat_B(:,:,:), cmat_C(:,:,:)
        real,    pointer :: rmat_C(:,:,:)
        integer :: ix, iy, ix_max, iy_max
        real    :: max_val
        call self%vis_imgs
        call self%imgA%fft
        call self%imgB%fft
        call self%imgC%fft
        call self%imgA%get_cmat_ptr(cmat_A)
        call self%imgB%get_cmat_ptr(cmat_B)
        call self%imgC%get_cmat_ptr(cmat_C)
        cmat_C = cmat_A * conjg(cmat_B)
        call self%imgC%ifft
        call self%imgC%vis
        call self%imgC%get_rmat_ptr(rmat_C)
        max_val = -HUGE(max_val)
        do ix = 1, Nx
            do iy = 1, Ny
                if (rmat_C(ix,iy,1) > max_val) then
                    ix_max = ix
                    iy_max = iy
                    max_val = rmat_C(ix,iy,1)
                end if
            end do
        end do
        write (*,*) 'max peak location: ', ix_max, iy_max
    end subroutine run

    subroutine kill( self )
        class(t_phasecorr), intent(inout) :: self
        call self%imgA%kill
        call self%imgB%kill
        call self%imgC%kill
    end subroutine kill

    subroutine alloc( self )
        class(t_phasecorr), intent(inout) :: self
        integer :: ldim(3)
        ldim(1) = Nx
        ldim(2) = Ny
        ldim(3) = 1
        call self%imgA%new(ldim, smpd)
        call self%imgB%new(ldim, smpd)
        call self%imgC%new(ldim, smpd)
    end subroutine alloc

    subroutine gen_images( self )
        class(t_phasecorr), intent(inout) :: self
        integer :: ix, iy
        real    :: distx, disty, dist
        real    :: tmp
        real, pointer :: rmat_A(:,:,:), rmat_B(:,:,:)
        call self%imgA%get_rmat_ptr(rmat_A)
        call self%imgB%get_rmat_ptr(rmat_B)
        rmat_A = 0.
        rmat_B = 0.
        do ix = 1, Nx
            do iy = 1, Ny
                distx = real(ix) - center_A_x
                disty = real(iy) - center_A_y
                dist  = sqrt(distx**2 + disty**2)
                if ((distx < 0.) .or. (distx > radius) .or.&
                                      (disty > radius)) then
                    tmp = 0.
                else
                    tmp = max(radius - dist, 0.)
                end if
                rmat_A(ix, iy, 1) = tmp
                distx = real(ix) - center_B_x
                disty = real(iy) - center_B_y
                dist  = sqrt(distx**2 + disty**2)
                if ((distx < 0.) .or. (distx > radius) .or. &
                                      (disty > radius)) then
                    tmp = 0.
                else
                    tmp = max(radius - dist, 0.)
                end if
                rmat_B(ix, iy, 1) = tmp
            end do
        end do
    end subroutine gen_images

    subroutine vis_imgs( self )
        class(t_phasecorr), intent(inout) :: self
        call self%imgA%vis
        call self%imgB%vis
    end subroutine vis_imgs
end module mod_phasecorr

program test_phasecorr
    use mod_phasecorr
    implicit none

    type(t_phasecorr) :: tphasecorr
    call tphasecorr%new
    call tphasecorr%run
    call tphasecorr%kill

end program test_phasecorr
