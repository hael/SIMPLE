program simple_test_calc_gradient
    include 'simple_lib.f08'
    use simple_image, only: image
    implicit none
    type(image) :: img
    integer     :: ldim(3)
    integer     :: k, l
    ldim(1) = 31
    ldim(2) = 31
    ldim(3) = 1
    call img%new(ldim, 1.)
    do k = -15, 15
    do l = -15, 15
        call img%set_rmat_at(k+16,l+16,1,(k/5.)**2 * (l/5.)**3)
    enddo
    enddo
    write(*, *) img%get_rmat_at(2, 3, 1)
end program simple_test_calc_gradient