program simple_test_binimage

    include 'simple_lib.f08'
    use simple_binimage

    type(binimage)       :: test_binimage, test_binimage_3D, ccimage
    integer, allocatable :: imat_ccs(:,:,:)
    integer :: imat(4,4,1), imat_3D(3,3,2), nccs, i
    real    :: dist, dist_3D

    !! 2D test
    print *, '2D TEST'
    call test_binimage%new_bimg(ldim=[4,4,1], smpd=1.)
    imat(1,:,1) = [1,0,0,0]
    imat(2,:,1) = [1,0,1,0]
    imat(3,:,1) = [0,0,0,0]
    imat(4,:,1) = [0,0,0,1]
    call test_binimage%set_imat(imat)
    call test_binimage%max_dist(dist)
    print *, 'MAXIMUM DISTANCE FROM CENTER IS ', dist
    call test_binimage%find_ccs(ccimage)
    call ccimage%get_nccs(nccs)
    print *, nccs
    call ccimage%get_imat(imat_ccs)
    do i = 1,4
        print *, imat_ccs(i,:,1)
    enddo
    call test_binimage%find2_ccs(ccimage)
    print *, nccs
    call ccimage%get_imat(imat_ccs)
    do i = 1,4
        print *, imat_ccs(i,:,1)
    enddo

    print *, ' '

    !! 3D test
    print *, '3D TEST'
    call test_binimage_3D%new_bimg(ldim=[3,3,2],smpd=1.)
    imat_3D = reshape( (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1 /), shape(imat_3D))
    call test_binimage_3D%set_imat(imat_3D)
    call test_binimage_3D%max_dist(dist_3D)
    print *, 'MAXIMUM DISTANCE FROM CENTER IS ', dist_3D


end program simple_test_binimage