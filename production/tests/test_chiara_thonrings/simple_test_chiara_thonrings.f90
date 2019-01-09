module simple_test_chiara_thonrings_mod
    use simple_powerspec_analysis
    use simple_image
    use simple_math
    use simple_segmentation
implicit none
contains

    !This function takes in input an image (suppsed to be
    !a connected components image) and deletes (sets to 0)
    !the internal shells (up tp 5A resolution) and the
    !external frame (up to 3. A).
    !It is meant to discard the connected components which
    !are not in the range to represent ice contaminations.
    subroutine delete_shells(img, box, smpd)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: box
        real,        intent(in)    :: smpd
        integer :: x, ldim(3), i, j
        real, allocatable :: rmat(:,:,:)
        x = calc_fourier_index(5.,box,smpd)
        ldim = img%get_ldim()
        rmat = img%get_rmat()
        !Discard internal shell
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if((real(i)-real(box/2))**2/(real(x)**2) + (real(j)-real(box/2))**2/(real(x)**2) - 1 < TINY) rmat(i,j,1) = 0.
            enddo
        enddo
        !Discard external frame
        x = calc_fourier_index(3.,box,smpd)
        rmat(:box/2-x,:,1) = 0.
        rmat(:,:box/2-x,1) = 0.
        rmat(box/2+x:,:,1) = 0.
        rmat(:,box/2+x:,1) = 0.
        call img%set_rmat(rmat)
        call img%write('Deleted.mrc')
    end subroutine delete_shells

    subroutine extract_3points(img, label, points)
        type(image), intent(inout) :: img
        integer,     intent(in)    :: label
        integer,     intent(out)   :: points(3,2)
        real, allocatable :: rmat(:,:,:), rmat_2d(:,:)
        integer :: i, ldim(3)
        rmat = img%get_rmat()
        where(abs(rmat-real(label)) > TINY) rmat = 0.
        ldim = img%get_ldim()
        allocate(rmat_2d(ldim(1), ldim(2)), source = 0.)
        rmat_2d(:,:) = rmat(:,:,1)
        do i = 1, 3
            points(i,:) = minloc(abs(rmat_2d(:,:)-real(label)))
            if(points(i,1)-1 > 0 .and. points(i,1)+1< ldim(1)) then
                rmat_2d(points(i,1)-1:points(i,1)+1,points(i,2)) = 0.
                rmat_2d(points(i,1)-1:points(i,1)+1,points(i,2)) = 0.
            else
                rmat_2d(points(i,1),points(i,2)) = 0.
            endif
            if(points(i,2)-1 > 0 .and. points(i,2)+1< ldim(2)) then
                rmat_2d(points(i,1),points(i,2)-1:points(i,2)+1) = 0.
                rmat_2d(points(i,1),points(i,2)-1:points(i,2)+1) = 0.
            else
                rmat_2d(points(i,1),points(i,2)) = 0.
            endif
        enddo
        if(any(points(:,:) .eq. 1))then
            rmat_2d(:,:) = rmat(:,:,1) !restore
            do i =1, 3
                points(i,:) = minloc(abs(rmat_2d(:,:)-real(label)))
                rmat_2d(points(i,1),points(i,2)) = 0.
            enddo
        endif
        deallocate(rmat, rmat_2d)
    end subroutine extract_3points

    subroutine identify_center_radius(points, center, radius)
        integer,     intent(in)    :: points(3,2)
        integer,     intent(out)   :: center(2)
        real,        intent(out)   :: radius
        real :: x, y
        if(   (2*(points(1,1)*(points(2,2)-points(3,2)) - &
            &     points(1,2)*(points(2,1)-points(3,1)) + &
            &     points(2,1)*points(3,2)-points(3,1)*points(2,2))) &
            &      .ne. 0) then
            center(1) = ((points(1,1)**2+points(1,2)**2)*(points(2,2)-points(3,2)) + &
            &            (points(2,1)**2+points(2,2)**2)*(points(3,2)-points(1,2)) + &
            &            (points(3,1)**2+points(3,2)**2)*(points(1,2)-points(2,2)))/ &
            & (2*(points(1,1)*(points(2,2)-points(3,2)) - &
            &     points(1,2)*(points(2,1)-points(3,1)) + &
            &     points(2,1)*points(3,2)-points(3,1)*points(2,2)))

            center(2) = ((points(1,1)**2+points(1,2)**2)*(points(3,1)-points(2,1)) + &
            &            (points(2,1)**2+points(2,2)**2)*(points(1,1)-points(3,1)) + &
            &            (points(3,1)**2+points(3,2)**2)*(points(2,1)-points(1,1)))/ &
            & (2*(points(1,1)*(points(2,2)-points(3,2)) - &
            &     points(1,2)*(points(2,1)-points(3,1)) + &
            &     points(2,1)*points(3,2)-points(3,1)*points(2,2)))
        else
            center(:) = points(1,:)
            print *, 'erroneous'
        endif

        x = real(center(1))
        y = real(center(2))
        radius = sqrt((x-real(points(1,1)))**2+(y-real(points(1,2)))**2)

        print *, 'CENTER = ', center, 'RADIUS = ', radius
    end subroutine identify_center_radius

    !This function  takes in input the box sz and the smpd
    !of a power spectra image and consequentially builds an
    !ice template (circular)
    function build_ice_template(box, smpd, radius) result(img_templ)
        integer, intent(in) :: box
        real,    intent(in) :: smpd
        real,    intent(in) :: radius
        type(image) :: img_templ
        call img_templ%new([box,box,1], smpd)
        call img_templ%ellipse([box/2,box/2], [radius,radius], 'yes')
    end function build_ice_template
end module simple_test_chiara_thonrings_mod

program simple_test_chiara_thonrings
    include 'simple_lib.f08'
    use simple_test_chiara_thonrings_mod
    use simple_powerspec_analysis
    use gnufor2
    use simple_ctf
    use simple_micops
    use simple_image
    use simple_math
    use simple_picker_chiara
    use simple_segmentation
    use simple_parameters, only: parameters
    use simple_cmdline,    only: cmdline
    type(image)       :: img, img_cc
    real, allocatable :: rmat(:,:,:)
    real, allocatable :: rmat_aux(:,:,:)
    integer :: i, ldim(3), box
    real :: label_mirror
    logical :: yes_no, discard
    real :: smpd
    integer :: points(3,2), px(2), j
    integer :: center(2)
    real :: radius
    logical :: outside
    type(image) :: img_win, img_templ, img_try
    real :: r
    integer :: h, k, m(3), sh, cnt, sz
    box = 512
    smpd = 1.
    cnt = 0
    call img%new([512,512,1],1.)
    ldim = img%get_ldim()
    call img_cc%new([512,512,1],1.)
    call img_try%new([512,512,1],1.)
    call img%read('/home/chiara/Desktop/Chiara/ThonRings/IceContaminants/SAGAWithICEBinCC.mrc')
    call img%find_connected_comps(img_cc)
    call delete_shells(img_cc, 512, 1.)
    call img_win%new([box/8, box/8, 1], smpd)
    rmat = img_cc%get_rmat()
    allocate(rmat_aux(ldim(1),ldim(2),1), source = 0.)
    do i = int(minval(rmat,rmat > 0.5)), int(maxval(rmat))
        sz = count(abs(rmat - real(i)) < TINY)
        if(sz < 5 ) cycle !discard too small cc
        print *,'Checking ', i, 'sz ', sz
        yes_no = is_symmetric(img_cc, i, discard, label_mirror)
        if(.not. discard) then !discard the cc that are not in the I or II quadrant
            if(yes_no) then    !if it is symmetric
                print *, 'CC number ', i, 'is symm with ', label_mirror
                call extract_3points(img_cc,i,points)
                print *, 'points = '
                call vis_mat(points)
                call identify_center_radius(points,center,radius)
                where(rmat-real(i) < TINY)
                     rmat_aux = 0.
                 elsewhere
                     rmat_aux = 1.
                 endwhere
                call img_try%set_rmat(rmat_aux)
                call img_try%window_slim(center-box/8/2, box/8, img_win, outside)
                print *, 'cc ', i, 'outside ', outside
                if(.not. outside) then
                    cnt = cnt + 1
                    call img_win%write('img_win.mrc', cnt)
                    img_templ = build_ice_template(box/8, smpd, radius)
                    call img_templ%write('img_templ.mrc', cnt)
                    r = img_win%real_corr(img_templ)
                    print *, 'correlation = ', r
                    rmat = img_cc%get_rmat()
                    m(:) =  minloc(abs(rmat- real(i)))
                    if(r > 0.4) print *, 'DETECTED ICE! At px ', m
                    h   = -box/2 + m(1) - 1
                    k   = -box/2 + m(2) - 1
                    sh  = nint(hyp(real(h),real(k)))
                    print *, 'shell = ', sh, ' Resolution = ', calc_fourier_index(real(sh),box,smpd)
                endif
            endif
        endif
    enddo
end program simple_test_chiara_thonrings
! call process_ps_stack('pspecs_saga_polii.mrc', 'analisedSAGA.mrc', 1.14, 35., 2, 10)
! call process_ps_stack('pspecs_sphire_tstdat.mrc', 'analisedSPHIRE.mrc', 1.41, 20.,1, 10)
