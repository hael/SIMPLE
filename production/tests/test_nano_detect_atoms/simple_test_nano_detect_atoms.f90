program simple_test_nano_detect_atoms

    include 'simple_lib.f08'
    use simple_image
    use simple_nanoparticle
    use simple_nanoparticle_utils
    use simple_parameters
    use simple_defs_atoms
    use simple_atoms
    use simple_parameters
    use simple_math
    use simple_linalg
    use simple_strings

    character(len=100) :: filename, pdb_filename
    integer :: iostat, Z, ldim(3), boxsize, ldim_box(3)
    integer :: offset, nxyz(3), nxyz_offset(3), nboxes, xind, yind, zind, pos(3)
    integer :: xoff, yoff, zoff, npeaks, nbox, ibox, jbox, loc, ipeak, iimg
    real :: radius, smpd, sxx, thres, dist, dthres, length_diffs
    type(image) :: simulated_atom, nano_img, simulated_NP, boximg, test_NP
    type(image), allocatable :: atms_array(:)
    type(atoms) :: atom
    type(nanoparticle) :: nano
    type(parameters), target :: params
    integer, allocatable :: positions(:,:), inds_offset(:,:,:), pos_inds(:)
    real, allocatable :: box_scores(:,:,:), tmp(:), scores_cen(:,:,:), pos_scores(:), coords(:,:)
    real, allocatable :: rmat_sub(:,:,:), diffs(:)
    logical :: l_err_atom, l_err_box, is_peak
    logical, allocatable :: mask(:), selected_pos(:)

    real, parameter :: box_expansion_factor = 0.111

    filename = 'recvol_state01_iter005.mrc' ! first draft of 3D reconstruction 
    pdb_filename = 'ATMS.pdb' ! first draft of pdb file of atomic positions based on identify_atomic_positions

    ! checking that files can be opened without issue
    open(unit=35, file=trim(filename), iostat=iostat, status='old')
    if (iostat .eq. 0) then
        print *, 'File ' // trim(filename) // ' exists'
    else
        print *, 'An error occured when trying to open the file ' // trim(filename)
    end if
    close(35)

    open(unit=36, file=trim(pdb_filename), iostat=iostat, status='old')
    if (iostat .eq. 0) then
        print *, 'File ' // trim(pdb_filename) // ' exists'
    else
        print *, 'An error occured when trying to open the file ' // trim(pdb_filename)
    end if
    close(36)

    ! get atomic radius for platinum
    call get_element_Z_and_radius('PT', Z, radius)

    ! simulate one atom
    call atom%new(1)
    call atom%set_element(1,'PT')
    ldim = [120,120,120]
    smpd = 0.358
    boxsize = round2even(radius  / smpd) * 2
    ldim_box = [boxsize,boxsize,boxsize]
    call atom%set_coord(1,(smpd*real(ldim_box)/2.)) ! make sure atom is in center of box
    call simulated_atom%new(ldim_box,smpd)
    call atom%convolve(simulated_atom, cutoff=8*smpd)
    call simulated_atom%write('simulated_atom.mrc')
    call simulated_atom%prenorm4real_corr(l_err_atom)

    ! create new nanoparticle from image stored at filename
    params_glob => params
    params_glob%element = 'PT'
    params_glob%smpd = smpd
    call nano%new(filename)
    !call nano%get_img(nano_img)

    ! simulating noise free NP map (for testing)
    call nano%set_atomic_coords(trim(pdb_filename))
    call nano%simulate_atoms(simatms=simulated_NP)
    call simulated_NP%write('simulated_NP.mrc')

    ! set up picking infrastructure
    offset = 2
    nxyz = simulated_NP%get_ldim() - boxsize ! number of pixels in each direction in simulated image, with boxsize subtracted
    nxyz_offset = 0
    ! find number of boxes
    do xind = 0, nxyz(1), offset
        nxyz_offset(1) = nxyz_offset(1) + 1
        nxyz_offset(2) = 0
        do yind = 0, nxyz(2), offset
            nxyz_offset(2) = nxyz_offset(2) + 1
            nxyz_offset(3) = 0
            do zind = 0, nxyz(3), offset
                nxyz_offset(3) = nxyz_offset(3) + 1
                nboxes = nboxes + 1
            end do
        end do
    end do
    ! set up positions and inds_offset
    allocate(positions(nboxes,3), inds_offset(nxyz_offset(1),nxyz_offset(2),nxyz_offset(3)), source=0)
    nxyz_offset = 0
    nboxes = 0
    do xind = 0, nxyz(1), offset
        nxyz_offset(1) = nxyz_offset(1) + 1
        nxyz_offset(2) = 0
        do yind = 0, nxyz(2), offset
            nxyz_offset(2) = nxyz_offset(2) + 1
            nxyz_offset(3) = 0
            do zind = 0, nxyz(3), offset
                nxyz_offset(3) = nxyz_offset(3) + 1
                nboxes = nboxes + 1
                positions(nboxes,:) = [xind,yind,zind]
                inds_offset(nxyz_offset(1),nxyz_offset(2),nxyz_offset(3)) = nboxes
            end do
        end do
    end do
    allocate(box_scores(nxyz_offset(1),nxyz_offset(2),nxyz_offset(3)), source = -1.)

    ! iterate through positions in simulated image, compare to simulated atom 
    do xoff = 1, nxyz_offset(1)
        do yoff = 1, nxyz_offset(2)
            do zoff = 1, nxyz_offset(3)
                call boximg%new([boxsize,boxsize,boxsize],smpd)
                pos = positions(inds_offset(xoff,yoff,zoff),:)
                call window_slim_3D(simulated_NP, pos, boxsize, boximg)
                call boximg%prenorm4real_corr(l_err_box)
                box_scores(xoff,yoff,zoff) = simulated_atom%real_corr_prenorm(boximg)
                call boximg%kill
            end do 
        end do 
    end do

    ! find peak thresholding value
    tmp = pack(box_scores, mask=(box_scores > -1 + 1e-10))
    call detect_peak_thres(size(tmp), nboxes, tmp, thres)
    print *, 'Peak threshold is ', thres

    ! center filter
    allocate(scores_cen(nxyz_offset(1), nxyz_offset(2), nxyz_offset(3)))
    do xoff = 1, nxyz_offset(1)
        do yoff = 1, nxyz_offset(2)
            do zoff = 1, nxyz_offset(3)
            if( box_scores(xoff,yoff,zoff) >= thres )then
                call boximg%new([boxsize,boxsize,boxsize],smpd)
                pos  = positions(inds_offset(xoff,yoff,zoff),:)
                call window_slim_3D(simulated_NP, pos, boxsize, boximg)
                scores_cen(xoff, yoff, zoff) = boximg%box_cen_arg(boximg)
                call boximg%kill
            else
                scores_cen(xoff,yoff,zoff) = real(offset) + 1.
            endif
            end do
        end do
    end do
    print *, 'NPEAKS BEFORE CENTER FILTER = ', count(box_scores >= thres)
    npeaks = count(scores_cen <= real(offset))
    where( scores_cen <= real(offset) )
            ! there's a peak
    elsewhere
            box_scores = -1.
    endwhere
    print *, 'NPEAKS AFTER CENTER FILTER = ', npeaks
    
    ! distance filter
    dthres = real(offset)
    pos_inds   = pack(inds_offset(:,:,:),  mask=box_scores(:,:,:) >= thres)
    pos_scores = pack( box_scores(:,:,:),  mask=box_scores(:,:,:) >= thres)
    nbox       = size(pos_inds)
    allocate(mask(nbox),         source=.false.)
    allocate(selected_pos(nbox), source=.true. )
    do ibox = 1, nbox
        mask = .false.
        ! identify boxes in neighborhood
        do jbox = 1, nbox
            dist = euclid(real(positions(pos_inds(ibox),:)),real(positions(pos_inds(jbox),:)))
            if( dist <= dthres ) mask(jbox) = .true.
        end do
        ! find highest correlation score in neighborhood
        loc = maxloc(pos_scores, mask=mask, dim=1)
        ! eliminate all but the best
        mask(loc) = .false.
        where( mask ) selected_pos = .false.
    end do
    npeaks = count(selected_pos)
    print *, 'NPEAKS BEFORE DISTANCE FILTER = ', nbox
    print *, 'NPEAKS AFTER DISTANCE FILTER = ', npeaks
    ! update packed arrays
    pos_inds   = pack(pos_inds,   mask=selected_pos)
    pos_scores = pack(pos_scores, mask=selected_pos)
    ! update box scores
    do xoff = 1, nxyz_offset(1)
        do yoff = 1, nxyz_offset(2)
            do zoff = 1, nxyz_offset(3)
                is_peak = .false.
                do ipeak = 1,npeaks
                    if( pos_inds(ipeak) == inds_offset(xoff,yoff,zoff) )then
                        is_peak = .true.
                        exit
                    endif
                end do
                if( .not. is_peak ) box_scores(xoff,yoff,zoff) = -1.
            end do
        end do
    end do

    ! remove outliers
    
    ! make array of images containing the images of identified atoms and extract coordinates of peaks
    pos_inds   = pack(inds_offset(:,:,:),  mask=box_scores(:,:,:) >= thres)
    nbox = size(pos_inds, dim=1)
    allocate(coords(nbox,3))
    allocate(atms_array(nbox))
    do iimg = 1, nbox
        pos = positions(pos_inds(iimg),:)
        call atms_array(iimg)%new(ldim_box,smpd)
        call window_slim_3D(simulated_NP, pos, boxsize, atms_array(iimg))
        call atms_array(iimg)%write('boximgs/boximg_'//trim(int2str(iimg))//'.mrc')
        ! want coordinates of atoms to be at the center of the images
        call atms_array(iimg)%masscen(coords(iimg,:)) 
        coords(iimg,:) = coords(iimg,:) + real(atms_array(iimg)%get_ldim())/2. + pos !adjust center by size and position of box
    end do
    call write_centers('test_atomic_centers_masscen',coords)
    deallocate(coords)

    ! simulate nanoparticle using new coordinates
    ! and compare to original (simulated) NP?
    call nano%set_atomic_coords('test_atomic_centers_masscen.pdb')
    call nano%simulate_atoms(simatms=test_NP)
    call test_NP%write('test_NP.mrc')

    ! subtract extracted atom images from centered reference image
    ! one example image
    allocate(rmat_sub(boxsize,boxsize,boxsize), diffs(boxsize**3))
    rmat_sub = simulated_atom%get_rmat() - atms_array(1)%get_rmat()
    diffs = pack(rmat_sub,mask=.true.)
    !print *, diffs
    length_diffs = norm_2_sp(diffs)
    print *, 'DIFFERENCE BETWEEN ACTUAL AND REFERENCE IS ', length_diffs


    deallocate(diffs)
    deallocate(rmat_sub)
    deallocate(atms_array)
    deallocate(selected_pos)
    deallocate(mask)
    deallocate(scores_cen)
    deallocate(box_scores)
    deallocate(positions)
    deallocate(inds_offset)


    contains

    subroutine window_slim_3D( img_in, coord, box, img_out)
        class(image), intent(in)    :: img_in
        integer,      intent(in)    :: coord(3), box
        class(image), intent(inout) :: img_out
        integer :: fromc(3), toc(3)
        real, allocatable :: img_in_rmat(:,:,:)
        allocate(img_in_rmat,source=img_in%get_rmat())
        fromc = coord + 1
        toc   = fromc + box - 1 
        call img_out%set_rmat(img_in_rmat(fromc(1):toc(1), fromc(2):toc(2), fromc(3):toc(3)), ft=.false.)
        deallocate(img_in_rmat) 
    end subroutine window_slim_3D

    function avg_loc_sdev_3D( img_in, winsz ) result( asdev )
        class(image), intent(in) :: img_in
        integer,      intent(in) :: winsz
        real    :: avg, asdev
        real, allocatable :: sdevs(:,:,:), rmat(:,:,:)
        integer :: i, j, k, ir(2), jr(2), kr(2), isz, jsz, ksz, npix, ldim(3)
        ldim = img_in%get_ldim()
        allocate(sdevs(ldim(1),ldim(2),ldim(3)))
        allocate(rmat( ldim(1),ldim(2),ldim(3)))
        rmat = img_in%get_rmat()
        do i = 1,ldim(1)
            ir(1) = max(1,       i - winsz)
            ir(2) = min(ldim(1), i + winsz)
            isz   = ir(2) - ir(1) + 1
            do j = 1,ldim(2)
                jr(1) = max(1,       j - winsz)
                jr(2) = min(ldim(2), j + winsz)
                jsz   = jr(2) - jr(1) + 1
                do k = 1, ldim(3)
                    kr(1) = max(1,       k - winsz)
                    kr(2) = min(ldim(3), k + winsz)
                    ksz   = kr(2) - kr(1) + 1
                    npix         = isz * jsz + ksz
                    avg          = sum(rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2))) / real(npix)
                    sdevs(i,j,k) = sqrt(sum((rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2)) - avg)**2.0) / real(npix - 1))
                end do
            end do
        end do
        asdev = sum(sdevs) / real(ldim(1) * ldim(2) * ldim(3))
        deallocate(rmat)
        deallocate(sdevs)
    end function avg_loc_sdev_3D

    subroutine write_centers(fname, coords )
        character(len=*),           intent(in)    :: fname
        real,                       intent(in)    :: coords(:,:)
        type(atoms) :: centers_pdb
        integer     :: cc
        call centers_pdb%new(size(coords, dim = 1), dummy=.true.)
        do cc = 1, size(coords, dim = 1)
            call centers_pdb%set_name(cc,params_glob%element)
            call centers_pdb%set_element(cc,params_glob%element)
            call centers_pdb%set_coord(cc,coords(cc,:)*smpd)
            !call centers_pdb%set_beta(cc,nano%atominfo(cc)%valid_corr) ! use per atom valid corr
            call centers_pdb%set_resnum(cc,cc)
        enddo
        call centers_pdb%writepdb(fname)

    end subroutine write_centers


end program simple_test_nano_detect_atoms