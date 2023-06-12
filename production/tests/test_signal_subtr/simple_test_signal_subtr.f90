! Subtracts reprojections from class averages.  
! Will want to also compare ptcl reprojections to ptcls
program simple_test_signal_subtr

    include 'simple_lib.f08'
    use simple_image,        only: image
    implicit none
    
    type(image),    allocatable :: cavgs(:), reprojs(:), diff(:), shells(:)
    real, allocatable   :: rmat(:,:,:), rmat_shells(:,:,:)
    real        :: smpd=0.358, nprad_A=15, nprad_vox, shell_size, max_rad, r, mean, var, com(2), tot, bckgrnd, sdev_noise, center(3)
    integer     :: iref, nrefs, ldim1(3), ldim2(3), ldim_refs(3), ifoo, nshells, i, j, n
    logical, allocatable    :: mask(:,:,:)
    character(*), parameter    :: fn_cavgs='cavgs.mrc', fn_reprojs='reprojections/reprojs.mrc'
    character(*), parameter    :: fn_diff='cavgs_minus_reprojections.mrc', fn_shells='shells.mrc'
    character(*), parameter    :: cavgs_out='cavgs_out.mrc', reprojs_out='reprojs_out.mrc'

    ! Read in stacks
    call find_ldim_nptcls(fn_cavgs, ldim1, ifoo)
    call find_ldim_nptcls(fn_reprojs, ldim2, ifoo)
    if (ldim1(1) /= ldim2(1) .or. ldim1(2) /= ldim2(2)) then
        print *, "CAVGS AND REPROJS OF DIFFERENT DIMENSIONS"
        stop
    else if (ldim1(3) /= ldim2(3)) then
        print *, "CAVGS AND REPROJS OF DIFFERENT NUMBER"
        stop
    end if
    ldim_refs = [ldim1(1), ldim1(2), 1]
    nrefs = ldim1(3)
    allocate(cavgs(nrefs), reprojs(nrefs), diff(nrefs), shells(nrefs))
    do iref = 1,nrefs
        call cavgs(iref)%new(ldim_refs, SMPD)
        call reprojs(iref)%new(ldim_refs, SMPD)
        call shells(iref)%new(ldim_refs, smpd)
        call cavgs(iref)%read(fn_cavgs, iref)
        call reprojs(iref)%read(fn_reprojs, iref)
    end do

    ! Normalize for comparison
    allocate(mask(ldim_refs(1), ldim_refs(2), ldim_refs(3)), source=.true.)
    do iref=1, nrefs
         call cavgs(iref)%norm_bin()
         call reprojs(iref)%norm_bin()
     end do

    ! Subtract reprojs from cavgs
    ! We want the magnitude of the fractional difference since some regions 
    ! may have higher signal than other regions.
    do iref=1,nrefs
        call diff(iref)%new(ldim_refs, SMPD)
        diff(iref) = (cavgs(iref) - reprojs(iref)) / ((cavgs(iref) + reprojs(iref)))
        call diff(iref)%set_rmat(2*abs(diff(iref)%get_rmat()), ft=.false.)
    end do

    ! Take averages of shells out to the NP radius
    shell_size = 2 / 0.358 ! Size of shell in voxel lengths
    max_rad = 18 / 0.358 ! Maximum NP radius in voxel lengths
    nshells = max_rad / shell_size
    center = ldim_refs / 2 + 0.5
    allocate(rmat_shells(ldim_refs(1), ldim_refs(2), ldim_refs(3)), source=0.)
    do iref=1,nrefs
        print *, "iref", iref - 1
        rmat = diff(iref)%get_rmat()
        rmat_shells = 0.
        do n=0, nshells
            mask = .false.
            do i=1, ldim_refs(1)
                do j=1, ldim_refs(2)
                    r = sqrt(real((i - center(1))**2 + (j - center(2))**2))
                    if (r > n*shell_size .and. r < (n+1)*shell_size) then
                        mask(i,j,1) = .true.
                        rmat_shells(i,j,1) = n + 1
                    end if
                end do
            end do
            ! Calculate variance within mask
            mean = sum(rmat, mask=mask) / count(mask)
            print *, count(mask), abs(mean), maxval(rmat, mask=mask), n*shell_size*smpd, (n+1)*shell_size*smpd
            call shells(iref)%set_rmat(rmat_shells, ft=.false.)
        end do
    end do

    ! Write output
    do iref=1, nrefs
        call cavgs(iref)%write(cavgs_out, iref)
        call reprojs(iref)%write(reprojs_out, iref)
        call diff(iref)%write(fn_diff, iref)
        call shells(iref)%write(fn_shells, iref)
    end do

end program simple_test_signal_subtr

