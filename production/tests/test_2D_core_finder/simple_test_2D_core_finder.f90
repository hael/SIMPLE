program simple_test_2D_core_finder

    include 'simple_lib.f08'
    use simple_image,        only: image
    implicit none

    type(image),    allocatable :: cavgs(:)
    real        :: smpd=0.358, shell_size, r, max_rad, mean, var
    real, allocatable      :: rmat(:,:,:)
    integer     :: i, j, iref, n, nshells, nrefs, ldim(3), ldim_refs(3), ifoo
    logical, allocatable    :: mask(:,:,:)
    character(*), parameter    :: stks='cavgs.mrc'

    ! Read in stacks
    call find_ldim_nptcls(stks, ldim, ifoo)
    ldim_refs = [ldim(1), ldim(2), 1]
    nrefs = ldim(3)
    allocate(cavgs(nrefs))
    do iref = 1,nrefs
        call cavgs(iref)%new(ldim_refs, SMPD)
        call cavgs(iref)%read(stks, iref)
    end do

    shell_size = 4 / 0.358 ! Size of shell in voxel lengths
    max_rad = 15 / 0.358 ! Maximum NP radius in voxel lengths
    nshells = max_rad / shell_size
    allocate(mask(ldim_refs(1), ldim_refs(2), ldim_refs(3)))

    do iref=1,nrefs
        print *, "iref", iref
        rmat = cavgs(iref)%get_rmat()
        do n=0, nshells
            mask = .false.
            do i=1, ldim_refs(1)
                do j=1, ldim_refs(2)
                    r = sqrt(real(i**2 + j**2))
                    if (r > n*shell_size .and. r < (n+1)*shell_size) then
                        mask(i,j,1) = .true.
                    end if
                end do
            end do
            ! Calculate variance within mask
            mean = sum(rmat, mask=mask) / count(mask)
            rmat = rmat - mean
            var = sum(rmat**2) / (count(mask)-1)
            ! Calculate variance within mask
            print *, count(mask), var, mean, n*shell_size*smpd, (n+1)*shell_size*smpd
        end do
    end do


end program simple_test_2D_core_finder