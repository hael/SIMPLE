! This programs the surface-core NanoParticle (NP) structure from 2D NP class averages
! and 2D reprojections of the 3D volume. A stack containing the fractional
! differences between the class averages and reprojections is output.
! Average intensities over shells starting from the NP center out to the 
! NP surface are reported in a .csv file. If the surface-core model is correct,
! the fractional differences should be higher towards the surface, where
! the structure is more flexible and dynamic.
program simple_test_2D_core_finder

    include 'simple_lib.f08'
    use simple_image,        only: image
    implicit none
    
    type(image),    allocatable :: cavgs(:), reprojs(:), diff(:), shells(:)
    type(image)                 :: temp
    real, parameter     :: smpd=0.358, nprad_A=14.0, shell_size_A = 2.0
    real                :: nprad_vox=nprad_A/smpd, shell_size_vox=shell_size_A/smpd
    real                :: r, mean, center(2), minmax(2)
    real, allocatable   :: rmat(:,:,:), rmat_shells(:,:,:)
    integer             :: iref, nrefs, ldim1(3), ldim2(3), ldim_refs(3), ifoo, nshells, i, j, n, funit
    logical, allocatable       :: is_populated(:), mask(:,:,:)
    character(len=256)         :: fn_cavgs, fn_reprojs
    character(*), parameter    :: fn_diff='cavgs_minus_reprojections.mrc', fn_shells='shells.mrc'
    character(*), parameter    :: fn_cavgs_out='cavgs_out.mrc', fn_reprojs_out='reprojs_out.mrc'
    character(*), parameter    :: fn_results='results.csv'
    character(*), parameter    :: results_header='INDEX'//CSV_DELIM//'RMIN'//CSV_DELIM//'RMAX'//CSV_DELIM//&
                                  &'COUNT'//CSV_DELIM//'MEAN'//CSV_DELIM//'MAX_INT'

    ! Handle command-line input
    if( command_argument_count() /= 2 )then
        write(logfhandle,'(a)')  'Usage: simple_test_2D_core_finder cavgs.mrc reprojs.mrc'
        write(logfhandle,'(a)')  'cavgs.mrc: stack of 2D class averages'
        write(logfhandle,'(2a)') 'reprojs.mrc: stack of 2D reprojections of 3D volume along '//&
                                &'same orientations as cavgs.mrc', NEW_LINE('a')
        stop
    else
        call get_command_argument(1, fn_cavgs)
        call get_command_argument(2, fn_reprojs)
    end if

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

    ! Identify which class averages are populated so we can ignore unpopulated ones
    allocate(is_populated(ldim1(3)), source=.false.)
    call temp%new(ldim_refs, smpd)
    do i=1, ldim1(3)
        call temp%read(fn_cavgs, i)
        minmax = temp%minmax()
        if (minmax(1) /= 0. .or. minmax(2) /= 0.) then
            is_populated(i) = .true.
        end if
    end do
    call temp%kill()
    nrefs = count(is_populated)
    allocate(cavgs(nrefs), reprojs(nrefs), diff(nrefs), shells(nrefs))

    ! Read in the populated class averages and corresponding reprojections
    iref = 1
    do i=1,ldim1(3)
        if (is_populated(i)) then
            call cavgs(iref)%new(ldim_refs, smpd)
            call reprojs(iref)%new(ldim_refs, smpd)
            call shells(iref)%new(ldim_refs, smpd)
            call cavgs(iref)%read(fn_cavgs, i)
            call reprojs(iref)%read(fn_reprojs, i)
            iref = iref + 1
        end if
    end do

    ! Normalize for comparison
    ! Unsure if this is a proper way to normalize
    allocate(mask(ldim_refs(1), ldim_refs(2), ldim_refs(3)), source=.true.)
    do iref=1, nrefs
         call cavgs(iref)%norm_bin()
         call reprojs(iref)%norm_bin()
     end do

    ! Subtract reprojs from cavgs
    ! We want the magnitude of the fractional difference since some regions 
    ! may have higher signal than other regions.
    do iref=1,nrefs
        call diff(iref)%new(ldim_refs, smpd)
        diff(iref) = (cavgs(iref) - reprojs(iref)) / ((cavgs(iref) + reprojs(iref)))
        call diff(iref)%set_rmat(2*abs(diff(iref)%get_rmat()), ft=.false.)
    end do

    ! Take averages of shells out to the NP radius
    nshells = nprad_vox / shell_size_vox
    center = ldim_refs(1:2) / 2 + 0.5
    allocate(rmat_shells(ldim_refs(1), ldim_refs(2), ldim_refs(3)), source=0.)
    call fopen(funit, FILE=trim(fn_results), STATUS='REPLACE', action='WRITE')
    write(funit, '(A)') results_header
    601 format(F12.4,A2)
    602 format(F12.4)
    do iref=1,nrefs
        rmat = diff(iref)%get_rmat()
        rmat_shells = 0.
        do n=0, nshells
            mask = .false. ! For now use mask in case we want to calculate other stats in the future
            do i=1, ldim_refs(1)
                do j=1, ldim_refs(2)
                    r = sqrt(real((i - center(1))**2 + (j - center(2))**2))
                    if (r > n*shell_size_vox .and. r < (n+1)*shell_size_vox) then
                        mask(i,j,1) = .true.
                        rmat_shells(i,j,1) = n + 1
                    end if
                end do
            end do
            mean = sum(rmat, mask=mask) / count(mask)
            ! Write csv file containing statistics
            write(funit,601,advance='no') real(iref-1),                CSV_DELIM ! INDEX
            write(funit,601,advance='no') n*shell_size_A,              CSV_DELIM ! RMIN
            write(funit,601,advance='no') (n+1)*shell_size_A,          CSV_DELIM ! RMAX
            write(funit,601,advance='no') real(count(mask)),           CSV_DELIM ! COUNT
            write(funit,601,advance='no') mean,                        CSV_DELIM ! MEAN
            write(funit,602)              maxval(rmat, mask=mask)                ! MAX_INT
            call shells(iref)%set_rmat(rmat_shells, ft=.false.)
        end do
    end do
    call fclose(funit)

    ! Write output
    do iref=1, nrefs
        call cavgs(iref)%write(fn_cavgs_out, iref)
        call reprojs(iref)%write(fn_reprojs_out, iref)
        call diff(iref)%write(fn_diff, iref)
        call shells(iref)%write(fn_shells, iref)
    end do

end program simple_test_2D_core_finder

