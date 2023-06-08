! Subtracts reprojections from class averages.  
! Will want to also compare ptcl reprojections to ptcls
program simple_test_signal_subtr

    include 'simple_lib.f08'
    use simple_image,        only: image
    implicit none
    
    type(image),    allocatable :: cavgs(:), reprojs(:), out(:)
    real        :: smpd=0.358, nprad_A=15, nprad_vox
    integer     :: iref, nrefs, ldim1(3), ldim2(3), ldim_refs(3), ifoo
    character(*), parameter    :: fn_cavgs='cavgs.mrc', fn_reprojs='reprojs.mrcs'
    character(*), parameter    :: fn_out='cavgs_minus_reprojections.mrc'
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
    allocate(cavgs(nrefs), reprojs(nrefs), out(nrefs))
    do iref = 1,nrefs
        call cavgs(iref)%new(ldim_refs, SMPD)
        call reprojs(iref)%new(ldim_refs, SMPD)
        call cavgs(iref)%read(fn_cavgs, iref)
        call reprojs(iref)%read(fn_reprojs, iref)
    end do

    ! Ensure proper normalizations of cavgs
    ! This is wrong
    do iref=1, nrefs
         call cavgs(iref)%norm()
         call reprojs(iref)%norm()
     end do

    ! Subtract reprojs from cavgs
    do iref=1,nrefs
        call out(iref)%new(ldim_refs, SMPD)
        out(iref) = cavgs(iref) - reprojs(iref)
    end do

    ! Take averages of shells out to the NP radius
    !--INSERT HERE---

    ! Write output
    do iref=1, nrefs
        call cavgs(iref)%write(cavgs_out, iref)
        call reprojs(iref)%write(reprojs_out, iref)
        call out(iref)%write(fn_out, iref)
    end do

end program simple_test_signal_subtr

