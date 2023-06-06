! Subtracts reprojections from class averages.  
! Will want to also compare ptcl reprojections to ptcls
program simple_test_signal_subtr

    include 'simple_lib.f08'
    use simple_image,        only: image
    implicit none
    
    type(image),    allocatable :: cavgs(:), reprojs(:), out(:)
    real        :: smpd=0.358
    integer     :: iref, nrefs, ldim(3), ldim_refs(3), ifoo
    character(*), parameter    :: stks='cavgs_vs_reprojections_rec_and_sim.mrc'
    character(*), parameter    :: fn_out='cavgs_minus_reprojections.mrc'
    character(*), parameter    :: cavgs_out='cavgs_out.mrc', reprojs_out='reprojs_out.mrc'

    ! Read in stacks
    call find_ldim_nptcls(stks, ldim, ifoo)
    ldim_refs = [ldim(1), ldim(2), 1]
    nrefs = ldim(3) / 3
    allocate(cavgs(nrefs), reprojs(nrefs), out(nrefs))
    do iref = 1,nrefs
        call cavgs(iref)%new(ldim_refs, SMPD)
        call reprojs(iref)%new(ldim_refs, SMPD)
        call cavgs(iref)%read(stks, 3*iref-2)
        call reprojs(iref)%read(stks, 3*iref-1)
    end do

    ! Ensure proper normalizations
    ! I'm not sure if this is the correct way to normalize
    do iref=1, nrefs
        call cavgs(iref)%norm()
        call reprojs(iref)%norm()
    end do

    ! Subtract reprojs from cavgs
    do iref=1,nrefs
        call out(iref)%new(ldim_refs, SMPD)
        out(iref) = cavgs(iref) - reprojs(iref)
    end do

    ! Write output
    do iref=1, nrefs
        call cavgs(iref)%write(cavgs_out, iref)
        call reprojs(iref)%write(reprojs_out, iref)
        call out(iref)%write(fn_out, iref)
    end do

end program simple_test_signal_subtr

