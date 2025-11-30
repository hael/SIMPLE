program simple_test_nanoparticle
    
    ! include 'simple_lib.f08'
    ! use simple_nanoparticle
    ! use simple_parameters
    ! use simple_image

    ! implicit none

    ! type(nanoparticle) :: nano
    ! character(len=100) :: img_name, pdb_name_shared, pdb_name_sim, pdb_name_henry, pdb_name_sim_not_henry, pdb_name_henry_not_sim
    ! character(len=2)   :: element
    ! real :: smpd
    ! type(parameters), target :: params
    ! type(image) :: simatms

    ! element = 'Pt'
    ! smpd = 0.358
    ! params_glob => params
    ! params_glob%element = element
    ! params_glob%smpd = smpd
    ! img_name = 'rec_merged.mrc'
    ! pdb_name_shared = 'common_atoms_01in02.pdb'
    ! pdb_name_sim = 'experimental_centers_6.pdb'
    ! pdb_name_henry = 'rec_merged_ATMS.pdb'
    ! pdb_name_sim_not_henry = 'different_atoms_01not_in02.pdb'
    ! pdb_name_henry_not_sim = 'different_atoms_02not_in01.pdb'

    ! print *, 'Shared atoms'
    ! call nano%new(trim(img_name))
    ! call nano%set_atomic_coords(trim(pdb_name_shared))
    ! call nano%simulate_atoms(simatms)
    ! call nano%validate_atoms(simatms, l_print=.true.)
    ! call nano%kill
    ! print *, 'New method'
    ! call nano%new(trim(img_name))
    ! call nano%set_atomic_coords(trim(pdb_name_sim))
    ! call nano%simulate_atoms(simatms)
    ! call nano%validate_atoms(simatms, l_print=.true.)
    ! call nano%kill
    ! print *, 'Old method'
    ! call nano%new(trim(img_name))
    ! call nano%set_atomic_coords(trim(pdb_name_henry))
    ! call nano%simulate_atoms(simatms)
    ! call nano%validate_atoms(simatms, l_print=.true.)
    ! call nano%kill
    ! print *, 'Atoms in old method and not new'
    ! call nano%new(trim(img_name))
    ! call nano%set_atomic_coords(trim(pdb_name_henry_not_sim))
    ! call nano%simulate_atoms(simatms)
    ! call nano%validate_atoms(simatms, l_print=.true.)
    ! call nano%kill
    ! print *, 'Atoms in new method and not old'
    ! call nano%new(trim(img_name))
    ! call nano%set_atomic_coords(trim(pdb_name_sim_not_henry))
    ! call nano%simulate_atoms(simatms)
    ! call nano%validate_atoms(simatms, l_print=.true.)
    ! call nano%kill

end program simple_test_nanoparticle