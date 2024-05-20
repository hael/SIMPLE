program simple_test_nanoparticle
    
    include 'simple_lib.f08'
    use simple_nanoparticle
    use simple_parameters
    use simple_image

    implicit none

    type(nanoparticle) :: nano
    character(len=100) :: img_name, pdb_name_shared, pdb_name_sim, pdb_name_henry
    character(len=2)   :: element
    real :: smpd
    type(parameters), target :: params
    type(image) :: simatms

    element = 'Pt'
    smpd = 0.358
    params_glob => params
    params_glob%element = element
    params_glob%smpd = smpd
    img_name = 'rec_merged.mrc'
    pdb_name_shared = 'common_atoms_01in02.pdb'
    pdb_name_sim = 'experimental_centers_5.pdb'
    pdb_name_henry = 'rec_merged_ATMS.pdb'

    print *, 'Shared atoms'
    call nano%new(trim(img_name))
    call nano%set_atomic_coords(trim(pdb_name_shared))
    call nano%simulate_atoms(simatms=simatms)
    call nano%validate_atoms(simatms=simatms)
    print *, 'New method'
    call nano%new(trim(img_name))
    call nano%set_atomic_coords(trim(pdb_name_sim))
    call nano%simulate_atoms(simatms=simatms)
    call nano%validate_atoms(simatms=simatms)
    print *, 'Old method'
    call nano%new(trim(img_name))
    call nano%set_atomic_coords(trim(pdb_name_henry))
    call nano%simulate_atoms(simatms=simatms)
    call nano%validate_atoms(simatms=simatms)

end program simple_test_nanoparticle