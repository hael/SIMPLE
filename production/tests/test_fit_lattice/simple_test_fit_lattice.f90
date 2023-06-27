! This program tests the calculations of lattice displacements (DISPL in atoms_stats.csv) output
! by atoms_stats. It takes in a pdb file which is used to find an ideal lattice and generate
! an ideal lattice with random displacements.  An integer nthr is taken as a command line arg
! for passing to atoms_stats.  The error in lattice displacement calculations is written to the file
! test_lattice_displ.csv.  Author: Henry Wietfeldt
program simple_test_fit_lattice
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters, params_glob
    use simple_atoms,              only: atoms
    use simple_nanoparticle_utils, only: fit_lattice
    implicit none

    type(cmdline)                 :: cline_atoms_stats
    type(atoms)                 :: atomsin
    real, allocatable       :: displ(:,:), displ_out(:), centers(:,:), ideal(:,:)
    real, parameter         :: displ_mean=0., limits(2)=[-0.7, 0.7], aMin=3.4, aMax=5.0, aStep=0.2
    real                    :: a0, a(3), a_fit(3), distTemp, dist, cMin(3), cMax(3), cMid(3), displ_std
    integer, allocatable    :: lat_pos(:,:)
    integer                 :: funit, natoms, cc, i, j, k, icenter(3), centerAtom, nsteps
    character(len=2)        :: element = 'Pt', el_ucase
    character(len=8)        :: crystal_system, strDispl
    character(len=256)      :: line, pdbin
    
    write(logfhandle,'(a)') '>>> Dev Note: This test is under development'
    
    if( command_argument_count() /= 1 )then
        write(logfhandle,'(a)') '>>> Usage: simple_test_fit_lattice file.pdb'
        write(logfhandle,'(2a)') 'file.pdb contains the input atomic elements and coordinates.', NEW_LINE('a')
        stop
    else
        call get_command_argument(1, pdbin)
    endif

    ! Read info from PDB (This allows us to test lattice fitting on lattices with gaps)
    write(logfhandle,'(a)') '>>> FINDING IDEAL LATTICE'
    call atomsin%new(pdbin)
    natoms = atomsin%get_n()
    write(logfhandle,'(a, i5, a)') '>>> Lattice contains ', natoms, ' atoms.'
    allocate(centers(3, natoms), source=0.)
    el_ucase = uppercase(trim(adjustl(atomsin%get_element(1))))
    call get_lattice_params(el_ucase, crystal_system, a0)
    do cc=1, natoms
        centers(1:3, cc) = atomsin%get_coord(cc)
    end do
    ! Get center atom
    centerAtom = 0
    distTemp = HUGE(distTemp)
    cMin       = minval(centers,dim=2)
    cMax       = maxval(centers,dim=2)
    cMid       = (cMax-cMin)/2.+cMin
    do cc=1, natoms
        dist = euclid(cMid, centers(1:3,cc))
        if (dist < distTemp) then
            centerAtom = cc
            distTemp = dist
        end if
    end do

    ! Generate ideal lattice
    allocate(lat_pos(3, natoms), source=0)
    do cc=1, natoms
        lat_pos(1:3,cc) = nint((centers(1:3,cc) - centers(1:3,centerAtom)) / (a0 / 2.))
        ! Get real positions
    end do

    ! Do this for a range of lattice params (3.5 to 4.5 in steps of 0.2)
    allocate(ideal(3, natoms), displ(3, natoms), source=0.)
    a = aMin
    do while(a(1) <= aMax)

        ! Get ideal positions from lattice params
        write(logfhandle,'(a)') '>>> FITTING IDEAL LATTICE'
        do cc=1, natoms
            ideal(1:3,cc) = centers(centerAtom,1:3) + a/2.*lat_pos(1:3,cc)
        end do

        ! Fit lattice (call directly)
        call fit_lattice(el_ucase, ideal, a_fit)
        print *,     "a: ", a
        print *, " a_fit:", a_fit

        ! Simulate random disiplacements
        displ_std = a(1) / 20. ! So different mean(a) with displacements are comparable
        write(logfhandle,'(a)') '>>> FITTING LATTICE WITH DISPLACEMENTS'
        do cc=1, natoms
            displ(1:3,cc) = ideal(1:3,cc) + gasdev(displ_mean, displ_std)
        end do

        ! Fit lattice again
        call fit_lattice(el_ucase, displ, a_fit)
        print *,     "a: ", a
        print *, " a_fit:", a_fit

        ! Increase lattice parameter
        a = a + aStep

    end do  

end program simple_test_fit_lattice