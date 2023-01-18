! This program tests the calculations of lattice displacements (DISPL in atoms_stats.csv) output
! by atoms_stats. It takes in a pdb file which is used to find an ideal lattice and generate
! an ideal lattice with random displacements.  An integer nthr is taken as a command line arg
! for passing to atoms_stats.  The error in lattice displacement calculations is written to the file
! test_lattice_displ.csv.  Author: Henry Wietfeldt
program simple_test_lattice_displ
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters, params_glob
    use simple_commander_atoms,    only: atoms_stats_commander
    use simple_binimage,           only: binimage
    use simple_image,              only: image
    use simple_atoms,              only: atoms
    implicit none

    type(cmdline)                 :: cline_atoms_stats
    type(atoms_stats_commander) :: xatoms_stats
    type(atoms)                 :: atomsin
    type(binimage)              :: ccmap
    type(image)                 :: vol
    real, allocatable       :: displ(:,:), displ_out(:), rmat(:,:,:), centers(:,:)
    real, parameter         :: displ_mean=0., displ_std=0.3, limits(2)=[-0.7, 0.7], smpd=0.358
    real                    :: a0, rad=0.6/smpd, distTemp, dist, cMin(3), cMax(3), cMid(3)
    integer, allocatable    :: imat(:,:,:), lat_pos(:,:)
    integer, parameter      :: window=10
    integer                 :: funit, natoms, cc, i, j, k, icenter(3), ldim(3)=[160,160,160], centerAtom, nthr
    character(len=2)        :: element = 'Pt', el_ucase
    character(len=8)        :: crystal_system, strDispl, strNthr
    character(len=256)      :: line, pdbin
    character(*), parameter :: pdbsim='sim_displacements', fnvol='dummy_vol.mrc', fncc='dummy_cc.mrc', &
                               &fnas='atoms_stats.csv', fnresults='test_lattice_displ.csv'
    character(*), parameter :: header='INDEX'//CSV_DELIM//'DISPLIN_A'//CSV_DELIM//'DISPLIN_VOX'//CSV_DELIM//&
                               &'DISPLOUT_A'//CSV_DELIM//'DISPLOUT_VOX'//CSV_DELIM//'DIFF_A'//CSV_DELIM//&
                               &'DIFF_VOX'
    
    write(logfhandle,'(a)') '>>> Dev Note: This test is under development'
    
    if( command_argument_count() /= 2 )then
        write(logfhandle,'(a)') 'Usage: simple_test_lattice_displ file.pdb nthr'
        write(logfhandle,'(a)') 'file.pdb contains the input atomic elements and coordinates.'
        write(logfhandle,'(a)') 'nthr [integer] for atoms_stats'
        stop
    else
        call get_command_argument(1, pdbin)
        call get_command_argument(2, strNthr)
        read (strNthr, '(i3)') nthr
    endif

    ! Read info from PDB
    write(logfhandle,'(a)') '>>> FINDING IDEAL LATTICE'
    call atomsin%new(pdbin)
    natoms = atomsin%get_n()
    print *, natoms
    allocate(centers(natoms, 3), source=0.)
    el_ucase = uppercase(trim(adjustl(atomsin%get_element(1))))
    call get_lattice_params(el_ucase, crystal_system, a0)
    do cc=1, natoms
        centers(cc, 1:3) = atomsin%get_coord(cc)
    end do
    ! Get center atom
    centerAtom = 0
    distTemp = HUGE(distTemp)
    cMin       = minval(centers,dim=1)
    cMax       = maxval(centers,dim=1)
    cMid       = (cMax-cMin)/2.+cMin
    do cc=1, natoms
        dist = euclid(cMid, centers(cc,1:3))
        if (dist < distTemp) then
            centerAtom = cc
            distTemp = dist
        end if
    end do
    ! Generate ideal lattice
    allocate(lat_pos(natoms, 3), source=0)
    do cc=1, natoms
        lat_pos(cc, 1:3) = nint((centers(cc,1:3) - centers(centerAtom,1:3)) / (a0 / 2.))
    end do

    ! Simulate random disiplacements
    write(logfhandle,'(a)') '>>> SIMULATING DISPLACEMENTS'
    allocate(displ(natoms, 3), source=0.)
    do cc=1, natoms
        do i=1,3
            ! Comment out the following line to ensure an ideal lattice was generated
            displ(cc, i) = gasdev(displ_mean, displ_std, limits)
        end do
        call atomsin%set_coord(cc, a0/2.*lat_pos(cc,1:3) + centers(centerAtom,1:3) + displ(cc,1:3))
    end do

    write(logfhandle,'(a)') '>>> GENERATING SIMULATED VOLUMES'
    ! Generate dummy volume and CC map (spherical atoms)
    allocate(rmat(ldim(1), ldim(2), ldim(3)), source=0.)
    allocate(imat(ldim(1), ldim(2), ldim(3)), source=0)
    do cc=1, natoms
        icenter = nint(atomsin%get_coord(cc) / smpd) + 1
        do k=icenter(3)-window, icenter(3)+window
            do j=icenter(2)-window, icenter(2)+window
                do i=icenter(1)-window, icenter(1)+window
                    if (euclid(1.*[i,j,k], 1.*icenter) < rad) then
                        imat(i,j,k) = cc
                        rmat(i,j,k) = 1.
                    end if
                end do
            end do
        end do
    end do

    ! Output volumes and pdb file
    call vol%new(ldim, smpd)
    call vol%set_rmat(rmat, .false.)
    call vol%write(fnvol)
    call ccmap%new_bimg(ldim, smpd)
    call ccmap%set_imat(imat)
    call ccmap%write_bimg(fncc)
    call atomsin%writepdb(pdbsim)

    ! Input pdb file and dummy volumes into atoms_stats
    write(logfhandle,'(a)') '>>> PASSING TO ATOMS_STATS'
    call cline_atoms_stats%set('vol1', fnvol)
    call cline_atoms_stats%set('vol2', fncc)
    call cline_atoms_stats%set('pdbfile', pdbsim//'.pdb')
    call cline_atoms_stats%set('smpd', smpd)
    call cline_atoms_stats%set('element', element)
    call cline_atoms_stats%set('mskdiam', 40.)
    call cline_atoms_stats%set('nthr', 1.*nthr)
    call xatoms_stats%execute(cline_atoms_stats)
    deallocate(rmat, imat)
    call vol%kill
    call ccmap%kill

    ! Read DISPL from atoms_stats CSV file
    allocate(displ_out(natoms), source=0.)
    call fopen(funit, FILE=trim(fnas), STATUS='OLD', action='READ')
    read(funit, *) line ! Header
    do cc=1, natoms
        read(funit, '(a)') line
        strDispl = trim(line(101:108))
        read(strDispl, '(f8.4)') displ_out(cc)
    end do
    call fclose(funit)

    ! Output the sim, fit, and differences
    write(logfhandle,'(a)') '>>> CALCULATING ERROR'
    call fopen(funit, FILE=trim(fnresults), STATUS='NEW', action='WRITE')
    write(funit, '(a)') header
    601 format(F7.3,A2)
    602 format(F7.3)
    do cc=1, natoms
        write(funit, 601, advance='no') real(cc),                                   CSV_DELIM   ! INDEX
        write(funit, 601, advance='no') norm_2(displ(cc,1:3)),                      CSV_DELIM   ! DISPLIN_A
        write(funit, 601, advance='no') norm_2(displ(cc,1:3))/smpd,                 CSV_DELIM   ! DISPLIN_VOX
        write(funit, 601, advance='no') displ_out(cc),                              CSV_DELIM   ! DISPLOUT_A
        write(funit, 601, advance='no') displ_out(cc)/smpd,                         CSV_DELIM   ! DISPLOUT_VOX
        write(funit, 601, advance='no') (norm_2(displ(cc,1:3))-displ_out(cc)),      CSV_DELIM   ! DIFF_A
        write(funit, 602)               (norm_2(displ(cc,1:3))-displ_out(cc))/smpd              ! DIFF_VOX
    end do
    call fclose(funit)

    deallocate(displ, displ_out)
end program simple_test_lattice_displ