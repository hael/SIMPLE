program simple_test_3D_opt_filt
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters
    use simple_commander_atoms,    only: atoms_stats_commander
    use simple_image,              only: image
    use simple_binimage,           only: binimage
    use simple_atoms,              only: atoms
    implicit none

    type (cmdline)              :: cline, cline_atoms_stats
    type (parameters)           :: p
    type(binimage)              :: simcc
    type(image)                 :: simvol
    type(atoms)                 :: simatoms
    type(atoms_stats_commander) :: xatoms_stats
    character(*), parameter     :: cc_out='sim_cc.mrc', vol_out='sim_vol.mrc', aniso_pdb='sim_aniso'
    character(len=256)          :: pdbin, nthrChar
    real, parameter             :: smpd = 0.358
    real, allocatable           :: ellipsoids(:,:), centers(:,:), sim_rmat(:,:,:), sim_aniso(:,:,:)
    real                        :: eigenvecs(3,3), eigenvecs_inv(3,3), u, v, w, lhs
    integer, parameter          :: ldim(3)=(/160,160,160/), window=30
    integer, allocatable        :: sim_imat(:,:,:)
    integer                     :: i, j, k, cc, funit, natoms, icenter(3), errflg, nthr

    if( command_argument_count() /= 2 )then
        write(logfhandle,'(a)') 'Usage: simple_test_calc_aniso_shell file.pdb nthr'
        write(logfhandle,'(a)') 'file.pdb contains the input atomic elements and coordinates.'
        write(logfhandle,'(a)') 'nthr [integer] for atoms_stats'
        stop
    else
        call get_command_argument(1, pdbin)
        call get_command_argument(2, nthrChar)
        read (nthrChar, '(i3)') nthr
    endif

    ! Read in centers from pdbfile
    call simatoms%new(pdbin)
    natoms = simatoms%get_n()

    ! Create simulated CC eigenvalues and eigenvectors
    allocate(ellipsoids(natoms, 6), source=0.)
    do cc=1, 20
        ellipsoids(cc,1) = 1.0
        ellipsoids(cc,2) = 1.0
        ellipsoids(cc,3) = 1.0
        ellipsoids(cc,4) = PI/2
        ellipsoids(cc,5:6) = 0.
    end do
    do cc=21, 40
        ellipsoids(cc,1) = 0.8
        ellipsoids(cc,2) = 0.8
        ellipsoids(cc,3) = 0.8
        ellipsoids(cc,4) = 0.0
        ellipsoids(cc,5) = PI/2
        ellipsoids(cc,6) = 0.0
    end do
    do cc=41, 60
        ellipsoids(cc,1) = 0.6
        ellipsoids(cc,2) = 0.6
        ellipsoids(cc,3) = 0.6
        ellipsoids(cc,4) = 0.0
        ellipsoids(cc,5) = 0.0
        ellipsoids(cc,6) = PI/2
    end do
    do cc=61, natoms
        ellipsoids(cc,1) = 0.4
        ellipsoids(cc,2) = 0.4
        ellipsoids(cc,3) = 0.4
        ellipsoids(cc,4:6) = 0.
    end do

    ! Create a simulated CC map, volume, and ANISOU pdb file
    allocate(sim_imat(ldim(1), ldim(2), ldim(3)), source=0)
    allocate(sim_rmat(ldim(1), ldim(2), ldim(3)), source=0.)
    allocate(sim_aniso(3,3,natoms), source=0.) ! 3x3 anisotropic displacment matrix
    do cc=1, natoms
        icenter=nint(simatoms%get_coord(cc)/smpd) + 1 ! Voxel indeces start at 1 while Angstroms start at 0
        eigenvecs = find_eigenvecs(ellipsoids(cc,4), ellipsoids(cc,5), ellipsoids(cc,6))
        ! Identify all voxels within the ellipsoid associated with this CC.
        do k=icenter(3)-window, icenter(3)+window
            do j=icenter(2)-window, icenter(2)+window
                do i=icenter(1)-window, icenter(1)+window
                    u = dot_product((/i,j,k/)-icenter, eigenvecs(:,1))*smpd
                    v = dot_product((/i,j,k/)-icenter, eigenvecs(:,2))*smpd
                    w = dot_product((/i,j,k/)-icenter, eigenvecs(:,3))*smpd
                    ! Conditions of ellipsoids
                    if ( sum(((/u,v,w/)/ellipsoids(cc,1:3))**2) <= 1. ) then
                        sim_imat(i,j,k) = cc
                        sim_rmat(i,j,k) = 1.
                    end if
                end do
            end do
        end do
        ! Calculated simulated anisotropic displacement matrix to view thermal ellipsoids in Chimera
        call matinv(eigenvecs, eigenvecs_inv, 3, errflg)
        do i=1,3
            sim_aniso(i,i,cc) = ellipsoids(cc,i)**2 ! Aniso format uses sqared displacments
        end do
        sim_aniso(:,:,cc) = matmul(matmul(eigenvecs, sim_aniso(:,:,cc)), eigenvecs_inv) ! (u,v,w)->(x,y,z)
    end do
    call simcc%new_bimg(ldim, smpd)
    call simcc%set_imat(sim_imat)
    call simcc%write_bimg(cc_out)
    write(logfhandle,'(a)') '>>> SIMULATED CC MAP: '//cc_out
    call simvol%new(ldim, smpd)
    call simvol%set_rmat(sim_rmat, .false.)
    call simvol%write(vol_out)
    write(logfhandle,'(a)') '>>> SIMULATED VOLUME: '//vol_out
    write(logfhandle,'(a)') '>>> OUTPUT SIMULATED ANISOU PDB: '//aniso_pdb//'.pdb'
    call simatoms%writepdb_aniso(aniso_pdb, sim_aniso)
    deallocate(ellipsoids, sim_imat, sim_rmat, sim_aniso)

    ! Pass to atoms stats
    write(logfhandle,'(a)') '>>> PASSING ELLIPTICAL ATOMS TO ATOMS_STATS'
    call cline_atoms_stats%set('vol1', vol_out)
    call cline_atoms_stats%set('vol2', cc_out)
    call cline_atoms_stats%set('pdbfile', pdbin)
    call cline_atoms_stats%set('smpd', smpd)
    call cline_atoms_stats%set('element', simatoms%get_element(1))
    call cline_atoms_stats%set('mskdiam', 40.)
    call cline_atoms_stats%set('nthr', 1.*nthr)
    call xatoms_stats%execute(cline_atoms_stats)

    contains

        pure function find_eigenvecs(a,b,c) result(eigenvecs)
            real, intent(in)    :: a, b, c ! Rotation angles in radians
            real                :: rotx(3,3), roty(3,3), rotz(3,3), eigenvecs(3,3)
            integer             :: i

            ! Rotation about x-axis
            rotx = 0.
            rotx(1,1) = 1.
            rotx(2,2) = cos(a)
            rotx(3,3) = cos(a)
            rotx(2,3) = -sin(a)
            rotx(3,2) = sin(a)

            ! Rotation about y-axis
            roty=0.
            roty(1,1) = cos(b)
            roty(2,2) = 1.
            roty(3,3) = cos(b)
            roty(1,3) = sin(b)
            roty(3,1) = -sin(b)

            ! Rotation about z-axis
            rotz = 0.
            rotz(1,1) = cos(c)
            rotz(2,2) = cos(c)
            rotz(1,2) = -sin(c)
            rotz(2,1) = sin(c)
            rotz(3,3) = 1.

            eigenvecs = 0.
            do i=1,3
                eigenvecs(i,i) = 1.
            end do
            eigenvecs = matmul(rotz, matmul(roty, matmul(rotx, eigenvecs)))
        end function find_eigenvecs

end program simple_test_3D_opt_filt