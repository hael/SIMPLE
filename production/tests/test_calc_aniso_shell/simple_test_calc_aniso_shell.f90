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
    real, parameter             :: smpd=0.358, min_axis=0.38, max_axis=1.00
    real, allocatable           :: ellipsoids(:,:), centers(:,:), sim_rmat(:,:,:), sim_aniso(:,:,:)
    real                        :: eigenvecs(3,3), eigenvecs_inv(3,3), u, v, w, lhs, fit_aniso(3,3), &
                                    &fit_egnvals(3), fit_egnvecs(3,3), axes_error(3), angular_error(3), step
    integer, parameter          :: ldim(3)=(/160,160,160/), window=30
    integer, allocatable        :: sim_imat(:,:,:), nvox(:)
    integer                     :: i, j, k, cc, fu_fit, fu_results, natoms, icenter(3), errflg, nthr, &
                                    &aniso_in_pdb(3,3), ifoo
    character(len=256)          :: pdbin, nthrChar, junk1, junk2
    character(*), parameter     :: cc_out='sim_cc.mrc', vol_out='sim_vol.mrc', aniso_pdb='sim_aniso', &
        &fn_results='test_calc_aniso_shell_results.csv', fn_fit='aniso_bfac_field.pdb', &
        &aniso_fmt='(a28,6i7,a10)'

    character(len=*), parameter :: RESULTS_HEAD = 'INDEX'//CSV_DELIM//'NVOX'//CSV_DELIM//'EVAL_RESULT_1'//&
    &CSV_DELIM//'EVAL_RESULT_2'//CSV_DELIM//'EVAL_RESULT_3'//CSV_DELIM//'EVAL_TRUTH_1'//CSV_DELIM//'EVAL_TRUTH_2'&
    &//CSV_DELIM//'EVAL_TRUTH_3'//CSV_DELIM//'PCT_ERR_1'//CSV_DELIM//'PCT_ERR_2'//CSV_DELIM//'PCT_ERR_3'&
    &//CSV_DELIM//'EVEC_RESULT_1X'//CSV_DELIM//'EVEC_RESULT_1Y'//CSV_DELIM//'EVEC_RESULT_1Z'&
    &//CSV_DELIM//'EVEC_RESULT_2X'//CSV_DELIM//'EVEC_RESULT_2Y'//CSV_DELIM//'EVEC_RESULT_2Z'&
    &//CSV_DELIM//'EVEC_RESULT_3X'//CSV_DELIM//'EVEC_RESULT_3Y'//CSV_DELIM//'EVEC_RESULT_3Z'&
    &//CSV_DELIM//'EVEC_TRUTH_1X'//CSV_DELIM//'EVEC_TRUTH_1Y'//CSV_DELIM//'EVEC_TRUTH_1Z'&
    &//CSV_DELIM//'EVEC_TRUTH_2X'//CSV_DELIM//'EVEC_TRUTH_2Y'//CSV_DELIM//'EVEC_TRUTH_2Z'&
    &//CSV_DELIM//'EVEC_TRUTH_3X'//CSV_DELIM//'EVEC_TRUTH_3Y'//CSV_DELIM//'EVEC_TRUTH_3Z'&
    &//CSV_DELIM//'ANGLE_ERR_1'//CSV_DELIM//'ANGLE_ERR_2'//CSV_DELIM//'ANGLE_ERR_3'

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
    allocate(ellipsoids(natoms, 6), source=0.) ! 1:3 are major axes, 4:6 are rotation angles about x,y,z
    ! Test spheres
    do cc=1, natoms/3
        ellipsoids(cc,1:3) = min_axis + (max_axis-min_axis)*cc/(natoms/3)
        ellipsoids(cc,4:6) = 0.
    end do
    ! Test ellipses with principal axes along x,y,z
    step = (max_axis-min_axis)/3. ! Eval1 > Eval2 > Eval3
    do cc = 1, natoms/3
        ellipsoids(cc+natoms/3,1) = min_axis + (step)*cc/(natoms/3) + 2*step
        ellipsoids(cc+natoms/3,2) = min_axis + (step)*cc/(natoms/3) + step
        ellipsoids(cc+natoms/3,3) = min_axis + (step)*cc/(natoms/3)
        ellipsoids(cc+natoms/3,4:6) = 0.
    end do
    ! Test ellipses with various orientations
    do cc = 1, natoms/3
        ellipsoids(cc+2*natoms/3,1) = min_axis + (step)*cc/(natoms/3) + 2*step
        ellipsoids(cc+2*natoms/3,2) = min_axis + (step)*cc/(natoms/3) + step
        ellipsoids(cc+2*natoms/3,3) = min_axis + (step)*cc/(natoms/3)
        ellipsoids(cc+2*natoms/3,4) = PI/2*cos(cc*PI/(natoms/3))
        ellipsoids(cc+2*natoms/3,5) = PI/2*sin(cc*3*PI/(natoms/3))
        ellipsoids(cc+2*natoms/3,6) = -PI/2*cos(cc*5*PI/(natoms/3))
    end do

    ! Create a simulated CC map, volume, and ANISOU pdb file
    allocate(nvox(natoms), source=0)
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
                        nvox(cc) = nvox(cc) + 1
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
    deallocate(sim_imat, sim_rmat, sim_aniso)

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

    ! Calculate and output error
    write(logfhandle,'(a)') '>>> CALCULATING ERROR OF ELLIPTICAL FITS'
    call fopen(fu_fit, FILE=trim(fn_fit), STATUS='OLD', action='READ')
    call fopen(fu_results, FILE=trim(fn_results), STATUS='REPLACE', action='WRITE')
    write (fu_results, '(a)') RESULTS_HEAD
    aniso_in_pdb = 0
    do cc=1, natoms
        ! Read in aniso matrix
        read (fu_fit, *)  ! Skip the 'ATOM' lines.
        read (fu_fit, aniso_fmt) junk1, aniso_in_pdb(1,1), aniso_in_pdb(2,2), aniso_in_pdb(3,3),&
            &aniso_in_pdb(1,2), aniso_in_pdb(1,3), aniso_in_pdb(2,3), junk2
        fit_aniso = aniso_in_pdb/10000.
        fit_aniso(2,1) = fit_aniso(1,2)
        fit_aniso(3,1) = fit_aniso(1,3)
        fit_aniso(3,2) = fit_aniso(2,3)
        ! Calculate fit eigenvalues and eigenvecs
        call jacobi(fit_aniso, 3, 3, fit_egnvals, fit_egnvecs, ifoo)
        call eigsrt(fit_egnvals, fit_egnvecs, 3, 3)
        fit_egnvals = sqrt(fit_egnvals) ! ANISOU format has squared aniso vals
        ! Calculate percent error in eigenvalues and eigenvecs
        eigenvecs = find_eigenvecs(ellipsoids(cc,4), ellipsoids(cc,5), ellipsoids(cc,6))
        do i=1,3
            axes_error(i) = 100. * (fit_egnvals(i) - ellipsoids(cc, i)) / ellipsoids(cc, i)
            angular_error(i) = abs(dot_product(fit_egnvecs(:,i),eigenvecs(:,i)))
            angular_error(i) = angular_error(i) / norm_2(fit_egnvecs(:,i)) / norm_2(eigenvecs(:,i))
            ! Floating point arithmetic errors can cause values slightly greater than 1 to be passed to acos() 
            angular_error(i) = acos(min(1.,angular_error(i))) * 180. / PI ! Output in degrees
        end do
        ! Output into results CSV file
        call write_results(cc, fu_results, nvox(cc), fit_egnvals, ellipsoids(cc, 1:3), axes_error, fit_egnvecs, &
            &eigenvecs, angular_error)
    end do
    call fclose(fu_fit)
    call fclose(fu_results)

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

        subroutine write_results(cc, funit, size, eval_result, eval_truth, eval_err, evecs_result, &
                &evecs_truth, evecs_err)
            integer, intent(in)     :: cc, funit, size
            real, intent(in)        :: eval_result(3), eval_truth(3), eval_err(3), evecs_result(3,3), &
                                        &evecs_truth(3,3), evecs_err(3)
            
            601 format(F7.3,A2)
            602 format(F7.3)

            write(funit,601,advance='no') real(cc),                 CSV_DELIM ! INDEX
            write(funit,601,advance='no') real(size),               CSV_DELIM ! NVOX
            ! Eigenvalues measured by atoms_stats
            write(funit,601,advance='no') eval_result(1),           CSV_DELIM ! EVAL_RESULT_1
            write(funit,601,advance='no') eval_result(2),           CSV_DELIM ! EVAL_RESULT_2
            write(funit,601,advance='no') eval_result(3),           CSV_DELIM ! EVAL_RESULT_3
            ! Input eigenvalues
            write(funit,601,advance='no') eval_truth(1),            CSV_DELIM ! EVAL_TRUTH_1
            write(funit,601,advance='no') eval_truth(2),            CSV_DELIM ! EVAL_TRUTH_2
            write(funit,601,advance='no') eval_truth(3),            CSV_DELIM ! EVAL_TRUTH_3
            ! Pct error in eigenvalues
            write(funit,601,advance='no') eval_err(1),              CSV_DELIM ! EVAL_ERR_1
            write(funit,601,advance='no') eval_err(2),              CSV_DELIM ! EVAL_ERR_2
            write(funit,601,advance='no') eval_err(3),              CSV_DELIM ! EVAL_ERR_3
            ! Eigenvectors measured by atoms_stats
            write(funit,601,advance='no') evecs_result(1,1),        CSV_DELIM ! EVEC_RESULT_1X
            write(funit,601,advance='no') evecs_result(2,1),        CSV_DELIM ! EVEC_RESULT_1Y
            write(funit,601,advance='no') evecs_result(3,1),        CSV_DELIM ! EVEC_RESULT_1Z
            write(funit,601,advance='no') evecs_result(1,2),        CSV_DELIM ! EVEC_RESULT_2X
            write(funit,601,advance='no') evecs_result(2,2),        CSV_DELIM ! EVEC_RESULT_2Y
            write(funit,601,advance='no') evecs_result(3,2),        CSV_DELIM ! EVEC_RESULT_2Z
            write(funit,601,advance='no') evecs_result(1,3),        CSV_DELIM ! EVEC_RESULT_3X
            write(funit,601,advance='no') evecs_result(2,3),        CSV_DELIM ! EVEC_RESULT_3Y
            write(funit,601,advance='no') evecs_result(3,3),        CSV_DELIM ! EVEC_RESULT_3Z
            ! Input eigenvectors
            write(funit,601,advance='no') evecs_truth(1,1),         CSV_DELIM ! EVEC_TRUTH_1X
            write(funit,601,advance='no') evecs_truth(2,1),         CSV_DELIM ! EVEC_TRUTH_1Y
            write(funit,601,advance='no') evecs_truth(3,1),         CSV_DELIM ! EVEC_TRUTH_1Z
            write(funit,601,advance='no') evecs_truth(1,2),         CSV_DELIM ! EVEC_TRUTH_2X
            write(funit,601,advance='no') evecs_truth(2,2),         CSV_DELIM ! EVEC_TRUTH_2Y
            write(funit,601,advance='no') evecs_truth(3,2),         CSV_DELIM ! EVEC_TRUTH_2Z
            write(funit,601,advance='no') evecs_truth(1,3),         CSV_DELIM ! EVEC_TRUTH_3X
            write(funit,601,advance='no') evecs_truth(2,3),         CSV_DELIM ! EVEC_TRUTH_3Y
            write(funit,601,advance='no') evecs_truth(3,3),         CSV_DELIM ! EVEC_TRUTH_3Z
            ! Angular 
            write(funit,601,advance='no') evecs_err(1),             CSV_DELIM ! EVEC_ERR_1
            write(funit,601,advance='no') evecs_err(2),             CSV_DELIM ! EVEC_ERR_2
            write(funit,602)              evecs_err(3)                        ! EVEC_ERR_3
        end subroutine

end program simple_test_3D_opt_filt