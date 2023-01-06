program simple_test_calc_aniso_shell
    include 'simple_lib.f08'
    use simple_cmdline,            only: cmdline
    use simple_builder,            only: builder
    use simple_parameters,         only: parameters, params_glob
    use simple_commander_atoms,    only: atoms_stats_commander
    use simple_image,              only: image
    use simple_binimage,           only: binimage
    use simple_atoms,              only: atoms
    implicit none

    type (cmdline)              :: cline, cline_atoms_stats, cline_atoms_stats_ccheck
    type (parameters)           :: p
    type(binimage)              :: simcc, fitcc
    type(image)                 :: simvol, fitvol
    type(atoms)                 :: simatoms
    type(atoms_stats_commander) :: xatoms_stats, xatoms_stats_ccheck
    real, parameter             :: smpd=0.358, min_axis=0.38, max_axis=1.00
    real, allocatable           :: ellipsoids(:,:), centers(:,:), rmat(:,:,:), sim_aniso(:,:,:)
    real                        :: eigenvecs(3,3), eigenvecs_inv(3,3), u, v, w, lhs,&
                                    &fit_egnvals(3), fit_egnvecs(3,3), twice_fit_evals(3), twice_fit_evecs(3,3), &
                                    &axes_error(3), angular_error(3), xyzdisp(3), fit_xyzdisp(3), twice_fit_xyzdisp(3), step
    integer, parameter          :: ldim(3)=(/160,160,160/), window=30
    integer, allocatable        :: imat(:,:,:), nvox(:), fit_nvox(:)
    integer                     :: i, j, k, cc, fu_fit, fu_twice_fit, fu_results, natoms, icenter(3), errflg, nthr, &
                                    &ifoo
    logical, parameter          :: check_consistency = .true.
    character(len=256)          :: pdbin, nthrChar
    character(*), parameter     :: cc_out='sim_cc.mrc', vol_out='sim_vol.mrc', aniso_pdb='sim_aniso', &
        &fn_results='test_calc_aniso_shell_results.csv', fn_fit='aniso_bfac_field.pdb', &
        &aniso_fmt='(a28,6i7,a10)', fn_fitcc='fit_cc.mrc', fn_fitvol='fit_vol.mrc', &
        &consistency_dir='test_consistency/'

    character(len=*), parameter :: RESULTS_HEAD = 'INDEX'//CSV_DELIM//'NVOX'//CSV_DELIM//'EVAL_RESULT_1'//&
    &CSV_DELIM//'EVAL_RESULT_2'//CSV_DELIM//'EVAL_RESULT_3'//CSV_DELIM//'EVAL_TRUTH_1'//CSV_DELIM//'EVAL_TRUTH_2'&
    &//CSV_DELIM//'EVAL_TRUTH_3'//CSV_DELIM//'PCT_ERR_1'//CSV_DELIM//'PCT_ERR_2'//CSV_DELIM//'PCT_ERR_3'&
    &//CSV_DELIM//'EVEC_RESULT_1X'//CSV_DELIM//'EVEC_RESULT_1Y'//CSV_DELIM//'EVEC_RESULT_1Z'&
    &//CSV_DELIM//'EVEC_RESULT_2X'//CSV_DELIM//'EVEC_RESULT_2Y'//CSV_DELIM//'EVEC_RESULT_2Z'&
    &//CSV_DELIM//'EVEC_RESULT_3X'//CSV_DELIM//'EVEC_RESULT_3Y'//CSV_DELIM//'EVEC_RESULT_3Z'&
    &//CSV_DELIM//'EVEC_TRUTH_1X'//CSV_DELIM//'EVEC_TRUTH_1Y'//CSV_DELIM//'EVEC_TRUTH_1Z'&
    &//CSV_DELIM//'EVEC_TRUTH_2X'//CSV_DELIM//'EVEC_TRUTH_2Y'//CSV_DELIM//'EVEC_TRUTH_2Z'&
    &//CSV_DELIM//'EVEC_TRUTH_3X'//CSV_DELIM//'EVEC_TRUTH_3Y'//CSV_DELIM//'EVEC_TRUTH_3Z'&
    &//CSV_DELIM//'ANGLE_ERR_1'//CSV_DELIM//'ANGLE_ERR_2'//CSV_DELIM//'ANGLE_ERR_3'&
    &//CSV_DELIM//'XRESULT'//CSV_DELIM//'YRESULT'//CSV_DELIM//'ZRESULT'&
    &//CSV_DELIM//'XTRUTH'//CSV_DELIM//'YTRUTH'//CSV_DELIM//'ZTRUTH'

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
    allocate(imat(ldim(1), ldim(2), ldim(3)), source=0)
    allocate(rmat(ldim(1), ldim(2), ldim(3)), source=0.)
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
                        imat(i,j,k) = cc
                        rmat(i,j,k) = 1.
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
    call simcc%set_imat(imat)
    call simcc%write_bimg(cc_out)
    write(logfhandle,'(a)') '>>> SIMULATED CC MAP: '//cc_out
    call simvol%new(ldim, smpd)
    call simvol%set_rmat(rmat, .false.)
    call simvol%write(vol_out)
    write(logfhandle,'(a)') '>>> SIMULATED VOLUME: '//vol_out
    write(logfhandle,'(a)') '>>> OUTPUT SIMULATED ANISOU PDB: '//aniso_pdb//'.pdb'
    call simatoms%writepdb_aniso(aniso_pdb, sim_aniso)
    deallocate(imat, rmat, sim_aniso)

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

    ! Calculate and output error.  Also prepare for consistency check if applicable
    write(logfhandle,'(a)') '>>> CALCULATING ERROR OF ELLIPTICAL FITS'
    call fopen(fu_fit, FILE=trim(fn_fit), STATUS='OLD', action='READ')
    call fopen(fu_results, FILE=trim(fn_results), STATUS='REPLACE', action='WRITE')
    write (fu_results, '(a)') RESULTS_HEAD
    if (check_consistency) then
        allocate(fit_nvox(natoms), source=0)
        allocate(imat(ldim(1), ldim(2), ldim(3)), source=0)
        allocate(rmat(ldim(1), ldim(2), ldim(3)), source=0.)
    end if
    do cc=1, natoms
        call read_from_aniso_pdb(fu_fit, fit_egnvals, fit_egnvecs)
        eigenvecs = find_eigenvecs(ellipsoids(cc,4), ellipsoids(cc,5), ellipsoids(cc,6))
        ! Get truth displacement in x, y, z direction from truth eigenvals and eigenvecs
        xyzdisp = get_xyzdisp(eigenvecs, ellipsoids(cc, 1:3))
        fit_xyzdisp = get_xyzdisp(fit_egnvecs, fit_egnvals)
        ! Calculate eigenvalue and eigevector error
        do i=1,3
            axes_error(i) = 100. * (fit_egnvals(i) - ellipsoids(cc, i)) / ellipsoids(cc, i) ! Signed PCT error
            angular_error(i) = angular_diff(fit_egnvecs(:,i), eigenvecs(:,i))
        end do
        ! Output into results CSV file
        call write_results(cc, fu_results, nvox(cc), fit_egnvals, ellipsoids(cc, 1:3), axes_error, fit_egnvecs, &
            &eigenvecs, angular_error, fit_xyzdisp, xyzdisp)
        if (check_consistency) then
            ! Create the fit imat and rmat to pass to atoms_stats and second time
            icenter=nint(simatoms%get_coord(cc)/smpd) + 1
            do k=icenter(3)-window, icenter(3)+window
                do j=icenter(2)-window, icenter(2)+window
                    do i=icenter(1)-window, icenter(1)+window
                        u = dot_product((/i,j,k/)-icenter, fit_egnvecs(:,1))*smpd
                        v = dot_product((/i,j,k/)-icenter, fit_egnvecs(:,2))*smpd
                        w = dot_product((/i,j,k/)-icenter, fit_egnvecs(:,3))*smpd
                        ! Conditions of ellipsoids
                        if ( sum(((/u,v,w/)/fit_egnvals(:))**2) <= 1. ) then
                            imat(i,j,k) = cc
                            rmat(i,j,k) = 1.
                            fit_nvox(cc) = fit_nvox(cc) + 1
                        end if
                    end do
                end do
            end do
        end if
    end do
    call fclose(fu_fit)
    call fclose(fu_results)

    ! Check how well atoms_stats can reproduce the fit ellipsoids from the above call to atoms_stats
    if (check_consistency) then
        write(logfhandle,'(a)') '>>> CHECKING CONSISTENCY OF ELLIPTICAL FITTING'
        ! Create fit volume and CC maps
        call simple_mkdir(consistency_dir)
        call simple_chdir(consistency_dir)
        write(logfhandle,'(a)') '>>> DIRECTORY FOR CONSISTENCY TEST: '
        call exec_cmdline('pwd')
        call fitcc%new_bimg(ldim, smpd)
        call fitcc%set_imat(imat)
        call fitcc%write_bimg(fn_fitcc)
        call fitvol%new(ldim, smpd)
        call fitvol%set_rmat(rmat, .false.)
        call fitvol%write(fn_fitvol)
        ! Pass to atoms_stats a second time
        write(logfhandle,'(a)') '>>> CONSISTENCY TEST: PASSING RESULTS TO ATOMS_STATS'
        call cline_atoms_stats%set('vol1', fn_fitvol)
        call cline_atoms_stats%set('vol2', fn_fitcc)
        call cline_atoms_stats%set('pdbfile', '../'//pdbin)
        call cline_atoms_stats%printline()
        params_glob => null() ! So that new global params can be set from cline
        call xatoms_stats%execute(cline_atoms_stats)

        ! Calculate Error and output in CSV file
        write(logfhandle,'(a)') '>>> CONSISTENCY TEST: CALCULATING ERROR OF ELLIPTICAL FITS'
        call fopen(fu_fit, FILE=trim('../'//fn_fit), STATUS='OLD', action='READ')
        call fopen(fu_twice_fit, FILE=trim(fn_fit), STATUS='OLD', action='READ')
        call fopen(fu_results, FILE=trim('consistency_'//fn_results), STATUS='REPLACE', action='WRITE')
        write (fu_results, '(a)') RESULTS_HEAD
        do cc=1, natoms
            call read_from_aniso_pdb(fu_fit, fit_egnvals, fit_egnvecs)
            call read_from_aniso_pdb(fu_twice_fit, twice_fit_evals, twice_fit_evecs)
            fit_xyzdisp = get_xyzdisp(fit_egnvecs, fit_egnvals)
            twice_fit_xyzdisp = get_xyzdisp(twice_fit_evecs, twice_fit_evals)
            do i=1,3
                axes_error(i) = 100. * (twice_fit_evals(i) - fit_egnvals(i)) / twice_fit_evals(i) ! PCT error
                angular_error(i) = angular_diff(twice_fit_evecs(:,i), fit_egnvecs(:,i))
            end do
            call write_results(cc, fu_results, nvox(cc), twice_fit_evals, fit_egnvals, axes_error, twice_fit_evecs, &
                fit_egnvecs, angular_error, fit_xyzdisp, xyzdisp)
        end do
        call fclose(fu_fit)
        call fclose(fu_twice_fit)
        call fclose(fu_results)
    end if
    write(logfhandle,'(a)') '****SIMPLE_TEST_CALC_ANISO_SHELL COMPLETE****'

    contains

        ! From an ANISOU pdb file, calculated the eigenvalues and eigenvectors of each atoms
        ! anisotropic displacement matrix
        subroutine read_from_aniso_pdb(funit, eigenvals, eigenvecs)
            integer, intent(in)         :: funit
            real,         intent(inout) :: eigenvals(3), eigenvecs(3,3)
            character(len=100)          :: junk1, junk2
            character(*), parameter     :: aniso_fmt='(a28,6i7,a10)'
            real                        :: aniso(3,3)
            integer                     :: aniso_in_pdb(3,3), ifoo

            aniso = 0.
            eigenvals = 0.
            eigenvecs = 0.
            aniso_in_pdb = 0
            read (funit, *)  ! Skip the 'ATOM' lines.
            read (funit, aniso_fmt) junk1, aniso_in_pdb(1,1), aniso_in_pdb(2,2), aniso_in_pdb(3,3),&
                &aniso_in_pdb(1,2), aniso_in_pdb(1,3), aniso_in_pdb(2,3), junk2
            aniso = aniso_in_pdb/10000.
            aniso(2,1) = aniso(1,2)
            aniso(3,1) = aniso(1,3)
            aniso(3,2) = aniso(2,3)
            call jacobi(aniso, 3, 3, eigenvals, eigenvecs, ifoo)
            call eigsrt(eigenvals, eigenvecs, 3, 3) ! Largest to smallest like simulated eigenvalues
            eigenvals = sqrt(eigenvals) ! ANISOU format has squared aniso valuess
        end subroutine read_from_aniso_pdb

        subroutine write_results(cc, funit, size, eval_result, eval_truth, eval_err, evecs_result, &
                &evecs_truth, evecs_err, xyz_result, xyz_truth)
            integer, intent(in)     :: cc, funit, size
            real, intent(in)        :: eval_result(3), eval_truth(3), eval_err(3), evecs_result(3,3), &
                                        &evecs_truth(3,3), evecs_err(3), xyz_result(3), xyz_truth(3)
            
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
            write(funit,601,advance='no') evecs_err(3)                        ! EVEC_ERR_3
            ! Input displacement in x,y,z directions
            write(funit,601,advance='no') xyz_result(1),            CSV_DELIM ! XRESULT
            write(funit,601,advance='no') xyz_result(2),            CSV_DELIM ! YRESULT
            write(funit,601,advance='no') xyz_result(3),            CSV_DELIM ! ZRESULT
            ! Output displacement in x,y,z directions
            write(funit,601,advance='no') xyz_truth(1),             CSV_DELIM ! XTRUTH
            write(funit,601,advance='no') xyz_truth(2),             CSV_DELIM ! YTRUTH
            write(funit,602)              xyz_truth(3)                        ! ZTRUTH
        end subroutine write_results

        ! Finds the eigenvectors by rotating the x,y,z basis by angles a,b,c about the x,y,z 
        ! axes, respectively.  Outputs a 3x3 matrix contatining the normalized eigenvectors.
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

        ! Returns the angle (in degrees) between two 3D vectors vec1 and vec2
        pure function angular_diff(vec1, vec2) result(theta)
            real, intent(in) :: vec1(3), vec2(3)
            real             :: theta
            theta = abs(dot_product(vec1,vec2))
            theta = theta / norm_2(vec1) / norm_2(vec2)
            ! Floating point arithmetic errors can cause values slightly greater than 1 to be passed to acos()
            ! so in that case replace with 1.
            theta = acos(min(1.,theta)) * 180. / PI ! Output in degrees
        end function angular_diff

        ! Returns the atomic displacements in the x,y,z directions as an array
        ! Given an input matrix of eigenvectors and array of eigenvals.
        function get_xyzdisp(eigenvecs, eigenvals) result(xyzdisp)
            real, intent(in)    :: eigenvecs(3,3), eigenvals(3)
            real                :: aniso(3,3), eigenvecs_inv(3,3), xyzdisp(3)
            real(kind=8)        :: a(3,3), b(3,3)
            integer             :: i, errflg, errflg1

            ! Build aniso matrix in principal basis
            aniso = 0.  
            aniso(1,1) = eigenvals(1)
            aniso(2,2) = eigenvals(2)
            aniso(3,3) = eigenvals(3)
            eigenvecs_inv = 0.
            ! Transform to x,y,z basis
            call matinv(eigenvecs, eigenvecs_inv, 3, errflg)
            aniso = matmul(matmul(eigenvecs, aniso), eigenvecs_inv)
            if (errflg /= 0) then
                print *, cc, eigenvecs(:,1), eigenvecs(:,2), eigenvecs(:,3)
                print *, aniso
            end if
            do i=1,3
                xyzdisp(i) = aniso(i,i)
            end do
        end function get_xyzdisp

end program simple_test_calc_aniso_shell