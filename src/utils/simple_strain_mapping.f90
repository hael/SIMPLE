! Fit the Pt NP model and reconstruction to a lattice
! Find displacements from fit lattice
! Find strain from displacements
! Ideal Pt lattice:
! [100] = 3.9242 A
! [110] = 2.77
! [111] = 6.79
module simple_strain_mapping
include 'simple_lib.f08'
use simple_atoms, only: atoms
implicit none

public :: strain_analysis
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG = .false.
integer, parameter :: NSTRAIN_COMPS = 7

contains

    ! ATTENTION: input coords of model have to be in ANGSTROMS
    subroutine strain_analysis( model, a, strain_array )
        use simple_strings
        real,           intent(in)    :: model(:,:)
        real,           intent(inout) :: a(3) ! fitted lattice parameter
        real, optional, intent(inout) :: strain_array(:,:)
        real,    parameter     :: H     = 2.   ! for weighted KDE differentiation
        real,    parameter     :: Pt_a  = 3.92 ! ideal
        integer, parameter     :: NITER = 30
        integer, allocatable   :: modelCstart(:)
        real,    allocatable   :: abc(:,:), modelFromABC(:,:), modelC(:,:), modelCS(:,:), modelU(:,:), modelU_spherical(:,:)
        real,    allocatable   :: eXX(:), eYY(:), eZZ(:), eXY(:), eYZ(:), eXZ(:), list_displacement(:,:)
        real,    allocatable   :: list_eXX(:,:), list_eYY(:,:), list_eZZ(:,:), list_eXY(:,:), list_eYZ(:,:), list_eXZ(:,:)
        real,    allocatable   :: list_eRR(:,:), list_Ux(:,:), list_Uy(:,:), list_Uz(:,:)
        real,    allocatable   :: r_CS(:), r_C(:), theta_CS(:), theta_C(:), phi_CS(:), phi_C(:), eRR(:)
        real(kind=8) :: p0(3, size(model,dim=2))
        real         :: cMin(3), cMax(3), cMid(3), centerAtomCoords(3),origin0(3), origin(3)
        real         :: dist, dist_temp, avgNNRR, avgD, subK, k, rMax, rMin
        real         :: u0(3), v0(3), w0(3), u(3), v(3), w(3)
        real         :: dx, dy, dz, x_local, y_local, z_local, x2_local, y2_local, z2_local, Ux_local, Uy_local, Uz_local
        real         :: uvwN0(3,3), coord(3), ref(3), Strain_local, r_local
        real         :: plusx_local, plusy_local, plusz_local, minusx_local, minusy_local, minusz_local
        real         :: plusR_local, minusR_local, plust_local, minust_local, plusp_local, minusp_local
        real         :: Ux_plusx_local, Ux_minusx_local, Ux_plusy_local, Ux_minusy_local, Ux_plusz_local, Ux_minusz_local
        real         :: Uy_plusx_local, Uy_minusx_local, Uy_plusy_local, Uy_minusy_local, Uy_plusz_local, Uy_minusz_local
        real         :: Uz_plusx_local, Uz_minusx_local, Uz_plusy_local, Uz_minusy_local, Uz_plusz_local, Uz_minusz_local
        real         :: UR_plusR_local, UR_minusR_local, Ut_plusR_local, Ut_minusR_local, Up_plusR_local, Up_minusR_local
        real         :: K1, K2, K3, K4, K5, K6, dR, dR_final
        integer      :: natoms,iatom,centerAtom,i,ii,filnum,n,io_stat
        logical      :: areNearest(size(model,dim=2)), strain_arr_present
        logical      :: plusx_surface, minusx_surface, plusy_surface, minusy_surface, plusz_surface, minusz_surface
        type(atoms)  :: Exx_strain,Eyy_strain,Ezz_strain,Exy_strain,Eyz_strain,Exz_strain, Err_strain
        type(atoms)  :: Ux_atoms, Uy_atoms, Uz_atoms
        character(len=4), parameter :: ATOM_NAME = ' Pt '
        character(len=2), parameter :: ELEMENT   = 'Pt'
        write(logfhandle, '(A)') '>>> STRAIN ANALYSIS'
        ! sanity check
        if( size(model,dim=1 ) /= 3 ) THROW_HARD('Wrong input coordinates! strain_analysis')
        natoms = size(model, dim=2)
        strain_arr_present = present(strain_array)
        if( strain_arr_present )then
            if( size(strain_array,dim=1) /= natoms )        THROW_HARD('dim=1 of strain_array not conforming with model! strain_analysis')
            if( size(strain_array,dim=2) /= NSTRAIN_COMPS ) THROW_HARD('dim=2 of strain_array not conforming with NSTRAIN_COMPS! strain_analysis')
        endif
        ! supercell size
        cMin       = minval(model, dim=2)
        cMax       = maxval(model, dim=2)
        cMid       = (cMax - cMin) / 2. + cMin
        dist_temp  = HUGE(dist_temp)
        centerAtom = 0
        ! find the center atom
        do iatom=1,natoms
            dist = euclid(model(:,iatom), cMid)
            if( dist < dist_temp )then
                centerAtom = iatom
                dist_temp  = dist
            endif
        enddo
        centerAtomCoords = model(:,centerAtom)
        ! find the nearest neighbors of the center atom
        rMax       = 3.5
        rMin       = 0.
        areNearest = .false.
        avgNNRR    = 0.
        do iatom = 1, natoms
            dist = euclid(model(:,iatom), cMid(:))
            if( dist < rMax .and. dist > rMin )then ! > rMin not to count the atom itself
                areNearest(iatom) = .true.
                avgNNRR = avgNNRR + dist
            endif
        enddo
        avgNNRR = avgNNRR / real(count(areNearest))
        avgD    = sqrt(avgNNRR**2/2.)
        p0      = model ! copy of the original model
        origin0 = centerAtomCoords
        subK    = 0.5
        ! find the value to plug into u0,v0 and w0
        k  = avgD / subK
        u0 = [k, 0., 0.]
        v0 = [0., k, 0.]
        w0 = [0., 0., k]
        u0 = u0 * subK
        v0 = v0 * subK
        w0 = w0 * subK
        ! keep copies
        origin = origin0
        u = [a(1)/2., 0., 0.]
        v = [0., a(2)/2., 0.]
        w = [0., 0., a(3)/2.]
        allocate( abc(3,natoms), source=0. )
        abc(1,:) = (model(1,:) - origin(1)) / (a(1)/2.)
        abc(2,:) = (model(2,:) - origin(2)) / (a(2)/2.)
        abc(3,:) = (model(3,:) - origin(3)) / (a(3)/2.)
        dx = a(1)
        dy = a(2)
        dz = a(3)
        ! expected Pt FCC lattice parameters as uvw matrix
        uvwN0      = 0.
        uvwN0(1,1) = Pt_a/2.
        uvwN0(2,2) = Pt_a/2.
        uvwN0(3,3) = Pt_a/2.
        allocate( modelFromABC(natoms,3), modelC(natoms,3), modelCS(natoms,3), modelU(natoms,3) )
        ! use perfect Pt FCC lattice and remove bad points
        modelFromABC = matmul(transpose(nint(abc)), uvwN0) ! create real positions from the fit; centered at the center atom
        do iatom = 1, natoms
            modelCS(iatom,:) = model(:,iatom) - origin(:)  ! the found positions centered at the fitted origin (in A)
        enddo
        modelC = modelFromABC ! expected positions centered at the fitted origin
        ! calculate displacements
        modelU = modelCS - modelC ! displacements of measured - expected
        allocate( list_displacement(natoms,6), source=0. )
        do i = 1, natoms
            x_local  = modelCS(i,1) + origin(1)
            y_local  = modelCS(i,2) + origin(2)
            z_local  = modelCS(i,3) + origin(3)
            Ux_local = modelU(i,1)
            Uy_local = modelU(i,2)
            Uz_local = modelU(i,3)
            list_displacement(i,1) = x_local
            list_displacement(i,2) = y_local
            list_displacement(i,3) = z_local
            list_displacement(i,4) = Ux_local
            list_displacement(i,5) = Uy_local
            list_displacement(i,6) = Uz_local
        enddo
        if( DEBUG )then
            call fopen(filnum, file='displacement.txt')
            do i = 1, natoms
                write(filnum,*) list_displacement(i,:)
            enddo
            call fclose(filnum)
        endif
        ! start of lattice
        allocate( modelCstart(1), source=0 )
        modelCstart(:) = minloc(sqrt(sum(modelC**2, dim=2)))
        write(logfhandle,*) 'Origin for strain calculation, atom #', modelCstart, 'at ', modelC(modelCstart,:)
        write(logfhandle,*) 'Origin Strain ', modelU(modelCstart,:)
        ! allocation
        allocate( eXX(natoms), eYY(natoms), eZZ(natoms), eXY(natoms), eYZ(natoms), eXZ(natoms), source=0. )
        ! main loop
        do i = 1, natoms
            dx = a(1)
            dy = a(2)
            dz = a(3)
            x_local = modelCS(i,1)
            y_local = modelCS(i,2)
            z_local = modelCS(i,3)
            plusx_local  = 0.
            minusx_local = 0.
            plusy_local  = 0.
            minusy_local = 0.
            plusz_local  = 0.
            minusz_local = 0.
            ! in x
            Ux_plusx_local  = 0.
            Ux_minusx_local = 0.
            Ux_plusy_local  = 0.
            Ux_minusy_local = 0.
            Ux_plusz_local  = 0.
            Ux_minusz_local = 0.
            ! in y
            Uy_plusx_local  = 0.
            Uy_minusx_local = 0.
            Uy_plusy_local  = 0.
            Uy_minusy_local = 0.
            Uy_plusz_local  = 0.
            Uy_minusz_local = 0.
            ! in z
            Uz_plusx_local  = 0.
            Uz_minusx_local = 0.
            Uz_plusy_local  = 0.
            Uz_minusy_local = 0.
            Uz_plusz_local  = 0.
            Uz_minusz_local = 0.
            plusx_surface  = .true.
            minusx_surface = .true.
            plusy_surface  = .true.
            minusy_surface = .true.
            plusz_surface  = .true.
            minusz_surface = .true.
            do ii = 1, natoms
                x2_local        = modelCS(ii,1)
                y2_local        = modelCS(ii,2)
                z2_local        = modelCS(ii,3)
                coord(:)        = [x2_local,y2_local,z2_local]
                ref(:)          = [x_local+dx/2.,y_local,z_local]
                call gaussian_kernel(coord,ref,h,K1)
                Ux_plusx_local  = Ux_plusx_local + modelU(ii,1) * K1
                Uy_plusx_local  = Uy_plusx_local + modelU(ii,2) * K1
                Uz_plusx_local  = Uz_plusx_local + modelU(ii,3) * K1
                plusx_local     = plusx_local + K1
                ref(:)          = [x_local-dx/2.,y_local,z_local]
                call gaussian_kernel(coord,ref,h,K2)
                Ux_minusx_local = Ux_minusx_local + modelU(ii,1) * K2
                Uy_minusx_local = Uy_minusx_local + modelU(ii,2) * K2
                Uz_minusx_local = Uz_minusx_local + modelU(ii,3) * K2
                minusx_local    = minusx_local + K2
                ref(:)          = [x_local,y_local+dy/2.,z_local]
                call gaussian_kernel(coord,ref,h,K3)
                Ux_plusy_local  = Ux_plusy_local + modelU(ii,1) * K3
                Uy_plusy_local  = Uy_plusy_local + modelU(ii,2) * K3
                Uz_plusy_local  = Uz_plusy_local + modelU(ii,3) * K3
                plusy_local     = plusy_local + K3
                ref(:)          = [x_local,y_local-dy/2.,z_local]
                call gaussian_kernel(coord,ref,h,K4)
                Ux_minusy_local = Ux_minusy_local + modelU(ii,1) * K4
                Uy_minusy_local = Uy_minusy_local + modelU(ii,2) * K4
                Uz_minusy_local = Uz_minusy_local + modelU(ii,3) * K4
                minusy_local    = minusy_local + K4
                ref(:)          = [x_local,y_local,z_local+dz/2.]
                call gaussian_kernel(coord,ref,h,K5)
                Ux_plusz_local  = Ux_plusz_local + modelU(ii,1) * K5
                Uy_plusz_local  = Uy_plusz_local + modelU(ii,2) * K5
                Uz_plusz_local  = Uz_plusz_local + modelU(ii,3) * K5
                plusz_local     = plusz_local + K5
                ref(:)          = [x_local,y_local,z_local-dz/2.]
                call gaussian_kernel(coord,ref,h,K6)
                Ux_minusz_local = Ux_minusz_local + modelU(ii,1) * K6
                Uy_minusz_local = Uy_minusz_local + modelU(ii,2) * K6
                Uz_minusz_local = Uz_minusz_local + modelU(ii,3) * K6
                minusz_local    = minusz_local + K6
                if(((modelC(ii,1) - (modelC(i,1) - a(1)))**2  + (modelC(ii,2) - modelC(i,2))**2 +&
                &   (modelC(ii,3) -  modelC(i,3))**2) < 2.**2)         minusx_surface = .false.
                if(((modelC(ii,1) - (modelC(i,1) + a(1)))**2  + (modelC(ii,2) - modelC(i,2))**2 +&
                &   (modelC(ii,3) -  modelC(i,3))**2) < 2.**2)         plusx_surface  = .false.
                if(((modelC(ii,1) -  modelC(i,1))**2          + (modelC(ii,2) - (modelC(i,2)-a(2)))**2 +&
                &   (modelC(ii,3) -  modelC(i,3))**2) < 2.**2)         minusy_surface = .false.
                if(((modelC(ii,1) -  modelC(i,1))**2          + (modelC(ii,2) - (modelC(i,2)+a(2)))**2 +&
                &   (modelC(ii,3) -  modelC(i,3))**2) < 2.**2)         plusy_surface  = .false.
                if(((modelC(ii,1) -  modelC(i,1))**2          + (modelC(ii,2) - modelC(i,2))**2 +&
                &   (modelC(ii,3) - (modelC(i,3) - a(3)))**2) < 2.**2) minusz_surface = .false.
                if(((modelC(ii,1) -  modelC(i,1))**2          + (modelC(ii,2) - modelC(i,2))**2 +&
                &   (modelC(ii,3) - (modelC(i,3) + a(3)))**2) < 2.**2) plusz_surface  = .false.
            enddo
            if( plusx_surface )then
                Ux_plusx_local = modelU(i,1) * plusx_local
                Uy_plusx_local = modelU(i,2) * plusx_local
                Uz_plusx_local = modelU(i,3) * plusx_local
                dx = dx / 2.
            endif
            if( plusy_surface )then
                Ux_plusy_local = modelU(i,1) * plusy_local
                Uy_plusy_local = modelU(i,2) * plusy_local
                Uz_plusy_local = modelU(i,3) * plusy_local
                dy = dy / 2.
            endif
            if( plusz_surface )then
                Ux_plusz_local = modelU(i,1) * plusz_local
                Uy_plusz_local = modelU(i,2) * plusz_local
                Uz_plusz_local = modelU(i,3) * plusz_local
                dz = dz / 2.
            endif
            if( minusx_surface )then
                Ux_minusx_local = modelU(i,1) * minusx_local
                Uy_minusx_local = modelU(i,2) * minusx_local
                Uz_minusx_local = modelU(i,3) * minusx_local
                dx = dx / 2.
            endif
            if( minusy_surface )then
                Ux_minusy_local = modelU(i,1) * minusy_local
                Uy_minusy_local = modelU(i,2) * minusy_local
                Uz_minusy_local = modelU(i,3) * minusy_local
                dy = dy / 2.
            endif
            if( minusz_surface )then
                Ux_minusz_local = modelU(i,1) * minusz_local
                Uy_minusz_local = modelU(i,2) * minusz_local
                Uz_minusz_local = modelU(i,3) * minusz_local
                dz = dz / 2.
            endif
            eXX(i) =  (Ux_plusx_local / plusx_local - Ux_minusx_local / minusx_local) / dx
            eYY(i) =  (Uy_plusy_local / plusy_local - Uy_minusy_local / minusy_local) / dy
            eZZ(i) =  (Uz_plusz_local / plusz_local - Uz_minusz_local / minusz_local) / dz
            eXY(i) = ((Ux_plusy_local / plusy_local - Ux_minusy_local / minusy_local) / dy +&
            &         (Uy_plusx_local / plusx_local - Uy_minusx_local / minusx_local) / dx) / 2.
            eYZ(i) = ((Uy_plusz_local / plusz_local - Uy_minusz_local / minusz_local) / dz +&
            &         (Uz_plusy_local / plusy_local - Uz_minusy_local / minusy_local) / dy) / 2.
            eXZ(i) = ((Ux_plusz_local / plusz_local - Ux_minusz_local / minusz_local) / dz +&
            &         (Uz_plusx_local / plusx_local - Uz_minusx_local / minusx_local) / dx) / 2.
        enddo
        allocate( list_eXX(natoms,4), list_eYY(natoms,4), list_eZZ(natoms,4), list_eXY(natoms,4), list_eYZ(natoms,4), list_eXZ(natoms,4) )
        ! type atoms to write pdb files
        call Exx_strain%new(natoms, dummy=.true.)
        call Eyy_strain%new(natoms, dummy=.true.)
        call Ezz_strain%new(natoms, dummy=.true.)
        call Exy_strain%new(natoms, dummy=.true.)
        call Eyz_strain%new(natoms, dummy=.true.)
        call Exz_strain%new(natoms, dummy=.true.)
        n = 0
        do i = 1,natoms
            x_local = modelCS(i,1) + origin(1)
            y_local = modelCS(i,2) + origin(2)
            z_local = modelCS(i,3) + origin(3)
            ref = [x_local, y_local, z_local]
            do ii = 1, natoms
                coord = model(:,ii)
                if( euclid(coord,ref) < 0.1 )then
                    n = n + 1
                    Strain_local  = eXX(i) * 100.
                    list_eXX(i,:) = [x_local,y_local,z_local,Strain_local]
                    Strain_local  = eYY(i) * 100.
                    list_eYY(i,:) = [x_local,y_local,z_local,Strain_local]
                    Strain_local  = eZZ(i) * 100.
                    list_eZZ(i,:) = [x_local,y_local,z_local,Strain_local]
                    Strain_local  = eXY(i) * 100.
                    list_eXY(i,:) = [x_local,y_local,z_local,Strain_local]
                    Strain_local  = eYZ(i) * 100.
                    list_eYZ(i,:) = [x_local,y_local,z_local,Strain_local]
                    Strain_local  = eXZ(i) * 100.
                    list_eXZ(i,:) = [x_local,y_local,z_local,Strain_local]
                endif
            enddo
            ! Exx strain
            call Exx_strain%set_name(i,ATOM_NAME)
            call Exx_strain%set_element(i,ELEMENT)
            call Exx_strain%set_coord(i,list_eXX(i,1:3))
            call Exx_strain%set_beta(i,list_eXX(i,4))
            call Exx_strain%set_resnum(i,i)
            ! Eyy strain
            call Eyy_strain%set_name(i,ATOM_NAME)
            call Eyy_strain%set_element(i,ELEMENT)
            call Eyy_strain%set_coord(i,list_eYY(i,1:3))
            call Eyy_strain%set_beta(i,list_eYY(i,4))
            call Eyy_strain%set_resnum(i,i)
            ! Ezz strain
            call Ezz_strain%set_name(i,ATOM_NAME)
            call Ezz_strain%set_element(i,ELEMENT)
            call Ezz_strain%set_coord(i,list_eZZ(i,1:3))
            call Ezz_strain%set_beta(i,list_eZZ(i,4))
            call Ezz_strain%set_resnum(i,i)
            ! Exy strain
            call Exy_strain%set_name(i,ATOM_NAME)
            call Exy_strain%set_element(i,ELEMENT)
            call Exy_strain%set_coord(i,list_eXY(i,1:3))
            call Exy_strain%set_beta(i,list_eXY(i,4))
            call Exy_strain%set_resnum(i,i)
            ! Eyz strain
            call Eyz_strain%set_name(i,ATOM_NAME)
            call Eyz_strain%set_element(i,ELEMENT)
            call Eyz_strain%set_coord(i,list_eYZ(i,1:3))
            call Eyz_strain%set_beta(i,list_eYZ(i,4))
            call Eyz_strain%set_resnum(i,i)
            ! Exz strain
            call Exz_strain%set_name(i,ATOM_NAME)
            call Exz_strain%set_element(i,ELEMENT)
            call Exz_strain%set_coord(i,list_eXZ(i,1:3))
            call Exz_strain%set_beta(i,list_eXZ(i,4))
            call Exz_strain%set_resnum(i,i)
            if( strain_arr_present )then
                ! fill in strain array
                strain_array(i,1) = list_eXX(i,4)
                strain_array(i,2) = list_eYY(i,4)
                strain_array(i,3) = list_eZZ(i,4)
                strain_array(i,4) = list_eXY(i,4)
                strain_array(i,5) = list_eYZ(i,4)
                strain_array(i,6) = list_eXZ(i,4)
            endif
        enddo
        ! output PDB files
        call Exx_strain%writepdb('Exx_strain')
        call Eyy_strain%writepdb('Eyy_strain')
        call Ezz_strain%writepdb('Ezz_strain')
        call Exy_strain%writepdb('Exy_strain')
        call Eyz_strain%writepdb('Eyz_strain')
        call Exz_strain%writepdb('Exz_strain')
        ! kill
        call Exx_strain%kill
        call Eyy_strain%kill
        call Ezz_strain%kill
        call Exy_strain%kill
        call Eyz_strain%kill
        call Exz_strain%kill
        ! calculate displacements
        allocate( r_CS(natoms), r_C(natoms), theta_CS(natoms), theta_C(natoms), phi_CS(natoms), phi_C(natoms), source = 0.)
        r_CS     = ((modelCS(:,1)**2 + modelCS(:,2)**2 + modelCS(:,3)**2)**0.5) + 10.**(-10.)
        r_C      = ((modelC(:,1)**2  + modelC(:,2)**2  + modelC(:,3)**2)**0.5)  + 10.**(-10.)
        theta_CS = acos(modelCS(:,3) / r_CS)
        theta_C  = acos(modelC(:,3)  / r_C)
        phi_CS   = atan(modelCS(:,1), modelCS(:,2)) ! in radians
        phi_C    = atan(modelC(:,1),  modelC(:,2))
        allocate(modelU_spherical(natoms,3), source = 0.)
        modelU_spherical(:,1) = r_CS     - r_C
        modelU_spherical(:,2) = theta_CS - theta_C
        modelU_spherical(:,3) = phi_CS   - phi_C
        ! this uses the change in u and distance
        dR = 2.
        allocate( eRR(natoms), source=0. )
        do i = 1, natoms
            x_local         = modelCS(i,1)
            y_local         = modelCS(i,2)
            z_local         = modelCS(i,3)
            r_local         = r_CS(i)
            plusR_local     = 0.
            minusR_local    = 0.
            plust_local     = 0.
            minust_local    = 0.
            plusp_local     = 0.
            minusp_local    = 0.
            UR_plusR_local  = 0.
            UR_minusR_local = 0.
            do ii = 1, natoms
                x2_local        = modelCS(ii,1)
                y2_local        = modelCS(ii,2)
                z2_local        = modelCS(ii,3)
                coord(:)        = [x2_local,y2_local,z2_local]
                ref(:)          = [x_local*(1.+dR/(2.*r_local)),y_local*(1.+dR/(2.*r_local)),z_local*(1.+dR/(2.*r_local))]
                call gaussian_kernel(coord,ref,h,K1)
                UR_plusR_local  = UR_plusR_local + modelU_spherical(ii,1) * K1
                Ut_plusR_local  = Ut_plusR_local + modelU_spherical(ii,2) * K1
                Up_plusR_local  = Up_plusR_local + modelU_spherical(ii,3) * K1
                plusR_local     = plusR_local + K1
                ref(:)          = [x_local*(1.-dR/(2.*r_local)),y_local*(1.-dR/(2.*r_local)),z_local*(1.-dR/(2.*r_local))]
                call gaussian_kernel(coord,ref,h,K2)
                UR_minusR_local = UR_minusR_local + modelU_spherical(ii,1) * K2
                Ut_minusR_local = Ut_minusR_local + modelU_spherical(ii,2) * K2
                Up_minusR_local = Up_minusR_local + modelU_spherical(ii,3) * K2
                minusR_local    = minusR_local + K2
            enddo
            dR_final = dR
            eRR(i)   = ((UR_plusR_local / plusR_local - UR_minusR_local / minusR_local) / dR_final)
        enddo
        allocate(list_eRR(natoms,4), source = 0.)
        call Err_strain%new(natoms, dummy=.true.)
        n = 0
        do i = 1,natoms
            x_local = modelCS(i,1) + origin(1)
            y_local = modelCS(i,2) + origin(2)
            z_local = modelCS(i,3) + origin(3)
            ref = [x_local, y_local, z_local]
            do ii = 1, natoms
                coord = model(:,ii)
                if(euclid(coord,ref) <0.1) then
                    n = n + 1
                    Strain_local   = eRR(i)*100.
                    list_eRR(i,:) = [x_local,y_local,z_local,Strain_local]
                endif
            enddo
            call Err_strain%set_name(i,ATOM_NAME)
            call Err_strain%set_element(i,ELEMENT)
            call Err_strain%set_coord(i,list_eRR(i,1:3))
            call Err_strain%set_beta(i,list_eRR(i,4))
            call Err_strain%set_resnum(i,i)
            ! fill in strain array
            if( strain_arr_present ) strain_array(i,7) = list_eRR(i,4)
        enddo
        ! output PDB file
        call Err_strain%writepdb('radial_strain')
        call Err_strain%kill
        allocate(list_Ux(natoms,4), list_Uy(natoms,4), list_Uz(natoms,4), source = 0.)
        call Ux_atoms%new(natoms, dummy=.true.)
        call Uy_atoms%new(natoms, dummy=.true.)
        call Uz_atoms%new(natoms, dummy=.true.)
        n = 0
        do i = 1, natoms
            x_local = modelCS(i,1) + origin(1)
            y_local = modelCS(i,2) + origin(2)
            z_local = modelCS(i,3) + origin(3)
            ref = [x_local, y_local, z_local]
            do ii = 1, natoms
                coord = model(:,ii)
                if(euclid(coord,ref) <0.1) then
                    n = n + 1
                    list_Ux(i,:) = [x_local,y_local,z_local,modelU(i,1)]
                    list_Uy(i,:) = [x_local,y_local,z_local,modelU(i,2)]
                    list_Uz(i,:) = [x_local,y_local,z_local,modelU(i,3)]
                endif
            enddo
            ! Ux
            call Ux_atoms%set_name(i,ATOM_NAME)
            call Ux_atoms%set_element(i,ELEMENT)
            call Ux_atoms%set_coord(i,list_Ux(i,1:3))
            call Ux_atoms%set_beta(i,list_Ux(i,4))
            call Ux_atoms%set_resnum(i,i)
            ! Uy
            call Uy_atoms%set_name(i,ATOM_NAME)
            call Uy_atoms%set_element(i,ELEMENT)
            call Uy_atoms%set_coord(i,list_Uy(i,1:3))
            call Uy_atoms%set_beta(i,list_Uy(i,4))
            call Uy_atoms%set_resnum(i,i)
            ! Uz
            call Uz_atoms%set_name(i,ATOM_NAME)
            call Uz_atoms%set_element(i,ELEMENT)
            call Uz_atoms%set_coord(i,list_Uz(i,1:3))
            call Uz_atoms%set_beta(i,list_Uz(i,4))
            call Uz_atoms%set_resnum(i,i)
        enddo
        if( DEBUG )then
            ! PDB files
            call Ux_atoms%writepdb('Ux')
            call Uy_atoms%writepdb('Uy')
            call Uz_atoms%writepdb('Uz')
           ! CSV files
           call fopen(filnum, file='Ux.csv', iostat=io_stat)
           write (filnum,*) 'ux'
           do i = 1, natoms
               write (filnum,'(A)', advance='yes') trim(real2str(list_Ux(i,4)))
           end do
           call fclose(filnum)
           call fopen(filnum, file='Uy.csv', iostat=io_stat)
           write (filnum,*) 'uy'
           do i = 1, natoms
               write (filnum,'(A)', advance='yes') trim(real2str(list_Uy(i,4)))
           end do
           call fclose(filnum)
           call fopen(filnum, file='Uz.csv', iostat=io_stat)
           write (filnum,*) 'uz'
           do i = 1, natoms
               write (filnum,'(A)', advance='yes') trim(real2str(list_Uz(i,4)))
           end do
           call fclose(filnum)
       endif
       ! kill
       call Ux_atoms%kill
       call Uy_atoms%kill
       call Uz_atoms%kill
       write(logfhandle, '(A)') '>>> STRAIN ANALYSIS, COMPLETED'

       contains

           subroutine gaussian_kernel( coord, ref, h, K )
               real, intent(in)  :: coord(3), ref(3), h
               real, intent(out) :: K
               real :: u
               u = euclid(coord, ref) / h
               K = exp( -0.5 * (u**2))
           end subroutine gaussian_kernel

end subroutine strain_analysis

end module simple_strain_mapping
