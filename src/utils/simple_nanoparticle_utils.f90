! Ideal Pt lattice:
! [100] = 3.9242A
! [110] = 2.77
! [111] = 6.79
module simple_nanoparticle_utils
include 'simple_lib.f08'
use simple_qr_solve
use simple_atoms, only: atoms
use simple_defs_atoms
implicit none

public :: fit_lattice, run_coord_number_analysis, strain_analysis, atoms_mask, dock_nanosPDB
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG = .false.
integer, parameter :: NSTRAIN_COMPS = 7

contains

    ! Identify the bound for defining the neighbourhood in
    ! fit_lattice and strain_analysis routines below
    function find_rMax( element ) result( rMax )
        character(len=2), intent(in) :: element
        character(len=2) :: el_ucase
        character(len=8) :: crystal_system
        real :: a_0, rMax
        el_ucase = uppercase(trim(adjustl(element)))
        call get_lattice_params(el_ucase, crystal_system, a_0)
        select case(trim(adjustl(crystal_system)))
            case('rocksalt')
                rMax = a_0 * ((1. / 2. + 1. / sqrt(2.)) / 2.)
            case('bcc')
                rMax = a_0 * ((1. + sqrt(3.) / 2.) / 2.)
            case DEFAULT ! FCC by default
                rMax = a_0 * ((1. + 1. / sqrt(2.)) / 2.)
        end select
        if( DEBUG ) write(logfhandle,*) 'rMax identified as ', rMax
    end function find_rMax

    ! ATTENTION: input coords of model have to be in ANGSTROMS
    subroutine fit_lattice( element, model, a )
        character(len=2), intent(in)    :: element
        real,             intent(inout) :: model(:,:)
        real,             intent(inout) :: a(3) ! fitted lattice parameter
        real(kind=8) , allocatable :: A_matrix(:,:), x2(:),x2_ref(:)
        real(kind=8) :: p0(3,size(model,dim=2))
        real(kind=8) :: b(3), uvw(3,3)
        real         :: cMin(3), cMax(3), cMid(3), centerAtomCoords(3), origin0(3), origin(3)
        real         :: dist, dist_temp, avgNNRR, avgD, subK, k
        real         :: rMaxsq, rMax, rMin, angleUV, angleUW, angleVW
        real         :: u0(3), v0(3), w0(3), u(3), v(3), w(3), uN(3), vN(3), wN(3), xyzbeta(4,3)
        integer      :: natoms, iatom, centerAtom,i
        logical      :: areNearest(size(model,dim=2))
        ! sanity check
        if( size(model,dim=1) /=3 ) then
            write(logfhandle,*) 'Nonconforming input coordinates! fit_lattice'
            stop
        endif
        natoms     = size  (model,dim=2)
        cMin       = minval(model,dim=2)
        cMax       = maxval(model,dim=2)
        cMid       = (cMax-cMin)/2.+cMin
        if( DEBUG ) write(logfhandle,*) 'cMid ', cMid
        dist_temp  = HUGE(dist_temp)
        centerAtom = 0
        do iatom = 1, natoms
            dist = euclid(model(:,iatom), cMid)
            if( dist < dist_temp )then
                if( DEBUG ) write(logfhandle,*) 'iatom ', iatom
                if( DEBUG ) write(logfhandle,*) 'dist ', dist
                centerAtom  = iatom
                dist_temp = dist
            endif
        enddo
        if( DEBUG ) write(logfhandle,*) 'dist_temp ', dist_temp
        if( DEBUG ) write(logfhandle,*) 'centerAtom ', centerAtom
        centerAtomCoords = model(:,centerAtom)
        if( DEBUG ) write(logfhandle,*) 'centerAtomCoords ', centerAtomCoords
        rMax       = find_rMax(element)
        rMin       = TINY
        rMaxsq     = rMax * rMax
        areNearest = .false.
        avgNNRR    = 0.
        do iatom = 1 ,natoms
            dist = euclid(model(:,iatom),cMid(:))
            if( dist <= rMax .and. dist > rMin )then ! > rMin not to count the atom itself
                areNearest(iatom) = .true.
                avgNNRR = avgNNRR + dist
            endif
        enddo
        avgNNRR = avgNNRR / real(count(areNearest))
        avgD    = sqrt(avgNNRR**2/2.)
        if( DEBUG ) write(logfhandle,*) 'avgNNRR ', avgNNRR, 'avgD ', avgD
        p0      = model ! copy of the original model
        origin0 = centerAtomCoords
        subK    = 0.5
        k       = avgD / subK
        u0      = [k,0.,0.]
        v0      = [0.,k,0.]
        w0      = [0.,0.,k]
        u0      = u0 * subK
        v0      = v0 * subK
        w0      = w0 * subK
        ! keep copies
        u       = u0
        v       = v0
        w       = w0
        origin  = origin0
        allocate( A_matrix(natoms,4), source=1.0_dp )
        if( DEBUG ) write(logfhandle,*) 'Before Loop: '
        if( DEBUG ) write(logfhandle,*) 'u v w ', u, v, w
        if( DEBUG ) write(logfhandle,*) 'origin ', origin
        uvw(1,:) = u
        uvw(2,:) = v
        uvw(3,:) = w
        do iatom = 1, natoms
            b(:) = (p0(:,iatom) - origin(:)) / subK
            call qr(uvw, b, x2)
            x2 = x2 * subK
            ! A_matrix(iatom,1) = 1. for initialization
            A_matrix(iatom,2) = nint(x2(1))
            A_matrix(iatom,3) = nint(x2(2))
            A_matrix(iatom,4) = nint(x2(3))
        enddo
        ! Refine lattice
        do i = 1, size(p0,1) ! 3
            call qr(A_matrix,p0(i,:),x2_ref)
            xyzbeta(:,i) = x2_ref
        enddo
        origin  = xyzbeta(1,:)
        u       = xyzbeta(2,:)
        v       = xyzbeta(3,:)
        w       = xyzbeta(4,:)
        ! normalize vectors
        uN      = u / norm2(u)
        vN      = v / norm2(v)
        wN      = w / norm2(w)
        ! find the angles between vectors
        angleUV = acos( dot_product(uN, vN) ) * 180. / PI
        angleUW = acos( dot_product(uN, wN) ) * 180. / PI
        angleVW = acos( dot_product(vN, wN) ) * 180. / PI
        if( DEBUG ) write(logfhandle,*) 'origin ', origin
        write(logfhandle,*) 'vector uN', uN
        write(logfhandle,*) 'vector vN', vN
        write(logfhandle,*) 'vector wN', wN
        write(logfhandle,*) '|u|', norm2(u) ! it should be in the order of 1.9
        write(logfhandle,*) '|v|', norm2(v)
        write(logfhandle,*) '|w|', norm2(w)
        write(logfhandle,*) 'angle UV', angleUV
        write(logfhandle,*) 'angle UW', angleUW
        write(logfhandle,*) 'angle VW', angleVW
        write(logfhandle,*) 'NP FCC a: ', 2. * norm2(u),2. * norm2(v),2. * norm2(w)
        ! Return calculated fitted lattice parameter
        a(1) = 2. * norm2(u)
        a(2) = 2. * norm2(v)
        a(3) = 2. * norm2(w)

        contains

            ! QR is a function with a simple interface which
            ! solves a linear system A*x = b in the least squares sense.
            ! It makes uses of the modules simple_r8lib and simple_qr_solve
            subroutine qr( a, b, x2 )
                real(kind=8) :: a(:,:)
                real(kind=8) :: b(:)
                real(kind=8), allocatable :: x2(:)
                integer :: m, n
                m = size(a,1)
                n = size(a,2)
                if( allocated(x2) ) deallocate(x2)
                allocate ( x2(n) )
                call qr_solve ( m, n, a, b, x2 )
            end subroutine qr

    end subroutine fit_lattice

    ! COORDINATION NUMBER ANALYSIS

    ! This function calculates the coordination number for each atom
    ! ATTENTION: input coords of model have to be in ANGSTROMS.
    subroutine run_coord_number_analysis( element, model, a, coord_nums_std, coord_nums_gen )
        character(len=2),  intent(in)    :: element
        real, allocatable, intent(in)    :: model(:,:)
        real,              intent(in)    :: a(3) ! lattice parameters
        integer,           intent(inout) :: coord_nums_std(size(model,2))
        real,              intent(inout) :: coord_nums_gen(size(model,2))
        character(len=2) :: el_ucase
        character(len=8) :: crystal_system
        integer :: filnum, io_stat, natoms, iatom, jatom, cnt, cn_max(size(model,2))
        real    :: dist, d, a0, foo
        ! Identify the bound for defining the neighbourhood
        a0 = sum(a)/real(size(a)) ! geometric mean of fitted lattice parameters
        el_ucase = uppercase(trim(adjustl(element)))
        call get_lattice_params(el_ucase, crystal_system, foo)
        select case(trim(adjustl(crystal_system)))
            case('rocksalt')
                d = a0 * ((1. / 2. + 1. / sqrt(2.)) / 2.)
            case('bcc')
                d = a0 * ((1. + sqrt(3.) / 2.) / 2.)
            case DEFAULT ! FCC by default
                d = a0 * ((1. + 1. / sqrt(2.)) / 2.)
        end select
        ! init
        natoms         = size(model,2)
        coord_nums_std = 0
        coord_nums_gen = 0
        cn_max         = 0
        ! cn_std
        do iatom = 1, natoms
            cnt = 0 ! reset counter, # neighbours
            do jatom = 1, natoms
                if( iatom /= jatom )then
                    dist = euclid(model(:,iatom), model(:,jatom))
                    if( dist < d ) cnt = cnt + 1
                endif
            enddo
            coord_nums_std(iatom) = cnt
        enddo
        ! cn_gen
        do iatom = 1, natoms
            cnt    = 0  ! sum of coordination numbers
            do jatom = 1, natoms
                if( iatom /= jatom )then
                    dist = euclid(model(:,iatom), model(:,jatom))
                    if( dist < d )then ! if they are neighbours
                        cnt = cnt + coord_nums_std(jatom)
                        if( coord_nums_std(jatom) > cn_max(iatom) )then
                            cn_max(iatom) = coord_nums_std(jatom)
                        endif
                    endif
                endif
            enddo
            if( cn_max(iatom) > 0 )then
                coord_nums_gen(iatom) = real(cnt) / real(cn_max(iatom))
            else
                coord_nums_gen(iatom) = 0.
            endif
        enddo
    end subroutine run_coord_number_analysis

    ! ATTENTION: input coords of model have to be in ANGSTROMS
    subroutine strain_analysis( element, model, a, strain_array )
        character(len=2), intent(in)    :: element
        real,             intent(in)    :: model(:,:)
        real,             intent(inout) :: a(3) ! fitted lattice parameter
        real,             intent(inout) :: strain_array(:,:)
        real,    parameter   :: H     = 2.      ! for weighted KDE differentiation
        integer, parameter   :: NITER = 30
        integer, allocatable :: modelCstart(:)
        real,    allocatable :: abc(:,:), modelFromABC(:,:), modelC(:,:), modelCS(:,:), modelU(:,:), modelU_spherical(:,:)
        real,    allocatable :: eXX(:), eYY(:), eZZ(:), eXY(:), eYZ(:), eXZ(:), list_displacement(:,:)
        real,    allocatable :: list_eXX(:,:), list_eYY(:,:), list_eZZ(:,:), list_eXY(:,:), list_eYZ(:,:), list_eXZ(:,:)
        real,    allocatable :: list_eRR(:,:), list_Ux(:,:), list_Uy(:,:), list_Uz(:,:)
        real,    allocatable :: r_CS(:), r_C(:), theta_CS(:), theta_C(:), phi_CS(:), phi_C(:), eRR(:)
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
        real         :: K1, K2, K3, K4, K5, K6, dR, dR_final, atm_a
        integer      :: natoms,iatom,centerAtom,i,ii,filnum,n,io_stat
        logical      :: areNearest(size(model,dim=2))
        logical      :: plusx_surface, minusx_surface, plusy_surface, minusy_surface, plusz_surface, minusz_surface
        type(atoms)  :: Exx_strain,Eyy_strain,Ezz_strain,Exy_strain,Eyz_strain,Exz_strain, Err_strain
        type(atoms)  :: Ux_atoms, Uy_atoms, Uz_atoms
        real(kind=8) :: p0(3, size(model,dim=2))
        character(len=2) :: el_ucase
        character(len=4) :: atom_name
        character(len=8) :: crystal_system
        write(logfhandle, '(A)') '>>> STRAIN ANALYSIS'
        ! sanity check
        if( size(model,dim=1 ) /= 3 ) THROW_HARD('Wrong input coordinates! strain_analysis')
        natoms = size(model, dim=2)
        if( size(strain_array,dim=1) /= natoms )        THROW_HARD('dim=1 of strain_array not conforming with model! strain_analysis')
        if( size(strain_array,dim=2) /= NSTRAIN_COMPS ) THROW_HARD('dim=2 of strain_array not conforming with NSTRAIN_COMPS! strain_analysis')
        ! naming convention (simple_atoms)
        atom_name = ' '//trim(element)//' '
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
        rMax       = find_rMax(element)
        rMin       = TINY
        areNearest = .false.
        avgNNRR    = 0.
        do iatom = 1, natoms
            dist = euclid(model(:,iatom), cMid(:))
            if( dist <= rMax .and. dist > rMin )then ! > rMin not to count the atom itself
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
        ! expected lattice parameters as uvw matrix
        el_ucase = uppercase(trim(adjustl(element)))
        call get_lattice_params(el_ucase, crystal_system, atm_a)
        uvwN0      = 0.
        uvwN0(1,1) = atm_a/2.
        uvwN0(2,2) = atm_a/2.
        uvwN0(3,3) = atm_a/2.
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
            plusx_surface   = .true.
            minusx_surface  = .true.
            plusy_surface   = .true.
            minusy_surface  = .true.
            plusz_surface   = .true.
            minusz_surface  = .true.
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
            call Exx_strain%set_name(i,atom_name)
            call Exx_strain%set_element(i,element)
            call Exx_strain%set_coord(i,list_eXX(i,1:3))
            call Exx_strain%set_beta(i,list_eXX(i,4))
            call Exx_strain%set_resnum(i,i)
            ! Eyy strain
            call Eyy_strain%set_name(i,atom_name)
            call Eyy_strain%set_element(i,element)
            call Eyy_strain%set_coord(i,list_eYY(i,1:3))
            call Eyy_strain%set_beta(i,list_eYY(i,4))
            call Eyy_strain%set_resnum(i,i)
            ! Ezz strain
            call Ezz_strain%set_name(i,atom_name)
            call Ezz_strain%set_element(i,element)
            call Ezz_strain%set_coord(i,list_eZZ(i,1:3))
            call Ezz_strain%set_beta(i,list_eZZ(i,4))
            call Ezz_strain%set_resnum(i,i)
            ! Exy strain
            call Exy_strain%set_name(i,atom_name)
            call Exy_strain%set_element(i,element)
            call Exy_strain%set_coord(i,list_eXY(i,1:3))
            call Exy_strain%set_beta(i,list_eXY(i,4))
            call Exy_strain%set_resnum(i,i)
            ! Eyz strain
            call Eyz_strain%set_name(i,atom_name)
            call Eyz_strain%set_element(i,element)
            call Eyz_strain%set_coord(i,list_eYZ(i,1:3))
            call Eyz_strain%set_beta(i,list_eYZ(i,4))
            call Eyz_strain%set_resnum(i,i)
            ! Exz strain
            call Exz_strain%set_name(i,atom_name)
            call Exz_strain%set_element(i,element)
            call Exz_strain%set_coord(i,list_eXZ(i,1:3))
            call Exz_strain%set_beta(i,list_eXZ(i,4))
            call Exz_strain%set_resnum(i,i)
            ! fill in strain array
            strain_array(i,1) = list_eXX(i,4)
            strain_array(i,2) = list_eYY(i,4)
            strain_array(i,3) = list_eZZ(i,4)
            strain_array(i,4) = list_eXY(i,4)
            strain_array(i,5) = list_eYZ(i,4)
            strain_array(i,6) = list_eXZ(i,4)
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
                    Strain_local  = eRR(i) * 100.
                    list_eRR(i,:) = [x_local,y_local,z_local,Strain_local]
                endif
            enddo
            call Err_strain%set_name(i,atom_name)
            call Err_strain%set_element(i,element)
            call Err_strain%set_coord(i,list_eRR(i,1:3))
            call Err_strain%set_beta(i,list_eRR(i,4))
            call Err_strain%set_resnum(i,i)
            ! fill in strain array
            strain_array(i,7) = list_eRR(i,4)
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
            call Ux_atoms%set_name(i,atom_name)
            call Ux_atoms%set_element(i,element)
            call Ux_atoms%set_coord(i,list_Ux(i,1:3))
            call Ux_atoms%set_beta(i,list_Ux(i,4))
            call Ux_atoms%set_resnum(i,i)
            ! Uy
            call Uy_atoms%set_name(i,atom_name)
            call Uy_atoms%set_element(i,element)
            call Uy_atoms%set_coord(i,list_Uy(i,1:3))
            call Uy_atoms%set_beta(i,list_Uy(i,4))
            call Uy_atoms%set_resnum(i,i)
            ! Uz
            call Uz_atoms%set_name(i,atom_name)
            call Uz_atoms%set_element(i,element)
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

    ! NANOPARTICLE PDB UTILS

    subroutine dock_nanosPDB( pdbfile1, pdbfile2, pdbfile1_out, pdbfile2_out )
        character(len=*), intent(inout) :: pdbfile1, pdbfile2         ! containing detected atomic positions
        character(len=*), intent(inout) :: pdbfile1_out, pdbfile2_out ! docked
        real, allocatable    :: model1(:,:), model2(:,:)              ! coords of the 2 models
        real, allocatable    :: model1_coupled(:,:),model2_coupled(:,:)
        logical, allocatable :: mask(:)
        type(atoms)          :: atom1, atom2, atom1_coupled, atom2_coupled
        character(len=2)     :: element
        character(len=4)     :: atom_name
        character(len=20)    :: atom1_cut_pdb, atom2_cut_pdb
        real    :: a1(3), a2(3) ! parameters of the fitted lattices
        real    :: geom_center1(3), geom_center2(3), origin(3), t
        real    :: radius1, radius2, radius, tmp_rad, d_max, ha, avg_params
        real    :: U(3,3), r(3), lrms ! kabsch
        integer :: iatom, n1, n2, i, nremoved
        element   = 'Pt'
        atom_name = 'Pt  '
        ! Find radii of the 2 nanos
        call atom1%new(pdbfile1)
        n1 = atom1%get_n() ! number of atoms in the model
        call atom2%new(pdbfile2)
        n2 = atom2%get_n() ! number of atoms in the model
        tmp_rad = 0.
        call read_pdb2matrix(pdbfile1, model1)
        call read_pdb2matrix(pdbfile2, model2)
        allocate(mask(n1), source = .true.)
        do iatom = 1, n1
            d_max = pixels_dist(model1(:,iatom), model1(:,:), 'max', mask)
            if( d_max > tmp_rad ) tmp_rad = d_max
        enddo
        radius1 = tmp_rad
        tmp_rad = 0.
        deallocate(mask)
        allocate( mask(n2), source=.true. )
        do iatom = 1, n2
            d_max = pixels_dist(model2(:,iatom), model2(:,:), 'max', mask)
            if( d_max > tmp_rad ) tmp_rad = d_max
        enddo
        radius2 = tmp_rad
        deallocate(mask)
        radius        = min(radius1, radius2)
        atom1_cut_pdb = 'Atom1_Cut.pdb'
        atom2_cut_pdb = 'Atom2_Cut.pdb'
        ! Remove atoms outside bigger radius to simplify couples identification
        call atoms_mask(pdbfile1, radius, atom1_cut_pdb, nremoved)
        call atoms_mask(pdbfile2, radius, atom2_cut_pdb, nremoved)
        call read_pdb2matrix(atom1_cut_pdb, model1)
        call read_pdb2matrix(atom2_cut_pdb, model2)
        call find_couples(model1, model2, element, model1_coupled, model2_coupled)
        ! Use Kabsch to register
        call kabsch(model1_coupled, model2_coupled, U, r, lrms)
        ! Write identified transformation
        write(logfhandle,*) 'Identified Rotation Matrix'
        call vis_mat(U)
        write(logfhandle,*) 'Identified Translation vector ', r
        ! Apply transformations
        call atom2%translate(-r)
        call atom2%rotate(transpose(U))
        ! Generate Output
        call atom1%writePDB(pdbfile1_out)
        call atom2%writePDB(pdbfile2_out)
        ! Delete useless files
        call del_file('Atom1_Cut.pdb')
        call del_file('Atom2_Cut.pdb')
        ! kill
        call atom1%kill
        call atom2%kill
        call atom1_coupled%kill
        call atom2_coupled%kill
    end subroutine dock_nanosPDB

    ! My implementation of the Matlab Kabsch algorithm
    ! source: https://au.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm
    ! Kabsch is for registering two pairs of 3D coords such that the RMSD is minimised
    subroutine kabsch(points_P, points_Q, U, r, lrms)
        real(sp), intent(inout) :: points_P(:,:), points_Q(:,:) ! set of points to register
        integer, parameter      :: D = 3  ! space dimension
        real(sp), intent(inout) :: U(D,D) ! rotation matrix
        real(sp), intent(inout) :: r(D)   ! translation
        real(sp), intent(inout) :: lrms   ! RMSD of registered points
        real(sp), allocatable   :: m(:), Pdm(:,:), Diff(:,:), Diff_divided(:,:)
        integer  :: i, j, N
        real(sp) :: centroid_P(D), centroid_Q(D), C(D,D), w(D),v(D,D), Id(D,D)
        if(size(points_P,dim=1) /= D  .or. size(points_Q,dim=1) /= D) then
            write(logfhandle,*) 'WORNG input dims', size(points_P,dim=1), 'and ', size(points_Q,dim=1), 'but D ', D
            stop
        elseif(size(points_P, dim = 2) /= size(points_Q, dim =2)) then
            write(logfhandle,*) 'Have to have same nb of points!'
            stop
        endif
        N = size(points_P, dim = 2)
        allocate(m(N), source = 1./real(N))
        ! translation
        centroid_P = sum(points_P(:,:), dim = 2)/real(N)
        centroid_Q = sum(points_Q(:,:), dim = 2)/real(N)
        do i = 1, N
            points_P(:,i) = points_P(:,i) - centroid_P(:)
        enddo
        do i = 1, N
            points_Q(:,i) = points_Q(:,i) - centroid_Q(:)
        enddo
        ! computation of the covariance matrix
        allocate(Pdm(D,N), source = 0.)
         do i = 1, N
             Pdm(:,i) = m(i) * points_P(:,i)
         enddo
         C = matmul(Pdm, transpose(points_Q)) ! covariance matrix of the coords
         ! svd
         call svdcmp(C,w,v)
         ! identity matrix
         Id(:,:) = 0.
         do i = 1, D
            Id(i,i) = 1.
         enddo
         if(det(matmul(C,transpose(v)),D,-1) < 0.) Id(2,2) = -1.
         U = matmul(matmul(v,Id),transpose(C))
         r = centroid_Q - matmul(U,centroid_P)
         allocate(Diff(D,N), source = 0.)
         Diff = matmul(U,points_P)-points_Q
         lrms = 0.
         do i = 1, N
             lrms = lrms + dot_product(m(i)*Diff(:,i),Diff(:,i))
         enddo
         lrms = sqrt(lrms)

         contains

             ! For calculating the determinant of a matrix
             ! source: https://rosettacode.org/wiki/Determinant_and_permanent#Fortran
             recursive function det(a,n,permanent) result(accumulation)
                 ! setting permanent to 1 computes the permanent.
                 ! setting permanent to -1 computes the determinant.
                 integer, intent(in) :: n, permanent
                 real,    intent(in) :: a(n,n)
                 real    :: b(n-1,n-1),accumulation
                 integer :: i, sgn
                 if (n .eq. 1) then
                     accumulation = a(1,1)
                 else
                     accumulation = 0
                     sgn = 1
                     do i=1, n
                         b(:, :(i-1)) = a(2:, :i-1)
                         b(:, i:) = a(2:, i+1:)
                         accumulation = accumulation + sgn * a(1, i) * det(b, n-1, permanent)
                         sgn = sgn * permanent
                     enddo
                 endif
             end function det

    end subroutine kabsch

    ! This subroutine takes in input a pdbfile and saves its
    ! coordinates in a matrix
    subroutine read_pdb2matrix( pdbfile, matrix )
        character(len=*),  intent(in)    :: pdbfile
        real, allocatable, intent(inout) :: matrix(:,:)
        integer     :: n, i
        type(atoms) :: a
        ! check that the file extension is .pdb
        if(fname2ext(pdbfile) .ne. 'pdb') THROW_HARD('Wrong extension input file! Should be pdb')
        call a%new(pdbfile)
        n = a%get_n()
        if( allocated(matrix) ) deallocate(matrix) ! for overwriting
        allocate(matrix(3,n), source = 0.) ! 3D points
        do i = 1, n
            matrix(:,i) = a%get_coord(i)
        enddo
        call a%kill
    end subroutine read_pdb2matrix

    subroutine find_couples( points_P, points_Q, element, P, Q )
        real,             intent(inout) :: points_P(:,:), points_Q(:,:)
        character(len=2), intent(inout) :: element
        real,allocatable, intent(inout) :: P(:,:), Q(:,:) ! just the couples of points
        real    :: theoretical_radius ! for threshold selection
        integer :: n   ! max number of couples
        integer :: cnt ! actual number of couples
        integer :: i, location(1), cnt_P, cnt_Q, z
        real    :: d
        character(len=2)     :: el_ucase
        logical, allocatable :: mask(:)
        real,    allocatable :: points_P_out(:,:), points_Q_out(:,:)
        real,    parameter   :: ABSURD = -10000.
        ! if(size(points_P) == size(points_Q))then
        !     allocate(P(size(points_P, dim=1),size(points_P, dim=2)),source = points_P)
        !     allocate(Q(size(points_Q, dim=1),size(points_Q, dim=2)),source = points_Q)
        !     return ! they are already coupled
        ! endif
        el_ucase = upperCase(trim(adjustl(element)))
        call get_element_Z_and_radius(el_ucase, z, theoretical_radius)
        if( z == 0 ) THROW_HARD('Unknown element: '//el_ucase)
        if( size(points_P, dim=2) <= size(points_Q, dim=2) )then
            n = size(points_P, dim=2)
            allocate(mask(n), source=.true.)
            allocate(points_P_out(3,n), points_Q_out(3,n), source = ABSURD) ! init to absurd value
            cnt = 0 ! counts the number of identified couples
            do i = 1, size(points_Q, dim=2)
                d = pixels_dist(points_Q(:,i),points_P(:,:),'min',mask,location,keep_zero=.true.)
                if( d <= 2.*theoretical_radius )then ! has a couple
                    cnt = cnt + 1
                    points_P_out(:,cnt) = points_P(:,location(1))
                    points_Q_out(:,cnt) = points_Q(:,i)
                    mask(location(1))   = .false. ! not to consider the same atom more than once
                endif
            enddo
            allocate(P(3,cnt), Q(3,cnt), source = 0.) ! correct size
            cnt_P = 0
            cnt_Q = 0
            do i = 1, n
                if( abs(points_P_out(1,i)-ABSURD) > TINY )then
                    cnt_P = cnt_P + 1
                    P(:3,cnt_P) = points_P_out(:3,i)
                endif
                if( abs(points_Q_out(1,i)-ABSURD) > TINY )then
                    cnt_Q = cnt_Q + 1
                    Q(:3,cnt_Q) = points_Q_out(:3,i)
                endif
            enddo
        else ! size(points_P, dim=2) > size(points_Q, dim=2)
            n = size(points_Q, dim=2)
            allocate(mask(n), source = .true.)
            allocate(points_P_out(3,n), points_Q_out(3,n), source = ABSURD) ! init to absurd value
            cnt = 0 ! counts the number of identified couples
            do i = 1, size(points_P, dim=2)
                d = pixels_dist(points_P(:,i),points_Q(:,:),'min',mask,location,keep_zero=.true.)
                if( d <= 2.*theoretical_radius )then ! has couple
                    cnt                 = cnt + 1
                    points_Q_out(:,cnt) = points_Q(:,location(1))
                    points_P_out(:,cnt) = points_P(:,i)
                    mask(location(1))   = .false. ! not to consider the same atom more than once
                endif
            enddo
            allocate(P(3,cnt), Q(3,cnt), source=0.)
            cnt_P = 0
            cnt_Q = 0
            do i = 1, n
                if( abs(points_P_out(1,i)-ABSURD) > TINY )then
                    cnt_P       = cnt_P + 1
                    P(:3,cnt_P) = points_P_out(:3,i)
                endif
                if( abs(points_Q_out(1,i)-ABSURD) > TINY )then
                    cnt_Q       = cnt_Q + 1
                    Q(:3,cnt_Q) = points_Q_out(:3,i)
                endif
            enddo
        endif
        deallocate(mask, points_P_out, points_Q_out)
    end subroutine find_couples

    ! This subrouine takes a pdb file input, removes all atoms beyond
    ! a given radius, outputs a new pdb file with the coordinates
    ! removed and reports how many atoms were removed
    subroutine atoms_mask( pdb_file_in, max_rad, pdb_file_out, nremoved )
        character(len=*), intent(in)    :: pdb_file_in
        character(len=*), intent(inout) :: pdb_file_out
        real,             intent(in)    :: max_rad
        integer,          intent(inout) :: nremoved
        real,    allocatable :: points_in(:,:)
        logical, allocatable :: mask(:)
        type(atoms) :: atoms_in, atoms_out
        real        :: m(3), dist
        integer     :: n, i, cnt
        character(len=STDLEN)         :: fbody
        character(len=:), allocatable :: ext
        call atoms_in%new(pdb_file_in)
        n = atoms_in%get_n()
        allocate(points_in(3,n), source = 0.)
        allocate(mask(n), source = .false.)
        ! center of mass of the input atoms
        m = 0.
        do i = 1,n
            points_in(:,i) = atoms_in%get_coord(i)
            m = m + points_in(:,i)
        enddo
        m = m/real(n)
        ! distance to the center of gravity
        do i = 1, n
            dist = euclid(m,points_in(:,i))
            if(dist <= max_rad) then
                mask(i) = .true.
            endif
        enddo
        ! save in the output pdb file
        call atoms_out%new(count(mask),dummy=.true.)
        cnt = 0
        do i = 1, n
            if(mask(i)) then
                cnt = cnt + 1
                call atoms_out%set_coord(cnt, points_in(:,i))
                ! set element and name one by one in case of biatomic nanos
                call atoms_out%set_element(cnt, atoms_in%get_element(i))
                call atoms_out%set_name(cnt, atoms_in%get_name(i))
            endif
        enddo
        ! remove the extension if present
        ext   = fname2ext(trim(pdb_file_out))
        fbody = get_fbody(trim(pdb_file_out), ext)
        call atoms_out%writePDB(fbody)
        nremoved = n - cnt
        ! kill
        deallocate(mask, points_in)
        call atoms_in%kill
        call atoms_out%kill
    end subroutine atoms_mask

end module simple_nanoparticle_utils
