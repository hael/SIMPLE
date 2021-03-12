! Fit the Pt NP model and reconstruction to a lattice
! Find displacements from fit lattice
! Find strain from displacements
! Ideal Pt lattice:
! [100] = 3.9242A
! [110] = 2.77
! [111] = 6.79
module simple_strain_mapping
  use simple_fileio
  use simple_math
  use simple_atoms,    only : atoms
  use simple_lattice_fitting
  implicit none

  public :: test_strain_analysis, strain_analysis
  private

contains

! ATTENTION: input coords of model have to be in ANGSTROMS.
subroutine strain_analysis(model,a)
  real, intent(inout)  :: model(:,:)
  real, intent(inout)  :: a(3)        ! fitted lattice parameter
  real, parameter      :: H = 2.      ! for Weighted KDE differentiation
  real, parameter      :: Pt_a = 3.92 ! ideal
  integer, parameter   :: NITER = 30
  integer, allocatable :: modelCstart(:)
  real, allocatable :: abc(:,:), modelFromABC(:,:), modelC(:,:), modelCS(:,:), modelU(:,:), modelU_spherical(:,:)
  real, allocatable :: eXX(:),eYY(:),eZZ(:),eXY(:),eYZ(:),eXZ(:),list_displacement(:,:)
  real, allocatable :: list_eXX(:,:),list_eYY(:,:),list_eZZ(:,:),list_eXY(:,:),list_eYZ(:,:),list_eXZ(:,:)
  real, allocatable :: list_eRR(:,:), list_Ux(:,:), list_Uy(:,:), list_Uz(:,:)
  real, allocatable :: r_CS(:),r_C(:), theta_CS(:),theta_C(:),phi_CS(:),phi_C(:), eRR(:)
  real (kind = 8) p0(3, size(model,dim=2))
  real    :: cMin(3), cMax(3), cMid(3), centerAtomCoords(3),origin0(3), origin(3)
  real    :: dist, dist_temp, avgNNRR, avgD, subK, k, rMax, rMin
  real    :: u0(3),v0(3),w0(3),u(3),v(3),w(3)
  real    :: dx, dy, dz, x_local, y_local, z_local, x2_local, y2_local, z2_local, Ux_local, Uy_local, Uz_local
  real    :: uvwN0(3,3),coord(3), ref(3), Strain_local, r_local
  real    :: plusx_local, plusy_local, plusz_local, minusx_local, minusy_local, minusz_local
  real    :: plusR_local, minusR_local, plust_local, minust_local, plusp_local, minusp_local
  real    :: Ux_plusx_local, Ux_minusx_local, Ux_plusy_local, Ux_minusy_local, Ux_plusz_local, Ux_minusz_local
  real    :: Uy_plusx_local, Uy_minusx_local, Uy_plusy_local, Uy_minusy_local, Uy_plusz_local, Uy_minusz_local
  real    :: Uz_plusx_local, Uz_minusx_local, Uz_plusy_local, Uz_minusy_local, Uz_plusz_local, Uz_minusz_local
  real    :: UR_plusR_local, UR_minusR_local, Ut_plusR_local, Ut_minusR_local, Up_plusR_local, Up_minusR_local
  real    :: K1, K2, K3, K4, K5, K6, dR, dR_final
  integer :: natoms,iatom,centerAtom,i,ii,filnum,n
  logical :: areNearest(size(model,dim=2))
  logical :: plusx_surface, minusx_surface, plusy_surface, minusy_surface, plusz_surface, minusz_surface
  ! sanity check
  if(size(model,dim=1) /=3 ) then
    write(logfhandle,*) 'Wrong input coordinates! fit_lattice'
    stop
  endif
  natoms=size(model,dim=2)
  !Supercell size
  cMin = minval(model,dim=2)
  cMax = maxval(model,dim=2)
  cMid = (cMax-cMin)/2.+cMin
  dist_temp = HUGE(dist_temp)
  centerAtom = 0
  !Find the center atom
  do iatom=1,natoms
      dist = euclid(model(:,iatom), cMid)
      if(dist < dist_temp) then
          centerAtom  = iatom
          dist_temp = dist
      endif
  enddo
  centerAtomCoords =  model(:,centerAtom)
  ! Find the nearest neighbors of the center atom
  rMax = 3.5
  rMin = 0.
  areNearest = .false.
  avgNNRR = 0.
  do iatom=1,natoms
    dist = euclid(model(:,iatom),cMid(:))
      if(dist < rMax .and. dist > rMin) then ! >rMin not to count the atom itself
        areNearest(iatom) = .true.
        avgNNRR = avgNNRR + dist
      endif
  enddo
  avgNNRR = avgNNRR/real(count(areNearest))
  avgD = sqrt(avgNNRR**2/2.)
  p0      = model ! copy of the original model
  origin0 = centerAtomCoords
  subK    = 0.5
  !Find the value to plug into u0,v0 and w0
  k       = avgD/subK
  u0 = [k,0.,0.]
  v0 = [0.,k,0.]
  w0 = [0.,0.,k]
  u0 = u0*subK
  v0 = v0*subK
  w0 = w0*subK

  ! keep copies
  origin = origin0
  u = [a(1)/2.,0.,0.]
  v = [0.,a(2)/2.,0.]
  w = [0.,0.,a(3)/2.]

  allocate(abc(3,natoms), source = 0.)
  abc(1,:) = (model(1,:)-origin(1))/(a(1)/2.)
  abc(2,:) = (model(2,:)-origin(2))/(a(2)/2.)
  abc(3,:) = (model(3,:)-origin(3))/(a(3)/2.)

  dx = a(1)
  dy = a(2)
  dz = a(3)

  ! Expected Pt FCC lattice parameters as uvw matrix
  uvwN0 = 0.
  uvwN0(1,1) = Pt_a/2.
  uvwN0(2,2) = Pt_a/2.
  uvwN0(3,3) = Pt_a/2.
  allocate(modelFromABC(natoms,3), modelC(natoms,3), modelCS(natoms,3), modelU(natoms,3))
  ! Use perfect Pt FCC lattice and remove bad points
  modelFromABC = matmul(transpose(nint(abc)),uvwN0) ! create real positions from the fit; centered at the center atom
  do iatom = 1, natoms
    modelCS(iatom,:) = model(:,iatom) - origin(:) ! the found positions centered at the fitted origin (in Ang)
  enddo
  modelC = modelFromABC     ! expected positions centered at the fitted origin

  ! Calculate displacements
  modelU = modelCS - modelC ! displacements of measured - expected
  allocate(list_displacement(natoms,6), source=0.)
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

  call fopen(filnum, file='displacement.txt')
  do i = 1, natoms
    write(filnum,*) list_displacement(i,:)
  enddo
  call fclose(filnum)
  ! Start of lattice
  allocate(modelCstart(1), source = 0)
  modelCstart(:) = minloc(sqrt(sum(modelC**2, dim=2)))
  write(logfhandle,*) 'Origin for strain calculation, atom #', modelCstart, 'at ', modelC(modelCstart,:)
  write(logfhandle,*) 'Origin Strain ', modelU(modelCstart, :)
  !allocation
  allocate(eXX(natoms), eYY(natoms), eZZ(natoms), eXY(natoms), eYZ(natoms), eXZ(natoms), source = 0.)
  ! main loop
  do i = 1, natoms
    dx=a(1)
    dy=a(2)
    dz=a(3)

    x_local = modelCS(i,1)
    y_local = modelCS(i,2)
    z_local = modelCS(i,3)

    plusx_local  = 0.
    minusx_local = 0.
    plusy_local  = 0.
    minusy_local = 0.
    plusz_local  = 0.
    minusz_local = 0.

    !In x
    Ux_plusx_local  = 0.
    Ux_minusx_local = 0.
    Ux_plusy_local  = 0.
    Ux_minusy_local = 0.
    Ux_plusz_local  = 0.
    Ux_minusz_local = 0.
    !In y
    Uy_plusx_local  = 0.
    Uy_minusx_local = 0.
    Uy_plusy_local  = 0.
    Uy_minusy_local = 0.
    Uy_plusz_local  = 0.
    Uy_minusz_local = 0.
    !In z
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
      x2_local = modelCS(ii,1)
      y2_local = modelCS(ii,2)
      z2_local = modelCS(ii,3)

      coord(:) = [x2_local,y2_local,z2_local]
      ref(:)   = [x_local+dx/2.,y_local,z_local]
      call gaussian_kernel(coord,ref,h,K1)
      Ux_plusx_local = Ux_plusx_local + modelU(ii,1)*K1
      Uy_plusx_local = Uy_plusx_local + modelU(ii,2)*K1
      Uz_plusx_local = Uz_plusx_local + modelU(ii,3)*K1
      plusx_local    = plusx_local + K1

      ref(:)   = [x_local-dx/2.,y_local,z_local]
      call gaussian_kernel(coord,ref,h,K2)
      Ux_minusx_local = Ux_minusx_local + modelU(ii,1)*K2
      Uy_minusx_local = Uy_minusx_local + modelU(ii,2)*K2
      Uz_minusx_local = Uz_minusx_local + modelU(ii,3)*K2
      minusx_local    = minusx_local + K2

      ref(:)   = [x_local,y_local+dy/2.,z_local]
      call gaussian_kernel(coord,ref,h,K3)
      Ux_plusy_local = Ux_plusy_local + modelU(ii,1)*K3
      Uy_plusy_local = Uy_plusy_local + modelU(ii,2)*K3
      Uz_plusy_local = Uz_plusy_local + modelU(ii,3)*K3
      plusy_local    = plusy_local + K3

      ref(:)   = [x_local,y_local-dy/2.,z_local]
      call gaussian_kernel(coord,ref,h,K4)
      Ux_minusy_local = Ux_minusy_local + modelU(ii,1)*K4
      Uy_minusy_local = Uy_minusy_local + modelU(ii,2)*K4
      Uz_minusy_local = Uz_minusy_local + modelU(ii,3)*K4
      minusy_local    = minusy_local + K4

      ref(:)   = [x_local,y_local,z_local+dz/2.]
      call gaussian_kernel(coord,ref,h,K5)
      Ux_plusz_local = Ux_plusz_local + modelU(ii,1)*K5
      Uy_plusz_local = Uy_plusz_local + modelU(ii,2)*K5
      Uz_plusz_local = Uz_plusz_local + modelU(ii,3)*K5
      plusz_local    = plusz_local + K5

      ref(:)   = [x_local,y_local,z_local-dz/2.]
      call gaussian_kernel(coord,ref,h,K6)
      Ux_minusz_local = Ux_minusz_local + modelU(ii,1)*K6
      Uy_minusz_local = Uy_minusz_local + modelU(ii,2)*K6
      Uz_minusz_local = Uz_minusz_local + modelU(ii,3)*K6
      minusz_local    = minusz_local + K6

      if(((modelC(ii,1)-(modelC(i,1)-a(1)))**2 + (modelC(ii,2)-modelC(i,2))**2 + (modelC(ii,3)-modelC(i,3))**2) < 2.**2) minusx_surface = .false.
      if(((modelC(ii,1)-(modelC(i,1)+a(1)))**2 + (modelC(ii,2)-modelC(i,2))**2 + (modelC(ii,3)-modelC(i,3))**2) < 2.**2) plusx_surface  = .false.
      if(((modelC(ii,1)-modelC(i,1))**2 + (modelC(ii,2)-(modelC(i,2)-a(2)))**2 + (modelC(ii,3)-modelC(i,3))**2) < 2.**2) minusy_surface = .false.
      if(((modelC(ii,1)-modelC(i,1))**2 + (modelC(ii,2)-(modelC(i,2)+a(2)))**2 + (modelC(ii,3)-modelC(i,3))**2) < 2.**2) plusy_surface  = .false.
      if(((modelC(ii,1)-modelC(i,1))**2 + (modelC(ii,2)-modelC(i,2))**2 + (modelC(ii,3)-(modelC(i,3)-a(3)))**2) < 2.**2) minusz_surface = .false.
      if(((modelC(ii,1)-modelC(i,1))**2 + (modelC(ii,2)-modelC(i,2))**2 + (modelC(ii,3)-(modelC(i,3)+a(3)))**2) < 2.**2) plusz_surface  = .false.
    enddo
    if(plusx_surface) then
      Ux_plusx_local = modelU(i,1)*plusx_local
      Uy_plusx_local = modelU(i,2)*plusx_local
      Uz_plusx_local = modelU(i,3)*plusx_local
      dx = dx/2.
    endif
    if(plusy_surface) then
      Ux_plusy_local = modelU(i,1)*plusy_local
      Uy_plusy_local = modelU(i,2)*plusy_local
      Uz_plusy_local = modelU(i,3)*plusy_local
      dy = dy/2.
    endif
    if(plusz_surface) then
      Ux_plusz_local = modelU(i,1)*plusz_local
      Uy_plusz_local = modelU(i,2)*plusz_local
      Uz_plusz_local = modelU(i,3)*plusz_local
      dz = dz/2.
    endif
    if(minusx_surface) then
      Ux_minusx_local = modelU(i,1)*minusx_local
      Uy_minusx_local = modelU(i,2)*minusx_local
      Uz_minusx_local = modelU(i,3)*minusx_local
      dx = dx/2.
    endif
    if(minusy_surface) then
      Ux_minusy_local = modelU(i,1)*minusy_local
      Uy_minusy_local = modelU(i,2)*minusy_local
      Uz_minusy_local = modelU(i,3)*minusy_local
      dy = dy/2.
    endif
    if(minusz_surface) then
      Ux_minusz_local = modelU(i,1)*minusz_local
      Uy_minusz_local = modelU(i,2)*minusz_local
      Uz_minusz_local = modelU(i,3)*minusz_local
      dz = dz/2.
    endif
    eXX(i) = (Ux_plusx_local/plusx_local  - Ux_minusx_local/minusx_local)/dx
    eYY(i) = (Uy_plusy_local/plusy_local  - Uy_minusy_local/minusy_local)/dy
    eZZ(i) = (Uz_plusz_local/plusz_local  - Uz_minusz_local/minusz_local)/dz
    eXY(i) = ((Ux_plusy_local/plusy_local - Ux_minusy_local/minusy_local)/dy + (Uy_plusx_local/plusx_local - Uy_minusx_local/minusx_local)/dx)/2.
    eYZ(i) = ((Uy_plusz_local/plusz_local - Uy_minusz_local/minusz_local)/dz + (Uz_plusy_local/plusy_local - Uz_minusy_local/minusy_local)/dy)/2.
    eXZ(i) = ((Ux_plusz_local/plusz_local - Ux_minusz_local/minusz_local)/dz + (Uz_plusx_local/plusx_local - Uz_minusx_local/minusx_local)/dx)/2.
  enddo

  allocate(list_eXX(natoms,4),list_eYY(natoms,4),list_eZZ(natoms,4),list_eXY(natoms,4),list_eYZ(natoms,4),list_eXZ(natoms,4))
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
        Strain_local = eXX(i)*100.
        list_eXX(i,:)  = [x_local,y_local,z_local,Strain_local]
        Strain_local = eYY(i)*100.
        list_eYY(i,:)  = [x_local,y_local,z_local,Strain_local]
        Strain_local = eZZ(i)*100.
        list_eZZ(i,:)  = [x_local,y_local,z_local,Strain_local]
        Strain_local = eXY(i)*100.
        list_eXY(i,:)  = [x_local,y_local,z_local,Strain_local]
        Strain_local = eYZ(i)*100.
        list_eYZ(i,:)  = [x_local,y_local,z_local,Strain_local]
        Strain_local = eXZ(i)*100.
        list_eXZ(i,:)  = [x_local,y_local,z_local,Strain_local]
      endif
    enddo
  enddo
  ! Calculate displacements
  allocate(r_CS(natoms),r_C(natoms),theta_CS(natoms),theta_C(natoms),phi_CS(natoms),phi_C(natoms),source = 0.)
  r_CS     = ((modelCS(:,1)**2+modelCS(:,2)**2+modelCS(:,3)**2)**0.5)+10.**(-10.)
  r_C      = ((modelC(:,1)**2+modelC(:,2)**2+modelC(:,3)**2)**0.5)+10.**(-10.)
  theta_CS = acos(modelCS(:,3)/r_CS)
  theta_C  = acos(modelC(:,3)/r_C)
  phi_CS   = atan(modelCS(:,1),modelCS(:,2)) ! in radians
  phi_C    = atan(modelC(:,1),modelC(:,2))

  allocate(modelU_spherical(natoms,3), source = 0.)
  modelU_spherical(:,1) = r_CS-r_C
  modelU_spherical(:,2) = theta_CS-theta_C
  modelU_spherical(:,3) = phi_CS-phi_C
  ! This uses the change in u and distance
  dR = 2.
  allocate(eRR(natoms), source = 0.)
  do i = 1, natoms
    x_local     = modelCS(i,1)
    y_local     = modelCS(i,2)
    z_local     = modelCS(i,3)
    r_local     = r_CS(i)
    plusR_local = 0.
    minusR_local = 0.
    plust_local = 0.
    minust_local = 0.
    plusp_local = 0.
    minusp_local = 0.
    UR_plusR_local = 0.
    UR_minusR_local = 0.
    do ii = 1, natoms
      x2_local = modelCS(ii,1)
      y2_local = modelCS(ii,2)
      z2_local = modelCS(ii,3)

      coord(:) = [x2_local,y2_local,z2_local]
      ref(:)   = [x_local*(1.+dR/(2.*r_local)),y_local*(1.+dR/(2.*r_local)),z_local*(1.+dR/(2.*r_local))]
      call gaussian_kernel(coord,ref,h,K1)
      UR_plusR_local = UR_plusR_local + modelU_spherical(ii,1)*K1
      Ut_plusR_local = Ut_plusR_local + modelU_spherical(ii,2)*K1
      Up_plusR_local = Up_plusR_local + modelU_spherical(ii,3)*K1
      plusR_local    = plusR_local + K1

      ref(:)   = [x_local*(1.-dR/(2.*r_local)),y_local*(1.-dR/(2.*r_local)),z_local*(1.-dR/(2.*r_local))]
      call gaussian_kernel(coord,ref,h,K2)
      UR_minusR_local = UR_minusR_local + modelU_spherical(ii,1)*K2
      Ut_minusR_local = Ut_minusR_local + modelU_spherical(ii,2)*K2
      Up_minusR_local = Up_minusR_local + modelU_spherical(ii,3)*K2
      minusR_local    = minusR_local + K2
    enddo
    dR_final = dR
    eRR(i) = ((UR_plusR_local/plusR_local-UR_minusR_local/minusR_local)/dR_final)
  enddo

  allocate(list_eRR(natoms,4), source = 0.)
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
  enddo

  allocate(list_Ux(natoms,4), list_Uy(natoms,4), list_Uz(natoms,4), source = 0.)
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
  enddo

contains

  subroutine gaussian_kernel(coord, ref, h, K)
    real, intent(inout) :: coord(3), ref(3)  !input
    real, intent(in)    :: h
    real, intent(inout) :: K                 !output
    real :: u
    u = euclid(coord,ref)/h
    K = exp(-0.5*(u**2))
  end subroutine gaussian_kernel
end subroutine strain_analysis


subroutine test_strain_analysis()
  real, allocatable  :: model(:,:)
  character(len=100) :: fname
  real :: a(3)

  ! In my laptop's folder Desktop/Nanoparticles/StrainAnalysis/Code_Packages_for_Chiara
  fname='model_coordinates_strain_mapping.txt'
  call read_3Dcoord(fname,model)
  ! just for testing purposes
  a(1) = 3.990
  a(2) = 3.970
  a(3) = 3.968
  call strain_analysis(model,a)

contains

  subroutine read_3Dcoord( fname, vals )
         character(len=*),  intent(in)  :: fname    !< input filename
         real, allocatable, intent(out) :: vals(:,:) !< array of values
         integer :: nl, funit, iline,io_stat
         nl = nlines(trim(fname))
         call fopen(funit,fname,'old','unknown',io_stat)
         call fileiochk("read_table failed to open file "//trim(fname),io_stat )
         allocate( vals(3,nl), stat=alloc_stat )
         if(alloc_stat /= 0) call allocchk ('In: read_filetable; simple_fileio  ', alloc_stat)
         do iline=1,nl
             read(funit,*) vals(1,iline), vals(2,iline), vals(3,iline)
         end do
         call fclose(funit,io_stat)
         call fileiochk("read_filetable failed to close",io_stat)
     end subroutine read_3Dcoord
end subroutine test_strain_analysis
end module simple_strain_mapping
