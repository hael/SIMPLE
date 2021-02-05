! Fit the Pt NP model and reconstruction to a lattice
! Find displacements from fit lattice
! Find strain from displacements
! Ideal Pt lattice:
! [100] = 3.9242A
! [110] = 2.77
! [111] = 6.79
module simple_lattice_fitting
  use simple_image,    only : image
  use simple_binimage, only : binimage
  use simple_fileio
  use simple_math
  use simple_atoms,    only : atoms
  implicit none

  public :: fit_lattice, test_lattice_fit
  private

contains

subroutine fit_lattice(model)
  ! use simple_lapacksgels, only : sgels
  real, intent(inout) :: model(:,:)
  real, parameter     :: xRes0=1e6, yRes0=1e6,zRes0=1e6
  integer, parameter  :: NITER = 30
  real    :: p0(3, size(model,dim=2))
  real    :: cMin(3), cMax(3), cMid(3), centerAtomCoords(3),origin0(3), origin(3)
  real    :: distsq_temp, distsq, dist
  real    :: avgNNRR, avgD, subK, k
  real    :: rMaxsq, rMax, rMin
  real    :: u0(3),v0(3),w0(3),u(3),v(3),w(3), uvw(3,3)
  integer :: iloop,natoms,iatom, centerAtom
  logical :: areNearest(size(model,dim=2))
  ! sgels
  real, allocatable :: work(:)
  integer :: info, lwork
  character(len=1) :: trans
  ! ! sanity check
  ! if(size(model,dim=1) /=3 ) then
  !   write(logfhandle,*) 'Wrong input coordinates! fit_lattice'
  !   stop
  ! endif
  ! natoms=size(model,dim=2)
  ! cMin = minval(model,dim=2)
  ! cMax = maxval(model,dim=2)
  ! cMid = (cMax-cMin)/2.+cMin
  ! print *, 'cMin, cMax, cMid', cMin, cMax,cMid
  ! distsq_temp = HUGE(distsq_temp)
  ! centerAtom = 0
  ! do iatom=1,natoms
  !     distsq = sum((model(:,iatom)-cMid(:))**2)
  !     if(distsq < distsq_temp) then
  !         centerAtom = iatom
  !         distsq_temp = distsq
  !     endif
  ! enddo
  ! centerAtomCoords =  model(:,centerAtom)
  ! print *, 'closest atom to the center: ',centerAtomCoords
  ! rMaxsq= 3.5*3.5
  ! areNearest = .false.
  ! avgNNRR = 0.
  ! do iatom=1,natoms
  !     distsq = sum((model(:,iatom)-centerAtomCoords(:))**2)
  !     if(distsq < rMaxsq) then
  !       areNearest(iatom) = .true.
  !       print *, 'CoordNearest ', model(:,iatom)
  !       avgNNRR = avgNNRR + sqrt(distsq)
  !     endif
  ! enddo
  ! areNearest(centerAtom) = .false. ! do not count center Atom
  ! avgNNRR = avgNNRR/real(count(areNearest))
  ! print *, 'avgNNRR', avgNNRR
  ! avgD = sqrt(avgNNRR**2/2.)
  ! print *, 'avgD ', avgD
  ! p0      = model ! copy of the original model
  ! origin0 = centerAtomCoords
  ! subK    = 0.5
  ! k       = avgD/subK
  ! print *, 'k: ', k
  ! u0 = [k,0.,0.]
  ! v0 = [0.,k,0.]
  ! w0 = [0.,0.,k]
  ! print *, 'u0,v0,w0 before: ', u0,v0,w0
  ! u0 = u0*subK
  ! v0 = v0*subK
  ! w0 = w0*subK
  ! print *, 'u0,v0,w0 after: ', u0,v0,w0
  ! ! keep copies
  ! u      = u0
  ! v      = v0
  ! w      = w0
  ! origin = origin0
  ! print *, 'u,v,w before loop', u,v,w
  ! print *, 'origin before loop', origin
  ! ! translate and scale now to facilitate the call in sgels
  ! do iatom = 1,natoms
  !   p0(:,iatom) = (p0(:,iatom)-origin(:))/subK
  ! enddo
  ! ! variables to be used in sgels
  ! allocate(work(natoms), source=-1.)
  ! trans = 'T'
  ! lwork =3*3+natoms
  ! ! LOOP
  !  !do iloop =1, NITER
  !    ! compute a,b,c values by finding points p0 in basis {u,v,w} from {i,j,k}
  !     uvw(1,:) = u
  !     uvw(2,:) = v
  !     uvw(3,:) = w
  !
  !     call sgels(trans,3,3,natoms,uvw,3,p0,natoms,w,lwork,info)
  !     ! call sgels(	TRANS,M,N,NRHS,A,LDA,B,LDB,WORK,LWORK,INFO)
  !
  ! !   abc = abc*subK
  ! !   ! Refine lattice
  ! !   A      = 1.
  ! !   A(:,1) = round(abc(1,:))
  ! !   A(:,1) = round(abc(1,:))
  ! !   A(:,1) = round(abc(1,:))
  ! ! enddo
end subroutine fit_lattice

subroutine test_lattice_fit
  real, allocatable  :: model(:,:)
  integer            :: i
  character(len=100) :: fname

  fname='model_coordinates.txt'
  call read_3Dcoord(fname,model)
  call fit_lattice(model)

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
end subroutine test_lattice_fit

end module simple_lattice_fitting
