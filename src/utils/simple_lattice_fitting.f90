! Fit the Pt NP model and reconstruction to a lattice
! Find displacements from fit lattice
! Find strain from displacements
! Ideal Pt lattice:
! [100] = 3.9242A
! [110] = 2.77
! [111] = 6.79
module simple_lattice_fitting
  use simple_fileio
  use simple_math
  use simple_atoms,    only : atoms
  use simple_qr_solve
  use simple_r8lib
  use simple_np_coordination_number, only : run_coord_number_analysis
  implicit none

  public :: test_lattice_fit, find_radius_for_coord_number, fit_lattice
  private

  logical, parameter :: DEBUG = .false.

contains

  ! QR is a function with a simple interface which
  ! solves a linear system A*x = b in the least squares sense.
  ! It makes uses of the modules simple_r8lib and simple_qr_solve
  subroutine qr ( a,b,x2 )
    real ( kind = 8 ) :: a(:,:)
    real ( kind = 8 ) :: b(:)
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n
    real ( kind = 8 ), allocatable :: x2(:)
    m = size(a,1)
    n = size(a,2)
    !  Allocate space.
    if(allocated(x2)) deallocate(x2)
    allocate ( x2(1:n) )
    !  Use QR_SOLVE on the problem.
      call qr_solve ( m, n, a, b, x2 )
    return
  end subroutine qr

! ATTENTION: input coords of model have to be in ANGSTROMS.
subroutine fit_lattice(model,a)
  real, intent(inout) :: model(:,:)
  real, intent(inout) :: a(3) ! fitted lattice parameter
  integer, parameter  :: NITER = 30
  real (kind = 8) p0(3, size(model,dim=2))
  real    :: cMin(3), cMax(3), cMid(3), centerAtomCoords(3),origin0(3), origin(3)
  real    :: dist, dist_temp, avgNNRR, avgD, subK, k
  real    :: rMaxsq, rMax, rMin, angleUV, angleUW, angleVW
  real    :: u0(3),v0(3),w0(3),u(3),v(3),w(3),uN(3),vN(3),wN(3),xyzbeta(4,3)
  integer :: iloop,natoms,iatom, centerAtom,i
  logical :: areNearest(size(model,dim=2))
  ! qr solve
  real (kind = 8) , allocatable :: A_matrix(:,:), x2(:),x2_ref(:)
  real (kind = 8) b(3), uvw(3,3)
  ! sanity check
  if(size(model,dim=1) /=3 ) then
    write(logfhandle,*) 'Wrong input coordinates! fit_lattice'
    stop
  endif
  natoms = size  (model,dim=2)
  cMin   = minval(model,dim=2)
  cMax   = maxval(model,dim=2)
  cMid   = (cMax-cMin)/2.+cMin
  if(DEBUG) write(logfhandle,*) 'cMid: ', cMid
  dist_temp = HUGE(dist_temp)
  centerAtom = 0
  do iatom=1,natoms
      dist = euclid(model(:,iatom), cMid)
      if(dist < dist_temp) then
        if(DEBUG) write(logfhandle,*)  'iatom ', iatom
        if(DEBUG) write(logfhandle,*)  'dist ', dist
          centerAtom  = iatom
          dist_temp = dist
      endif
  enddo
  if(DEBUG) write(logfhandle,*) 'dist_temp', dist_temp
  if(DEBUG) write(logfhandle,*) 'centerAtom ', centerAtom
  centerAtomCoords =  model(:,centerAtom)
  if(DEBUG) write(logfhandle,*) '  centerAtomCoords ',   centerAtomCoords
  rMax = 3.5
  rMin = 0.
  rMaxsq= rMax*rMax
  areNearest = .false.
  avgNNRR = 0.
  do iatom=1,natoms
    ! dist = euclid(model(:,iatom),centerAtomCoords(:))
      dist = euclid(model(:,iatom),cMid(:)) ! ALREADY CHANGED
      if(dist < rMax .and. dist > rMin) then ! >rMin not to count the atom itself
        areNearest(iatom) = .true.
        avgNNRR = avgNNRR + dist
      endif
  enddo
  avgNNRR = avgNNRR/real(count(areNearest))
  avgD = sqrt(avgNNRR**2/2.)
  if(DEBUG) write(logfhandle,*) 'avgNNRR', avgNNRR, 'avgD ', avgD
  p0      = model ! copy of the original model
  origin0 = centerAtomCoords
  subK    = 0.5
  k       = avgD/subK
  u0 = [k,0.,0.]
  v0 = [0.,k,0.]
  w0 = [0.,0.,k]
  u0 = u0*subK
  v0 = v0*subK
  w0 = w0*subK
  ! keep copies
  u      = u0
  v      = v0
  w      = w0
  origin = origin0
  allocate(A_matrix(natoms,4), source = 1.0_dp)
  if(DEBUG) write(logfhandle,*)  'Before Loop: '
  if(DEBUG) write(logfhandle,*) 'u v w ', u, v, w
  if(DEBUG) write(logfhandle,*) 'origin ', origin
  ! LOOP
   do iloop =1,NITER+1 !indexes start from 0 in Python
      uvw(1,:) = u
      uvw(2,:) = v
      uvw(3,:) = w
      do iatom = 1,natoms
        b(:) = (p0(:,iatom)-origin(:))/subK
        call qr(uvw,b,x2)
        x2 = x2 * subK
        !A_matrix(iatom,1)= 1. for initiation
        A_matrix(iatom,2) = nint(x2(1))
        A_matrix(iatom,3) = nint(x2(2))
        A_matrix(iatom,4) = nint(x2(3))
      enddo
      ! Refine lattice
      do i = 1, size(p0,1) ! 3
          call qr(A_matrix,p0(i,:),x2_ref)
          xyzbeta(:,i) = x2_ref
      enddo
      origin = xyzbeta(1,:)
      u = xyzbeta(2,:)
      v = xyzbeta(3,:)
      w = xyzbeta(4,:)
   enddo
   ! normalize vectors
   uN = u/norm2(u)
   vN = v/norm2(v)
   wN = w/norm2(w)
   ! find the angles between vectors
   angleUV =  acos( dot_product(uN, vN) )*180./PI
   angleUW =  acos( dot_product(uN, wN) )*180./PI
   angleVW =  acos( dot_product(vN, wN) )*180./PI
   if(DEBUG) write(logfhandle,*) 'origin ', origin
   write(logfhandle,*) 'vector uN', uN
   write(logfhandle,*) 'vector vN', vN
   write(logfhandle,*) 'vector wN', wN
   write(logfhandle,*) '|u|', norm2(u) ! it should be in the order of 1.9
   write(logfhandle,*) '|v|', norm2(v)
   write(logfhandle,*) '|w|', norm2(w)
   write(logfhandle,*) 'angle UV', angleUV
   write(logfhandle,*) 'angle UW', angleUW
   write(logfhandle,*) 'angle VW', angleVW
   write(logfhandle,*) 'PtNP FCC a: ', 2.*norm2(u),2.*norm2(v),2.*norm2(w)
   ! Return calculated fitted lattice parameter
   a(1) = 2.*norm2(u)
   a(2) = 2.*norm2(v)
   a(3) = 2.*norm2(w)
end subroutine fit_lattice

subroutine test_lattice_fit()
  real, allocatable  :: model(:,:)
  integer            :: i
  character(len=100) :: fname
  real :: a(3), d

  fname='model_coordinates.txt'
  call read_3Dcoord(fname,model)
  call fit_lattice(model,a)
  call find_radius_for_coord_number(a,d)
  call run_coord_number_analysis(model, d)

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

! CN of an atom is calculated as the number of neighboring atoms within specific distance d
! 	d=(d_1+d_2)/2, where
!   d_1=(a_0)/√2
! 	d_2=a_0	where
!   a_0 is fitted lattice constant of corresponding fcc nanocrystal.
!   a_0 is geometric mean of three components of the lattice parameter
subroutine find_radius_for_coord_number(a,d)
  real, intent(in)    :: a(3)  ! fitted lattice parameter
  real, intent(inout) :: d     ! radius for coordination number calculation
  real :: a0, d1, d2
  a0 = sum(a)/real(size(a)) ! geometric mean
  d1 = a0/sqrt(2.)
  d2 = a0
  d = (d1+d2)/2.
  if(DEBUG) write(logfhandle,*) 'Identified radius for coord numbe calculation ', d
end subroutine find_radius_for_coord_number
end module simple_lattice_fitting
