module simple_nano_utils
  use simple_image,    only : image
  use simple_binimage, only : binimage
  use simple_rnd
  use simple_math
  use simple_defs
  use simple_atoms,        only : atoms
  use simple_nanoparticle, only :nanoparticle
  implicit none

  public :: kabsch, read_pdb2matrix, write_matrix2pdb, find_couples, atoms_mask, dock_nanosPDB, test_dockPDB
  private


  interface kabsch
      module procedure kabsch_sp, kabsch_dp
  end interface kabsch

  interface read_pdb2matrix
      module procedure read_pdb2matrix_sp, read_pdb2matrix_dp
  end interface read_pdb2matrix

contains

  subroutine test_dockPDB(pdbfile,theta,t)
    character(len=*), intent(inout) :: pdbfile
    real, intent(inout) :: theta ! rotation angle
    real, intent(inout) :: t(3)  ! translation vector
    type(atoms)         :: a
    real                :: U(3,3), r(3) ! identified rotation and translation
    character(len=20)   :: pdbfile2,pdbfile1_out,pdbfile2_out
    ! Definition of rotation matrix
    U(1,1) = cos(theta)
    U(1,2) = -sin(theta)
    U(1,3) = 0.
    U(2,1) = sin(theta)
    U(2,2) = cos(theta)
    U(2,3) = 0.
    U(3,1) = 0.
    U(3,2) = 0.
    U(3,3) = 1.
    write(logfhandle,*) 'Set Rotation Matrix to '
    call vis_mat(U)
    write(logfhandle,*) 'Set translation Vector to ', t
    call a%new(pdbfile)
    call a%rotate(U)
    call a%translate(t)
    pdbfile2 = 'Transformed'
    call a%writePDB(pdbfile2)
    pdbfile1_out = 'Pdbfile1_out'
    pdbfile2_out = 'Pdbfile2_out'
    pdbfile2 = 'Transformed.pdb'
    call dock_nanosPDB(pdbfile,pdbfile2,pdbfile1_out,pdbfile2_out)
  end subroutine test_dockPDB

  subroutine dock_nanosPDB(pdbfile1,pdbfile2,pdbfile1_out,pdbfile2_out)
    use simple_syslib
    character(len=*), intent(inout) :: pdbfile1,pdbfile2         ! containing detected atomic positions
    character(len=*), intent(inout) :: pdbfile1_out,pdbfile2_out ! docked
    real, allocatable    :: model1(:,:), model2(:,:) ! coords of the 2 models
    real, allocatable    :: model1_coupled(:,:),model2_coupled(:,:)
    logical, allocatable :: mask(:)
    type(atoms)          :: atom1,atom2,atom1_coupled,atom2_coupled
    character(len=2)     :: element
    character(len=4)     :: atom_name
    character(len=20)    :: atom1_cut_pdb,atom2_cut_pdb
    real    :: a1(3), a2(3)  ! parameters of the fitted lattice
    real    :: geom_center1(3),geom_center2(3),origin(3),t
    real    :: radius1,radius2,radius,tmp_rad,d_max,ha,avg_params
    real    :: U(3,3), r(3), lrms ! kabsch
    integer :: iatom,n1,n2,i,nremoved
    element   = 'Pt'
    atom_name = 'Pt  '
    ! Find radii of the 2 nanos
    call atom1%new(pdbfile1)
    n1 = atom1%get_n() ! number of atoms in the model
    call atom2%new(pdbfile2)
    n2 = atom2%get_n() ! number of atoms in the model
    tmp_rad = 0.
    call read_pdb2matrix(pdbfile1,model1)
    call read_pdb2matrix(pdbfile2,model2)
    allocate(mask(n1), source = .true.)
    do iatom = 1, n1
      d_max = pixels_dist(model1(:,iatom),model1(:,:),'max',mask)
      if(d_max > tmp_rad) tmp_rad = d_max
    enddo
    radius1 = tmp_rad
    tmp_rad = 0.
    deallocate(mask)
    allocate(mask(n2), source = .true.)
    do iatom = 1, n2
      d_max = pixels_dist(model2(:,iatom),model2(:,:),'max',mask)
      if(d_max > tmp_rad) tmp_rad = d_max
    enddo
    radius2 = tmp_rad
    deallocate(mask)
    radius   = min(radius1,radius2)
    atom1_cut_pdb = 'Atom1_Cut.pdb'
    atom2_cut_pdb = 'Atom2_Cut.pdb'
    ! Remove atoms outside bigger radius to simplify couples identification
    call atoms_mask(pdbfile1,radius,atom1_cut_pdb,nremoved)
    call atoms_mask(pdbfile2,radius,atom2_cut_pdb,nremoved)
    call read_pdb2matrix(atom1_cut_pdb,model1)
    call read_pdb2matrix(atom2_cut_pdb,model2)
    call find_couples(model1,model2,element,model1_coupled,model2_coupled)
    ! Use Kabsch to register
    call kabsch(model1_coupled,model2_coupled,U,r,lrms)
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
  ! Kabsch is for registering two pairs of 3D coords such that the RMSD
  ! is minimised
  subroutine kabsch_sp(points_P, points_Q, U, r, lrms)
     real(sp), intent(inout) :: points_P(:,:), points_Q(:,:) ! set of points to register
     integer, parameter      :: D = 3  !space dimension
     real(sp), intent(inout) :: U(D,D) !rotation matrix
     real(sp), intent(inout) :: r(D)   !translation
     real(sp), intent(inout) :: lrms   !RMSD of registered points
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
         points_P(:,i) = points_P(:,i)-centroid_P(:)
     enddo
     do i = 1, N
         points_Q(:,i) = points_Q(:,i)-centroid_Q(:)
     enddo
     ! computation of the covariance matrix
     allocate(Pdm(D,N), source = 0.)
     do i = 1, N
         Pdm(:,i) = m(i)*points_P(:,i)
     enddo
     C = matmul(Pdm, transpose(points_Q)) !covariance matrix of the coords
     ! svd
     call svdcmp(C,w,v)
     ! identity matrix
     Id(:,:) = 0.
     do i = 1, D
         do j = 1, D
           if(i==j)  Id(i,j) = 1.
         enddo
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
  end subroutine kabsch_sp

  ! My implementation of the Matlab Kabsch algorithm
  ! source: https://au.mathworks.com/matlabcentral/fileexchange/25746-kabsch-algorithm
  ! Kabsch is for registering two pairs of 3D coords such that the RMSD
  ! is minimised
  subroutine kabsch_dp(points_P, points_Q, U, r, lrms)
     real(dp), intent(inout) :: points_P(:,:), points_Q(:,:) ! set of points to register
     integer,  parameter     :: D = 3  !space dimension
     real(dp), intent(inout) :: U(D,D) !rotation matrix
     real(dp), intent(inout) :: r(D)   !translation
     real(dp), intent(inout) :: lrms   !RMSD of registered points
     real(dp), allocatable   :: m(:), Pdm(:,:), Diff(:,:), Diff_divided(:,:)
     integer  :: i, j, N
     real(dp) :: centroid_P(D), centroid_Q(D), C(D,D), w(D),v(D,D), Id(D,D)
     if(size(points_P,dim=1) /= D  .or. size(points_Q,dim=1) /= D) then
       write(logfhandle,*) 'WORNG input dims'
       stop
     elseif(size(points_P, dim = 2) /= size(points_Q, dim =2)) then
       write(logfhandle,*) 'Have to have same nb of points!'
       stop
     endif
     N = size(points_P, dim = 2)
     allocate(m(N), source = real(1./real(N),dp))
     ! translation
     centroid_P = sum(points_P(:,:), dim = 2)/real(N)
     centroid_Q = sum(points_Q(:,:), dim = 2)/real(N)
     do i = 1, N
         points_P(:,i) = points_P(:,i)-centroid_P(:)
     enddo
     do i = 1, N
         points_Q(:,i) = points_Q(:,i)-centroid_Q(:)
     enddo
     ! computation of the covariance matrix
     allocate(Pdm(D,N), source = 0.d0)
     do i = 1, N
         Pdm(:,i) = m(i)*points_P(:,i)
     enddo
     C = matmul(Pdm, transpose(points_Q)) !covariance matrix of the coords
     ! svd
     call svdcmp(C,w,v)
     ! identity matrix
     Id(:,:) = 0.
     do i = 1, D
         do j = 1, D
           if(i==j)  Id(i,j) = 1.
         enddo
     enddo
     if(det(real(matmul(C,transpose(v)),sp),D,-1) < 0.) Id(2,2) = -1.
     U = matmul(matmul(v,Id),transpose(C))
     r = centroid_Q - matmul(U,centroid_P)
     allocate(Diff(D,N), source = 0.d0)
     Diff = matmul(U,points_P)-points_Q
     lrms = 0.d0
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
  end subroutine kabsch_dp

  ! This subroutine takes in input a pdbfile and saves its
  ! coordinates onto matrix
  subroutine read_pdb2matrix_sp(pdbfile, matrix)
    use simple_fileio, only: fname2ext
    character(len=*),  intent(in)    :: pdbfile
    real, allocatable, intent(inout) :: matrix(:,:)
    integer     :: n, i
    type(atoms) :: a
    ! check that the file extension is .pdb
    if(fname2ext(pdbfile) .ne. 'pdb') then
        write(logfhandle, *) 'Wrong extension input file! Should be pdb'
        stop
    endif
    call a%new(pdbfile)
    n = a%get_n()
    if(allocated(matrix)) deallocate(matrix) ! for overwriting
    allocate(matrix(3,n), source = 0.) ! 3D points
    do i = 1, n
        matrix(:3,i) = a%get_coord(i)
    enddo
    call a%kill
  end subroutine read_pdb2matrix_sp

  ! This subroutine takes in input a pdbfile and saves its
  ! coordinates onto matrix
  subroutine read_pdb2matrix_dp(pdbfile, matrix)
    use simple_fileio, only: fname2ext
    character(len=*),    intent(in)    :: pdbfile
    real(dp),allocatable,intent(inout) :: matrix(:,:)
    integer     :: n, i
    type(atoms) :: a
    ! check that the file extension is .pdb
    if(fname2ext(pdbfile) .ne. 'pdb') then
        write(logfhandle, *) 'Wrong extension input file! Should be pdb'
        stop
    endif
    call a%new(pdbfile)
    n = a%get_n()
    allocate(matrix(3,n), source = 0.d0) ! 3D points
    do i = 1, n
        matrix(:3,i) = real(a%get_coord(i),dp)
    enddo
    call a%kill
  end subroutine read_pdb2matrix_dp

  ! This subroutine takes in input a matrix and saves its
  ! entries onto a pdbfile
  subroutine write_matrix2pdb(matrix, pdbfile)
    use simple_fileio, only: fname2ext
    character(len=*),  intent(in) :: pdbfile
    real,              intent(in) :: matrix(:,:)
    integer     :: n, i
    type(atoms) :: a
    if(size(matrix, dim=1) /= 3) then
        write(logfhandle,*) 'Wrong input matrix dimension!'
        stop
    endif
    n = size(matrix,dim=2)
    call a%new(n, dummy = .true.)
    do i = 1, n
        call a%set_coord(i,matrix(:3,i))
    enddo
    call a%writePDB(pdbfile)
    call a%kill
  end subroutine write_matrix2pdb

  subroutine find_couples(points_P, points_Q, element, P, Q)
      use simple_strings, only: upperCase
      real,             intent(inout) :: points_P(:,:), points_Q(:,:)
      character(len=2), intent(inout) :: element
      real,allocatable, intent(inout) :: P(:,:), Q(:,:) ! just the couples of points
      real    :: theoretical_radius ! for threshold selection
      integer :: n   ! max number of couples
      integer :: cnt ! actual number of couples
      integer :: i, location(1), cnt_P, cnt_Q
      real    :: d
      logical, allocatable :: mask(:)
      real,    allocatable :: points_P_out(:,:), points_Q_out(:,:)
      real,    parameter   :: ABSURD = -10000.
      if(size(points_P) == size(points_Q))then
        allocate(P(size(points_P, dim=1),size(points_P, dim=2)),source = points_P)
        allocate(Q(size(points_Q, dim=1),size(points_Q, dim=2)),source = points_Q)
        return ! they are already coupled
      endif
      ! To modify when we make a module for definitions
      element = upperCase(element)
      select case(element)
          case('PT')
              theoretical_radius = 1.1
              ! thoretical radius is already set
          case('PD')
              theoretical_radius = 1.12
          case('FE')
              theoretical_radius = 1.02
          case('AU')
              theoretical_radius = 1.23
          case default
              write(logfhandle, *)('Unknown atom element; find_couples')
              stop
       end select
       if(size(points_P, dim=2) < size(points_Q, dim=2)) then
         n  = size(points_P, dim=2)
         allocate(mask(n), source = .true.)
         allocate(points_P_out(3,n), points_Q_out(3,n), source = ABSURD) ! init to absurd value
         cnt = 0 ! counts the number of identified couples
         do i = 1, size(points_Q, dim=2)
            d = pixels_dist(points_Q(:,i),points_P(:,:),'min',mask,location,keep_zero=.true.)
            if(d <= 2.*theoretical_radius)then ! has a couple
                cnt = cnt + 1
                points_P_out(:,cnt) = points_P(:,location(1))
                points_Q_out(:,cnt) = points_Q(:,i)
                mask(location(1)) = .false. ! not to consider the same atom more than once
            endif
         enddo
         allocate(P(3,cnt), Q(3,cnt), source = 0.) ! correct size
         print *, 'n: ', n, 'cnt: ', cnt
         cnt_P = 0
         cnt_Q = 0
         do i = 1, n
            if(abs(points_P_out(1,i)-ABSURD) > TINY) then ! check just firts component
                cnt_P = cnt_P + 1
                P(:3,cnt_P) = points_P_out(:3,i)
            endif
            if(abs(points_Q_out(1,i)-ABSURD) > TINY) then
                cnt_Q = cnt_Q + 1
                Q(:3,cnt_Q) = points_Q_out(:3,i)
            endif
         enddo
       else ! size(points_P, dim=2) > size(points_Q, dim=2)
         n  = size(points_Q, dim=2)
         allocate(mask(n), source = .true.)
         allocate(points_P_out(3,n), points_Q_out(3,n), source = ABSURD) ! init to absurd value
         cnt = 0 ! counts the number of identified couples
         do i = 1, size(points_P, dim=2)
            d = pixels_dist(points_P(:,i),points_Q(:,:),'min',mask,location, keep_zero=.true.)
            if(d <= 2.*theoretical_radius)then ! has couple
                cnt = cnt + 1
                points_Q_out(:,cnt) = points_Q(:,location(1))
                points_P_out(:,cnt) = points_P(:,i)
                mask(location(1)) = .false. ! not to consider the same atom more than once
            endif
         enddo
         allocate(P(3,cnt), Q(3,cnt), source = 0.)
         cnt_P = 0
         cnt_Q = 0
         do i = 1, n
            if(abs(points_P_out(1,i)-ABSURD) > TINY) then
                cnt_P = cnt_P + 1
                P(:3,cnt_P) = points_P_out(:3,i)
            endif
            if(abs(points_Q_out(1,i)-ABSURD) > TINY) then
                cnt_Q = cnt_Q + 1
                Q(:3,cnt_Q) = points_Q_out(:3,i)
            endif
         enddo
       endif
       deallocate(mask, points_P_out, points_Q_out)
  end subroutine find_couples
  ! subroutine find_couples_dp(points_P, points_Q, element, P, Q)
  !     use simple_strings, only: upperCase
  !     real(dp),             intent(inout) :: points_P(:,:), points_Q(:,:)
  !     character(len=2),     intent(inout) :: element
  !     real(dp),allocatable, intent(inout) :: P(:,:), Q(:,:) ! just the couples of points
  !     real(dp)    :: theoretical_radius ! for threshold selection
  !     integer :: n   ! max number of couples
  !     integer :: cnt ! actual number of couples
  !     integer :: i, location(1), cnt_P, cnt_Q
  !     real    :: d
  !     logical,  allocatable :: mask(:)
  !     real(dp), allocatable :: points_P_out(:,:), points_Q_out(:,:)
  !     real(dp), parameter   :: ABSURD = -10000.d0
  !     if(size(points_P) == size(points_Q)) return ! they are already coupled
  !     ! To modify when we make a module for definitions
  !     element = upperCase(element)
  !     select case(element)
  !         case('PT')
  !             theoretical_radius = 1.1d0
  !             ! thoretical radius is already set
  !         case('PD')
  !             theoretical_radius = 1.12d0
  !         case('FE')
  !             theoretical_radius = 1.02d0
  !         case('AU')
  !             theoretical_radius = 1.23d0
  !         case default
  !             write(logfhandle, *)('Unknown atom element; find_couples')
  !             stop
  !      end select
  !      if(size(points_P, dim=2) < size(points_Q, dim=2)) then
  !        n  = size(points_P, dim=2)
  !        allocate(mask(n), source = .true.)
  !        allocate(points_P_out(3,n), points_Q_out(3,n), source = ABSURD) ! init to absurd value
  !        cnt = 0 ! counts the number of identified couples
  !        do i = 1, size(points_Q, dim=2)
  !           d = pixels_dist(real(points_Q(:,i),sp),real(points_P(:,:),sp),'min',mask,location,keep_zero=.true.)
  !           if(d <= 2.d0*theoretical_radius)then ! has a couple
  !               cnt = cnt + 1
  !               points_P_out(:,cnt) = points_P(:,location(1))
  !               points_Q_out(:,cnt) = points_Q(:,i)
  !               mask(location(1)) = .false. ! not to consider the same atom more than once
  !           endif
  !        enddo
  !        allocate(P(3,cnt), Q(3,cnt), source = 0.d0) ! correct size
  !        print *, 'n: ', n, 'cnt: ', cnt
  !        cnt_P = 0
  !        cnt_Q = 0
  !        do i = 1, n
  !           if(abs(points_P_out(1,i)-ABSURD) > TINY) then ! check just firts component
  !               cnt_P = cnt_P + 1
  !               P(:3,cnt_P) = points_P_out(:3,i)
  !           endif
  !           if(abs(points_Q_out(1,i)-ABSURD) > TINY) then
  !               cnt_Q = cnt_Q + 1
  !               Q(:3,cnt_Q) = points_Q_out(:3,i)
  !           endif
  !        enddo
  !      else ! size(points_P, dim=2) > size(points_Q, dim=2)
  !        n  = size(points_Q, dim=2)
  !        allocate(mask(n), source = .true.)
  !        allocate(points_P_out(3,n), points_Q_out(3,n), source = ABSURD) ! init to absurd value
  !        cnt = 0 ! counts the number of identified couples
  !        do i = 1, size(points_P, dim=2)
  !           d = pixels_dist(real(points_P(:,i),sp),real(points_Q(:,:),sp),'min',mask,location, keep_zero=.true.)
  !           if(d <= 2.d0*theoretical_radius)then ! has couple
  !               cnt = cnt + 1
  !               points_Q_out(:,cnt) = points_Q(:,location(1))
  !               points_P_out(:,cnt) = points_P(:,i)
  !               mask(location(1)) = .false. ! not to consider the same atom more than once
  !           endif
  !        enddo
  !        allocate(P(3,cnt), Q(3,cnt), source = 0.d0)
  !        cnt_P = 0
  !        cnt_Q = 0
  !        print *, 'n: ', n, 'cnt: ', cnt
  !        do i = 1, n
  !           if(abs(points_P_out(1,i)-ABSURD) > TINY) then
  !               cnt_P = cnt_P + 1
  !               P(:3,cnt_P) = points_P_out(:3,i)
  !           endif
  !           if(abs(points_Q_out(1,i)-ABSURD) > TINY) then
  !               cnt_Q = cnt_Q + 1
  !               Q(:3,cnt_Q) = points_Q_out(:3,i)
  !           endif
  !        enddo
  !      endif
  !      deallocate(mask, points_P_out, points_Q_out)
  ! end subroutine find_couples_dp

  ! This subrouine takes a pdb file input, removes all atoms beyond
  ! a given diameter, outputs a new pdb file with the coordinates
  ! removed and reports how many atoms were removed
  subroutine atoms_mask(pdb_file_in, max_rad, pdb_file_out, nremoved)
    use simple_fileio
    character(len=*), intent(in)    :: pdb_file_in
    character(len=*), intent(inout) :: pdb_file_out
    real,             intent(in)    :: max_rad
    integer,          intent(inout) :: nremoved
    real, allocatable    :: points_in(:,:)
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

end module simple_nano_utils
