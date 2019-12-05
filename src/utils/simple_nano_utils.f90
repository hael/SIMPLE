module simple_nano_utils
  use simple_image,    only : image
  use simple_binimage, only : binimage
  use simple_rnd
  use simple_math
  use simple_defs ! singleton
  use simple_atoms,    only : atoms
  implicit none

  public :: remove_graphene_peaks, kabsch, read_pdb2matrix, write_matrix2pdb, find_couples
  private

  interface kabsch
      module procedure kabsch_sp, kabsch_dp
  end interface kabsch

  interface read_pdb2matrix
      module procedure read_pdb2matrix_sp, read_pdb2matrix_dp
  end interface read_pdb2matrix

  interface find_couples
      module procedure find_couples_sp, find_couples_dp
  end interface find_couples

contains

  ! This subroutine generates a mask that identifies the graphene
  ! peaks in the spectrum.
  subroutine generate_peakmask(img_spec, msk)
    type(image),          intent(inout) :: img_spec
    integer, allocatable, intent(inout) :: msk(:,:,:)
    type(binimage)       :: img_spec_roavg
    integer, parameter   :: WIDTH   =2 ! width of the shell for peak identification
    integer, parameter   :: N_LAYERS=1 ! nb of layers to grow
    real,    pointer     :: rmat(:,:,:)
    logical, allocatable :: lmsk(:)
    real, allocatable    :: x(:),rmat_aux(:,:,:)
    real                 :: thresh, smpd
    integer              :: ldim(3),i,j,h,k,sh,nframe,cnt,UFlim,LFlim,lims(2),nyq
    type(binimage)       :: img_mask
    nthr_glob=8 !! TO REMOVE AFTER TESTING
    if(allocated(msk)) deallocate(msk)
    ldim = img_spec%get_ldim()
    ldim(3) = 1 ! stack
    smpd = img_spec%get_smpd()
    ! calculate rotational avg of the spectrum
    call img_spec_roavg%new_bimg(ldim,smpd)
    call img_spec%roavg(60,img_spec_roavg)
    call img_spec_roavg%get_rmat_ptr(rmat)
    allocate(rmat_aux(1:ldim(1),1:ldim(2),1), source= 0.)
    rmat_aux(1:ldim(1),1:ldim(2),1)=rmat(1:ldim(1),1:ldim(2),1)
    UFlim   = calc_fourier_index(max(2.*smpd,GRAPHENE_BAND1),ldim(1),smpd)
    lims(1) = UFlim-WIDTH
    lims(2) = UFlim+WIDTH
    cnt = 0
    do i = 1, ldim(1)
      do j = 1, ldim(2)
        h   = -int(ldim(1)/2) + i - 1
        k   = -int(ldim(2)/2) + j - 1
        sh  =  nint(hyp(real(h),real(k)))
        if(sh < lims(2) .and. sh > lims(1)) then
          !do nothing
        else
          rmat_aux(i,j,1) = 0. ! set to zero everything outside the considered shell
        endif
      enddo
    enddo
    x = pack(rmat_aux, abs(rmat_aux) > TINY)
    call otsu(x, thresh)
    deallocate(x)
    allocate(msk(ldim(1),ldim(2),1), source =0)
    where(rmat(1:ldim(1),1:ldim(2),1)>=thresh .and. abs(rmat_aux(:,:,1)) > TINY)  msk(:,:,1) = 1
    ! Do the same thing for the other band
    ! reset
    rmat_aux(:,:,1) = rmat(1:ldim(1),1:ldim(2),1)
    LFlim   = calc_fourier_index(max(2.*smpd,GRAPHENE_BAND2),ldim(1),smpd)
    lims(1) = LFlim-WIDTH
    lims(2) = LFlim+WIDTH
    do i = 1, ldim(1)
      do j = 1, ldim(2)
        h   = -int(ldim(1)/2) + i - 1
        k   = -int(ldim(2)/2) + j - 1
        sh  =  nint(hyp(real(h),real(k)))
        if(sh < lims(2) .and. sh > lims(1)) then
          ! do nothing
        else
          rmat_aux(i,j,1) = 0.
        endif
      enddo
    enddo
    x = pack(rmat_aux, abs(rmat_aux) > 0.) ! pack where it's not 0
    call otsu(x, thresh)
    where(rmat(1:ldim(1),1:ldim(2),1)>=thresh  .and. abs(rmat_aux(:,:,1)) > TINY) msk(:,:,1) = 1
    call img_mask%new_bimg(ldim, smpd)
    call img_mask%set_imat(msk)
    call img_mask%update_img_rmat()
    call img_mask%get_rmat_ptr(rmat)
    where(rmat<TINY)
      rmat(:,:,:) = 0.
    elsewhere
      rmat(:,:,:) = 1.
    endwhere
    call img_mask%set_imat()
    call img_mask%grow_bins(N_LAYERS)
    call img_mask%get_imat(msk)
    ! filter results with graphene-bands mask
    lmsk = calc_graphene_mask( ldim(1), smpd )
    nyq  = size(lmsk)
    do i = 1, ldim(1)
      do j = 1, ldim(2)
        h   =  i - (ldim(1)/2 + 1)
        k   =  j - (ldim(1)/2 + 1)
        sh  =  nint(hyp(real(h),real(k)))
        if( sh > nyq .or. sh == 0 ) then
          cycle
        else
          if(lmsk(sh)) msk(i,j,1) = 0 ! reset
        endif
      enddo
    enddo
    call img_mask%kill_bimg
    call img_spec_roavg%kill_bimg
  end subroutine generate_peakmask


  ! This subroutine operates on raw_img on the basis of the
  ! provided logical mask. For each pxl in the msk, it subsitutes
  ! its value with the avg value in the neighbours of the pxl that
  ! do not belong to the msk. The neighbours are in a shpere of radius RAD.
  subroutine filter_peaks(raw_img,msk)
    type(image), intent(inout) :: raw_img
    integer,     intent(in)    :: msk(:,:,:)
    complex, pointer   :: cmat(:,:,:)
    real,    parameter :: RAD = 4.
    real        :: rabs, ang
    integer     :: nptcls,i,j,ii,jj,ldim(3),cnt,hp,h,k,kp,phys(3)
    ldim  = raw_img%get_ldim()
    call raw_img%norm()
    call raw_img%fft()
    call raw_img%get_cmat_ptr(cmat)
    do i = ldim(1)/2+1, ldim(1)
      do j = 1, ldim(2)
        if(msk(i,j,1)>0) then
          cnt = 0 ! number of pxls in the sphere that are not flagged by msk
          rabs = 0.
          do ii = 1, ldim(1)
            h  = ii-(ldim(1)/2+1)
            do jj = 1, ldim(2)
              if(ii/=i .and. jj/=j) then
                if(euclid(real([i,j]), real([ii,jj])) <= RAD) then
                  if(msk(ii,jj,1)<1) then
                    cnt = cnt + 1
                    k = jj-(ldim(2)/2+1)
                    phys = raw_img%comp_addr_phys(h,k,0)
                    rabs = rabs + sqrt(csq(cmat(phys(1),phys(2),1))) ! magnitude
                  endif
                endif
              endif
            enddo
          enddo
          h = i-(ldim(1)/2+1)
          k = j-(ldim(2)/2+1)
          phys = raw_img%comp_addr_phys(h,k,0)
          if(cnt == 0) then
            ! random number
            rabs = ran3()
          else
            rabs = rabs/real(cnt)
          endif
          ang = TWOPI*ran3()
          cmat(phys(1),phys(2),1) = 0.! rabs * cmplx(cos(ang), sin(ang))
        endif
      enddo
    enddo
    call raw_img%ifft()
  end subroutine filter_peaks

  subroutine remove_graphene_peaks(raw_img, spectrum)
    type(image), intent(inout) :: raw_img
    type(image), intent(inout) :: spectrum
    integer, allocatable       :: msk(:,:,:)
    call generate_peakmask(spectrum, msk)
    call filter_peaks(raw_img,msk)
  end subroutine remove_graphene_peaks

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
       write(logfhandle,*) 'WORNG input dims'
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

  subroutine find_couples_sp(points_P, points_Q, element, P, Q)
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
      if(size(points_P) == size(points_Q)) return ! they are already coupled
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
         print *, 'n: ', n, 'cnt: ', cnt
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
  end subroutine find_couples_sp


  subroutine find_couples_dp(points_P, points_Q, element, P, Q)
      use simple_strings, only: upperCase
      real(dp),             intent(inout) :: points_P(:,:), points_Q(:,:)
      character(len=2),     intent(inout) :: element
      real(dp),allocatable, intent(inout) :: P(:,:), Q(:,:) ! just the couples of points
      real(dp)    :: theoretical_radius ! for threshold selection
      integer :: n   ! max number of couples
      integer :: cnt ! actual number of couples
      integer :: i, location(1), cnt_P, cnt_Q
      real    :: d
      logical,  allocatable :: mask(:)
      real(dp), allocatable :: points_P_out(:,:), points_Q_out(:,:)
      real(dp), parameter   :: ABSURD = -10000.d0
      if(size(points_P) == size(points_Q)) return ! they are already coupled
      ! To modify when we make a module for definitions
      element = upperCase(element)
      select case(element)
          case('PT')
              theoretical_radius = 1.1d0
              ! thoretical radius is already set
          case('PD')
              theoretical_radius = 1.12d0
          case('FE')
              theoretical_radius = 1.02d0
          case('AU')
              theoretical_radius = 1.23d0
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
            d = pixels_dist(real(points_Q(:,i),sp),real(points_P(:,:),sp),'min',mask,location,keep_zero=.true.)
            if(d <= 2.d0*theoretical_radius)then ! has a couple
                cnt = cnt + 1
                points_P_out(:,cnt) = points_P(:,location(1))
                points_Q_out(:,cnt) = points_Q(:,i)
                mask(location(1)) = .false. ! not to consider the same atom more than once
            endif
         enddo
         allocate(P(3,cnt), Q(3,cnt), source = 0.d0) ! correct size
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
            d = pixels_dist(real(points_P(:,i),sp),real(points_Q(:,:),sp),'min',mask,location, keep_zero=.true.)
            if(d <= 2.d0*theoretical_radius)then ! has couple
                cnt = cnt + 1
                points_Q_out(:,cnt) = points_Q(:,location(1))
                points_P_out(:,cnt) = points_P(:,i)
                mask(location(1)) = .false. ! not to consider the same atom more than once
            endif
         enddo
         allocate(P(3,cnt), Q(3,cnt), source = 0.d0)
         cnt_P = 0
         cnt_Q = 0
         print *, 'n: ', n, 'cnt: ', cnt
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
  end subroutine find_couples_dp
end module simple_nano_utils
