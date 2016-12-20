  integer                                          :: n=5600, m=5600
  integer                                          :: k
  integer                                          :: lda, ldb, ldc
  integer                                          :: iseed
  double precision,dimension(:,:),allocatable      :: matDA
  double precision,dimension(:,:),allocatable      :: matDB
  double precision,dimension(:,:),allocatable      :: matDC
  double precision,dimension(:,:),allocatable      :: matDA_mpi
  double precision,dimension(:,:),allocatable      :: matDB_mpi
  double precision,dimension(:,:),allocatable      :: matDC_mpi
  double precision,dimension(nx,ny,nz,nt,mu,nc,nc) :: u1r_rnd,u1i_rnd
  double precision,dimension(nx,ny,nz,nt,mu,nc,nc) :: u2r_rnd,u2i_rnd
  double precision,dimension(nx,ny,nz,nt,mu,nc,nc) :: u3r,u3i

