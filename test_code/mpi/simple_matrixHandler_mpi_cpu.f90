!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 22th of October 2013.
!
! Name:
! matrixMultiplier_cpu - Various utilities and matrix multiplier to test on CPU
!
! Description:
! tester module provides test code.
! Subroutine: my_dgemm_cpu(m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
!
! \param[in] M
! \verbatim
!          M is INTEGER
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
! \endverbatim
!
! \param[in] N
! \verbatim
!          N is INTEGER
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
! \endverbatim
!
! \param[in] K
! \verbatim
!          K is INTEGER
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
! \endverbatim
!
! \param[in] A
! \verbatim
!          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
! \endverbatim
!
! \param[in] LDA
! \verbatim
!          LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
! \endverbatim
!
! \param[in] B
! \verbatim
!          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
! \endverbatim
!
! \param[in] LDB
! \verbatim
!          LDB is INTEGER
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
! \endverbatim
!
! \param[in,out] C
! \verbatim
!          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
! \endverbatim
!
! \param[in] LDC
! \verbatim
!          LDC is INTEGER
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
! \endverbatim
!*******************************************************************************
!
module simple_matrixHandler_cpu
!  use mpi
  use simple_defs
  use simple_timing
  use simple_random
  use matrixGetter

implicit none
#if defined (SIMPLE_MPI)
include 'mpif.h'
#endif 
contains
!*******************************************************************************
! DESCRIPTION
! subroutine to multiply dense mxn matrix
!*******************************************************************************
! SOURCE
  subroutine my_dgemm_cpu(m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
    implicit none
    !global variables
    integer                                    :: lda,ldb,ldc
    integer                                    :: n,m,k
    double precision                           :: matDA(lda,*)
    double precision                           :: matDB(ldb,*)
    double precision                           :: matDC(ldc,*)
    !local variables
    integer                                    :: in,im,ik
    !start execution commands

    matDC(1:m,1:n) = 0.0d0
    do im = 1,m
       do in = 1,n
          do ik = 1,k
             matDC(im,in) = matDC(im,in) + ( matDA(im,ik) * matDB(ik,in) )
          end do
       end do
    end do

    return
  end subroutine my_dgemm_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to multiply two SU(3) random gauge
!*******************************************************************************
! SOURCE
  subroutine su3randomlinks_mul_cpu(u1r_rnd,u1i_rnd, u2r_rnd,u2i_rnd, u3r,u3i)
    implicit none

    !global variables
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u1r_rnd,u1i_rnd
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u2r_rnd,u2i_rnd
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: u3r,u3i

    !local variables

    integer                         :: ic,jc,kc  !loop indices
    !start of the excution commands

    !serial multiplication
    u3r = 0.0d0
    u3i = 0.0d0
    do ic=1,nc
       do jc=1,nc
          do kc=1,nc
             u3r(:,:,:,:,:,ic,jc) = u3r(:,:,:,:,:,ic,jc) + &
                  ( u1r_rnd(:,:,:,:,:,ic,kc) * u2r_rnd(:,:,:,:,:,kc,jc) - &
                    u1i_rnd(:,:,:,:,:,ic,kc) * u2i_rnd(:,:,:,:,:,kc,jc) )
             u3i(:,:,:,:,:,ic,jc) = u3i(:,:,:,:,:,ic,jc) + &
                  ( u1r_rnd(:,:,:,:,:,ic,kc) * u2i_rnd(:,:,:,:,:,kc,jc) + &
                    u1i_rnd(:,:,:,:,:,ic,kc) * u2r_rnd(:,:,:,:,:,kc,jc) )
          end do
       end do
    end do

    return
  end subroutine su3randomlinks_mul_cpu

!*******************************************************************************
! DESCRIPTION
! subroutine to multiply dense mxn matrix using mpi distributed memory
!*******************************************************************************
! SOURCE
  subroutine my_dgemm_mpi_1(my_rank,n_proc,mpierr,m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
    !    use physsym_mpi_defs
    implicit none
#if defined (SIMPLE_MPI)
    include 'mpif.h'
#endif
    !global variables
    integer                                    :: my_rank      !rank of process
    integer                                    :: n_proc       !number of processes
    integer                                    :: mpierr 
    integer                                    :: lda,ldb,ldc
    integer                                    :: n,m,k
    double precision,intent(in)                :: matDA(lda,*)
    double precision,intent(in)                :: matDB(ldb,*)
    double precision,intent(inout)             :: matDC(ldc,*)
    !local variables
    integer                                    :: nworkers
    integer                                    :: nra, nca, ncb

    integer                                    :: extra, avecol, offset, cols
    !MPI variables
#if defined (SIMPLE_MPI)
    integer,parameter                          :: MASTER = 0
    integer,parameter                          :: FROM_MASTER = 1
    integer,parameter                          :: FROM_WORKER = 2
    integer                                    :: status(MPI_STATUS_SIZE)
#endif
    !The counters
    integer                                    :: im,in,ik
    integer                                    :: iworkers, isource
    !start execution commands
    nra=m
    nca=k
    ncb=n
#ifdef SIMPLE_MPI
    !initializing the matrix C
    !matDC(1:m,1:n) = 0.0d0
    !write(*,*) "Rank my_rank = ",my_rank,"n_proc = ",n_proc
    nworkers = n_proc-1
    !*************************  Master task **************************
    if ( my_rank .eq. MASTER ) then
       !send matrix data to the work tasks
       avecol = ncb / nworkers
       extra = mod(ncb, nworkers)
       offset = 1
       do iworkers = 1, nworkers
          if (iworkers .le. extra ) then
             cols = avecol + 1
          else
             cols = avecol
          end if
          !write(*,*)'   sending',cols,' cols to task',iworkers
          call MPI_SEND(offset         ,        1,          MPI_INTEGER, iworkers,FROM_MASTER, MPI_COMM_WORLD, mpierr)
          call MPI_SEND(cols           ,        1,          MPI_INTEGER, iworkers,FROM_MASTER, MPI_COMM_WORLD, mpierr)
          call MPI_SEND(matDA          ,  nra*nca, MPI_DOUBLE_PRECISION, iworkers,FROM_MASTER, MPI_COMM_WORLD, mpierr)
          call MPI_SEND(matDB(1,offset), cols*nca, MPI_DOUBLE_PRECISION, iworkers,FROM_MASTER, MPI_COMM_WORLD, mpierr)
          offset = offset + cols
       end do
       !write(*,*) "after sending data from master task"
       do isource = 1, nworkers
          !isource = i
          call MPI_RECV(offset,                 1,          MPI_INTEGER, isource, FROM_WORKER, MPI_COMM_WORLD, status, mpierr)
          call MPI_RECV(cols,                   1,          MPI_INTEGER, isource, FROM_WORKER, MPI_COMM_WORLD, status, mpierr)
          call MPI_RECV(matDC(1,offset), cols*nra, MPI_DOUBLE_PRECISION, isource, FROM_WORKER, MPI_COMM_WORLD, status, mpierr)
       end do

    endif

    !*************************  Worker task **************************
    if ( my_rank > MASTER ) then
       !receive matrix data from master task
       call MPI_RECV(offset,       1,          MPI_INTEGER, MASTER, FROM_MASTER, MPI_COMM_WORLD, status, mpierr)
       call MPI_RECV(cols,         1,          MPI_INTEGER, MASTER, FROM_MASTER, MPI_COMM_WORLD, status, mpierr)
       call MPI_RECV(matDA,  nra*nca, MPI_DOUBLE_PRECISION, MASTER, FROM_MASTER, MPI_COMM_WORLD, status, mpierr)
       call MPI_RECV(matDB, cols*nca, MPI_DOUBLE_PRECISION, MASTER, FROM_MASTER, MPI_COMM_WORLD, status, mpierr)

       do ik=1, cols
          do im=1, NRA
             matDC(im,ik) = 0.0
             do in=1, NCA
                matDC(im,ik) = matDC(im,ik) + matDA(im,in) * matDB(in,ik)
             end do
          end do
       end do
       !send result back to master task
       call MPI_SEND(offset,       1,          MPI_INTEGER, MASTER, FROM_WORKER, MPI_COMM_WORLD, mpierr)
       call MPI_SEND(cols,         1,          MPI_INTEGER, MASTER, FROM_WORKER, MPI_COMM_WORLD, mpierr)
       call MPI_SEND(matDC, cols*nra, MPI_DOUBLE_PRECISION, MASTER, FROM_WORKER, MPI_COMM_WORLD, mpierr)
    end if
#else 
    write(*,*)"***************************WARNING************************************"
    write(*,*)"You need to compile with -DSIMPLE_MPI to acces the mpi computtation  "
    write(*,*)"Proceeding with standard 1-CPU calculation using my_degemm_cpu routine"
    write(*,*)"**********************************************************************"
    call my_dgemm_cpu(m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
#endif
    return
  end subroutine my_dgemm_mpi_1
!*******************************************************************************
! DESCRIPTION
! subroutine to multiply dense mxn matrix using mpi distributed memory
!*******************************************************************************
! SOURCE
  subroutine my_dgemm_mpi_2(my_rank,n_proc,mpierr,m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
    !use mpi
    !use physsym_mpi_defs
    implicit none 
#if defined (SIMPLE_MPI)
    include 'mpif.h'
#endif    
    !global variables
    integer                                    :: my_rank      !rank of process
    integer                                    :: n_proc       !number of processes
    integer                                    :: mpierr 
    integer                                    :: lda,ldb,ldc
    integer                                    :: n,m,k
    double precision,intent(in)                :: matDA(lda,*)
    double precision,intent(in)                :: matDB(ldb,*)
    double precision,intent(inout)             :: matDC(ldc,*)
    !local variables
    integer                                    :: nworkers
    integer                                    :: numset
    integer                                    :: sender, anstype, row
    integer                                    :: matDA_rows, matDA_cols
    integer                                    :: matDB_rows, matDB_cols
    integer                                    :: matDC_rows, matDC_cols
    double precision,dimension(:),allocatable  :: buffer, ans
    !MPI variables
#if defined (SIMPLE_MPI)
    integer                                    :: status(MPI_STATUS_SIZE)
    integer,parameter                          :: MASTER = 0
    integer,parameter                          :: FROM_MASTER = 1
    integer,parameter                          :: FROM_WORKER = 2
#endif
    !The counters
    integer                                    :: im,in,ik
    integer                                    :: i_proc
    !start of the execution commands
    matDA_rows = m
    matDA_cols = k
    matDB_rows = k
    matDB_cols = n
    matDC_rows = m
    matDC_cols = n
    allocate(buffer(m))
    allocate(ans(n))

#ifdef SIMPLE_MPI

    nworkers = n_proc - 1
    !*************************  Master task **************************
    if ( my_rank .eq. MASTER ) then
       !send matrix data to the work tasks
       numset = 0
       do in=1,matDB_cols
          call mpi_bcast(matDB(1,in),matDB_rows,MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD,mpierr)
       end do

       !send a row of matDA to each of the processors
       do i_proc=1,nworkers
          do im=1,matDA_cols
             buffer(im) = matDA(i_proc,im)
          end do
          call mpi_send(buffer, matDA_cols, MPI_DOUBLE_PRECISION, i_proc, i_proc, MPI_COMM_WORLD, mpierr )
          numset = numset + 1
       end do

       do im=1,matDC_rows
          call mpi_recv(ans, matDC_cols, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
          sender = status(MPI_SOURCE)
          anstype = status(MPI_TAG)
          do in=1,matDC_cols
             matDC(anstype,in) = ans(in)
          end do

          if (numset .lt. matDA_rows ) then
             do ik = 1,matDA_cols
                buffer(ik) = matDA(numset+1,ik)
             end do
             call mpi_send(buffer, matDA_cols, MPI_DOUBLE_PRECISION, sender, numset+1, MPI_COMM_WORLD, mpierr)
             numset = numset + 1
          else
             call mpi_send(1.0, 1, MPI_DOUBLE_PRECISION ,sender, 0, MPI_COMM_WORLD, mpierr)
          end if

       end do
    endif
    !*************************  Worker task **************************
    if ( my_rank > MASTER ) then
       !workers receives B and computes rows of C until message
       do in=1,matDB_cols
          call mpi_bcast(matDB(1,in), matDB_rows, MPI_DOUBLE_PRECISION, MASTER, MPI_COMM_WORLD, mpierr )
       end do
90     call MPI_RECV(buffer, matDA_cols, MPI_DOUBLE_PRECISION, master, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
       if (status(MPI_TAG) .eq. 0) then
          go to 200
       else
          row = status(MPI_TAG)
          do in = 1,matDB_cols
             ans(in) = 0.0
             do ik = 1,matDA_cols
                ans(in) = ans(in) + buffer(ik)*matDB(ik,in)
             end do
          end do
          call MPI_SEND(ans, matDB_cols, MPI_DOUBLE_PRECISION, master, row, MPI_COMM_WORLD, mpierr)
          go to 90
       endif
200    continue

    end if

    !freeing the memory
    deallocate(ans)
    deallocate(buffer)
#else 
    write(*,*)"***************************WARNING************************************"
    write(*,*)"You need to compile with -DSIMPLE_MPI to acces the mpi computtation  "
    write(*,*)"Proceeding with standard 1-CPU calculation using my_degemm_cpu routine"
    write(*,*)"**********************************************************************"
    call my_dgemm_cpu(m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
#endif

  return
end subroutine my_dgemm_mpi_2

end module simple_matrixHandler_cpu
