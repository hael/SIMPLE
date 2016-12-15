!==module simple_math_gpu
!
! simple_math_gpu contains various mathematical subroutines and functions.
! the code is distributed with the hope that it will be useful,
! but _without_ _any_ _warranty_. redistribution or modification is regulated
! by the ! gnu general public license. *author:* Frederic Bonnet, 2015-04-22.
! 
!==changes are documented below
!* incorporated in the _simple_ library, the 2015-04-22
!
!> \brief \b get_matmul
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
! 
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*op( A )*op( B ) + beta*C,
!>
!> where  op( X ) is one of
!>
!>    op( X ) = X   or   op( X ) = X**T,
!>
!> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n',  op( A ) = A.
!>
!>              TRANSA = 'T' or 't',  op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c',  op( A ) = A**T.
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER*1
!>           On entry, TRANSB specifies the form of op( B ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSB = 'N' or 'n',  op( B ) = B.
!>
!>              TRANSB = 'T' or 't',  op( B ) = B**T.
!>
!>              TRANSB = 'C' or 'c',  op( B ) = B**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies  the number  of rows  of the  matrix
!>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N  specifies the number  of columns of the matrix
!>           op( B ) and the number of columns of the matrix C. N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry,  K  specifies  the number of columns of the matrix
!>           op( A ) and the number of rows of the matrix op( B ). K must
!>           be at least  zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by m  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  n by k  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!>           least  max( 1, n ).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
!>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!>           supplied as zero then C need not be set on input.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!>           Before entry, the leading  m by n  part of the array  C must
!>           contain the matrix  C,  except when  beta  is zero, in which
!>           case C need not be set on entry.
!>           On exit, the array  C  is overwritten by the  m by n  matrix
!>           ( alpha*op( A )*op( B ) + beta*C ).
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!
!  =====================================================================
module simple_math_gpu
  use simple_defs
  use simple_cuda_defs

  implicit none
#define devptr_t integer*8

contains

  !> \brief  generates polar coordinates
  subroutine gen_polar_coords_gpu( kfromto, ring2, coords, angtab )
    integer, intent(in)            :: kfromto(2), ring2
    real, allocatable, intent(out) :: coords(:,:,:), angtab(:) 

#if defined (CUDA)
    call gen_polar_coords_cuda_gpu( kfromto, ring2, coords, angtab )
#else
    call get_invert_gpu_utils_warning()
#endif

    return
  end subroutine gen_polar_coords_gpu

! LINEAR ALGEBRA STUFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix vector product of a Single precision 2d array
! of size y(n) = alpha * A(m,n)*x(n) + beta*y(n)
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_SMatvec(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    use matvec_cpu_utils
    use matvec_gpu_utils
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    real(sp)                      :: alpha,beta

    real(sp)                      :: A(lda,*), x(*)
    real(sp)                      :: y(*)
    character                     :: trans

#if defined (OPENMP) && !defined (CUDA)
    !$omp parallel default(shared)
    !$omp workshare
    y = matmul(A,x)
    !$omp end workshare nowait
    !$omp end parallel
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    call get_SMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    ! else defautlting on the standard lapack
    call get_SMatvec_cpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#endif    

    return
  end subroutine my_SMatvec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix vector product of a Single precision 2d array
! of size y(n) = alpha * A(m,n)*x(n) + beta*y(n)
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_DMatvec(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    use matvec_cpu_utils
    use matvec_gpu_utils
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    real(dp)                      :: alpha,beta

    real(dp)                      :: A(lda,*), x(*)
    real(dp)                      :: y(*)
    character                     :: trans

#if defined (OPENMP) && !defined (CUDA)
    !$omp parallel default(shared)
    !$omp workshare
    y = matmul(A,x)
    !$omp end workshare nowait
    !$omp end parallel
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    call get_DMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    ! else defautlting on the standard lapack
    call get_DMatvec_cpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#endif    

    return
  end subroutine my_DMatvec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix vector product of a Single complex precision 2d
! array
! of size y(n) = alpha * A(m,n)*x(n) + beta*y(n)
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_CMatvec(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    use matvec_gpu_utils
    use matvec_cpu_utils
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    complex(sp)                   :: alpha,beta

    complex(sp)                   :: A(lda,*), x(*)
    complex(sp)                   :: y(*)
    character                     :: trans

#if defined (OPENMP) && !defined (CUDA)
    !$omp parallel default(shared)
    !$omp workshare
    y = matmul(A,x)
    !$omp end workshare nowait
    !$omp end parallel
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    call get_CMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    ! else defautlting on the standard lapack
    call get_CMatvec_cpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#endif    

    return
  end subroutine my_CMatvec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix vector product of a Double complex precision 2d
! array
! of size y(n) = alpha * A(m,n)*x(n) + beta*y(n)
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_ZMatvec(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    use matvec_cpu_utils
    use matvec_gpu_utils
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    complex(dp)                   :: alpha,beta

    complex(dp)                   :: A(lda,*), x(*)
    complex(dp)                   :: y(*)
    character                     :: trans

#if defined (OPENMP) && !defined (CUDA)
    !$omp parallel default(shared)
    !$omp workshare
    y = matmul(A,x)
    !$omp end workshare nowait
    !$omp end parallel
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    call get_ZMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    ! else defautlting on the standard lapack
    call get_ZMatvec_cpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#endif    

    return
  end subroutine my_ZMatvec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single precision 2d array of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_SMatmul(transa, transb, m, n, k, &
                        alpha,                   &
                        A, lda,                  &
                        B, ldb,                  &
                        beta,                    &
                        prod_gpu,ldc             )
    use matmul_gpu_utils
    use matmul_cpu_utils
    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    real(sp)                      :: alpha,beta
    
    real(sp)                      :: A(lda,*)
    real(sp)                      :: B(ldb,*)
    real(sp)                      :: prod_gpu(ldc,*)
    character                     :: transa,transb
    
#if defined (OPENMP) && !defined (CUDA)
    !TODO: here put the openMP code for the matmul
    !$omp parallel default(shared)
    !$omp workshare
    prod_gpu = matmul(A,B)
    !$omp end workshare nowait
    !$omp end parallel
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    write(*,*)"call get_SMatmul_gpu with m: ",m," n: ",n," and k: ",k
    call get_SMatmul_gpu(transa,transb,m,n,k,    &
                        alpha,                   &
                        A,lda,                   &
                        B,ldb,                   &
                        beta,                    &
                        prod_gpu,ldc             )
    write(*,*)"end of get_SMatmul_gpu"
#else
    ! else defautlting on the standard lapack
    call get_SMatmul_cpu(transa,transb,m,n,k,&
                         alpha,              &
                         A,lda,              &
                         B,ldb,              &
                         beta,               &
                         prod_gpu,ldc)
#endif    
    
    return
  end subroutine my_SMatmul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single precision 2d array of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_DMatmul(transa, transb, m, n, k, &
                        alpha,                   &
                        A, lda,                  &
                        B, ldb,                  &
                        beta,                    &
                        prod_gpu,ldc             )
    use matmul_gpu_utils
    use matmul_cpu_utils
    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    real(dp)                      :: alpha,beta
    
    real(dp)                      :: A(lda,*)
    real(dp)                      :: B(ldb,*)
    real(dp)                      :: prod_gpu(ldc,*)
    character                     :: transa,transb
    
#if defined (OPENMP) && !defined (CUDA)
    !TODO: here put the openMP code for the matmul
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    call get_DMatmul_gpu(transa,transb,m,n,k,    &
                        alpha,                   &
                        A,lda,                   &
                        B,ldb,                   &
                        beta,                    &
                        prod_gpu,ldc             )
#else
    ! else defautlting on the standard lapack
    call get_DMatmul_cpu(transa,transb,m,n,k,&
                         alpha,              &
                         A,lda,              &
                         B,ldb,              &
                         beta,               &
                         prod_gpu,ldc)
#endif    
    
    return
  end subroutine my_DMatmul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single complex precision 2d array
! of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_CMatmul(transa, transb, m, n, k, &
                        alpha,                   &
                        A, lda,                  &
                        B, ldb,                  &
                        beta,                    &
                        prod_gpu,ldc             )
    use matmul_gpu_utils
    use matmul_cpu_utils
    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    complex(sp)                   :: alpha,beta
    
    complex(sp)                   :: A(lda,*)
    complex(sp)                   :: B(ldb,*)
    complex(sp)                   :: prod_gpu(ldc,*)
    character                     :: transa,transb
    
#if defined (OPENMP) && !defined (CUDA)
    !TODO: here put the openMP code for the matmul
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    call get_CMatmul_gpu(transa,transb,m,n,k,    &
                        alpha,                   &
                        A,lda,                   &
                        B,ldb,                   &
                        beta,                    &
                        prod_gpu,ldc             )
#else
    ! else defautlting on the standard lapack
    call get_CMatmul_cpu(transa,transb,m,n,k,&
                         alpha,              &
                         A,lda,              &
                         B,ldb,              &
                         beta,               &
                         prod_gpu,ldc)
#endif    
    
    return
  end subroutine my_CMatmul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Double complex precision 2d array
! of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine my_ZMatmul(transa, transb, m, n, k, &
                        alpha,                   &
                        A, lda,                  &
                        B, ldb,                  &
                        beta,                    &
                        prod_gpu,ldc             )
    use matmul_gpu_utils
    use matmul_cpu_utils
    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    complex(dp)                   :: alpha,beta
    
    complex(dp)                   :: A(lda,*)
    complex(dp)                   :: B(ldb,*)
    complex(dp)                   :: prod_gpu(ldc,*)
    character                     :: transa,transb
    
#if defined (OPENMP) && !defined (CUDA)
    !TODO: here put the openMP code for the matmul
#elif defined (CUDA)
    !calling the cublas for thge matrix product
    call get_ZMatmul_gpu(transa,transb,m,n,k,    &
                        alpha,                   &
                        A,lda,                   &
                        B,ldb,                   &
                        beta,                    &
                        prod_gpu,ldc             )
#else
    ! else defautlting on the standard lapack
    call get_ZMatmul_cpu(transa,transb,m,n,k,&
                         alpha,              &
                         A,lda,              &
                         B,ldb,              &
                         beta,               &
                         prod_gpu,ldc)
#endif    
    
    return
  end subroutine my_ZMatmul
!*******************************************************************************
! DESCRIPTION
! subroutine to take the determinante of a single precision 2d array of size
! (n)x(n) on GPU 
!
!*******************************************************************************
! SOURCE
  !> \brief  subroutine to find the determinant of a Single precision
  !!         square matrix on GPU
  !!         author : Frederic Bonnet
  function det_gpu(matrix, n, errflg) result (d)
    use invert_gpu_utils
    use invert_cpu_utils
    use simple_math
    integer, intent(out)  :: errflg  !return code. -1 for error, 0 for normal
    integer, intent(in)   :: n
    real(sp), dimension(n,n)  :: matrix
    real(sp)                  :: d
    !local variable

    devptr_t              :: devPtrA_new
    integer               :: lda
    integer               :: rc = RC_SUCCESS
    !counters 
    integer               :: i

    lda = n

#if defined (OPENMP) && !defined (CUDA)
    !TODO: here put the openMP code for the matinv
     d = det(matrix, n)
#elif defined (CUDA) && defined (MAGMA)
     rc = cublas_alloc(n*n, size_of_double, devPtrA_new)
     if (rc .ne. 0 ) call simple_cudblas_stat_return(rc)

     rc = cublas_set_matrix (n,n,size_of_double,dble(matrix),n,devPtrA_new,n)
     if ( rc .ne. 0 ) call simple_cudblas_stat_return(rc)

     call getLU_DMatrix_gpu(n,n,lda,size_of_double,n,devPtrA_new,dble(matrix))

     d = 0.0 !initializing the determinant of matrix
     do i=1,n
        d = d*matrix(i,i) !returning the determinant of the triangular matrix
     end do

     !freeing ressources on the device
     rc = cublas_free(devPtrA_new)

     errflg = rc
#else
     call getLU_SMatrix_cpu(n,n,lda,size_of_double,n,matrix)
     d = 0.0 !initializing the determinant of matrix
     do i=1,n
        d = d*matrix(i,i) !returning the determinant of the triangular matrix
     end do
#endif

  end function det_gpu
!*******************************************************************************
! DESCRIPTION
! subroutine to take the determinante of a single precision 2d array of size
! (n)x(n) on GPU 
!
!*******************************************************************************
! SOURCE
  !> \brief  subroutine to find the determinant of a Single precision
  !!         square matrix on GPU
  !!         author : Frederic Bonnet
  function det_D_gpu(matrix, n, errflg) result (d)
    use invert_gpu_utils
    use invert_cpu_utils
    use simple_math
    integer, intent(out)  :: errflg  !return code. -1 for error, 0 for normal
    integer, intent(in)   :: n
    real(dp), dimension(n,n)  :: matrix
    real(dp)                  :: d
    !local variable

    devptr_t              :: devPtrA_new
    integer               :: lda
    integer               :: rc = RC_SUCCESS
    !counters 
    integer               :: i

    lda = n

#if defined (OPENMP) && !defined (CUDA)
    !TODO: here put the openMP code for the matinv
     d = det_D(matrix, n)
#elif defined (CUDA) && defined (MAGMA)
     rc = cublas_alloc(n*n, size_of_double, devPtrA_new)
     if (rc .ne. 0 ) call simple_cudblas_stat_return(rc)

     rc = cublas_set_matrix (n,n,size_of_double,dble(matrix),n,devPtrA_new,n)
     if ( rc .ne. 0 ) call simple_cudblas_stat_return(rc)

     call getLU_DMatrix_gpu(n,n,lda,size_of_double,n,devPtrA_new,matrix)

     d = 0.0 !initializing the determinant of matrix
     do i=1,n
        d = d*matrix(i,i) !returning the determinant of the triangular matrix
     end do

     !freeing ressources on the device
     rc = cublas_free(devPtrA_new)

     errflg = rc
#else
     call getLU_DMatrix_cpu(n,n,lda,size_of_double,n,matrix)
     d = 0.0 !initializing the determinant of matrix
     do i=1,n
        d = d*matrix(i,i) !returning the determinant of the triangular matrix
     end do
#endif

   end function det_D_gpu
!*******************************************************************************
! DESCRIPTION
! subroutine to take the inverse of a real precision 2d array of size (m)x(n)
! on GPU 
!
!*******************************************************************************
! SOURCE
  !> \brief  subroutine to find the inverse of a Single precision
  !!         square matrix on GPU
  !!         author : Frederic Bonnet
  subroutine matinv_S_gpu(matrix, inverse, n, errflg)
    use invert_gpu_utils
    use invert_cpu_utils
    use simple_math
    integer, intent(in)  :: n
    integer, intent(out) :: errflg  !return code. -1 for error, 0 for normal
    real, intent(in), dimension(n,n)  :: matrix  !input matrix
    real, intent(out), dimension(n,n) :: inverse !inverted matrix
    !local variables

    devptr_t                        :: devPtrA_new
    integer                         :: lda
    integer                         :: rc = RC_SUCCESS

    !start of the execution commands

    ! now using the magma to get iunverse matrix

    lda = n

#if defined (OPENMP) && !defined (CUDA)
    !$omp parallel default(shared)
    !$omp workshare
    call matinv(matrix, inverse, n, errflg)
    !$omp end workshare nowait
    !$omp end parallel
#elif defined (CUDA) && defined (MAGMA)
    rc = cublas_alloc(n*n, size_of_double, devPtrA_new)
    if (rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    rc = cublas_set_matrix(n,n,size_of_double,dble(matrix),n,devPtrA_new,n)
    if ( rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    inverse = 0.0d0 !initializing the inverted matrix
    call getBlockInverseDMatrix_gpu(n,n,lda,size_of_double,n,devPtrA_new,dble(inverse))

    rc = cublas_free(devPtrA_new)

    errflg = rc
#else
    call getBlockInverseSMatrix_cpu(n,n,lda,size_of_double,n,matrix,inverse)
#endif
 
    return
  end subroutine matinv_S_gpu
!*******************************************************************************
! DESCRIPTION
! subroutine to take the inverse of a Double precision 2d array of size (m)x(n)
! on GPU
!
!*******************************************************************************
! SOURCE
  !> \brief  subroutine to find the inverse of a Double precision
  !!         square matrix on GPU
  !!         author : Frederic Bonnet
  subroutine matinv_D_gpu(matrix, inverse, n, errflg)
    use invert_gpu_utils
    use invert_cpu_utils
    use simple_math
    integer, intent(in)  :: n
    integer, intent(out) :: errflg  !return code. -1 for error, 0 for normal
    real(dp), intent(in), dimension(n,n)  :: matrix  !input matrix
    real(dp), intent(out), dimension(n,n) :: inverse !inverted matrix
    !local variables

    devptr_t                        :: devPtrA_new
    integer                         :: lda
    integer                         :: rc = RC_SUCCESS

    !start of the execution commands

    lda = n

#if defined (OPENMP) && !defined (CUDA)
    !$omp parallel default(shared)
    !$omp workshare
    call matinv_D(matrix,inverse,n,errflg)
    !$omp end workshare nowait
    !$omp end parallel
#elif defined (CUDA) && defined (MAGMA)
    ! now using the magma to get inverse matrix
    rc = cublas_alloc(n*n, size_of_double, devPtrA_new)
    if (rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    rc = cublas_set_matrix (n,n,size_of_double,matrix,n,devPtrA_new,n)
    if ( rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    inverse = 0.0d0 !initializing the inverted matrix
    call getBlockInverseDMatrix_gpu(n,n,lda,size_of_double,n,devPtrA_new, &
                                    inverse)

    rc = cublas_free(devPtrA_new)

    errflg = rc
#else
    call getBlockInverseDMatrix_cpu(n,n,lda,size_of_double,n,matrix,inverse)
#endif
 
    return
  end subroutine matinv_D_gpu
!*******************************************************************************
! DESCRIPTION
! subroutine to take the inverse of a DoubleComplex precision 2d array of
! size (m)x(n) on GPU 
!
!*******************************************************************************
! SOURCE
  !> \brief  subroutine to find the inverse of a DoubleComplex precision
  !!         square matrix on GPU
  !!         author : Frederic Bonnet
  subroutine matinv_Z_gpu(matrix, inverse, n, errflg)
    use invert_gpu_utils
    use invert_cpu_utils
    integer, intent(in)  :: n
    integer, intent(out) :: errflg !return code. -1 for error, 0 for normal
    complex(dp), intent(in), dimension(n,n)  :: matrix  !input matrix
    complex(dp), intent(out), dimension(n,n) :: inverse !inverted matrix
    !local variables

    devptr_t                        :: devPtrA_new
    integer                         :: lda
    integer                         :: rc = RC_SUCCESS

    !start of the execution commands

    lda = n

#if defined (OPENMP) && !defined (CUDA)
    !TODO: here put the openMP code for the matinv
    write(*,*)"Need to have the matinv_Z for the case OPENMP and not CUDA" 
#elif defined (CUDA) && defined (MAGMA)
    ! now using the magma to get inverse matrix
    rc = cublas_alloc(n*n, size_of_double_complex, devPtrA_new)
    if (rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    rc = cublas_set_matrix (n, n, size_of_double_complex , matrix, n, &
         devPtrA_new, n )
    if ( rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    inverse = 0.0d0 !initializing the inverted matrix
    call getBlockInverseZMatrix_gpu(n, n, lda, size_of_double_complex, n, &
                                    devPtrA_new, inverse)

    rc = cublas_free(devPtrA_new)

    errflg = rc
#else
    call getBlockInverseZMatrix_cpu(n,n,lda,size_of_double_complex,n,matrix,&
                                    inverse)
!     call get_invert_gpu_utils_warning()
!     errflg = -1
!     stop
#endif
 
    return
  end subroutine matinv_Z_gpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to warns on the compilation derectives
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_invert_gpu_utils_warning()
    implicit none
     write(*,*)"**************************WARNING******************************"
     write(*,*)"You need to compile with -DCUDA and -DMAGMA                    "
     write(*,*)"to acces the CUDA environment computation using GPU            "
     write(*,*)"switching to the CPU version of matinv                         "
     write(*,*)"***************************************************************"
    return
  end subroutine get_invert_gpu_utils_warning

end module simple_math_gpu
!> \brief \b SGEMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at 
!            http://www.netlib.org/lapack/explore-html/ 
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
! 
!       .. Scalar Arguments ..
!       REAL ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),X(*),Y(*)
!       ..
!  
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!>
!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!>
!>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array of DIMENSION ( LDA, n ).
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array of DIMENSION at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!>           Before entry, the incremented array X must contain the
!>           vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array of DIMENSION at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!>           Before entry with BETA non-zero, the incremented array Y
!>           must contain the vector y. On exit, Y is overwritten by the
!>           updated vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee 
!> \author Univ. of California Berkeley 
!> \author Univ. of Colorado Denver 
!> \author NAG Ltd. 
!
!> \date November 2011
!
!> \ingroup single_blas_level2
!
!> \par Further Details:
!
!  =====================================================================
