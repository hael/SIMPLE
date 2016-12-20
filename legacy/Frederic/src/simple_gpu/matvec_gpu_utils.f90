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
module matvec_gpu_utils

  use simple_defs
  use simple_cuda_defs
  use simple_timing

  implicit none

#define size_t integer*8
#define devptr_t integer*8

  logical,external :: LSAME

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single precision 2d array of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_SMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    real(sp)                      :: alpha,beta

    real(sp)                      :: A(lda,*), x(*)
    real(sp)                      :: y(*)
    character                     :: trans

    !local variables
    integer :: rc = RC_SUCCESS
    !function calls
    integer :: smatvec_tgpu

#if defined (CUDA)
    rc= smatvec_tgpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    call get_matvec_gpu_utils_warning()
    rc = RC_FAIL
    stop
#endif
    
    return
  end subroutine get_SMatvec_gpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single precision 2d array of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_DMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    real(dp)                      :: alpha,beta

    real(dp)                      :: A(lda,*), x(*)
    real(dp)                      :: y(*)
    character                     :: trans

    !local variables
    integer :: rc = RC_SUCCESS
    !function calls
    integer :: dmatvec_tgpu

#if defined (CUDA)
    rc= dmatvec_tgpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    call get_matvec_gpu_utils_warning()
    rc = RC_FAIL
    stop
#endif
    
    return
  end subroutine get_DMatvec_gpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single precision 2d array of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_CMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    complex(sp)                      :: alpha,beta

    complex(sp)                      :: A(lda,*), x(*)
    complex(sp)                      :: y(*)
    character                     :: trans

    !local variables
    integer :: rc = RC_SUCCESS
    !function calls
    integer :: cmatvec_tgpu

#if defined (CUDA)
    rc= cmatvec_tgpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    call get_matvec_gpu_utils_warning()
    rc = RC_FAIL
    stop
#endif
    
    return
  end subroutine get_CMatvec_gpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single precision 2d array of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_ZMatvec_gpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
    implicit none
    integer                       :: n,m
    integer                       :: lda,incx,incy

    complex(dp)                      :: alpha,beta

    complex(dp)                      :: A(lda,*), x(*)
    complex(dp)                      :: y(*)
    character                     :: trans

    !local variables
    integer :: rc = RC_SUCCESS
    !function calls
    integer :: zmatvec_tgpu

#if defined (CUDA)
    rc= zmatvec_tgpu(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
#else
    call get_matvec_gpu_utils_warning()
    rc = RC_FAIL
    stop
#endif
    
    return
  end subroutine get_ZMatvec_gpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to warns on the compilation derectives
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_matvec_gpu_utils_warning()
    implicit none
     write(*,*)"**************************WARNING******************************"
     write(*,*)"You need to compile with -DCUDA                                "
     write(*,*)"to acces the CUDA environment computation using GPU            "
     write(*,*)"Try to recompile witht he proper compilation directives        "
     write(*,*)"***************************************************************"
     stop
    return
  end subroutine get_matvec_gpu_utils_warning

end module matvec_gpu_utils

  
