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
module matmul_cpu_utils

  use simple_defs
  use simple_timing

  implicit none

  logical,external :: LSAME

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single precision 2d array of size
! (n)x(m) on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_SMatmul_cpu(transa,transb,m,n,k, &
                             alpha,               &
                             A,lda,               &
                             B,ldb,               &
                             beta,                &
                             prod_cpu,ldc         )

    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    real(sp)                      :: alpha,beta
    real(sp)                      :: A(lda,*), B(ldb,*)
    real(sp)                      :: prod_cpu(ldc,*)
    character                     :: transa,transb

    call sgemm(transa,transb,m,n,k,&
               alpha,              &
               A,lda,              &
               B,ldb,              &
               beta,               &
               prod_cpu,ldc)
    return
  end subroutine get_SMatmul_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Double precision 2d array of size
! (n)x(m) on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_DMatmul_cpu(transa,transb,m,n,k, &
                             alpha,               &
                             A,lda,               &
                             B,ldb,               &
                             beta,                &
                             prod_cpu,ldc         )

    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    real(dp)                      :: alpha,beta
    real(dp)                      :: A(lda,*), B(ldb,*)
    real(dp)                      :: prod_cpu(ldc,*)
    character                     :: transa,transb

    call dgemm(transa,transb,m,n,k,&
               alpha,              &
               A,lda,              &
               B,ldb,              &
               beta,               &
               prod_cpu,ldc)
    return
  end subroutine get_DMatmul_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Single complex precision 2d array
! of size
! (n)x(m) on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_CMatmul_cpu(transa,transb,m,n,k, &
                             alpha,               &
                             A,lda,               &
                             B,ldb,               &
                             beta,                &
                             prod_cpu,ldc         )

    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    complex(sp)                   :: alpha,beta
    complex(sp)                   :: A(lda,*), B(ldb,*)
    complex(sp)                   :: prod_cpu(ldc,*)
    character                     :: transa,transb

    call cgemm(transa,transb,m,n,k,&
               alpha,              &
               A,lda,              &
               B,ldb,              &
               beta,               &
               prod_cpu,ldc)
    return
  end subroutine get_CMatmul_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the matrix product of a Double complex precision 2d array
! of size
! (n)x(m) on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine get_ZMatmul_cpu(transa,transb,m,n,k, &
                             alpha,               &
                             A,lda,               &
                             B,ldb,               &
                             beta,                &
                             prod_cpu,ldc         )

    implicit none
    integer                       :: n,m,k
    integer                       :: lda,ldb,ldc
    complex(dp)                   :: alpha,beta
    complex(dp)                   :: A(lda,*), B(ldb,*)
    complex(dp)                   :: prod_cpu(ldc,*)
    character                     :: transa,transb

    call zgemm(transa,transb,m,n,k,&
               alpha,              &
               A,lda,              &
               B,ldb,              &
               beta,               &
               prod_cpu,ldc)
    
    return
  end subroutine get_ZMatmul_cpu

end module matmul_cpu_utils

  
