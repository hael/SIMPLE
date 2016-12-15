!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
   INTEGER, INTENT(IN), OPTIONAL :: ITYPE
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
   REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!        LA_SYGV, LA_SYGVD, LA_HEGV and LA_HEGVD compute all eigenvalues
! and, optionally, all eigenvectors of generalized eigenvalue problems of
! the form A*z = lambda*B*z, A*B*z = lambda*z, and B*A*z = lambda*z,
! where A and B are real symmetric in the cases of LA_SYGV and LA_SYGVD 
! and complex Hermitian in the cases of LA_HEGV and LA_HEGVD. In all four
! cases B is positive deffinite.
!        LA_SYGVD and LA_HEGVD use a divide and conquer algorithm. If 
! eigenvectors are desired, they can be much faster than LA_SYGV and
! LA_HEGV for large matrices but use more workspace.
! 
! =========
! 
!       SUBROUTINE LA_SYGV / LA_SYGVD / LA_HEGV / LA_HEGVD( A, B, &
!                  W, ITYPE=itype, JOBZ=jobz, UPLO=uplo, INFO=info )
!            <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
!            REAL(<wp>), INTENT(OUT) :: W(:)
!            INTEGER, INTENT(IN), OPTIONAL :: ITYPE
!            CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!            INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!            <type> ::= REAL | COMPLEX
!            <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX square array, shape (:,:).
!         On entry, the matrix A.
!         If UPLO = 'U', the upper triangular part of A contains the 
!         upper triangular part of matrix A. If UPLO = 'L', the lower 
!         triangular part of A contains the lower triangular part of 
!         matrix A.
!         On exit, if JOBZ = 'V', then the columns of A contain the 
!         eigenvectors, normalized as follows:
!            if ITYPE = 1 or 2: Z^H*B*Z = I ,
!            if ITYPE = 3: Z^H*B^-1*Z = I .
!         If JOBZ = 'N', then the upper triangle (if UPLO = 'U') or the
!         lower triangle (if UPLO = 'L') of A, including the diagonal, 
!         is destroyed. 
! B       (input/output) REAL or COMPLEX square array, shape (:,:) with 
!         size(B,1) = size(A,1).
!         On entry, the matrix B. If UPLO = 'U', the upper triangular 
!         part of B contains the upper triangular part of matrix B. If 
! 	  UPLO = 'L', the lower triangular part of B contains the lower
!         triangular part of matrix B. 
!         On exit, if the part of B containing the matrix is overwritten 
!         by the triangular factor U or L of the Cholesky factorization 
!         B = U^H*U or B = L*L^H , respectively.
! W       (output) REAL array, shape (:) with size(W) = size(A,1).
!         The eigenvalues in ascending order.
! ITYPE   Optional (input) INTEGER.
!         Specifies the problem type to be solved:
!           = 1: A*z = lambda*B*z
!           = 2: A*B*z = lambda*z
!           = 3: B*A*z = lambda*z
!         Default value: 1.
! JOBZ    Optional (input) CHARACTER(LEN=1).
!           = 'N': Compute eigenvalues only;
!           = 'V': Compute eigenvalues and eigenvectors.
!         Default value: 'N'.
! UPLO    Optional (input) CHARACTER(LEN=1).
!           = 'U': Upper triangles of A and B are stored;
!           = 'L': Lower triangles of A and B are stored.
!         Default value: 'U'.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: the algorithm failed to converge or matrix B is not 
! 	  positive deffinite:
!             <= n: if INFO = i, i off-diagonal elements of an 
!                   intermediate tridiagonal form did not converge to 
!                   zero.
!             > n: if INFO = n+i, for 1 <= i <= n, then the leading minor
!                   of order i of B is not positive deffinite. The 
!                   factorization of B could not be completed and no
!                   eigenvalues or eigenvectors were computed.
!             n is the order of A.
!         If INFO is not present and an error occurs, then the program is 
!         terminated with an error message.
!------------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBZ, LUPLO
!  .. LOCAL ARRAYS ..
   COMPLEX(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC SIZE, MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; N = SIZE(A,1); LD = MAX(1,N); ISTAT = 0
   IF( PRESENT(ITYPE) )THEN
      LITYPE = ITYPE
   ELSE
      LITYPE = 1
   END IF
   IF( PRESENT(JOBZ) ) THEN
      LJOBZ = JOBZ
   ELSE
      LJOBZ = 'N'
   END IF
   IF( PRESENT(UPLO) ) THEN
      LUPLO = UPLO
   ELSE
      LUPLO = 'U'
   END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN
      LINFO = -1
   ELSE IF( SIZE( B, 1 ) /= N .OR. SIZE( B, 2 ) /= N  )THEN
      LINFO = -2
   ELSE IF( SIZE( W ) /= N )THEN
      LINFO = -3
   ELSE IF( LITYPE < 1 .OR. LITYPE > 3 )THEN
      LINFO = -4
   ELSE IF( .NOT.LSAME(LJOBZ,'N') .AND. .NOT.LSAME(LJOBZ,'V') )THEN
      LINFO = -5
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN
      LINFO = -6
   ELSE IF( N > 0 )THEN
!     .. DETERMINE THE WORKSPACE
      NB = ILAENV_F77( 1, BSNAME, LUPLO, N, -1, -1, -1 )
      IF( NB <= 1 .OR. NB >= N )THEN
         NB = 1
      ENDIF
