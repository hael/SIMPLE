SUBROUTINE DGERFS_F95(A, AF, IPIV, B, X, TRANS, FERR, BERR, INFO)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: GERFS_F77 => LA_GERFS
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   INTEGER, INTENT(IN), TARGET :: IPIV(:)
   REAL(WP), INTENT(IN) :: A(:,:), AF(:,:), B(:,:)
   REAL(WP), INTENT(INOUT) :: X(:,:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: FERR(:), BERR(:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_GERFS improves the computed solution X of a system of linear 
! equations   A X = B  or  A^T X = B
! and provides error bounds and backward error estimates for
! the solution. LA_GERFS uses the LU factors computed by LA_GETRF.
!
! Arguments
! =========
! SUBROUTINE LA_GERFS (A, AF, IPIV, B, X, TRANS, FERR, BERR, INFO)
!    <type>(<wp>), INTENT(IN)  :: A(:,:), AF(:,:), <rhs>
!    INTEGER, INTENT(IN) :: IPIV(:)
!    <type>(<wp>), INTENT(INOUT) :: <sol>
!    REAL(<wp>), INTENT(OUT), OPTIONAL :: <err>
!    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
!    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!    <type> ::= REAL | COMPLEX
!    <wp>   ::= KIND(1.0) | KIND(1.0D0)
!    <rhs>  ::= B(:,:) | B(:)
!    <sol>  ::= X(:,:) | X(:)
!    <err>  ::= FERR(:), BERR(:) | FERR, BERR
!
! =====================
!
! A     (input) either REAL or COMPLEX square array,
!       shape (:,:), size(A,1) == size(A,2).
!       The original matrix A.
!
! AF    (input) either REAL or COMPLEX square array,
!       shape (:,:), size(AF,1) == size(AF,2) == size(A,1).
!       The factors L and U from the factorization A = PLU
!       as computed by LA_GETRF.
!
! IPIV  (input) INTEGER array, shape (:), size(IPIV) == size(A,1).
!       The pivot indices from LA_GETRF; for 1<=i<=size(A,1), row i
!       of the matrix was interchanged with row IPIV(i).
!
! B     (input) either REAL or COMPLEX rectangular array,
!       shape either (:,:) or (:), size(B,1) or size(B) == size(A,1).
!       The right hand side vector(s) of matrix B for
!       the system of equations AX = B.
!
! X     (input/output) either REAL or COMPLEX rectangular array,
!       shape either (:,:) or (:), size(X,1) or size(X) == size(A,1).
!       On entry, the solution matrix X, as computed by LA_GETRS.
!       On exit, the improved solution matrix X.
!
! TRANS Optional (input) CHARACTER*1
!       If TRANS is present, it specifies the form of the system 
!       of equations:
!          = 'N':  A X = B    (No transpose)
!          = 'T':  A^T X = B  (Transpose)
!          = 'C':  A^H X = B  (Conjugate transpose = Transpose)
!       otherwise TRANS = 'N' is assumed.
!
! FERR  Optional (output) either REAL array of shape (:) or REAL
!       scalar. If it is an array, size(FERR) == size(X,2).
!       The estimated forward error bound for each solution vector
!       X(j) (the j-th column of the solution matrix X).
!       If XTRUE is the true solution corresponding to X(j), FERR(j)
!       is an estimated upper bound for the magnitude of the largest
!       element in (X(j) - XTRUE) divided by the magnitude of the
!       largest element in X(j).  The estimate is as reliable as
!       the estimate for RCOND, and is almost always a slight
!       overestimate of the true error.
!
! BERR  Optional (output) either REAL array of shape (:) or REAL
!       scalar. If it is an array, size(BERR) == size(X,2).
!       The componentwise relative backward error of each solution
!       vector X(j) (i.e., the smallest relative change in
!       any element of A or B that makes X(j) an exact solution).
!
! INFO  Optional (output) INTEGER.
!       If INFO is present
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!       If INFO is not present and an error occurs, then the program is
!          terminated with an error message.
!
! Internal Parameters
! ===================
!
! ITMAX is the maximum number of steps of iterative refinement.
! It is set to 5 in the LAPACK77 subroutines
! -----------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GERFS'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LTRANS
   INTEGER :: LINFO, NRHS, N, ISTAT, ISTAT1, SFERR, SBERR
!  .. LOCAL ARRAYS, POINTERS ..
   INTEGER, POINTER :: IWORK(:)
   REAL(WP), POINTER :: LFERR(:), LBERR(:)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, SIZE, MAX
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A, 1); NRHS = SIZE(B, 2)
   IF( PRESENT(FERR) )THEN; SFERR = SIZE(FERR); ELSE; SFERR = NRHS; END IF
   IF( PRESENT(BERR) )THEN; SBERR = SIZE(BERR); ELSE; SBERR = NRHS; END IF
   IF(PRESENT(TRANS))THEN; LTRANS = TRANS; ELSE; LTRANS='N'
   END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE(A, 2) /= N .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE(AF, 1) /= N .OR. SIZE(AF, 2) /= N )THEN; LINFO = -2
   ELSE IF ( SIZE( IPIV ) /= N ) THEN; LINFO = -3
   ELSE IF ( SIZE(B, 1) /= N ) THEN; LINFO = -4
   ELSE IF ( SIZE(X, 1) /= N .OR. SIZE(X, 2) /= NRHS ) THEN; LINFO = -5
   ELSE IF( .NOT.( LSAME(LTRANS,'N') .OR.  LSAME(LTRANS,'T') .OR. &
                  LSAME(LTRANS,'C') ) )THEN; LINFO = -6
   ELSE IF( SFERR /= NRHS )THEN; LINFO = -7
   ELSE IF( SBERR /= NRHS )THEN; LINFO = -8
   ELSE IF ( N > 0 ) THEN
      IF( .NOT.PRESENT(FERR) ) THEN
         ALLOCATE( LFERR(NRHS), STAT=ISTAT )
      ELSE; LFERR => FERR; END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(BERR) ) THEN
            ALLOCATE( LBERR(NRHS), STAT=ISTAT )
         ELSE; LBERR => BERR; END IF
      END IF
      IF( ISTAT == 0 )THEN; ALLOCATE( WORK(3*N), IWORK(N), STAT=ISTAT ); END IF
      IF( ISTAT == 0 )THEN
!        .. CALL LAPACK77 ROUTINE
         CALL GERFS_F77( LTRANS, N, NRHS, A, MAX(1,N), AF, MAX(1,N), &
                         IPIV, B, MAX(1,N), X, MAX(1,N), LFERR, &
                         LBERR, WORK, IWORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(FERR) ) DEALLOCATE( LFERR, STAT=ISTAT1 )
      IF( .NOT.PRESENT(BERR) ) DEALLOCATE( LBERR, STAT=ISTAT1 )
      DEALLOCATE(WORK, IWORK, STAT=ISTAT1 )
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
END SUBROUTINE DGERFS_F95
