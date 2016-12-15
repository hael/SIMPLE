SUBROUTINE DGETRS_F95( A, IPIV, B, TRANS, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: LSAME, ERINFO
   USE F77_LAPACK, ONLY: GETRS_F77 => LA_GETRS
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   INTEGER, INTENT(IN) :: IPIV(:)
   REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_GETRS solves a system of linear equations
!    A X = B, A^T X = B or  A^H X = B
! with a general square matrix A using the LU factorization computed
! by LA_GETRF.
!
! Arguments
! =========
! SUBROUTINE LA_GETRS (A, IPIV, B, TRANS, INFO)
!    <type>(<wp>), INTENT(IN)  :: A(:,:)
!    <type>(<wp>), INTENT(INOUT) :: <rhs>
!    INTEGER, INTENT(IN) :: IPIV(:)
!    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
!    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!    <type> ::= REAL | COMPLEX
!    <wp>   ::= KIND(1.0) | KIND(1.0D0)
!    <rhs>  ::= B(:,:) | B(:)
!
! =====================
!
! A     (input) either REAL or COMPLEX square array,
!       shape (:,:), size(A,1) == size(A,2).
!       The factors L and U from the factorization A = PLU as computed 
!       by LA_GETRF.
!
! IPIV  (input) INTEGER array, shape (:), size(IPIV) == size(A,1).
!       The pivot indices from LA_GETRF; for 1<=i<=size(A,1), row i
!       of the matrix was interchanged with row IPIV(i).
!
! B     (input/output) either REAL or COMPLEX rectangular array,
!       shape either (:,:) or (:), size(B,1) or size(B) == size(A,1).
!       On entry, the right hand side vector(s) of matrix B for
!          the system of equations AX = B.
!       On exit, if there is no error, the matrix of solution
!          vector(s) X.
!
! TRANS Optional (input) CHARACTER*1
!       If TRANS is present, it specifies the form of the system 
!          of equations:
!          = 'N':  A X = B    (No transpose)
!          = 'T':  A^T X = B  (Transpose)
!          = 'C':  A^H X = B  (Conjugate transpose = Transpose)
!       otherwise TRANS = 'N' is assumed.
!
! INFO  Optional (output) INTEGER.
!       If INFO is present
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!       If INFO is not present and an error occurs, then the program is
!          terminated with an error message.
!---------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GETRS'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LTRANS
   INTEGER    :: LINFO, NRHS, N, LD
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC SIZE, MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; N = SIZE(A, 1); NRHS = SIZE(B,2); LD = MAX(1,N)
   IF(PRESENT(TRANS))THEN; LTRANS = TRANS; ELSE; LTRANS='N'; END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 ) THEN; LINFO = -1
   ELSE IF( SIZE( IPIV ) /= N ) THEN; LINFO = -2
   ELSE IF( SIZE( B, 1 ) /= N ) THEN; LINFO = -3
   ELSE IF(.NOT.LSAME(LTRANS,'N') .AND. .NOT.LSAME(LTRANS,'T').AND. &
           .NOT.LSAME(LTRANS,'C'))THEN; LINFO = -4
   ELSE
!  .. CALL LAPACK77 ROUTINE
      CALL GETRS_F77( LTRANS, N, NRHS, A, LD, IPIV, B, LD, LINFO )
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO )
END SUBROUTINE DGETRS_F95
