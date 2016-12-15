SUBROUTINE ZGETRI_F95( A, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO
   USE F77_LAPACK, ONLY: GETRI_F77 => LA_GETRI, ILAENV_F77 => ILAENV
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   INTEGER, INTENT(IN) :: IPIV(:)
   COMPLEX(WP), INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_GETRI computes the inverse of a matrix using the LU factorization 
! computed by LA_GETRF.
!
! Arguments
! =========
! SUBROUTINE LA_GETRI (A, IPIV, INFO)
!    <type>(<wp>), INTENT(INOUT)  :: A(:,:)
!    INTEGER, INTENT(IN) :: IPIV(:)
!    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!    <type> ::= REAL | COMPLEX
!    <wp>   ::= KIND(1.0) | KIND(1.0D0)
!
! =====================
!
! A      (input/output) either REAL or COMPLEX square array, shape (:,:),
!        size(A,1) == size(A,2).
!        On entry contains the factors L and U from the factorization
!           A = PLU as computed by LA_GETRF.
!        On exit, if INFO = 0, the inverse of the original matrix A.
!
! IPIV   (input) INTEGER array, shape (:), size(IPIV) == size(A,1).
!        The pivot indices from LA_GETRF; for 1<=i<=size(A,1), row i of
!        the matrix was interchanged with row IPIV(i).
!
! INFO   Optional (output) INTEGER.
!        If INFO is present
!           = 0: successful exit
!           < 0: if INFO = -k, the k-th argument had an illegal value
!           > 0: if INFO = k, U(k,k) is exactly zero.  The matrix is
!               singular and its inverse could not be computed.
!        If INFO is not present and an error occurs, then the program is
!           terminated with an error message.
!-----------------------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GETRI'
   CHARACTER(LEN=6), PARAMETER :: BSNAME = 'ZGETRI'
!  .. LOCAL SCALARS ..
   INTEGER    :: LINFO, N, LD, LWORK, ISTAT, ISTAT1, NB
!  .. LOCAL ARRAY ..
   COMPLEX(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC SIZE, MAX
!  .. EXECUTABLE STATEMENTS ..
   N = SIZE(A,1); LINFO = 0; LD = MAX(1,N); ISTAT = 0
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE( IPIV ) /= N )THEN; LINFO = -2
   ELSE IF( N > 0 )THEN
!     DETERMINE THE WORK SPACE.
      NB = ILAENV_F77( 1, BSNAME, ' ', N, -1, -1, -1 )
      IF( NB < 1 .OR. NB >= N )THEN; NB = 1; END IF
      LWORK = MAX( N*NB, 1 )
      ALLOCATE(WORK(LWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK, STAT=ISTAT1)
         LWORK = MAX(1,N); ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT == 0 ) CALL ERINFO( -200, SRNAME, LINFO )
      END IF
      IF( LINFO == 0 )THEN
         CALL GETRI_F77( N, A, LD, IPIV, WORK, LWORK, LINFO )
      ELSE; LINFO = -100; END IF
      DEALLOCATE(WORK, STAT=ISTAT1)
   END IF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZGETRI_F95
