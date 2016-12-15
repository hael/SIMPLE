      SUBROUTINE ZGESV1_F95( A, B, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!     .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO
      USE F77_LAPACK, ONLY: GESV_F77 => LA_GESV
!     .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!     .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. ARRAY ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
      COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:)
!     .. PARAMETERS ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GESV'
!     .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV
!     .. LOCAL POINTERS ..
      INTEGER, POINTER :: LPIV(:)
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. EXECUTABLE STATEMENTS ..
      LINFO = 0; ISTAT = 0
      IF( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = SIZE(A,1)
      END IF
!     .. TEST THE ARGUMENTS
      IF( SIZE( A, 2 ) /= SIZE(A,1) .OR. SIZE(A,1) < 0 ) THEN; LINFO = -1
      ELSE IF( SIZE( B ) /= SIZE(A,1) ) THEN; LINFO = -2
      ELSE IF( SIPIV /= SIZE(A,1) )THEN; LINFO = -3
      ELSE IF ( SIZE(A,1) > 0 ) THEN
         IF( PRESENT(IPIV) )THEN; LPIV => IPIV
         ELSE; ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT ); END IF
         IF ( ISTAT == 0 ) THEN
!        .. CALL LAPACK77 ROUTINE ..
            CALL GESV_F77( SIZE(A,1), 1, A, MAX(1,SIZE(A,1)), LPIV, B, MAX(1,SIZE(A,1)),         &
                           LINFO )
         ELSE; LINFO = -100; END IF
         IF( .NOT.PRESENT(IPIV) ) DEALLOCATE(LPIV, STAT = ISTAT1 )
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
      END SUBROUTINE ZGESV1_F95
