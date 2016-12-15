      FUNCTION ZLANGE1_F95( A, NORM, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!     .. "Use Statements" ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: LANGE_F77 => LA_LANGE
!     .. "Implicit Statement" ..
      IMPLICIT NONE
!     .. "Function Type" ..
      REAL(WP) :: ZLANGE1_F95
!     .. "Scalar Arguments" ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      COMPLEX(WP), INTENT(IN) :: A(:)
!     .. "Parameters" ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_LANGE'
!     .. "Local Scalars" ..
      CHARACTER(LEN=1) :: LNORM
      INTEGER :: ISTAT, ISTAT1, LDA, LINFO, M, N
      REAL(WP), TARGET :: LLWORK(1)
!     .. "Local Pointers" ..
      REAL(WP), POINTER :: WORK(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. "Executable Statements" ..
      LINFO = 0; M = SIZE(A); N = 1; LDA = MAX(1,M)
      ISTAT = 0
      IF( PRESENT(NORM) )THEN; LNORM = NORM; ELSE; LNORM = '1'; ENDIF
!     .. "Testing The Arguments" ..
      IF( M < 0 )THEN; LINFO = -1
      ELSE IF( .NOT. ( LSAME(LNORM,'M') .OR. LSAME(LNORM,'1') .OR. &
         LSAME(LNORM,'I') .OR. LSAME(LNORM,'F') .OR. &
         LSAME(LNORM,'E') ) )THEN; LINFO = -2
      ELSE
         IF( LSAME(LNORM,'I') )THEN; ALLOCATE( WORK( M), STAT=ISTAT )
         ELSE; WORK => LLWORK; ENDIF
         IF( ISTAT == 0 ) &
           ZLANGE1_F95 = LANGE_F77( LNORM, M, N, A, LDA, WORK )
         IF( LSAME(LNORM,'I') )DEALLOCATE( WORK, STAT=ISTAT1 )
      ENDIF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
      END FUNCTION ZLANGE1_F95
