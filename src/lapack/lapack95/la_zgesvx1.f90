      SUBROUTINE ZGESVX1_F95( A, B, X, AF, IPIV, FACT, TRANS, EQUED, R, C, FERR, BERR, RCOND,    &
                              RPVGRW, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!     .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: LSAME, ERINFO
      USE F77_LAPACK, ONLY: GESVX_F77 => LA_GESVX
!     .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!     .. SCALAR ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS, FACT
      CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: EQUED
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
      REAL(WP), INTENT(OUT), OPTIONAL :: RCOND, RPVGRW
!     .. ARRAY ARGUMENTS ..
      COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:)
      COMPLEX(WP), INTENT(OUT) :: X(:)
      INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: IPIV(:)
      REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: C(:), R(:)
      COMPLEX(WP), INTENT(INOUT), OPTIONAL, TARGET :: AF(:,:)
      REAL(WP), INTENT(OUT), OPTIONAL :: FERR, BERR
!     .. PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GESVX'
!     .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LFACT, LTRANS, LEQUED
      INTEGER :: ISTAT, ISTAT1, LD, LINFO, N, NRHS, S1AF, S2AF, SC, SIPIV, SR
      REAL(WP) :: LRCOND, MVR, MVC, LFERR, LBERR
!     .. LOCAL POINTERS ..
      INTEGER, POINTER :: LPIV(:)
      REAL(WP),  POINTER :: LC(:), LR(:), RWORK(:)
      COMPLEX(WP),  POINTER :: WORK(:), LAF(:, :)
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX, PRESENT, SIZE, MINVAL, TINY
!     .. EXECUTABLE STATEMENTS ..
      LINFO = 0; ISTAT = 0; N = SIZE(A, 1); NRHS = 1; LD = MAX(1,N)
      IF( PRESENT(RCOND) ) RCOND = 1.0_WP
      IF( PRESENT(RPVGRW) ) RPVGRW = 1.0_WP
      IF( PRESENT(FACT) )THEN; LFACT = FACT; ELSE; LFACT='N'; END IF
      IF( PRESENT(EQUED) .AND. LSAME(LFACT,'F') )THEN; LEQUED = EQUED
      ELSE; LEQUED='N'; END IF
      IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = N; END IF
      IF( PRESENT(AF) )THEN; S1AF = SIZE(AF,1); S2AF = SIZE(AF,2)
      ELSE; S1AF = N; S2AF = N; END IF
      IF( ( PRESENT(C) ) )THEN; SC = SIZE(C); ELSE; SC = N; END IF
      IF( ( PRESENT(C) .AND. LSAME(LFACT,'F') ) .AND.                   &
     &    ( LSAME(LEQUED,'C') .OR. LSAME(LEQUED,'B') ) )THEN; MVC = MINVAL(C)
      ELSE; MVC = TINY(1.0_WP); END IF
      IF( PRESENT(R) )THEN
         SR = SIZE(R)
      ELSE
         SR = N
      END IF
      IF( ( PRESENT(R) .AND. LSAME(LFACT,'F') ) .AND.                   &
     &    ( LSAME(LEQUED,'R') .OR. LSAME(LEQUED,'B') ) )THEN
         MVR = MINVAL(R)
      ELSE
         MVR = TINY(1.0_WP)
      END IF
      IF(PRESENT(TRANS))THEN
         LTRANS = TRANS
      ELSE
         LTRANS='N'
      END IF
!     .. TEST THE ARGUMENTS
      IF( SIZE(A, 2) /= N .OR. N < 0 )THEN
         LINFO = -1
      ELSE IF( SIZE(B) /= N )THEN
         LINFO = -2
      ELSE IF( SIZE(X) /= N )THEN
         LINFO = -3
      ELSE IF( S1AF /= N .OR. S2AF /= N ) THEN
         LINFO = -4
      ELSE IF( SIPIV /= N )THEN
         LINFO = -5
      ELSE IF( SR /= N .OR. MVR <= 0.0_WP )THEN
         LINFO = -9
      ELSE IF( SC /= N .OR. MVC <= 0.0_WP )THEN
         LINFO = -10
      ELSE IF( ( .NOT. ( LSAME(LFACT,'F') .OR. LSAME(LFACT,'N') .OR.    &
     &                 LSAME(LFACT,'E') ) ) .OR.                        &
     &    ( LSAME(LFACT,'F') .AND. .NOT.( PRESENT(AF) .AND.             &
     &      PRESENT(IPIV) ) ) )THEN
         LINFO = -6
      ELSE IF( .NOT.( LSAME(LTRANS,'N') .OR.  LSAME(LTRANS,'T') .OR.    &
     &               LSAME(LTRANS,'C') ) )THEN
         LINFO = -7
      ELSE IF( ( .NOT.( LSAME(LEQUED,'N') .OR. LSAME(LEQUED,'R') .OR.   &
     &       LSAME(LEQUED,'C') .OR. LSAME(LEQUED,'B') )                 &
     &           .AND. LSAME(LFACT,'F') ) .OR.                          &
     &      ( ( LSAME(LEQUED,'R') .OR. LSAME(LEQUED,'B') ) .AND.        &
     &           .NOT.PRESENT(R) ) .OR.                                 &
     &      ( ( LSAME(LEQUED,'C') .OR. LSAME(LEQUED,'B') ) .AND.        &
     &           .NOT.PRESENT(C) ) )THEN
         LINFO = -8
      ELSE IF ( N > 0 )THEN
         IF( .NOT.PRESENT(AF) ) THEN
            ALLOCATE( LAF(LD,N), STAT=ISTAT )
         ELSE
            LAF => AF
         END IF
         IF( .NOT.PRESENT(IPIV) )THEN
            ALLOCATE( LPIV(N), STAT=ISTAT )
         ELSE
            LPIV => IPIV
         END IF
         IF( .NOT.PRESENT(R) )THEN
            ALLOCATE( LR(N), STAT=ISTAT )
         ELSE
            LR => R
         END IF
         IF( .NOT.PRESENT(C) )THEN
            ALLOCATE( LC(N), STAT=ISTAT )
         ELSE
            LC => C
         END IF
            ALLOCATE(WORK(2*N), RWORK(2*N), STAT=ISTAT )
         IF( ISTAT == 0 )THEN
!           .. CALL LAPACK77 ROUTINE
            CALL GESVX_F77( LFACT, LTRANS, N, NRHS, A, LD, LAF, LD, LPIV, LEQUED, LR, LC, B, LD, &
                            X, LD, LRCOND, LFERR, LBERR, WORK, RWORK, LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT.PRESENT(R) ) DEALLOCATE( LR, STAT=ISTAT1 )
         IF( .NOT.PRESENT(C) ) DEALLOCATE( LC, STAT=ISTAT1 )
         IF( .NOT.PRESENT(AF) ) DEALLOCATE( LAF, STAT=ISTAT1 )
         IF( .NOT.PRESENT(IPIV) ) DEALLOCATE( LPIV, STAT=ISTAT1 )
         IF( PRESENT(FERR) ) FERR = LFERR
         IF( PRESENT(BERR) ) BERR = LBERR
         IF( PRESENT(EQUED) .AND. .NOT.LSAME(LFACT,'F') ) EQUED=LEQUED
         IF( PRESENT(RCOND) ) RCOND=LRCOND
         IF( PRESENT(RPVGRW) ) RPVGRW=RWORK(1)
         DEALLOCATE( WORK, RWORK, STAT=ISTAT1 )
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
      END SUBROUTINE ZGESVX1_F95
