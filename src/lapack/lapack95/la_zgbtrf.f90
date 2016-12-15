      SUBROUTINE ZGBTRF_F95( A, K, M, IPIV, RCOND, NORM, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: LSAME, ERINFO
      USE F77_LAPACK, ONLY: GBTRF_F77 => LA_GBTRF, LANGB_F77 => LA_LANGB, &
                                          GBCON_F77 => LA_GBCON
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
      INTEGER, INTENT(IN), OPTIONAL :: K, M
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
      REAL(WP), INTENT( OUT ), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
      INTEGER, INTENT( OUT ), OPTIONAL, TARGET :: IPIV( : )
      COMPLEX(WP), INTENT( INOUT ) :: A( :, : )
!  .. PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GBTRF'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LNORM
      INTEGER :: LINFO, N, LD, ISTAT, ISTAT1, MINMN, LWORK, SIPIV, &
                 LK, KU, LM
      REAL(WP) :: LANORM
!  .. LOCAL POINTERS ..
      INTEGER, POINTER :: LIPIV(:)
      REAL(WP), POINTER :: RWORK(:)
      COMPLEX(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC PRESENT, MIN, SIZE, TINY
!  .. EXECUTABLE STATEMENTS ..
      LD = SIZE(A,1); N = SIZE(A,2); LINFO = 0; ISTAT = 0
      IF( PRESENT(K) ) THEN; LK = K; ELSE; LK = (LD-1)/3; ENDIF
      IF( PRESENT(M) ) THEN; LM = M; ELSE; LM = N; ENDIF; MINMN = MIN(LM,N)
      IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = MINMN; ENDIF
      IF ( PRESENT(NORM) ) THEN; LNORM = NORM; ELSE; LNORM = '1'; END IF
!  .. TEST THE ARGUMENTS
      IF( N < 0 .OR. LD < 0 )THEN; LINFO = -1
      ELSE IF( LD - 2*LK -1 < 0 .OR. LK < 0 ) THEN; LINFO = -2
      ELSE IF( LM < 0 )THEN; LINFO = -3
      ELSE IF( SIPIV /= MINMN )THEN; LINFO = -4
      ELSE IF( PRESENT(RCOND) .AND. M /= N )THEN; LINFO = -5
      ELSE IF( ( .NOT.PRESENT(RCOND) .AND. PRESENT(NORM) ) .OR. &
               ( .NOT.LSAME(LNORM,'I') .AND. .NOT.LSAME(LNORM,'O') &
                                       .AND. LNORM /= '1' ) ) THEN; LINFO = -6
      ELSE IF( M > 0 .AND. N > 0 ) THEN
         IF( PRESENT(RCOND) .AND. M == N ) THEN
!        .. COMPUTE THE NORM OF THE MATRIX A
            IF( LNORM == 'I' ) THEN; LWORK = MINMN; ELSE; LWORK = 1; END IF
            ALLOCATE( RWORK(LWORK), STAT=ISTAT )
            IF( ISTAT == 0 )THEN
               LANORM = LANGB_F77( LNORM, MINMN, LK, KU, A, LD, RWORK )
            ELSE
               LINFO = -100
            END IF
            DEALLOCATE(RWORK, STAT=ISTAT1)
         END IF
         IF( LINFO == 0 ) THEN
            IF( PRESENT(IPIV) )THEN; LIPIV => IPIV
            ELSE; ALLOCATE( LIPIV( MINMN ), STAT=ISTAT ); ENDIF
            IF( ISTAT /= 0 )LINFO = -100
         END IF
         IF( LINFO == 0 ) THEN
!           .. COMPUTE THE LU FACTORS OF THE MATRIX A
            CALL GBTRF_F77( LM, N, LK, KU, A, LD, LIPIV, LINFO )
            IF( .NOT. PRESENT(IPIV) )DEALLOCATE( LIPIV, STAT=ISTAT )
            IF( PRESENT(RCOND) ) THEN
!              .. COMPUTE THE RECIPROCAL OF THE CONDITION NUMBER OF A
               IF( LANORM <= TINY(1.0_WP) .OR. M /= N .OR. LINFO /= 0 ) THEN
                  RCOND = 0.0_WP
               ELSE 
                  ALLOCATE(WORK(2*MINMN), RWORK(MINMN), STAT=ISTAT)
                  IF( ISTAT == 0 )THEN
                     CALL GBCON_F77( LNORM, MINMN, LK, KU, A, LD, LIPIV, &
                                     LANORM, RCOND, WORK, RWORK, LINFO )
                  ELSE
                     LINFO = -100
                  END IF
                  DEALLOCATE( WORK, RWORK, STAT=ISTAT1 )
               END IF
            END IF
         END IF
      ELSE IF( PRESENT(RCOND) ) THEN
         IF( M == N )THEN
            RCOND = 1.0_WP
         ELSE
            RCOND = 0.0_WP
         END IF
      END IF
      CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
      END SUBROUTINE ZGBTRF_F95
