      INTEGER FUNCTION LA_WS_GELS( VER, M, N, NRHS, TRANS )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
      USE LA_AUXMOD, ONLY: LSAME
      USE F77_LAPACK, ONLY: ILAENV_F77 => LA_ILAENV
!     .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!     .. SCALAR ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN) :: TRANS, VER
      INTEGER, INTENT(IN) :: M, N, NRHS
!     .. PARAMETERS ..
      CHARACTER(LEN=5), PARAMETER :: NAME1='GEQRF', NAME2='ORMQR', NAME3='ORMLQ', NAME4='GELQF'
!     .. LOCAL SCALARS ..
      INTEGER :: NB, MN
      LOGICAL :: TPSD
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX, MIN
!     .. EXECUTABLE STATEMENTS ..
      MN = MIN( M, N )
      IF( LSAME( TRANS, 'N' ) )THEN
        TPSD = .FALSE.
      ELSE
        TPSD = .TRUE.
      ENDIF
      IF( M.GE.N ) THEN
        NB = ILAENV_F77( 1, VER//NAME1, ' ', M, N, -1, -1 )
        IF( TPSD ) THEN
          NB = MAX( NB,ILAENV_F77( 1,VER//NAME2,'LN',M,NRHS,N,-1 ) )
        ELSE
          NB = MAX( NB, ILAENV_F77( 1,VER//NAME2,'LT',M,NRHS,N,-1 ) )
        END IF
      ELSE
        NB = ILAENV_F77( 1, VER//NAME4, ' ', M, N, -1, -1 )
        IF( TPSD ) THEN
          NB = MAX( NB, ILAENV_F77( 1,VER//NAME3,'LT',N,NRHS,M,-1 ) )
        ELSE
          NB = MAX( NB, ILAENV_F77( 1,VER//NAME3,'LN',N,NRHS,M,-1 ) )
        END IF
      END IF
      LA_WS_GELS = MN + MAX( M, N, NRHS )*NB
      END FUNCTION LA_WS_GELS
