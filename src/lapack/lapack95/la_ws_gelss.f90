      INTEGER FUNCTION LA_WS_GELSS( VER, M, N, NRHS )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
      USE F77_LAPACK, ONLY: ILAENV_F77 => LA_ILAENV
!     .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!     .. SCALAR ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN) :: VER
      INTEGER, INTENT(IN) :: M, N, NRHS
!     .. PARAMETERS ..
      CHARACTER(LEN=5), PARAMETER :: NAME1='GELSS', NAME2='GEQRF', NAME3='ORMQR', NAME4='GEBRD', &
                                     NAME5='ORMBR', NAME6='ORGBR', NAME7='GELQF', NAME8='ORMLQ'
!     .. LOCAL SCALARS ..
      INTEGER :: MNTHR, MINWRK, MAXWRK, MM, BDSPAC
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX
!     .. EXECUTABLE STATEMENTS ..
      MNTHR = ILAENV_F77( 6, VER//NAME1, ' ', M, N, NRHS, -1 )
!
!     COMPUTE WORKSPACE
!      (NOTE: COMMENTS IN THE CODE BEGINNING "Workspace:" DESCRIBE THE
!       MINIMAL AMOUNT OF WORKSPACE NEEDED AT THAT POINT IN THE CODE,
!       AS WELL AS THE PREFERRED AMOUNT FOR GOOD PERFORMANCE.
!       NB REFERS TO THE OPTIMAL BLOCK SIZE FOR THE IMMEDIATELY
!       FOLLOWING SUBROUTINE, AS RETURNED BY ILAENV.)
!
      MINWRK = 1
      MAXWRK = 0
      MM = M
      IF( M.GE.N .AND. M.GE.MNTHR ) THEN
!
!        PATH 1A - OVERDETERMINED, WITH MANY MORE ROWS THAN COLUMNS
!
         MM = N
         MAXWRK = MAX(MAXWRK,N+N*ILAENV_F77(1,VER//NAME2,' ',M,N,-1,-1))
         MAXWRK = MAX( MAXWRK, N+NRHS*                                  &
     &            ILAENV_F77( 1, VER//NAME3, 'LT', M, NRHS, N, -1 ) )
      END IF
      IF( M.GE.N ) THEN
!
!        PATH 1 - OVERDETERMINED OR EXACTLY DETERMINED
!
!        COMPUTE WORKSPACE NEEDE FOR DBDSQR
!
         BDSPAC = MAX( 1, 5*N-4 )
         MAXWRK = MAX( MAXWRK, 3*N+( MM+N )*                            &
     &            ILAENV_F77( 1, VER//NAME4, ' ', MM, N, -1, -1 ) )
         MAXWRK = MAX( MAXWRK, 3*N+NRHS*                                &
     &            ILAENV_F77( 1, VER//NAME5, 'QLT', MM, NRHS, N, -1 ) )
         MAXWRK = MAX( MAXWRK, 3*N+( N-1 )*                             &
     &            ILAENV_F77( 1, VER//NAME6, 'P', N, N, N, -1 ) )
         MAXWRK = MAX( MAXWRK, BDSPAC )
         MAXWRK = MAX( MAXWRK, N*NRHS )
         MINWRK = MAX( 3*N+MM, 3*N+NRHS, BDSPAC )
         MAXWRK = MAX( MINWRK, MAXWRK )
      END IF
      IF( N.GT.M ) THEN
!
!        COMPUTE WORKSPACE NEEDE FOR DBDSQR
!
         BDSPAC = MAX( 1, 5*M-4 )
         MINWRK = MAX( 3*M+NRHS, 3*M+N, BDSPAC )
         IF( N.GE.MNTHR ) THEN
!
!           PATH 2A - UNDERDETERMINED, WITH MANY MORE COLUMNS
!           THAN ROWS
!
            MAXWRK = M + M*ILAENV_F77( 1,VER//NAME7,' ',M,N,-1,-1 )
            MAXWRK = MAX( MAXWRK, M*M+4*M+2*M*                          &
     &               ILAENV_F77( 1, VER//NAME4, ' ', M, M, -1, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS*                         &
     &               ILAENV_F77( 1,VER//NAME5,'QLT',M,NRHS,M,-1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )*                      &
     &               ILAENV_F77( 1, VER//NAME6, 'P', M, M, M, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+M+BDSPAC )
            IF( NRHS.GT.1 ) THEN
               MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
            ELSE
               MAXWRK = MAX( MAXWRK, M*M+2*M )
            END IF
            MAXWRK = MAX( MAXWRK, M+NRHS*                               &
     &               ILAENV_F77( 1, VER//NAME8, 'LT', N, NRHS, M, -1 ) )
         ELSE
!
!           PATH 2 - UNDERDETERMINED
!
            MAXWRK = 3*M+(N+M)*ILAENV_F77(1,VER//NAME4,' ',M,N,-1,-1)
            MAXWRK = MAX( MAXWRK, 3*M+NRHS*                             &
     &               ILAENV_F77( 1,VER//NAME5,'QLT',M,NRHS,M,-1 ) )
            MAXWRK = MAX( MAXWRK, 3*M+M*                                &
     &               ILAENV_F77( 1, VER//NAME6, 'P', M, N, M, -1 ) )
            MAXWRK = MAX( MAXWRK, BDSPAC )
            MAXWRK = MAX( MAXWRK, N*NRHS )
         END IF
      END IF
      LA_WS_GELSS = MAX( MINWRK, MAXWRK )
      END FUNCTION LA_WS_GELSS
