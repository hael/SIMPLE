! subroutines imported from blas, lapack from netlib.org, manually made consistent with fortran 90
module simple_lapackblas
implicit none
public :: SGEEV

private

contains

    SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        REAL SA
        INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
        REAL SX(*),SY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MOD
        !     ..
        IF (N.LE.0) RETURN
        IF (SA.EQ.0.0) RETURN
        IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
            !
            !        code for both increments equal to 1
            !
            !
            !        clean-up loop
            !
            M = MOD(N,4)
            IF (M.NE.0) THEN
                DO I = 1,M
                    SY(I) = SY(I) + SA*SX(I)
                END DO
            END IF
            IF (N.LT.4) RETURN
            MP1 = M + 1
            DO I = MP1,N,4
                SY(I) = SY(I) + SA*SX(I)
                SY(I+1) = SY(I+1) + SA*SX(I+1)
                SY(I+2) = SY(I+2) + SA*SX(I+2)
                SY(I+3) = SY(I+3) + SA*SX(I+3)
            END DO
        ELSE
            !
            !        code for unequal increments or equal increments
            !          not equal to 1
            !
            IX = 1
            IY = 1
            IF (INCX.LT.0) IX = (-N+1)*INCX + 1
            IF (INCY.LT.0) IY = (-N+1)*INCY + 1
            DO I = 1,N
                SY(IY) = SY(IY) + SA*SX(IX)
                IX = IX + INCX
                IY = IY + INCY
            END DO
        END IF
        RETURN
    END SUBROUTINE SAXPY

    SUBROUTINE SCOPY(N,SX,INCX,SY,INCY)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
        REAL SX(*),SY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MOD
        !     ..
        IF (N.LE.0) RETURN
        IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
            !
            !        code for both increments equal to 1
            !
            !
            !        clean-up loop
            !
            M = MOD(N,7)
            IF (M.NE.0) THEN
                DO I = 1,M
                    SY(I) = SX(I)
                END DO
                IF (N.LT.7) RETURN
            END IF
            MP1 = M + 1
            DO I = MP1,N,7
                SY(I) = SX(I)
                SY(I+1) = SX(I+1)
                SY(I+2) = SX(I+2)
                SY(I+3) = SX(I+3)
                SY(I+4) = SX(I+4)
                SY(I+5) = SX(I+5)
                SY(I+6) = SX(I+6)
            END DO
        ELSE
            !
            !        code for unequal increments or equal increments
            !          not equal to 1
            !
            IX = 1
            IY = 1
            IF (INCX.LT.0) IX = (-N+1)*INCX + 1
            IF (INCY.LT.0) IY = (-N+1)*INCY + 1
            DO I = 1,N
                SY(IY) = SX(IX)
                IX = IX + INCX
                IY = IY + INCY
            END DO
        END IF
        RETURN
    END SUBROUTINE SCOPY

    SUBROUTINE SGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, &
        INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          JOB, SIDE
        INTEGER            IHI, ILO, INFO, LDV, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               V( LDV, * ), SCALE( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE
        PARAMETER          ( ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            LEFTV, RIGHTV
        INTEGER            I, II, K
        REAL               S
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Decode and Test the input parameters
        !
        RIGHTV = LSAME( SIDE, 'R' )
        LEFTV = LSAME( SIDE, 'L' )
        !
        INFO = 0
        IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
            .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
            INFO = -4
        ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
            INFO = -5
        ELSE IF( M.LT.0 ) THEN
            INFO = -7
        ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
            INFO = -9
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SGEBAK', -INFO )
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.EQ.0 ) &
            RETURN
        IF( M.EQ.0 ) &
            RETURN
        IF( LSAME( JOB, 'N' ) ) &
            RETURN
        !
        IF( ILO.EQ.IHI ) &
            GO TO 30
        !
        !     Backward balance
        !
        IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
            !
            IF( RIGHTV ) THEN
                DO  I = ILO, IHI
                    S = SCALE( I )
                    CALL SSCAL( M, S, V( I, 1 ), LDV )
                END DO
            END IF
            !
            IF( LEFTV ) THEN
                DO  I = ILO, IHI
                    S = ONE / SCALE( I )
                    CALL SSCAL( M, S, V( I, 1 ), LDV )
                END DO
            END IF
            !
        END IF
        !
        !     Backward permutation
        !
        !     For  I = ILO-1 step -1 until 1,
        !              IHI+1 step 1 until N do --
        !
30      CONTINUE
        IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
            IF( RIGHTV ) THEN
                DO  II = 1, N
                    I = II
                    IF( I.GE.ILO .AND. I.LE.IHI ) &
                        GO TO 40
                    IF( I.LT.ILO ) &
                        I = ILO - II
                    K = SCALE( I )
                    IF( K.EQ.I ) &
                        GO TO 40
                    CALL SSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
40              END DO
            END IF
            !
            IF( LEFTV ) THEN
                DO  II = 1, N
                    I = II
                    IF( I.GE.ILO .AND. I.LE.IHI ) &
                        GO TO 50
                    IF( I.LT.ILO ) &
                        I = ILO - II
                    K = SCALE( I )
                    IF( K.EQ.I ) &
                        GO TO 50
                    CALL SSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
50              END DO
            END IF
        END IF
        !
        RETURN
        !
        !     End of SGEBAK
        !
    END SUBROUTINE SGEBAK

    SUBROUTINE SGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          JOB
        INTEGER            IHI, ILO, INFO, LDA, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), SCALE( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
        REAL               SCLFAC
        PARAMETER          ( SCLFAC = 2.0E+0 )
        REAL               FACTOR
        PARAMETER          ( FACTOR = 0.95E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            NOCONV
        INTEGER            I, ICA, IEXC, IRA, J, K, L, M
        REAL               C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1, &
            SFMIN2
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX, MIN
        !
        !     Test the input parameters
        !
        INFO = 0
        IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND. &
            .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
            INFO = -1
        ELSE IF( N.LT.0 ) THEN
            INFO = -2
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -4
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SGEBAL', -INFO )
            RETURN
        END IF
        !
        K = 1
        L = N
        !
        IF( N.EQ.0 ) &
            GO TO 210
        !
        IF( LSAME( JOB, 'N' ) ) THEN
            DO  I = 1, N
                SCALE( I ) = ONE
            END DO
            GO TO 210
        END IF
        !
        IF( LSAME( JOB, 'S' ) ) &
            GO TO 120
        !
        !     Permutation to isolate eigenvalues if possible
        !
        GO TO 50
        !
        !     Row and column exchange.
        !
20      CONTINUE
        SCALE( M ) = J
        IF( J.EQ.M ) &
            GO TO 30
        !
        CALL SSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
        CALL SSWAP( N-K+1, A( J, K ), LDA, A( M, K ), LDA )
        !
30      CONTINUE
        !!! GO TO ( 40, 80 )IEXC -> obsolecent feature; changed manually
        if ( IEXC == 1) THEN
            GO TO 40
        ELSE
            GO TO 80
        END if
        !
        !     Search for rows isolating an eigenvalue and push them down.
        !
40      CONTINUE
        IF( L.EQ.1 ) &
            GO TO 210
        L = L - 1
        !
50      CONTINUE
        DO  J = L, 1, -1
            !
            DO  I = 1, L
                IF( I.EQ.J ) &
                    GO TO 60
                IF( A( J, I ).NE.ZERO ) &
                    GO TO 70
60          END DO
            !
            M = L
            IEXC = 1
            GO TO 20
70      END DO
        !
        GO TO 90
        !
        !     Search for columns isolating an eigenvalue and push them left.
        !
80      CONTINUE
        K = K + 1
        !
90      CONTINUE
        DO  J = K, L
            !
            DO  I = K, L
                IF( I.EQ.J ) &
                    GO TO 100
                IF( A( I, J ).NE.ZERO ) &
                    GO TO 110
100         END DO
            !
            M = K
            IEXC = 2
            GO TO 20
110     END DO
        !
120     CONTINUE
        DO  I = K, L
            SCALE( I ) = ONE
        END DO
        !
        IF( LSAME( JOB, 'P' ) ) &
            GO TO 210
        !
        !     Balance the submatrix in rows K to L.
        !
        !     Iterative loop for norm reduction
        !
        SFMIN1 = SLAMCH( 'S' ) / SLAMCH( 'P' )
        SFMAX1 = ONE / SFMIN1
        SFMIN2 = SFMIN1*SCLFAC
        SFMAX2 = ONE / SFMIN2
140     CONTINUE
        NOCONV = .FALSE.
        !
        DO  I = K, L
            !
            C = SNRM2( L-K+1, A( K, I ), 1 )
            R = SNRM2( L-K+1, A( I, K ), LDA )
            ICA = ISAMAX( L, A( 1, I ), 1 )
            CA = ABS( A( ICA, I ) )
            IRA = ISAMAX( N-K+1, A( I, K ), LDA )
            RA = ABS( A( I, IRA+K-1 ) )
            !
            !        Guard against zero C or R due to underflow.
            !
            IF( C.EQ.ZERO .OR. R.EQ.ZERO ) &
                GO TO 200
            G = R / SCLFAC
            F = ONE
            S = C + R
160         CONTINUE
            IF( C.GE.G .OR. MAX( F, C, CA ).GE.SFMAX2 .OR. &
                MIN( R, G, RA ).LE.SFMIN2 )GO TO 170
            F = F*SCLFAC
            C = C*SCLFAC
            CA = CA*SCLFAC
            R = R / SCLFAC
            G = G / SCLFAC
            RA = RA / SCLFAC
            GO TO 160
            !
170         CONTINUE
            G = C / SCLFAC
180         CONTINUE
            IF( G.LT.R .OR. MAX( R, RA ).GE.SFMAX2 .OR. &
                MIN( F, C, G, CA ).LE.SFMIN2 )GO TO 190
            IF( SISNAN( C+F+CA+R+G+RA ) ) THEN
                !
                !           Exit if NaN to avoid infinite loop
                !
                INFO = -3
                CALL XERBLA( 'SGEBAL', -INFO )
                RETURN
            END IF
            F = F / SCLFAC
            C = C / SCLFAC
            G = G / SCLFAC
            CA = CA / SCLFAC
            R = R*SCLFAC
            RA = RA*SCLFAC
            GO TO 180
            !
            !        Now balance.
            !
190         CONTINUE
            IF( ( C+R ).GE.FACTOR*S ) &
                GO TO 200
            IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
                IF( F*SCALE( I ).LE.SFMIN1 ) &
                    GO TO 200
            END IF
            IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
                IF( SCALE( I ).GE.SFMAX1 / F ) &
                    GO TO 200
            END IF
            G = ONE / F
            SCALE( I ) = SCALE( I )*F
            NOCONV = .TRUE.
            !
            CALL SSCAL( N-K+1, G, A( I, K ), LDA )
            CALL SSCAL( L, F, A( 1, I ), 1 )
            !
200     END DO
        !
        IF( NOCONV ) &
            GO TO 140
        !
210     CONTINUE
        ILO = K
        IHI = L
        !
        RETURN
        !
        !     End of SGEBAL
        !
    END SUBROUTINE SGEBAL

    SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
        LDVR, WORK, LWORK, INFO )
        implicit none
        !
        !  -- LAPACK driver routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          JOBVL, JOBVR
        INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
        !     ..
        !     .. Array Arguments ..
        REAL   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
            WI( * ), WORK( * ), WR( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL   ZERO, ONE
        PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
        CHARACTER          SIDE
        INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K, &
            LWORK_TREVC, MAXWRK, MINWRK, NOUT
        REAL   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, &
            SN
        !     ..
        !     .. Local Arrays ..
        LOGICAL            SELECT( 1 )
        REAL   DUM( 1 )
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        INFO = 0
        LQUERY = ( LWORK.EQ.-1 )
        WANTVL = LSAME( JOBVL, 'V' )
        WANTVR = LSAME( JOBVR, 'V' )
        IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
            INFO = -1
        ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
            INFO = -9
        ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
            INFO = -11
        END IF
        !
        !     Compute workspace
        !      (Note: Comments in the code beginning "Workspace:" describe the
        !       minimal amount of workspace needed at that point in the code,
        !       as well as the preferred amount for good performance.
        !       NB refers to the optimal block size for the immediately
        !       following subroutine, as returned by ILAENV.
        !       HSWORK refers to the workspace preferred by SHSEQR, as
        !       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        !       the worst case.)
        !
        IF( INFO.EQ.0 ) THEN
            IF( N.EQ.0 ) THEN
                MINWRK = 1
                MAXWRK = 1
            ELSE
                MAXWRK = 2*N + N*ILAENV( 1, 'SGEHRD', ' ', N, 1, N, 0 )
                IF( WANTVL ) THEN
                    MINWRK = 4*N
                    MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, &
                        'SORGHR', ' ', N, 1, N, -1 ) )
                    CALL SHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, &
                        WORK, -1, INFO )
                    HSWORK = INT( WORK(1) )
                    MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
                    CALL STREVC3( 'L', 'B', SELECT, N, A, LDA, &
                        VL, LDVL, VR, LDVR, N, NOUT, &
                        WORK, -1, IERR )
                    LWORK_TREVC = INT( WORK(1) )
                    MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
                    MAXWRK = MAX( MAXWRK, 4*N )
                ELSE IF( WANTVR ) THEN
                    MINWRK = 4*N
                    MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, &
                        'SORGHR', ' ', N, 1, N, -1 ) )
                    CALL SHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, &
                        WORK, -1, INFO )
                    HSWORK = INT( WORK(1) )
                    MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
                    CALL STREVC3( 'R', 'B', SELECT, N, A, LDA, &
                        VL, LDVL, VR, LDVR, N, NOUT, &
                        WORK, -1, IERR )
                    LWORK_TREVC = INT( WORK(1) )
                    MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
                    MAXWRK = MAX( MAXWRK, 4*N )
                ELSE
                    MINWRK = 3*N
                    CALL SHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, &
                        WORK, -1, INFO )
                    HSWORK = INT( WORK(1) )
                    MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
                END IF
                MAXWRK = MAX( MAXWRK, MINWRK )
            END IF
            WORK( 1 ) = MAXWRK
            !
            IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
                INFO = -13
            END IF
        END IF
        !
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SGEEV ', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.EQ.0 ) &
            RETURN
        !
        !     Get machine constants
        !
        EPS = SLAMCH( 'P' )
        SMLNUM = SLAMCH( 'S' )
        BIGNUM = ONE / SMLNUM
        CALL SLABAD( SMLNUM, BIGNUM )
        SMLNUM = SQRT( SMLNUM ) / EPS
        BIGNUM = ONE / SMLNUM
        !
        !     Scale A if max element outside range [SMLNUM,BIGNUM]
        !
        ANRM = SLANGE( 'M', N, N, A, LDA, DUM )
        SCALEA = .FALSE.
        IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
            SCALEA = .TRUE.
            CSCALE = SMLNUM
        ELSE IF( ANRM.GT.BIGNUM ) THEN
            SCALEA = .TRUE.
            CSCALE = BIGNUM
        END IF
        IF( SCALEA ) &
            CALL SLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
        !
        !     Balance the matrix
        !     (Workspace: need N)
        !
        IBAL = 1
        CALL SGEBAL( 'B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )
        !
        !     Reduce to upper Hessenberg form
        !     (Workspace: need 3*N, prefer 2*N+N*NB)
        !
        ITAU = IBAL + N
        IWRK = ITAU + N
        CALL SGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), &
            LWORK-IWRK+1, IERR )
        !
        IF( WANTVL ) THEN
            !
            !        Want left eigenvectors
            !        Copy Householder vectors to VL
            !
            SIDE = 'L'
            CALL SLACPY( 'L', N, N, A, LDA, VL, LDVL )
            !
            !        Generate orthogonal matrix in VL
            !        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
            !
            CALL SORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), &
                LWORK-IWRK+1, IERR )
            !
            !        Perform QR iteration, accumulating Schur vectors in VL
            !        (Workspace: need N+1, prefer N+HSWORK (see comments) )
            !
            IWRK = ITAU
            CALL SHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, &
                WORK( IWRK ), LWORK-IWRK+1, INFO )
            !
            IF( WANTVR ) THEN
                !
                !           Want left and right eigenvectors
                !           Copy Schur vectors to VR
                !
                SIDE = 'B'
                CALL SLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
            END IF
            !
        ELSE IF( WANTVR ) THEN
            !
            !        Want right eigenvectors
            !        Copy Householder vectors to VR
            !
            SIDE = 'R'
            CALL SLACPY( 'L', N, N, A, LDA, VR, LDVR )
            !
            !        Generate orthogonal matrix in VR
            !        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
            !
            CALL SORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), &
                LWORK-IWRK+1, IERR )
            !
            !        Perform QR iteration, accumulating Schur vectors in VR
            !        (Workspace: need N+1, prefer N+HSWORK (see comments) )
            !
            IWRK = ITAU
            CALL SHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                WORK( IWRK ), LWORK-IWRK+1, INFO )
            !
        ELSE
            !
            !        Compute eigenvalues only
            !        (Workspace: need N+1, prefer N+HSWORK (see comments) )
            !
            IWRK = ITAU
            CALL SHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, &
                WORK( IWRK ), LWORK-IWRK+1, INFO )
        END IF
        !
        !     If INFO .NE. 0 from SHSEQR, then quit
        !
        IF( INFO.NE.0 ) &
            GO TO 50
        !
        IF( WANTVL .OR. WANTVR ) THEN
            !
            !        Compute left and/or right eigenvectors
            !        (Workspace: need 4*N, prefer N + N + 2*N*NB)
            !
            CALL STREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, &
                N, NOUT, WORK( IWRK ), LWORK-IWRK+1, IERR )
        END IF
        !
        IF( WANTVL ) THEN
            !
            !        Undo balancing of left eigenvectors
            !        (Workspace: need N)
            !
            CALL SGEBAK( 'B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL, &
                IERR )
            !
            !        Normalize left eigenvectors and make largest component real
            !
            DO  I = 1, N
                IF( WI( I ).EQ.ZERO ) THEN
                    SCL = ONE / SNRM2( N, VL( 1, I ), 1 )
                    CALL SSCAL( N, SCL, VL( 1, I ), 1 )
                ELSE IF( WI( I ).GT.ZERO ) THEN
                    SCL = ONE / SLAPY2( SNRM2( N, VL( 1, I ), 1 ), &
                        SNRM2( N, VL( 1, I+1 ), 1 ) )
                    CALL SSCAL( N, SCL, VL( 1, I ), 1 )
                    CALL SSCAL( N, SCL, VL( 1, I+1 ), 1 )
                    DO  K = 1, N
                        WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2
                    END DO
                    K = ISAMAX( N, WORK( IWRK ), 1 )
                    CALL SLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
                    CALL SROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
                    VL( K, I+1 ) = ZERO
                END IF
            END DO
        END IF
        !
        IF( WANTVR ) THEN
            !
            !        Undo balancing of right eigenvectors
            !        (Workspace: need N)
            !
            CALL SGEBAK( 'B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR, &
                IERR )
            !
            !        Normalize right eigenvectors and make largest component real
            !
            DO  I = 1, N
                IF( WI( I ).EQ.ZERO ) THEN
                    SCL = ONE / SNRM2( N, VR( 1, I ), 1 )
                    CALL SSCAL( N, SCL, VR( 1, I ), 1 )
                ELSE IF( WI( I ).GT.ZERO ) THEN
                    SCL = ONE / SLAPY2( SNRM2( N, VR( 1, I ), 1 ), &
                        SNRM2( N, VR( 1, I+1 ), 1 ) )
                    CALL SSCAL( N, SCL, VR( 1, I ), 1 )
                    CALL SSCAL( N, SCL, VR( 1, I+1 ), 1 )
                    DO  K = 1, N
                        WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2
                    END DO
                    K = ISAMAX( N, WORK( IWRK ), 1 )
                    CALL SLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
                    CALL SROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
                    VR( K, I+1 ) = ZERO
                END IF
            END DO
        END IF
        !
        !     Undo scaling if necessary
        !
50      CONTINUE
        IF( SCALEA ) THEN
            CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), &
                MAX( N-INFO, 1 ), IERR )
            CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), &
                MAX( N-INFO, 1 ), IERR )
            IF( INFO.GT.0 ) THEN
                CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, &
                    IERR )
                CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, &
                    IERR )
            END IF
        END IF
        !
        WORK( 1 ) = MAXWRK
        RETURN
        !
        !     End of SGEEV
        !
    END SUBROUTINE SGEEV

    SUBROUTINE SGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, ILO, INFO, LDA, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE
        PARAMETER          ( ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            I
        REAL               AII
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters
        !
        INFO = 0
        IF( N.LT.0 ) THEN
            INFO = -1
        ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
            INFO = -2
        ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SGEHD2', -INFO )
            RETURN
        END IF
        !
        DO  I = ILO, IHI - 1
            !
            !        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
            !
            CALL SLARFG( IHI-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, &
                TAU( I ) )
            AII = A( I+1, I )
            A( I+1, I ) = ONE
            !
            !        Apply H(i) to A(1:ihi,i+1:ihi) from the right
            !
            CALL SLARF( 'Right', IHI, IHI-I, A( I+1, I ), 1, TAU( I ), &
                A( 1, I+1 ), LDA, WORK )
            !
            !        Apply H(i) to A(i+1:ihi,i+1:n) from the left
            !
            CALL SLARF( 'Left', IHI-I, N-I, A( I+1, I ), 1, TAU( I ), &
                A( I+1, I+1 ), LDA, WORK )
            !
            A( I+1, I ) = AII
        END DO
        !
        RETURN
        !
        !     End of SGEHD2
        !
    END SUBROUTINE SGEHD2

    SUBROUTINE SGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, ILO, INFO, LDA, LWORK, N
        !     ..
        !     .. Array Arguments ..
        REAL              A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        INTEGER            NBMAX, LDT, TSIZE
        PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, &
            TSIZE = LDT*NBMAX )
        REAL              ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, &
            ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            LQUERY
        INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB, &
            NBMIN, NH, NX
        REAL              EI
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters
        !
        INFO = 0
        LQUERY = ( LWORK.EQ.-1 )
        IF( N.LT.0 ) THEN
            INFO = -1
        ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
            INFO = -2
        ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
        END IF
        !
        IF( INFO.EQ.0 ) THEN
            !
            !       Compute the workspace requirements
            !
            NB = MIN( NBMAX, ILAENV( 1, 'SGEHRD', ' ', N, ILO, IHI, -1 ) )
            LWKOPT = N*NB + TSIZE
            WORK( 1 ) = LWKOPT
        END IF
        !
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SGEHRD', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        !
        !     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
        !
        DO  I = 1, ILO - 1
            TAU( I ) = ZERO
        END DO
        DO  I = MAX( 1, IHI ), N - 1
            TAU( I ) = ZERO
        END DO
        !
        !     Quick return if possible
        !
        NH = IHI - ILO + 1
        IF( NH.LE.1 ) THEN
            WORK( 1 ) = 1
            RETURN
        END IF
        !
        !     Determine the block size
        !
        NB = MIN( NBMAX, ILAENV( 1, 'SGEHRD', ' ', N, ILO, IHI, -1 ) )
        NBMIN = 2
        IF( NB.GT.1 .AND. NB.LT.NH ) THEN
            !
            !        Determine when to cross over from blocked to unblocked code
            !        (last block is always handled by unblocked code)
            !
            NX = MAX( NB, ILAENV( 3, 'SGEHRD', ' ', N, ILO, IHI, -1 ) )
            IF( NX.LT.NH ) THEN
                !
                !           Determine if workspace is large enough for blocked code
                !
                IF( LWORK.LT.N*NB+TSIZE ) THEN
                    !
                    !              Not enough workspace to use optimal NB:  determine the
                    !              minimum value of NB, and reduce NB or force use of
                    !              unblocked code
                    !
                    NBMIN = MAX( 2, ILAENV( 2, 'SGEHRD', ' ', N, ILO, IHI, &
                        -1 ) )
                    IF( LWORK.GE.(N*NBMIN + TSIZE) ) THEN
                        NB = (LWORK-TSIZE) / N
                    ELSE
                        NB = 1
                    END IF
                END IF
            END IF
        END IF
        LDWORK = N
        !
        IF( NB.LT.NBMIN .OR. NB.GE.NH ) THEN
            !
            !        Use unblocked code below
            !
            I = ILO
            !
        ELSE
            !
            !        Use blocked code
            !
            IWT = 1 + N*NB
            DO  I = ILO, IHI - 1 - NX, NB
                IB = MIN( NB, IHI-I )
                !
                !           Reduce columns i:i+ib-1 to Hessenberg form, returning the
                !           matrices V and T of the block reflector H = I - V*T*V**T
                !           which performs the reduction, and also the matrix Y = A*V*T
                !
                CALL SLAHR2( IHI, I, IB, A( 1, I ), LDA, TAU( I ), &
                    WORK( IWT ), LDT, WORK, LDWORK )
                !
                !           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
                !           right, computing  A := A - Y * V**T. V(i+ib,ib-1) must be set
                !           to 1
                !
                EI = A( I+IB, I+IB-1 )
                A( I+IB, I+IB-1 ) = ONE
                CALL SGEMM( 'No transpose', 'Transpose', &
                    IHI, IHI-I-IB+1, &
                    IB, -ONE, WORK, LDWORK, A( I+IB, I ), LDA, ONE, &
                    A( 1, I+IB ), LDA )
                A( I+IB, I+IB-1 ) = EI
                !
                !           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
                !           right
                !
                CALL STRMM( 'Right', 'Lower', 'Transpose', &
                    'Unit', I, IB-1, &
                    ONE, A( I+1, I ), LDA, WORK, LDWORK )
                DO  J = 0, IB-2
                    CALL SAXPY( I, -ONE, WORK( LDWORK*J+1 ), 1, &
                        A( 1, I+J+1 ), 1 )
                END DO
                !
                !           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
                !           left
                !
                CALL SLARFB( 'Left', 'Transpose', 'Forward', &
                    'Columnwise', &
                    IHI-I, N-I-IB+1, IB, A( I+1, I ), LDA, &
                    WORK( IWT ), LDT, A( I+1, I+IB ), LDA, &
                    WORK, LDWORK )
            END DO
        END IF
        !
        !     Use unblocked code to reduce the rest of the matrix
        !
        CALL SGEHD2( N, I, IHI, A, LDA, TAU, WORK, IINFO )
        WORK( 1 ) = LWKOPT
        !
        RETURN
        !
        !     End of SGEHRD
        !
    END SUBROUTINE SGEHRD

    SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        !
        !  -- Reference BLAS level3 routine (version 3.7.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        REAL ALPHA,BETA
        INTEGER K,LDA,LDB,LDC,M,N
        CHARACTER TRANSA,TRANSB
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),B(LDB,*),C(LDC,*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !     .. Local Scalars ..
        REAL TEMP
        INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
        LOGICAL NOTA,NOTB
        !     ..
        !     .. Parameters ..
        REAL ONE,ZERO
        PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
        !     ..
        !
        !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
        !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
        !     and  columns of  A  and the  number of  rows  of  B  respectively.
        !
        NOTA = LSAME(TRANSA,'N')
        NOTB = LSAME(TRANSB,'N')
        IF (NOTA) THEN
            NROWA = M
            NCOLA = K
        ELSE
            NROWA = K
            NCOLA = M
        END IF
        IF (NOTB) THEN
            NROWB = K
        ELSE
            NROWB = N
        END IF
        !
        !     Test the input parameters.
        !
        INFO = 0
        IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
            (.NOT.LSAME(TRANSA,'T'))) THEN
            INFO = 1
        ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
            (.NOT.LSAME(TRANSB,'T'))) THEN
            INFO = 2
        ELSE IF (M.LT.0) THEN
            INFO = 3
        ELSE IF (N.LT.0) THEN
            INFO = 4
        ELSE IF (K.LT.0) THEN
            INFO = 5
        ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
            INFO = 8
        ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
            INFO = 10
        ELSE IF (LDC.LT.MAX(1,M)) THEN
            INFO = 13
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('SGEMM ',INFO)
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
            (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
        !
        !     And if  alpha.eq.zero.
        !
        IF (ALPHA.EQ.ZERO) THEN
            IF (BETA.EQ.ZERO) THEN
                DO  J = 1,N
                    DO  I = 1,M
                        C(I,J) = ZERO
                    END DO
                END DO
            ELSE
                DO  J = 1,N
                    DO  I = 1,M
                        C(I,J) = BETA*C(I,J)
                    END DO
                END DO
            END IF
            RETURN
        END IF
        !
        !     Start the operations.
        !
        IF (NOTB) THEN
            IF (NOTA) THEN
                !
                !           Form  C := alpha*A*B + beta*C.
                !
                DO  J = 1,N
                    IF (BETA.EQ.ZERO) THEN
                        DO  I = 1,M
                            C(I,J) = ZERO
                        END DO
                    ELSE IF (BETA.NE.ONE) THEN
                        DO  I = 1,M
                            C(I,J) = BETA*C(I,J)
                        END DO
                    END IF
                    DO  L = 1,K
                        TEMP = ALPHA*B(L,J)
                        DO  I = 1,M
                            C(I,J) = C(I,J) + TEMP*A(I,L)
                        END DO
                    END DO
                END DO
            ELSE
                !
                !           Form  C := alpha*A**T*B + beta*C
                !
                DO  J = 1,N
                    DO  I = 1,M
                        TEMP = ZERO
                        DO  L = 1,K
                            TEMP = TEMP + A(L,I)*B(L,J)
                        END DO
                        IF (BETA.EQ.ZERO) THEN
                            C(I,J) = ALPHA*TEMP
                        ELSE
                            C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                        END IF
                    END DO
                END DO
            END IF
        ELSE
            IF (NOTA) THEN
                !
                !           Form  C := alpha*A*B**T + beta*C
                !
                DO  J = 1,N
                    IF (BETA.EQ.ZERO) THEN
                        DO  I = 1,M
                            C(I,J) = ZERO
                        END DO
                    ELSE IF (BETA.NE.ONE) THEN
                        DO  I = 1,M
                            C(I,J) = BETA*C(I,J)
                        END DO
                    END IF
                    DO  L = 1,K
                        TEMP = ALPHA*B(J,L)
                        DO  I = 1,M
                            C(I,J) = C(I,J) + TEMP*A(I,L)
                        END DO
                    END DO
                END DO
            ELSE
                !
                !           Form  C := alpha*A**T*B**T + beta*C
                !
                DO  J = 1,N
                    DO  I = 1,M
                        TEMP = ZERO
                        DO  L = 1,K
                            TEMP = TEMP + A(L,I)*B(J,L)
                        END DO
                        IF (BETA.EQ.ZERO) THEN
                            C(I,J) = ALPHA*TEMP
                        ELSE
                            C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                        END IF
                    END DO
                END DO
            END IF
        END IF
        !
        RETURN
        !
        !     End of SGEMM .
        !
    END SUBROUTINE SGEMM

    SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        !
        !  -- Reference BLAS level2 routine (version 3.7.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        REAL ALPHA,BETA
        INTEGER INCX,INCY,LDA,M,N
        CHARACTER TRANS
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL ONE,ZERO
        PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
        !     ..
        !     .. Local Scalars ..
        REAL TEMP
        INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
        INFO = 0
        IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
            .NOT.LSAME(TRANS,'C')) THEN
            INFO = 1
        ELSE IF (M.LT.0) THEN
            INFO = 2
        ELSE IF (N.LT.0) THEN
            INFO = 3
        ELSE IF (LDA.LT.MAX(1,M)) THEN
            INFO = 6
        ELSE IF (INCX.EQ.0) THEN
            INFO = 8
        ELSE IF (INCY.EQ.0) THEN
            INFO = 11
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('SGEMV ',INFO)
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
            ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
        !
        !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
        !     up the start points in  X  and  Y.
        !
        IF (LSAME(TRANS,'N')) THEN
            LENX = N
            LENY = M
        ELSE
            LENX = M
            LENY = N
        END IF
        IF (INCX.GT.0) THEN
            KX = 1
        ELSE
            KX = 1 - (LENX-1)*INCX
        END IF
        IF (INCY.GT.0) THEN
            KY = 1
        ELSE
            KY = 1 - (LENY-1)*INCY
        END IF
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through A.
        !
        !     First form  y := beta*y.
        !
        IF (BETA.NE.ONE) THEN
            IF (INCY.EQ.1) THEN
                IF (BETA.EQ.ZERO) THEN
                    DO  I = 1,LENY
                        Y(I) = ZERO
                    END DO
                ELSE
                    DO  I = 1,LENY
                        Y(I) = BETA*Y(I)
                    END DO
                END IF
            ELSE
                IY = KY
                IF (BETA.EQ.ZERO) THEN
                    DO  I = 1,LENY
                        Y(IY) = ZERO
                        IY = IY + INCY
                    END DO
                ELSE
                    DO  I = 1,LENY
                        Y(IY) = BETA*Y(IY)
                        IY = IY + INCY
                    END DO
                END IF
            END IF
        END IF
        IF (ALPHA.EQ.ZERO) RETURN
        IF (LSAME(TRANS,'N')) THEN
            !
            !        Form  y := alpha*A*x + y.
            !
            JX = KX
            IF (INCY.EQ.1) THEN
                DO  J = 1,N
                    TEMP = ALPHA*X(JX)
                    DO  I = 1,M
                        Y(I) = Y(I) + TEMP*A(I,J)
                    END DO
                    JX = JX + INCX
                END DO
            ELSE
                DO  J = 1,N
                    TEMP = ALPHA*X(JX)
                    IY = KY
                    DO  I = 1,M
                        Y(IY) = Y(IY) + TEMP*A(I,J)
                        IY = IY + INCY
                    END DO
                    JX = JX + INCX
                END DO
            END IF
        ELSE
            !
            !        Form  y := alpha*A**T*x + y.
            !
            JY = KY
            IF (INCX.EQ.1) THEN
                DO  J = 1,N
                    TEMP = ZERO
                    DO  I = 1,M
                        TEMP = TEMP + A(I,J)*X(I)
                    END DO
                    Y(JY) = Y(JY) + ALPHA*TEMP
                    JY = JY + INCY
                END DO
            ELSE
                DO  J = 1,N
                    TEMP = ZERO
                    IX = KX
                    DO  I = 1,M
                        TEMP = TEMP + A(I,J)*X(IX)
                        IX = IX + INCX
                    END DO
                    Y(JY) = Y(JY) + ALPHA*TEMP
                    JY = JY + INCY
                END DO
            END IF
        END IF
        !
        RETURN
        !
        !     End of SGEMV .
        !
    END SUBROUTINE SGEMV

    SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        !
        !  -- Reference BLAS level2 routine (version 3.7.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        REAL ALPHA
        INTEGER INCX,INCY,LDA,M,N
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL ZERO
        PARAMETER (ZERO=0.0E+0)
        !     ..
        !     .. Local Scalars ..
        REAL TEMP
        INTEGER I,INFO,IX,J,JY,KX
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
        INFO = 0
        IF (M.LT.0) THEN
            INFO = 1
        ELSE IF (N.LT.0) THEN
            INFO = 2
        ELSE IF (INCX.EQ.0) THEN
            INFO = 5
        ELSE IF (INCY.EQ.0) THEN
            INFO = 7
        ELSE IF (LDA.LT.MAX(1,M)) THEN
            INFO = 9
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('SGER  ',INFO)
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through A.
        !
        IF (INCY.GT.0) THEN
            JY = 1
        ELSE
            JY = 1 - (N-1)*INCY
        END IF
        IF (INCX.EQ.1) THEN
            DO  J = 1,N
                IF (Y(JY).NE.ZERO) THEN
                    TEMP = ALPHA*Y(JY)
                    DO  I = 1,M
                        A(I,J) = A(I,J) + X(I)*TEMP
                    END DO
                END IF
                JY = JY + INCY
            END DO
        ELSE
            IF (INCX.GT.0) THEN
                KX = 1
            ELSE
                KX = 1 - (M-1)*INCX
            END IF
            DO  J = 1,N
                IF (Y(JY).NE.ZERO) THEN
                    TEMP = ALPHA*Y(JY)
                    IX = KX
                    DO  I = 1,M
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
                    END DO
                END IF
                JY = JY + INCY
            END DO
        END IF
        !
        RETURN
        !
        !     End of SGER  .
        !
    END SUBROUTINE SGER

    SUBROUTINE SHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z, &
        LDZ, WORK, LWORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
        CHARACTER          COMPZ, JOB
        !     ..
        !     .. Array Arguments ..
        REAL               H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
            Z( LDZ, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        !
        !     ==== Matrices of order NTINY or smaller must be processed by
        !     .    SLAHQR because of insufficient subdiagonal scratch space.
        !     .    (This is a hard limit.) ====
        INTEGER            NTINY
        PARAMETER          ( NTINY = 11 )
        !
        !     ==== NL allocates some local workspace to help small matrices
        !     .    through a rare SLAHQR failure.  NL .GT. NTINY = 11 is
        !     .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
        !     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
        !     .    allows up to six simultaneous shifts and a 16-by-16
        !     .    deflation window.  ====
        INTEGER            NL
        PARAMETER          ( NL = 49 )
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0e0, ONE = 1.0e0 )
        !     ..
        !     .. Local Arrays ..
        REAL               HL( NL, NL ), WORKL( NL )
        !     ..
        !     .. Local Scalars ..
        INTEGER            I, KBOT, NMIN
        LOGICAL            INITZ, LQUERY, WANTT, WANTZ
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN, REAL
        !     ..
        !     .. Executable Statements ..
        !
        !     ==== Decode and check the input parameters. ====
        !
        WANTT = LSAME( JOB, 'S' )
        INITZ = LSAME( COMPZ, 'I' )
        WANTZ = INITZ .OR. LSAME( COMPZ, 'V' )
        WORK( 1 ) = REAL( MAX( 1, N ) )
        LQUERY = LWORK.EQ.-1
        !
        INFO = 0
        IF( .NOT.LSAME( JOB, 'E' ) .AND. .NOT.WANTT ) THEN
            INFO = -1
        ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.WANTZ ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
            INFO = -4
        ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
            INFO = -5
        ELSE IF( LDH.LT.MAX( 1, N ) ) THEN
            INFO = -7
        ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
            INFO = -11
        ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -13
        END IF
        !
        IF( INFO.NE.0 ) THEN
            !
            !        ==== Quick return in case of invalid argument. ====
            !
            CALL XERBLA( 'SHSEQR', -INFO )
            RETURN
            !
        ELSE IF( N.EQ.0 ) THEN
            !
            !        ==== Quick return in case N = 0; nothing to do. ====
            !
            RETURN
            !
        ELSE IF( LQUERY ) THEN
            !
            !        ==== Quick return in case of a workspace query ====
            !
            CALL SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, &
                IHI, Z, LDZ, WORK, LWORK, INFO )
            !        ==== Ensure reported workspace size is backward-compatible with
            !        .    previous LAPACK versions. ====
            WORK( 1 ) = MAX( REAL( MAX( 1, N ) ), WORK( 1 ) )
            RETURN
            !
        ELSE
            !
            !        ==== copy eigenvalues isolated by SGEBAL ====
            !
            DO  I = 1, ILO - 1
                WR( I ) = H( I, I )
                WI( I ) = ZERO
            END DO
            DO  I = IHI + 1, N
                WR( I ) = H( I, I )
                WI( I ) = ZERO
            END DO
            !
            !        ==== Initialize Z, if requested ====
            !
            IF( INITZ ) &
                CALL SLASET( 'A', N, N, ZERO, ONE, Z, LDZ )
            !
            !        ==== Quick return if possible ====
            !
            IF( ILO.EQ.IHI ) THEN
                WR( ILO ) = H( ILO, ILO )
                WI( ILO ) = ZERO
                RETURN
            END IF
            !
            !        ==== SLAHQR/SLAQR0 crossover point ====
            !
            NMIN = ILAENV( 12, 'SHSEQR', JOB( : 1 ) // COMPZ( : 1 ), N, &
                ILO, IHI, LWORK )
            NMIN = MAX( NTINY, NMIN )
            !
            !        ==== SLAQR0 for big matrices; SLAHQR for small ones ====
            !
            IF( N.GT.NMIN ) THEN
                CALL SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, &
                    IHI, Z, LDZ, WORK, LWORK, INFO )
            ELSE
                !
                !           ==== Small matrix ====
                !
                CALL SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, &
                    IHI, Z, LDZ, INFO )
                !
                IF( INFO.GT.0 ) THEN
                    !
                    !              ==== A rare SLAHQR failure!  SLAQR0 sometimes succeeds
                    !              .    when SLAHQR fails. ====
                    !
                    KBOT = INFO
                    !
                    IF( N.GE.NL ) THEN
                        !
                        !                 ==== Larger matrices have enough subdiagonal scratch
                        !                 .    space to call SLAQR0 directly. ====
                        !
                        CALL SLAQR0( WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, &
                            WI, ILO, IHI, Z, LDZ, WORK, LWORK, INFO )
                        !
                    ELSE
                        !
                        !                 ==== Tiny matrices don't have enough subdiagonal
                        !                 .    scratch space to benefit from SLAQR0.  Hence,
                        !                 .    tiny matrices must be copied into a larger
                        !                 .    array before calling SLAQR0. ====
                        !
                        CALL SLACPY( 'A', N, N, H, LDH, HL, NL )
                        HL( N+1, N ) = ZERO
                        CALL SLASET( 'A', NL, NL-N, ZERO, ZERO, HL( 1, N+1 ), &
                            NL )
                        CALL SLAQR0( WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, &
                            WI, ILO, IHI, Z, LDZ, WORKL, NL, INFO )
                        IF( WANTT .OR. INFO.NE.0 ) &
                            CALL SLACPY( 'A', N, N, HL, NL, H, LDH )
                    END IF
                END IF
            END IF
            !
            !        ==== Clear out the trash, if necessary. ====
            !
            IF( ( WANTT .OR. INFO.NE.0 ) .AND. N.GT.2 ) &
                CALL SLASET( 'L', N-2, N-2, ZERO, ZERO, H( 3, 1 ), LDH )
            !
            !        ==== Ensure reported workspace size is backward-compatible with
            !        .    previous LAPACK versions. ====
            !
            WORK( 1 ) = MAX( REAL( MAX( 1, N ) ), WORK( 1 ) )
        END IF
        !
        !     ==== End of SHSEQR ====
        !
    END SUBROUTINE SHSEQR

    SUBROUTINE SLABAD( SMALL, LARGE )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        REAL               LARGE, SMALL
        !     ..
        !
        !  =====================================================================
        !
        !     .. Intrinsic Functions ..
        INTRINSIC          LOG10, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     If it looks like we're on a Cray, take the square root of
        !     SMALL and LARGE to avoid overflow and underflow problems.
        !
        IF( LOG10( LARGE ).GT.2000. ) THEN
            SMALL = SQRT( SMALL )
            LARGE = SQRT( LARGE )
        END IF
        !
        RETURN
        !
        !     End of SLABAD
        !
    END SUBROUTINE SLABAD

    SUBROUTINE SLACPY( UPLO, M, N, A, LDA, B, LDB )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          UPLO
        INTEGER            LDA, LDB, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), B( LDB, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        INTEGER            I, J
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MIN
        !     ..
        !     .. Executable Statements ..
        !
        IF( LSAME( UPLO, 'U' ) ) THEN
            DO  J = 1, N
                DO  I = 1, MIN( J, M )
                    B( I, J ) = A( I, J )
                END DO
            END DO
        ELSE IF( LSAME( UPLO, 'L' ) ) THEN
            DO  J = 1, N
                DO  I = J, M
                    B( I, J ) = A( I, J )
                END DO
            END DO
        ELSE
            DO  J = 1, N
                DO  I = 1, M
                    B( I, J ) = A( I, J )
                END DO
            END DO
        END IF
        RETURN
        !
        !     End of SLACPY
        !
    END SUBROUTINE SLACPY

    SUBROUTINE SLADIV( A, B, C, D, P, Q )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     January 2013
        !
        !     .. Scalar Arguments ..
        REAL               A, B, C, D, P, Q
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               BS
        PARAMETER          ( BS = 2.0E0 )
        REAL               HALF
        PARAMETER          ( HALF = 0.5E0 )
        REAL               TWO
        PARAMETER          ( TWO = 2.0E0 )
        !
        !     .. Local Scalars ..
        REAL               AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX
        !     ..
        !     .. Executable Statements ..
        !
        AA = A
        BB = B
        CC = C
        DD = D
        AB = MAX( ABS(A), ABS(B) )
        CD = MAX( ABS(C), ABS(D) )
        S = 1.0E0

        OV = SLAMCH( 'Overflow threshold' )
        UN = SLAMCH( 'Safe minimum' )
        EPS = SLAMCH( 'Epsilon' )
        BE = BS / (EPS*EPS)

        IF( AB >= HALF*OV ) THEN
            AA = HALF * AA
            BB = HALF * BB
            S  = TWO * S
        END IF
        IF( CD >= HALF*OV ) THEN
            CC = HALF * CC
            DD = HALF * DD
            S  = HALF * S
        END IF
        IF( AB <= UN*BS/EPS ) THEN
            AA = AA * BE
            BB = BB * BE
            S  = S / BE
        END IF
        IF( CD <= UN*BS/EPS ) THEN
            CC = CC * BE
            DD = DD * BE
            S  = S * BE
        END IF
        IF( ABS( D ).LE.ABS( C ) ) THEN
            CALL SLADIV1(AA, BB, CC, DD, P, Q)
        ELSE
            CALL SLADIV1(BB, AA, DD, CC, P, Q)
            Q = -Q
        END IF
        P = P * S
        Q = Q * S
        !
        RETURN
        !
        !     End of SLADIV
        !
    END SUBROUTINE SLADIV

    !> \ingroup realOTHERauxiliary


    SUBROUTINE SLADIV1( A, B, C, D, P, Q )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     January 2013
        !
        !     .. Scalar Arguments ..
        REAL               A, B, C, D, P, Q
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE
        PARAMETER          ( ONE = 1.0E0 )
        !
        !     .. Local Scalars ..
        REAL               R, T
        !     ..
        !     .. Executable Statements ..
        !
        R = D / C
        T = ONE / (C + D * R)
        P = SLADIV2(A, B, C, D, R, T)
        A = -A
        Q = SLADIV2(B, A, C, D, R, T)
        !
        RETURN
        !
        !     End of SLADIV1
        !
    END SUBROUTINE SLADIV1

    !> \ingroup realOTHERauxiliary

    REAL FUNCTION SLADIV2( A, B, C, D, R, T )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     January 2013
        !
        !     .. Scalar Arguments ..
        REAL               A, B, C, D, R, T
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO
        PARAMETER          ( ZERO = 0.0E0 )
        !
        !     .. Local Scalars ..
        REAL               BR
        !     ..
        !     .. Executable Statements ..
        !
        IF( R.NE.ZERO ) THEN
            BR = B * R
            if( BR.NE.ZERO ) THEN
                SLADIV2 = (A + BR) * T
            ELSE
                SLADIV2 = A * T + (B * T) * R
            END IF
        ELSE
            SLADIV2 = (A + D * (B / C)) * T
        END IF
        !
        RETURN
        !
        !     End of SLADIV
        !
    END FUNCTION SLADIV2

    SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
        ILOZ, IHIZ, Z, LDZ, INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
        LOGICAL            WANTT, WANTZ
        !     ..
        !     .. Array Arguments ..
        REAL               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
        !     ..
        !
        !  =========================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE, TWO
        PARAMETER          ( ZERO = 0.0e0, ONE = 1.0e0, TWO = 2.0e0 )
        REAL               DAT1, DAT2
        PARAMETER          ( DAT1 = 3.0e0 / 4.0e0, DAT2 = -0.4375e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               AA, AB, BA, BB, CS, DET, H11, H12, H21, H21S, &
            H22, RT1I, RT1R, RT2I, RT2R, RTDISC, S, SAFMAX, &
            SAFMIN, SMLNUM, SN, SUM, T1, T2, T3, TR, TST, &
            ULP, V2, V3
        INTEGER            I, I1, I2, ITS, ITMAX, J, K, L, M, NH, NR, NZ
        !     ..
        !     .. Local Arrays ..
        REAL               V( 3 )
        !     ..
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX, MIN, REAL, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        INFO = 0
        !
        !     Quick return if possible
        !
        IF( N.EQ.0 ) &
            RETURN
        IF( ILO.EQ.IHI ) THEN
            WR( ILO ) = H( ILO, ILO )
            WI( ILO ) = ZERO
            RETURN
        END IF
        !
        !     ==== clear out the trash ====
        DO  J = ILO, IHI - 3
            H( J+2, J ) = ZERO
            H( J+3, J ) = ZERO
        END DO
        IF( ILO.LE.IHI-2 ) &
            H( IHI, IHI-2 ) = ZERO
        !
        NH = IHI - ILO + 1
        NZ = IHIZ - ILOZ + 1
        !
        !     Set machine-dependent constants for the stopping criterion.
        !
        SAFMIN = SLAMCH( 'SAFE MINIMUM' )
        SAFMAX = ONE / SAFMIN
        CALL SLABAD( SAFMIN, SAFMAX )
        ULP = SLAMCH( 'PRECISION' )
        SMLNUM = SAFMIN*( REAL( NH ) / ULP )
        !
        !     I1 and I2 are the indices of the first row and last column of H
        !     to which transformations must be applied. If eigenvalues only are
        !     being computed, I1 and I2 are set inside the main loop.
        !
        IF( WANTT ) THEN
            I1 = 1
            I2 = N
        END IF
        !
        !     ITMAX is the total number of QR iterations allowed.
        !
        ITMAX = 30 * MAX( 10, NH )
        !
        !     The main loop begins here. I is the loop index and decreases from
        !     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
        !     with the active submatrix in rows and columns L to I.
        !     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
        !     H(L,L-1) is negligible so that the matrix splits.
        !
        I = IHI
20      CONTINUE
        L = ILO
        IF( I.LT.ILO ) &
            GO TO 160
        !
        !     Perform QR iterations on rows and columns ILO to I until a
        !     submatrix of order 1 or 2 splits off at the bottom because a
        !     subdiagonal element has become negligible.
        !
        DO  ITS = 0, ITMAX
            !
            !        Look for a single small subdiagonal element.
            !
            DO  K = I, L + 1, -1
                IF( ABS( H( K, K-1 ) ).LE.SMLNUM ) &
                    GO TO 40
                TST = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
                IF( TST.EQ.ZERO ) THEN
                    IF( K-2.GE.ILO ) &
                        TST = TST + ABS( H( K-1, K-2 ) )
                    IF( K+1.LE.IHI ) &
                        TST = TST + ABS( H( K+1, K ) )
                END IF
                !           ==== The following is a conservative small subdiagonal
                !           .    deflation  criterion due to Ahues & Tisseur (LAWN 122,
                !           .    1997). It has better mathematical foundation and
                !           .    improves accuracy in some cases.  ====
                IF( ABS( H( K, K-1 ) ).LE.ULP*TST ) THEN
                    AB = MAX( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
                    BA = MIN( ABS( H( K, K-1 ) ), ABS( H( K-1, K ) ) )
                    AA = MAX( ABS( H( K, K ) ), &
                        ABS( H( K-1, K-1 )-H( K, K ) ) )
                    BB = MIN( ABS( H( K, K ) ), &
                        ABS( H( K-1, K-1 )-H( K, K ) ) )
                    S = AA + AB
                    IF( BA*( AB / S ).LE.MAX( SMLNUM, &
                        ULP*( BB*( AA / S ) ) ) )GO TO 40
                END IF
            END DO
40          CONTINUE
            L = K
            IF( L.GT.ILO ) THEN
                !
                !           H(L,L-1) is negligible
                !
                H( L, L-1 ) = ZERO
            END IF
            !
            !        Exit from loop if a submatrix of order 1 or 2 has split off.
            !
            IF( L.GE.I-1 ) &
                GO TO 150
            !
            !        Now the active submatrix is in rows and columns L to I. If
            !        eigenvalues only are being computed, only the active submatrix
            !        need be transformed.
            !
            IF( .NOT.WANTT ) THEN
                I1 = L
                I2 = I
            END IF
            !
            IF( ITS.EQ.10 ) THEN
                !
                !           Exceptional shift.
                !
                S = ABS( H( L+1, L ) ) + ABS( H( L+2, L+1 ) )
                H11 = DAT1*S + H( L, L )
                H12 = DAT2*S
                H21 = S
                H22 = H11
            ELSE IF( ITS.EQ.20 ) THEN
                !
                !           Exceptional shift.
                !
                S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
                H11 = DAT1*S + H( I, I )
                H12 = DAT2*S
                H21 = S
                H22 = H11
            ELSE
                !
                !           Prepare to use Francis' double shift
                !           (i.e. 2nd degree generalized Rayleigh quotient)
                !
                H11 = H( I-1, I-1 )
                H21 = H( I, I-1 )
                H12 = H( I-1, I )
                H22 = H( I, I )
            END IF
            S = ABS( H11 ) + ABS( H12 ) + ABS( H21 ) + ABS( H22 )
            IF( S.EQ.ZERO ) THEN
                RT1R = ZERO
                RT1I = ZERO
                RT2R = ZERO
                RT2I = ZERO
            ELSE
                H11 = H11 / S
                H21 = H21 / S
                H12 = H12 / S
                H22 = H22 / S
                TR = ( H11+H22 ) / TWO
                DET = ( H11-TR )*( H22-TR ) - H12*H21
                RTDISC = SQRT( ABS( DET ) )
                IF( DET.GE.ZERO ) THEN
                    !
                    !              ==== complex conjugate shifts ====
                    !
                    RT1R = TR*S
                    RT2R = RT1R
                    RT1I = RTDISC*S
                    RT2I = -RT1I
                ELSE
                    !
                    !              ==== real shifts (use only one of them)  ====
                    !
                    RT1R = TR + RTDISC
                    RT2R = TR - RTDISC
                    IF( ABS( RT1R-H22 ).LE.ABS( RT2R-H22 ) ) THEN
                        RT1R = RT1R*S
                        RT2R = RT1R
                    ELSE
                        RT2R = RT2R*S
                        RT1R = RT2R
                    END IF
                    RT1I = ZERO
                    RT2I = ZERO
                END IF
            END IF
            !
            !        Look for two consecutive small subdiagonal elements.
            !
            DO  M = I - 2, L, -1
                !           Determine the effect of starting the double-shift QR
                !           iteration at row M, and see if this would make H(M,M-1)
                !           negligible.  (The following uses scaling to avoid
                !           overflows and most underflows.)
                !
                H21S = H( M+1, M )
                S = ABS( H( M, M )-RT2R ) + ABS( RT2I ) + ABS( H21S )
                H21S = H( M+1, M ) / S
                V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )* &
                    ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S )
                V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R )
                V( 3 ) = H21S*H( M+2, M+1 )
                S = ABS( V( 1 ) ) + ABS( V( 2 ) ) + ABS( V( 3 ) )
                V( 1 ) = V( 1 ) / S
                V( 2 ) = V( 2 ) / S
                V( 3 ) = V( 3 ) / S
                IF( M.EQ.L ) &
                    GO TO 60
                IF( ABS( H( M, M-1 ) )*( ABS( V( 2 ) )+ABS( V( 3 ) ) ).LE. &
                    ULP*ABS( V( 1 ) )*( ABS( H( M-1, M-1 ) )+ABS( H( M, &
                    M ) )+ABS( H( M+1, M+1 ) ) ) )GO TO 60
            END DO
60          CONTINUE
            !
            !        Double-shift QR step
            !
            DO  K = M, I - 1
                !
                !           The first iteration of this loop determines a reflection G
                !           from the vector V and applies it from left and right to H,
                !           thus creating a nonzero bulge below the subdiagonal.
                !
                !           Each subsequent iteration determines a reflection G to
                !           restore the Hessenberg form in the (K-1)th column, and thus
                !           chases the bulge one step toward the bottom of the active
                !           submatrix. NR is the order of G.
                !
                NR = MIN( 3, I-K+1 )
                IF( K.GT.M ) &
                    CALL SCOPY( NR, H( K, K-1 ), 1, V, 1 )
                CALL SLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
                IF( K.GT.M ) THEN
                    H( K, K-1 ) = V( 1 )
                    H( K+1, K-1 ) = ZERO
                    IF( K.LT.I-1 ) &
                        H( K+2, K-1 ) = ZERO
                ELSE IF( M.GT.L ) THEN
                    !               ==== Use the following instead of
                    !               .    H( K, K-1 ) = -H( K, K-1 ) to
                    !               .    avoid a bug when v(2) and v(3)
                    !               .    underflow. ====
                    H( K, K-1 ) = H( K, K-1 )*( ONE-T1 )
                END IF
                V2 = V( 2 )
                T2 = T1*V2
                IF( NR.EQ.3 ) THEN
                    V3 = V( 3 )
                    T3 = T1*V3
                    !
                    !              Apply G from the left to transform the rows of the matrix
                    !              in columns K to I2.
                    !
                    DO  J = K, I2
                        SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                        H( K, J ) = H( K, J ) - SUM*T1
                        H( K+1, J ) = H( K+1, J ) - SUM*T2
                        H( K+2, J ) = H( K+2, J ) - SUM*T3
                    END DO
                    !
                    !              Apply G from the right to transform the columns of the
                    !              matrix in rows I1 to min(K+3,I).
                    !
                    DO  J = I1, MIN( K+3, I )
                        SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                        H( J, K ) = H( J, K ) - SUM*T1
                        H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                        H( J, K+2 ) = H( J, K+2 ) - SUM*T3
                    END DO
                    !
                    IF( WANTZ ) THEN
                        !
                        !                 Accumulate transformations in the matrix Z
                        !
                        DO  J = ILOZ, IHIZ
                            SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                            Z( J, K ) = Z( J, K ) - SUM*T1
                            Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                            Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
                        END DO
                    END IF
                ELSE IF( NR.EQ.2 ) THEN
                    !
                    !              Apply G from the left to transform the rows of the matrix
                    !              in columns K to I2.
                    !
                    DO  J = K, I2
                        SUM = H( K, J ) + V2*H( K+1, J )
                        H( K, J ) = H( K, J ) - SUM*T1
                        H( K+1, J ) = H( K+1, J ) - SUM*T2
                    END DO
                    !
                    !              Apply G from the right to transform the columns of the
                    !              matrix in rows I1 to min(K+3,I).
                    !
                    DO  J = I1, I
                        SUM = H( J, K ) + V2*H( J, K+1 )
                        H( J, K ) = H( J, K ) - SUM*T1
                        H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                    END DO
                    !
                    IF( WANTZ ) THEN
                        !
                        !                 Accumulate transformations in the matrix Z
                        !
                        DO  J = ILOZ, IHIZ
                            SUM = Z( J, K ) + V2*Z( J, K+1 )
                            Z( J, K ) = Z( J, K ) - SUM*T1
                            Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                        END DO
                    END IF
                END IF
            END DO
            !
        END DO
        !
        !     Failure to converge in remaining number of iterations
        !
        INFO = I
        RETURN
        !
150     CONTINUE
        !
        IF( L.EQ.I ) THEN
            !
            !        H(I,I-1) is negligible: one eigenvalue has converged.
            !
            WR( I ) = H( I, I )
            WI( I ) = ZERO
        ELSE IF( L.EQ.I-1 ) THEN
            !
            !        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
            !
            !        Transform the 2-by-2 submatrix to standard Schur form,
            !        and compute and store the eigenvalues.
            !
            CALL SLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ), &
                H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ), &
                CS, SN )
            !
            IF( WANTT ) THEN
                !
                !           Apply the transformation to the rest of H.
                !
                IF( I2.GT.I ) &
                    CALL SROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH, &
                    CS, SN )
                CALL SROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
            END IF
            IF( WANTZ ) THEN
                !
                !           Apply the transformation to Z.
                !
                CALL SROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
            END IF
        END IF
        !
        !     return to start of the main loop with new value of I.
        !
        I = L - 1
        GO TO 20
        !
160     CONTINUE
        RETURN
        !
        !     End of SLAHQR
        !
    END SUBROUTINE SLAHQR

    SUBROUTINE SLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            K, LDA, LDT, LDY, N, NB
        !     ..
        !     .. Array Arguments ..
        REAL              A( LDA, * ), T( LDT, NB ), TAU( NB ), &
            Y( LDY, NB )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL              ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, &
            ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            I
        REAL              EI
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick return if possible
        !
        IF( N.LE.1 ) &
            RETURN
        !
        DO  I = 1, NB
            IF( I.GT.1 ) THEN
                !
                !           Update A(K+1:N,I)
                !
                !           Update I-th column of A - Y * V**T
                !
                CALL SGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, Y(K+1,1), LDY, &
                    A( K+I-1, 1 ), LDA, ONE, A( K+1, I ), 1 )
                !
                !           Apply I - V * T**T * V**T to this column (call it b) from the
                !           left, using the last column of T as workspace
                !
                !           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
                !                    ( V2 )             ( b2 )
                !
                !           where V1 is unit lower triangular
                !
                !           w := V1**T * b1
                !
                CALL SCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
                CALL STRMV( 'Lower', 'Transpose', 'UNIT', &
                    I-1, A( K+1, 1 ), &
                    LDA, T( 1, NB ), 1 )
                !
                !           w := w + V2**T * b2
                !
                CALL SGEMV( 'Transpose', N-K-I+1, I-1, &
                    ONE, A( K+I, 1 ), &
                    LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
                !
                !           w := T**T * w
                !
                CALL STRMV( 'Upper', 'Transpose', 'NON-UNIT', &
                    I-1, T, LDT, &
                    T( 1, NB ), 1 )
                !
                !           b2 := b2 - V2*w
                !
                CALL SGEMV( 'NO TRANSPOSE', N-K-I+1, I-1, -ONE, &
                    A( K+I, 1 ), &
                    LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
                !
                !           b1 := b1 - V1*w
                !
                CALL STRMV( 'Lower', 'NO TRANSPOSE', &
                    'UNIT', I-1, &
                    A( K+1, 1 ), LDA, T( 1, NB ), 1 )
                CALL SAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
                !
                A( K+I-1, I-1 ) = EI
            END IF
            !
            !        Generate the elementary reflector H(I) to annihilate
            !        A(K+I+1:N,I)
            !
            CALL SLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1, &
                TAU( I ) )
            EI = A( K+I, I )
            A( K+I, I ) = ONE
            !
            !        Compute  Y(K+1:N,I)
            !
            CALL SGEMV( 'NO TRANSPOSE', N-K, N-K-I+1, &
                ONE, A( K+1, I+1 ), &
                LDA, A( K+I, I ), 1, ZERO, Y( K+1, I ), 1 )
            CALL SGEMV( 'Transpose', N-K-I+1, I-1, &
                ONE, A( K+I, 1 ), LDA, &
                A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
            CALL SGEMV( 'NO TRANSPOSE', N-K, I-1, -ONE, &
                Y( K+1, 1 ), LDY, &
                T( 1, I ), 1, ONE, Y( K+1, I ), 1 )
            CALL SSCAL( N-K, TAU( I ), Y( K+1, I ), 1 )
            !
            !        Compute T(1:I,I)
            !
            CALL SSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
            CALL STRMV( 'Upper', 'No Transpose', 'NON-UNIT', &
                I-1, T, LDT, &
                T( 1, I ), 1 )
            T( I, I ) = TAU( I )
            !
        END DO
        A( K+NB, NB ) = EI
        !
        !     Compute Y(1:K,1:NB)
        !
        CALL SLACPY( 'ALL', K, NB, A( 1, 2 ), LDA, Y, LDY )
        CALL STRMM( 'RIGHT', 'Lower', 'NO TRANSPOSE', &
            'UNIT', K, NB, &
            ONE, A( K+1, 1 ), LDA, Y, LDY )
        IF( N.GT.K+NB ) &
            CALL SGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', K, &
            NB, N-K-NB, ONE, &
            A( 1, 2+NB ), LDA, A( K+1+NB, 1 ), LDA, ONE, Y, &
            LDY )
        CALL STRMM( 'RIGHT', 'Upper', 'NO TRANSPOSE', &
            'NON-UNIT', K, NB, &
            ONE, T, LDT, Y, LDY )
        !
        RETURN
        !
        !     End of SLAHR2
        !
    END SUBROUTINE SLAHR2

    SUBROUTINE SLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, &
        LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        LOGICAL            LTRANS
        INTEGER            INFO, LDA, LDB, LDX, NA, NW
        REAL               CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), B( LDB, * ), X( LDX, * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
        REAL               TWO
        PARAMETER          ( TWO = 2.0E0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            ICMAX, J
        REAL               BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21, &
            CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21, &
            LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R, &
            UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S, &
            UR22, XI1, XI2, XR1, XR2
        !     ..
        !     .. Local Arrays ..
        LOGICAL            CSWAP( 4 ), RSWAP( 4 )
        INTEGER            IPIVOT( 4, 4 )
        REAL               CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX
        !     ..
        !     .. Equivalences ..
        EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ), &
            ( CR( 1, 1 ), CRV( 1 ) )
        !     ..
        !     .. Data statements ..
        DATA               CSWAP / .FALSE., .FALSE., .TRUE., .TRUE. /
        DATA               RSWAP / .FALSE., .TRUE., .FALSE., .TRUE. /
        DATA               IPIVOT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, &
            3, 2, 1 /
        !     ..
        !     .. Executable Statements ..
        !
        !     Compute BIGNUM
        !
        SMLNUM = TWO*SLAMCH( 'Safe minimum' )
        BIGNUM = ONE / SMLNUM
        SMINI = MAX( SMIN, SMLNUM )
        !
        !     Don't check for input errors
        !
        INFO = 0
        !
        !     Standard Initializations
        !
        SCALE = ONE
        !
        IF( NA.EQ.1 ) THEN
            !
            !        1 x 1  (i.e., scalar) system   C X = B
            !
            IF( NW.EQ.1 ) THEN
                !
                !           Real 1x1 system.
                !
                !           C = ca A - w D
                !
                CSR = CA*A( 1, 1 ) - WR*D1
                CNORM = ABS( CSR )
                !
                !           If | C | < SMINI, use C = SMINI
                !
                IF( CNORM.LT.SMINI ) THEN
                    CSR = SMINI
                    CNORM = SMINI
                    INFO = 1
                END IF
                !
                !           Check scaling for  X = B / C
                !
                BNORM = ABS( B( 1, 1 ) )
                IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
                    IF( BNORM.GT.BIGNUM*CNORM ) &
                        SCALE = ONE / BNORM
                END IF
                !
                !           Compute X
                !
                X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / CSR
                XNORM = ABS( X( 1, 1 ) )
            ELSE
                !
                !           Complex 1x1 system (w is complex)
                !
                !           C = ca A - w D
                !
                CSR = CA*A( 1, 1 ) - WR*D1
                CSI = -WI*D1
                CNORM = ABS( CSR ) + ABS( CSI )
                !
                !           If | C | < SMINI, use C = SMINI
                !
                IF( CNORM.LT.SMINI ) THEN
                    CSR = SMINI
                    CSI = ZERO
                    CNORM = SMINI
                    INFO = 1
                END IF
                !
                !           Check scaling for  X = B / C
                !
                BNORM = ABS( B( 1, 1 ) ) + ABS( B( 1, 2 ) )
                IF( CNORM.LT.ONE .AND. BNORM.GT.ONE ) THEN
                    IF( BNORM.GT.BIGNUM*CNORM ) &
                        SCALE = ONE / BNORM
                END IF
                !
                !           Compute X
                !
                CALL SLADIV( SCALE*B( 1, 1 ), SCALE*B( 1, 2 ), CSR, CSI, &
                    X( 1, 1 ), X( 1, 2 ) )
                XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
            END IF
            !
        ELSE
            !
            !        2x2 System
            !
            !        Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
            !
            CR( 1, 1 ) = CA*A( 1, 1 ) - WR*D1
            CR( 2, 2 ) = CA*A( 2, 2 ) - WR*D2
            IF( LTRANS ) THEN
                CR( 1, 2 ) = CA*A( 2, 1 )
                CR( 2, 1 ) = CA*A( 1, 2 )
            ELSE
                CR( 2, 1 ) = CA*A( 2, 1 )
                CR( 1, 2 ) = CA*A( 1, 2 )
            END IF
            !
            IF( NW.EQ.1 ) THEN
                !
                !           Real 2x2 system  (w is real)
                !
                !           Find the largest element in C
                !
                CMAX = ZERO
                ICMAX = 0
                !
                DO  J = 1, 4
                    IF( ABS( CRV( J ) ).GT.CMAX ) THEN
                        CMAX = ABS( CRV( J ) )
                        ICMAX = J
                    END IF
                END DO
                !
                !           If norm(C) < SMINI, use SMINI*identity.
                !
                IF( CMAX.LT.SMINI ) THEN
                    BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 2, 1 ) ) )
                    IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                        IF( BNORM.GT.BIGNUM*SMINI ) &
                            SCALE = ONE / BNORM
                    END IF
                    TEMP = SCALE / SMINI
                    X( 1, 1 ) = TEMP*B( 1, 1 )
                    X( 2, 1 ) = TEMP*B( 2, 1 )
                    XNORM = TEMP*BNORM
                    INFO = 1
                    RETURN
                END IF
                !
                !           Gaussian elimination with complete pivoting.
                !
                UR11 = CRV( ICMAX )
                CR21 = CRV( IPIVOT( 2, ICMAX ) )
                UR12 = CRV( IPIVOT( 3, ICMAX ) )
                CR22 = CRV( IPIVOT( 4, ICMAX ) )
                UR11R = ONE / UR11
                LR21 = UR11R*CR21
                UR22 = CR22 - UR12*LR21
                !
                !           If smaller pivot < SMINI, use SMINI
                !
                IF( ABS( UR22 ).LT.SMINI ) THEN
                    UR22 = SMINI
                    INFO = 1
                END IF
                IF( RSWAP( ICMAX ) ) THEN
                    BR1 = B( 2, 1 )
                    BR2 = B( 1, 1 )
                ELSE
                    BR1 = B( 1, 1 )
                    BR2 = B( 2, 1 )
                END IF
                BR2 = BR2 - LR21*BR1
                BBND = MAX( ABS( BR1*( UR22*UR11R ) ), ABS( BR2 ) )
                IF( BBND.GT.ONE .AND. ABS( UR22 ).LT.ONE ) THEN
                    IF( BBND.GE.BIGNUM*ABS( UR22 ) ) &
                        SCALE = ONE / BBND
                END IF
                !
                XR2 = ( BR2*SCALE ) / UR22
                XR1 = ( SCALE*BR1 )*UR11R - XR2*( UR11R*UR12 )
                IF( CSWAP( ICMAX ) ) THEN
                    X( 1, 1 ) = XR2
                    X( 2, 1 ) = XR1
                ELSE
                    X( 1, 1 ) = XR1
                    X( 2, 1 ) = XR2
                END IF
                XNORM = MAX( ABS( XR1 ), ABS( XR2 ) )
                !
                !           Further scaling if  norm(A) norm(X) > overflow
                !
                IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
                    IF( XNORM.GT.BIGNUM / CMAX ) THEN
                        TEMP = CMAX / BIGNUM
                        X( 1, 1 ) = TEMP*X( 1, 1 )
                        X( 2, 1 ) = TEMP*X( 2, 1 )
                        XNORM = TEMP*XNORM
                        SCALE = TEMP*SCALE
                    END IF
                END IF
            ELSE
                !
                !           Complex 2x2 system  (w is complex)
                !
                !           Find the largest element in C
                !
                CI( 1, 1 ) = -WI*D1
                CI( 2, 1 ) = ZERO
                CI( 1, 2 ) = ZERO
                CI( 2, 2 ) = -WI*D2
                CMAX = ZERO
                ICMAX = 0
                !
                DO  J = 1, 4
                    IF( ABS( CRV( J ) )+ABS( CIV( J ) ).GT.CMAX ) THEN
                        CMAX = ABS( CRV( J ) ) + ABS( CIV( J ) )
                        ICMAX = J
                    END IF
                END DO
                !
                !           If norm(C) < SMINI, use SMINI*identity.
                !
                IF( CMAX.LT.SMINI ) THEN
                    BNORM = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ), &
                        ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) )
                    IF( SMINI.LT.ONE .AND. BNORM.GT.ONE ) THEN
                        IF( BNORM.GT.BIGNUM*SMINI ) &
                            SCALE = ONE / BNORM
                    END IF
                    TEMP = SCALE / SMINI
                    X( 1, 1 ) = TEMP*B( 1, 1 )
                    X( 2, 1 ) = TEMP*B( 2, 1 )
                    X( 1, 2 ) = TEMP*B( 1, 2 )
                    X( 2, 2 ) = TEMP*B( 2, 2 )
                    XNORM = TEMP*BNORM
                    INFO = 1
                    RETURN
                END IF
                !
                !           Gaussian elimination with complete pivoting.
                !
                UR11 = CRV( ICMAX )
                UI11 = CIV( ICMAX )
                CR21 = CRV( IPIVOT( 2, ICMAX ) )
                CI21 = CIV( IPIVOT( 2, ICMAX ) )
                UR12 = CRV( IPIVOT( 3, ICMAX ) )
                UI12 = CIV( IPIVOT( 3, ICMAX ) )
                CR22 = CRV( IPIVOT( 4, ICMAX ) )
                CI22 = CIV( IPIVOT( 4, ICMAX ) )
                IF( ICMAX.EQ.1 .OR. ICMAX.EQ.4 ) THEN
                    !
                    !              Code when off-diagonals of pivoted C are real
                    !
                    IF( ABS( UR11 ).GT.ABS( UI11 ) ) THEN
                        TEMP = UI11 / UR11
                        UR11R = ONE / ( UR11*( ONE+TEMP**2 ) )
                        UI11R = -TEMP*UR11R
                    ELSE
                        TEMP = UR11 / UI11
                        UI11R = -ONE / ( UI11*( ONE+TEMP**2 ) )
                        UR11R = -TEMP*UI11R
                    END IF
                    LR21 = CR21*UR11R
                    LI21 = CR21*UI11R
                    UR12S = UR12*UR11R
                    UI12S = UR12*UI11R
                    UR22 = CR22 - UR12*LR21
                    UI22 = CI22 - UR12*LI21
                ELSE
                    !
                    !              Code when diagonals of pivoted C are real
                    !
                    UR11R = ONE / UR11
                    UI11R = ZERO
                    LR21 = CR21*UR11R
                    LI21 = CI21*UR11R
                    UR12S = UR12*UR11R
                    UI12S = UI12*UR11R
                    UR22 = CR22 - UR12*LR21 + UI12*LI21
                    UI22 = -UR12*LI21 - UI12*LR21
                END IF
                U22ABS = ABS( UR22 ) + ABS( UI22 )
                !
                !           If smaller pivot < SMINI, use SMINI
                !
                IF( U22ABS.LT.SMINI ) THEN
                    UR22 = SMINI
                    UI22 = ZERO
                    INFO = 1
                END IF
                IF( RSWAP( ICMAX ) ) THEN
                    BR2 = B( 1, 1 )
                    BR1 = B( 2, 1 )
                    BI2 = B( 1, 2 )
                    BI1 = B( 2, 2 )
                ELSE
                    BR1 = B( 1, 1 )
                    BR2 = B( 2, 1 )
                    BI1 = B( 1, 2 )
                    BI2 = B( 2, 2 )
                END IF
                BR2 = BR2 - LR21*BR1 + LI21*BI1
                BI2 = BI2 - LI21*BR1 - LR21*BI1
                BBND = MAX( ( ABS( BR1 )+ABS( BI1 ) )* &
                    ( U22ABS*( ABS( UR11R )+ABS( UI11R ) ) ), &
                    ABS( BR2 )+ABS( BI2 ) )
                IF( BBND.GT.ONE .AND. U22ABS.LT.ONE ) THEN
                    IF( BBND.GE.BIGNUM*U22ABS ) THEN
                        SCALE = ONE / BBND
                        BR1 = SCALE*BR1
                        BI1 = SCALE*BI1
                        BR2 = SCALE*BR2
                        BI2 = SCALE*BI2
                    END IF
                END IF
                !
                CALL SLADIV( BR2, BI2, UR22, UI22, XR2, XI2 )
                XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2
                XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2
                IF( CSWAP( ICMAX ) ) THEN
                    X( 1, 1 ) = XR2
                    X( 2, 1 ) = XR1
                    X( 1, 2 ) = XI2
                    X( 2, 2 ) = XI1
                ELSE
                    X( 1, 1 ) = XR1
                    X( 2, 1 ) = XR2
                    X( 1, 2 ) = XI1
                    X( 2, 2 ) = XI2
                END IF
                XNORM = MAX( ABS( XR1 )+ABS( XI1 ), ABS( XR2 )+ABS( XI2 ) )
                !
                !           Further scaling if  norm(A) norm(X) > overflow
                !
                IF( XNORM.GT.ONE .AND. CMAX.GT.ONE ) THEN
                    IF( XNORM.GT.BIGNUM / CMAX ) THEN
                        TEMP = CMAX / BIGNUM
                        X( 1, 1 ) = TEMP*X( 1, 1 )
                        X( 2, 1 ) = TEMP*X( 2, 1 )
                        X( 1, 2 ) = TEMP*X( 1, 2 )
                        X( 2, 2 ) = TEMP*X( 2, 2 )
                        XNORM = TEMP*XNORM
                        SCALE = TEMP*SCALE
                    END IF
                END IF
            END IF
        END IF
        !
        RETURN
        !
        !     End of SLALN2
        !
    END SUBROUTINE SLALN2

    SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        REAL               A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, HALF, ONE
        PARAMETER          ( ZERO = 0.0E+0, HALF = 0.5E+0, ONE = 1.0E+0 )
        REAL               MULTPL
        PARAMETER          ( MULTPL = 4.0E+0 )
        !     ..
        !     .. Local Scalars ..
        REAL               AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB, &
            SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX, MIN, SIGN, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        EPS = SLAMCH( 'P' )
        IF( C.EQ.ZERO ) THEN
            CS = ONE
            SN = ZERO
            GO TO 10
            !
        ELSE IF( B.EQ.ZERO ) THEN
            !
            !        Swap rows and columns
            !
            CS = ZERO
            SN = ONE
            TEMP = D
            D = A
            A = TEMP
            B = -C
            C = ZERO
            GO TO 10
        ELSE IF( (A-D).EQ.ZERO .AND. SIGN( ONE, B ).NE. &
            SIGN( ONE, C ) ) THEN
            CS = ONE
            SN = ZERO
            GO TO 10
        ELSE
            !
            TEMP = A - D
            P = HALF*TEMP
            BCMAX = MAX( ABS( B ), ABS( C ) )
            BCMIS = MIN( ABS( B ), ABS( C ) )*SIGN( ONE, B )*SIGN( ONE, C )
            SCALE = MAX( ABS( P ), BCMAX )
            Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS
            !
            !        If Z is of the order of the machine accuracy, postpone the
            !        decision on the nature of eigenvalues
            !
            IF( Z.GE.MULTPL*EPS ) THEN
                !
                !           Real eigenvalues. Compute A and D.
                !
                Z = P + SIGN( SQRT( SCALE )*SQRT( Z ), P )
                A = D + Z
                D = D - ( BCMAX / Z )*BCMIS
                !
                !           Compute B and the rotation matrix
                !
                TAU = SLAPY2( C, Z )
                CS = Z / TAU
                SN = C / TAU
                B = B - C
                C = ZERO
            ELSE
                !
                !           Complex eigenvalues, or real (almost) equal eigenvalues.
                !           Make diagonal elements equal.
                !
                SIGMA = B + C
                TAU = SLAPY2( SIGMA, TEMP )
                CS = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
                SN = -( P / ( TAU*CS ) )*SIGN( ONE, SIGMA )
                !
                !           Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
                !                   [ CC  DD ]   [ C  D ] [ SN  CS ]
                !
                AA = A*CS + B*SN
                BB = -A*SN + B*CS
                CC = C*CS + D*SN
                DD = -C*SN + D*CS
                !
                !           Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
                !                   [ C  D ]   [-SN  CS ] [ CC  DD ]
                !
                A = AA*CS + CC*SN
                B = BB*CS + DD*SN
                C = -AA*SN + CC*CS
                D = -BB*SN + DD*CS
                !
                TEMP = HALF*( A+D )
                A = TEMP
                D = TEMP
                !
                IF( C.NE.ZERO ) THEN
                    IF( B.NE.ZERO ) THEN
                        IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN
                            !
                            !                    Real eigenvalues: reduce to upper triangular form
                            !
                            SAB = SQRT( ABS( B ) )
                            SAC = SQRT( ABS( C ) )
                            P = SIGN( SAB*SAC, C )
                            TAU = ONE / SQRT( ABS( B+C ) )
                            A = TEMP + P
                            D = TEMP - P
                            B = B - C
                            C = ZERO
                            CS1 = SAB*TAU
                            SN1 = SAC*TAU
                            TEMP = CS*CS1 - SN*SN1
                            SN = CS*SN1 + SN*CS1
                            CS = TEMP
                        END IF
                    ELSE
                        B = -C
                        C = ZERO
                        TEMP = CS
                        CS = -SN
                        SN = TEMP
                    END IF
                END IF
            END IF
            !
        END IF
        !
10      CONTINUE
        !
        !     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
        !
        RT1R = A
        RT2R = D
        IF( C.EQ.ZERO ) THEN
            RT1I = ZERO
            RT2I = ZERO
        ELSE
            RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
            RT2I = -RT1I
        END IF
        RETURN
        !
        !     End of SLANV2
        !
    END SUBROUTINE SLANV2

    SUBROUTINE SLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
        ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
        LOGICAL            WANTT, WANTZ
        !     ..
        !     .. Array Arguments ..
        REAL               H( LDH, * ), WI( * ), WORK( * ), WR( * ), &
            Z( LDZ, * )
        !     ..
        !
        !  ================================================================
        !     .. Parameters ..
        !
        !     ==== Matrices of order NTINY or smaller must be processed by
        !     .    SLAHQR because of insufficient subdiagonal scratch space.
        !     .    (This is a hard limit.) ====
        INTEGER            NTINY
        PARAMETER          ( NTINY = 11 )
        !
        !     ==== Exceptional deflation windows:  try to cure rare
        !     .    slow convergence by varying the size of the
        !     .    deflation window after KEXNW iterations. ====
        INTEGER            KEXNW
        PARAMETER          ( KEXNW = 5 )
        !
        !     ==== Exceptional shifts: try to cure rare slow convergence
        !     .    with ad-hoc exceptional shifts every KEXSH iterations.
        !     .    ====
        INTEGER            KEXSH
        PARAMETER          ( KEXSH = 6 )
        !
        !     ==== The constants WILK1 and WILK2 are used to form the
        !     .    exceptional shifts. ====
        REAL               WILK1, WILK2
        PARAMETER          ( WILK1 = 0.75e0, WILK2 = -0.4375e0 )
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0e0, ONE = 1.0e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               AA, BB, CC, CS, DD, SN, SS, SWAP
        INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
            KT, KTOP, KU, KV, KWH, KWTOP, KWV, LD, LS, &
            LWKOPT, NDEC, NDFL, NH, NHO, NIBBLE, NMIN, NS, &
            NSMAX, NSR, NVE, NW, NWMAX, NWR, NWUPBD
        LOGICAL            SORTED
        CHARACTER          JBCMPZ*2
        !     ..
        !     .. Local Arrays ..
        REAL               ZDUM( 1, 1 )
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, INT, MAX, MIN, MOD, REAL
        !     ..
        !     .. Executable Statements ..
        INFO = 0
        !
        !     ==== Quick return for N = 0: nothing to do. ====
        !
        IF( N.EQ.0 ) THEN
            WORK( 1 ) = ONE
            RETURN
        END IF
        !
        IF( N.LE.NTINY ) THEN
            !
            !        ==== Tiny matrices must use SLAHQR. ====
            !
            LWKOPT = 1
            IF( LWORK.NE.-1 ) &
                CALL SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                ILOZ, IHIZ, Z, LDZ, INFO )
        ELSE
            !
            !        ==== Use small bulge multi-shift QR with aggressive early
            !        .    deflation on larger-than-tiny matrices. ====
            !
            !        ==== Hope for the best. ====
            !
            INFO = 0
            !
            !        ==== Set up job flags for ILAENV. ====
            !
            IF( WANTT ) THEN
                JBCMPZ( 1: 1 ) = 'S'
            ELSE
                JBCMPZ( 1: 1 ) = 'E'
            END IF
            IF( WANTZ ) THEN
                JBCMPZ( 2: 2 ) = 'V'
            ELSE
                JBCMPZ( 2: 2 ) = 'N'
            END IF
            !
            !        ==== NWR = recommended deflation window size.  At this
            !        .    point,  N .GT. NTINY = 11, so there is enough
            !        .    subdiagonal workspace for NWR.GE.2 as required.
            !        .    (In fact, there is enough subdiagonal space for
            !        .    NWR.GE.3.) ====
            !
            NWR = ILAENV( 13, 'SLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
            NWR = MAX( 2, NWR )
            NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )
            !
            !        ==== NSR = recommended number of simultaneous shifts.
            !        .    At this point N .GT. NTINY = 11, so there is at
            !        .    enough subdiagonal workspace for NSR to be even
            !        .    and greater than or equal to two as required. ====
            !
            NSR = ILAENV( 15, 'SLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
            NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
            NSR = MAX( 2, NSR-MOD( NSR, 2 ) )
            !
            !        ==== Estimate optimal workspace ====
            !
            !        ==== Workspace query call to SLAQR3 ====
            !
            CALL SLAQR3( WANTT, WANTZ, N, ILO, IHI, NWR+1, H, LDH, ILOZ, &
                IHIZ, Z, LDZ, LS, LD, WR, WI, H, LDH, N, H, LDH, &
                N, H, LDH, WORK, -1 )
            !
            !        ==== Optimal workspace = MAX(SLAQR5, SLAQR3) ====
            !
            LWKOPT = MAX( 3*NSR / 2, INT( WORK( 1 ) ) )
            !
            !        ==== Quick return in case of workspace query. ====
            !
            IF( LWORK.EQ.-1 ) THEN
                WORK( 1 ) = REAL( LWKOPT )
                RETURN
            END IF
            !
            !        ==== SLAHQR/SLAQR0 crossover point ====
            !
            NMIN = ILAENV( 12, 'SLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
            NMIN = MAX( NTINY, NMIN )
            !
            !        ==== Nibble crossover point ====
            !
            NIBBLE = ILAENV( 14, 'SLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
            NIBBLE = MAX( 0, NIBBLE )
            !
            !        ==== Accumulate reflections during ttswp?  Use block
            !        .    2-by-2 structure during matrix-matrix multiply? ====
            !
            KACC22 = ILAENV( 16, 'SLAQR0', JBCMPZ, N, ILO, IHI, LWORK )
            KACC22 = MAX( 0, KACC22 )
            KACC22 = MIN( 2, KACC22 )
            !
            !        ==== NWMAX = the largest possible deflation window for
            !        .    which there is sufficient workspace. ====
            !
            NWMAX = MIN( ( N-1 ) / 3, LWORK / 2 )
            NW = NWMAX
            !
            !        ==== NSMAX = the Largest number of simultaneous shifts
            !        .    for which there is sufficient workspace. ====
            !
            NSMAX = MIN( ( N+6 ) / 9, 2*LWORK / 3 )
            NSMAX = NSMAX - MOD( NSMAX, 2 )
            !
            !        ==== NDFL: an iteration count restarted at deflation. ====
            !
            NDFL = 1
            !
            !        ==== ITMAX = iteration limit ====
            !
            ITMAX = MAX( 30, 2*KEXSH )*MAX( 10, ( IHI-ILO+1 ) )
            !
            !        ==== Last row and column in the active block ====
            !
            KBOT = IHI
            !
            !        ==== Main Loop ====
            !
            DO  IT = 1, ITMAX
                !
                !           ==== Done when KBOT falls below ILO ====
                !
                IF( KBOT.LT.ILO ) &
                    GO TO 90
                !
                !           ==== Locate active block ====
                !
                DO  K = KBOT, ILO + 1, -1
                    IF( H( K, K-1 ).EQ.ZERO ) &
                        GO TO 20
                END DO
                K = ILO
20              CONTINUE
                KTOP = K
                !
                !           ==== Select deflation window size:
                !           .    Typical Case:
                !           .      If possible and advisable, nibble the entire
                !           .      active block.  If not, use size MIN(NWR,NWMAX)
                !           .      or MIN(NWR+1,NWMAX) depending upon which has
                !           .      the smaller corresponding subdiagonal entry
                !           .      (a heuristic).
                !           .
                !           .    Exceptional Case:
                !           .      If there have been no deflations in KEXNW or
                !           .      more iterations, then vary the deflation window
                !           .      size.   At first, because, larger windows are,
                !           .      in general, more powerful than smaller ones,
                !           .      rapidly increase the window to the maximum possible.
                !           .      Then, gradually reduce the window size. ====
                !
                NH = KBOT - KTOP + 1
                NWUPBD = MIN( NH, NWMAX )
                IF( NDFL.LT.KEXNW ) THEN
                    NW = MIN( NWUPBD, NWR )
                ELSE
                    NW = MIN( NWUPBD, 2*NW )
                END IF
                IF( NW.LT.NWMAX ) THEN
                    IF( NW.GE.NH-1 ) THEN
                        NW = NH
                    ELSE
                        KWTOP = KBOT - NW + 1
                        IF( ABS( H( KWTOP, KWTOP-1 ) ).GT. &
                            ABS( H( KWTOP-1, KWTOP-2 ) ) )NW = NW + 1
                    END IF
                END IF
                IF( NDFL.LT.KEXNW ) THEN
                    NDEC = -1
                ELSE IF( NDEC.GE.0 .OR. NW.GE.NWUPBD ) THEN
                    NDEC = NDEC + 1
                    IF( NW-NDEC.LT.2 ) &
                        NDEC = 0
                    NW = NW - NDEC
                END IF
                !
                !           ==== Aggressive early deflation:
                !           .    split workspace under the subdiagonal into
                !           .      - an nw-by-nw work array V in the lower
                !           .        left-hand-corner,
                !           .      - an NW-by-at-least-NW-but-more-is-better
                !           .        (NW-by-NHO) horizontal work array along
                !           .        the bottom edge,
                !           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
                !           .        vertical work array along the left-hand-edge.
                !           .        ====
                !
                KV = N - NW + 1
                KT = NW + 1
                NHO = ( N-NW-1 ) - KT + 1
                KWV = NW + 2
                NVE = ( N-NW ) - KWV + 1
                !
                !           ==== Aggressive early deflation ====
                !
                CALL SLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
                    IHIZ, Z, LDZ, LS, LD, WR, WI, H( KV, 1 ), LDH, &
                    NHO, H( KV, KT ), LDH, NVE, H( KWV, 1 ), LDH, &
                    WORK, LWORK )
                !
                !           ==== Adjust KBOT accounting for new deflations. ====
                !
                KBOT = KBOT - LD
                !
                !           ==== KS points to the shifts. ====
                !
                KS = KBOT - LS + 1
                !
                !           ==== Skip an expensive QR sweep if there is a (partly
                !           .    heuristic) reason to expect that many eigenvalues
                !           .    will deflate without it.  Here, the QR sweep is
                !           .    skipped if many eigenvalues have just been deflated
                !           .    or if the remaining active block is small.
                !
                IF( ( LD.EQ.0 ) .OR. ( ( 100*LD.LE.NW*NIBBLE ) .AND. ( KBOT- &
                    KTOP+1.GT.MIN( NMIN, NWMAX ) ) ) ) THEN
                    !
                    !              ==== NS = nominal number of simultaneous shifts.
                    !              .    This may be lowered (slightly) if SLAQR3
                    !              .    did not provide that many shifts. ====
                    !
                    NS = MIN( NSMAX, NSR, MAX( 2, KBOT-KTOP ) )
                    NS = NS - MOD( NS, 2 )
                    !
                    !              ==== If there have been no deflations
                    !              .    in a multiple of KEXSH iterations,
                    !              .    then try exceptional shifts.
                    !              .    Otherwise use shifts provided by
                    !              .    SLAQR3 above or from the eigenvalues
                    !              .    of a trailing principal submatrix. ====
                    !
                    IF( MOD( NDFL, KEXSH ).EQ.0 ) THEN
                        KS = KBOT - NS + 1
                        DO  I = KBOT, MAX( KS+1, KTOP+2 ), -2
                            SS = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
                            AA = WILK1*SS + H( I, I )
                            BB = SS
                            CC = WILK2*SS
                            DD = AA
                            CALL SLANV2( AA, BB, CC, DD, WR( I-1 ), WI( I-1 ), &
                                WR( I ), WI( I ), CS, SN )
                        END DO
                        IF( KS.EQ.KTOP ) THEN
                            WR( KS+1 ) = H( KS+1, KS+1 )
                            WI( KS+1 ) = ZERO
                            WR( KS ) = WR( KS+1 )
                            WI( KS ) = WI( KS+1 )
                        END IF
                    ELSE
                        !
                        !                 ==== Got NS/2 or fewer shifts? Use SLAQR4 or
                        !                 .    SLAHQR on a trailing principal submatrix to
                        !                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
                        !                 .    there is enough space below the subdiagonal
                        !                 .    to fit an NS-by-NS scratch array.) ====
                        !
                        IF( KBOT-KS+1.LE.NS / 2 ) THEN
                            KS = KBOT - NS + 1
                            KT = N - NS + 1
                            CALL SLACPY( 'A', NS, NS, H( KS, KS ), LDH, &
                                H( KT, 1 ), LDH )
                            IF( NS.GT.NMIN ) THEN
                                CALL SLAQR4( .false., .false., NS, 1, NS, &
                                    H( KT, 1 ), LDH, WR( KS ), &
                                    WI( KS ), 1, 1, ZDUM, 1, WORK, &
                                    LWORK, INF )
                            ELSE
                                CALL SLAHQR( .false., .false., NS, 1, NS, &
                                    H( KT, 1 ), LDH, WR( KS ), &
                                    WI( KS ), 1, 1, ZDUM, 1, INF )
                            END IF
                            KS = KS + INF
                            !
                            !                    ==== In case of a rare QR failure use
                            !                    .    eigenvalues of the trailing 2-by-2
                            !                    .    principal submatrix.  ====
                            !
                            IF( KS.GE.KBOT ) THEN
                                AA = H( KBOT-1, KBOT-1 )
                                CC = H( KBOT, KBOT-1 )
                                BB = H( KBOT-1, KBOT )
                                DD = H( KBOT, KBOT )
                                CALL SLANV2( AA, BB, CC, DD, WR( KBOT-1 ), &
                                    WI( KBOT-1 ), WR( KBOT ), &
                                    WI( KBOT ), CS, SN )
                                KS = KBOT - 1
                            END IF
                        END IF
                        !
                        IF( KBOT-KS+1.GT.NS ) THEN
                            !
                            !                    ==== Sort the shifts (Helps a little)
                            !                    .    Bubble sort keeps complex conjugate
                            !                    .    pairs together. ====
                            !
                            SORTED = .false.
                            DO  K = KBOT, KS + 1, -1
                                IF( SORTED ) &
                                    GO TO 60
                                SORTED = .true.
                                DO  I = KS, K - 1
                                    IF( ABS( WR( I ) )+ABS( WI( I ) ).LT. &
                                        ABS( WR( I+1 ) )+ABS( WI( I+1 ) ) ) THEN
                                        SORTED = .false.
                                        !
                                        SWAP = WR( I )
                                        WR( I ) = WR( I+1 )
                                        WR( I+1 ) = SWAP
                                        !
                                        SWAP = WI( I )
                                        WI( I ) = WI( I+1 )
                                        WI( I+1 ) = SWAP
                                    END IF
                                END DO
                            END DO
60                          CONTINUE
                        END IF
                        !
                        !                 ==== Shuffle shifts into pairs of real shifts
                        !                 .    and pairs of complex conjugate shifts
                        !                 .    assuming complex conjugate shifts are
                        !                 .    already adjacent to one another. (Yes,
                        !                 .    they are.)  ====
                        !
                        DO  I = KBOT, KS + 2, -2
                            IF( WI( I ).NE.-WI( I-1 ) ) THEN
                                !
                                SWAP = WR( I )
                                WR( I ) = WR( I-1 )
                                WR( I-1 ) = WR( I-2 )
                                WR( I-2 ) = SWAP
                                !
                                SWAP = WI( I )
                                WI( I ) = WI( I-1 )
                                WI( I-1 ) = WI( I-2 )
                                WI( I-2 ) = SWAP
                            END IF
                        END DO
                    END IF
                    !
                    !              ==== If there are only two shifts and both are
                    !              .    real, then use only one.  ====
                    !
                    IF( KBOT-KS+1.EQ.2 ) THEN
                        IF( WI( KBOT ).EQ.ZERO ) THEN
                            IF( ABS( WR( KBOT )-H( KBOT, KBOT ) ).LT. &
                                ABS( WR( KBOT-1 )-H( KBOT, KBOT ) ) ) THEN
                                WR( KBOT-1 ) = WR( KBOT )
                            ELSE
                                WR( KBOT ) = WR( KBOT-1 )
                            END IF
                        END IF
                    END IF
                    !
                    !              ==== Use up to NS of the the smallest magnatiude
                    !              .    shifts.  If there aren't NS shifts available,
                    !              .    then use them all, possibly dropping one to
                    !              .    make the number of shifts even. ====
                    !
                    NS = MIN( NS, KBOT-KS+1 )
                    NS = NS - MOD( NS, 2 )
                    KS = KBOT - NS + 1
                    !
                    !              ==== Small-bulge multi-shift QR sweep:
                    !              .    split workspace under the subdiagonal into
                    !              .    - a KDU-by-KDU work array U in the lower
                    !              .      left-hand-corner,
                    !              .    - a KDU-by-at-least-KDU-but-more-is-better
                    !              .      (KDU-by-NHo) horizontal work array WH along
                    !              .      the bottom edge,
                    !              .    - and an at-least-KDU-but-more-is-better-by-KDU
                    !              .      (NVE-by-KDU) vertical work WV arrow along
                    !              .      the left-hand-edge. ====
                    !
                    KDU = 3*NS - 3
                    KU = N - KDU + 1
                    KWH = KDU + 1
                    NHO = ( N-KDU+1-4 ) - ( KDU+1 ) + 1
                    KWV = KDU + 4
                    NVE = N - KDU - KWV + 1
                    !
                    !              ==== Small-bulge multi-shift QR sweep ====
                    !
                    CALL SLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NS, &
                        WR( KS ), WI( KS ), H, LDH, ILOZ, IHIZ, Z, &
                        LDZ, WORK, 3, H( KU, 1 ), LDH, NVE, &
                        H( KWV, 1 ), LDH, NHO, H( KU, KWH ), LDH )
                END IF
                !
                !           ==== Note progress (or the lack of it). ====
                !
                IF( LD.GT.0 ) THEN
                    NDFL = 1
                ELSE
                    NDFL = NDFL + 1
                END IF
                !
                !           ==== End of main loop ====
            END DO
            !
            !        ==== Iteration limit exceeded.  Set INFO to show where
            !        .    the problem occurred and exit. ====
            !
            INFO = KBOT
90          CONTINUE
        END IF
        !
        !     ==== Return the optimal value of LWORK. ====
        !
        WORK( 1 ) = REAL( LWKOPT )
        !
        !     ==== End of SLAQR0 ====
        !
    END SUBROUTINE SLAQR0

    SUBROUTINE SLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
        IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, &
        LDT, NV, WV, LDWV, WORK, LWORK )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
            LDZ, LWORK, N, ND, NH, NS, NV, NW
        LOGICAL            WANTT, WANTZ
        !     ..
        !     .. Array Arguments ..
        REAL               H( LDH, * ), SI( * ), SR( * ), T( LDT, * ), &
            V( LDV, * ), WORK( * ), WV( LDWV, * ), &
            Z( LDZ, * )
        !     ..
        !
        !  ================================================================
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0e0, ONE = 1.0e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, &
            SAFMAX, SAFMIN, SMLNUM, SN, TAU, ULP
        INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, &
            KEND, KLN, KROW, KWTOP, LTOP, LWK1, LWK2, LWK3, &
            LWKOPT, NMIN
        LOGICAL            BULGE, SORTED
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, INT, MAX, MIN, REAL, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     ==== Estimate optimal workspace. ====
        !
        JW = MIN( NW, KBOT-KTOP+1 )
        IF( JW.LE.2 ) THEN
            LWKOPT = 1
        ELSE
            !
            !        ==== Workspace query call to SGEHRD ====
            !
            CALL SGEHRD( JW, 1, JW-1, T, LDT, WORK, WORK, -1, INFO )
            LWK1 = INT( WORK( 1 ) )
            !
            !        ==== Workspace query call to SORMHR ====
            !
            CALL SORMHR( 'R', 'N', JW, JW, 1, JW-1, T, LDT, WORK, V, LDV, &
                WORK, -1, INFO )
            LWK2 = INT( WORK( 1 ) )
            !
            !        ==== Workspace query call to SLAQR4 ====
            !
            CALL SLAQR4( .true., .true., JW, 1, JW, T, LDT, SR, SI, 1, JW, &
                V, LDV, WORK, -1, INFQR )
            LWK3 = INT( WORK( 1 ) )
            !
            !        ==== Optimal workspace ====
            !
            LWKOPT = MAX( JW+MAX( LWK1, LWK2 ), LWK3 )
        END IF
        !
        !     ==== Quick return in case of workspace query. ====
        !
        IF( LWORK.EQ.-1 ) THEN
            WORK( 1 ) = REAL( LWKOPT )
            RETURN
        END IF
        !
        !     ==== Nothing to do ...
        !     ... for an empty active block ... ====
        NS = 0
        ND = 0
        WORK( 1 ) = ONE
        IF( KTOP.GT.KBOT ) &
            RETURN
        !     ... nor for an empty deflation window. ====
        IF( NW.LT.1 ) &
            RETURN
        !
        !     ==== Machine constants ====
        !
        SAFMIN = SLAMCH( 'SAFE MINIMUM' )
        SAFMAX = ONE / SAFMIN
        CALL SLABAD( SAFMIN, SAFMAX )
        ULP = SLAMCH( 'PRECISION' )
        SMLNUM = SAFMIN*( REAL( N ) / ULP )
        !
        !     ==== Setup deflation window ====
        !
        JW = MIN( NW, KBOT-KTOP+1 )
        KWTOP = KBOT - JW + 1
        IF( KWTOP.EQ.KTOP ) THEN
            S = ZERO
        ELSE
            S = H( KWTOP, KWTOP-1 )
        END IF
        !
        IF( KBOT.EQ.KWTOP ) THEN
            !
            !        ==== 1-by-1 deflation window: not much to do ====
            !
            SR( KWTOP ) = H( KWTOP, KWTOP )
            SI( KWTOP ) = ZERO
            NS = 1
            ND = 0
            IF( ABS( S ).LE.MAX( SMLNUM, ULP*ABS( H( KWTOP, KWTOP ) ) ) ) &
                THEN
                NS = 0
                ND = 1
                IF( KWTOP.GT.KTOP ) &
                    H( KWTOP, KWTOP-1 ) = ZERO
            END IF
            WORK( 1 ) = ONE
            RETURN
        END IF
        !
        !     ==== Convert to spike-triangular form.  (In case of a
        !     .    rare QR failure, this routine continues to do
        !     .    aggressive early deflation using that part of
        !     .    the deflation window that converged using INFQR
        !     .    here and there to keep track.) ====
        !
        CALL SLACPY( 'U', JW, JW, H( KWTOP, KWTOP ), LDH, T, LDT )
        CALL SCOPY( JW-1, H( KWTOP+1, KWTOP ), LDH+1, T( 2, 1 ), LDT+1 )
        !
        CALL SLASET( 'A', JW, JW, ZERO, ONE, V, LDV )
        NMIN = ILAENV( 12, 'SLAQR3', 'SV', JW, 1, JW, LWORK )
        IF( JW.GT.NMIN ) THEN
            CALL SLAQR4( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ), &
                SI( KWTOP ), 1, JW, V, LDV, WORK, LWORK, INFQR )
        ELSE
            CALL SLAHQR( .true., .true., JW, 1, JW, T, LDT, SR( KWTOP ), &
                SI( KWTOP ), 1, JW, V, LDV, INFQR )
        END IF
        !
        !     ==== STREXC needs a clean margin near the diagonal ====
        !
        DO  J = 1, JW - 3
            T( J+2, J ) = ZERO
            T( J+3, J ) = ZERO
        END DO
        IF( JW.GT.2 ) &
            T( JW, JW-2 ) = ZERO
        !
        !     ==== Deflation detection loop ====
        !
        NS = JW
        ILST = INFQR + 1
20      CONTINUE
        IF( ILST.LE.NS ) THEN
            IF( NS.EQ.1 ) THEN
                BULGE = .FALSE.
            ELSE
                BULGE = T( NS, NS-1 ).NE.ZERO
            END IF
            !
            !        ==== Small spike tip test for deflation ====
            !
            IF( .NOT. BULGE ) THEN
                !
                !           ==== Real eigenvalue ====
                !
                FOO = ABS( T( NS, NS ) )
                IF( FOO.EQ.ZERO ) &
                    FOO = ABS( S )
                IF( ABS( S*V( 1, NS ) ).LE.MAX( SMLNUM, ULP*FOO ) ) THEN
                    !
                    !              ==== Deflatable ====
                    !
                    NS = NS - 1
                ELSE
                    !
                    !              ==== Undeflatable.   Move it up out of the way.
                    !              .    (STREXC can not fail in this case.) ====
                    !
                    IFST = NS
                    CALL STREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                        INFO )
                    ILST = ILST + 1
                END IF
            ELSE
                !
                !           ==== Complex conjugate pair ====
                !
                FOO = ABS( T( NS, NS ) ) + SQRT( ABS( T( NS, NS-1 ) ) )* &
                    SQRT( ABS( T( NS-1, NS ) ) )
                IF( FOO.EQ.ZERO ) &
                    FOO = ABS( S )
                IF( MAX( ABS( S*V( 1, NS ) ), ABS( S*V( 1, NS-1 ) ) ).LE. &
                    MAX( SMLNUM, ULP*FOO ) ) THEN
                    !
                    !              ==== Deflatable ====
                    !
                    NS = NS - 2
                ELSE
                    !
                    !              ==== Undeflatable. Move them up out of the way.
                    !              .    Fortunately, STREXC does the right thing with
                    !              .    ILST in case of a rare exchange failure. ====
                    !
                    IFST = NS
                    CALL STREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                        INFO )
                    ILST = ILST + 2
                END IF
            END IF
            !
            !        ==== End deflation detection loop ====
            !
            GO TO 20
        END IF
        !
        !        ==== Return to Hessenberg form ====
        !
        IF( NS.EQ.0 ) &
            S = ZERO
        !
        IF( NS.LT.JW ) THEN
            !
            !        ==== sorting diagonal blocks of T improves accuracy for
            !        .    graded matrices.  Bubble sort deals well with
            !        .    exchange failures. ====
            !
            SORTED = .false.
            I = NS + 1
30          CONTINUE
            IF( SORTED ) &
                GO TO 50
            SORTED = .true.
            !
            KEND = I - 1
            I = INFQR + 1
            IF( I.EQ.NS ) THEN
                K = I + 1
            ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
                K = I + 1
            ELSE
                K = I + 2
            END IF
40          CONTINUE
            IF( K.LE.KEND ) THEN
                IF( K.EQ.I+1 ) THEN
                    EVI = ABS( T( I, I ) )
                ELSE
                    EVI = ABS( T( I, I ) ) + SQRT( ABS( T( I+1, I ) ) )* &
                        SQRT( ABS( T( I, I+1 ) ) )
                END IF
                !
                IF( K.EQ.KEND ) THEN
                    EVK = ABS( T( K, K ) )
                ELSE IF( T( K+1, K ).EQ.ZERO ) THEN
                    EVK = ABS( T( K, K ) )
                ELSE
                    EVK = ABS( T( K, K ) ) + SQRT( ABS( T( K+1, K ) ) )* &
                        SQRT( ABS( T( K, K+1 ) ) )
                END IF
                !
                IF( EVI.GE.EVK ) THEN
                    I = K
                ELSE
                    SORTED = .false.
                    IFST = I
                    ILST = K
                    CALL STREXC( 'V', JW, T, LDT, V, LDV, IFST, ILST, WORK, &
                        INFO )
                    IF( INFO.EQ.0 ) THEN
                        I = ILST
                    ELSE
                        I = K
                    END IF
                END IF
                IF( I.EQ.KEND ) THEN
                    K = I + 1
                ELSE IF( T( I+1, I ).EQ.ZERO ) THEN
                    K = I + 1
                ELSE
                    K = I + 2
                END IF
                GO TO 40
            END IF
            GO TO 30
50          CONTINUE
        END IF
        !
        !     ==== Restore shift/eigenvalue array from T ====
        !
        I = JW
60      CONTINUE
        IF( I.GE.INFQR+1 ) THEN
            IF( I.EQ.INFQR+1 ) THEN
                SR( KWTOP+I-1 ) = T( I, I )
                SI( KWTOP+I-1 ) = ZERO
                I = I - 1
            ELSE IF( T( I, I-1 ).EQ.ZERO ) THEN
                SR( KWTOP+I-1 ) = T( I, I )
                SI( KWTOP+I-1 ) = ZERO
                I = I - 1
            ELSE
                AA = T( I-1, I-1 )
                CC = T( I, I-1 )
                BB = T( I-1, I )
                DD = T( I, I )
                CALL SLANV2( AA, BB, CC, DD, SR( KWTOP+I-2 ), &
                    SI( KWTOP+I-2 ), SR( KWTOP+I-1 ), &
                    SI( KWTOP+I-1 ), CS, SN )
                I = I - 2
            END IF
            GO TO 60
        END IF
        !
        IF( NS.LT.JW .OR. S.EQ.ZERO ) THEN
            IF( NS.GT.1 .AND. S.NE.ZERO ) THEN
                !
                !           ==== Reflect spike back into lower triangle ====
                !
                CALL SCOPY( NS, V, LDV, WORK, 1 )
                BETA = WORK( 1 )
                CALL SLARFG( NS, BETA, WORK( 2 ), 1, TAU )
                WORK( 1 ) = ONE
                !
                CALL SLASET( 'L', JW-2, JW-2, ZERO, ZERO, T( 3, 1 ), LDT )
                !
                CALL SLARF( 'L', NS, JW, WORK, 1, TAU, T, LDT, &
                    WORK( JW+1 ) )
                CALL SLARF( 'R', NS, NS, WORK, 1, TAU, T, LDT, &
                    WORK( JW+1 ) )
                CALL SLARF( 'R', JW, NS, WORK, 1, TAU, V, LDV, &
                    WORK( JW+1 ) )
                !
                CALL SGEHRD( JW, 1, NS, T, LDT, WORK, WORK( JW+1 ), &
                    LWORK-JW, INFO )
            END IF
            !
            !        ==== Copy updated reduced window into place ====
            !
            IF( KWTOP.GT.1 ) &
                H( KWTOP, KWTOP-1 ) = S*V( 1, 1 )
            CALL SLACPY( 'U', JW, JW, T, LDT, H( KWTOP, KWTOP ), LDH )
            CALL SCOPY( JW-1, T( 2, 1 ), LDT+1, H( KWTOP+1, KWTOP ), &
                LDH+1 )
            !
            !        ==== Accumulate orthogonal matrix in order update
            !        .    H and Z, if requested.  ====
            !
            IF( NS.GT.1 .AND. S.NE.ZERO ) &
                CALL SORMHR( 'R', 'N', JW, NS, 1, NS, T, LDT, WORK, V, LDV, &
                WORK( JW+1 ), LWORK-JW, INFO )
            !
            !        ==== Update vertical slab in H ====
            !
            IF( WANTT ) THEN
                LTOP = 1
            ELSE
                LTOP = KTOP
            END IF
            DO  KROW = LTOP, KWTOP - 1, NV
                KLN = MIN( NV, KWTOP-KROW )
                CALL SGEMM( 'N', 'N', KLN, JW, JW, ONE, H( KROW, KWTOP ), &
                    LDH, V, LDV, ZERO, WV, LDWV )
                CALL SLACPY( 'A', KLN, JW, WV, LDWV, H( KROW, KWTOP ), LDH )
            END DO
            !
            !        ==== Update horizontal slab in H ====
            !
            IF( WANTT ) THEN
                DO  KCOL = KBOT + 1, N, NH
                    KLN = MIN( NH, N-KCOL+1 )
                    CALL SGEMM( 'C', 'N', JW, KLN, JW, ONE, V, LDV, &
                        H( KWTOP, KCOL ), LDH, ZERO, T, LDT )
                    CALL SLACPY( 'A', JW, KLN, T, LDT, H( KWTOP, KCOL ), &
                        LDH )
                END DO
            END IF
            !
            !        ==== Update vertical slab in Z ====
            !
            IF( WANTZ ) THEN
                DO  KROW = ILOZ, IHIZ, NV
                    KLN = MIN( NV, IHIZ-KROW+1 )
                    CALL SGEMM( 'N', 'N', KLN, JW, JW, ONE, Z( KROW, KWTOP ), &
                        LDZ, V, LDV, ZERO, WV, LDWV )
                    CALL SLACPY( 'A', KLN, JW, WV, LDWV, Z( KROW, KWTOP ), &
                        LDZ )
                END DO
            END IF
        END IF
        !
        !     ==== Return the number of deflations ... ====
        !
        ND = JW - NS
        !
        !     ==== ... and the number of shifts. (Subtracting
        !     .    INFQR from the spike length takes care
        !     .    of the case of a rare QR failure while
        !     .    calculating eigenvalues of the deflation
        !     .    window.)  ====
        !
        NS = NS - INFQR
        !
        !      ==== Return optimal workspace. ====
        !
        WORK( 1 ) = REAL( LWKOPT )
        !
        !     ==== End of SLAQR3 ====
        !
    END SUBROUTINE SLAQR3

    SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
        T, LDT, C, LDC, WORK, LDWORK )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2013
        !
        !     .. Scalar Arguments ..
        CHARACTER          DIRECT, SIDE, STOREV, TRANS
        INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ), &
            WORK( LDWORK, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE
        PARAMETER          ( ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        CHARACTER          TRANST
        INTEGER            I, J
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick return if possible
        !
        IF( M.LE.0 .OR. N.LE.0 ) &
            RETURN
        !
        IF( LSAME( TRANS, 'N' ) ) THEN
            TRANST = 'T'
        ELSE
            TRANST = 'N'
        END IF
        !
        IF( LSAME( STOREV, 'C' ) ) THEN
            !
            IF( LSAME( DIRECT, 'F' ) ) THEN
                !
                !           Let  V =  ( V1 )    (first K rows)
                !                     ( V2 )
                !           where  V1  is unit lower triangular.
                !
                IF( LSAME( SIDE, 'L' ) ) THEN
                    !
                    !              Form  H * C  or  H**T * C  where  C = ( C1 )
                    !                                                    ( C2 )
                    !
                    !              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
                    !
                    !              W := C1**T
                    !
                    DO  J = 1, K
                        CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V1
                    !
                    CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                        K, ONE, V, LDV, WORK, LDWORK )
                    IF( M.GT.K ) THEN
                        !
                        !                 W := W + C2**T * V2
                        !
                        CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                            ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T**T  or  W * T
                    !
                    CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - V * W**T
                    !
                    IF( M.GT.K ) THEN
                        !
                        !                 C2 := C2 - V2 * W**T
                        !
                        CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                            -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, &
                            C( K+1, 1 ), LDC )
                    END IF
                    !
                    !              W := W * V1**T
                    !
                    CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                        ONE, V, LDV, WORK, LDWORK )
                    !
                    !              C1 := C1 - W**T
                    !
                    DO  J = 1, K
                        DO  I = 1, N
                            C( J, I ) = C( J, I ) - WORK( I, J )
                        END DO
                    END DO
                    !
                ELSE IF( LSAME( SIDE, 'R' ) ) THEN
                    !
                    !              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
                    !
                    !              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
                    !
                    !              W := C1
                    !
                    DO  J = 1, K
                        CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V1
                    !
                    CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                        K, ONE, V, LDV, WORK, LDWORK )
                    IF( N.GT.K ) THEN
                        !
                        !                 W := W + C2 * V2
                        !
                        CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                            ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T  or  W * T**T
                    !
                    CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - W * V**T
                    !
                    IF( N.GT.K ) THEN
                        !
                        !                 C2 := C2 - W * V2**T
                        !
                        CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                            C( 1, K+1 ), LDC )
                    END IF
                    !
                    !              W := W * V1**T
                    !
                    CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                        ONE, V, LDV, WORK, LDWORK )
                    !
                    !              C1 := C1 - W
                    !
                    DO  J = 1, K
                        DO  I = 1, M
                            C( I, J ) = C( I, J ) - WORK( I, J )
                        END DO
                    END DO
                END IF
                !
            ELSE
                !
                !           Let  V =  ( V1 )
                !                     ( V2 )    (last K rows)
                !           where  V2  is unit upper triangular.
                !
                IF( LSAME( SIDE, 'L' ) ) THEN
                    !
                    !              Form  H * C  or  H**T * C  where  C = ( C1 )
                    !                                                    ( C2 )
                    !
                    !              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
                    !
                    !              W := C2**T
                    !
                    DO  J = 1, K
                        CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V2
                    !
                    CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                        K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
                    IF( M.GT.K ) THEN
                        !
                        !                 W := W + C1**T * V1
                        !
                        CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T**T  or  W * T
                    !
                    CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - V * W**T
                    !
                    IF( M.GT.K ) THEN
                        !
                        !                 C1 := C1 - V1 * W**T
                        !
                        CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                            -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
                    END IF
                    !
                    !              W := W * V2**T
                    !
                    CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                        ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
                    !
                    !              C2 := C2 - W**T
                    !
                    DO  J = 1, K
                        DO  I = 1, N
                            C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
                        END DO
                    END DO
                    !
                ELSE IF( LSAME( SIDE, 'R' ) ) THEN
                    !
                    !              Form  C * H  or  C * H'  where  C = ( C1  C2 )
                    !
                    !              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
                    !
                    !              W := C2
                    !
                    DO  J = 1, K
                        CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V2
                    !
                    CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                        K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
                    IF( N.GT.K ) THEN
                        !
                        !                 W := W + C1 * V1
                        !
                        CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T  or  W * T**T
                    !
                    CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - W * V**T
                    !
                    IF( N.GT.K ) THEN
                        !
                        !                 C1 := C1 - W * V1**T
                        !
                        CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
                    END IF
                    !
                    !              W := W * V2**T
                    !
                    CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                        ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
                    !
                    !              C2 := C2 - W
                    !
                    DO  J = 1, K
                        DO  I = 1, M
                            C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
                        END DO
                    END DO
                END IF
            END IF
            !
        ELSE IF( LSAME( STOREV, 'R' ) ) THEN
            !
            IF( LSAME( DIRECT, 'F' ) ) THEN
                !
                !           Let  V =  ( V1  V2 )    (V1: first K columns)
                !           where  V1  is unit upper triangular.
                !
                IF( LSAME( SIDE, 'L' ) ) THEN
                    !
                    !              Form  H * C  or  H**T * C  where  C = ( C1 )
                    !                                                    ( C2 )
                    !
                    !              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
                    !
                    !              W := C1**T
                    !
                    DO  J = 1, K
                        CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V1**T
                    !
                    CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                        ONE, V, LDV, WORK, LDWORK )
                    IF( M.GT.K ) THEN
                        !
                        !                 W := W + C2**T * V2**T
                        !
                        CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                            C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, &
                            WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T**T  or  W * T
                    !
                    CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - V**T * W**T
                    !
                    IF( M.GT.K ) THEN
                        !
                        !                 C2 := C2 - V2**T * W**T
                        !
                        CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                            V( 1, K+1 ), LDV, WORK, LDWORK, ONE, &
                            C( K+1, 1 ), LDC )
                    END IF
                    !
                    !              W := W * V1
                    !
                    CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                        K, ONE, V, LDV, WORK, LDWORK )
                    !
                    !              C1 := C1 - W**T
                    !
                    DO  J = 1, K
                        DO  I = 1, N
                            C( J, I ) = C( J, I ) - WORK( I, J )
                        END DO
                    END DO
                    !
                ELSE IF( LSAME( SIDE, 'R' ) ) THEN
                    !
                    !              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
                    !
                    !              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
                    !
                    !              W := C1
                    !
                    DO  J = 1, K
                        CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V1**T
                    !
                    CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                        ONE, V, LDV, WORK, LDWORK )
                    IF( N.GT.K ) THEN
                        !
                        !                 W := W + C2 * V2**T
                        !
                        CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                            ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, &
                            ONE, WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T  or  W * T**T
                    !
                    CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - W * V
                    !
                    IF( N.GT.K ) THEN
                        !
                        !                 C2 := C2 - W * V2
                        !
                        CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, &
                            C( 1, K+1 ), LDC )
                    END IF
                    !
                    !              W := W * V1
                    !
                    CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                        K, ONE, V, LDV, WORK, LDWORK )
                    !
                    !              C1 := C1 - W
                    !
                    DO  J = 1, K
                        DO  I = 1, M
                            C( I, J ) = C( I, J ) - WORK( I, J )
                        END DO
                    END DO
                    !
                END IF
                !
            ELSE
                !
                !           Let  V =  ( V1  V2 )    (V2: last K columns)
                !           where  V2  is unit lower triangular.
                !
                IF( LSAME( SIDE, 'L' ) ) THEN
                    !
                    !              Form  H * C  or  H**T * C  where  C = ( C1 )
                    !                                                    ( C2 )
                    !
                    !              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
                    !
                    !              W := C2**T
                    !
                    DO  J = 1, K
                        CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V2**T
                    !
                    CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
                        ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
                    IF( M.GT.K ) THEN
                        !
                        !                 W := W + C1**T * V1**T
                        !
                        CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                            C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T**T  or  W * T
                    !
                    CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - V**T * W**T
                    !
                    IF( M.GT.K ) THEN
                        !
                        !                 C1 := C1 - V1**T * W**T
                        !
                        CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                            V, LDV, WORK, LDWORK, ONE, C, LDC )
                    END IF
                    !
                    !              W := W * V2
                    !
                    CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                        K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
                    !
                    !              C2 := C2 - W**T
                    !
                    DO  J = 1, K
                        DO  I = 1, N
                            C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
                        END DO
                    END DO
                    !
                ELSE IF( LSAME( SIDE, 'R' ) ) THEN
                    !
                    !              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
                    !
                    !              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
                    !
                    !              W := C2
                    !
                    DO  J = 1, K
                        CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
                    END DO
                    !
                    !              W := W * V2**T
                    !
                    CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                        ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
                    IF( N.GT.K ) THEN
                        !
                        !                 W := W + C1 * V1**T
                        !
                        CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                            ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
                    END IF
                    !
                    !              W := W * T  or  W * T**T
                    !
                    CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
                        ONE, T, LDT, WORK, LDWORK )
                    !
                    !              C := C - W * V
                    !
                    IF( N.GT.K ) THEN
                        !
                        !                 C1 := C1 - W * V1
                        !
                        CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
                    END IF
                    !
                    !              W := W * V2
                    !
                    CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
                        K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
                    !
                    !              C1 := C1 - W
                    !
                    DO  J = 1, K
                        DO  I = 1, M
                            C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
                        END DO
                    END DO
                    !
                END IF
                !
            END IF
        END IF
        !
        RETURN
        !
        !     End of SLARFB
        !
    END SUBROUTINE SLARFB

    SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          SIDE
        INTEGER            INCV, LDC, M, N
        REAL               TAU
        !     ..
        !     .. Array Arguments ..
        REAL               C( LDC, * ), V( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE, ZERO
        PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            APPLYLEFT
        INTEGER            I, LASTV, LASTC
        !     ..
        !     .. Executable Statements ..
        !
        APPLYLEFT = LSAME( SIDE, 'L' )
        LASTV = 0
        LASTC = 0
        IF( TAU.NE.ZERO ) THEN
            !     Set up variables for scanning V.  LASTV begins pointing to the end
            !     of V.
            IF( APPLYLEFT ) THEN
                LASTV = M
            ELSE
                LASTV = N
            END IF
            IF( INCV.GT.0 ) THEN
                I = 1 + (LASTV-1) * INCV
            ELSE
                I = 1
            END IF
            !     Look for the last non-zero row in V.
            DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
                LASTV = LASTV - 1
                I = I - INCV
            END DO
            IF( APPLYLEFT ) THEN
                !     Scan for the last non-zero column in C(1:lastv,:).
                LASTC = ILASLC(LASTV, N, C, LDC)
            ELSE
                !     Scan for the last non-zero row in C(:,1:lastv).
                LASTC = ILASLR(M, LASTV, C, LDC)
            END IF
        END IF
        !     Note that lastc.eq.0 renders the BLAS operations null; no special
        !     case is needed at this level.
        IF( APPLYLEFT ) THEN
            !
            !        Form  H * C
            !
            IF( LASTV.GT.0 ) THEN
                !
                !           w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
                !
                CALL SGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, &
                    ZERO, WORK, 1 )
                !
                !           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
                !
                CALL SGER( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
            END IF
        ELSE
            !
            !        Form  C * H
            !
            IF( LASTV.GT.0 ) THEN
                !
                !           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
                !
                CALL SGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, &
                    V, INCV, ZERO, WORK, 1 )
                !
                !           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
                !
                CALL SGER( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
            END IF
        END IF
        RETURN
        !
        !     End of SLARF
        !
    END SUBROUTINE SLARF

    SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )
        !
        !  -- LAPACK auxiliary routine (version 3.8.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        INTEGER            INCX, N
        REAL               ALPHA, TAU
        !     ..
        !     .. Array Arguments ..
        REAL               X( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE, ZERO
        PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            J, KNT
        REAL               BETA, RSAFMN, SAFMIN, XNORM
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, SIGN
        !     ..
        !     .. Executable Statements ..
        !
        IF( N.LE.1 ) THEN
            TAU = ZERO
            RETURN
        END IF
        !
        XNORM = SNRM2( N-1, X, INCX )
        !
        IF( XNORM.EQ.ZERO ) THEN
            !
            !        H  =  I
            !
            TAU = ZERO
        ELSE
            !
            !        general case
            !
            BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
            SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
            KNT = 0
            IF( ABS( BETA ).LT.SAFMIN ) THEN
                !
                !           XNORM, BETA may be inaccurate; scale X and recompute them
                !
                RSAFMN = ONE / SAFMIN
10              CONTINUE
                KNT = KNT + 1
                CALL SSCAL( N-1, RSAFMN, X, INCX )
                BETA = BETA*RSAFMN
                ALPHA = ALPHA*RSAFMN
                IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) &
                    GO TO 10
                !
                !           New BETA is at most 1, at least SAFMIN
                !
                XNORM = SNRM2( N-1, X, INCX )
                BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
            END IF
            TAU = ( BETA-ALPHA ) / BETA
            CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            !
            !        If ALPHA is subnormal, it may lose relative accuracy
            !
            DO  J = 1, KNT
                BETA = BETA*SAFMIN
            END DO
            ALPHA = BETA
        END IF
        !
        RETURN
        !
        !     End of SLARFG
        !
    END SUBROUTINE SLARFG

    SUBROUTINE SLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          DIRECT, STOREV
        INTEGER            K, LDT, LDV, N
        !     ..
        !     .. Array Arguments ..
        REAL               T( LDT, * ), TAU( * ), V( LDV, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE, ZERO
        PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            I, J, PREVLASTV, LASTV
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick return if possible
        !
        IF( N.EQ.0 ) &
            RETURN
        !
        IF( LSAME( DIRECT, 'F' ) ) THEN
            PREVLASTV = N
            DO I = 1, K
                PREVLASTV = MAX( I, PREVLASTV )
                IF( TAU( I ).EQ.ZERO ) THEN
                    !
                    !              H(i)  =  I
                    !
                    DO J = 1, I
                        T( J, I ) = ZERO
                    END DO
                ELSE
                    !
                    !              general case
                    !
                    IF( LSAME( STOREV, 'C' ) ) THEN
                        !                 Skip any trailing zeros.
                        DO LASTV = N, I+1, -1
                            IF( V( LASTV, I ).NE.ZERO ) EXIT
                        END DO
                        DO J = 1, I-1
                            T( J, I ) = -TAU( I ) * V( I , J )
                        END DO
                        J = MIN( LASTV, PREVLASTV )
                        !
                        !                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**T * V(i:j,i)
                        !
                        CALL SGEMV( 'Transpose', J-I, I-1, -TAU( I ), &
                            V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE, &
                            T( 1, I ), 1 )
                    ELSE
                        !                 Skip any trailing zeros.
                        DO LASTV = N, I+1, -1
                            IF( V( I, LASTV ).NE.ZERO ) EXIT
                        END DO
                        DO J = 1, I-1
                            T( J, I ) = -TAU( I ) * V( J , I )
                        END DO
                        J = MIN( LASTV, PREVLASTV )
                        !
                        !                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**T
                        !
                        CALL SGEMV( 'No transpose', I-1, J-I, -TAU( I ), &
                            V( 1, I+1 ), LDV, V( I, I+1 ), LDV, &
                            ONE, T( 1, I ), 1 )
                    END IF
                    !
                    !              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
                    !
                    CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
                        LDT, T( 1, I ), 1 )
                    T( I, I ) = TAU( I )
                    IF( I.GT.1 ) THEN
                        PREVLASTV = MAX( PREVLASTV, LASTV )
                    ELSE
                        PREVLASTV = LASTV
                    END IF
                END IF
            END DO
        ELSE
            PREVLASTV = 1
            DO I = K, 1, -1
                IF( TAU( I ).EQ.ZERO ) THEN
                    !
                    !              H(i)  =  I
                    !
                    DO J = I, K
                        T( J, I ) = ZERO
                    END DO
                ELSE
                    !
                    !              general case
                    !
                    IF( I.LT.K ) THEN
                        IF( LSAME( STOREV, 'C' ) ) THEN
                            !                    Skip any leading zeros.
                            DO LASTV = 1, I-1
                                IF( V( LASTV, I ).NE.ZERO ) EXIT
                            END DO
                            DO J = I+1, K
                                T( J, I ) = -TAU( I ) * V( N-K+I , J )
                            END DO
                            J = MAX( LASTV, PREVLASTV )
                            !
                            !                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
                            !
                            CALL SGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ), &
                                V( J, I+1 ), LDV, V( J, I ), 1, ONE, &
                                T( I+1, I ), 1 )
                        ELSE
                            !                    Skip any leading zeros.
                            DO LASTV = 1, I-1
                                IF( V( I, LASTV ).NE.ZERO ) EXIT
                            END DO
                            DO J = I+1, K
                                T( J, I ) = -TAU( I ) * V( J, N-K+I )
                            END DO
                            J = MAX( LASTV, PREVLASTV )
                            !
                            !                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**T
                            !
                            CALL SGEMV( 'No transpose', K-I, N-K+I-J, &
                                -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, &
                                ONE, T( I+1, I ), 1 )
                        END IF
                        !
                        !                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
                        !
                        CALL STRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                            T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
                        IF( I.GT.1 ) THEN
                            PREVLASTV = MIN( PREVLASTV, LASTV )
                        ELSE
                            PREVLASTV = LASTV
                        END IF
                    END IF
                    T( I, I ) = TAU( I )
                END IF
            END DO
        END IF
        RETURN
        !
        !     End of SLARFT
        !
    END SUBROUTINE SLARFT

    SUBROUTINE SLARTG( F, G, CS, SN, R )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        REAL               CS, F, G, R, SN
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO
        PARAMETER          ( ZERO = 0.0E0 )
        REAL               ONE
        PARAMETER          ( ONE = 1.0E0 )
        REAL               TWO
        PARAMETER          ( TWO = 2.0E0 )
        !     ..
        !     .. Local Scalars ..
        !     LOGICAL            FIRST
        INTEGER            COUNT, I
        REAL               EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, INT, LOG, MAX, SQRT
        !     ..
        !     .. Save statement ..
        !     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
        !     ..
        !     .. Data statements ..
        !     DATA               FIRST / .TRUE. /
        !     ..
        !     .. Executable Statements ..
        !
        !     IF( FIRST ) THEN
        SAFMIN = SLAMCH( 'S' )
        EPS = SLAMCH( 'E' )
        SAFMN2 = SLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
            LOG( SLAMCH( 'B' ) ) / TWO )
        SAFMX2 = ONE / SAFMN2
        !        FIRST = .FALSE.
        !     END IF
        IF( G.EQ.ZERO ) THEN
            CS = ONE
            SN = ZERO
            R = F
        ELSE IF( F.EQ.ZERO ) THEN
            CS = ZERO
            SN = ONE
            R = G
        ELSE
            F1 = F
            G1 = G
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE.GE.SAFMX2 ) THEN
                COUNT = 0
10              CONTINUE
                COUNT = COUNT + 1
                F1 = F1*SAFMN2
                G1 = G1*SAFMN2
                SCALE = MAX( ABS( F1 ), ABS( G1 ) )
                IF( SCALE.GE.SAFMX2 ) &
                    GO TO 10
                R = SQRT( F1**2+G1**2 )
                CS = F1 / R
                SN = G1 / R
                DO  I = 1, COUNT
                    R = R*SAFMX2
                END DO
            ELSE IF( SCALE.LE.SAFMN2 ) THEN
                COUNT = 0
30              CONTINUE
                COUNT = COUNT + 1
                F1 = F1*SAFMX2
                G1 = G1*SAFMX2
                SCALE = MAX( ABS( F1 ), ABS( G1 ) )
                IF( SCALE.LE.SAFMN2 ) &
                    GO TO 30
                R = SQRT( F1**2+G1**2 )
                CS = F1 / R
                SN = G1 / R
                DO  I = 1, COUNT
                    R = R*SAFMN2
                END DO
            ELSE
                R = SQRT( F1**2+G1**2 )
                CS = F1 / R
                SN = G1 / R
            END IF
            IF( ABS( F ).GT.ABS( G ) .AND. CS.LT.ZERO ) THEN
                CS = -CS
                SN = -SN
                R = -R
            END IF
        END IF
        RETURN
        !
        !     End of SLARTG
        !
    END SUBROUTINE SLARTG

    SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          TYPE
        INTEGER            INFO, KL, KU, LDA, M, N
        REAL               CFROM, CTO
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            DONE
        INTEGER            I, ITYPE, J, K1, K2, K3, K4
        REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        INFO = 0
        !
        IF( LSAME( TYPE, 'G' ) ) THEN
            ITYPE = 0
        ELSE IF( LSAME( TYPE, 'L' ) ) THEN
            ITYPE = 1
        ELSE IF( LSAME( TYPE, 'U' ) ) THEN
            ITYPE = 2
        ELSE IF( LSAME( TYPE, 'H' ) ) THEN
            ITYPE = 3
        ELSE IF( LSAME( TYPE, 'B' ) ) THEN
            ITYPE = 4
        ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
            ITYPE = 5
        ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
            ITYPE = 6
        ELSE
            ITYPE = -1
        END IF
        !
        IF( ITYPE.EQ.-1 ) THEN
            INFO = -1
        ELSE IF( CFROM.EQ.ZERO .OR. SISNAN(CFROM) ) THEN
            INFO = -4
        ELSE IF( SISNAN(CTO) ) THEN
            INFO = -5
        ELSE IF( M.LT.0 ) THEN
            INFO = -6
        ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. &
            ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
            INFO = -7
        ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
            INFO = -9
        ELSE IF( ITYPE.GE.4 ) THEN
            IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
                INFO = -2
            ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. &
                ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) &
                THEN
                INFO = -3
            ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. &
                ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. &
                ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
                INFO = -9
            END IF
        END IF
        !
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SLASCL', -INFO )
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.EQ.0 .OR. M.EQ.0 ) &
            RETURN
        !
        !     Get machine parameters
        !
        SMLNUM = SLAMCH( 'S' )
        BIGNUM = ONE / SMLNUM
        !
        CFROMC = CFROM
        CTOC = CTO
        !
10      CONTINUE
        CFROM1 = CFROMC*SMLNUM
        IF( CFROM1.EQ.CFROMC ) THEN
            !        CFROMC is an inf.  Multiply by a correctly signed zero for
            !        finite CTOC, or a NaN if CTOC is infinite.
            MUL = CTOC / CFROMC
            DONE = .TRUE.
            CTO1 = CTOC
        ELSE
            CTO1 = CTOC / BIGNUM
            IF( CTO1.EQ.CTOC ) THEN
                !           CTOC is either 0 or an inf.  In both cases, CTOC itself
                !           serves as the correct multiplication factor.
                MUL = CTOC
                DONE = .TRUE.
                CFROMC = ONE
            ELSE IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
                MUL = SMLNUM
                DONE = .FALSE.
                CFROMC = CFROM1
            ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
                MUL = BIGNUM
                DONE = .FALSE.
                CTOC = CTO1
            ELSE
                MUL = CTOC / CFROMC
                DONE = .TRUE.
            END IF
        END IF
        !
        IF( ITYPE.EQ.0 ) THEN
            !
            !        Full matrix
            !
            DO  J = 1, N
                DO  I = 1, M
                    A( I, J ) = A( I, J )*MUL
                END DO
            END DO
            !
        ELSE IF( ITYPE.EQ.1 ) THEN
            !
            !        Lower triangular matrix
            !
            DO  J = 1, N
                DO  I = J, M
                    A( I, J ) = A( I, J )*MUL
                END DO
            END DO
            !
        ELSE IF( ITYPE.EQ.2 ) THEN
            !
            !        Upper triangular matrix
            !
            DO  J = 1, N
                DO  I = 1, MIN( J, M )
                    A( I, J ) = A( I, J )*MUL
                END DO
            END DO
            !
        ELSE IF( ITYPE.EQ.3 ) THEN
            !
            !        Upper Hessenberg matrix
            !
            DO  J = 1, N
                DO  I = 1, MIN( J+1, M )
                    A( I, J ) = A( I, J )*MUL
                END DO
            END DO
            !
        ELSE IF( ITYPE.EQ.4 ) THEN
            !
            !        Lower half of a symmetric band matrix
            !
            K3 = KL + 1
            K4 = N + 1
            DO  J = 1, N
                DO  I = 1, MIN( K3, K4-J )
                    A( I, J ) = A( I, J )*MUL
                END DO
            END DO
            !
        ELSE IF( ITYPE.EQ.5 ) THEN
            !
            !        Upper half of a symmetric band matrix
            !
            K1 = KU + 2
            K3 = KU + 1
            DO  J = 1, N
                DO  I = MAX( K1-J, 1 ), K3
                    A( I, J ) = A( I, J )*MUL
                END DO
            END DO
            !
        ELSE IF( ITYPE.EQ.6 ) THEN
            !
            !        Band matrix
            !
            K1 = KL + KU + 2
            K2 = KL + 1
            K3 = 2*KL + KU + 1
            K4 = KL + KU + 1 + M
            DO  J = 1, N
                DO  I = MAX( K1-J, K2 ), MIN( K3, K4-J )
                    A( I, J ) = A( I, J )*MUL
                END DO
            END DO
            !
        END IF
        !
        IF( .NOT.DONE ) &
            GO TO 10
        !
        RETURN
        !
        !     End of SLASCL
        !
    END SUBROUTINE SLASCL

    SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          UPLO
        INTEGER            LDA, M, N
        REAL               ALPHA, BETA
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. Local Scalars ..
        INTEGER            I, J
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MIN
        !     ..
        !     .. Executable Statements ..
        !
        IF( LSAME( UPLO, 'U' ) ) THEN
            !
            !        Set the strictly upper triangular or trapezoidal part of the
            !        array to ALPHA.
            !
            DO  J = 2, N
                DO  I = 1, MIN( J-1, M )
                    A( I, J ) = ALPHA
                END DO
            END DO
            !
        ELSE IF( LSAME( UPLO, 'L' ) ) THEN
            !
            !        Set the strictly lower triangular or trapezoidal part of the
            !        array to ALPHA.
            !
            DO  J = 1, MIN( M, N )
                DO  I = J + 1, M
                    A( I, J ) = ALPHA
                END DO
            END DO
            !
        ELSE
            !
            !        Set the leading m-by-n submatrix to ALPHA.
            !
            DO  J = 1, N
                DO  I = 1, M
                    A( I, J ) = ALPHA
                END DO
            END DO
        END IF
        !
        !     Set the first min(M,N) diagonal elements to BETA.
        !
        DO  I = 1, MIN( M, N )
            A( I, I ) = BETA
        END DO
        !
        RETURN
        !
        !     End of SLASET
        !
    END SUBROUTINE SLASET

    SUBROUTINE SORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            INFO, K, LDA, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE, ZERO
        PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            I, J, L
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        INFO = 0
        IF( M.LT.0 ) THEN
            INFO = -1
        ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
            INFO = -2
        ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
            INFO = -5
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SORG2R', -INFO )
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.LE.0 ) &
            RETURN
        !
        !     Initialise columns k+1:n to columns of the unit matrix
        !
        DO  J = K + 1, N
            DO  L = 1, M
                A( L, J ) = ZERO
            END DO
            A( J, J ) = ONE
        END DO
        !
        DO  I = K, 1, -1
            !
            !        Apply H(i) to A(i:m,i:n) from the left
            !
            IF( I.LT.N ) THEN
                A( I, I ) = ONE
                CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                    A( I, I+1 ), LDA, WORK )
            END IF
            IF( I.LT.M ) &
                CALL SSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
            A( I, I ) = ONE - TAU( I )
            !
            !        Set A(1:i-1,i) to zero
            !
            DO  L = 1, I - 1
                A( L, I ) = ZERO
            END DO
        END DO
        RETURN
        !
        !     End of SORG2R
        !
    END SUBROUTINE SORG2R

    SUBROUTINE SORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, ILO, INFO, LDA, LWORK, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            LQUERY
        INTEGER            I, IINFO, J, LWKOPT, NB, NH
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        INFO = 0
        NH = IHI - ILO
        LQUERY = ( LWORK.EQ.-1 )
        IF( N.LT.0 ) THEN
            INFO = -1
        ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
            INFO = -2
        ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        ELSE IF( LWORK.LT.MAX( 1, NH ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
        END IF
        !
        IF( INFO.EQ.0 ) THEN
            NB = ILAENV( 1, 'SORGQR', ' ', NH, NH, NH, -1 )
            LWKOPT = MAX( 1, NH )*NB
            WORK( 1 ) = LWKOPT
        END IF
        !
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SORGHR', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.EQ.0 ) THEN
            WORK( 1 ) = 1
            RETURN
        END IF
        !
        !     Shift the vectors which define the elementary reflectors one
        !     column to the right, and set the first ilo and the last n-ihi
        !     rows and columns to those of the unit matrix
        !
        DO  J = IHI, ILO + 1, -1
            DO  I = 1, J - 1
                A( I, J ) = ZERO
            END DO
            DO  I = J + 1, IHI
                A( I, J ) = A( I, J-1 )
            END DO
            DO  I = IHI + 1, N
                A( I, J ) = ZERO
            END DO
        END DO
        DO  J = 1, ILO
            DO  I = 1, N
                A( I, J ) = ZERO
            END DO
            A( J, J ) = ONE
        END DO
        DO  J = IHI + 1, N
            DO  I = 1, N
                A( I, J ) = ZERO
            END DO
            A( J, J ) = ONE
        END DO
        !
        IF( NH.GT.0 ) THEN
            !
            !        Generate Q(ilo+1:ihi,ilo+1:ihi)
            !
            CALL SORGQR( NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), &
                WORK, LWORK, IINFO )
        END IF
        WORK( 1 ) = LWKOPT
        RETURN
        !
        !     End of SORGHR
        !
    END SUBROUTINE SORGHR

    SUBROUTINE SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            INFO, K, LDA, LWORK, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO
        PARAMETER          ( ZERO = 0.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            LQUERY
        INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
            LWKOPT, NB, NBMIN, NX
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        INFO = 0
        NB = ILAENV( 1, 'SORGQR', ' ', M, N, K, -1 )
        LWKOPT = MAX( 1, N )*NB
        WORK( 1 ) = LWKOPT
        LQUERY = ( LWORK.EQ.-1 )
        IF( M.LT.0 ) THEN
            INFO = -1
        ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
            INFO = -2
        ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
            INFO = -5
        ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SORGQR', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.LE.0 ) THEN
            WORK( 1 ) = 1
            RETURN
        END IF
        !
        NBMIN = 2
        NX = 0
        IWS = N
        IF( NB.GT.1 .AND. NB.LT.K ) THEN
            !
            !        Determine when to cross over from blocked to unblocked code.
            !
            NX = MAX( 0, ILAENV( 3, 'SORGQR', ' ', M, N, K, -1 ) )
            IF( NX.LT.K ) THEN
                !
                !           Determine if workspace is large enough for blocked code.
                !
                LDWORK = N
                IWS = LDWORK*NB
                IF( LWORK.LT.IWS ) THEN
                    !
                    !              Not enough workspace to use optimal NB:  reduce NB and
                    !              determine the minimum value of NB.
                    !
                    NB = LWORK / LDWORK
                    NBMIN = MAX( 2, ILAENV( 2, 'SORGQR', ' ', M, N, K, -1 ) )
                END IF
            END IF
        END IF
        !
        IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
            !
            !        Use blocked code after the last block.
            !        The first kk columns are handled by the block method.
            !
            KI = ( ( K-NX-1 ) / NB )*NB
            KK = MIN( K, KI+NB )
            !
            !        Set A(1:kk,kk+1:n) to zero.
            !
            DO  J = KK + 1, N
                DO  I = 1, KK
                    A( I, J ) = ZERO
                END DO
            END DO
        ELSE
            KK = 0
        END IF
        !
        !     Use unblocked code for the last or only block.
        !
        IF( KK.LT.N ) &
            CALL SORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
            TAU( KK+1 ), WORK, IINFO )
        !
        IF( KK.GT.0 ) THEN
            !
            !        Use blocked code
            !
            DO  I = KI + 1, 1, -NB
                IB = MIN( NB, K-I+1 )
                IF( I+IB.LE.N ) THEN
                    !
                    !              Form the triangular factor of the block reflector
                    !              H = H(i) H(i+1) . . . H(i+ib-1)
                    !
                    CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                        A( I, I ), LDA, TAU( I ), WORK, LDWORK )
                    !
                    !              Apply H to A(i:m,i+ib:n) from the left
                    !
                    CALL SLARFB( 'Left', 'No transpose', 'Forward', &
                        'Columnwise', M-I+1, N-I-IB+1, IB, &
                        A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                        LDA, WORK( IB+1 ), LDWORK )
                END IF
                !
                !           Apply H to rows i:m of current block
                !
                CALL SORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, &
                    IINFO )
                !
                !           Set rows 1:i-1 of current block to zero
                !
                DO  J = I, I + IB - 1
                    DO  L = 1, I - 1
                        A( L, J ) = ZERO
                    END DO
                END DO
            END DO
        END IF
        !
        WORK( 1 ) = IWS
        RETURN
        !
        !     End of SORGQR
        !
    END SUBROUTINE SORGQR

    SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        REAL C,S
        INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
        REAL SX(*),SY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        REAL STEMP
        INTEGER I,IX,IY
        !     ..
        IF (N.LE.0) RETURN
        IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
            !
            !       code for both increments equal to 1
            !
            DO I = 1,N
                STEMP = C*SX(I) + S*SY(I)
                SY(I) = C*SY(I) - S*SX(I)
                SX(I) = STEMP
            END DO
        ELSE
            !
            !       code for unequal increments or equal increments not equal
            !         to 1
            !
            IX = 1
            IY = 1
            IF (INCX.LT.0) IX = (-N+1)*INCX + 1
            IF (INCY.LT.0) IY = (-N+1)*INCY + 1
            DO I = 1,N
                STEMP = C*SX(IX) + S*SY(IY)
                SY(IY) = C*SY(IY) - S*SX(IX)
                SX(IX) = STEMP
                IX = IX + INCX
                IY = IY + INCY
            END DO
        END IF
        RETURN
    END SUBROUTINE SROT

    SUBROUTINE SSCAL(N,SA,SX,INCX)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        REAL SA
        INTEGER INCX,N
        !     ..
        !     .. Array Arguments ..
        REAL SX(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        INTEGER I,M,MP1,NINCX
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MOD
        !     ..
        IF (N.LE.0 .OR. INCX.LE.0) RETURN
        IF (INCX.EQ.1) THEN
            !
            !        code for increment equal to 1
            !
            !
            !        clean-up loop
            !
            M = MOD(N,5)
            IF (M.NE.0) THEN
                DO I = 1,M
                    SX(I) = SA*SX(I)
                END DO
                IF (N.LT.5) RETURN
            END IF
            MP1 = M + 1
            DO I = MP1,N,5
                SX(I) = SA*SX(I)
                SX(I+1) = SA*SX(I+1)
                SX(I+2) = SA*SX(I+2)
                SX(I+3) = SA*SX(I+3)
                SX(I+4) = SA*SX(I+4)
            END DO
        ELSE
            !
            !        code for increment not equal to 1
            !
            NINCX = N*INCX
            DO I = 1,NINCX,INCX
                SX(I) = SA*SX(I)
            END DO
        END IF
        RETURN
    END SUBROUTINE SSCAL

    SUBROUTINE SSWAP(N,SX,INCX,SY,INCY)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
        REAL SX(*),SY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        REAL STEMP
        INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MOD
        !     ..
        IF (N.LE.0) RETURN
        IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
            !
            !       code for both increments equal to 1
            !
            !
            !       clean-up loop
            !
            M = MOD(N,3)
            IF (M.NE.0) THEN
                DO I = 1,M
                    STEMP = SX(I)
                    SX(I) = SY(I)
                    SY(I) = STEMP
                END DO
                IF (N.LT.3) RETURN
            END IF
            MP1 = M + 1
            DO I = MP1,N,3
                STEMP = SX(I)
                SX(I) = SY(I)
                SY(I) = STEMP
                STEMP = SX(I+1)
                SX(I+1) = SY(I+1)
                SY(I+1) = STEMP
                STEMP = SX(I+2)
                SX(I+2) = SY(I+2)
                SY(I+2) = STEMP
            END DO
        ELSE
            !
            !       code for unequal increments or equal increments not equal
            !         to 1
            !
            IX = 1
            IY = 1
            IF (INCX.LT.0) IX = (-N+1)*INCX + 1
            IF (INCY.LT.0) IY = (-N+1)*INCY + 1
            DO I = 1,N
                STEMP = SX(IX)
                SX(IX) = SY(IY)
                SY(IY) = STEMP
                IX = IX + INCX
                IY = IY + INCY
            END DO
        END IF
        RETURN
    END SUBROUTINE SSWAP

    SUBROUTINE STREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, &
        VR, LDVR, MM, M, WORK, LWORK, INFO )
        IMPLICIT NONE
        !
        !  -- LAPACK computational routine (version 3.8.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        CHARACTER          HOWMNY, SIDE
        INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
        !     ..
        !     .. Array Arguments ..
        LOGICAL            SELECT( * )
        REAL   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
            WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL   ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
        INTEGER            NBMIN, NBMAX
        PARAMETER          ( NBMIN = 8, NBMAX = 128 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, PAIR, &
            RIGHTV, SOMEV
        INTEGER            I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, &
            IV, MAXWRK, NB, KI2
        REAL   BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE, &
            SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR, &
            XNORM
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX, SQRT
        !     ..
        !     .. Local Arrays ..
        REAL   X( 2, 2 )
        INTEGER            ISCOMPLEX( NBMAX )
        !     ..
        !     .. Executable Statements ..
        !
        !     Decode and test the input parameters
        !
        BOTHV  = LSAME( SIDE, 'B' )
        RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
        LEFTV  = LSAME( SIDE, 'L' ) .OR. BOTHV
        !
        ALLV  = LSAME( HOWMNY, 'A' )
        OVER  = LSAME( HOWMNY, 'B' )
        SOMEV = LSAME( HOWMNY, 'S' )
        !
        INFO = 0
        NB = ILAENV( 1, 'STREVC', SIDE // HOWMNY, N, -1, -1, -1 )
        MAXWRK = N + 2*N*NB
        WORK(1) = MAXWRK
        LQUERY = ( LWORK.EQ.-1 )
        IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
            INFO = -1
        ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -4
        ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
            INFO = -6
        ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
            INFO = -8
        ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
            INFO = -10
        ELSE IF( LWORK.LT.MAX( 1, 3*N ) .AND. .NOT.LQUERY ) THEN
            INFO = -14
        ELSE
            !
            !        Set M to the number of columns required to store the selected
            !        eigenvectors, standardize the array SELECT if necessary, and
            !        test MM.
            !
            IF( SOMEV ) THEN
                M = 0
                PAIR = .FALSE.
                DO  J = 1, N
                    IF( PAIR ) THEN
                        PAIR = .FALSE.
                        SELECT( J ) = .FALSE.
                    ELSE
                        IF( J.LT.N ) THEN
                            IF( T( J+1, J ).EQ.ZERO ) THEN
                                IF( SELECT( J ) ) &
                                    M = M + 1
                            ELSE
                                PAIR = .TRUE.
                                IF( SELECT( J ) .OR. SELECT( J+1 ) ) THEN
                                    SELECT( J ) = .TRUE.
                                    M = M + 2
                                END IF
                            END IF
                        ELSE
                            IF( SELECT( N ) ) &
                                M = M + 1
                        END IF
                    END IF
                END DO
            ELSE
                M = N
            END IF
            !
            IF( MM.LT.M ) THEN
                INFO = -11
            END IF
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'STREVC3', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF( N.EQ.0 ) &
            RETURN
        !
        !     Use blocked version of back-transformation if sufficient workspace.
        !     Zero-out the workspace to avoid potential NaN propagation.
        !
        IF( OVER .AND. LWORK .GE. N + 2*N*NBMIN ) THEN
            NB = (LWORK - N) / (2*N)
            NB = MIN( NB, NBMAX )
            CALL SLASET( 'F', N, 1+2*NB, ZERO, ZERO, WORK, N )
        ELSE
            NB = 1
        END IF
        !
        !     Set the constants to control overflow.
        !
        UNFL = SLAMCH( 'Safe minimum' )
        OVFL = ONE / UNFL
        CALL SLABAD( UNFL, OVFL )
        ULP = SLAMCH( 'Precision' )
        SMLNUM = UNFL*( N / ULP )
        BIGNUM = ( ONE-ULP ) / SMLNUM
        !
        !     Compute 1-norm of each column of strictly upper triangular
        !     part of T to control overflow in triangular solver.
        !
        WORK( 1 ) = ZERO
        DO  J = 2, N
            WORK( J ) = ZERO
            DO  I = 1, J - 1
                WORK( J ) = WORK( J ) + ABS( T( I, J ) )
            END DO
        END DO
        !
        !     Index IP is used to specify the real or complex eigenvalue:
        !       IP = 0, real eigenvalue,
        !            1, first  of conjugate complex pair: (wr,wi)
        !           -1, second of conjugate complex pair: (wr,wi)
        !       ISCOMPLEX array stores IP for each column in current block.
        !
        IF( RIGHTV ) THEN
            !
            !        ============================================================
            !        Compute right eigenvectors.
            !
            !        IV is index of column in current block.
            !        For complex right vector, uses IV-1 for real part and IV for complex part.
            !        Non-blocked version always uses IV=2;
            !        blocked     version starts with IV=NB, goes down to 1 or 2.
            !        (Note the "0-th" column is used for 1-norms computed above.)
            IV = 2
            IF( NB.GT.2 ) THEN
                IV = NB
            END IF

            IP = 0
            IS = M
            DO  KI = N, 1, -1
                IF( IP.EQ.-1 ) THEN
                    !              previous iteration (ki+1) was second of conjugate pair,
                    !              so this ki is first of conjugate pair; skip to end of loop
                    IP = 1
                    GO TO 140
                ELSE IF( KI.EQ.1 ) THEN
                    !              last column, so this ki must be real eigenvalue
                    IP = 0
                ELSE IF( T( KI, KI-1 ).EQ.ZERO ) THEN
                    !              zero on sub-diagonal, so this ki is real eigenvalue
                    IP = 0
                ELSE
                    !              non-zero on sub-diagonal, so this ki is second of conjugate pair
                    IP = -1
                END IF

                IF( SOMEV ) THEN
                    IF( IP.EQ.0 ) THEN
                        IF( .NOT.SELECT( KI ) ) &
                            GO TO 140
                    ELSE
                        IF( .NOT.SELECT( KI-1 ) ) &
                            GO TO 140
                    END IF
                END IF
                !
                !           Compute the KI-th eigenvalue (WR,WI).
                !
                WR = T( KI, KI )
                WI = ZERO
                IF( IP.NE.0 ) &
                    WI = SQRT( ABS( T( KI, KI-1 ) ) )* &
                    SQRT( ABS( T( KI-1, KI ) ) )
                SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
                !
                IF( IP.EQ.0 ) THEN
                    !
                    !              --------------------------------------------------------
                    !              Real right eigenvector
                    !
                    WORK( KI + IV*N ) = ONE
                    !
                    !              Form right-hand side.
                    !
                    DO  K = 1, KI - 1
                        WORK( K + IV*N ) = -T( K, KI )
                    END DO
                    !
                    !              Solve upper quasi-triangular system:
                    !              [ T(1:KI-1,1:KI-1) - WR ]*X = SCALE*WORK.
                    !
                    JNXT = KI - 1
                    DO  J = KI - 1, 1, -1
                        IF( J.GT.JNXT ) &
                            GO TO 60
                        J1 = J
                        J2 = J
                        JNXT = J - 1
                        IF( J.GT.1 ) THEN
                            IF( T( J, J-1 ).NE.ZERO ) THEN
                                J1   = J - 1
                                JNXT = J - 2
                            END IF
                        END IF
                        !
                        IF( J1.EQ.J2 ) THEN
                            !
                            !                    1-by-1 diagonal block
                            !
                            CALL SLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                LDT, ONE, ONE, WORK( J+IV*N ), N, WR, &
                                ZERO, X, 2, SCALE, XNORM, IERR )
                            !
                            !                    Scale X(1,1) to avoid overflow when updating
                            !                    the right-hand side.
                            !
                            IF( XNORM.GT.ONE ) THEN
                                IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                                    X( 1, 1 ) = X( 1, 1 ) / XNORM
                                    SCALE = SCALE / XNORM
                                END IF
                            END IF
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) &
                                CALL SSCAL( KI, SCALE, WORK( 1+IV*N ), 1 )
                            WORK( J+IV*N ) = X( 1, 1 )
                            !
                            !                    Update right-hand side
                            !
                            CALL SAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                WORK( 1+IV*N ), 1 )
                            !
                        ELSE
                            !
                            !                    2-by-2 diagonal block
                            !
                            CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, &
                                T( J-1, J-1 ), LDT, ONE, ONE, &
                                WORK( J-1+IV*N ), N, WR, ZERO, X, 2, &
                                SCALE, XNORM, IERR )
                            !
                            !                    Scale X(1,1) and X(2,1) to avoid overflow when
                            !                    updating the right-hand side.
                            !
                            IF( XNORM.GT.ONE ) THEN
                                BETA = MAX( WORK( J-1 ), WORK( J ) )
                                IF( BETA.GT.BIGNUM / XNORM ) THEN
                                    X( 1, 1 ) = X( 1, 1 ) / XNORM
                                    X( 2, 1 ) = X( 2, 1 ) / XNORM
                                    SCALE = SCALE / XNORM
                                END IF
                            END IF
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) &
                                CALL SSCAL( KI, SCALE, WORK( 1+IV*N ), 1 )
                            WORK( J-1+IV*N ) = X( 1, 1 )
                            WORK( J  +IV*N ) = X( 2, 1 )
                            !
                            !                    Update right-hand side
                            !
                            CALL SAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                WORK( 1+IV*N ), 1 )
                            CALL SAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                WORK( 1+IV*N ), 1 )
                        END IF
60                  END DO
                    !
                    !              Copy the vector x or Q*x to VR and normalize.
                    !
                    IF( .NOT.OVER ) THEN
                        !                 ------------------------------
                        !                 no back-transform: copy x to VR and normalize.
                        CALL SCOPY( KI, WORK( 1 + IV*N ), 1, VR( 1, IS ), 1 )
                        !
                        II = ISAMAX( KI, VR( 1, IS ), 1 )
                        REMAX = ONE / ABS( VR( II, IS ) )
                        CALL SSCAL( KI, REMAX, VR( 1, IS ), 1 )
                        !
                        DO  K = KI + 1, N
                            VR( K, IS ) = ZERO
                        END DO
                        !
                    ELSE IF( NB.EQ.1 ) THEN
                        !                 ------------------------------
                        !                 version 1: back-transform each vector with GEMV, Q*x.
                        IF( KI.GT.1 ) &
                            CALL SGEMV( 'N', N, KI-1, ONE, VR, LDVR, &
                            WORK( 1 + IV*N ), 1, WORK( KI + IV*N ), &
                            VR( 1, KI ), 1 )
                        !
                        II = ISAMAX( N, VR( 1, KI ), 1 )
                        REMAX = ONE / ABS( VR( II, KI ) )
                        CALL SSCAL( N, REMAX, VR( 1, KI ), 1 )
                        !
                    ELSE
                        !                 ------------------------------
                        !                 version 2: back-transform block of vectors with GEMM
                        !                 zero out below vector
                        DO K = KI + 1, N
                            WORK( K + IV*N ) = ZERO
                        END DO
                        ISCOMPLEX( IV ) = IP
                        !                 back-transform and normalization is done below
                    END IF
                ELSE
                    !
                    !              --------------------------------------------------------
                    !              Complex right eigenvector.
                    !
                    !              Initial solve
                    !              [ ( T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I*WI) ]*X = 0.
                    !              [ ( T(KI,  KI-1) T(KI,  KI) )               ]
                    !
                    IF( ABS( T( KI-1, KI ) ).GE.ABS( T( KI, KI-1 ) ) ) THEN
                        WORK( KI-1 + (IV-1)*N ) = ONE
                        WORK( KI   + (IV  )*N ) = WI / T( KI-1, KI )
                    ELSE
                        WORK( KI-1 + (IV-1)*N ) = -WI / T( KI, KI-1 )
                        WORK( KI   + (IV  )*N ) = ONE
                    END IF
                    WORK( KI   + (IV-1)*N ) = ZERO
                    WORK( KI-1 + (IV  )*N ) = ZERO
                    !
                    !              Form right-hand side.
                    !
                    DO  K = 1, KI - 2
                        WORK( K+(IV-1)*N ) = -WORK( KI-1+(IV-1)*N )*T(K,KI-1)
                        WORK( K+(IV  )*N ) = -WORK( KI  +(IV  )*N )*T(K,KI  )
                    END DO
                    !
                    !              Solve upper quasi-triangular system:
                    !              [ T(1:KI-2,1:KI-2) - (WR+i*WI) ]*X = SCALE*(WORK+i*WORK2)
                    !
                    JNXT = KI - 2
                    DO  J = KI - 2, 1, -1
                        IF( J.GT.JNXT ) &
                            GO TO 90
                        J1 = J
                        J2 = J
                        JNXT = J - 1
                        IF( J.GT.1 ) THEN
                            IF( T( J, J-1 ).NE.ZERO ) THEN
                                J1   = J - 1
                                JNXT = J - 2
                            END IF
                        END IF
                        !
                        IF( J1.EQ.J2 ) THEN
                            !
                            !                    1-by-1 diagonal block
                            !
                            CALL SLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                LDT, ONE, ONE, WORK( J+(IV-1)*N ), N, &
                                WR, WI, X, 2, SCALE, XNORM, IERR )
                            !
                            !                    Scale X(1,1) and X(1,2) to avoid overflow when
                            !                    updating the right-hand side.
                            !
                            IF( XNORM.GT.ONE ) THEN
                                IF( WORK( J ).GT.BIGNUM / XNORM ) THEN
                                    X( 1, 1 ) = X( 1, 1 ) / XNORM
                                    X( 1, 2 ) = X( 1, 2 ) / XNORM
                                    SCALE = SCALE / XNORM
                                END IF
                            END IF
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) THEN
                                CALL SSCAL( KI, SCALE, WORK( 1+(IV-1)*N ), 1 )
                                CALL SSCAL( KI, SCALE, WORK( 1+(IV  )*N ), 1 )
                            END IF
                            WORK( J+(IV-1)*N ) = X( 1, 1 )
                            WORK( J+(IV  )*N ) = X( 1, 2 )
                            !
                            !                    Update the right-hand side
                            !
                            CALL SAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                WORK( 1+(IV-1)*N ), 1 )
                            CALL SAXPY( J-1, -X( 1, 2 ), T( 1, J ), 1, &
                                WORK( 1+(IV  )*N ), 1 )
                            !
                        ELSE
                            !
                            !                    2-by-2 diagonal block
                            !
                            CALL SLALN2( .FALSE., 2, 2, SMIN, ONE, &
                                T( J-1, J-1 ), LDT, ONE, ONE, &
                                WORK( J-1+(IV-1)*N ), N, WR, WI, X, 2, &
                                SCALE, XNORM, IERR )
                            !
                            !                    Scale X to avoid overflow when updating
                            !                    the right-hand side.
                            !
                            IF( XNORM.GT.ONE ) THEN
                                BETA = MAX( WORK( J-1 ), WORK( J ) )
                                IF( BETA.GT.BIGNUM / XNORM ) THEN
                                    REC = ONE / XNORM
                                    X( 1, 1 ) = X( 1, 1 )*REC
                                    X( 1, 2 ) = X( 1, 2 )*REC
                                    X( 2, 1 ) = X( 2, 1 )*REC
                                    X( 2, 2 ) = X( 2, 2 )*REC
                                    SCALE = SCALE*REC
                                END IF
                            END IF
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) THEN
                                CALL SSCAL( KI, SCALE, WORK( 1+(IV-1)*N ), 1 )
                                CALL SSCAL( KI, SCALE, WORK( 1+(IV  )*N ), 1 )
                            END IF
                            WORK( J-1+(IV-1)*N ) = X( 1, 1 )
                            WORK( J  +(IV-1)*N ) = X( 2, 1 )
                            WORK( J-1+(IV  )*N ) = X( 1, 2 )
                            WORK( J  +(IV  )*N ) = X( 2, 2 )
                            !
                            !                    Update the right-hand side
                            !
                            CALL SAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                WORK( 1+(IV-1)*N   ), 1 )
                            CALL SAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                WORK( 1+(IV-1)*N   ), 1 )
                            CALL SAXPY( J-2, -X( 1, 2 ), T( 1, J-1 ), 1, &
                                WORK( 1+(IV  )*N ), 1 )
                            CALL SAXPY( J-2, -X( 2, 2 ), T( 1, J ), 1, &
                                WORK( 1+(IV  )*N ), 1 )
                        END IF
90                  END DO
                    !
                    !              Copy the vector x or Q*x to VR and normalize.
                    !
                    IF( .NOT.OVER ) THEN
                        !                 ------------------------------
                        !                 no back-transform: copy x to VR and normalize.
                        CALL SCOPY( KI, WORK( 1+(IV-1)*N ), 1, VR(1,IS-1), 1 )
                        CALL SCOPY( KI, WORK( 1+(IV  )*N ), 1, VR(1,IS  ), 1 )
                        !
                        EMAX = ZERO
                        DO  K = 1, KI
                            EMAX = MAX( EMAX, ABS( VR( K, IS-1 ) )+ &
                                ABS( VR( K, IS   ) ) )
                        END DO
                        REMAX = ONE / EMAX
                        CALL SSCAL( KI, REMAX, VR( 1, IS-1 ), 1 )
                        CALL SSCAL( KI, REMAX, VR( 1, IS   ), 1 )
                        !
                        DO  K = KI + 1, N
                            VR( K, IS-1 ) = ZERO
                            VR( K, IS   ) = ZERO
                        END DO
                        !
                    ELSE IF( NB.EQ.1 ) THEN
                        !                 ------------------------------
                        !                 version 1: back-transform each vector with GEMV, Q*x.
                        IF( KI.GT.2 ) THEN
                            CALL SGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                WORK( 1    + (IV-1)*N ), 1, &
                                WORK( KI-1 + (IV-1)*N ), VR(1,KI-1), 1)
                            CALL SGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                WORK( 1  + (IV)*N ), 1, &
                                WORK( KI + (IV)*N ), VR( 1, KI ), 1 )
                        ELSE
                            CALL SSCAL( N, WORK(KI-1+(IV-1)*N), VR(1,KI-1), 1)
                            CALL SSCAL( N, WORK(KI  +(IV  )*N), VR(1,KI  ), 1)
                        END IF
                        !
                        EMAX = ZERO
                        DO  K = 1, N
                            EMAX = MAX( EMAX, ABS( VR( K, KI-1 ) )+ &
                                ABS( VR( K, KI   ) ) )
                        END DO
                        REMAX = ONE / EMAX
                        CALL SSCAL( N, REMAX, VR( 1, KI-1 ), 1 )
                        CALL SSCAL( N, REMAX, VR( 1, KI   ), 1 )
                        !
                    ELSE
                        !                 ------------------------------
                        !                 version 2: back-transform block of vectors with GEMM
                        !                 zero out below vector
                        DO K = KI + 1, N
                            WORK( K + (IV-1)*N ) = ZERO
                            WORK( K + (IV  )*N ) = ZERO
                        END DO
                        ISCOMPLEX( IV-1 ) = -IP
                        ISCOMPLEX( IV   ) =  IP
                        IV = IV - 1
                        !                 back-transform and normalization is done below
                    END IF
                END IF

                IF( NB.GT.1 ) THEN
                    !              --------------------------------------------------------
                    !              Blocked version of back-transform
                    !              For complex case, KI2 includes both vectors (KI-1 and KI)
                    IF( IP.EQ.0 ) THEN
                        KI2 = KI
                    ELSE
                        KI2 = KI - 1
                    END IF

                    !              Columns IV:NB of work are valid vectors.
                    !              When the number of vectors stored reaches NB-1 or NB,
                    !              or if this was last vector, do the GEMM
                    IF( (IV.LE.2) .OR. (KI2.EQ.1) ) THEN
                        CALL SGEMM( 'N', 'N', N, NB-IV+1, KI2+NB-IV, ONE, &
                            VR, LDVR, &
                            WORK( 1 + (IV)*N    ), N, &
                            ZERO, &
                            WORK( 1 + (NB+IV)*N ), N )
                        !                 normalize vectors
                        DO K = IV, NB
                            IF( ISCOMPLEX(K).EQ.0 ) THEN
                                !                       real eigenvector
                                II = ISAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                                REMAX = ONE / ABS( WORK( II + (NB+K)*N ) )
                            ELSE IF( ISCOMPLEX(K).EQ.1 ) THEN
                                !                       first eigenvector of conjugate pair
                                EMAX = ZERO
                                DO II = 1, N
                                    EMAX = MAX( EMAX, &
                                        ABS( WORK( II + (NB+K  )*N ) )+ &
                                        ABS( WORK( II + (NB+K+1)*N ) ) )
                                END DO
                                REMAX = ONE / EMAX
                                !                    else if ISCOMPLEX(K).EQ.-1
                                !                       second eigenvector of conjugate pair
                                !                       reuse same REMAX as previous K
                            END IF
                            CALL SSCAL( N, REMAX, WORK( 1 + (NB+K)*N ), 1 )
                        END DO
                        CALL SLACPY( 'F', N, NB-IV+1, &
                            WORK( 1 + (NB+IV)*N ), N, &
                            VR( 1, KI2 ), LDVR )
                        IV = NB
                    ELSE
                        IV = IV - 1
                    END IF
                END IF ! blocked back-transform
                !
                IS = IS - 1
                IF( IP.NE.0 ) &
                    IS = IS - 1
140         END DO
        END IF

        IF( LEFTV ) THEN
            !
            !        ============================================================
            !        Compute left eigenvectors.
            !
            !        IV is index of column in current block.
            !        For complex left vector, uses IV for real part and IV+1 for complex part.
            !        Non-blocked version always uses IV=1;
            !        blocked     version starts with IV=1, goes up to NB-1 or NB.
            !        (Note the "0-th" column is used for 1-norms computed above.)
            IV = 1
            IP = 0
            IS = 1
            DO  KI = 1, N
                IF( IP.EQ.1 ) THEN
                    !              previous iteration (ki-1) was first of conjugate pair,
                    !              so this ki is second of conjugate pair; skip to end of loop
                    IP = -1
                    GO TO 260
                ELSE IF( KI.EQ.N ) THEN
                    !              last column, so this ki must be real eigenvalue
                    IP = 0
                ELSE IF( T( KI+1, KI ).EQ.ZERO ) THEN
                    !              zero on sub-diagonal, so this ki is real eigenvalue
                    IP = 0
                ELSE
                    !              non-zero on sub-diagonal, so this ki is first of conjugate pair
                    IP = 1
                END IF
                !
                IF( SOMEV ) THEN
                    IF( .NOT.SELECT( KI ) ) &
                        GO TO 260
                END IF
                !
                !           Compute the KI-th eigenvalue (WR,WI).
                !
                WR = T( KI, KI )
                WI = ZERO
                IF( IP.NE.0 ) &
                    WI = SQRT( ABS( T( KI, KI+1 ) ) )* &
                    SQRT( ABS( T( KI+1, KI ) ) )
                SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
                !
                IF( IP.EQ.0 ) THEN
                    !
                    !              --------------------------------------------------------
                    !              Real left eigenvector
                    !
                    WORK( KI + IV*N ) = ONE
                    !
                    !              Form right-hand side.
                    !
                    DO  K = KI + 1, N
                        WORK( K + IV*N ) = -T( KI, K )
                    END DO
                    !
                    !              Solve transposed quasi-triangular system:
                    !              [ T(KI+1:N,KI+1:N) - WR ]**T * X = SCALE*WORK
                    !
                    VMAX = ONE
                    VCRIT = BIGNUM
                    !
                    JNXT = KI + 1
                    DO  J = KI + 1, N
                        IF( J.LT.JNXT ) &
                            GO TO 170
                        J1 = J
                        J2 = J
                        JNXT = J + 1
                        IF( J.LT.N ) THEN
                            IF( T( J+1, J ).NE.ZERO ) THEN
                                J2 = J + 1
                                JNXT = J + 2
                            END IF
                        END IF
                        !
                        IF( J1.EQ.J2 ) THEN
                            !
                            !                    1-by-1 diagonal block
                            !
                            !                    Scale if necessary to avoid overflow when forming
                            !                    the right-hand side.
                            !
                            IF( WORK( J ).GT.VCRIT ) THEN
                                REC = ONE / VMAX
                                CALL SSCAL( N-KI+1, REC, WORK( KI+IV*N ), 1 )
                                VMAX = ONE
                                VCRIT = BIGNUM
                            END IF
                            !
                            WORK( J+IV*N ) = WORK( J+IV*N ) - &
                                SDOT( J-KI-1, T( KI+1, J ), 1, &
                                WORK( KI+1+IV*N ), 1 )
                            !
                            !                    Solve [ T(J,J) - WR ]**T * X = WORK
                            !
                            CALL SLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                LDT, ONE, ONE, WORK( J+IV*N ), N, WR, &
                                ZERO, X, 2, SCALE, XNORM, IERR )
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) &
                                CALL SSCAL( N-KI+1, SCALE, WORK( KI+IV*N ), 1 )
                            WORK( J+IV*N ) = X( 1, 1 )
                            VMAX = MAX( ABS( WORK( J+IV*N ) ), VMAX )
                            VCRIT = BIGNUM / VMAX
                            !
                        ELSE
                            !
                            !                    2-by-2 diagonal block
                            !
                            !                    Scale if necessary to avoid overflow when forming
                            !                    the right-hand side.
                            !
                            BETA = MAX( WORK( J ), WORK( J+1 ) )
                            IF( BETA.GT.VCRIT ) THEN
                                REC = ONE / VMAX
                                CALL SSCAL( N-KI+1, REC, WORK( KI+IV*N ), 1 )
                                VMAX = ONE
                                VCRIT = BIGNUM
                            END IF
                            !
                            WORK( J+IV*N ) = WORK( J+IV*N ) - &
                                SDOT( J-KI-1, T( KI+1, J ), 1, &
                                WORK( KI+1+IV*N ), 1 )
                            !
                            WORK( J+1+IV*N ) = WORK( J+1+IV*N ) - &
                                SDOT( J-KI-1, T( KI+1, J+1 ), 1, &
                                WORK( KI+1+IV*N ), 1 )
                            !
                            !                    Solve
                            !                    [ T(J,J)-WR   T(J,J+1)      ]**T * X = SCALE*( WORK1 )
                            !                    [ T(J+1,J)    T(J+1,J+1)-WR ]                ( WORK2 )
                            !
                            CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, T( J, J ), &
                                LDT, ONE, ONE, WORK( J+IV*N ), N, WR, &
                                ZERO, X, 2, SCALE, XNORM, IERR )
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) &
                                CALL SSCAL( N-KI+1, SCALE, WORK( KI+IV*N ), 1 )
                            WORK( J  +IV*N ) = X( 1, 1 )
                            WORK( J+1+IV*N ) = X( 2, 1 )
                            !
                            VMAX = MAX( ABS( WORK( J  +IV*N ) ), &
                                ABS( WORK( J+1+IV*N ) ), VMAX )
                            VCRIT = BIGNUM / VMAX
                            !
                        END IF
170                 END DO
                    !
                    !              Copy the vector x or Q*x to VL and normalize.
                    !
                    IF( .NOT.OVER ) THEN
                        !                 ------------------------------
                        !                 no back-transform: copy x to VL and normalize.
                        CALL SCOPY( N-KI+1, WORK( KI + IV*N ), 1, &
                            VL( KI, IS ), 1 )
                        !
                        II = ISAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
                        REMAX = ONE / ABS( VL( II, IS ) )
                        CALL SSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
                        !
                        DO  K = 1, KI - 1
                            VL( K, IS ) = ZERO
                        END DO
                        !
                    ELSE IF( NB.EQ.1 ) THEN
                        !                 ------------------------------
                        !                 version 1: back-transform each vector with GEMV, Q*x.
                        IF( KI.LT.N ) &
                            CALL SGEMV( 'N', N, N-KI, ONE, &
                            VL( 1, KI+1 ), LDVL, &
                            WORK( KI+1 + IV*N ), 1, &
                            WORK( KI   + IV*N ), VL( 1, KI ), 1 )
                        !
                        II = ISAMAX( N, VL( 1, KI ), 1 )
                        REMAX = ONE / ABS( VL( II, KI ) )
                        CALL SSCAL( N, REMAX, VL( 1, KI ), 1 )
                        !
                    ELSE
                        !                 ------------------------------
                        !                 version 2: back-transform block of vectors with GEMM
                        !                 zero out above vector
                        !                 could go from KI-NV+1 to KI-1
                        DO K = 1, KI - 1
                            WORK( K + IV*N ) = ZERO
                        END DO
                        ISCOMPLEX( IV ) = IP
                        !                 back-transform and normalization is done below
                    END IF
                ELSE
                    !
                    !              --------------------------------------------------------
                    !              Complex left eigenvector.
                    !
                    !              Initial solve:
                    !              [ ( T(KI,KI)    T(KI,KI+1)  )**T - (WR - I* WI) ]*X = 0.
                    !              [ ( T(KI+1,KI) T(KI+1,KI+1) )                   ]
                    !
                    IF( ABS( T( KI, KI+1 ) ).GE.ABS( T( KI+1, KI ) ) ) THEN
                        WORK( KI   + (IV  )*N ) = WI / T( KI, KI+1 )
                        WORK( KI+1 + (IV+1)*N ) = ONE
                    ELSE
                        WORK( KI   + (IV  )*N ) = ONE
                        WORK( KI+1 + (IV+1)*N ) = -WI / T( KI+1, KI )
                    END IF
                    WORK( KI+1 + (IV  )*N ) = ZERO
                    WORK( KI   + (IV+1)*N ) = ZERO
                    !
                    !              Form right-hand side.
                    !
                    DO  K = KI + 2, N
                        WORK( K+(IV  )*N ) = -WORK( KI  +(IV  )*N )*T(KI,  K)
                        WORK( K+(IV+1)*N ) = -WORK( KI+1+(IV+1)*N )*T(KI+1,K)
                    END DO
                    !
                    !              Solve transposed quasi-triangular system:
                    !              [ T(KI+2:N,KI+2:N)**T - (WR-i*WI) ]*X = WORK1+i*WORK2
                    !
                    VMAX = ONE
                    VCRIT = BIGNUM
                    !
                    JNXT = KI + 2
                    DO  J = KI + 2, N
                        IF( J.LT.JNXT ) &
                            GO TO 200
                        J1 = J
                        J2 = J
                        JNXT = J + 1
                        IF( J.LT.N ) THEN
                            IF( T( J+1, J ).NE.ZERO ) THEN
                                J2 = J + 1
                                JNXT = J + 2
                            END IF
                        END IF
                        !
                        IF( J1.EQ.J2 ) THEN
                            !
                            !                    1-by-1 diagonal block
                            !
                            !                    Scale if necessary to avoid overflow when
                            !                    forming the right-hand side elements.
                            !
                            IF( WORK( J ).GT.VCRIT ) THEN
                                REC = ONE / VMAX
                                CALL SSCAL( N-KI+1, REC, WORK(KI+(IV  )*N), 1 )
                                CALL SSCAL( N-KI+1, REC, WORK(KI+(IV+1)*N), 1 )
                                VMAX = ONE
                                VCRIT = BIGNUM
                            END IF
                            !
                            WORK( J+(IV  )*N ) = WORK( J+(IV)*N ) - &
                                SDOT( J-KI-2, T( KI+2, J ), 1, &
                                WORK( KI+2+(IV)*N ), 1 )
                            WORK( J+(IV+1)*N ) = WORK( J+(IV+1)*N ) - &
                                SDOT( J-KI-2, T( KI+2, J ), 1, &
                                WORK( KI+2+(IV+1)*N ), 1 )
                            !
                            !                    Solve [ T(J,J)-(WR-i*WI) ]*(X11+i*X12)= WK+I*WK2
                            !
                            CALL SLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                LDT, ONE, ONE, WORK( J+IV*N ), N, WR, &
                                -WI, X, 2, SCALE, XNORM, IERR )
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) THEN
                                CALL SSCAL( N-KI+1, SCALE, WORK(KI+(IV  )*N), 1)
                                CALL SSCAL( N-KI+1, SCALE, WORK(KI+(IV+1)*N), 1)
                            END IF
                            WORK( J+(IV  )*N ) = X( 1, 1 )
                            WORK( J+(IV+1)*N ) = X( 1, 2 )
                            VMAX = MAX( ABS( WORK( J+(IV  )*N ) ), &
                                ABS( WORK( J+(IV+1)*N ) ), VMAX )
                            VCRIT = BIGNUM / VMAX
                            !
                        ELSE
                            !
                            !                    2-by-2 diagonal block
                            !
                            !                    Scale if necessary to avoid overflow when forming
                            !                    the right-hand side elements.
                            !
                            BETA = MAX( WORK( J ), WORK( J+1 ) )
                            IF( BETA.GT.VCRIT ) THEN
                                REC = ONE / VMAX
                                CALL SSCAL( N-KI+1, REC, WORK(KI+(IV  )*N), 1 )
                                CALL SSCAL( N-KI+1, REC, WORK(KI+(IV+1)*N), 1 )
                                VMAX = ONE
                                VCRIT = BIGNUM
                            END IF
                            !
                            WORK( J  +(IV  )*N ) = WORK( J+(IV)*N ) - &
                                SDOT( J-KI-2, T( KI+2, J ), 1, &
                                WORK( KI+2+(IV)*N ), 1 )
                            !
                            WORK( J  +(IV+1)*N ) = WORK( J+(IV+1)*N ) - &
                                SDOT( J-KI-2, T( KI+2, J ), 1, &
                                WORK( KI+2+(IV+1)*N ), 1 )
                            !
                            WORK( J+1+(IV  )*N ) = WORK( J+1+(IV)*N ) - &
                                SDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                WORK( KI+2+(IV)*N ), 1 )
                            !
                            WORK( J+1+(IV+1)*N ) = WORK( J+1+(IV+1)*N ) - &
                                SDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                WORK( KI+2+(IV+1)*N ), 1 )
                            !
                            !                    Solve 2-by-2 complex linear equation
                            !                    [ (T(j,j)   T(j,j+1)  )**T - (wr-i*wi)*I ]*X = SCALE*B
                            !                    [ (T(j+1,j) T(j+1,j+1))                  ]
                            !
                            CALL SLALN2( .TRUE., 2, 2, SMIN, ONE, T( J, J ), &
                                LDT, ONE, ONE, WORK( J+IV*N ), N, WR, &
                                -WI, X, 2, SCALE, XNORM, IERR )
                            !
                            !                    Scale if necessary
                            !
                            IF( SCALE.NE.ONE ) THEN
                                CALL SSCAL( N-KI+1, SCALE, WORK(KI+(IV  )*N), 1)
                                CALL SSCAL( N-KI+1, SCALE, WORK(KI+(IV+1)*N), 1)
                            END IF
                            WORK( J  +(IV  )*N ) = X( 1, 1 )
                            WORK( J  +(IV+1)*N ) = X( 1, 2 )
                            WORK( J+1+(IV  )*N ) = X( 2, 1 )
                            WORK( J+1+(IV+1)*N ) = X( 2, 2 )
                            VMAX = MAX( ABS( X( 1, 1 ) ), ABS( X( 1, 2 ) ), &
                                ABS( X( 2, 1 ) ), ABS( X( 2, 2 ) ), &
                                VMAX )
                            VCRIT = BIGNUM / VMAX
                            !
                        END IF
200                 END DO
                    !
                    !              Copy the vector x or Q*x to VL and normalize.
                    !
                    IF( .NOT.OVER ) THEN
                        !                 ------------------------------
                        !                 no back-transform: copy x to VL and normalize.
                        CALL SCOPY( N-KI+1, WORK( KI + (IV  )*N ), 1, &
                            VL( KI, IS   ), 1 )
                        CALL SCOPY( N-KI+1, WORK( KI + (IV+1)*N ), 1, &
                            VL( KI, IS+1 ), 1 )
                        !
                        EMAX = ZERO
                        DO  K = KI, N
                            EMAX = MAX( EMAX, ABS( VL( K, IS   ) )+ &
                                ABS( VL( K, IS+1 ) ) )
                        END DO
                        REMAX = ONE / EMAX
                        CALL SSCAL( N-KI+1, REMAX, VL( KI, IS   ), 1 )
                        CALL SSCAL( N-KI+1, REMAX, VL( KI, IS+1 ), 1 )
                        !
                        DO  K = 1, KI - 1
                            VL( K, IS   ) = ZERO
                            VL( K, IS+1 ) = ZERO
                        END DO
                        !
                    ELSE IF( NB.EQ.1 ) THEN
                        !                 ------------------------------
                        !                 version 1: back-transform each vector with GEMV, Q*x.
                        IF( KI.LT.N-1 ) THEN
                            CALL SGEMV( 'N', N, N-KI-1, ONE, &
                                VL( 1, KI+2 ), LDVL, &
                                WORK( KI+2 + (IV)*N ), 1, &
                                WORK( KI   + (IV)*N ), &
                                VL( 1, KI ), 1 )
                            CALL SGEMV( 'N', N, N-KI-1, ONE, &
                                VL( 1, KI+2 ), LDVL, &
                                WORK( KI+2 + (IV+1)*N ), 1, &
                                WORK( KI+1 + (IV+1)*N ), &
                                VL( 1, KI+1 ), 1 )
                        ELSE
                            CALL SSCAL( N, WORK(KI+  (IV  )*N), VL(1, KI  ), 1)
                            CALL SSCAL( N, WORK(KI+1+(IV+1)*N), VL(1, KI+1), 1)
                        END IF
                        !
                        EMAX = ZERO
                        DO  K = 1, N
                            EMAX = MAX( EMAX, ABS( VL( K, KI   ) )+ &
                                ABS( VL( K, KI+1 ) ) )
                        END DO
                        REMAX = ONE / EMAX
                        CALL SSCAL( N, REMAX, VL( 1, KI   ), 1 )
                        CALL SSCAL( N, REMAX, VL( 1, KI+1 ), 1 )
                        !
                    ELSE
                        !                 ------------------------------
                        !                 version 2: back-transform block of vectors with GEMM
                        !                 zero out above vector
                        !                 could go from KI-NV+1 to KI-1
                        DO K = 1, KI - 1
                            WORK( K + (IV  )*N ) = ZERO
                            WORK( K + (IV+1)*N ) = ZERO
                        END DO
                        ISCOMPLEX( IV   ) =  IP
                        ISCOMPLEX( IV+1 ) = -IP
                        IV = IV + 1
                        !                 back-transform and normalization is done below
                    END IF
                END IF

                IF( NB.GT.1 ) THEN
                    !              --------------------------------------------------------
                    !              Blocked version of back-transform
                    !              For complex case, KI2 includes both vectors (KI and KI+1)
                    IF( IP.EQ.0 ) THEN
                        KI2 = KI
                    ELSE
                        KI2 = KI + 1
                    END IF

                    !              Columns 1:IV of work are valid vectors.
                    !              When the number of vectors stored reaches NB-1 or NB,
                    !              or if this was last vector, do the GEMM
                    IF( (IV.GE.NB-1) .OR. (KI2.EQ.N) ) THEN
                        CALL SGEMM( 'N', 'N', N, IV, N-KI2+IV, ONE, &
                            VL( 1, KI2-IV+1 ), LDVL, &
                            WORK( KI2-IV+1 + (1)*N ), N, &
                            ZERO, &
                            WORK( 1 + (NB+1)*N ), N )
                        !                 normalize vectors
                        DO K = 1, IV
                            IF( ISCOMPLEX(K).EQ.0) THEN
                                !                       real eigenvector
                                II = ISAMAX( N, WORK( 1 + (NB+K)*N ), 1 )
                                REMAX = ONE / ABS( WORK( II + (NB+K)*N ) )
                            ELSE IF( ISCOMPLEX(K).EQ.1) THEN
                                !                       first eigenvector of conjugate pair
                                EMAX = ZERO
                                DO II = 1, N
                                    EMAX = MAX( EMAX, &
                                        ABS( WORK( II + (NB+K  )*N ) )+ &
                                        ABS( WORK( II + (NB+K+1)*N ) ) )
                                END DO
                                REMAX = ONE / EMAX
                                !                    else if ISCOMPLEX(K).EQ.-1
                                !                       second eigenvector of conjugate pair
                                !                       reuse same REMAX as previous K
                            END IF
                            CALL SSCAL( N, REMAX, WORK( 1 + (NB+K)*N ), 1 )
                        END DO
                        CALL SLACPY( 'F', N, IV, &
                            WORK( 1 + (NB+1)*N ), N, &
                            VL( 1, KI2-IV+1 ), LDVL )
                        IV = 1
                    ELSE
                        IV = IV + 1
                    END IF
                END IF ! blocked back-transform
                !
                IS = IS + 1
                IF( IP.NE.0 ) &
                    IS = IS + 1
260         END DO
        END IF
        !
        RETURN
        !
        !     End of STREVC3
        !
    END SUBROUTINE STREVC3

    SUBROUTINE STREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, &
        INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          COMPQ
        INTEGER            IFST, ILST, INFO, LDQ, LDT, N
        !     ..
        !     .. Array Arguments ..
        REAL               Q( LDQ, * ), T( LDT, * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO
        PARAMETER          ( ZERO = 0.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            WANTQ
        INTEGER            HERE, NBF, NBL, NBNEXT
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Decode and test the input arguments.
        !
        INFO = 0
        WANTQ = LSAME( COMPQ, 'V' )
        IF( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) THEN
            INFO = -1
        ELSE IF( N.LT.0 ) THEN
            INFO = -2
        ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
            INFO = -4
        ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) ) THEN
            INFO = -6
        ELSE IF(( IFST.LT.1 .OR. IFST.GT.N ).AND.( N.GT.0 )) THEN
            INFO = -7
        ELSE IF(( ILST.LT.1 .OR. ILST.GT.N ).AND.( N.GT.0 )) THEN
            INFO = -8
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'STREXC', -INFO )
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.LE.1 ) &
            RETURN
        !
        !     Determine the first row of specified block
        !     and find out it is 1 by 1 or 2 by 2.
        !
        IF( IFST.GT.1 ) THEN
            IF( T( IFST, IFST-1 ).NE.ZERO ) &
                IFST = IFST - 1
        END IF
        NBF = 1
        IF( IFST.LT.N ) THEN
            IF( T( IFST+1, IFST ).NE.ZERO ) &
                NBF = 2
        END IF
        !
        !     Determine the first row of the final block
        !     and find out it is 1 by 1 or 2 by 2.
        !
        IF( ILST.GT.1 ) THEN
            IF( T( ILST, ILST-1 ).NE.ZERO ) &
                ILST = ILST - 1
        END IF
        NBL = 1
        IF( ILST.LT.N ) THEN
            IF( T( ILST+1, ILST ).NE.ZERO ) &
                NBL = 2
        END IF
        !
        IF( IFST.EQ.ILST ) &
            RETURN
        !
        IF( IFST.LT.ILST ) THEN
            !
            !        Update ILST
            !
            IF( NBF.EQ.2 .AND. NBL.EQ.1 ) &
                ILST = ILST - 1
            IF( NBF.EQ.1 .AND. NBL.EQ.2 ) &
                ILST = ILST + 1
            !
            HERE = IFST
            !
10          CONTINUE
            !
            !        Swap block with next one below
            !
            IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
                !
                !           Current block either 1 by 1 or 2 by 2
                !
                NBNEXT = 1
                IF( HERE+NBF+1.LE.N ) THEN
                    IF( T( HERE+NBF+1, HERE+NBF ).NE.ZERO ) &
                        NBNEXT = 2
                END IF
                CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT, &
                    WORK, INFO )
                IF( INFO.NE.0 ) THEN
                    ILST = HERE
                    RETURN
                END IF
                HERE = HERE + NBNEXT
                !
                !           Test if 2 by 2 block breaks into two 1 by 1 blocks
                !
                IF( NBF.EQ.2 ) THEN
                    IF( T( HERE+1, HERE ).EQ.ZERO ) &
                        NBF = 3
                END IF
                !
            ELSE
                !
                !           Current block consists of two 1 by 1 blocks each of which
                !           must be swapped individually
                !
                NBNEXT = 1
                IF( HERE+3.LE.N ) THEN
                    IF( T( HERE+3, HERE+2 ).NE.ZERO ) &
                        NBNEXT = 2
                END IF
                CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT, &
                    WORK, INFO )
                IF( INFO.NE.0 ) THEN
                    ILST = HERE
                    RETURN
                END IF
                IF( NBNEXT.EQ.1 ) THEN
                    !
                    !              Swap two 1 by 1 blocks, no problems possible
                    !
                    CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, &
                        WORK, INFO )
                    HERE = HERE + 1
                ELSE
                    !
                    !              Recompute NBNEXT in case 2 by 2 split
                    !
                    IF( T( HERE+2, HERE+1 ).EQ.ZERO ) &
                        NBNEXT = 1
                    IF( NBNEXT.EQ.2 ) THEN
                        !
                        !                 2 by 2 Block did not split
                        !
                        CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, &
                            NBNEXT, WORK, INFO )
                        IF( INFO.NE.0 ) THEN
                            ILST = HERE
                            RETURN
                        END IF
                        HERE = HERE + 2
                    ELSE
                        !
                        !                 2 by 2 Block did split
                        !
                        CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                            WORK, INFO )
                        CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1, &
                            WORK, INFO )
                        HERE = HERE + 2
                    END IF
                END IF
            END IF
            IF( HERE.LT.ILST ) &
                GO TO 10
            !
        ELSE
            !
            HERE = IFST
20          CONTINUE
            !
            !        Swap block with next one above
            !
            IF( NBF.EQ.1 .OR. NBF.EQ.2 ) THEN
                !
                !           Current block either 1 by 1 or 2 by 2
                !
                NBNEXT = 1
                IF( HERE.GE.3 ) THEN
                    IF( T( HERE-1, HERE-2 ).NE.ZERO ) &
                        NBNEXT = 2
                END IF
                CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                    NBF, WORK, INFO )
                IF( INFO.NE.0 ) THEN
                    ILST = HERE
                    RETURN
                END IF
                HERE = HERE - NBNEXT
                !
                !           Test if 2 by 2 block breaks into two 1 by 1 blocks
                !
                IF( NBF.EQ.2 ) THEN
                    IF( T( HERE+1, HERE ).EQ.ZERO ) &
                        NBF = 3
                END IF
                !
            ELSE
                !
                !           Current block consists of two 1 by 1 blocks each of which
                !           must be swapped individually
                !
                NBNEXT = 1
                IF( HERE.GE.3 ) THEN
                    IF( T( HERE-1, HERE-2 ).NE.ZERO ) &
                        NBNEXT = 2
                END IF
                CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                    1, WORK, INFO )
                IF( INFO.NE.0 ) THEN
                    ILST = HERE
                    RETURN
                END IF
                IF( NBNEXT.EQ.1 ) THEN
                    !
                    !              Swap two 1 by 1 blocks, no problems possible
                    !
                    CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1, &
                        WORK, INFO )
                    HERE = HERE - 1
                ELSE
                    !
                    !              Recompute NBNEXT in case 2 by 2 split
                    !
                    IF( T( HERE, HERE-1 ).EQ.ZERO ) &
                        NBNEXT = 1
                    IF( NBNEXT.EQ.2 ) THEN
                        !
                        !                 2 by 2 Block did not split
                        !
                        CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1, &
                            WORK, INFO )
                        IF( INFO.NE.0 ) THEN
                            ILST = HERE
                            RETURN
                        END IF
                        HERE = HERE - 2
                    ELSE
                        !
                        !                 2 by 2 Block did split
                        !
                        CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                            WORK, INFO )
                        CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1, &
                            WORK, INFO )
                        HERE = HERE - 2
                    END IF
                END IF
            END IF
            IF( HERE.GT.ILST ) &
                GO TO 20
        END IF
        ILST = HERE
        !
        RETURN
        !
        !     End of STREXC
        !
    END SUBROUTINE STREXC

    SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
        !
        !  -- Reference BLAS level3 routine (version 3.7.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        REAL ALPHA
        INTEGER LDA,LDB,M,N
        CHARACTER DIAG,SIDE,TRANSA,UPLO
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),B(LDB,*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !     .. Local Scalars ..
        REAL TEMP
        INTEGER I,INFO,J,K,NROWA
        LOGICAL LSIDE,NOUNIT,UPPER
        !     ..
        !     .. Parameters ..
        REAL ONE,ZERO
        PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
        !     ..
        !
        !     Test the input parameters.
        !
        LSIDE = LSAME(SIDE,'L')
        IF (LSIDE) THEN
            NROWA = M
        ELSE
            NROWA = N
        END IF
        NOUNIT = LSAME(DIAG,'N')
        UPPER = LSAME(UPLO,'U')
        !
        INFO = 0
        IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
            INFO = 1
        ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
            INFO = 2
        ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND. &
            (.NOT.LSAME(TRANSA,'T')) .AND. &
            (.NOT.LSAME(TRANSA,'C'))) THEN
            INFO = 3
        ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
            INFO = 4
        ELSE IF (M.LT.0) THEN
            INFO = 5
        ELSE IF (N.LT.0) THEN
            INFO = 6
        ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
            INFO = 9
        ELSE IF (LDB.LT.MAX(1,M)) THEN
            INFO = 11
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('STRMM ',INFO)
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF (M.EQ.0 .OR. N.EQ.0) RETURN
        !
        !     And when  alpha.eq.zero.
        !
        IF (ALPHA.EQ.ZERO) THEN
            DO  J = 1,N
                DO  I = 1,M
                    B(I,J) = ZERO
                END DO
            END DO
            RETURN
        END IF
        !
        !     Start the operations.
        !
        IF (LSIDE) THEN
            IF (LSAME(TRANSA,'N')) THEN
                !
                !           Form  B := alpha*A*B.
                !
                IF (UPPER) THEN
                    DO  J = 1,N
                        DO  K = 1,M
                            IF (B(K,J).NE.ZERO) THEN
                                TEMP = ALPHA*B(K,J)
                                DO  I = 1,K - 1
                                    B(I,J) = B(I,J) + TEMP*A(I,K)
                                END DO
                                IF (NOUNIT) TEMP = TEMP*A(K,K)
                                B(K,J) = TEMP
                            END IF
                        END DO
                    END DO
                ELSE
                    DO  J = 1,N
                        DO  K = M,1,-1
                            IF (B(K,J).NE.ZERO) THEN
                                TEMP = ALPHA*B(K,J)
                                B(K,J) = TEMP
                                IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                                DO  I = K + 1,M
                                    B(I,J) = B(I,J) + TEMP*A(I,K)
                                END DO
                            END IF
                        END DO
                    END DO
                END IF
            ELSE
                !
                !           Form  B := alpha*A**T*B.
                !
                IF (UPPER) THEN
                    DO  J = 1,N
                        DO  I = M,1,-1
                            TEMP = B(I,J)
                            IF (NOUNIT) TEMP = TEMP*A(I,I)
                            DO  K = 1,I - 1
                                TEMP = TEMP + A(K,I)*B(K,J)
                            END DO
                            B(I,J) = ALPHA*TEMP
                        END DO
                    END DO
                ELSE
                    DO  J = 1,N
                        DO  I = 1,M
                            TEMP = B(I,J)
                            IF (NOUNIT) TEMP = TEMP*A(I,I)
                            DO  K = I + 1,M
                                TEMP = TEMP + A(K,I)*B(K,J)
                            END DO
                            B(I,J) = ALPHA*TEMP
                        END DO
                    END DO
                END IF
            END IF
        ELSE
            IF (LSAME(TRANSA,'N')) THEN
                !
                !           Form  B := alpha*B*A.
                !
                IF (UPPER) THEN
                    DO  J = N,1,-1
                        TEMP = ALPHA
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO  I = 1,M
                            B(I,J) = TEMP*B(I,J)
                        END DO
                        DO  K = 1,J - 1
                            IF (A(K,J).NE.ZERO) THEN
                                TEMP = ALPHA*A(K,J)
                                DO  I = 1,M
                                    B(I,J) = B(I,J) + TEMP*B(I,K)
                                END DO
                            END IF
                        END DO
                    END DO
                ELSE
                    DO  J = 1,N
                        TEMP = ALPHA
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO  I = 1,M
                            B(I,J) = TEMP*B(I,J)
                        END DO
                        DO  K = J + 1,N
                            IF (A(K,J).NE.ZERO) THEN
                                TEMP = ALPHA*A(K,J)
                                DO  I = 1,M
                                    B(I,J) = B(I,J) + TEMP*B(I,K)
                                END DO
                            END IF
                        END DO
                    END DO
                END IF
            ELSE
                !
                !           Form  B := alpha*B*A**T.
                !
                IF (UPPER) THEN
                    DO  K = 1,N
                        DO  J = 1,K - 1
                            IF (A(J,K).NE.ZERO) THEN
                                TEMP = ALPHA*A(J,K)
                                DO  I = 1,M
                                    B(I,J) = B(I,J) + TEMP*B(I,K)
                                END DO
                            END IF
                        END DO
                        TEMP = ALPHA
                        IF (NOUNIT) TEMP = TEMP*A(K,K)
                        IF (TEMP.NE.ONE) THEN
                            DO  I = 1,M
                                B(I,K) = TEMP*B(I,K)
                            END DO
                        END IF
                    END DO
                ELSE
                    DO  K = N,1,-1
                        DO  J = K + 1,N
                            IF (A(J,K).NE.ZERO) THEN
                                TEMP = ALPHA*A(J,K)
                                DO  I = 1,M
                                    B(I,J) = B(I,J) + TEMP*B(I,K)
                                END DO
                            END IF
                        END DO
                        TEMP = ALPHA
                        IF (NOUNIT) TEMP = TEMP*A(K,K)
                        IF (TEMP.NE.ONE) THEN
                            DO  I = 1,M
                                B(I,K) = TEMP*B(I,K)
                            END DO
                        END IF
                    END DO
                END IF
            END IF
        END IF
        !
        RETURN
        !
        !     End of STRMM .
        !
    END SUBROUTINE STRMM

    SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
        !
        !  -- Reference BLAS level2 routine (version 3.7.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER INCX,LDA,N
        CHARACTER DIAG,TRANS,UPLO
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),X(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL ZERO
        PARAMETER (ZERO=0.0E+0)
        !     ..
        !     .. Local Scalars ..
        REAL TEMP
        INTEGER I,INFO,IX,J,JX,KX
        LOGICAL NOUNIT
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
        INFO = 0
        IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
            INFO = 1
        ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. &
            .NOT.LSAME(TRANS,'C')) THEN
            INFO = 2
        ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
            INFO = 3
        ELSE IF (N.LT.0) THEN
            INFO = 4
        ELSE IF (LDA.LT.MAX(1,N)) THEN
            INFO = 6
        ELSE IF (INCX.EQ.0) THEN
            INFO = 8
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('STRMV ',INFO)
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF (N.EQ.0) RETURN
        !
        NOUNIT = LSAME(DIAG,'N')
        !
        !     Set up the start point in X if the increment is not unity. This
        !     will be  ( N - 1 )*INCX  too small for descending loops.
        !
        IF (INCX.LE.0) THEN
            KX = 1 - (N-1)*INCX
        ELSE IF (INCX.NE.1) THEN
            KX = 1
        END IF
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through A.
        !
        IF (LSAME(TRANS,'N')) THEN
            !
            !        Form  x := A*x.
            !
            IF (LSAME(UPLO,'U')) THEN
                IF (INCX.EQ.1) THEN
                    DO  J = 1,N
                        IF (X(J).NE.ZERO) THEN
                            TEMP = X(J)
                            DO  I = 1,J - 1
                                X(I) = X(I) + TEMP*A(I,J)
                            END DO
                            IF (NOUNIT) X(J) = X(J)*A(J,J)
                        END IF
                    END DO
                ELSE
                    JX = KX
                    DO  J = 1,N
                        IF (X(JX).NE.ZERO) THEN
                            TEMP = X(JX)
                            IX = KX
                            DO  I = 1,J - 1
                                X(IX) = X(IX) + TEMP*A(I,J)
                                IX = IX + INCX
                            END DO
                            IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                        END IF
                        JX = JX + INCX
                    END DO
                END IF
            ELSE
                IF (INCX.EQ.1) THEN
                    DO  J = N,1,-1
                        IF (X(J).NE.ZERO) THEN
                            TEMP = X(J)
                            DO  I = N,J + 1,-1
                                X(I) = X(I) + TEMP*A(I,J)
                            END DO
                            IF (NOUNIT) X(J) = X(J)*A(J,J)
                        END IF
                    END DO
                ELSE
                    KX = KX + (N-1)*INCX
                    JX = KX
                    DO  J = N,1,-1
                        IF (X(JX).NE.ZERO) THEN
                            TEMP = X(JX)
                            IX = KX
                            DO  I = N,J + 1,-1
                                X(IX) = X(IX) + TEMP*A(I,J)
                                IX = IX - INCX
                            END DO
                            IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                        END IF
                        JX = JX - INCX
                    END DO
                END IF
            END IF
        ELSE
            !
            !        Form  x := A**T*x.
            !
            IF (LSAME(UPLO,'U')) THEN
                IF (INCX.EQ.1) THEN
                    DO  J = N,1,-1
                        TEMP = X(J)
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO  I = J - 1,1,-1
                            TEMP = TEMP + A(I,J)*X(I)
                        END DO
                        X(J) = TEMP
                    END DO
                ELSE
                    JX = KX + (N-1)*INCX
                    DO  J = N,1,-1
                        TEMP = X(JX)
                        IX = JX
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO  I = J - 1,1,-1
                            IX = IX - INCX
                            TEMP = TEMP + A(I,J)*X(IX)
                        END DO
                        X(JX) = TEMP
                        JX = JX - INCX
                    END DO
                END IF
            ELSE
                IF (INCX.EQ.1) THEN
                    DO  J = 1,N
                        TEMP = X(J)
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO  I = J + 1,N
                            TEMP = TEMP + A(I,J)*X(I)
                        END DO
                        X(J) = TEMP
                    END DO
                ELSE
                    JX = KX
                    DO  J = 1,N
                        TEMP = X(JX)
                        IX = JX
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO  I = J + 1,N
                            IX = IX + INCX
                            TEMP = TEMP + A(I,J)*X(IX)
                        END DO
                        X(JX) = TEMP
                        JX = JX + INCX
                    END DO
                END IF
            END IF
        END IF
        !
        RETURN
        !
        !     End of STRMV .
        !
    END SUBROUTINE STRMV

    SUBROUTINE XERBLA( SRNAME, INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER*(*)      SRNAME
        INTEGER            INFO
        !     ..
        !
        ! =====================================================================
        !
        !     .. Intrinsic Functions ..
        INTRINSIC          LEN_TRIM
        !     ..
        !     .. Executable Statements ..
        !
        WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
        !
        STOP
        !
9999    FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ', &
            'an illegal value' )
        !
        !     End of XERBLA
        !
    END SUBROUTINE XERBLA

    LOGICAL FUNCTION lsame(CA,CB)
        !
        !  -- Reference BLAS level1 routine (version 3.1) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER CA,CB
        !     ..
        !
        ! =====================================================================
        !
        !     .. Intrinsic Functions ..
        INTRINSIC ichar
        !     ..
        !     .. Local Scalars ..
        INTEGER INTA,INTB,ZCODE
        !     ..
        !
        !     Test if the characters are equal
        !
        lsame = ca .EQ. cb
        IF (lsame) RETURN
        !
        !     Now test for equivalence if both characters are alphabetic.
        !
        zcode = ichar('Z')
        !
        !     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
        !     machines, on which ICHAR returns a value with bit 8 set.
        !     ICHAR('A') on Prime machines returns 193 which is the same as
        !     ICHAR('A') on an EBCDIC machine.
        !
        inta = ichar(ca)
        intb = ichar(cb)
        !
        IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN
            !
            !        ASCII is assumed - ZCODE is the ASCII code of either lower or
            !        upper case 'Z'.
            !
            IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
            IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
            !
        ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN
            !
            !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
            !        upper case 'Z'.
            !
            IF (inta.GE.129 .AND. inta.LE.137 .OR. &
                +        inta.GE.145 .AND. inta.LE.153 .OR. &
                +        inta.GE.162 .AND. inta.LE.169) inta = inta + 64
            IF (intb.GE.129 .AND. intb.LE.137 .OR. &
                +        intb.GE.145 .AND. intb.LE.153 .OR. &
                +        intb.GE.162 .AND. intb.LE.169) intb = intb + 64
            !
        ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN
            !
            !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
            !        plus 128 of either lower or upper case 'Z'.
            !
            IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
            IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
        END IF
        lsame = inta .EQ. intb
        !
        !     RETURN
        !
        !     End of LSAME
        !
    END FUNCTION lsame

    INTEGER FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
        !
        !  -- LAPACK auxiliary routine (version 3.8.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        CHARACTER*( * )    NAME, OPTS
        INTEGER            ISPEC, N1, N2, N3, N4
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        INTEGER            I, IC, IZ, NB, NBMIN, NX
        LOGICAL            CNAME, SNAME, TWOSTAGE
        CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*16
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          char, ichar, int, min, real
        !     ..
        !     ..
        !     .. Executable Statements ..
        !
        !!! GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
        !!!        130, 140, 150, 160, 160, 160, 160, 160)ispec
        !!! -> obsolecent feature; changed manually
        IF ((ISPEC == 1).OR.(ISPEC == 2).OR.(ISPEC == 3)) THEN
            GO TO 10
        ELSE IF (ISPEC == 4) THEN
            GO TO 80
        ELSE IF (ISPEC == 5) THEN
            GO TO 90
        ELSE IF (ISPEC == 6) THEN
            GO TO 100
        ELSE IF (ISPEC == 7) THEN
            GO TO 110
        ELSE IF (ISPEC == 8) THEN
            GO TO 120
        ELSE IF (ISPEC == 9) THEN
            GO TO 130
        ELSE IF (ISPEC == 10) THEN
            GO TO 140
        ELSE IF (ISPEC == 11) THEN
            GO TO 150
        ELSE
            GO TO 160
        END IF
        !
        !     Invalid value for ISPEC
        !
        ilaenv = -1
        RETURN
        !
10      CONTINUE
        !
        !     Convert NAME to upper case if the first character is lower case.
        !
        ilaenv = 1
        subnam = name
        ic = ichar( subnam( 1: 1 ) )
        iz = ichar( 'Z' )
        IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
            !
            !        ASCII character set
            !
            IF( ic.GE.97 .AND. ic.LE.122 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                    ic = ichar( subnam( i: i ) )
                    IF( ic.GE.97 .AND. ic.LE.122 ) &
                        subnam( i: i ) = char( ic-32 )
                END DO
            END IF
            !
        ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
            !
            !        EBCDIC character set
            !
            IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
                ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
                ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
                subnam( 1: 1 ) = char( ic+64 )
                DO i = 2, 6
                    ic = ichar( subnam( i: i ) )
                    IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
                        ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
                        ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i: &
                        i ) = char( ic+64 )
                END DO
            END IF
            !
        ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
            !
            !        Prime machines:  ASCII+128
            !
            IF( ic.GE.225 .AND. ic.LE.250 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                    ic = ichar( subnam( i: i ) )
                    IF( ic.GE.225 .AND. ic.LE.250 ) &
                        subnam( i: i ) = char( ic-32 )
                END DO
            END IF
        END IF
        !
        c1 = subnam( 1: 1 )
        sname = c1.EQ.'S' .OR. c1.EQ.'D'
        cname = c1.EQ.'C' .OR. c1.EQ.'Z'
        IF( .NOT.( cname .OR. sname ) ) &
            RETURN
        c2 = subnam( 2: 3 )
        c3 = subnam( 4: 6 )
        c4 = c3( 2: 3 )
        twostage = len( subnam ).GE.11 &
            .AND. subnam( 11: 11 ).EQ.'2'
        !
        !!! GO TO ( 50, 60, 70 )ispec -> obsolecent feature; changed manually
        IF (ISPEC == 1) THEN
            GO TO 50
        ELSE IF (ISPEC == 2) THEN
            GO TO 60
        ELSE
            GO TO 70
        END IF
        !
50      CONTINUE
        !
        !     ISPEC = 1:  block size
        !
        !     In these examples, separate code is provided for setting NB for
        !     real and complex.  We assume that NB will take the same value in
        !     single or double precision.
        !
        nb = 1
        !
        IF( c2.EQ.'GE' ) THEN
            IF( c3.EQ.'TRF' ) THEN
                IF( sname ) THEN
                    nb = 64
                ELSE
                    nb = 64
                END IF
            ELSE IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. &
                c3.EQ.'QLF' ) THEN
                IF( sname ) THEN
                    nb = 32
                ELSE
                    nb = 32
                END IF
            ELSE IF( c3.EQ.'QR ') THEN
                IF( n3 .EQ. 1) THEN
                    IF( sname ) THEN
                        !     M*N
                        IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                            nb = n1
                        ELSE
                            nb = 32768/n2
                        END IF
                    ELSE
                        IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                            nb = n1
                        ELSE
                            nb = 32768/n2
                        END IF
                    END IF
                ELSE
                    IF( sname ) THEN
                        nb = 1
                    ELSE
                        nb = 1
                    END IF
                END IF
            ELSE IF( c3.EQ.'LQ ') THEN
                IF( n3 .EQ. 2) THEN
                    IF( sname ) THEN
                        !     M*N
                        IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                            nb = n1
                        ELSE
                            nb = 32768/n2
                        END IF
                    ELSE
                        IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                            nb = n1
                        ELSE
                            nb = 32768/n2
                        END IF
                    END IF
                ELSE
                    IF( sname ) THEN
                        nb = 1
                    ELSE
                        nb = 1
                    END IF
                END IF
            ELSE IF( c3.EQ.'HRD' ) THEN
                IF( sname ) THEN
                    nb = 32
                ELSE
                    nb = 32
                END IF
            ELSE IF( c3.EQ.'BRD' ) THEN
                IF( sname ) THEN
                    nb = 32
                ELSE
                    nb = 32
                END IF
            ELSE IF( c3.EQ.'TRI' ) THEN
                IF( sname ) THEN
                    nb = 64
                ELSE
                    nb = 64
                END IF
            END IF
        ELSE IF( c2.EQ.'PO' ) THEN
            IF( c3.EQ.'TRF' ) THEN
                IF( sname ) THEN
                    nb = 64
                ELSE
                    nb = 64
                END IF
            END IF
        ELSE IF( c2.EQ.'SY' ) THEN
            IF( c3.EQ.'TRF' ) THEN
                IF( sname ) THEN
                    IF( twostage ) THEN
                        nb = 192
                    ELSE
                        nb = 64
                    END IF
                ELSE
                    IF( twostage ) THEN
                        nb = 192
                    ELSE
                        nb = 64
                    END IF
                END IF
            ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
                nb = 32
            ELSE IF( sname .AND. c3.EQ.'GST' ) THEN
                nb = 64
            END IF
        ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
            IF( c3.EQ.'TRF' ) THEN
                IF( twostage ) THEN
                    nb = 192
                ELSE
                    nb = 64
                END IF
            ELSE IF( c3.EQ.'TRD' ) THEN
                nb = 32
            ELSE IF( c3.EQ.'GST' ) THEN
                nb = 64
            END IF
        ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
            IF( c3( 1: 1 ).EQ.'G' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nb = 32
                END IF
            ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nb = 32
                END IF
            END IF
        ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
            IF( c3( 1: 1 ).EQ.'G' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nb = 32
                END IF
            ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nb = 32
                END IF
            END IF
        ELSE IF( c2.EQ.'GB' ) THEN
            IF( c3.EQ.'TRF' ) THEN
                IF( sname ) THEN
                    IF( n4.LE.64 ) THEN
                        nb = 1
                    ELSE
                        nb = 32
                    END IF
                ELSE
                    IF( n4.LE.64 ) THEN
                        nb = 1
                    ELSE
                        nb = 32
                    END IF
                END IF
            END IF
        ELSE IF( c2.EQ.'PB' ) THEN
            IF( c3.EQ.'TRF' ) THEN
                IF( sname ) THEN
                    IF( n2.LE.64 ) THEN
                        nb = 1
                    ELSE
                        nb = 32
                    END IF
                ELSE
                    IF( n2.LE.64 ) THEN
                        nb = 1
                    ELSE
                        nb = 32
                    END IF
                END IF
            END IF
        ELSE IF( c2.EQ.'TR' ) THEN
            IF( c3.EQ.'TRI' ) THEN
                IF( sname ) THEN
                    nb = 64
                ELSE
                    nb = 64
                END IF
            ELSE IF ( c3.EQ.'EVC' ) THEN
                IF( sname ) THEN
                    nb = 64
                ELSE
                    nb = 64
                END IF
            END IF
        ELSE IF( c2.EQ.'LA' ) THEN
            IF( c3.EQ.'UUM' ) THEN
                IF( sname ) THEN
                    nb = 64
                ELSE
                    nb = 64
                END IF
            END IF
        ELSE IF( sname .AND. c2.EQ.'ST' ) THEN
            IF( c3.EQ.'EBZ' ) THEN
                nb = 1
            END IF
        ELSE IF( c2.EQ.'GG' ) THEN
            nb = 32
            IF( c3.EQ.'HD3' ) THEN
                IF( sname ) THEN
                    nb = 32
                ELSE
                    nb = 32
                END IF
            END IF
        END IF
        ilaenv = nb
        RETURN
        !
60      CONTINUE
        !
        !     ISPEC = 2:  minimum block size
        !
        nbmin = 2
        IF( c2.EQ.'GE' ) THEN
            IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. &
                'QLF' ) THEN
                IF( sname ) THEN
                    nbmin = 2
                ELSE
                    nbmin = 2
                END IF
            ELSE IF( c3.EQ.'HRD' ) THEN
                IF( sname ) THEN
                    nbmin = 2
                ELSE
                    nbmin = 2
                END IF
            ELSE IF( c3.EQ.'BRD' ) THEN
                IF( sname ) THEN
                    nbmin = 2
                ELSE
                    nbmin = 2
                END IF
            ELSE IF( c3.EQ.'TRI' ) THEN
                IF( sname ) THEN
                    nbmin = 2
                ELSE
                    nbmin = 2
                END IF
            END IF
        ELSE IF( c2.EQ.'SY' ) THEN
            IF( c3.EQ.'TRF' ) THEN
                IF( sname ) THEN
                    nbmin = 8
                ELSE
                    nbmin = 8
                END IF
            ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
                nbmin = 2
            END IF
        ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
            IF( c3.EQ.'TRD' ) THEN
                nbmin = 2
            END IF
        ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
            IF( c3( 1: 1 ).EQ.'G' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nbmin = 2
                END IF
            ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nbmin = 2
                END IF
            END IF
        ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
            IF( c3( 1: 1 ).EQ.'G' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nbmin = 2
                END IF
            ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nbmin = 2
                END IF
            END IF
        ELSE IF( c2.EQ.'GG' ) THEN
            nbmin = 2
            IF( c3.EQ.'HD3' ) THEN
                nbmin = 2
            END IF
        END IF
        ilaenv = nbmin
        RETURN
        !
70      CONTINUE
        !
        !     ISPEC = 3:  crossover point
        !
        nx = 0
        IF( c2.EQ.'GE' ) THEN
            IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. &
                'QLF' ) THEN
                IF( sname ) THEN
                    nx = 128
                ELSE
                    nx = 128
                END IF
            ELSE IF( c3.EQ.'HRD' ) THEN
                IF( sname ) THEN
                    nx = 128
                ELSE
                    nx = 128
                END IF
            ELSE IF( c3.EQ.'BRD' ) THEN
                IF( sname ) THEN
                    nx = 128
                ELSE
                    nx = 128
                END IF
            END IF
        ELSE IF( c2.EQ.'SY' ) THEN
            IF( sname .AND. c3.EQ.'TRD' ) THEN
                nx = 32
            END IF
        ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
            IF( c3.EQ.'TRD' ) THEN
                nx = 32
            END IF
        ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
            IF( c3( 1: 1 ).EQ.'G' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nx = 128
                END IF
            END IF
        ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
            IF( c3( 1: 1 ).EQ.'G' ) THEN
                IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                    'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                    THEN
                    nx = 128
                END IF
            END IF
        ELSE IF( c2.EQ.'GG' ) THEN
            nx = 128
            IF( c3.EQ.'HD3' ) THEN
                nx = 128
            END IF
        END IF
        ilaenv = nx
        RETURN
        !
80      CONTINUE
        !
        !     ISPEC = 4:  number of shifts (used by xHSEQR)
        !
        ilaenv = 6
        RETURN
        !
90      CONTINUE
        !
        !     ISPEC = 5:  minimum column dimension (not used)
        !
        ilaenv = 2
        RETURN
        !
100     CONTINUE
        !
        !     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
        !
        ilaenv = int( REAL( MIN( N1, N2 ) )*1.6e0 )
        RETURN
        !
110     CONTINUE
        !
        !     ISPEC = 7:  number of processors (not used)
        !
        ilaenv = 1
        RETURN
        !
120     CONTINUE
        !
        !     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
        !
        ilaenv = 50
        RETURN
        !
130     CONTINUE
        !
        !     ISPEC = 9:  maximum size of the subproblems at the bottom of the
        !                 computation tree in the divide-and-conquer algorithm
        !                 (used by xGELSD and xGESDD)
        !
        ilaenv = 25
        RETURN
        !
140     CONTINUE
        !
        !     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
        !
        !     ILAENV = 0
        ilaenv = 1
        IF( ilaenv.EQ.1 ) THEN
            ilaenv = ieeeck( 1, 0.0, 1.0 )
        END IF
        RETURN
        !
150     CONTINUE
        !
        !     ISPEC = 11: infinity arithmetic can be trusted not to trap
        !
        !     ILAENV = 0
        ilaenv = 1
        IF( ilaenv.EQ.1 ) THEN
            ilaenv = ieeeck( 0, 0.0, 1.0 )
        END IF
        RETURN
        !
160     CONTINUE
        !
        !     12 <= ISPEC <= 16: xHSEQR or related subroutines.
        !
        ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
        RETURN
        !
        !     End of ILAENV
        !
    END FUNCTION ilaenv
    INTEGER          FUNCTION ieeeck( ISPEC, ZERO, ONE )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            ISPEC
        REAL               ONE, ZERO
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
            negzro, newzro, posinf
        !     ..
        !     .. Executable Statements ..
        ieeeck = 1
        !
        posinf = one / zero
        IF( posinf.LE.one ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        neginf = -one / zero
        IF( neginf.GE.zero ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        negzro = one / ( neginf+one )
        IF( negzro.NE.zero ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        neginf = one / negzro
        IF( neginf.GE.zero ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        newzro = negzro + zero
        IF( newzro.NE.zero ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        posinf = one / newzro
        IF( posinf.LE.one ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        neginf = neginf*posinf
        IF( neginf.GE.zero ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        posinf = posinf*posinf
        IF( posinf.LE.one ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        !
        !
        !
        !     Return if we were only asked to check infinity arithmetic
        !
        IF( ispec.EQ.0 ) &
            RETURN
        !
        nan1 = posinf + neginf
        !
        nan2 = posinf / neginf
        !
        nan3 = posinf / posinf
        !
        nan4 = posinf*zero
        !
        nan5 = neginf*negzro
        !
        nan6 = nan5*zero
        !
        IF( nan1.EQ.nan1 ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        IF( nan2.EQ.nan2 ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        IF( nan3.EQ.nan3 ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        IF( nan4.EQ.nan4 ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        IF( nan5.EQ.nan5 ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        IF( nan6.EQ.nan6 ) THEN
            ieeeck = 0
            RETURN
        END IF
        !
        RETURN
    END FUNCTION ieeeck

    INTEGER FUNCTION iparmq( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2017
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, ILO, ISPEC, LWORK, N
        CHARACTER          NAME*( * ), OPTS*( * )
        !
        !  ================================================================
        !     .. Parameters ..
        INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
        parameter( inmin = 12, inwin = 13, inibl = 14, &
            ishfts = 15, iacc22 = 16 )
        INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
        parameter( nmin = 75, k22min = 14, kacmin = 14, &
            nibble = 14, knwswp = 500 )
        REAL               TWO
        parameter( two = 2.0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            NH, NS
        INTEGER            I, IC, IZ
        CHARACTER          SUBNAM*6
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          log, max, mod, nint, real
        !     ..
        !     .. Executable Statements ..
        IF( ( ispec.EQ.ishfts ) .OR. ( ispec.EQ.inwin ) .OR. &
            ( ispec.EQ.iacc22 ) ) THEN
            !
            !        ==== Set the number simultaneous shifts ====
            !
            nh = ihi - ilo + 1
            ns = 2
            IF( nh.GE.30 ) &
                ns = 4
            IF( nh.GE.60 ) &
                ns = 10
            IF( nh.GE.150 ) &
                ns = max( 10, nh / nint( log( REAL( NH ) ) / log( TWO ) ) )
            IF( nh.GE.590 ) &
                ns = 64
            IF( nh.GE.3000 ) &
                ns = 128
            IF( nh.GE.6000 ) &
                ns = 256
            ns = max( 2, ns-mod( ns, 2 ) )
        END IF
        !
        IF( ispec.EQ.inmin ) THEN
            !
            !
            !        ===== Matrices of order smaller than NMIN get sent
            !        .     to xLAHQR, the classic double shift algorithm.
            !        .     This must be at least 11. ====
            !
            iparmq = nmin
            !
        ELSE IF( ispec.EQ.inibl ) THEN
            !
            !        ==== INIBL: skip a multi-shift qr iteration and
            !        .    whenever aggressive early deflation finds
            !        .    at least (NIBBLE*(window size)/100) deflations. ====
            !
            iparmq = nibble
            !
        ELSE IF( ispec.EQ.ishfts ) THEN
            !
            !        ==== NSHFTS: The number of simultaneous shifts =====
            !
            iparmq = ns
            !
        ELSE IF( ispec.EQ.inwin ) THEN
            !
            !        ==== NW: deflation window size.  ====
            !
            IF( nh.LE.knwswp ) THEN
                iparmq = ns
            ELSE
                iparmq = 3*ns / 2
            END IF
            !
        ELSE IF( ispec.EQ.iacc22 ) THEN
            !
            !        ==== IACC22: Whether to accumulate reflections
            !        .     before updating the far-from-diagonal elements
            !        .     and whether to use 2-by-2 block structure while
            !        .     doing it.  A small amount of work could be saved
            !        .     by making this choice dependent also upon the
            !        .     NH=IHI-ILO+1.
            !
            !
            !        Convert NAME to upper case if the first character is lower case.
            !
            iparmq = 0
            subnam = name
            ic = ichar( subnam( 1: 1 ) )
            iz = ichar( 'Z' )
            IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
                !
                !           ASCII character set
                !
                IF( ic.GE.97 .AND. ic.LE.122 ) THEN
                    subnam( 1: 1 ) = char( ic-32 )
                    DO i = 2, 6
                        ic = ichar( subnam( i: i ) )
                        IF( ic.GE.97 .AND. ic.LE.122 ) &
                            subnam( i: i ) = char( ic-32 )
                    END DO
                END IF
                !
            ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
                !
                !           EBCDIC character set
                !
                IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
                    ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
                    ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
                    subnam( 1: 1 ) = char( ic+64 )
                    DO i = 2, 6
                        ic = ichar( subnam( i: i ) )
                        IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
                            ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
                            ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i: &
                            i ) = char( ic+64 )
                    END DO
                END IF
                !
            ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
                !
                !           Prime machines:  ASCII+128
                !
                IF( ic.GE.225 .AND. ic.LE.250 ) THEN
                    subnam( 1: 1 ) = char( ic-32 )
                    DO i = 2, 6
                        ic = ichar( subnam( i: i ) )
                        IF( ic.GE.225 .AND. ic.LE.250 ) &
                            subnam( i: i ) = char( ic-32 )
                    END DO
                END IF
            END IF
            !
            IF( subnam( 2:6 ).EQ.'GGHRD' .OR. &
                subnam( 2:6 ).EQ.'GGHD3' ) THEN
                iparmq = 1
                IF( nh.GE.k22min ) &
                    iparmq = 2
            ELSE IF ( subnam( 4:6 ).EQ.'EXC' ) THEN
                IF( nh.GE.kacmin ) &
                    iparmq = 1
                IF( nh.GE.k22min ) &
                    iparmq = 2
            ELSE IF ( subnam( 2:6 ).EQ.'HSEQR' .OR. &
                subnam( 2:5 ).EQ.'LAQR' ) THEN
                IF( ns.GE.kacmin ) &
                    iparmq = 1
                IF( ns.GE.k22min ) &
                    iparmq = 2
            END IF
            !
        ELSE
            !        ===== invalid value of ispec =====
            iparmq = -1
            !
        END IF
        !
        !     ==== End of IPARMQ ====
        !
    END FUNCTION iparmq

    REAL             FUNCTION slamch( CMACH )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          CMACH
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE, ZERO
        parameter( one = 1.0e+0, zero = 0.0e+0 )
        !     ..
        !     .. Local Scalars ..
        REAL               RND, EPS, SFMIN, SMALL, RMACH

        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          digits, epsilon, huge, maxexponent, &
            minexponent, radix, tiny
        !     ..
        !     .. Executable Statements ..
        !
        !
        !     Assume rounding, not chopping. Always.
        !
        rnd = one
        !
        IF( one.EQ.rnd ) THEN
            eps = epsilon(zero) * 0.5
        ELSE
            eps = epsilon(zero)
        END IF
        !
        IF( lsame( cmach, 'E' ) ) THEN
            rmach = eps
        ELSE IF( lsame( cmach, 'S' ) ) THEN
            sfmin = tiny(zero)
            small = one / huge(zero)
            IF( small.GE.sfmin ) THEN
                !
                !           Use SMALL plus a bit, to avoid the possibility of rounding
                !           causing overflow when computing  1/sfmin.
                !
                sfmin = small*( one+eps )
            END IF
            rmach = sfmin
        ELSE IF( lsame( cmach, 'B' ) ) THEN
            rmach = radix(zero)
        ELSE IF( lsame( cmach, 'P' ) ) THEN
            rmach = eps * radix(zero)
        ELSE IF( lsame( cmach, 'N' ) ) THEN
            rmach = digits(zero)
        ELSE IF( lsame( cmach, 'R' ) ) THEN
            rmach = rnd
        ELSE IF( lsame( cmach, 'M' ) ) THEN
            rmach = minexponent(zero)
        ELSE IF( lsame( cmach, 'U' ) ) THEN
            rmach = tiny(zero)
        ELSE IF( lsame( cmach, 'L' ) ) THEN
            rmach = maxexponent(zero)
        ELSE IF( lsame( cmach, 'O' ) ) THEN
            rmach = huge(zero)
        ELSE
            rmach = zero
        END IF
        !
        slamch = rmach
        RETURN
        !
        !     End of SLAMCH
        !
    END FUNCTION slamch

    INTEGER FUNCTION isamax(N,SX,INCX)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        INTEGER INCX,N
        !     ..
        !     .. Array Arguments ..
        REAL SX(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        REAL SMAX
        INTEGER I,IX
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC abs
        !     ..
        isamax = 0
        IF (n.LT.1 .OR. incx.LE.0) RETURN
        isamax = 1
        IF (n.EQ.1) RETURN
        IF (incx.EQ.1) THEN
            !
            !        code for increment equal to 1
            !
            smax = abs(sx(1))
            DO i = 2,n
                IF (abs(sx(i)).GT.smax) THEN
                    isamax = i
                    smax = abs(sx(i))
                END IF
            END DO
        ELSE
            !
            !        code for increment not equal to 1
            !
            ix = 1
            smax = abs(sx(1))
            ix = ix + incx
            DO i = 2,n
                IF (abs(sx(ix)).GT.smax) THEN
                    isamax = i
                    smax = abs(sx(ix))
                END IF
                ix = ix + incx
            END DO
        END IF
        RETURN
    END FUNCTION isamax

    REAL FUNCTION sdot(N,SX,INCX,SY,INCY)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        INTEGER INCX,INCY,N
        !     ..
        !     .. Array Arguments ..
        REAL SX(*),SY(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        REAL STEMP
        INTEGER I,IX,IY,M,MP1
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC mod
        !     ..
        stemp = 0.0e0
        sdot = 0.0e0
        IF (n.LE.0) RETURN
        IF (incx.EQ.1 .AND. incy.EQ.1) THEN
            !
            !        code for both increments equal to 1
            !
            !
            !        clean-up loop
            !
            m = mod(n,5)
            IF (m.NE.0) THEN
                DO i = 1,m
                    stemp = stemp + sx(i)*sy(i)
                END DO
                IF (n.LT.5) THEN
                    sdot=stemp
                    RETURN
                END IF
            END IF
            mp1 = m + 1
            DO i = mp1,n,5
                stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) + &
                    sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
            END DO
        ELSE
            !
            !        code for unequal increments or equal increments
            !          not equal to 1
            !
            ix = 1
            iy = 1
            IF (incx.LT.0) ix = (-n+1)*incx + 1
            IF (incy.LT.0) iy = (-n+1)*incy + 1
            DO i = 1,n
                stemp = stemp + sx(ix)*sy(iy)
                ix = ix + incx
                iy = iy + incy
            END DO
        END IF
        sdot = stemp
        RETURN
    END FUNCTION sdot

    INTEGER FUNCTION ilaslc( M, N, A, LDA )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2017
        !
        !     .. Scalar Arguments ..
        INTEGER            M, N, LDA
        !     ..
        !     .. Array Arguments ..
        REAL               A( lda, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL             ZERO
        parameter( zero = 0.0e+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER I
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick test for the common case where one corner is non-zero.
        IF( n.EQ.0 ) THEN
            ilaslc = n
        ELSE IF( a(1, n).NE.zero .OR. a(m, n).NE.zero ) THEN
            ilaslc = n
        ELSE
            !     Now scan each column from the end, returning with the first non-zero.
            DO ilaslc = n, 1, -1
                DO i = 1, m
                    IF( a(i, ilaslc).NE.zero ) RETURN
                END DO
            END DO
        END IF
        RETURN
    END FUNCTION ilaslc

    REAL FUNCTION snrm2(N,X,INCX)
        !
        !  -- Reference BLAS level1 routine (version 3.8.0) --
        !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2017
        !
        !     .. Scalar Arguments ..
        INTEGER INCX,N
        !     ..
        !     .. Array Arguments ..
        REAL X(*)
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL ONE,ZERO
        parameter(one=1.0e+0,zero=0.0e+0)
        !     ..
        !     .. Local Scalars ..
        REAL ABSXI,NORM,SCALE,SSQ
        INTEGER IX
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC abs,sqrt
        !     ..
        IF (n.LT.1 .OR. incx.LT.1) THEN
            norm = zero
        ELSE IF (n.EQ.1) THEN
            norm = abs(x(1))
        ELSE
            scale = zero
            ssq = one
            !        The following loop is equivalent to this call to the LAPACK
            !        auxiliary routine:
            !        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
            !
            DO ix = 1,1 + (n-1)*incx,incx
                IF (x(ix).NE.zero) THEN
                    absxi = abs(x(ix))
                    IF (scale.LT.absxi) THEN
                        ssq = one + ssq* (scale/absxi)**2
                        scale = absxi
                    ELSE
                        ssq = ssq + (absxi/scale)**2
                    END IF
                END IF
            END DO
            norm = scale*sqrt(ssq)
        END IF
        !
        snrm2 = norm
        RETURN
        !
        !     End of SNRM2.
        !
    END FUNCTION snrm2

    LOGICAL FUNCTION sisnan( SIN )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2017
        !
        !     .. Scalar Arguments ..
        REAL, INTENT(IN) :: SIN
        !     ..
        !
        !  =====================================================================
        !
        !  .. Executable Statements ..
        sisnan = slaisnan(sin,sin)
        RETURN
    END FUNCTION sisnan

    LOGICAL FUNCTION slaisnan( SIN1, SIN2 )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2017
        !
        !     .. Scalar Arguments ..
        REAL, INTENT(IN) :: SIN1, SIN2
        !     ..
        !
        !  =====================================================================
        !
        !  .. Executable Statements ..
        slaisnan = (sin1.NE.sin2)
        RETURN
    END FUNCTION slaisnan

    REAL             FUNCTION slapy2( X, Y )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2017
        !
        !     .. Scalar Arguments ..
        REAL               X, Y
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO
        parameter( zero = 0.0e0 )
        REAL               ONE
        parameter( one = 1.0e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               W, XABS, YABS, Z
        LOGICAL            X_IS_NAN, Y_IS_NAN
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          abs, max, min, sqrt
        !     ..
        !     .. Executable Statements ..
        !
        !     ..
        !     .. Executable Statements ..
        !
        x_is_nan = sisnan( x )
        y_is_nan = sisnan( y )
        IF ( x_is_nan ) slapy2 = x
        IF ( y_is_nan ) slapy2 = y
        !
        IF ( .NOT.( x_is_nan.OR.y_is_nan ) ) THEN
            xabs = abs( x )
            yabs = abs( y )
            w = max( xabs, yabs )
            z = min( xabs, yabs )
            IF( z.EQ.zero ) THEN
                slapy2 = w
            ELSE
                slapy2 = w*sqrt( one+( z / w )**2 )
            END IF
        END IF
        RETURN
        !
        !     End of SLAPY2
        !
    END FUNCTION slapy2

    REAL             FUNCTION slange( NORM, M, N, A, LDA, WORK )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          NORM
        INTEGER            LDA, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( lda, * ), WORK( * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE, ZERO
        parameter( one = 1.0e+0, zero = 0.0e+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            I, J
        REAL               SCALE, SUM, VALUE, TEMP
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          abs, min, sqrt
        !     ..
        !     .. Executable Statements ..
        !
        IF( min( m, n ).EQ.0 ) THEN
            VALUE = zero
        ELSE IF( lsame( norm, 'M' ) ) THEN
            !
            !        Find max(abs(A(i,j))).
            !
            VALUE = zero
            DO j = 1, n
                DO i = 1, m
                    temp = abs( a( i, j ) )
                    IF( VALUE.LT.temp .OR. sisnan( temp ) ) VALUE = temp
                END DO
            END DO
        ELSE IF( ( lsame( norm, 'O' ) ) .OR. ( norm.EQ.'1' ) ) THEN
            !
            !        Find norm1(A).
            !
            VALUE = zero
            DO j = 1, n
                sum = zero
                DO i = 1, m
                    sum = sum + abs( a( i, j ) )
                END DO
                IF( VALUE.LT.sum .OR. sisnan( sum ) ) VALUE = sum
            END DO
        ELSE IF( lsame( norm, 'I' ) ) THEN
            !
            !        Find normI(A).
            !
            DO i = 1, m
                work( i ) = zero
            END DO
            DO  j = 1, n
                DO  i = 1, m
                    work( i ) = work( i ) + abs( a( i, j ) )
                END DO
            END DO
            VALUE = zero
            DO i = 1, m
                temp = work( i )
                IF( VALUE.LT.temp .OR. sisnan( temp ) ) VALUE = temp
            END DO
        ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN
            !
            !        Find normF(A).
            !
            scale = zero
            sum = one
            DO j = 1, n
                CALL slassq( m, a( 1, j ), 1, scale, sum )
            END DO
            VALUE = scale*sqrt( sum )
        END IF
        !
        slange = VALUE
        RETURN
        !
        !     End of SLANGE
        !
    END FUNCTION slange

    SUBROUTINE slassq( N, X, INCX, SCALE, SUMSQ )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            INCX, N
        REAL               SCALE, SUMSQ
        !     ..
        !     .. Array Arguments ..
        REAL               X( * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO
        parameter( zero = 0.0e+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            IX
        REAL               ABSXI
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          abs
        !     ..
        !     .. Executable Statements ..
        !
        IF( n.GT.0 ) THEN
            DO ix = 1, 1 + ( n-1 )*incx, incx
                absxi = abs( x( ix ) )
                IF( absxi.GT.zero.OR.sisnan( absxi ) ) THEN
                    IF( scale.LT.absxi ) THEN
                        sumsq = 1 + sumsq*( scale / absxi )**2
                        scale = absxi
                    ELSE
                        sumsq = sumsq + ( absxi / scale )**2
                    END IF
                END IF
            END DO
        END IF
        RETURN
        !
        !     End of SLASSQ
        !
    END SUBROUTINE slassq

    INTEGER FUNCTION ilaslr( M, N, A, LDA )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            M, N, LDA
        !     ..
        !     .. Array Arguments ..
        REAL               A( lda, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL             ZERO
        parameter( zero = 0.0e+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER I, J
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick test for the common case where one corner is non-zero.
        IF( m.EQ.0 ) THEN
            ilaslr = m
        ELSEIF( a(m, 1).NE.zero .OR. a(m, n).NE.zero ) THEN
            ilaslr = m
        ELSE
            !     Scan up each column tracking the last zero row seen.
            ilaslr = 0
            DO j = 1, n
                i=m
                DO WHILE((a(max(i,1),j).EQ.zero).AND.(i.GE.1))
                    i=i-1
                ENDDO
                ilaslr = max( ilaslr, i )
            END DO
        END IF
        RETURN
    END FUNCTION ilaslr

    SUBROUTINE slaexc( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, &
        INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        LOGICAL            WANTQ
        INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
        !     ..
        !     .. Array Arguments ..
        REAL               Q( ldq, * ), T( ldt, * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        parameter( zero = 0.0e+0, one = 1.0e+0 )
        REAL               TEN
        parameter( ten = 1.0e+1 )
        INTEGER            LDD, LDX
        parameter( ldd = 4, ldx = 2 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            IERR, J2, J3, J4, K, ND
        REAL               CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22, &
            t33, tau, tau1, tau2, temp, thresh, wi1, wi2, &
            wr1, wr2, xnorm
        !     ..
        !     .. Local Arrays ..
        REAL               D( ldd, 4 ), U( 3 ), U1( 3 ), U2( 3 ), &
            x( ldx, 2 )
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          abs, max
        !     ..
        !     .. Executable Statements ..
        !
        info = 0
        !
        !     Quick return if possible
        !
        IF( n.EQ.0 .OR. n1.EQ.0 .OR. n2.EQ.0 ) &
            RETURN
        IF( j1+n1.GT.n ) &
            RETURN
        !
        j2 = j1 + 1
        j3 = j1 + 2
        j4 = j1 + 3
        !
        IF( n1.EQ.1 .AND. n2.EQ.1 ) THEN
            !
            !        Swap two 1-by-1 blocks.
            !
            t11 = t( j1, j1 )
            t22 = t( j2, j2 )
            !
            !        Determine the transformation to perform the interchange.
            !
            CALL slartg( t( j1, j2 ), t22-t11, cs, sn, temp )
            !
            !        Apply transformation to the matrix T.
            !
            IF( j3.LE.n ) &
                CALL srot( n-j1-1, t( j1, j3 ), ldt, t( j2, j3 ), ldt, cs, &
                sn )
            CALL srot( j1-1, t( 1, j1 ), 1, t( 1, j2 ), 1, cs, sn )
            !
            t( j1, j1 ) = t22
            t( j2, j2 ) = t11
            !
            IF( wantq ) THEN
                !
                !           Accumulate transformation in the matrix Q.
                !
                CALL srot( n, q( 1, j1 ), 1, q( 1, j2 ), 1, cs, sn )
            END IF
            !
        ELSE
            !
            !        Swapping involves at least one 2-by-2 block.
            !
            !        Copy the diagonal block of order N1+N2 to the local array D
            !        and compute its norm.
            !
            nd = n1 + n2
            CALL slacpy( 'Full', nd, nd, t( j1, j1 ), ldt, d, ldd )
            dnorm = slange( 'Max', nd, nd, d, ldd, work )
            !
            !        Compute machine-dependent threshold for test for accepting
            !        swap.
            !
            eps = slamch( 'P' )
            smlnum = slamch( 'S' ) / eps
            thresh = max( ten*eps*dnorm, smlnum )
            !
            !        Solve T11*X - X*T22 = scale*T12 for X.
            !
            CALL slasy2( .false., .false., -1, n1, n2, d, ldd, &
                d( n1+1, n1+1 ), ldd, d( 1, n1+1 ), ldd, scale, x, &
                ldx, xnorm, ierr )
            !
            !        Swap the adjacent diagonal blocks.
            !
            k = n1 + n1 + n2 - 3
            !!! GO TO ( 10, 20, 30 )k -> obsolecent feature; changed manually
            IF (k == 1) THEN
                GO TO 10
            ELSE IF (k == 2) THEN
                GO TO 20
            ELSE
                GO TO 30
            END IF
            !
10          CONTINUE
            !
            !        N1 = 1, N2 = 2: generate elementary reflector H so that:
            !
            !        ( scale, X11, X12 ) H = ( 0, 0, * )
            !
            u( 1 ) = scale
            u( 2 ) = x( 1, 1 )
            u( 3 ) = x( 1, 2 )
            CALL slarfg( 3, u( 3 ), u, 1, tau )
            u( 3 ) = one
            t11 = t( j1, j1 )
            !
            !        Perform swap provisionally on diagonal block in D.
            !
            CALL slarfx( 'L', 3, 3, u, tau, d, ldd, work )
            CALL slarfx( 'R', 3, 3, u, tau, d, ldd, work )
            !
            !        Test whether to reject swap.
            !
            IF( max( abs( d( 3, 1 ) ), abs( d( 3, 2 ) ), abs( d( 3, &
                3 )-t11 ) ).GT.thresh )GO TO 50
            !
            !        Accept swap: apply transformation to the entire matrix T.
            !
            CALL slarfx( 'L', 3, n-j1+1, u, tau, t( j1, j1 ), ldt, work )
            CALL slarfx( 'R', j2, 3, u, tau, t( 1, j1 ), ldt, work )
            !
            t( j3, j1 ) = zero
            t( j3, j2 ) = zero
            t( j3, j3 ) = t11
            !
            IF( wantq ) THEN
                !
                !           Accumulate transformation in the matrix Q.
                !
                CALL slarfx( 'R', n, 3, u, tau, q( 1, j1 ), ldq, work )
            END IF
            GO TO 40
            !
20          CONTINUE
            !
            !        N1 = 2, N2 = 1: generate elementary reflector H so that:
            !
            !        H (  -X11 ) = ( * )
            !          (  -X21 ) = ( 0 )
            !          ( scale ) = ( 0 )
            !
            u( 1 ) = -x( 1, 1 )
            u( 2 ) = -x( 2, 1 )
            u( 3 ) = scale
            CALL slarfg( 3, u( 1 ), u( 2 ), 1, tau )
            u( 1 ) = one
            t33 = t( j3, j3 )
            !
            !        Perform swap provisionally on diagonal block in D.
            !
            CALL slarfx( 'L', 3, 3, u, tau, d, ldd, work )
            CALL slarfx( 'R', 3, 3, u, tau, d, ldd, work )
            !
            !        Test whether to reject swap.
            !
            IF( max( abs( d( 2, 1 ) ), abs( d( 3, 1 ) ), abs( d( 1, &
                1 )-t33 ) ).GT.thresh )GO TO 50
            !
            !        Accept swap: apply transformation to the entire matrix T.
            !
            CALL slarfx( 'R', j3, 3, u, tau, t( 1, j1 ), ldt, work )
            CALL slarfx( 'L', 3, n-j1, u, tau, t( j1, j2 ), ldt, work )
            !
            t( j1, j1 ) = t33
            t( j2, j1 ) = zero
            t( j3, j1 ) = zero
            !
            IF( wantq ) THEN
                !
                !           Accumulate transformation in the matrix Q.
                !
                CALL slarfx( 'R', n, 3, u, tau, q( 1, j1 ), ldq, work )
            END IF
            GO TO 40
            !
30          CONTINUE
            !
            !        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
            !        that:
            !
            !        H(2) H(1) (  -X11  -X12 ) = (  *  * )
            !                  (  -X21  -X22 )   (  0  * )
            !                  ( scale    0  )   (  0  0 )
            !                  (    0  scale )   (  0  0 )
            !
            u1( 1 ) = -x( 1, 1 )
            u1( 2 ) = -x( 2, 1 )
            u1( 3 ) = scale
            CALL slarfg( 3, u1( 1 ), u1( 2 ), 1, tau1 )
            u1( 1 ) = one
            !
            temp = -tau1*( x( 1, 2 )+u1( 2 )*x( 2, 2 ) )
            u2( 1 ) = -temp*u1( 2 ) - x( 2, 2 )
            u2( 2 ) = -temp*u1( 3 )
            u2( 3 ) = scale
            CALL slarfg( 3, u2( 1 ), u2( 2 ), 1, tau2 )
            u2( 1 ) = one
            !
            !        Perform swap provisionally on diagonal block in D.
            !
            CALL slarfx( 'L', 3, 4, u1, tau1, d, ldd, work )
            CALL slarfx( 'R', 4, 3, u1, tau1, d, ldd, work )
            CALL slarfx( 'L', 3, 4, u2, tau2, d( 2, 1 ), ldd, work )
            CALL slarfx( 'R', 4, 3, u2, tau2, d( 1, 2 ), ldd, work )
            !
            !        Test whether to reject swap.
            !
            IF( max( abs( d( 3, 1 ) ), abs( d( 3, 2 ) ), abs( d( 4, 1 ) ), &
                abs( d( 4, 2 ) ) ).GT.thresh )GO TO 50
            !
            !        Accept swap: apply transformation to the entire matrix T.
            !
            CALL slarfx( 'L', 3, n-j1+1, u1, tau1, t( j1, j1 ), ldt, work )
            CALL slarfx( 'R', j4, 3, u1, tau1, t( 1, j1 ), ldt, work )
            CALL slarfx( 'L', 3, n-j1+1, u2, tau2, t( j2, j1 ), ldt, work )
            CALL slarfx( 'R', j4, 3, u2, tau2, t( 1, j2 ), ldt, work )
            !
            t( j3, j1 ) = zero
            t( j3, j2 ) = zero
            t( j4, j1 ) = zero
            t( j4, j2 ) = zero
            !
            IF( wantq ) THEN
                !
                !           Accumulate transformation in the matrix Q.
                !
                CALL slarfx( 'R', n, 3, u1, tau1, q( 1, j1 ), ldq, work )
                CALL slarfx( 'R', n, 3, u2, tau2, q( 1, j2 ), ldq, work )
            END IF
            !
40          CONTINUE
            !
            IF( n2.EQ.2 ) THEN
                !
                !           Standardize new 2-by-2 block T11
                !
                CALL slanv2( t( j1, j1 ), t( j1, j2 ), t( j2, j1 ), &
                    t( j2, j2 ), wr1, wi1, wr2, wi2, cs, sn )
                CALL srot( n-j1-1, t( j1, j1+2 ), ldt, t( j2, j1+2 ), ldt, &
                    cs, sn )
                CALL srot( j1-1, t( 1, j1 ), 1, t( 1, j2 ), 1, cs, sn )
                IF( wantq ) &
                    CALL srot( n, q( 1, j1 ), 1, q( 1, j2 ), 1, cs, sn )
            END IF
            !
            IF( n1.EQ.2 ) THEN
                !
                !           Standardize new 2-by-2 block T22
                !
                j3 = j1 + n2
                j4 = j3 + 1
                CALL slanv2( t( j3, j3 ), t( j3, j4 ), t( j4, j3 ), &
                    t( j4, j4 ), wr1, wi1, wr2, wi2, cs, sn )
                IF( j3+2.LE.n ) &
                    CALL srot( n-j3-1, t( j3, j3+2 ), ldt, t( j4, j3+2 ), &
                    ldt, cs, sn )
                CALL srot( j3-1, t( 1, j3 ), 1, t( 1, j4 ), 1, cs, sn )
                IF( wantq ) &
                    CALL srot( n, q( 1, j3 ), 1, q( 1, j4 ), 1, cs, sn )
            END IF
            !
        END IF
        RETURN
        !
        !     Exit with INFO = 1 if swap was rejected.
        !
50      info = 1
        RETURN
        !
        !     End of SLAEXC
        !
    END SUBROUTINE slaexc

    SUBROUTINE SLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, &
        LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2016
        !
        !     .. Scalar Arguments ..
        LOGICAL            LTRANL, LTRANR
        INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
        REAL               SCALE, XNORM
        !     ..
        !     .. Array Arguments ..
        REAL               B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ), &
            X( LDX, * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
        REAL               TWO, HALF, EIGHT
        PARAMETER          ( TWO = 2.0E+0, HALF = 0.5E+0, EIGHT = 8.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            BSWAP, XSWAP
        INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
        REAL               BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1, &
            TEMP, U11, U12, U22, XMAX
        !     ..
        !     .. Local Arrays ..
        LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
        INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ), &
            LOCU22( 4 )
        REAL               BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX
        !     ..
        !     .. Data statements ..
        DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / , &
            LOCU22 / 4, 3, 2, 1 /
        DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
        DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
        !     ..
        !     .. Executable Statements ..
        !
        !     Do not check the input parameters for errors
        !
        INFO = 0
        !
        !     Quick return if possible
        !
        IF( N1.EQ.0 .OR. N2.EQ.0 ) &
            RETURN
        !
        !     Set constants to control overflow
        !
        EPS = SLAMCH( 'P' )
        SMLNUM = SLAMCH( 'S' ) / EPS
        SGN = ISGN
        !
        K = N1 + N1 + N2 - 2
        !!! GO TO ( 10, 20, 30, 50 )K -> obsolecent feature; changed manually
        IF (K == 1) THEN
            GO TO 10
        ELSE IF (K == 2) THEN
            GO To 20
        ELSE IF (K == 3) THEN
            GO TO 30
        ELSE
            GO TO 50
        END IF
        !
        !     1 by 1: TL11*X + SGN*X*TR11 = B11
        !
10      CONTINUE
        TAU1 = TL( 1, 1 ) + SGN*TR( 1, 1 )
        BET = ABS( TAU1 )
        IF( BET.LE.SMLNUM ) THEN
            TAU1 = SMLNUM
            BET = SMLNUM
            INFO = 1
        END IF
        !
        SCALE = ONE
        GAM = ABS( B( 1, 1 ) )
        IF( SMLNUM*GAM.GT.BET ) &
            SCALE = ONE / GAM
        !
        X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
        XNORM = ABS( X( 1, 1 ) )
        RETURN
        !
        !     1 by 2:
        !     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
        !                                       [TR21 TR22]
        !
20      CONTINUE
        !
        SMIN = MAX( EPS*MAX( ABS( TL( 1, 1 ) ), ABS( TR( 1, 1 ) ), &
            ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) ), &
            SMLNUM )
        TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
        TMP( 4 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
        IF( LTRANR ) THEN
            TMP( 2 ) = SGN*TR( 2, 1 )
            TMP( 3 ) = SGN*TR( 1, 2 )
        ELSE
            TMP( 2 ) = SGN*TR( 1, 2 )
            TMP( 3 ) = SGN*TR( 2, 1 )
        END IF
        BTMP( 1 ) = B( 1, 1 )
        BTMP( 2 ) = B( 1, 2 )
        GO TO 40
        !
        !     2 by 1:
        !          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
        !            [TL21 TL22] [X21]         [X21]         [B21]
        !
30      CONTINUE
        SMIN = MAX( EPS*MAX( ABS( TR( 1, 1 ) ), ABS( TL( 1, 1 ) ), &
            ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) ), &
            SMLNUM )
        TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
        TMP( 4 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
        IF( LTRANL ) THEN
            TMP( 2 ) = TL( 1, 2 )
            TMP( 3 ) = TL( 2, 1 )
        ELSE
            TMP( 2 ) = TL( 2, 1 )
            TMP( 3 ) = TL( 1, 2 )
        END IF
        BTMP( 1 ) = B( 1, 1 )
        BTMP( 2 ) = B( 2, 1 )
40      CONTINUE
        !
        !     Solve 2 by 2 system using complete pivoting.
        !     Set pivots less than SMIN to SMIN.
        !
        IPIV = ISAMAX( 4, TMP, 1 )
        U11 = TMP( IPIV )
        IF( ABS( U11 ).LE.SMIN ) THEN
            INFO = 1
            U11 = SMIN
        END IF
        U12 = TMP( LOCU12( IPIV ) )
        L21 = TMP( LOCL21( IPIV ) ) / U11
        U22 = TMP( LOCU22( IPIV ) ) - U12*L21
        XSWAP = XSWPIV( IPIV )
        BSWAP = BSWPIV( IPIV )
        IF( ABS( U22 ).LE.SMIN ) THEN
            INFO = 1
            U22 = SMIN
        END IF
        IF( BSWAP ) THEN
            TEMP = BTMP( 2 )
            BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
            BTMP( 1 ) = TEMP
        ELSE
            BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
        END IF
        SCALE = ONE
        IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( U22 ) .OR. &
            ( TWO*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( U11 ) ) THEN
            SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
            BTMP( 1 ) = BTMP( 1 )*SCALE
            BTMP( 2 ) = BTMP( 2 )*SCALE
        END IF
        X2( 2 ) = BTMP( 2 ) / U22
        X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
        IF( XSWAP ) THEN
            TEMP = X2( 2 )
            X2( 2 ) = X2( 1 )
            X2( 1 ) = TEMP
        END IF
        X( 1, 1 ) = X2( 1 )
        IF( N1.EQ.1 ) THEN
            X( 1, 2 ) = X2( 2 )
            XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
        ELSE
            X( 2, 1 ) = X2( 2 )
            XNORM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
        END IF
        RETURN
        !
        !     2 by 2:
        !     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
        !       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
        !
        !     Solve equivalent 4 by 4 system using complete pivoting.
        !     Set pivots less than SMIN to SMIN.
        !
50      CONTINUE
        SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ), &
            ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
        SMIN = MAX( SMIN, ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ), &
            ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
        SMIN = MAX( EPS*SMIN, SMLNUM )
        BTMP( 1 ) = ZERO
        CALL SCOPY( 16, BTMP, 0, T16, 1 )
        T16( 1, 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
        T16( 2, 2 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
        T16( 3, 3 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
        T16( 4, 4 ) = TL( 2, 2 ) + SGN*TR( 2, 2 )
        IF( LTRANL ) THEN
            T16( 1, 2 ) = TL( 2, 1 )
            T16( 2, 1 ) = TL( 1, 2 )
            T16( 3, 4 ) = TL( 2, 1 )
            T16( 4, 3 ) = TL( 1, 2 )
        ELSE
            T16( 1, 2 ) = TL( 1, 2 )
            T16( 2, 1 ) = TL( 2, 1 )
            T16( 3, 4 ) = TL( 1, 2 )
            T16( 4, 3 ) = TL( 2, 1 )
        END IF
        IF( LTRANR ) THEN
            T16( 1, 3 ) = SGN*TR( 1, 2 )
            T16( 2, 4 ) = SGN*TR( 1, 2 )
            T16( 3, 1 ) = SGN*TR( 2, 1 )
            T16( 4, 2 ) = SGN*TR( 2, 1 )
        ELSE
            T16( 1, 3 ) = SGN*TR( 2, 1 )
            T16( 2, 4 ) = SGN*TR( 2, 1 )
            T16( 3, 1 ) = SGN*TR( 1, 2 )
            T16( 4, 2 ) = SGN*TR( 1, 2 )
        END IF
        BTMP( 1 ) = B( 1, 1 )
        BTMP( 2 ) = B( 2, 1 )
        BTMP( 3 ) = B( 1, 2 )
        BTMP( 4 ) = B( 2, 2 )
        !
        !     Perform elimination
        !
        DO  I = 1, 3
            XMAX = ZERO
            DO  IP = I, 4
                DO  JP = I, 4
                    IF( ABS( T16( IP, JP ) ).GE.XMAX ) THEN
                        XMAX = ABS( T16( IP, JP ) )
                        IPSV = IP
                        JPSV = JP
                    END IF
                END DO
            END DO
            IF( IPSV.NE.I ) THEN
                CALL SSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
                TEMP = BTMP( I )
                BTMP( I ) = BTMP( IPSV )
                BTMP( IPSV ) = TEMP
            END IF
            IF( JPSV.NE.I ) &
                CALL SSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
            JPIV( I ) = JPSV
            IF( ABS( T16( I, I ) ).LT.SMIN ) THEN
                INFO = 1
                T16( I, I ) = SMIN
            END IF
            DO  J = I + 1, 4
                T16( J, I ) = T16( J, I ) / T16( I, I )
                BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
                DO  K = I + 1, 4
                    T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
                END DO
            END DO
        END DO
        IF( ABS( T16( 4, 4 ) ).LT.SMIN ) THEN
            INFO = 1
            T16( 4, 4 ) = SMIN
        END IF
        SCALE = ONE
        IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) ).GT.ABS( T16( 1, 1 ) ) .OR. &
            ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) ).GT.ABS( T16( 2, 2 ) ) .OR. &
            ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) ).GT.ABS( T16( 3, 3 ) ) .OR. &
            ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) ).GT.ABS( T16( 4, 4 ) ) ) THEN
            SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ), &
                ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ), ABS( BTMP( 4 ) ) )
            BTMP( 1 ) = BTMP( 1 )*SCALE
            BTMP( 2 ) = BTMP( 2 )*SCALE
            BTMP( 3 ) = BTMP( 3 )*SCALE
            BTMP( 4 ) = BTMP( 4 )*SCALE
        END IF
        DO  I = 1, 4
            K = 5 - I
            TEMP = ONE / T16( K, K )
            TMP( K ) = BTMP( K )*TEMP
            DO  J = K + 1, 4
                TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
            END DO
        END DO
        DO  I = 1, 3
            IF( JPIV( 4-I ).NE.4-I ) THEN
                TEMP = TMP( 4-I )
                TMP( 4-I ) = TMP( JPIV( 4-I ) )
                TMP( JPIV( 4-I ) ) = TEMP
            END IF
        END DO
        X( 1, 1 ) = TMP( 1 )
        X( 2, 1 ) = TMP( 2 )
        X( 1, 2 ) = TMP( 3 )
        X( 2, 2 ) = TMP( 4 )
        XNORM = MAX( ABS( TMP( 1 ) )+ABS( TMP( 3 ) ), &
            ABS( TMP( 2 ) )+ABS( TMP( 4 ) ) )
        RETURN
        !
        !     End of SLASY2
        !
    END SUBROUTINE SLASY2

    SUBROUTINE SLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          SIDE
        INTEGER            LDC, M, N
        REAL               TAU
        !     ..
        !     .. Array Arguments ..
        REAL               C( LDC, * ), V( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        INTEGER            J
        REAL               SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, &
            V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
        !     ..
        !     .. Executable Statements ..
        !
        IF( TAU.EQ.ZERO ) &
            RETURN
        IF( LSAME( SIDE, 'L' ) ) THEN
            !
            !        Form  H * C, where H has order m.
            !
            !!!GO TO ( 10, 30, 50, 70, 90, 110, 130, 150, &
            !!!    170, 190 )M
            !!! -> obsolecent feature; changed manually
            IF (M == 1) THEN
                GO TO 10
            ELSE IF (M == 2) THEN
                GO TO 30                
            ELSE IF (M == 3) THEN
                GO TO 50                
            ELSE IF (M == 4) THEN
                GO TO 70                
            ELSE IF (M == 5) THEN
                GO TO 90                
            ELSE IF (M == 6) THEN
                GO TO 110                
            ELSE IF (M == 7) THEN
                GO TO 130                
            ELSE IF (M == 8) THEN
                GO TO 150                
            ELSE IF (M == 9) THEN
                GO TO 170                
            ELSE
                GO TO 190
            END IF
                
            !
            !        Code for general M
            !
            CALL SLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
            GO TO 410
10          CONTINUE
            !
            !        Special code for 1 x 1 Householder
            !
            T1 = ONE - TAU*V( 1 )*V( 1 )
            DO  J = 1, N
                C( 1, J ) = T1*C( 1, J )
            END DO
            GO TO 410
30          CONTINUE
            !
            !        Special code for 2 x 2 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
            END DO
            GO TO 410
50          CONTINUE
            !
            !        Special code for 3 x 3 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
            END DO
            GO TO 410
70          CONTINUE
            !
            !        Special code for 4 x 4 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                    V4*C( 4, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
                C( 4, J ) = C( 4, J ) - SUM*T4
            END DO
            GO TO 410
90          CONTINUE
            !
            !        Special code for 5 x 5 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                    V4*C( 4, J ) + V5*C( 5, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
                C( 4, J ) = C( 4, J ) - SUM*T4
                C( 5, J ) = C( 5, J ) - SUM*T5
            END DO
            GO TO 410
110         CONTINUE
            !
            !        Special code for 6 x 6 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                    V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
                C( 4, J ) = C( 4, J ) - SUM*T4
                C( 5, J ) = C( 5, J ) - SUM*T5
                C( 6, J ) = C( 6, J ) - SUM*T6
            END DO
            GO TO 410
130         CONTINUE
            !
            !        Special code for 7 x 7 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                    V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                    V7*C( 7, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
                C( 4, J ) = C( 4, J ) - SUM*T4
                C( 5, J ) = C( 5, J ) - SUM*T5
                C( 6, J ) = C( 6, J ) - SUM*T6
                C( 7, J ) = C( 7, J ) - SUM*T7
            END DO
            GO TO 410
150         CONTINUE
            !
            !        Special code for 8 x 8 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            V8 = V( 8 )
            T8 = TAU*V8
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                    V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                    V7*C( 7, J ) + V8*C( 8, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
                C( 4, J ) = C( 4, J ) - SUM*T4
                C( 5, J ) = C( 5, J ) - SUM*T5
                C( 6, J ) = C( 6, J ) - SUM*T6
                C( 7, J ) = C( 7, J ) - SUM*T7
                C( 8, J ) = C( 8, J ) - SUM*T8
            END DO
            GO TO 410
170         CONTINUE
            !
            !        Special code for 9 x 9 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            V8 = V( 8 )
            T8 = TAU*V8
            V9 = V( 9 )
            T9 = TAU*V9
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                    V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                    V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
                C( 4, J ) = C( 4, J ) - SUM*T4
                C( 5, J ) = C( 5, J ) - SUM*T5
                C( 6, J ) = C( 6, J ) - SUM*T6
                C( 7, J ) = C( 7, J ) - SUM*T7
                C( 8, J ) = C( 8, J ) - SUM*T8
                C( 9, J ) = C( 9, J ) - SUM*T9
            END DO
            GO TO 410
190         CONTINUE
            !
            !        Special code for 10 x 10 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            V8 = V( 8 )
            T8 = TAU*V8
            V9 = V( 9 )
            T9 = TAU*V9
            V10 = V( 10 )
            T10 = TAU*V10
            DO  J = 1, N
                SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                    V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                    V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) + &
                    V10*C( 10, J )
                C( 1, J ) = C( 1, J ) - SUM*T1
                C( 2, J ) = C( 2, J ) - SUM*T2
                C( 3, J ) = C( 3, J ) - SUM*T3
                C( 4, J ) = C( 4, J ) - SUM*T4
                C( 5, J ) = C( 5, J ) - SUM*T5
                C( 6, J ) = C( 6, J ) - SUM*T6
                C( 7, J ) = C( 7, J ) - SUM*T7
                C( 8, J ) = C( 8, J ) - SUM*T8
                C( 9, J ) = C( 9, J ) - SUM*T9
                C( 10, J ) = C( 10, J ) - SUM*T10
            END DO
            GO TO 410
        ELSE
            !
            !        Form  C * H, where H has order n.
            !
            !!! GO TO ( 210, 230, 250, 270, 290, 310, 330, 350, &
            !!!    370, 390 )N
            !!! -> obsolecent feature; changed manually
            If (N==1) THEN
                GO TO 210
            ELSE IF (N==2) THEN
                GO TO 230                
            ELSE IF (N==3) THEN
                GO TO 250                
            ELSE IF (N==4) THEN
                GO TO 270                
            ELSE IF (N==5) THEN
                GO TO 290                
            ELSE IF (N==6) THEN
                GO TO 310                
            ELSE IF (N==7) THEN
                GO TO 330                
            ELSE IF (N==8) THEN
                GO TO 350                
            ELSE IF (N==9) THEN
                GO TO 370                
            ELSE
                GO TO 390                
            END IF
                
            !        Code for general N
            !
            CALL SLARF( SIDE, M, N, V, 1, TAU, C, LDC, WORK )
            GO TO 410
210         CONTINUE
            !
            !        Special code for 1 x 1 Householder
            !
            T1 = ONE - TAU*V( 1 )*V( 1 )
            DO  J = 1, M
                C( J, 1 ) = T1*C( J, 1 )
            END DO
            GO TO 410
230         CONTINUE
            !
            !        Special code for 2 x 2 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
            END DO
            GO TO 410
250         CONTINUE
            !
            !        Special code for 3 x 3 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
            END DO
            GO TO 410
270         CONTINUE
            !
            !        Special code for 4 x 4 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                    V4*C( J, 4 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
                C( J, 4 ) = C( J, 4 ) - SUM*T4
            END DO
            GO TO 410
290         CONTINUE
            !
            !        Special code for 5 x 5 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                    V4*C( J, 4 ) + V5*C( J, 5 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
                C( J, 4 ) = C( J, 4 ) - SUM*T4
                C( J, 5 ) = C( J, 5 ) - SUM*T5
            END DO
            GO TO 410
310         CONTINUE
            !
            !        Special code for 6 x 6 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                    V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
                C( J, 4 ) = C( J, 4 ) - SUM*T4
                C( J, 5 ) = C( J, 5 ) - SUM*T5
                C( J, 6 ) = C( J, 6 ) - SUM*T6
            END DO
            GO TO 410
330         CONTINUE
            !
            !        Special code for 7 x 7 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                    V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                    V7*C( J, 7 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
                C( J, 4 ) = C( J, 4 ) - SUM*T4
                C( J, 5 ) = C( J, 5 ) - SUM*T5
                C( J, 6 ) = C( J, 6 ) - SUM*T6
                C( J, 7 ) = C( J, 7 ) - SUM*T7
            END DO
            GO TO 410
350         CONTINUE
            !
            !        Special code for 8 x 8 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            V8 = V( 8 )
            T8 = TAU*V8
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                    V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                    V7*C( J, 7 ) + V8*C( J, 8 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
                C( J, 4 ) = C( J, 4 ) - SUM*T4
                C( J, 5 ) = C( J, 5 ) - SUM*T5
                C( J, 6 ) = C( J, 6 ) - SUM*T6
                C( J, 7 ) = C( J, 7 ) - SUM*T7
                C( J, 8 ) = C( J, 8 ) - SUM*T8
            END DO
            GO TO 410
370         CONTINUE
            !
            !        Special code for 9 x 9 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            V8 = V( 8 )
            T8 = TAU*V8
            V9 = V( 9 )
            T9 = TAU*V9
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                    V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                    V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
                C( J, 4 ) = C( J, 4 ) - SUM*T4
                C( J, 5 ) = C( J, 5 ) - SUM*T5
                C( J, 6 ) = C( J, 6 ) - SUM*T6
                C( J, 7 ) = C( J, 7 ) - SUM*T7
                C( J, 8 ) = C( J, 8 ) - SUM*T8
                C( J, 9 ) = C( J, 9 ) - SUM*T9
            END DO
            GO TO 410
390         CONTINUE
            !
            !        Special code for 10 x 10 Householder
            !
            V1 = V( 1 )
            T1 = TAU*V1
            V2 = V( 2 )
            T2 = TAU*V2
            V3 = V( 3 )
            T3 = TAU*V3
            V4 = V( 4 )
            T4 = TAU*V4
            V5 = V( 5 )
            T5 = TAU*V5
            V6 = V( 6 )
            T6 = TAU*V6
            V7 = V( 7 )
            T7 = TAU*V7
            V8 = V( 8 )
            T8 = TAU*V8
            V9 = V( 9 )
            T9 = TAU*V9
            V10 = V( 10 )
            T10 = TAU*V10
            DO  J = 1, M
                SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                    V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                    V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) + &
                    V10*C( J, 10 )
                C( J, 1 ) = C( J, 1 ) - SUM*T1
                C( J, 2 ) = C( J, 2 ) - SUM*T2
                C( J, 3 ) = C( J, 3 ) - SUM*T3
                C( J, 4 ) = C( J, 4 ) - SUM*T4
                C( J, 5 ) = C( J, 5 ) - SUM*T5
                C( J, 6 ) = C( J, 6 ) - SUM*T6
                C( J, 7 ) = C( J, 7 ) - SUM*T7
                C( J, 8 ) = C( J, 8 ) - SUM*T8
                C( J, 9 ) = C( J, 9 ) - SUM*T9
                C( J, 10 ) = C( J, 10 ) - SUM*T10
            END DO
            GO TO 410
        END IF
410     RETURN
        !
        !     End of SLARFX
        !
    END SUBROUTINE SLARFX

    SUBROUTINE sormhr( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, &
        LDC, WORK, LWORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          SIDE, TRANS
        INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( lda, * ), C( ldc, * ), TAU( * ), &
            work( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        LOGICAL            LEFT, LQUERY
        INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          max, min
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        info = 0
        nh = ihi - ilo
        left = lsame( side, 'L' )
        lquery = ( lwork.EQ.-1 )
        !
        !     NQ is the order of Q and NW is the minimum dimension of WORK
        !
        IF( left ) THEN
            nq = m
            nw = n
        ELSE
            nq = n
            nw = m
        END IF
        IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
            info = -1
        ELSE IF( .NOT.lsame( trans, 'N' ) .AND. .NOT.lsame( trans, 'T' ) ) &
            THEN
            info = -2
        ELSE IF( m.LT.0 ) THEN
            info = -3
        ELSE IF( n.LT.0 ) THEN
            info = -4
        ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, nq ) ) THEN
            info = -5
        ELSE IF( ihi.LT.min( ilo, nq ) .OR. ihi.GT.nq ) THEN
            info = -6
        ELSE IF( lda.LT.max( 1, nq ) ) THEN
            info = -8
        ELSE IF( ldc.LT.max( 1, m ) ) THEN
            info = -11
        ELSE IF( lwork.LT.max( 1, nw ) .AND. .NOT.lquery ) THEN
            info = -13
        END IF
        !
        IF( info.EQ.0 ) THEN
            IF( left ) THEN
                nb = ilaenv( 1, 'SORMQR', side // trans, nh, n, nh, -1 )
            ELSE
                nb = ilaenv( 1, 'SORMQR', side // trans, m, nh, nh, -1 )
            END IF
            lwkopt = max( 1, nw )*nb
            work( 1 ) = lwkopt
        END IF
        !
        IF( info.NE.0 ) THEN
            CALL xerbla( 'SORMHR', -info )
            RETURN
        ELSE IF( lquery ) THEN
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( m.EQ.0 .OR. n.EQ.0 .OR. nh.EQ.0 ) THEN
            work( 1 ) = 1
            RETURN
        END IF
        !
        IF( left ) THEN
            mi = nh
            ni = n
            i1 = ilo + 1
            i2 = 1
        ELSE
            mi = m
            ni = nh
            i1 = 1
            i2 = ilo + 1
        END IF
        !
        CALL sormqr( side, trans, mi, ni, nh, a( ilo+1, ilo ), lda, &
            tau( ilo ), c( i1, i2 ), ldc, work, lwork, iinfo )
        !
        work( 1 ) = lwkopt
        RETURN
        !
        !     End of SORMHR
        !
    END SUBROUTINE sormhr

    SUBROUTINE SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
        WORK, LWORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          SIDE, TRANS
        INTEGER            INFO, K, LDA, LDC, LWORK, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), C( LDC, * ), TAU( * ), &
            WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        INTEGER            NBMAX, LDT, TSIZE
        PARAMETER          ( NBMAX = 64, LDT = NBMAX+1, &
            TSIZE = LDT*NBMAX )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            LEFT, LQUERY, NOTRAN
        INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK, &
            LWKOPT, MI, NB, NBMIN, NI, NQ, NW
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        INFO = 0
        LEFT = LSAME( SIDE, 'L' )
        NOTRAN = LSAME( TRANS, 'N' )
        LQUERY = ( LWORK.EQ.-1 )
        !
        !     NQ is the order of Q and NW is the minimum dimension of WORK
        !
        IF( LEFT ) THEN
            NQ = M
            NW = N
        ELSE
            NQ = N
            NW = M
        END IF
        IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
            INFO = -2
        ELSE IF( M.LT.0 ) THEN
            INFO = -3
        ELSE IF( N.LT.0 ) THEN
            INFO = -4
        ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
            INFO = -5
        ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
            INFO = -7
        ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
            INFO = -10
        ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
            INFO = -12
        END IF
        !
        IF( INFO.EQ.0 ) THEN
            !
            !        Compute the workspace requirements
            !
            NB = MIN( NBMAX, ILAENV( 1, 'SORMQR', SIDE // TRANS, M, N, K, &
                -1 ) )
            LWKOPT = MAX( 1, NW )*NB + TSIZE
            WORK( 1 ) = LWKOPT
        END IF
        !
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SORMQR', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
            WORK( 1 ) = 1
            RETURN
        END IF
        !
        NBMIN = 2
        LDWORK = NW
        IF( NB.GT.1 .AND. NB.LT.K ) THEN
            IF( LWORK.LT.NW*NB+TSIZE ) THEN
                NB = (LWORK-TSIZE) / LDWORK
                NBMIN = MAX( 2, ILAENV( 2, 'SORMQR', SIDE // TRANS, M, N, K, &
                    -1 ) )
            END IF
        END IF
        !
        IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
            !
            !        Use unblocked code
            !
            CALL SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, &
                IINFO )
        ELSE
            !
            !        Use blocked code
            !
            IWT = 1 + NW*NB
            IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. &
                ( .NOT.LEFT .AND. NOTRAN ) ) THEN
                I1 = 1
                I2 = K
                I3 = NB
            ELSE
                I1 = ( ( K-1 ) / NB )*NB + 1
                I2 = 1
                I3 = -NB
            END IF
            !
            IF( LEFT ) THEN
                NI = N
                JC = 1
            ELSE
                MI = M
                IC = 1
            END IF
            !
            DO  I = I1, I2, I3
                IB = MIN( NB, K-I+1 )
                !
                !           Form the triangular factor of the block reflector
                !           H = H(i) H(i+1) . . . H(i+ib-1)
                !
                CALL SLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), &
                    LDA, TAU( I ), WORK( IWT ), LDT )
                IF( LEFT ) THEN
                    !
                    !              H or H**T is applied to C(i:m,1:n)
                    !
                    MI = M - I + 1
                    IC = I
                ELSE
                    !
                    !              H or H**T is applied to C(1:m,i:n)
                    !
                    NI = N - I + 1
                    JC = I
                END IF
                !
                !           Apply H or H**T
                !
                CALL SLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
                    IB, A( I, I ), LDA, WORK( IWT ), LDT, &
                    C( IC, JC ), LDC, WORK, LDWORK )
            END DO
        END IF
        WORK( 1 ) = LWKOPT
        RETURN
        !
        !     End of SORMQR
        !
    END SUBROUTINE SORMQR

    SUBROUTINE SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
        WORK, INFO )
        !
        !  -- LAPACK computational routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        CHARACTER          SIDE, TRANS
        INTEGER            INFO, K, LDA, LDC, M, N
        !     ..
        !     .. Array Arguments ..
        REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE
        PARAMETER          ( ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            LEFT, NOTRAN
        INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
        REAL               AII
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input arguments
        !
        INFO = 0
        LEFT = LSAME( SIDE, 'L' )
        NOTRAN = LSAME( TRANS, 'N' )
        !
        !     NQ is the order of Q
        !
        IF( LEFT ) THEN
            NQ = M
        ELSE
            NQ = N
        END IF
        IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
            INFO = -2
        ELSE IF( M.LT.0 ) THEN
            INFO = -3
        ELSE IF( N.LT.0 ) THEN
            INFO = -4
        ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
            INFO = -5
        ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
            INFO = -7
        ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
            INFO = -10
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SORM2R', -INFO )
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) &
            RETURN
        !
        IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) &
            THEN
            I1 = 1
            I2 = K
            I3 = 1
        ELSE
            I1 = K
            I2 = 1
            I3 = -1
        END IF
        !
        IF( LEFT ) THEN
            NI = N
            JC = 1
        ELSE
            MI = M
            IC = 1
        END IF
        !
        DO  I = I1, I2, I3
            IF( LEFT ) THEN
                !
                !           H(i) is applied to C(i:m,1:n)
                !
                MI = M - I + 1
                IC = I
            ELSE
                !
                !           H(i) is applied to C(1:m,i:n)
                !
                NI = N - I + 1
                JC = I
            END IF
            !
            !        Apply H(i)
            !
            AII = A( I, I )
            A( I, I ) = ONE
            CALL SLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), &
                LDC, WORK )
            A( I, I ) = AII
        END DO
        RETURN
        !
        !     End of SORM2R
        !
    END SUBROUTINE SORM2R

    SUBROUTINE slaqr4( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
        ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
        !
        !  -- LAPACK auxiliary routine (version 3.7.0) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     December 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
        LOGICAL            WANTT, WANTZ
        !     ..
        !     .. Array Arguments ..
        REAL               H( ldh, * ), WI( * ), WORK( * ), WR( * ), &
            z( ldz, * )
        !     ..
        !
        !  ================================================================
        !
        !     .. Parameters ..
        !
        !     ==== Matrices of order NTINY or smaller must be processed by
        !     .    SLAHQR because of insufficient subdiagonal scratch space.
        !     .    (This is a hard limit.) ====
        INTEGER            NTINY
        parameter( ntiny = 11 )
        !
        !     ==== Exceptional deflation windows:  try to cure rare
        !     .    slow convergence by varying the size of the
        !     .    deflation window after KEXNW iterations. ====
        INTEGER            KEXNW
        parameter( kexnw = 5 )
        !
        !     ==== Exceptional shifts: try to cure rare slow convergence
        !     .    with ad-hoc exceptional shifts every KEXSH iterations.
        !     .    ====
        INTEGER            KEXSH
        parameter( kexsh = 6 )
        !
        !     ==== The constants WILK1 and WILK2 are used to form the
        !     .    exceptional shifts. ====
        REAL               WILK1, WILK2
        parameter( wilk1 = 0.75e0, wilk2 = -0.4375e0 )
        REAL               ZERO, ONE
        parameter( zero = 0.0e0, one = 1.0e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               AA, BB, CC, CS, DD, SN, SS, SWAP
        INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS, &
            kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls, &
            lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns, &
            nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
        LOGICAL            SORTED
        CHARACTER          JBCMPZ*2
        !     ..
        !     .. Local Arrays ..
        REAL               ZDUM( 1, 1 )
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          abs, int, max, min, mod, real
        !     ..
        !     .. Executable Statements ..
        info = 0
        !
        !     ==== Quick return for N = 0: nothing to do. ====
        !
        IF( n.EQ.0 ) THEN
            work( 1 ) = one
            RETURN
        END IF
        !
        IF( n.LE.ntiny ) THEN
            !
            !        ==== Tiny matrices must use SLAHQR. ====
            !
            lwkopt = 1
            IF( lwork.NE.-1 ) &
                CALL slahqr( wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, &
                iloz, ihiz, z, ldz, info )
        ELSE
            !
            !        ==== Use small bulge multi-shift QR with aggressive early
            !        .    deflation on larger-than-tiny matrices. ====
            !
            !        ==== Hope for the best. ====
            !
            info = 0
            !
            !        ==== Set up job flags for ILAENV. ====
            !
            IF( wantt ) THEN
                jbcmpz( 1: 1 ) = 'S'
            ELSE
                jbcmpz( 1: 1 ) = 'E'
            END IF
            IF( wantz ) THEN
                jbcmpz( 2: 2 ) = 'V'
            ELSE
                jbcmpz( 2: 2 ) = 'N'
            END IF
            !
            !        ==== NWR = recommended deflation window size.  At this
            !        .    point,  N .GT. NTINY = 11, so there is enough
            !        .    subdiagonal workspace for NWR.GE.2 as required.
            !        .    (In fact, there is enough subdiagonal space for
            !        .    NWR.GE.3.) ====
            !
            nwr = ilaenv( 13, 'SLAQR4', jbcmpz, n, ilo, ihi, lwork )
            nwr = max( 2, nwr )
            nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
            !
            !        ==== NSR = recommended number of simultaneous shifts.
            !        .    At this point N .GT. NTINY = 11, so there is at
            !        .    enough subdiagonal workspace for NSR to be even
            !        .    and greater than or equal to two as required. ====
            !
            nsr = ilaenv( 15, 'SLAQR4', jbcmpz, n, ilo, ihi, lwork )
            nsr = min( nsr, ( n+6 ) / 9, ihi-ilo )
            nsr = max( 2, nsr-mod( nsr, 2 ) )
            !
            !        ==== Estimate optimal workspace ====
            !
            !        ==== Workspace query call to SLAQR2 ====
            !
            CALL slaqr2( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz, &
                ihiz, z, ldz, ls, ld, wr, wi, h, ldh, n, h, ldh, &
                n, h, ldh, work, -1 )
            !
            !        ==== Optimal workspace = MAX(SLAQR5, SLAQR2) ====
            !
            lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
            !
            !        ==== Quick return in case of workspace query. ====
            !
            IF( lwork.EQ.-1 ) THEN
                work( 1 ) = REAL( lwkopt )
                RETURN
            END IF
            !
            !        ==== SLAHQR/SLAQR0 crossover point ====
            !
            nmin = ilaenv( 12, 'SLAQR4', jbcmpz, n, ilo, ihi, lwork )
            nmin = max( ntiny, nmin )
            !
            !        ==== Nibble crossover point ====
            !
            nibble = ilaenv( 14, 'SLAQR4', jbcmpz, n, ilo, ihi, lwork )
            nibble = max( 0, nibble )
            !
            !        ==== Accumulate reflections during ttswp?  Use block
            !        .    2-by-2 structure during matrix-matrix multiply? ====
            !
            kacc22 = ilaenv( 16, 'SLAQR4', jbcmpz, n, ilo, ihi, lwork )
            kacc22 = max( 0, kacc22 )
            kacc22 = min( 2, kacc22 )
            !
            !        ==== NWMAX = the largest possible deflation window for
            !        .    which there is sufficient workspace. ====
            !
            nwmax = min( ( n-1 ) / 3, lwork / 2 )
            nw = nwmax
            !
            !        ==== NSMAX = the Largest number of simultaneous shifts
            !        .    for which there is sufficient workspace. ====
            !
            nsmax = min( ( n+6 ) / 9, 2*lwork / 3 )
            nsmax = nsmax - mod( nsmax, 2 )
            !
            !        ==== NDFL: an iteration count restarted at deflation. ====
            !
            ndfl = 1
            !
            !        ==== ITMAX = iteration limit ====
            !
            itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
            !
            !        ==== Last row and column in the active block ====
            !
            kbot = ihi
            !
            !        ==== Main Loop ====
            !
            DO it = 1, itmax
                !
                !           ==== Done when KBOT falls below ILO ====
                !
                IF( kbot.LT.ilo ) &
                    GO TO 90
                !
                !           ==== Locate active block ====
                !
                DO k = kbot, ilo + 1, -1
                    IF( h( k, k-1 ).EQ.zero ) &
                        GO TO 20
                END DO
                k = ilo
20              CONTINUE
                ktop = k
                !
                !           ==== Select deflation window size:
                !           .    Typical Case:
                !           .      If possible and advisable, nibble the entire
                !           .      active block.  If not, use size MIN(NWR,NWMAX)
                !           .      or MIN(NWR+1,NWMAX) depending upon which has
                !           .      the smaller corresponding subdiagonal entry
                !           .      (a heuristic).
                !           .
                !           .    Exceptional Case:
                !           .      If there have been no deflations in KEXNW or
                !           .      more iterations, then vary the deflation window
                !           .      size.   At first, because, larger windows are,
                !           .      in general, more powerful than smaller ones,
                !           .      rapidly increase the window to the maximum possible.
                !           .      Then, gradually reduce the window size. ====
                !
                nh = kbot - ktop + 1
                nwupbd = min( nh, nwmax )
                IF( ndfl.LT.kexnw ) THEN
                    nw = min( nwupbd, nwr )
                ELSE
                    nw = min( nwupbd, 2*nw )
                END IF
                IF( nw.LT.nwmax ) THEN
                    IF( nw.GE.nh-1 ) THEN
                        nw = nh
                    ELSE
                        kwtop = kbot - nw + 1
                        IF( abs( h( kwtop, kwtop-1 ) ).GT. &
                            abs( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
                    END IF
                END IF
                IF( ndfl.LT.kexnw ) THEN
                    ndec = -1
                ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
                    ndec = ndec + 1
                    IF( nw-ndec.LT.2 ) &
                        ndec = 0
                    nw = nw - ndec
                END IF
                !
                !           ==== Aggressive early deflation:
                !           .    split workspace under the subdiagonal into
                !           .      - an nw-by-nw work array V in the lower
                !           .        left-hand-corner,
                !           .      - an NW-by-at-least-NW-but-more-is-better
                !           .        (NW-by-NHO) horizontal work array along
                !           .        the bottom edge,
                !           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
                !           .        vertical work array along the left-hand-edge.
                !           .        ====
                !
                kv = n - nw + 1
                kt = nw + 1
                nho = ( n-nw-1 ) - kt + 1
                kwv = nw + 2
                nve = ( n-nw ) - kwv + 1
                !
                !           ==== Aggressive early deflation ====
                !
                CALL slaqr2( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, &
                    ihiz, z, ldz, ls, ld, wr, wi, h( kv, 1 ), ldh, &
                    nho, h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh, &
                    work, lwork )
                !
                !           ==== Adjust KBOT accounting for new deflations. ====
                !
                kbot = kbot - ld
                !
                !           ==== KS points to the shifts. ====
                !
                ks = kbot - ls + 1
                !
                !           ==== Skip an expensive QR sweep if there is a (partly
                !           .    heuristic) reason to expect that many eigenvalues
                !           .    will deflate without it.  Here, the QR sweep is
                !           .    skipped if many eigenvalues have just been deflated
                !           .    or if the remaining active block is small.
                !
                IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot- &
                    ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
                    !
                    !              ==== NS = nominal number of simultaneous shifts.
                    !              .    This may be lowered (slightly) if SLAQR2
                    !              .    did not provide that many shifts. ====
                    !
                    ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
                    ns = ns - mod( ns, 2 )
                    !
                    !              ==== If there have been no deflations
                    !              .    in a multiple of KEXSH iterations,
                    !              .    then try exceptional shifts.
                    !              .    Otherwise use shifts provided by
                    !              .    SLAQR2 above or from the eigenvalues
                    !              .    of a trailing principal submatrix. ====
                    !
                    IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                        ks = kbot - ns + 1
                        DO i = kbot, max( ks+1, ktop+2 ), -2
                            ss = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
                            aa = wilk1*ss + h( i, i )
                            bb = ss
                            cc = wilk2*ss
                            dd = aa
                            CALL slanv2( aa, bb, cc, dd, wr( i-1 ), wi( i-1 ), &
                                wr( i ), wi( i ), cs, sn )
                        END DO
                        IF( ks.EQ.ktop ) THEN
                            wr( ks+1 ) = h( ks+1, ks+1 )
                            wi( ks+1 ) = zero
                            wr( ks ) = wr( ks+1 )
                            wi( ks ) = wi( ks+1 )
                        END IF
                    ELSE
                        !
                        !                 ==== Got NS/2 or fewer shifts? Use SLAHQR
                        !                 .    on a trailing principal submatrix to
                        !                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
                        !                 .    there is enough space below the subdiagonal
                        !                 .    to fit an NS-by-NS scratch array.) ====
                        !
                        IF( kbot-ks+1.LE.ns / 2 ) THEN
                            ks = kbot - ns + 1
                            kt = n - ns + 1
                            CALL slacpy( 'A', ns, ns, h( ks, ks ), ldh, &
                                h( kt, 1 ), ldh )
                            CALL slahqr( .false., .false., ns, 1, ns, &
                                h( kt, 1 ), ldh, wr( ks ), wi( ks ), &
                                1, 1, zdum, 1, inf )
                            ks = ks + inf
                            !
                            !                    ==== In case of a rare QR failure use
                            !                    .    eigenvalues of the trailing 2-by-2
                            !                    .    principal submatrix.  ====
                            !
                            IF( ks.GE.kbot ) THEN
                                aa = h( kbot-1, kbot-1 )
                                cc = h( kbot, kbot-1 )
                                bb = h( kbot-1, kbot )
                                dd = h( kbot, kbot )
                                CALL slanv2( aa, bb, cc, dd, wr( kbot-1 ), &
                                    wi( kbot-1 ), wr( kbot ), &
                                    wi( kbot ), cs, sn )
                                ks = kbot - 1
                            END IF
                        END IF
                        !
                        IF( kbot-ks+1.GT.ns ) THEN
                            !
                            !                    ==== Sort the shifts (Helps a little)
                            !                    .    Bubble sort keeps complex conjugate
                            !                    .    pairs together. ====
                            !
                            sorted = .false.
                            DO k = kbot, ks + 1, -1
                                IF( sorted ) &
                                    GO TO 60
                                sorted = .true.
                                DO i = ks, k - 1
                                    IF( abs( wr( i ) )+abs( wi( i ) ).LT. &
                                        abs( wr( i+1 ) )+abs( wi( i+1 ) ) ) THEN
                                        sorted = .false.
                                        !
                                        swap = wr( i )
                                        wr( i ) = wr( i+1 )
                                        wr( i+1 ) = swap
                                        !
                                        swap = wi( i )
                                        wi( i ) = wi( i+1 )
                                        wi( i+1 ) = swap
                                    END IF
                                END DO
                            END DO
60                          CONTINUE
                        END IF
                        !
                        !                 ==== Shuffle shifts into pairs of real shifts
                        !                 .    and pairs of complex conjugate shifts
                        !                 .    assuming complex conjugate shifts are
                        !                 .    already adjacent to one another. (Yes,
                        !                 .    they are.)  ====
                        !
                        DO i = kbot, ks + 2, -2
                            IF( wi( i ).NE.-wi( i-1 ) ) THEN
                                !
                                swap = wr( i )
                                wr( i ) = wr( i-1 )
                                wr( i-1 ) = wr( i-2 )
                                wr( i-2 ) = swap
                                !
                                swap = wi( i )
                                wi( i ) = wi( i-1 )
                                wi( i-1 ) = wi( i-2 )
                                wi( i-2 ) = swap
                            END IF
                        END DO
                    END IF
                    !
                    !              ==== If there are only two shifts and both are
                    !              .    real, then use only one.  ====
                    !
                    IF( kbot-ks+1.EQ.2 ) THEN
                        IF( wi( kbot ).EQ.zero ) THEN
                            IF( abs( wr( kbot )-h( kbot, kbot ) ).LT. &
                                abs( wr( kbot-1 )-h( kbot, kbot ) ) ) THEN
                                wr( kbot-1 ) = wr( kbot )
                            ELSE
                                wr( kbot ) = wr( kbot-1 )
                            END IF
                        END IF
                    END IF
                    !
                    !              ==== Use up to NS of the the smallest magnatiude
                    !              .    shifts.  If there aren't NS shifts available,
                    !              .    then use them all, possibly dropping one to
                    !              .    make the number of shifts even. ====
                    !
                    ns = min( ns, kbot-ks+1 )
                    ns = ns - mod( ns, 2 )
                    ks = kbot - ns + 1
                    !
                    !              ==== Small-bulge multi-shift QR sweep:
                    !              .    split workspace under the subdiagonal into
                    !              .    - a KDU-by-KDU work array U in the lower
                    !              .      left-hand-corner,
                    !              .    - a KDU-by-at-least-KDU-but-more-is-better
                    !              .      (KDU-by-NHo) horizontal work array WH along
                    !              .      the bottom edge,
                    !              .    - and an at-least-KDU-but-more-is-better-by-KDU
                    !              .      (NVE-by-KDU) vertical work WV arrow along
                    !              .      the left-hand-edge. ====
                    !
                    kdu = 3*ns - 3
                    ku = n - kdu + 1
                    kwh = kdu + 1
                    nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
                    kwv = kdu + 4
                    nve = n - kdu - kwv + 1
                    !
                    !              ==== Small-bulge multi-shift QR sweep ====
                    !
                    CALL slaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns, &
                        wr( ks ), wi( ks ), h, ldh, iloz, ihiz, z, &
                        ldz, work, 3, h( ku, 1 ), ldh, nve, &
                        h( kwv, 1 ), ldh, nho, h( ku, kwh ), ldh )
                END IF
                !
                !           ==== Note progress (or the lack of it). ====
                !
                IF( ld.GT.0 ) THEN
                    ndfl = 1
                ELSE
                    ndfl = ndfl + 1
                END IF
                !
                !           ==== End of main loop ====
            END DO
            !
            !        ==== Iteration limit exceeded.  Set INFO to show where
            !        .    the problem occurred and exit. ====
            !
            info = kbot
90          CONTINUE
        END IF
        !
        !     ==== Return the optimal value of LWORK. ====
        !
        work( 1 ) = REAL( lwkopt )
        !
        !     ==== End of SLAQR4 ====
        !
    END SUBROUTINE slaqr4

    SUBROUTINE slaqr2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ, &
        IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T, &
        LDT, NV, WV, LDWV, WORK, LWORK )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2017
        !
        !     .. Scalar Arguments ..
        INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV, &
            ldz, lwork, n, nd, nh, ns, nv, nw
        LOGICAL            WANTT, WANTZ
        !     ..
        !     .. Array Arguments ..
        REAL               H( ldh, * ), SI( * ), SR( * ), T( ldt, * ), &
            v( ldv, * ), work( * ), wv( ldwv, * ), &
            z( ldz, * )
        !     ..
        !
        !  ================================================================
        !     .. Parameters ..
        REAL               ZERO, ONE
        parameter( zero = 0.0e0, one = 1.0e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               AA, BB, BETA, CC, CS, DD, EVI, EVK, FOO, S, &
            safmax, safmin, smlnum, sn, tau, ulp
        INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, K, KCOL, &
            kend, kln, krow, kwtop, ltop, lwk1, lwk2, &
            lwkopt
        LOGICAL            BULGE, SORTED
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          abs, int, max, min, REAL, SQRT
        !     ..
        !     .. Executable Statements ..
        !
        !     ==== Estimate optimal workspace. ====
        !
        jw = min( nw, kbot-ktop+1 )
        IF( jw.LE.2 ) THEN
            lwkopt = 1
        ELSE
            !
            !        ==== Workspace query call to SGEHRD ====
            !
            CALL sgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
            lwk1 = int( work( 1 ) )
            !
            !        ==== Workspace query call to SORMHR ====
            !
            CALL sormhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, ldv, &
                work, -1, info )
            lwk2 = int( work( 1 ) )
            !
            !        ==== Optimal workspace ====
            !
            lwkopt = jw + max( lwk1, lwk2 )
        END IF
        !
        !     ==== Quick return in case of workspace query. ====
        !
        IF( lwork.EQ.-1 ) THEN
            work( 1 ) = REAL( lwkopt )
            RETURN
        END IF
        !
        !     ==== Nothing to do ...
        !     ... for an empty active block ... ====
        ns = 0
        nd = 0
        work( 1 ) = one
        IF( ktop.GT.kbot ) &
            RETURN
        !     ... nor for an empty deflation window. ====
        IF( nw.LT.1 ) &
            RETURN
        !
        !     ==== Machine constants ====
        !
        safmin = slamch( 'SAFE MINIMUM' )
        safmax = one / safmin
        CALL slabad( safmin, safmax )
        ulp = slamch( 'PRECISION' )
        smlnum = safmin*( REAL( N ) / ULP )
        !
        !     ==== Setup deflation window ====
        !
        jw = min( nw, kbot-ktop+1 )
        kwtop = kbot - jw + 1
        IF( kwtop.EQ.ktop ) THEN
            s = zero
        ELSE
            s = h( kwtop, kwtop-1 )
        END IF
        !
        IF( kbot.EQ.kwtop ) THEN
            !
            !        ==== 1-by-1 deflation window: not much to do ====
            !
            sr( kwtop ) = h( kwtop, kwtop )
            si( kwtop ) = zero
            ns = 1
            nd = 0
            IF( abs( s ).LE.max( smlnum, ulp*abs( h( kwtop, kwtop ) ) ) ) &
                THEN
                ns = 0
                nd = 1
                IF( kwtop.GT.ktop ) &
                    h( kwtop, kwtop-1 ) = zero
            END IF
            work( 1 ) = one
            RETURN
        END IF
        !
        !     ==== Convert to spike-triangular form.  (In case of a
        !     .    rare QR failure, this routine continues to do
        !     .    aggressive early deflation using that part of
        !     .    the deflation window that converged using INFQR
        !     .    here and there to keep track.) ====
        !
        CALL slacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
        CALL scopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), ldt+1 )
        !
        CALL slaset( 'A', jw, jw, zero, one, v, ldv )
        CALL slahqr( .true., .true., jw, 1, jw, t, ldt, sr( kwtop ), &
            si( kwtop ), 1, jw, v, ldv, infqr )
        !
        !     ==== STREXC needs a clean margin near the diagonal ====
        !
        DO j = 1, jw - 3
            t( j+2, j ) = zero
            t( j+3, j ) = zero
        END DO
        IF( jw.GT.2 ) &
            t( jw, jw-2 ) = zero
        !
        !     ==== Deflation detection loop ====
        !
        ns = jw
        ilst = infqr + 1
20      CONTINUE
        IF( ilst.LE.ns ) THEN
            IF( ns.EQ.1 ) THEN
                bulge = .false.
            ELSE
                bulge = t( ns, ns-1 ).NE.zero
            END IF
            !
            !        ==== Small spike tip test for deflation ====
            !
            IF( .NOT.bulge ) THEN
                !
                !           ==== Real eigenvalue ====
                !
                foo = abs( t( ns, ns ) )
                IF( foo.EQ.zero ) &
                    foo = abs( s )
                IF( abs( s*v( 1, ns ) ).LE.max( smlnum, ulp*foo ) ) THEN
                    !
                    !              ==== Deflatable ====
                    !
                    ns = ns - 1
                ELSE
                    !
                    !              ==== Undeflatable.   Move it up out of the way.
                    !              .    (STREXC can not fail in this case.) ====
                    !
                    ifst = ns
                    CALL strexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work, &
                        info )
                    ilst = ilst + 1
                END IF
            ELSE
                !
                !           ==== Complex conjugate pair ====
                !
                foo = abs( t( ns, ns ) ) + sqrt( abs( t( ns, ns-1 ) ) )* &
                    sqrt( abs( t( ns-1, ns ) ) )
                IF( foo.EQ.zero ) &
                    foo = abs( s )
                IF( max( abs( s*v( 1, ns ) ), abs( s*v( 1, ns-1 ) ) ).LE. &
                    max( smlnum, ulp*foo ) ) THEN
                    !
                    !              ==== Deflatable ====
                    !
                    ns = ns - 2
                ELSE
                    !
                    !              ==== Undeflatable. Move them up out of the way.
                    !              .    Fortunately, STREXC does the right thing with
                    !              .    ILST in case of a rare exchange failure. ====
                    !
                    ifst = ns
                    CALL strexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work, &
                        info )
                    ilst = ilst + 2
                END IF
            END IF
            !
            !        ==== End deflation detection loop ====
            !
            GO TO 20
        END IF
        !
        !        ==== Return to Hessenberg form ====
        !
        IF( ns.EQ.0 ) &
            s = zero
        !
        IF( ns.LT.jw ) THEN
            !
            !        ==== sorting diagonal blocks of T improves accuracy for
            !        .    graded matrices.  Bubble sort deals well with
            !        .    exchange failures. ====
            !
            sorted = .false.
            i = ns + 1
30          CONTINUE
            IF( sorted ) &
                GO TO 50
            sorted = .true.
            !
            kend = i - 1
            i = infqr + 1
            IF( i.EQ.ns ) THEN
                k = i + 1
            ELSE IF( t( i+1, i ).EQ.zero ) THEN
                k = i + 1
            ELSE
                k = i + 2
            END IF
40          CONTINUE
            IF( k.LE.kend ) THEN
                IF( k.EQ.i+1 ) THEN
                    evi = abs( t( i, i ) )
                ELSE
                    evi = abs( t( i, i ) ) + sqrt( abs( t( i+1, i ) ) )* &
                        sqrt( abs( t( i, i+1 ) ) )
                END IF
                !
                IF( k.EQ.kend ) THEN
                    evk = abs( t( k, k ) )
                ELSE IF( t( k+1, k ).EQ.zero ) THEN
                    evk = abs( t( k, k ) )
                ELSE
                    evk = abs( t( k, k ) ) + sqrt( abs( t( k+1, k ) ) )* &
                        sqrt( abs( t( k, k+1 ) ) )
                END IF
                !
                IF( evi.GE.evk ) THEN
                    i = k
                ELSE
                    sorted = .false.
                    ifst = i
                    ilst = k
                    CALL strexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, work, &
                        info )
                    IF( info.EQ.0 ) THEN
                        i = ilst
                    ELSE
                        i = k
                    END IF
                END IF
                IF( i.EQ.kend ) THEN
                    k = i + 1
                ELSE IF( t( i+1, i ).EQ.zero ) THEN
                    k = i + 1
                ELSE
                    k = i + 2
                END IF
                GO TO 40
            END IF
            GO TO 30
50          CONTINUE
        END IF
        !
        !     ==== Restore shift/eigenvalue array from T ====
        !
        i = jw
60      CONTINUE
        IF( i.GE.infqr+1 ) THEN
            IF( i.EQ.infqr+1 ) THEN
                sr( kwtop+i-1 ) = t( i, i )
                si( kwtop+i-1 ) = zero
                i = i - 1
            ELSE IF( t( i, i-1 ).EQ.zero ) THEN
                sr( kwtop+i-1 ) = t( i, i )
                si( kwtop+i-1 ) = zero
                i = i - 1
            ELSE
                aa = t( i-1, i-1 )
                cc = t( i, i-1 )
                bb = t( i-1, i )
                dd = t( i, i )
                CALL slanv2( aa, bb, cc, dd, sr( kwtop+i-2 ), &
                    si( kwtop+i-2 ), sr( kwtop+i-1 ), &
                    si( kwtop+i-1 ), cs, sn )
                i = i - 2
            END IF
            GO TO 60
        END IF
        !
        IF( ns.LT.jw .OR. s.EQ.zero ) THEN
            IF( ns.GT.1 .AND. s.NE.zero ) THEN
                !
                !           ==== Reflect spike back into lower triangle ====
                !
                CALL scopy( ns, v, ldv, work, 1 )
                beta = work( 1 )
                CALL slarfg( ns, beta, work( 2 ), 1, tau )
                work( 1 ) = one
                !
                CALL slaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), ldt )
                !
                CALL slarf( 'L', ns, jw, work, 1, tau, t, ldt, &
                    work( jw+1 ) )
                CALL slarf( 'R', ns, ns, work, 1, tau, t, ldt, &
                    work( jw+1 ) )
                CALL slarf( 'R', jw, ns, work, 1, tau, v, ldv, &
                    work( jw+1 ) )
                !
                CALL sgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ), &
                    lwork-jw, info )
            END IF
            !
            !        ==== Copy updated reduced window into place ====
            !
            IF( kwtop.GT.1 ) &
                h( kwtop, kwtop-1 ) = s*v( 1, 1 )
            CALL slacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
            CALL scopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ), &
                ldh+1 )
            !
            !        ==== Accumulate orthogonal matrix in order update
            !        .    H and Z, if requested.  ====
            !
            IF( ns.GT.1 .AND. s.NE.zero ) &
                CALL sormhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, ldv, &
                work( jw+1 ), lwork-jw, info )
            !
            !        ==== Update vertical slab in H ====
            !
            IF( wantt ) THEN
                ltop = 1
            ELSE
                ltop = ktop
            END IF
            DO krow = ltop, kwtop - 1, nv
                kln = min( nv, kwtop-krow )
                CALL sgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ), &
                    ldh, v, ldv, zero, wv, ldwv )
                CALL slacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), ldh )
            END DO
            !
            !        ==== Update horizontal slab in H ====
            !
            IF( wantt ) THEN
                DO kcol = kbot + 1, n, nh
                    kln = min( nh, n-kcol+1 )
                    CALL sgemm( 'C', 'N', jw, kln, jw, one, v, ldv, &
                        h( kwtop, kcol ), ldh, zero, t, ldt )
                    CALL slacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ), &
                        ldh )
                END DO
            END IF
            !
            !        ==== Update vertical slab in Z ====
            !
            IF( wantz ) THEN
                DO krow = iloz, ihiz, nv
                    kln = min( nv, ihiz-krow+1 )
                    CALL sgemm( 'N', 'N', kln, jw, jw, one, z( krow, kwtop ), &
                        ldz, v, ldv, zero, wv, ldwv )
                    CALL slacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ), &
                        ldz )
                END DO
            END IF
        END IF

    END SUBROUTINE slaqr2

    SUBROUTINE slaqr5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, &
        SR, SI, H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, &
        LDU, NV, WV, LDWV, NH, WH, LDWH )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2016
        !
        !     .. Scalar Arguments ..
        INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV, &
            ldwh, ldwv, ldz, n, nh, nshfts, nv
        LOGICAL            WANTT, WANTZ
        !     ..
        !     .. Array Arguments ..
        REAL               H( ldh, * ), SI( * ), SR( * ), U( ldu, * ), &
            v( ldv, * ), wh( ldwh, * ), wv( ldwv, * ), &
            z( ldz, * )
        !     ..
        !
        !  ================================================================
        !     .. Parameters ..
        REAL               ZERO, ONE
        parameter( zero = 0.0e0, one = 1.0e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               ALPHA, BETA, H11, H12, H21, H22, REFSUM, &
            safmax, safmin, scl, smlnum, swap, tst1, tst2, &
            ulp
        INTEGER            I, I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN, &
            jrow, jtop, k, k1, kdu, kms, knz, krcol, kzs, &
            m, m22, mbot, mend, mstart, mtop, nbmps, ndcol, &
            ns, nu
        LOGICAL            ACCUM, BLK22, BMP22
        !     ..
        !     .. Intrinsic Functions ..
        !
        INTRINSIC          abs, max, min, mod, real
        !     ..
        !     .. Local Arrays ..
        REAL               VT( 3 )
        !     ..
        !     .. Executable Statements ..
        !
        !     ==== If there are no shifts, then there is nothing to do. ====
        !
        IF( nshfts.LT.2 ) &
            RETURN
        !
        !     ==== If the active block is empty or 1-by-1, then there
        !     .    is nothing to do. ====
        !
        IF( ktop.GE.kbot ) &
            RETURN
        !
        !     ==== Shuffle shifts into pairs of real shifts and pairs
        !     .    of complex conjugate shifts assuming complex
        !     .    conjugate shifts are already adjacent to one
        !     .    another. ====
        !
        DO i = 1, nshfts - 2, 2
            IF( si( i ).NE.-si( i+1 ) ) THEN
                !
                swap = sr( i )
                sr( i ) = sr( i+1 )
                sr( i+1 ) = sr( i+2 )
                sr( i+2 ) = swap
                !
                swap = si( i )
                si( i ) = si( i+1 )
                si( i+1 ) = si( i+2 )
                si( i+2 ) = swap
            END IF
        END DO
        !
        !     ==== NSHFTS is supposed to be even, but if it is odd,
        !     .    then simply reduce it by one.  The shuffle above
        !     .    ensures that the dropped shift is real and that
        !     .    the remaining shifts are paired. ====
        !
        ns = nshfts - mod( nshfts, 2 )
        !
        !     ==== Machine constants for deflation ====
        !
        safmin = slamch( 'SAFE MINIMUM' )
        safmax = one / safmin
        CALL slabad( safmin, safmax )
        ulp = slamch( 'PRECISION' )
        smlnum = safmin*( REAL( N ) / ULP )
        !
        !     ==== Use accumulated reflections to update far-from-diagonal
        !     .    entries ? ====
        !
        accum = ( kacc22.EQ.1 ) .OR. ( kacc22.EQ.2 )
        !
        !     ==== If so, exploit the 2-by-2 block structure? ====
        !
        blk22 = ( ns.GT.2 ) .AND. ( kacc22.EQ.2 )
        !
        !     ==== clear trash ====
        !
        IF( ktop+2.LE.kbot ) &
            h( ktop+2, ktop ) = zero
        !
        !     ==== NBMPS = number of 2-shift bulges in the chain ====
        !
        nbmps = ns / 2
        !
        !     ==== KDU = width of slab ====
        !
        kdu = 6*nbmps - 3
        !
        !     ==== Create and chase chains of NBMPS bulges ====
        !
        DO incol = 3*( 1-nbmps ) + ktop - 1, kbot - 2, 3*nbmps - 2
            ndcol = incol + kdu
            IF( accum ) &
                CALL slaset( 'ALL', kdu, kdu, zero, one, u, ldu )
            !
            !        ==== Near-the-diagonal bulge chase.  The following loop
            !        .    performs the near-the-diagonal part of a small bulge
            !        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
            !        .    chunk extends from column INCOL to column NDCOL
            !        .    (including both column INCOL and column NDCOL). The
            !        .    following loop chases a 3*NBMPS column long chain of
            !        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
            !        .    may be less than KTOP and and NDCOL may be greater than
            !        .    KBOT indicating phantom columns from which to chase
            !        .    bulges before they are actually introduced or to which
            !        .    to chase bulges beyond column KBOT.)  ====
            !
            DO  krcol = incol, min( incol+3*nbmps-3, kbot-2 )
                !
                !           ==== Bulges number MTOP to MBOT are active double implicit
                !           .    shift bulges.  There may or may not also be small
                !           .    2-by-2 bulge, if there is room.  The inactive bulges
                !           .    (if any) must wait until the active bulges have moved
                !           .    down the diagonal to make room.  The phantom matrix
                !           .    paradigm described above helps keep track.  ====
                !
                mtop = max( 1, ( ( ktop-1 )-krcol+2 ) / 3+1 )
                mbot = min( nbmps, ( kbot-krcol ) / 3 )
                m22 = mbot + 1
                bmp22 = ( mbot.LT.nbmps ) .AND. ( krcol+3*( m22-1 ) ).EQ. &
                    ( kbot-2 )
                !
                !           ==== Generate reflections to chase the chain right
                !           .    one column.  (The minimum value of K is KTOP-1.) ====
                !
                DO m = mtop, mbot
                    k = krcol + 3*( m-1 )
                    IF( k.EQ.ktop-1 ) THEN
                        CALL slaqr1( 3, h( ktop, ktop ), ldh, sr( 2*m-1 ), &
                            si( 2*m-1 ), sr( 2*m ), si( 2*m ), &
                            v( 1, m ) )
                        alpha = v( 1, m )
                        CALL slarfg( 3, alpha, v( 2, m ), 1, v( 1, m ) )
                    ELSE
                        beta = h( k+1, k )
                        v( 2, m ) = h( k+2, k )
                        v( 3, m ) = h( k+3, k )
                        CALL slarfg( 3, beta, v( 2, m ), 1, v( 1, m ) )
                        !
                        !                 ==== A Bulge may collapse because of vigilant
                        !                 .    deflation or destructive underflow.  In the
                        !                 .    underflow case, try the two-small-subdiagonals
                        !                 .    trick to try to reinflate the bulge.  ====
                        !
                        IF( h( k+3, k ).NE.zero .OR. h( k+3, k+1 ).NE. &
                            zero .OR. h( k+3, k+2 ).EQ.zero ) THEN
                            !
                            !                    ==== Typical case: not collapsed (yet). ====
                            !
                            h( k+1, k ) = beta
                            h( k+2, k ) = zero
                            h( k+3, k ) = zero
                        ELSE
                            !
                            !                    ==== Atypical case: collapsed.  Attempt to
                            !                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
                            !                    .    If the fill resulting from the new
                            !                    .    reflector is too large, then abandon it.
                            !                    .    Otherwise, use the new one. ====
                            !
                            CALL slaqr1( 3, h( k+1, k+1 ), ldh, sr( 2*m-1 ), &
                                si( 2*m-1 ), sr( 2*m ), si( 2*m ), &
                                vt )
                            alpha = vt( 1 )
                            CALL slarfg( 3, alpha, vt( 2 ), 1, vt( 1 ) )
                            refsum = vt( 1 )*( h( k+1, k )+vt( 2 )* &
                                h( k+2, k ) )
                            !
                            IF( abs( h( k+2, k )-refsum*vt( 2 ) )+ &
                                abs( refsum*vt( 3 ) ).GT.ulp* &
                                ( abs( h( k, k ) )+abs( h( k+1, &
                                k+1 ) )+abs( h( k+2, k+2 ) ) ) ) THEN
                                !
                                !                       ==== Starting a new bulge here would
                                !                       .    create non-negligible fill.  Use
                                !                       .    the old one with trepidation. ====
                                !
                                h( k+1, k ) = beta
                                h( k+2, k ) = zero
                                h( k+3, k ) = zero
                            ELSE
                                !
                                !                       ==== Stating a new bulge here would
                                !                       .    create only negligible fill.
                                !                       .    Replace the old reflector with
                                !                       .    the new one. ====
                                !
                                h( k+1, k ) = h( k+1, k ) - refsum
                                h( k+2, k ) = zero
                                h( k+3, k ) = zero
                                v( 1, m ) = vt( 1 )
                                v( 2, m ) = vt( 2 )
                                v( 3, m ) = vt( 3 )
                            END IF
                        END IF
                    END IF
                END DO
                !
                !           ==== Generate a 2-by-2 reflection, if needed. ====
                !
                k = krcol + 3*( m22-1 )
                IF( bmp22 ) THEN
                    IF( k.EQ.ktop-1 ) THEN
                        CALL slaqr1( 2, h( k+1, k+1 ), ldh, sr( 2*m22-1 ), &
                            si( 2*m22-1 ), sr( 2*m22 ), si( 2*m22 ), &
                            v( 1, m22 ) )
                        beta = v( 1, m22 )
                        CALL slarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
                    ELSE
                        beta = h( k+1, k )
                        v( 2, m22 ) = h( k+2, k )
                        CALL slarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
                        h( k+1, k ) = beta
                        h( k+2, k ) = zero
                    END IF
                END IF
                !
                !           ==== Multiply H by reflections from the left ====
                !
                IF( accum ) THEN
                    jbot = min( ndcol, kbot )
                ELSE IF( wantt ) THEN
                    jbot = n
                ELSE
                    jbot = kbot
                END IF
                DO j = max( ktop, krcol ), jbot
                    mend = min( mbot, ( j-krcol+2 ) / 3 )
                    DO m = mtop, mend
                        k = krcol + 3*( m-1 )
                        refsum = v( 1, m )*( h( k+1, j )+v( 2, m )* &
                            h( k+2, j )+v( 3, m )*h( k+3, j ) )
                        h( k+1, j ) = h( k+1, j ) - refsum
                        h( k+2, j ) = h( k+2, j ) - refsum*v( 2, m )
                        h( k+3, j ) = h( k+3, j ) - refsum*v( 3, m )
                    END DO
                END DO
                IF( bmp22 ) THEN
                    k = krcol + 3*( m22-1 )
                    DO j = max( k+1, ktop ), jbot
                        refsum = v( 1, m22 )*( h( k+1, j )+v( 2, m22 )* &
                            h( k+2, j ) )
                        h( k+1, j ) = h( k+1, j ) - refsum
                        h( k+2, j ) = h( k+2, j ) - refsum*v( 2, m22 )
                    END DO
                END IF
                !
                !           ==== Multiply H by reflections from the right.
                !           .    Delay filling in the last row until the
                !           .    vigilant deflation check is complete. ====
                !
                IF( accum ) THEN
                    jtop = max( ktop, incol )
                ELSE IF( wantt ) THEN
                    jtop = 1
                ELSE
                    jtop = ktop
                END IF
                DO m = mtop, mbot
                    IF( v( 1, m ).NE.zero ) THEN
                        k = krcol + 3*( m-1 )
                        DO j = jtop, min( kbot, k+3 )
                            refsum = v( 1, m )*( h( j, k+1 )+v( 2, m )* &
                                h( j, k+2 )+v( 3, m )*h( j, k+3 ) )
                            h( j, k+1 ) = h( j, k+1 ) - refsum
                            h( j, k+2 ) = h( j, k+2 ) - refsum*v( 2, m )
                            h( j, k+3 ) = h( j, k+3 ) - refsum*v( 3, m )
                        END DO
                        !
                        IF( accum ) THEN
                            !
                            !                    ==== Accumulate U. (If necessary, update Z later
                            !                    .    with with an efficient matrix-matrix
                            !                    .    multiply.) ====
                            !
                            kms = k - incol
                            DO j = max( 1, ktop-incol ), kdu
                                refsum = v( 1, m )*( u( j, kms+1 )+v( 2, m )* &
                                    u( j, kms+2 )+v( 3, m )*u( j, kms+3 ) )
                                u( j, kms+1 ) = u( j, kms+1 ) - refsum
                                u( j, kms+2 ) = u( j, kms+2 ) - refsum*v( 2, m )
                                u( j, kms+3 ) = u( j, kms+3 ) - refsum*v( 3, m )
                            END DO
                        ELSE IF( wantz ) THEN
                            !
                            !                    ==== U is not accumulated, so update Z
                            !                    .    now by multiplying by reflections
                            !                    .    from the right. ====
                            !
                            DO j = iloz, ihiz
                                refsum = v( 1, m )*( z( j, k+1 )+v( 2, m )* &
                                    z( j, k+2 )+v( 3, m )*z( j, k+3 ) )
                                z( j, k+1 ) = z( j, k+1 ) - refsum
                                z( j, k+2 ) = z( j, k+2 ) - refsum*v( 2, m )
                                z( j, k+3 ) = z( j, k+3 ) - refsum*v( 3, m )
                            END DO
                        END IF
                    END IF
                END DO
                !
                !           ==== Special case: 2-by-2 reflection (if needed) ====
                !
                k = krcol + 3*( m22-1 )
                IF( bmp22 ) THEN
                    IF ( v( 1, m22 ).NE.zero ) THEN
                        DO j = jtop, min( kbot, k+3 )
                            refsum = v( 1, m22 )*( h( j, k+1 )+v( 2, m22 )* &
                                h( j, k+2 ) )
                            h( j, k+1 ) = h( j, k+1 ) - refsum
                            h( j, k+2 ) = h( j, k+2 ) - refsum*v( 2, m22 )
                        END DO
                        !
                        IF( accum ) THEN
                            kms = k - incol
                            DO j = max( 1, ktop-incol ), kdu
                                refsum = v( 1, m22 )*( u( j, kms+1 )+ &
                                    v( 2, m22 )*u( j, kms+2 ) )
                                u( j, kms+1 ) = u( j, kms+1 ) - refsum
                                u( j, kms+2 ) = u( j, kms+2 ) - refsum* &
                                    v( 2, m22 )
                            END DO
                        ELSE IF( wantz ) THEN
                            DO j = iloz, ihiz
                                refsum = v( 1, m22 )*( z( j, k+1 )+v( 2, m22 )* &
                                    z( j, k+2 ) )
                                z( j, k+1 ) = z( j, k+1 ) - refsum
                                z( j, k+2 ) = z( j, k+2 ) - refsum*v( 2, m22 )
                            END DO
                        END IF
                    END IF
                END IF
                !
                !           ==== Vigilant deflation check ====
                !
                mstart = mtop
                IF( krcol+3*( mstart-1 ).LT.ktop ) &
                    mstart = mstart + 1
                mend = mbot
                IF( bmp22 ) &
                    mend = mend + 1
                IF( krcol.EQ.kbot-2 ) &
                    mend = mend + 1
                DO m = mstart, mend
                    k = min( kbot-1, krcol+3*( m-1 ) )
                    !
                    !              ==== The following convergence test requires that
                    !              .    the tradition small-compared-to-nearby-diagonals
                    !              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
                    !              .    criteria both be satisfied.  The latter improves
                    !              .    accuracy in some examples. Falling back on an
                    !              .    alternate convergence criterion when TST1 or TST2
                    !              .    is zero (as done here) is traditional but probably
                    !              .    unnecessary. ====
                    !
                    IF( h( k+1, k ).NE.zero ) THEN
                        tst1 = abs( h( k, k ) ) + abs( h( k+1, k+1 ) )
                        IF( tst1.EQ.zero ) THEN
                            IF( k.GE.ktop+1 ) &
                                tst1 = tst1 + abs( h( k, k-1 ) )
                            IF( k.GE.ktop+2 ) &
                                tst1 = tst1 + abs( h( k, k-2 ) )
                            IF( k.GE.ktop+3 ) &
                                tst1 = tst1 + abs( h( k, k-3 ) )
                            IF( k.LE.kbot-2 ) &
                                tst1 = tst1 + abs( h( k+2, k+1 ) )
                            IF( k.LE.kbot-3 ) &
                                tst1 = tst1 + abs( h( k+3, k+1 ) )
                            IF( k.LE.kbot-4 ) &
                                tst1 = tst1 + abs( h( k+4, k+1 ) )
                        END IF
                        IF( abs( h( k+1, k ) ).LE.max( smlnum, ulp*tst1 ) ) &
                            THEN
                            h12 = max( abs( h( k+1, k ) ), abs( h( k, k+1 ) ) )
                            h21 = min( abs( h( k+1, k ) ), abs( h( k, k+1 ) ) )
                            h11 = max( abs( h( k+1, k+1 ) ), &
                                abs( h( k, k )-h( k+1, k+1 ) ) )
                            h22 = min( abs( h( k+1, k+1 ) ), &
                                abs( h( k, k )-h( k+1, k+1 ) ) )
                            scl = h11 + h12
                            tst2 = h22*( h11 / scl )
                            !
                            IF( tst2.EQ.zero .OR. h21*( h12 / scl ).LE. &
                                max( smlnum, ulp*tst2 ) )h( k+1, k ) = zero
                        END IF
                    END IF
                END DO
                !
                !           ==== Fill in the last row of each bulge. ====
                !
                mend = min( nbmps, ( kbot-krcol-1 ) / 3 )
                DO m = mtop, mend
                    k = krcol + 3*( m-1 )
                    refsum = v( 1, m )*v( 3, m )*h( k+4, k+3 )
                    h( k+4, k+1 ) = -refsum
                    h( k+4, k+2 ) = -refsum*v( 2, m )
                    h( k+4, k+3 ) = h( k+4, k+3 ) - refsum*v( 3, m )
                END DO
                !
                !           ==== End of near-the-diagonal bulge chase. ====
                !
            END DO
            !
            !        ==== Use U (if accumulated) to update far-from-diagonal
            !        .    entries in H.  If required, use U to update Z as
            !        .    well. ====
            !
            IF( accum ) THEN
                IF( wantt ) THEN
                    jtop = 1
                    jbot = n
                ELSE
                    jtop = ktop
                    jbot = kbot
                END IF
                IF( ( .NOT.blk22 ) .OR. ( incol.LT.ktop ) .OR. &
                    ( ndcol.GT.kbot ) .OR. ( ns.LE.2 ) ) THEN
                    !
                    !              ==== Updates not exploiting the 2-by-2 block
                    !              .    structure of U.  K1 and NU keep track of
                    !              .    the location and size of U in the special
                    !              .    cases of introducing bulges and chasing
                    !              .    bulges off the bottom.  In these special
                    !              .    cases and in case the number of shifts
                    !              .    is NS = 2, there is no 2-by-2 block
                    !              .    structure to exploit.  ====
                    !
                    k1 = max( 1, ktop-incol )
                    nu = ( kdu-max( 0, ndcol-kbot ) ) - k1 + 1
                    !
                    !              ==== Horizontal Multiply ====
                    !
                    DO jcol = min( ndcol, kbot ) + 1, jbot, nh
                        jlen = min( nh, jbot-jcol+1 )
                        CALL sgemm( 'C', 'N', nu, jlen, nu, one, u( k1, k1 ), &
                            ldu, h( incol+k1, jcol ), ldh, zero, wh, &
                            ldwh )
                        CALL slacpy( 'ALL', nu, jlen, wh, ldwh, &
                            h( incol+k1, jcol ), ldh )
                    END DO
                    !
                    !              ==== Vertical multiply ====
                    !
                    DO jrow = jtop, max( ktop, incol ) - 1, nv
                        jlen = min( nv, max( ktop, incol )-jrow )
                        CALL sgemm( 'N', 'N', jlen, nu, nu, one, &
                            h( jrow, incol+k1 ), ldh, u( k1, k1 ), &
                            ldu, zero, wv, ldwv )
                        CALL slacpy( 'ALL', jlen, nu, wv, ldwv, &
                            h( jrow, incol+k1 ), ldh )
                    END DO
                    !
                    !              ==== Z multiply (also vertical) ====
                    !
                    IF( wantz ) THEN
                        DO jrow = iloz, ihiz, nv
                            jlen = min( nv, ihiz-jrow+1 )
                            CALL sgemm( 'N', 'N', jlen, nu, nu, one, &
                                z( jrow, incol+k1 ), ldz, u( k1, k1 ), &
                                ldu, zero, wv, ldwv )
                            CALL slacpy( 'ALL', jlen, nu, wv, ldwv, &
                                z( jrow, incol+k1 ), ldz )
                        END DO
                    END IF
                ELSE
                    !
                    !              ==== Updates exploiting U's 2-by-2 block structure.
                    !              .    (I2, I4, J2, J4 are the last rows and columns
                    !              .    of the blocks.) ====
                    !
                    i2 = ( kdu+1 ) / 2
                    i4 = kdu
                    j2 = i4 - i2
                    j4 = kdu
                    !
                    !              ==== KZS and KNZ deal with the band of zeros
                    !              .    along the diagonal of one of the triangular
                    !              .    blocks. ====
                    !
                    kzs = ( j4-j2 ) - ( ns+1 )
                    knz = ns + 1
                    !
                    !              ==== Horizontal multiply ====
                    !
                    DO jcol = min( ndcol, kbot ) + 1, jbot, nh
                        jlen = min( nh, jbot-jcol+1 )
                        !
                        !                 ==== Copy bottom of H to top+KZS of scratch ====
                        !                  (The first KZS rows get multiplied by zero.) ====
                        !
                        CALL slacpy( 'ALL', knz, jlen, h( incol+1+j2, jcol ), &
                            ldh, wh( kzs+1, 1 ), ldwh )
                        !
                        !                 ==== Multiply by U21**T ====
                        !
                        CALL slaset( 'ALL', kzs, jlen, zero, zero, wh, ldwh )
                        CALL strmm( 'L', 'U', 'C', 'N', knz, jlen, one, &
                            u( j2+1, 1+kzs ), ldu, wh( kzs+1, 1 ), &
                            ldwh )
                        !
                        !                 ==== Multiply top of H by U11**T ====
                        !
                        CALL sgemm( 'C', 'N', i2, jlen, j2, one, u, ldu, &
                            h( incol+1, jcol ), ldh, one, wh, ldwh )
                        !
                        !                 ==== Copy top of H to bottom of WH ====
                        !
                        CALL slacpy( 'ALL', j2, jlen, h( incol+1, jcol ), ldh, &
                            wh( i2+1, 1 ), ldwh )
                        !
                        !                 ==== Multiply by U21**T ====
                        !
                        CALL strmm( 'L', 'L', 'C', 'N', j2, jlen, one, &
                            u( 1, i2+1 ), ldu, wh( i2+1, 1 ), ldwh )
                        !
                        !                 ==== Multiply by U22 ====
                        !
                        CALL sgemm( 'C', 'N', i4-i2, jlen, j4-j2, one, &
                            u( j2+1, i2+1 ), ldu, &
                            h( incol+1+j2, jcol ), ldh, one, &
                            wh( i2+1, 1 ), ldwh )
                        !
                        !                 ==== Copy it back ====
                        !
                        CALL slacpy( 'ALL', kdu, jlen, wh, ldwh, &
                            h( incol+1, jcol ), ldh )
                    END DO
                    !
                    !              ==== Vertical multiply ====
                    !
                    DO jrow = jtop, max( incol, ktop ) - 1, nv
                        jlen = min( nv, max( incol, ktop )-jrow )
                        !
                        !                 ==== Copy right of H to scratch (the first KZS
                        !                 .    columns get multiplied by zero) ====
                        !
                        CALL slacpy( 'ALL', jlen, knz, h( jrow, incol+1+j2 ), &
                            ldh, wv( 1, 1+kzs ), ldwv )
                        !
                        !                 ==== Multiply by U21 ====
                        !
                        CALL slaset( 'ALL', jlen, kzs, zero, zero, wv, ldwv )
                        CALL strmm( 'R', 'U', 'N', 'N', jlen, knz, one, &
                            u( j2+1, 1+kzs ), ldu, wv( 1, 1+kzs ), &
                            ldwv )
                        !
                        !                 ==== Multiply by U11 ====
                        !
                        CALL sgemm( 'N', 'N', jlen, i2, j2, one, &
                            h( jrow, incol+1 ), ldh, u, ldu, one, wv, &
                            ldwv )
                        !
                        !                 ==== Copy left of H to right of scratch ====
                        !
                        CALL slacpy( 'ALL', jlen, j2, h( jrow, incol+1 ), ldh, &
                            wv( 1, 1+i2 ), ldwv )
                        !
                        !                 ==== Multiply by U21 ====
                        !
                        CALL strmm( 'R', 'L', 'N', 'N', jlen, i4-i2, one, &
                            u( 1, i2+1 ), ldu, wv( 1, 1+i2 ), ldwv )
                        !
                        !                 ==== Multiply by U22 ====
                        !
                        CALL sgemm( 'N', 'N', jlen, i4-i2, j4-j2, one, &
                            h( jrow, incol+1+j2 ), ldh, &
                            u( j2+1, i2+1 ), ldu, one, wv( 1, 1+i2 ), &
                            ldwv )
                        !
                        !                 ==== Copy it back ====
                        !
                        CALL slacpy( 'ALL', jlen, kdu, wv, ldwv, &
                            h( jrow, incol+1 ), ldh )
                    END DO
                    !
                    !              ==== Multiply Z (also vertical) ====
                    !
                    IF( wantz ) THEN
                        DO jrow = iloz, ihiz, nv
                            jlen = min( nv, ihiz-jrow+1 )
                            !
                            !                    ==== Copy right of Z to left of scratch (first
                            !                    .     KZS columns get multiplied by zero) ====
                            !
                            CALL slacpy( 'ALL', jlen, knz, &
                                z( jrow, incol+1+j2 ), ldz, &
                                wv( 1, 1+kzs ), ldwv )
                            !
                            !                    ==== Multiply by U12 ====
                            !
                            CALL slaset( 'ALL', jlen, kzs, zero, zero, wv, &
                                ldwv )
                            CALL strmm( 'R', 'U', 'N', 'N', jlen, knz, one, &
                                u( j2+1, 1+kzs ), ldu, wv( 1, 1+kzs ), &
                                ldwv )
                            !
                            !                    ==== Multiply by U11 ====
                            !
                            CALL sgemm( 'N', 'N', jlen, i2, j2, one, &
                                z( jrow, incol+1 ), ldz, u, ldu, one, &
                                wv, ldwv )
                            !
                            !                    ==== Copy left of Z to right of scratch ====
                            !
                            CALL slacpy( 'ALL', jlen, j2, z( jrow, incol+1 ), &
                                ldz, wv( 1, 1+i2 ), ldwv )
                            !
                            !                    ==== Multiply by U21 ====
                            !
                            CALL strmm( 'R', 'L', 'N', 'N', jlen, i4-i2, one, &
                                u( 1, i2+1 ), ldu, wv( 1, 1+i2 ), &
                                ldwv )
                            !
                            !                    ==== Multiply by U22 ====
                            !
                            CALL sgemm( 'N', 'N', jlen, i4-i2, j4-j2, one, &
                                z( jrow, incol+1+j2 ), ldz, &
                                u( j2+1, i2+1 ), ldu, one, &
                                wv( 1, 1+i2 ), ldwv )
                            !
                            !                    ==== Copy the result back to Z ====
                            !
                            CALL slacpy( 'ALL', jlen, kdu, wv, ldwv, &
                                z( jrow, incol+1 ), ldz )
                        END DO
                    END IF
                END IF
            END IF
        END DO
        !
        !     ==== End of SLAQR5 ====
        !
    END SUBROUTINE slaqr5

    SUBROUTINE slaqr1( N, H, LDH, SR1, SI1, SR2, SI2, V )
        !
        !  -- LAPACK auxiliary routine (version 3.7.1) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     June 2017
        !
        !     .. Scalar Arguments ..
        REAL               SI1, SI2, SR1, SR2
        INTEGER            LDH, N
        !     ..
        !     .. Array Arguments ..
        REAL               H( ldh, * ), V( * )
        !     ..
        !
        !  ================================================================
        !
        !     .. Parameters ..
        REAL               ZERO
        parameter( zero = 0.0e0 )
        !     ..
        !     .. Local Scalars ..
        REAL               H21S, H31S, S
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          abs
        !     ..
        !     .. Executable Statements ..
        IF( n.EQ.2 ) THEN
            s = abs( h( 1, 1 )-sr2 ) + abs( si2 ) + abs( h( 2, 1 ) )
            IF( s.EQ.zero ) THEN
                v( 1 ) = zero
                v( 2 ) = zero
            ELSE
                h21s = h( 2, 1 ) / s
                v( 1 ) = h21s*h( 1, 2 ) + ( h( 1, 1 )-sr1 )* &
                    ( ( h( 1, 1 )-sr2 ) / s ) - si1*( si2 / s )
                v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-sr1-sr2 )
            END IF
        ELSE
            s = abs( h( 1, 1 )-sr2 ) + abs( si2 ) + abs( h( 2, 1 ) ) + &
                abs( h( 3, 1 ) )
            IF( s.EQ.zero ) THEN
                v( 1 ) = zero
                v( 2 ) = zero
                v( 3 ) = zero
            ELSE
                h21s = h( 2, 1 ) / s
                h31s = h( 3, 1 ) / s
                v( 1 ) = ( h( 1, 1 )-sr1 )*( ( h( 1, 1 )-sr2 ) / s ) - &
                    si1*( si2 / s ) + h( 1, 2 )*h21s + h( 1, 3 )*h31s
                v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-sr1-sr2 ) + &
                    h( 2, 3 )*h31s
                v( 3 ) = h31s*( h( 1, 1 )+h( 3, 3 )-sr1-sr2 ) + &
                    h21s*h( 3, 2 )
            END IF
        END IF
    END SUBROUTINE slaqr1

end module simple_lapackblas

! Copyright (c) 1992-2013 The University of Tennessee and The University
!                         of Tennessee Research Foundation.  All rights
!                         reserved.
! Copyright (c) 2000-2013 The University of California Berkeley. All
!                         rights reserved.
! Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
!                         reserved.
!
! $COPYRIGHT$
!
! Additional copyrights may follow
!
! $HEADER$
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! - Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
!
! - Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer listed
!   in this license in the documentation and/or other materials
!   provided with the distribution.
!
! - Neither the name of the copyright holders nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
! The copyright holders provide no reassurances that the source code
! provided does not infringe any patent, copyright, or any other
! intellectual property rights of third parties.  The copyright holders
! disclaim any liability to any recipient for claims brought against
! recipient by any third party for infringement of that parties
! intellectual property rights.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
