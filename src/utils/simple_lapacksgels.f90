! subroutines imported from sgels, lapack from netlib.org, manually made consistent with fortran 90
module simple_lapacksgels
implicit none
! public :: SGELS
public
!
! private

contains

  !  =====================================================================
        INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
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
        REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, NEGZRO, NEWZRO, POSINF
  !     ..
  !     .. Executable Statements ..
        IEEECK = 1
  !
        POSINF = ONE / ZERO
        IF( POSINF.LE.ONE ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        NEGINF = -ONE / ZERO
        IF( NEGINF.GE.ZERO ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        NEGZRO = ONE / ( NEGINF+ONE )
        IF( NEGZRO.NE.ZERO ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        NEGINF = ONE / NEGZRO
        IF( NEGINF.GE.ZERO ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        NEWZRO = NEGZRO + ZERO
        IF( NEWZRO.NE.ZERO ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        POSINF = ONE / NEWZRO
        IF( POSINF.LE.ONE ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        NEGINF = NEGINF*POSINF
        IF( NEGINF.GE.ZERO ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        POSINF = POSINF*POSINF
        IF( POSINF.LE.ONE ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
  !
  !
  !
  !     Return if we were only asked to check infinity arithmetic
  !
        IF( ISPEC.EQ.0 ) RETURN
  !
        NAN1 = POSINF + NEGINF
  !
        NAN2 = POSINF / NEGINF
  !
        NAN3 = POSINF / POSINF
  !
        NAN4 = POSINF*ZERO
  !
        NAN5 = NEGINF*NEGZRO
  !
        NAN6 = NAN5*ZERO
  !
        IF( NAN1.EQ.NAN1 ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        IF( NAN2.EQ.NAN2 ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        IF( NAN3.EQ.NAN3 ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        IF( NAN4.EQ.NAN4 ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        IF( NAN5.EQ.NAN5 ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        IF( NAN6.EQ.NAN6 ) THEN
           IEEECK = 0
           RETURN
        END IF
  !
        RETURN
        END FUNCTION IEEECK

    !  =====================================================================

        INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
    !
    !  -- LAPACK auxiliary routine (version 3.9.0) --
    !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    !     November 2019
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
        INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
    !     ..
    !     .. External Functions ..
        INTEGER            IEEECK, IPARMQ, IPARAM2STAGE
        EXTERNAL           IEEECK, IPARMQ, IPARAM2STAGE
    !     ..
    !     .. Executable Statements ..
    !
        GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, 130, 140, 150, 160, 160, 160, 160, 160)ISPEC
    !
    !     Invalid value for ISPEC
    !
        ILAENV = -1
        RETURN
    !
     10 CONTINUE
    !
    !     Convert NAME to upper case if the first character is lower case.
    !
        ILAENV = 1
        SUBNAM = NAME
        IC = ICHAR( SUBNAM( 1: 1 ) )
        IZ = ICHAR( 'Z' )
        IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
    !
    !        ASCII character set
    !
           IF( IC.GE.97 .AND. IC.LE.122 ) THEN
              SUBNAM( 1: 1 ) = CHAR( IC-32 )
              DO 20 I = 2, 6
                 IC = ICHAR( SUBNAM( I: I ) )
                 IF( IC.GE.97 .AND. IC.LE.122 ) SUBNAM( I: I ) = CHAR( IC-32 )
     20       CONTINUE
           END IF
    !
        ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
    !
    !        EBCDIC character set
    !
           IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.( IC.GE.145 .AND. IC.LE.153 ).OR.( IC.GE.162 .AND. IC.LE.169 ) ) THEN
              SUBNAM( 1: 1 ) = CHAR( IC+64 )
              DO 30 I = 2, 6
                 IC = ICHAR( SUBNAM( I: I ) )
                 IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.( IC.GE.145 .AND. IC.LE.153 ) .OR.( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
     30       CONTINUE
           END IF
    !
        ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
    !
    !        Prime machines:  ASCII+128
    !
           IF( IC.GE.225 .AND. IC.LE.250 ) THEN
              SUBNAM( 1: 1 ) = CHAR( IC-32 )
              DO 40 I = 2, 6
                 IC = ICHAR( SUBNAM( I: I ) )
                 IF( IC.GE.225 .AND. IC.LE.250 ) SUBNAM( I: I ) = CHAR( IC-32 )
     40       CONTINUE
           END IF
        END IF
    !
        C1 = SUBNAM( 1: 1 )
        SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
        CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
        IF( .NOT.( CNAME .OR. SNAME ) )RETURN
        C2 = SUBNAM( 2: 3 )
        C3 = SUBNAM( 4: 6 )
        C4 = C3( 2: 3 )
        TWOSTAGE = LEN( SUBNAM ).GE.11.AND. SUBNAM( 11: 11 ).EQ.'2'
    !
        GO TO ( 50, 60, 70 )ISPEC
    !
     50 CONTINUE
    !
    !     ISPEC = 1:  block size
    !
    !     In these examples, separate code is provided for setting NB for
    !     real and complex.  We assume that NB will take the same value in
    !     single or double precision.
    !
        NB = 1
    !
        IF( SUBNAM(2:6).EQ.'LAORH' ) THEN
    !
    !        This is for *LAORHR_GETRFNP routine
    !
           IF( SNAME ) THEN
               NB = 32
           ELSE
               NB = 32
           END IF
        ELSE IF( C2.EQ.'GE' ) THEN
           IF( C3.EQ.'TRF' ) THEN
              IF( SNAME ) THEN
                 NB = 64
              ELSE
                 NB = 64
              END IF
           ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
              IF( SNAME ) THEN
                 NB = 32
              ELSE
                 NB = 32
              END IF
           ELSE IF( C3.EQ.'QR ') THEN
              IF( N3 .EQ. 1) THEN
                 IF( SNAME ) THEN
    !     M*N
                    IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                       NB = N1
                    ELSE
                       NB = 32768/N2
                    END IF
                 ELSE
                    IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                       NB = N1
                    ELSE
                       NB = 32768/N2
                    END IF
                 END IF
              ELSE
                 IF( SNAME ) THEN
                    NB = 1
                 ELSE
                    NB = 1
                 END IF
              END IF
           ELSE IF( C3.EQ.'LQ ') THEN
              IF( N3 .EQ. 2) THEN
                 IF( SNAME ) THEN
    !     M*N
                    IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                       NB = N1
                    ELSE
                       NB = 32768/N2
                    END IF
                 ELSE
                    IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                       NB = N1
                    ELSE
                       NB = 32768/N2
                    END IF
                 END IF
              ELSE
                 IF( SNAME ) THEN
                    NB = 1
                 ELSE
                    NB = 1
                 END IF
              END IF
           ELSE IF( C3.EQ.'HRD' ) THEN
              IF( SNAME ) THEN
                 NB = 32
              ELSE
                 NB = 32
              END IF
           ELSE IF( C3.EQ.'BRD' ) THEN
              IF( SNAME ) THEN
                 NB = 32
              ELSE
                 NB = 32
              END IF
           ELSE IF( C3.EQ.'TRI' ) THEN
              IF( SNAME ) THEN
                 NB = 64
              ELSE
                 NB = 64
              END IF
           END IF
        ELSE IF( C2.EQ.'PO' ) THEN
           IF( C3.EQ.'TRF' ) THEN
              IF( SNAME ) THEN
                 NB = 64
              ELSE
                 NB = 64
              END IF
           END IF
        ELSE IF( C2.EQ.'SY' ) THEN
           IF( C3.EQ.'TRF' ) THEN
              IF( SNAME ) THEN
                 IF( TWOSTAGE ) THEN
                    NB = 192
                 ELSE
                    NB = 64
                 END IF
              ELSE
                 IF( TWOSTAGE ) THEN
                    NB = 192
                 ELSE
                    NB = 64
                 END IF
              END IF
           ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
              NB = 32
           ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
              NB = 64
           END IF
        ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
           IF( C3.EQ.'TRF' ) THEN
              IF( TWOSTAGE ) THEN
                 NB = 192
              ELSE
                 NB = 64
              END IF
           ELSE IF( C3.EQ.'TRD' ) THEN
              NB = 32
           ELSE IF( C3.EQ.'GST' ) THEN
              NB = 64
           END IF
        ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
           IF( C3( 1: 1 ).EQ.'G' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )THEN
                 NB = 32
              END IF
           ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )THEN
                 NB = 32
              END IF
           END IF
        ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
           IF( C3( 1: 1 ).EQ.'G' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )  THEN
                 NB = 32
              END IF
           ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                 NB = 32
              END IF
           END IF
        ELSE IF( C2.EQ.'GB' ) THEN
           IF( C3.EQ.'TRF' ) THEN
              IF( SNAME ) THEN
                 IF( N4.LE.64 ) THEN
                    NB = 1
                 ELSE
                    NB = 32
                 END IF
              ELSE
                 IF( N4.LE.64 ) THEN
                    NB = 1
                 ELSE
                    NB = 32
                 END IF
              END IF
           END IF
        ELSE IF( C2.EQ.'PB' ) THEN
           IF( C3.EQ.'TRF' ) THEN
              IF( SNAME ) THEN
                 IF( N2.LE.64 ) THEN
                    NB = 1
                 ELSE
                    NB = 32
                 END IF
              ELSE
                 IF( N2.LE.64 ) THEN
                    NB = 1
                 ELSE
                    NB = 32
                 END IF
              END IF
           END IF
        ELSE IF( C2.EQ.'TR' ) THEN
           IF( C3.EQ.'TRI' ) THEN
              IF( SNAME ) THEN
                 NB = 64
              ELSE
                 NB = 64
              END IF
           ELSE IF ( C3.EQ.'EVC' ) THEN
              IF( SNAME ) THEN
                 NB = 64
              ELSE
                 NB = 64
              END IF
           END IF
        ELSE IF( C2.EQ.'LA' ) THEN
           IF( C3.EQ.'UUM' ) THEN
              IF( SNAME ) THEN
                 NB = 64
              ELSE
                 NB = 64
              END IF
           END IF
        ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
           IF( C3.EQ.'EBZ' ) THEN
              NB = 1
           END IF
        ELSE IF( C2.EQ.'GG' ) THEN
           NB = 32
           IF( C3.EQ.'HD3' ) THEN
              IF( SNAME ) THEN
                 NB = 32
              ELSE
                 NB = 32
              END IF
           END IF
        END IF
        ILAENV = NB
        RETURN
    !
     60 CONTINUE
    !
    !     ISPEC = 2:  minimum block size
    !
        NBMIN = 2
        IF( C2.EQ.'GE' ) THEN
           IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
              IF( SNAME ) THEN
                 NBMIN = 2
              ELSE
                 NBMIN = 2
              END IF
           ELSE IF( C3.EQ.'HRD' ) THEN
              IF( SNAME ) THEN
                 NBMIN = 2
              ELSE
                 NBMIN = 2
              END IF
           ELSE IF( C3.EQ.'BRD' ) THEN
              IF( SNAME ) THEN
                 NBMIN = 2
              ELSE
                 NBMIN = 2
              END IF
           ELSE IF( C3.EQ.'TRI' ) THEN
              IF( SNAME ) THEN
                 NBMIN = 2
              ELSE
                 NBMIN = 2
              END IF
           END IF
        ELSE IF( C2.EQ.'SY' ) THEN
           IF( C3.EQ.'TRF' ) THEN
              IF( SNAME ) THEN
                 NBMIN = 8
              ELSE
                 NBMIN = 8
              END IF
           ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
              NBMIN = 2
           END IF
        ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
           IF( C3.EQ.'TRD' ) THEN
              NBMIN = 2
           END IF
        ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
           IF( C3( 1: 1 ).EQ.'G' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )THEN
                 NBMIN = 2
              END IF
           ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                 NBMIN = 2
              END IF
           END IF
        ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
           IF( C3( 1: 1 ).EQ.'G' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                 NBMIN = 2
              END IF
           ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )THEN
                 NBMIN = 2
              END IF
           END IF
        ELSE IF( C2.EQ.'GG' ) THEN
           NBMIN = 2
           IF( C3.EQ.'HD3' ) THEN
              NBMIN = 2
           END IF
        END IF
        ILAENV = NBMIN
        RETURN
    !
     70 CONTINUE
    !
    !     ISPEC = 3:  crossover point
    !
        NX = 0
        IF( C2.EQ.'GE' ) THEN
           IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.'QLF' ) THEN
              IF( SNAME ) THEN
                 NX = 128
              ELSE
                 NX = 128
              END IF
           ELSE IF( C3.EQ.'HRD' ) THEN
              IF( SNAME ) THEN
                 NX = 128
              ELSE
                 NX = 128
              END IF
           ELSE IF( C3.EQ.'BRD' ) THEN
              IF( SNAME ) THEN
                 NX = 128
              ELSE
                 NX = 128
              END IF
           END IF
        ELSE IF( C2.EQ.'SY' ) THEN
           IF( SNAME .AND. C3.EQ.'TRD' ) THEN
              NX = 32
           END IF
        ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
           IF( C3.EQ.'TRD' ) THEN
              NX = 32
           END IF
        ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
           IF( C3( 1: 1 ).EQ.'G' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )  THEN
                 NX = 128
              END IF
           END IF
        ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
           IF( C3( 1: 1 ).EQ.'G' ) THEN
              IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ. 'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' ) THEN
                 NX = 128
              END IF
           END IF
        ELSE IF( C2.EQ.'GG' ) THEN
           NX = 128
           IF( C3.EQ.'HD3' ) THEN
              NX = 128
           END IF
        END IF
        ILAENV = NX
        RETURN
    !
     80 CONTINUE
    !
    !     ISPEC = 4:  number of shifts (used by xHSEQR)
    !
        ILAENV = 6
        RETURN
    !
     90 CONTINUE
    !
    !     ISPEC = 5:  minimum column dimension (not used)
    !
        ILAENV = 2
        RETURN
    !
    100 CONTINUE
    !
    !     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
    !
        ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
        RETURN
    !
    110 CONTINUE
    !
    !     ISPEC = 7:  number of processors (not used)
    !
        ILAENV = 1
        RETURN
    !
    120 CONTINUE
    !
    !     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
    !
        ILAENV = 50
        RETURN
    !
    130 CONTINUE
    !
    !     ISPEC = 9:  maximum size of the subproblems at the bottom of the
    !                 computation tree in the divide-and-conquer algorithm
    !                 (used by xGELSD and xGESDD)
    !
        ILAENV = 25
        RETURN
    !
    140 CONTINUE
    !
    !     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
    !
    !     ILAENV = 0
        ILAENV = 1
        IF( ILAENV.EQ.1 ) THEN
           ILAENV = IEEECK( 1, 0.0, 1.0 )
        END IF
        RETURN
    !
    150 CONTINUE
    !
    !     ISPEC = 11: infinity arithmetic can be trusted not to trap
    !
    !     ILAENV = 0
        ILAENV = 1
        IF( ILAENV.EQ.1 ) THEN
           ILAENV = IEEECK( 0, 0.0, 1.0 )
        END IF
        RETURN
    !
    160 CONTINUE
    !
    !     12 <= ISPEC <= 16: xHSEQR or related subroutines.
    !
        ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
        RETURN
    !
    !     End of ILAENV
    !
        END FUNCTION ILAENV
!  =====================================================================
    INTEGER FUNCTION ILASLC( M, N, A, LDA )
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
    REAL               A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
    REAL             ZERO
    PARAMETER ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
    INTEGER I
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
    IF( N.EQ.0 ) THEN
       ILASLC = N
    ELSE IF( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
       ILASLC = N
    ELSE
!     Now scan each column from the end, returning with the first non-zero.
       DO ILASLC = N, 1, -1
          DO I = 1, M
             IF( A(I, ILASLC).NE.ZERO ) RETURN
          END DO
       END DO
    END IF
    RETURN
  END FUNCTION ILASLC

!  =====================================================================
      INTEGER FUNCTION ILASLR( M, N, A, LDA )
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
      REAL               A( LDA, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL             ZERO
      PARAMETER ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER I, J
!     ..
!     .. Executable Statements ..
!
!     Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILASLR = M
      ELSEIF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILASLR = M
      ELSE
!     Scan up each column tracking the last zero row seen.
         ILASLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILASLR = MAX( ILASLR, I )
         END DO
      END IF
      RETURN
    END FUNCTION ILASLR

  !  =====================================================================
        INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
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
        PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,ISHFTS = 15, IACC22 = 16 )
        INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
        PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500 )
        REAL               TWO
        PARAMETER          ( TWO = 2.0 )
  !     ..
  !     .. Local Scalars ..
        INTEGER            NH, NS
        INTEGER            I, IC, IZ
        CHARACTER          SUBNAM*6
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          LOG, MAX, MOD, NINT, REAL
  !     ..
  !     .. Executable Statements ..
        IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR. ( ISPEC.EQ.IACC22 ) ) THEN
  !
  !        ==== Set the number simultaneous shifts ====
  !
           NH = IHI - ILO + 1
           NS = 2
           IF( NH.GE.30 )NS = 4
           IF( NH.GE.60 ) NS = 10
           IF( NH.GE.150 )NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
           IF( NH.GE.590 )NS = 64
           IF( NH.GE.3000 )NS = 128
           IF( NH.GE.6000 )NS = 256
           NS = MAX( 2, NS-MOD( NS, 2 ) )
        END IF
  !
        IF( ISPEC.EQ.INMIN ) THEN
  !
  !
  !        ===== Matrices of order smaller than NMIN get sent
  !        .     to xLAHQR, the classic double shift algorithm.
  !        .     This must be at least 11. ====
  !
           IPARMQ = NMIN
  !
        ELSE IF( ISPEC.EQ.INIBL ) THEN
  !
  !        ==== INIBL: skip a multi-shift qr iteration and
  !        .    whenever aggressive early deflation finds
  !        .    at least (NIBBLE*(window size)/100) deflations. ====
  !
           IPARMQ = NIBBLE
  !
        ELSE IF( ISPEC.EQ.ISHFTS ) THEN
  !
  !        ==== NSHFTS: The number of simultaneous shifts =====
  !
           IPARMQ = NS
  !
        ELSE IF( ISPEC.EQ.INWIN ) THEN
  !
  !        ==== NW: deflation window size.  ====
  !
           IF( NH.LE.KNWSWP ) THEN
              IPARMQ = NS
           ELSE
              IPARMQ = 3*NS / 2
           END IF
  !
        ELSE IF( ISPEC.EQ.IACC22 ) THEN
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
           IPARMQ = 0
           SUBNAM = NAME
           IC = ICHAR( SUBNAM( 1: 1 ) )
           IZ = ICHAR( 'Z' )
           IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
  !
  !           ASCII character set
  !
              IF( IC.GE.97 .AND. IC.LE.122 ) THEN
                 SUBNAM( 1: 1 ) = CHAR( IC-32 )
                 DO I = 2, 6
                    IC = ICHAR( SUBNAM( I: I ) )
                    IF( IC.GE.97 .AND. IC.LE.122 )SUBNAM( I: I ) = CHAR( IC-32 )
                 END DO
              END IF
  !
           ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
  !
  !           EBCDIC character set
  !
              IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.( IC.GE.145 .AND. IC.LE.153 ) .OR.( IC.GE.162 .AND. IC.LE.169 ) ) THEN
                 SUBNAM( 1: 1 ) = CHAR( IC+64 )
                 DO I = 2, 6
                    IC = ICHAR( SUBNAM( I: I ) )
                    IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.( IC.GE.145 .AND. IC.LE.153 ) .OR.( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
                 END DO
              END IF
  !
           ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
  !
  !           Prime machines:  ASCII+128
  !
              IF( IC.GE.225 .AND. IC.LE.250 ) THEN
                 SUBNAM( 1: 1 ) = CHAR( IC-32 )
                 DO I = 2, 6
                    IC = ICHAR( SUBNAM( I: I ) )
                    IF( IC.GE.225 .AND. IC.LE.250 ) SUBNAM( I: I ) = CHAR( IC-32 )
                 END DO
              END IF
           END IF
  !
           IF( SUBNAM( 2:6 ).EQ.'GGHRD' .OR.  SUBNAM( 2:6 ).EQ.'GGHD3' ) THEN
              IPARMQ = 1
              IF( NH.GE.K22MIN ) IPARMQ = 2
           ELSE IF ( SUBNAM( 4:6 ).EQ.'EXC' ) THEN
              IF( NH.GE.KACMIN )  IPARMQ = 1
              IF( NH.GE.K22MIN ) IPARMQ = 2
           ELSE IF ( SUBNAM( 2:6 ).EQ.'HSEQR' .OR.SUBNAM( 2:5 ).EQ.'LAQR' ) THEN
              IF( NS.GE.KACMIN ) IPARMQ = 1
              IF( NS.GE.K22MIN ) IPARMQ = 2
           END IF
  !
        ELSE
  !        ===== invalid value of ispec =====
           IPARMQ = -1
  !
        END IF
  !
  !     ==== End of IPARMQ ====
  !
END FUNCTION IPARMQ

!  =====================================================================
      LOGICAL FUNCTION LSAME(CA,CB)
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
      INTRINSIC ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
!     ..
!
!     Test if the characters are equal
!
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR('Z')
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
!
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
          IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
          IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
!
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
          IF (INTA.GE.129 .AND. INTA.LE.137 .OR. INTA.GE.145 .AND. INTA.LE.153 .OR.INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
          IF (INTB.GE.129 .AND. INTB.LE.137 .OR.INTB.GE.145 .AND. INTB.LE.153 .OR. INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
!
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
          IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
          IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
!
!     RETURN
!
!     End of LSAME
!
END FUNCTION LSAME

!  =====================================================================
      SUBROUTINE SCOMBSSQ( V1, V2 )
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2018
!
!     .. Array Arguments ..
      REAL               V1( 2 ), V2( 2 )
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
      IF( V1( 1 ).GE.V2( 1 ) ) THEN
         IF( V1( 1 ).NE.ZERO ) THEN
            V1( 2 ) = V1( 2 ) + ( V2( 1 ) / V1( 1 ) )**2 * V2( 2 )
         END IF
      ELSE
         V1( 2 ) = V2( 2 ) + ( V1( 1 ) / V2( 1 ) )**2 * V1( 2 )
         V1( 1 ) = V2( 1 )
      END IF
      RETURN
!
!     End of SCOMBSSQ
!
END SUBROUTINE SCOMBSSQ

!  =====================================================================
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


  !  =====================================================================
        SUBROUTINE SGELQ2( M, N, A, LDA, TAU, WORK, INFO )
  !
  !  -- LAPACK computational routine (version 3.9.0) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     November 2019
  !
  !     .. Scalar Arguments ..
        INTEGER            INFO, LDA, M, N
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
        INTEGER            I, K
        REAL               AII
  !     ..
  !     .. External Subroutines ..
        EXTERNAL           SLARF, SLARFG, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input arguments
  !
        INFO = 0
        IF( M.LT.0 ) THEN
           INFO = -1
        ELSE IF( N.LT.0 ) THEN
           INFO = -2
        ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
           INFO = -4
        END IF
        IF( INFO.NE.0 ) THEN
           CALL XERBLA( 'SGELQ2', -INFO )
           RETURN
        END IF
  !
        K = MIN( M, N )
  !
        DO 10 I = 1, K
  !
  !        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
  !
           CALL SLARFG( N-I+1, A( I, I ), A( I, MIN( I+1, N ) ), LDA, TAU( I ) )
           IF( I.LT.M ) THEN
  !
  !           Apply H(i) to A(i+1:m,i:n) from the right
  !
              AII = A( I, I )
              A( I, I ) = ONE
              CALL SLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ), A( I+1, I ), LDA, WORK )
              A( I, I ) = AII
           END IF
     10 CONTINUE
        RETURN
  !
  !     End of SGELQ2
  !
END SUBROUTINE SGELQ2

!  =====================================================================
      SUBROUTINE SGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB, NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGELQ2, SLARFB, SLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
      LWKOPT = M*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGELQF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'SGELQF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SGELQF', ' ', M, N, -1,-1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!
!           Compute the LQ factorization of the current block
!           A(i:i+ib-1,i:n)
!
            CALL SGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK,IINFO )
            IF( I+IB.LE.M ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL SLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ),LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i+ib:m,i:n) from the right
!
               CALL SLARFB( 'Right', 'No transpose', 'Forward','Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ), LDA, WORK, LDWORK, A( I+IB, I ), LDA,WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K )CALL SGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of SGELQF
!
END SUBROUTINE SGELQF

!  =====================================================================
      SUBROUTINE SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, TPSD
      INTEGER            BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE
      REAL               ANRM, BIGNUM, BNRM, SMLNUM
!     ..
!     .. Local Arrays ..
      REAL               RWORK( 1 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      REAL               SLAMCH, SLANGE
      EXTERNAL           LSAME, ILAENV, SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGELQF, SGEQRF, SLABAD, SLASCL, SLASET, SORMLQ,SORMQR, STRTRS, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'T' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, M, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, MN + MAX( MN, NRHS ) ) .AND..NOT.LQUERY ) THEN
         INFO = -10
      END IF
!
!     Figure out optimal block size
!
      IF( INFO.EQ.0 .OR. INFO.EQ.-10 ) THEN
!
         TPSD = .TRUE.
         IF( LSAME( TRANS, 'N' ) )TPSD = .FALSE.
!
         IF( M.GE.N ) THEN
            NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LN', M, NRHS, N,-1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LT', M, NRHS, N,-1 ) )
            END IF
         ELSE
            NB = ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
            IF( TPSD ) THEN
               NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LT', N, NRHS, M,-1 ) )
            ELSE
               NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LN', N, NRHS, M, -1 ) )
            END IF
         END IF
!
         WSIZE = MAX( 1, MN + MAX( MN, NRHS )*NB )
         WORK( 1 ) = REAL( WSIZE )
!
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGELS ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( MIN( M, N, NRHS ).EQ.0 ) THEN
         CALL SLASET( 'Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         RETURN
      END IF
!
!     Get machine parameters
!
      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
!
!     Scale A, B if max element outside range [SMLNUM,BIGNUM]
!
      ANRM = SLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
!
!        Matrix all zero. Return zero solution.
!
         CALL SLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         GO TO 50
      END IF
!
      BROW = M
      IF( TPSD )BROW = N
      BNRM = SLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
!
!        Scale matrix norm up to SMLNUM
!
         CALL SLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB,INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
!
!        Scale matrix norm down to BIGNUM
!
         CALL SLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB,INFO )
         IBSCL = 2
      END IF
!
      IF( M.GE.N ) THEN
!
!        compute QR factorization of A
!
         CALL SGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN,INFO )
!
!        workspace at least N, optimally N*NB
!
         IF( .NOT.TPSD ) THEN
!
!           Least-Squares Problem min || A * X - B ||
!
!           B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)
!
            CALL SORMQR( 'Left', 'Transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
!
            CALL STRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS,A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
            SCLLEN = N
!
         ELSE
!
!           Underdetermined system of equations A**T * X = B
!
!           B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)
!
            CALL STRTRS( 'Upper', 'Transpose', 'Non-unit', N, NRHS,A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
!           B(N+1:M,1:NRHS) = ZERO
!
            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
!
!           B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
!
            CALL SORMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
            SCLLEN = M
!
         END IF
!
      ELSE
!
!        Compute LQ factorization of A
!
         CALL SGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )
!
!        workspace at least M, optimally M*NB.
!
         IF( .NOT.TPSD ) THEN
!
!           underdetermined system of equations A * X = B
!
!           B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
!
            CALL STRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS,A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
!           B(M+1:N,1:NRHS) = 0
!
            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
!
!           B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)
!
            CALL SORMLQ( 'Left', 'Transpose', N, NRHS, M, A, LDA,WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
            SCLLEN = N
!
         ELSE
!
!           overdetermined system min || A**T * X - B ||
!
!           B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
!
            CALL SORMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA,WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN,INFO )
!
!           workspace at least NRHS, optimally NRHS*NB
!
!           B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)
!
            CALL STRTRS( 'Lower', 'Transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO )
!
            IF( INFO.GT.0 ) THEN
               RETURN
            END IF
!
            SCLLEN = M
!
         END IF
!
      END IF
!
!     Undo scaling
!
      IF( IASCL.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB,INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB,INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL SLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL SLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB,INFO )
      END IF
!
   50 CONTINUE
      WORK( 1 ) = REAL( WSIZE )
!
      RETURN
!
!     End of SGELS
!
END SUBROUTINE SGELS

!  =====================================================================
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
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
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
      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. (.NOT.LSAME(TRANSA,'T'))) THEN
          INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND.(.NOT.LSAME(TRANSB,'T'))) THEN
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
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.(((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
!
!     And if  alpha.eq.zero.
!
      IF (ALPHA.EQ.ZERO) THEN
          IF (BETA.EQ.ZERO) THEN
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      C(I,J) = ZERO
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 J = 1,N
                  DO 30 I = 1,M
                      C(I,J) = BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
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
              DO 90 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 50 I = 1,M
                          C(I,J) = ZERO
   50                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 60 I = 1,M
                          C(I,J) = BETA*C(I,J)
   60                 CONTINUE
                  END IF
                  DO 80 L = 1,K
                      TEMP = ALPHA*B(L,J)
                      DO 70 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B + beta*C
!
              DO 120 J = 1,N
                  DO 110 I = 1,M
                      TEMP = ZERO
                      DO 100 L = 1,K
                          TEMP = TEMP + A(L,I)*B(L,J)
  100                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  110             CONTINUE
  120         CONTINUE
          END IF
      ELSE
          IF (NOTA) THEN
!
!           Form  C := alpha*A*B**T + beta*C
!
              DO 170 J = 1,N
                  IF (BETA.EQ.ZERO) THEN
                      DO 130 I = 1,M
                          C(I,J) = ZERO
  130                 CONTINUE
                  ELSE IF (BETA.NE.ONE) THEN
                      DO 140 I = 1,M
                          C(I,J) = BETA*C(I,J)
  140                 CONTINUE
                  END IF
                  DO 160 L = 1,K
                      TEMP = ALPHA*B(J,L)
                      DO 150 I = 1,M
                          C(I,J) = C(I,J) + TEMP*A(I,L)
  150                 CONTINUE
  160             CONTINUE
  170         CONTINUE
          ELSE
!
!           Form  C := alpha*A**T*B**T + beta*C
!
              DO 200 J = 1,N
                  DO 190 I = 1,M
                      TEMP = ZERO
                      DO 180 L = 1,K
                          TEMP = TEMP + A(L,I)*B(J,L)
  180                 CONTINUE
                      IF (BETA.EQ.ZERO) THEN
                          C(I,J) = ALPHA*TEMP
                      ELSE
                          C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                      END IF
  190             CONTINUE
  200         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of SGEMM .
!
END SUBROUTINE SGEMM

!  =====================================================================
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
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.LSAME(TRANS,'C')) THEN
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
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
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
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
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
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
!
!        Form  y := alpha*A**T*x + y.
!
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
!
      RETURN
!
!     End of SGEMV .
!
END SUBROUTINE SGEMV

!  =====================================================================
      SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
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
      INTEGER            I, K
      REAL               AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, SLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQR2', -INFO )
         RETURN
      END IF
!
      K = MIN( M, N )
!
      DO 10 I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
         CALL SLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,TAU( I ) )
         IF( I.LT.N ) THEN
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
            AII = A( I, I )
            A( I, I ) = ONE
            CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
!
!     End of SGEQR2
!
END SUBROUTINE SGEQR2

!  =====================================================================
      SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEQR2, SLARFB, SLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
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
         NX = MAX( 0, ILAENV( 3, 'SGEQRF', ' ', M, N, -1, -1 ) )
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
               NBMIN = MAX( 2, ILAENV( 2, 'SGEQRF', ' ', M, N, -1,-1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!
!           Compute the QR factorization of the current block
!           A(i:m,i:i+ib-1)
!
            CALL SGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,IINFO )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB,A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H**T to A(i:m,i+ib:n) from the left
!
               CALL SLARFB( 'Left', 'Transpose', 'Forward', 'Columnwise', M-I+1, N-I-IB+1, IB,A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K )CALL SGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of SGEQRF
!
END SUBROUTINE SGEQRF

!  =====================================================================
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
!     .. External Subroutines ..
      EXTERNAL XERBLA
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
          DO 20 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  DO 10 I = 1,M
                      A(I,J) = A(I,J) + X(I)*TEMP
   10             CONTINUE
              END IF
              JY = JY + INCY
   20     CONTINUE
      ELSE
          IF (INCX.GT.0) THEN
              KX = 1
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
          DO 40 J = 1,N
              IF (Y(JY).NE.ZERO) THEN
                  TEMP = ALPHA*Y(JY)
                  IX = KX
                  DO 30 I = 1,M
                      A(I,J) = A(I,J) + X(IX)*TEMP
                      IX = IX + INCX
   30             CONTINUE
              END IF
              JY = JY + INCY
   40     CONTINUE
      END IF
!
      RETURN
!
!     End of SGER  .
!
END SUBROUTINE SGER

!  =====================================================================
      LOGICAL FUNCTION SISNAN( SIN )
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
!  .. External Functions ..
      LOGICAL SLAISNAN
      EXTERNAL SLAISNAN
!  ..
!  .. Executable Statements ..
      SISNAN = SLAISNAN(SIN,SIN)
      RETURN
    END FUNCTION SISNAN

  !  =====================================================================
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
  !     If it looks like we are on a Cray, take the square root of
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

!  =====================================================================
    LOGICAL FUNCTION SLAISNAN( SIN1, SIN2 )
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
    SLAISNAN = (SIN1.NE.SIN2)
    RETURN
  END FUNCTION SLAISNAN

  !  =====================================================================
        REAL             FUNCTION SLAMCH( CMACH )
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
        PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     ..
  !     .. Local Scalars ..
        REAL               RND, EPS, SFMIN, SMALL, RMACH
  !     ..
  !     .. External Functions ..
        LOGICAL            LSAME
        EXTERNAL           LSAME
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          DIGITS, EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, RADIX, TINY
  !     ..
  !     .. Executable Statements ..
  !
  !
  !     Assume rounding, not chopping. Always.
  !
        RND = ONE
  !
        IF( ONE.EQ.RND ) THEN
           EPS = EPSILON(ZERO) * 0.5
        ELSE
           EPS = EPSILON(ZERO)
        END IF
  !
        IF( LSAME( CMACH, 'E' ) ) THEN
           RMACH = EPS
        ELSE IF( LSAME( CMACH, 'S' ) ) THEN
           SFMIN = TINY(ZERO)
           SMALL = ONE / HUGE(ZERO)
           IF( SMALL.GE.SFMIN ) THEN
  !
  !           Use SMALL plus a bit, to avoid the possibility of rounding
  !           causing overflow when computing  1/sfmin.
  !
              SFMIN = SMALL*( ONE+EPS )
           END IF
           RMACH = SFMIN
        ELSE IF( LSAME( CMACH, 'B' ) ) THEN
           RMACH = RADIX(ZERO)
        ELSE IF( LSAME( CMACH, 'P' ) ) THEN
           RMACH = EPS * RADIX(ZERO)
        ELSE IF( LSAME( CMACH, 'N' ) ) THEN
           RMACH = DIGITS(ZERO)
        ELSE IF( LSAME( CMACH, 'R' ) ) THEN
           RMACH = RND
        ELSE IF( LSAME( CMACH, 'M' ) ) THEN
           RMACH = MINEXPONENT(ZERO)
        ELSE IF( LSAME( CMACH, 'U' ) ) THEN
           RMACH = tiny(zero)
        ELSE IF( LSAME( CMACH, 'L' ) ) THEN
           RMACH = MAXEXPONENT(ZERO)
        ELSE IF( LSAME( CMACH, 'O' ) ) THEN
           RMACH = HUGE(ZERO)
        ELSE
           RMACH = ZERO
        END IF
  !
        SLAMCH = RMACH
        RETURN
  !
  !     End of SLAMCH
  !
END FUNCTION SLAMCH

  !> \brief \b SLAMC3
  !> \details
  !> \b Purpose:
  !> \verbatim
  !> SLAMC3  is intended to force  A  and  B  to be stored prior to doing
  !> the addition of  A  and  B ,  for use in situations where optimizers
  !> might hold one of these in a register.
  !> \endverbatim
  !> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
  !> \date December 2016
  !> \ingroup auxOTHERauxiliary
  !>
  !> \param[in] A
  !> \verbatim
  !> \endverbatim
  !>
  !> \param[in] B
  !> \verbatim
  !>          The values A and B.
  !> \endverbatim
  !>
  !
        REAL             FUNCTION SLAMC3( A, B )
  !
  !  -- LAPACK auxiliary routine (version 3.7.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  !     November 2010
  !
  !     .. Scalar Arguments ..
        REAL               A, B
  !     ..
  ! =====================================================================
  !
  !     .. Executable Statements ..
  !
        SLAMC3 = A + B
  !
        RETURN
  !
  !     End of SLAMC3
  !
END FUNCTION SLAMC3
  !
  !  =====================================================================
        REAL             FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )
  !
  !  -- LAPACK auxiliary routine (version 3.7.0) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     December 2016
  !
        IMPLICIT NONE
  !     .. Scalar Arguments ..
        CHARACTER          NORM
        INTEGER            LDA, M, N
  !     ..
  !     .. Array Arguments ..
        REAL               A( LDA, * ), WORK( * )
  !     ..
  !
  ! =====================================================================
  !
  !     .. Parameters ..
        REAL               ONE, ZERO
        PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
  !     ..
  !     .. Local Scalars ..
        INTEGER            I, J
        REAL               SUM, VALUE, TEMP
  !     ..
  !     .. Local Arrays ..
        REAL               SSQ( 2 ), COLSSQ( 2 )
  !     ..
  !     .. External Subroutines ..
        EXTERNAL           SLASSQ, SCOMBSSQ
  !     ..
  !     .. External Functions ..
        LOGICAL            LSAME, SISNAN
        EXTERNAL           LSAME, SISNAN
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MIN, SQRT
  !     ..
  !     .. Executable Statements ..
  !
        IF( MIN( M, N ).EQ.0 ) THEN
           VALUE = ZERO
        ELSE IF( LSAME( NORM, 'M' ) ) THEN
  !
  !        Find max(abs(A(i,j))).
  !
           VALUE = ZERO
           DO 20 J = 1, N
              DO 10 I = 1, M
                 TEMP = ABS( A( I, J ) )
                 IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
     10       CONTINUE
     20    CONTINUE
        ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
  !
  !        Find norm1(A).
  !
           VALUE = ZERO
           DO 40 J = 1, N
              SUM = ZERO
              DO 30 I = 1, M
                 SUM = SUM + ABS( A( I, J ) )
     30       CONTINUE
              IF( VALUE.LT.SUM .OR. SISNAN( SUM ) ) VALUE = SUM
     40    CONTINUE
        ELSE IF( LSAME( NORM, 'I' ) ) THEN
  !
  !        Find normI(A).
  !
           DO 50 I = 1, M
              WORK( I ) = ZERO
     50    CONTINUE
           DO 70 J = 1, N
              DO 60 I = 1, M
                 WORK( I ) = WORK( I ) + ABS( A( I, J ) )
     60       CONTINUE
     70    CONTINUE
           VALUE = ZERO
           DO 80 I = 1, M
              TEMP = WORK( I )
              IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
     80    CONTINUE
        ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
  !
  !        Find normF(A).
  !        SSQ(1) is scale
  !        SSQ(2) is sum-of-squares
  !        For better accuracy, sum each column separately.
  !
           SSQ( 1 ) = ZERO
           SSQ( 2 ) = ONE
           DO 90 J = 1, N
              COLSSQ( 1 ) = ZERO
              COLSSQ( 2 ) = ONE
              CALL SLASSQ( M, A( 1, J ), 1, COLSSQ( 1 ), COLSSQ( 2 ) )
              CALL SCOMBSSQ( SSQ, COLSSQ )
     90    CONTINUE
           VALUE = SSQ( 1 )*SQRT( SSQ( 2 ) )
        END IF
  !
        SLANGE = VALUE
        RETURN
  !
  !     End of SLANGE
  !
END FUNCTION SLANGE


!  =====================================================================
      REAL             FUNCTION SLAPY2( X, Y )
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
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      REAL               W, XABS, YABS, Z
      LOGICAL            X_IS_NAN, Y_IS_NAN
!     ..
!     .. External Functions ..
      LOGICAL            SISNAN
      EXTERNAL           SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
!     ..
!     .. Executable Statements ..
!
      X_IS_NAN = SISNAN( X )
      Y_IS_NAN = SISNAN( Y )
      IF ( X_IS_NAN ) SLAPY2 = X
      IF ( Y_IS_NAN ) SLAPY2 = Y
!
      IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
         XABS = ABS( X )
         YABS = ABS( Y )
         W = MAX( XABS, YABS )
         Z = MIN( XABS, YABS )
         IF( Z.EQ.ZERO ) THEN
            SLAPY2 = W
         ELSE
            SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
         END IF
      END IF
      RETURN
!
!     End of SLAPY2
!
END FUNCTION SLAPY2

!  =====================================================================
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
!     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILASLR, ILASLC
      EXTERNAL           LSAME, ILASLR, ILASLC
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
            CALL SGEMV( 'Transpose', LASTV, LASTC, ONE, C, LDC, V, INCV,ZERO, WORK, 1 )
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
            CALL SGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
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

!  =====================================================================
      SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,T, LDT, C, LDC, WORK, LDWORK )
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
      REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ),WORK( LDWORK, * )
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
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SCOPY, SGEMM, STRMM
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( M.LE.0 .OR. N.LE.0 )RETURN
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
               DO 10 J = 1, K
                  CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2**T * V2
!
                  CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K,ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**T
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2 * W**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K,-ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**T
!
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C1
!
               DO 40 J = 1, K
                  CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K,ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**T
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K,-ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
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
               DO 70 J = 1, K
                  CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1**T * V1
!
                  CALL SGEMM( 'Transpose', 'No transpose', N, K, M-K,ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V * W**T
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1 * W**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', M-K, N, K,-ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W**T
!
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * Htranspose  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
!
!              W := C2
!
               DO 100 J = 1, K
                  CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, K, N-K,ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V**T
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, N-K, K,-ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W
!
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
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
               DO 130 J = 1, K
                  CALL SCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K,ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C2**T * V2**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**T * W**T
!
               IF( M.GT.K ) THEN
!
!                 C2 := C2 - V2**T * W**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,V( 1, K+1 ), LDV, WORK, LDWORK, ONE,C( K+1, 1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W**T
!
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C1
!
               DO 160 J = 1, K
                  CALL SCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
!
!              W := W * V1**T
!
               CALL STRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K,ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C2 * V2**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K,ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C2 := C2 - W * V2
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K, -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,C( 1, K+1 ), LDC )
               END IF
!
!              W := W * V1
!
               CALL STRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
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
               DO 190 J = 1, K
                  CALL SCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
!
!                 W := W + C1**T * V1**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE,C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T**T  or  W * T
!
               CALL STRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,ONE, T, LDT, WORK, LDWORK )
!
!              C := C - V**T * W**T
!
               IF( M.GT.K ) THEN
!
!                 C1 := C1 - V1**T * W**T
!
                  CALL SGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE,V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 := C2 - W**T
!
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
!
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
!
!              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
!
!              W := C2
!
               DO 220 J = 1, K
                  CALL SCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
!
!              W := W * V2**T
!
               CALL STRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K,ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
!
!                 W := W + C1 * V1**T
!
                  CALL SGEMM( 'No transpose', 'Transpose', M, K, N-K,ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
!
!              W := W * T  or  W * T**T
!
               CALL STRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,ONE, T, LDT, WORK, LDWORK )
!
!              C := C - W * V
!
               IF( N.GT.K ) THEN
!
!                 C1 := C1 - W * V1
!
                  CALL SGEMM( 'No transpose', 'No transpose', M, N-K, K,-ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
!
!              W := W * V2
!
               CALL STRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 := C1 - W
!
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
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

!  =====================================================================
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
!     .. External Functions ..
      REAL               SLAMCH, SLAPY2, SNRM2
      EXTERNAL           SLAMCH, SLAPY2, SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           SSCAL
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
   10       CONTINUE
            KNT = KNT + 1
            CALL SSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( (ABS( BETA ).LT.SAFMIN) .AND. (KNT .LT. 20) ) GO TO 10
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
         DO 20 J = 1, KNT
            BETA = BETA*SAFMIN
 20      CONTINUE
         ALPHA = BETA
      END IF
!
      RETURN
!
!     End of SLARFG
!
END SUBROUTINE SLARFG

!  =====================================================================
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
!     .. External Subroutines ..
      EXTERNAL           SGEMV, STRMV
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 )RETURN
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
                  CALL SGEMV( 'Transpose', J-I, I-1, -TAU( I ),V( I+1, 1 ), LDV, V( I+1, I ), 1, ONE,T( 1, I ), 1 )
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
                  CALL SGEMV( 'No transpose', I-1, J-I, -TAU( I ),V( 1, I+1 ), LDV, V( I, I+1 ), LDV,ONE, T( 1, I ), 1 )
               END IF
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,LDT, T( 1, I ), 1 )
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
                     CALL SGEMV( 'Transpose', N-K+I-J, K-I, -TAU( I ), V( J, I+1 ), LDV, V( J, I ), 1, ONE,T( I+1, I ), 1 )
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
                     CALL SGEMV( 'No transpose', K-I, N-K+I-J, -TAU( I ), V( I+1, J ), LDV, V( I, J ), LDV, ONE, T( I+1, I ), 1 )
                  END IF
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  CALL STRMV( 'Lower', 'No transpose', 'Non-unit', K-I, T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
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

!  =====================================================================
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
!     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH, SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
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
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
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
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN
!
!     Get machine parameters
!
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
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
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE.EQ.1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE.EQ.2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE.EQ.3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE.EQ.4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE.EQ.5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE.EQ.6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE )GO TO 10
!
      RETURN
!
!     End of SLASCL
!
END SUBROUTINE SLASCL

!  =====================================================================
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
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
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
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
!
      RETURN
!
!     End of SLASET
!
END SUBROUTINE SLASET

!  =====================================================================
      SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )
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
      PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      REAL               ABSXI
!     ..
!     .. External Functions ..
      LOGICAL            SISNAN
      EXTERNAL           SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            ABSXI = ABS( X( IX ) )
            IF( ABSXI.GT.ZERO.OR.SISNAN( ABSXI ) ) THEN
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of SLASSQ
!
END SUBROUTINE SLASSQ

!  =====================================================================
      REAL FUNCTION SNRM2(N,X,INCX)
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
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      REAL ABSXI,NORM,SCALE,SSQ
      INTEGER IX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
!     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
!
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
!
      SNRM2 = NORM
      RETURN
!
!     End of SNRM2.
!
END FUNCTION SNRM2

!  =====================================================================
      SUBROUTINE SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,WORK, INFO )
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
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, XERBLA
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
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )THEN
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
      DO 10 I = I1, I2, I3
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
         CALL SLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of SORM2R
!
END SUBROUTINE SORM2R

!  =====================================================================
      SUBROUTINE SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )
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
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, XERBLA
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
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORML2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )RETURN
!
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) )THEN
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
      DO 10 I = I1, I2, I3
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
         CALL SLARF( SIDE, MI, NI, A( I, I ), LDA, TAU( I ), C( IC, JC ), LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of SORML2
!
END SUBROUTINE SORML2

!  =====================================================================
      SUBROUTINE SORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,WORK, LWORK, INFO )
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
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1,TSIZE = LDT*NBMAX )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK,LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORML2, XERBLA
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
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
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
         NB = MIN( NBMAX, ILAENV( 1, 'SORMLQ', SIDE // TRANS, M, N, K, -1 ) )
         LWKOPT = MAX( 1, NW )*NB + TSIZE
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SORMLQ', -INFO )
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
            NBMIN = MAX( 2, ILAENV( 2, 'SORMLQ', SIDE // TRANS, M, N, K, -1 ) )
         END IF
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL SORML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,IINFO )
      ELSE
!
!        Use blocked code
!
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. NOTRAN ) .OR.( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
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
         IF( NOTRAN ) THEN
            TRANST = 'T'
         ELSE
            TRANST = 'N'
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL SLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ), LDA, TAU( I ), WORK( IWT ), LDT )
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
            CALL SLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB,A( I, I ), LDA, WORK( IWT ), LDT,C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of SORMLQ
!
END SUBROUTINE SORMLQ

!  =====================================================================
      SUBROUTINE SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
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
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ),WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1,TSIZE = LDT*NBMAX )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK,LWKOPT, MI, NB, NBMIN, NI, NQ, NW
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARFB, SLARFT, SORM2R, XERBLA
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
         NB = MIN( NBMAX, ILAENV( 1, 'SORMQR', SIDE // TRANS, M, N, K,-1 ) )
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
            NBMIN = MAX( 2, ILAENV( 2, 'SORMQR', SIDE // TRANS, M, N, K,-1 ) )
         END IF
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,IINFO )
      ELSE
!
!        Use blocked code
!
         IWT = 1 + NW*NB
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) THEN
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
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL SLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),LDA, TAU( I ), WORK( IWT ), LDT )
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
            CALL SLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,IB, A( I, I ), LDA, WORK( IWT ), LDT,C( IC, JC ), LDC, WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of SORMQR
!
END SUBROUTINE SORMQR

!  =====================================================================
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


    !  =====================================================================
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
    !     .. External Functions ..
          LOGICAL LSAME
          EXTERNAL LSAME
    !     ..
    !     .. External Subroutines ..
          EXTERNAL XERBLA
    !     ..
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
          ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.(.NOT.LSAME(TRANSA,'T')) .AND. (.NOT.LSAME(TRANSA,'C'))) THEN
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
              DO 20 J = 1,N
                  DO 10 I = 1,M
                      B(I,J) = ZERO
       10         CONTINUE
       20     CONTINUE
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
                      DO 50 J = 1,N
                          DO 40 K = 1,M
                              IF (B(K,J).NE.ZERO) THEN
                                  TEMP = ALPHA*B(K,J)
                                  DO 30 I = 1,K - 1
                                      B(I,J) = B(I,J) + TEMP*A(I,K)
       30                         CONTINUE
                                  IF (NOUNIT) TEMP = TEMP*A(K,K)
                                  B(K,J) = TEMP
                              END IF
       40                 CONTINUE
       50             CONTINUE
                  ELSE
                      DO 80 J = 1,N
                          DO 70 K = M,1,-1
                              IF (B(K,J).NE.ZERO) THEN
                                  TEMP = ALPHA*B(K,J)
                                  B(K,J) = TEMP
                                  IF (NOUNIT) B(K,J) = B(K,J)*A(K,K)
                                  DO 60 I = K + 1,M
                                      B(I,J) = B(I,J) + TEMP*A(I,K)
       60                         CONTINUE
                              END IF
       70                 CONTINUE
       80             CONTINUE
                  END IF
              ELSE
    !
    !           Form  B := alpha*A**T*B.
    !
                  IF (UPPER) THEN
                      DO 110 J = 1,N
                          DO 100 I = M,1,-1
                              TEMP = B(I,J)
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 90 K = 1,I - 1
                                  TEMP = TEMP + A(K,I)*B(K,J)
       90                     CONTINUE
                              B(I,J) = ALPHA*TEMP
      100                 CONTINUE
      110             CONTINUE
                  ELSE
                      DO 140 J = 1,N
                          DO 130 I = 1,M
                              TEMP = B(I,J)
                              IF (NOUNIT) TEMP = TEMP*A(I,I)
                              DO 120 K = I + 1,M
                                  TEMP = TEMP + A(K,I)*B(K,J)
      120                     CONTINUE
                              B(I,J) = ALPHA*TEMP
      130                 CONTINUE
      140             CONTINUE
                  END IF
              END IF
          ELSE
              IF (LSAME(TRANSA,'N')) THEN
    !
    !           Form  B := alpha*B*A.
    !
                  IF (UPPER) THEN
                      DO 180 J = N,1,-1
                          TEMP = ALPHA
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 150 I = 1,M
                              B(I,J) = TEMP*B(I,J)
      150                 CONTINUE
                          DO 170 K = 1,J - 1
                              IF (A(K,J).NE.ZERO) THEN
                                  TEMP = ALPHA*A(K,J)
                                  DO 160 I = 1,M
                                      B(I,J) = B(I,J) + TEMP*B(I,K)
      160                         CONTINUE
                              END IF
      170                 CONTINUE
      180             CONTINUE
                  ELSE
                      DO 220 J = 1,N
                          TEMP = ALPHA
                          IF (NOUNIT) TEMP = TEMP*A(J,J)
                          DO 190 I = 1,M
                              B(I,J) = TEMP*B(I,J)
      190                 CONTINUE
                          DO 210 K = J + 1,N
                              IF (A(K,J).NE.ZERO) THEN
                                  TEMP = ALPHA*A(K,J)
                                  DO 200 I = 1,M
                                      B(I,J) = B(I,J) + TEMP*B(I,K)
      200                         CONTINUE
                              END IF
      210                 CONTINUE
      220             CONTINUE
                  END IF
              ELSE
    !
    !           Form  B := alpha*B*A**T.
    !
                  IF (UPPER) THEN
                      DO 260 K = 1,N
                          DO 240 J = 1,K - 1
                              IF (A(J,K).NE.ZERO) THEN
                                  TEMP = ALPHA*A(J,K)
                                  DO 230 I = 1,M
                                      B(I,J) = B(I,J) + TEMP*B(I,K)
      230                         CONTINUE
                              END IF
      240                 CONTINUE
                          TEMP = ALPHA
                          IF (NOUNIT) TEMP = TEMP*A(K,K)
                          IF (TEMP.NE.ONE) THEN
                              DO 250 I = 1,M
                                  B(I,K) = TEMP*B(I,K)
      250                     CONTINUE
                          END IF
      260             CONTINUE
                  ELSE
                      DO 300 K = N,1,-1
                          DO 280 J = K + 1,N
                              IF (A(J,K).NE.ZERO) THEN
                                  TEMP = ALPHA*A(J,K)
                                  DO 270 I = 1,M
                                      B(I,J) = B(I,J) + TEMP*B(I,K)
      270                         CONTINUE
                              END IF
      280                 CONTINUE
                          TEMP = ALPHA
                          IF (NOUNIT) TEMP = TEMP*A(K,K)
                          IF (TEMP.NE.ONE) THEN
                              DO 290 I = 1,M
                                  B(I,K) = TEMP*B(I,K)
      290                     CONTINUE
                          END IF
      300             CONTINUE
                  END IF
              END IF
          END IF
    !
          RETURN
    !
    !     End of STRMM .
    !
  END SUBROUTINE STRMM

  !  =====================================================================
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
  !     .. External Functions ..
        LOGICAL LSAME
        EXTERNAL LSAME
  !     ..
  !     .. External Subroutines ..
        EXTERNAL XERBLA
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
        ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.LSAME(TRANS,'C')) THEN
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
                    DO 20 J = 1,N
                        IF (X(J).NE.ZERO) THEN
                            TEMP = X(J)
                            DO 10 I = 1,J - 1
                                X(I) = X(I) + TEMP*A(I,J)
     10                     CONTINUE
                            IF (NOUNIT) X(J) = X(J)*A(J,J)
                        END IF
     20             CONTINUE
                ELSE
                    JX = KX
                    DO 40 J = 1,N
                        IF (X(JX).NE.ZERO) THEN
                            TEMP = X(JX)
                            IX = KX
                            DO 30 I = 1,J - 1
                                X(IX) = X(IX) + TEMP*A(I,J)
                                IX = IX + INCX
     30                     CONTINUE
                            IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                        END IF
                        JX = JX + INCX
     40             CONTINUE
                END IF
            ELSE
                IF (INCX.EQ.1) THEN
                    DO 60 J = N,1,-1
                        IF (X(J).NE.ZERO) THEN
                            TEMP = X(J)
                            DO 50 I = N,J + 1,-1
                                X(I) = X(I) + TEMP*A(I,J)
     50                     CONTINUE
                            IF (NOUNIT) X(J) = X(J)*A(J,J)
                        END IF
     60             CONTINUE
                ELSE
                    KX = KX + (N-1)*INCX
                    JX = KX
                    DO 80 J = N,1,-1
                        IF (X(JX).NE.ZERO) THEN
                            TEMP = X(JX)
                            IX = KX
                            DO 70 I = N,J + 1,-1
                                X(IX) = X(IX) + TEMP*A(I,J)
                                IX = IX - INCX
     70                     CONTINUE
                            IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                        END IF
                        JX = JX - INCX
     80             CONTINUE
                END IF
            END IF
        ELSE
  !
  !        Form  x := A**T*x.
  !
            IF (LSAME(UPLO,'U')) THEN
                IF (INCX.EQ.1) THEN
                    DO 100 J = N,1,-1
                        TEMP = X(J)
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO 90 I = J - 1,1,-1
                            TEMP = TEMP + A(I,J)*X(I)
     90                 CONTINUE
                        X(J) = TEMP
    100             CONTINUE
                ELSE
                    JX = KX + (N-1)*INCX
                    DO 120 J = N,1,-1
                        TEMP = X(JX)
                        IX = JX
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO 110 I = J - 1,1,-1
                            IX = IX - INCX
                            TEMP = TEMP + A(I,J)*X(IX)
    110                 CONTINUE
                        X(JX) = TEMP
                        JX = JX - INCX
    120             CONTINUE
                END IF
            ELSE
                IF (INCX.EQ.1) THEN
                    DO 140 J = 1,N
                        TEMP = X(J)
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO 130 I = J + 1,N
                            TEMP = TEMP + A(I,J)*X(I)
    130                 CONTINUE
                        X(J) = TEMP
    140             CONTINUE
                ELSE
                    JX = KX
                    DO 160 J = 1,N
                        TEMP = X(JX)
                        IX = JX
                        IF (NOUNIT) TEMP = TEMP*A(J,J)
                        DO 150 I = J + 1,N
                            IX = IX + INCX
                            TEMP = TEMP + A(I,J)*X(IX)
    150                 CONTINUE
                        X(JX) = TEMP
                        JX = JX + INCX
    160             CONTINUE
                END IF
            END IF
        END IF
  !
        RETURN
  !
  !     End of STRMV .
  !
END SUBROUTINE STRMV


!  =====================================================================
      SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
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
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
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
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.(.NOT.LSAME(TRANSA,'T')) .AND.(.NOT.LSAME(TRANSA,'C'))) THEN
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
          CALL XERBLA('STRSM ',INFO)
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
          DO 20 J = 1,N
              DO 10 I = 1,M
                  B(I,J) = ZERO
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
!
!     Start the operations.
!
      IF (LSIDE) THEN
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*inv( A )*B.
!
              IF (UPPER) THEN
                  DO 60 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 30 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   30                     CONTINUE
                      END IF
                      DO 50 K = M,1,-1
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 40 I = 1,K - 1
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40                         CONTINUE
                          END IF
   50                 CONTINUE
   60             CONTINUE
              ELSE
                  DO 100 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 70 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
   70                     CONTINUE
                      END IF
                      DO 90 K = 1,M
                          IF (B(K,J).NE.ZERO) THEN
                              IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                              DO 80 I = K + 1,M
                                  B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80                         CONTINUE
                          END IF
   90                 CONTINUE
  100             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*inv( A**T )*B.
!
              IF (UPPER) THEN
                  DO 130 J = 1,N
                      DO 120 I = 1,M
                          TEMP = ALPHA*B(I,J)
                          DO 110 K = 1,I - 1
                              TEMP = TEMP - A(K,I)*B(K,J)
  110                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  120                 CONTINUE
  130             CONTINUE
              ELSE
                  DO 160 J = 1,N
                      DO 150 I = M,1,-1
                          TEMP = ALPHA*B(I,J)
                          DO 140 K = I + 1,M
                              TEMP = TEMP - A(K,I)*B(K,J)
  140                     CONTINUE
                          IF (NOUNIT) TEMP = TEMP/A(I,I)
                          B(I,J) = TEMP
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (LSAME(TRANSA,'N')) THEN
!
!           Form  B := alpha*B*inv( A ).
!
              IF (UPPER) THEN
                  DO 210 J = 1,N
                      IF (ALPHA.NE.ONE) THEN
                          DO 170 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  170                     CONTINUE
                      END IF
                      DO 190 K = 1,J - 1
                          IF (A(K,J).NE.ZERO) THEN
                              DO 180 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 200 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  200                     CONTINUE
                      END IF
  210             CONTINUE
              ELSE
                  DO 260 J = N,1,-1
                      IF (ALPHA.NE.ONE) THEN
                          DO 220 I = 1,M
                              B(I,J) = ALPHA*B(I,J)
  220                     CONTINUE
                      END IF
                      DO 240 K = J + 1,N
                          IF (A(K,J).NE.ZERO) THEN
                              DO 230 I = 1,M
                                  B(I,J) = B(I,J) - A(K,J)*B(I,K)
  230                         CONTINUE
                          END IF
  240                 CONTINUE
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(J,J)
                          DO 250 I = 1,M
                              B(I,J) = TEMP*B(I,J)
  250                     CONTINUE
                      END IF
  260             CONTINUE
              END IF
          ELSE
!
!           Form  B := alpha*B*inv( A**T ).
!
              IF (UPPER) THEN
                  DO 310 K = N,1,-1
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 270 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  270                     CONTINUE
                      END IF
                      DO 290 J = 1,K - 1
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 280 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  280                         CONTINUE
                          END IF
  290                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 300 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  300                     CONTINUE
                      END IF
  310             CONTINUE
              ELSE
                  DO 360 K = 1,N
                      IF (NOUNIT) THEN
                          TEMP = ONE/A(K,K)
                          DO 320 I = 1,M
                              B(I,K) = TEMP*B(I,K)
  320                     CONTINUE
                      END IF
                      DO 340 J = K + 1,N
                          IF (A(J,K).NE.ZERO) THEN
                              TEMP = A(J,K)
                              DO 330 I = 1,M
                                  B(I,J) = B(I,J) - TEMP*B(I,K)
  330                         CONTINUE
                          END IF
  340                 CONTINUE
                      IF (ALPHA.NE.ONE) THEN
                          DO 350 I = 1,M
                              B(I,K) = ALPHA*B(I,K)
  350                     CONTINUE
                      END IF
  360             CONTINUE
              END IF
          END IF
      END IF
!
      RETURN
!
!     End of STRSM .
!
END SUBROUTINE STRSM


!  =====================================================================
      SUBROUTINE STRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,INFO )
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           STRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT. LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STRTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 )RETURN
!
!     Check for singularity.
!
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO ) RETURN
   10    CONTINUE
      END IF
      INFO = 0
!
!     Solve A * x = b  or  A**T * x = b.
!
      CALL STRSM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB )
!
      RETURN
!
!     End of STRTRS
!
END SUBROUTINE STRTRS

!  =====================================================================
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
 9999 FORMAT( ' ** On entry to ', A, ' parameter number ', I2, ' had ','an illegal value' )
!
!     End of XERBLA
!
END SUBROUTINE XERBLA



























end module simple_lapacksgels

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
