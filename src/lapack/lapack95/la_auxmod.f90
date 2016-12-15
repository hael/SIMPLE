      MODULE LA_PRECISION
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  DEFINES SINGLE AND DOUBLE PRECISION PARAMETERS, SP AND DP.
!  THESE VALUES ARE COMPILER DEPENDENT.
!
      INTEGER, PARAMETER :: SP=KIND(1.0), DP=KIND(1.0D0)
!
      END MODULE LA_PRECISION

      MODULE LA_AUXMOD
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
      INTERFACE
         SUBROUTINE    ERINFO(LINFO, SRNAME, INFO, ISTAT)
            CHARACTER( LEN = * ), INTENT(IN)              :: SRNAME
            INTEGER             , INTENT(IN)              :: LINFO
            INTEGER             , INTENT(OUT), OPTIONAL   :: INFO
            INTEGER             , INTENT(IN), OPTIONAL    :: ISTAT
         END SUBROUTINE ERINFO
         INTEGER FUNCTION LA_WS_GELS( VER, M, N, NRHS, TRANS )
             CHARACTER( LEN=1 ), INTENT(IN) :: TRANS, VER
             INTEGER, INTENT(IN) :: M, N, NRHS
         END FUNCTION LA_WS_GELS
         INTEGER FUNCTION LA_WS_GELSS( VER, M, N, NRHS )
          CHARACTER(LEN=1), INTENT(IN) :: VER
          INTEGER, INTENT(IN) :: M, N, NRHS
         END FUNCTION LA_WS_GELSS
      END INTERFACE
!
      CONTAINS
!
      LOGICAL FUNCTION LSAME( CA, CB )
!
!  PURPOSE
!  =======
!
!  LSAME  TESTS IF CA IS THE SAME LETTER AS CB REGARDLESS OF CASE.
!
!  PARAMETERS
!  ==========
!
!  CA      (INPUT) CHARACTER*1
!  CB      (INPUT) CHARACTER*1
!          CHARACTERS TO BE COMPARED.
!
!  .. SCALAR ARGUMENTS ..
      CHARACTER*1, INTENT(IN) :: CA, CB
!  .. PARAMETERS ..
      INTEGER, PARAMETER      :: IOFF=32
!  .. LOCAL SCALARS ..
      INTEGER                 :: INTA, INTB, ZCODE
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC                  ICHAR
!
!  .. EXECUTABLE STATEMENTS ..
!
!  TEST IF THE CHARACTERS ARE EQUAL
!
      LSAME = CA == CB
!
!  NOW TEST FOR EQUIVALENCE
!
      IF( .NOT.LSAME )THEN
!
!     USE 'Z' RATHER THAN 'A' SO THAT ASCII CAN BE DETECTED ON PRIME
!     MACHINES, ON WHICH ICHAR RETURNS A VALUE WITH BIT 8 SET.
!     ICHAR('A') ON PRIME MACHINES RETURNS 193 WHICH IS THE SAME AS
!     ICHAR('A') ON AN EBCDIC MACHINE.
!
         ZCODE = ICHAR( 'Z' )
!
         INTA = ICHAR( CA )
         INTB = ICHAR( CB )
!
         IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 )THEN
!
!        ASCII IS ASSUMED - ZCODE IS THE ASCII CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
            IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
            IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
         ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 )THEN
!
!        EBCDIC IS ASSUMED - ZCODE IS THE EBCDIC CODE OF EITHER LOWER OR
!        UPPER CASE 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
!    &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
         ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 )THEN
!
!        ASCII IS ASSUMED, ON PRIME MACHINES - ZCODE IS THE ASCII CODE
!        PLUS 128 OF EITHER LOWER OR UPPER CASE 'Z'.
!
            IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
         ENDIF
         LSAME = INTA == INTB
      ENDIF
      END FUNCTION LSAME

      END MODULE LA_AUXMOD
