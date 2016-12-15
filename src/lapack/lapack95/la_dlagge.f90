      SUBROUTINE DLAGGE_F95( A, KL, KU, D, ISEED, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!     .. "Use Statements" ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO
      USE F77_LAPACK, ONLY: LAGGE_F77 => LA_LAGGE
!     .. "Implicit Statement" ..
      IMPLICIT NONE
!     .. "Scalar Arguments" ..
      INTEGER, INTENT(IN), OPTIONAL :: KL, KU
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: ISEED(4)
      REAL(WP), INTENT(IN), OPTIONAL, TARGET :: D(:)
      REAL(WP), INTENT(OUT) :: A(:,:)
!-----------------------------------------------------------------
!
!  Purpose
!  =======
!
!  LA_LAGGE generates a real general m by n matrix A, by pre- and post-
!  multiplying a real diagonal matrix D with random orthogonal matrices:
!  A = U*D*V. The lower and upper bandwidths may then be reduced to
!  kl and ku by additional orthogonal transformations.
!
!  Arguments
!  =========
!
!  SUBROUTINE LA_LAGGE( A, KL, KU, D, ISEED, INFO )
!     <type>(<wp>), INTENT(OUT) :: A(:,:)
!     INTEGER, INTENT(IN), OPTIONAL :: KL, KU
!     REAL(<wp>), INTENT(IN), OPTIONAL, TARGET :: D(:)
!     INTEGER, INTENT(INOUT), OPTIONAL :: ISEED(4)
!     INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     <type> ::= REAL | COMPLEX
!     <wp>   ::= KIND(1.0) | KIND(1.0D0)
!
!  =====================
!
!  A       (output) REAL array, shape (:,:), SIZE(A,1) == m,
!          SIZE(A,2) == n.
!          The generated m by n matrix A.
!
!  KL      (input) INTEGER
!          The number of nonzero subdiagonals within the band of A.
!          0 <= KL <= M-1.
!
!  KU      (input) INTEGER
!          The number of nonzero superdiagonals within the band of A.
!          0 <= KU <= N-1.
!
!  D       (input) REAL array, dimension (min(M,N))
!          The diagonal elements of the diagonal matrix D.
!
!  ISEED   Optional (input/output) INTEGER array, shape (:),
!          SIZE(ISEED) == 4.
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!      If INFO is not present and an error occurs, then the program is
!         terminated with an error message.
!
!  ---------------------------------------------------------------------
!     .. "Parameters" ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_LAGGE'
!     .. "LOCAL Scalars" ..
      INTEGER :: ISTAT, ISTAT1, LDA, LINFO, LKL, LKU, M, MN, N, SD, SISEED
!     .. "Local Arrays" ..
      INTEGER :: LISEED(4)
!     .. "Local Pointers" ..
      REAL(WP), POINTER :: LD(:)
      REAL(WP), POINTER :: WORK(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX, MIN
!     .. "Executable Statements" ..
      LINFO = 0; M = SIZE(A,1); N = SIZE(A,2); LDA = MAX(1,M)
      MN = MIN(M,N); ISTAT = 0
      IF( PRESENT(KL) )THEN; LKL = KL; ELSE; LKL = M-1; ENDIF
      IF( PRESENT(KU) )THEN; LKU = KU; ELSE; LKU = M-1; ENDIF
      IF( PRESENT(D) )THEN; SD = SIZE(D); ELSE; SD = MN; ENDIF
      IF( PRESENT(ISEED) )THEN; SISEED = SIZE(ISEED); LISEED = ISEED
      ELSE; SISEED = 4
            LISEED(1) = 15; LISEED(2) = 1926
            LISEED(3) = 16; LISEED(4) = 1931; ENDIF
!     .. "Test the arguments" ..
      IF( M < 0 .OR. N < 0 )THEN; LINFO = -1
      ELSE IF( LKL < 0 .OR. LKL > M-1 )THEN; LINFO = -2
      ELSE IF( LKU < 0 .OR. LKU > N-1 )THEN; LINFO = -3
      ELSE IF( SD /= MN )THEN; LINFO = -4
      ELSE IF( SISEED /= 4 .OR. PRESENT(ISEED) .AND. (                  &
     &   ( LISEED(1) < 0 .OR. LISEED(1) > 4095 ) .OR.                   &
     &   ( LISEED(2) < 0 .OR. LISEED(2) > 4095 ) .OR.                   &
     &   ( LISEED(3) < 0 .OR. LISEED(3) > 4095 ) .OR.                   &
     &   ( LISEED(4) < 0 .OR. LISEED(4) > 4095 ) .OR.                   &
     &   ( MOD( LISEED(4), 2 ) == 0 ) ) )THEN; LINFO = -5
      ELSE
         IF( PRESENT(D) )THEN; LD => D
         ELSE; ALLOCATE( LD( MN ), STAT=ISTAT); ENDIF
         IF( ISTAT == 0 )THEN
            ALLOCATE( WORK( M + N ), STAT=ISTAT )
            IF( ISTAT == 0 )THEN
               IF( .NOT. PRESENT(D) )THEN
                  LD(1:MN-1) = 1.0_WP; LD(MN) = 0.5_WP; ENDIF
               CALL LAGGE_F77( M, N, LKL, LKU, LD, A, LDA, LISEED, WORK, LINFO )
               IF( PRESENT(ISEED) )ISEED = LISEED
            ELSE; LINFO = -100; END IF
            DEALLOCATE( WORK, STAT=ISTAT1 )
         ENDIF
         IF( .NOT. PRESENT(D) )DEALLOCATE( LD, STAT=ISTAT1)
      ENDIF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
      END SUBROUTINE DLAGGE_F95
