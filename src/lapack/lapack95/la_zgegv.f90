    SUBROUTINE ZGEGV_F95( A, B, ALPHA, BETA, VL, VR, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO
   USE F77_LAPACK, ONLY: GEGV_F77 => LA_GEGV
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
   COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: ALPHA(:), BETA(:)
   COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
!  LA_GEGV computes for a pair of N-by-N complex nonsymmetric
!  matrices A and B, the generalized eigenvalues (alpha, beta),
!  and optionally, the left and/or right generalized eigenvectors
!  (VL and VR).
!  A generalized eigenvalue for a pair of matrices (A,B) is, roughly
!  speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B
!  is singular.  It is usually represented as the pair (alpha,beta),
!  as there is a reasonable interpretation for beta=0, and even for
!  both being zero.  A good beginning reference is the book, "Matrix
!  Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)
!  A right generalized eigenvector corresponding to a generalized
!  eigenvalue  w  for a pair of matrices (A,B) is a vector  r  such
!  that  (A - w B) r = 0 .  A left generalized eigenvector is a vector
!  l such that l**H * (A - w B) = 0, where l**H is the
!  conjugate-transpose of l.
!  Note: this routine performs "full balancing" on A and B -- see
!  "Further Details", below.
!
! =========
!
!   SUBROUTINE LA_GEGV( A, B, <alpha>, BETA, VL, VR, INFO )
!      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!      <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
!      <type>(<wp>), INTENT(OUT), OPTIONAL :: <alpha'>, BETA(:)
!      <type>(<wp>), INTENT(OUT), OPTIONAL :: VL(:,:), VR(:,:)
!   where
!      <type>   ::= REAL | COMPLEX
!      <wp>     ::= KIND(1.0) | KIND(1.0D0)
!      <alpha>  ::= ALPHAR, ALPHAI | ALPHA
!      <alpha'> ::= ALPHAR(:), ALPHAI(:) | ALPHA(:)
!
! Arguments
! =========
!
! A    (input/output) REAL / COMPLEX array, shape (:,:),
!      SIZE(A,1) == SIZE(A,2) == n.
!      On entry, the first of the pair of matrices whose
!      generalized eigenvalues and (optionally) generalized
!      eigenvectors are to be computed.
!      On exit, the contents will have been destroyed.  (For a
!      description of the contents of A on exit, see "Further
!      Details", below.)
!
! B    (input/output) REAL / COMPLEX array, shape (:,:),
!      SIZE(rBA,1) == SIZE(rBA,2) == n.
!      On entry, the second of the pair of matrices whose
!      generalized eigenvalues and (optionally) generalized
!      eigenvectors are to be computed.
!      On exit, the contents will have been destroyed.  (For a
!      description of the contents of B on exit, see "Further
!      Details", below.)
!
! ALPHAR  Only for the real case. Optional (output) REAL array,
! ALPHAI  shape (:), SIZE(ALPHA') == n. ALPHA ::= ALPHAR | ALPHAI
! ALPHA   Only for the complex case. Optional (output) COMPLEX
!         array, shape (:), SIZE(ALPHAI) == n.
! BETA Optional (output) REAL / COMPLEX array, shape (:),
!      SIZE(BETA) == n.
!      On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,n, will
!      be the generalized eigenvalues.  If ALPHAI(j) is zero, then
!      the j-th eigenvalue is real; if positive, then the j-th and
!      (j+1)-st eigenvalues are a complex conjugate pair, with
!      ALPHAI(j+1) negative.
!      Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
!      may easily over- or underflow, and BETA(j) may even be zero.
!      Thus, the user should avoid naively computing the ratio
!      alpha/beta.  However, ALPHAR and ALPHAI will be always less
!      than and usually comparable with norm(A) in magnitude, and
!      BETA always less than and usually comparable with norm(B).
!
! VL   Optional (output) REAL / COMPLEX array, shape (:,:),
!      SIZE(VL,1) == SIZE(VL,2) == n.
!      The left generalized eigenvectors.  (See "Purpose", above.)
!      Real eigenvectors take one column,
!      complex take two columns, the first for the real part and
!      the second for the imaginary part.  Complex eigenvectors
!      correspond to an eigenvalue with positive imaginary part.
!      Each eigenvector will be scaled so the largest component
!      will have abs(real part) + abs(imag. part) = 1, *except*
!      that for eigenvalues with alpha=beta=0, a zero vector will
!      be returned as the corresponding eigenvector.
!
! VR   Optional (output) REAL / COMPLEX array, shape (:,:),
!      SIZE(VR,1) == SIZE(VR,2) == n.
!      The right generalized eigenvectors.  (See "Purpose", above.)
!      Real eigenvectors take one column,
!      complex take two columns, the first for the real part and
!      the second for the imaginary part.  Complex eigenvectors
!      correspond to an eigenvalue with positive imaginary part.
!      Each eigenvector will be scaled so the largest component
!      will have abs(real part) + abs(imag. part) = 1, *except*
!      that for eigenvalues with alpha=beta=0, a zero vector will
!      be returned as the corresponding eigenvector.
!
! INFO (output) INTEGER
!      = 0:  successful exit
!      < 0:  if INFO = -i, the i-th argument had an illegal value.
!      = 1,...,n:
!            The QZ iteration failed.  No eigenvectors have been
!            calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
!            should be correct for j=INFO+1,...,n.
!      > n:  errors that usually indicate LAPACK problems:
!            =n+1: error return from LA_GGBAL
!            =n+2: error return from LA_GEQRF
!            =n+3: error return from LA_ORMQR
!            =n+4: error return from LA_ORGQR
!            =n+5: error return from LA_GGHRD
!            =n+6: error return from LA_HGEQZ (other than failed
!                                            iteration)
!            =n+7: error return from LA_TGEVC
!            =n+8: error return from LA_GGBAK (computing VL)
!            =n+9: error return from LA_GGBAK (computing VR)
!            =n+10: error return from LA_LASCL (various calls)
!      If INFO is not present and an error occurs, then the program is
!         terminated with an error message.
!
! Further Details
! ===============
!
! Balancing
! ---------
!
! This driver calls SGGBAL to both permute and scale rows and columns
! of A and B.  The permutations PL and PR are chosen so that PL*A*PR
! and PL*B*R will be upper triangular except for the diagonal blocks
! A(i:j,i:j) and B(i:j,i:j), with i and j as close together as
! possible.  The diagonal scaling matrices DL and DR are chosen so
! that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to
! one (except for the elements that start out zero.)
!
! After the eigenvalues and eigenvectors of the balanced matrices
! have been computed, SGGBAK transforms the eigenvectors back to what
! they would have been (in perfect arithmetic) if they had not been
! balanced.
!
! Contents of A and B on Exit
! -------- -- - --- - -- ----
!
! If any eigenvectors are computed (either VL or VR or both), then on
! exit the arrays A and B will contain the real Schur form[*] of the
! "balanced" versions of A and B. If no eigenvectors are computed,
! then only the diagonal blocks will be correct.
!
! [*] See LA_HGEQZ, LA_GEGS, or read the book "Matrix Computations",
!     by Golub & van Loan, pub. by Johns Hopkins U. Press.
!--------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GEGV'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBVL, LJOBVR
   INTEGER, SAVE :: LWORK = 0
   INTEGER :: N, LINFO, LD, ISTAT, ISTAT1, S1VL, S2VL, S1VR, S2VR, &
              SALPHA, SBETA
!  .. LOCAL ARRAYS ..
   COMPLEX(WP), TARGET :: LLVL(1,1), LLVR(1,1)
   COMPLEX(WP), POINTER :: WORK(:), LALPHA(:), LBETA(:)
   REAL(WP), POINTER :: RWORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
   IF( PRESENT(ALPHA) )THEN; SALPHA = SIZE(ALPHA); ELSE; SALPHA = N; ENDIF
   IF( PRESENT(BETA) )THEN; SBETA = SIZE(BETA); ELSE; SBETA = N; ENDIF
   IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
   ELSE; S1VL = 1; S2VL = 1; LJOBVL = 'N'; END IF
   IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
   ELSE; S1VR = 1; S2VR = 1; LJOBVR = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
   ELSE IF( SIZE(B,1) /= N .OR. SIZE(B,2) /= N )THEN; LINFO = -2
   ELSE IF( SALPHA /= N )THEN; LINFO = -3
   ELSE IF( SBETA /= N )THEN; LINFO = -4
   ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -5
   ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -6
   ELSE IF( N > 0 )THEN
      IF( PRESENT(ALPHA) )THEN; LALPHA => ALPHA
      ELSE; ALLOCATE(LALPHA(N),STAT=ISTAT); ENDIF
      IF( ISTAT == 0 )THEN
         IF( PRESENT(BETA) )THEN; LBETA => BETA
         ELSE; ALLOCATE(LBETA(N),STAT=ISTAT); ENDIF
      END IF
      IF( ISTAT == 0 )THEN
         ALLOCATE( RWORK(MAX(1, 8*N)), STAT=ISTAT )
      END IF
      IF( ISTAT == 0 )THEN
         LWORK = MAX( 1, 2*N, LWORK); ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
            LWORK = MAX( 1, 2*N ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
            IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
         END IF
      END IF
      IF( ISTAT == 0 ) THEN
         IF( PRESENT(VL) )THEN
          IF( PRESENT(VR) )THEN
               CALL GEGV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, LALPHA, &
                        LBETA, VL, S1VL, VR, S1VR, WORK, LWORK, RWORK, LINFO )
	     ELSE
	       CALL GEGV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, LALPHA, &
                        LBETA, VL, S1VL, LLVR, S1VR, WORK, LWORK, RWORK, LINFO )
	     ENDIF
	    ELSE
	     IF( PRESENT(VR) )THEN 
	       CALL GEGV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, LALPHA, &
                        LBETA, LLVL, S1VL, VR, S1VR, WORK, LWORK, RWORK, LINFO )
	     ELSE
	       CALL GEGV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, LALPHA, &
                        LBETA, LLVL, S1VL, LLVR, S1VR, WORK, LWORK, RWORK, LINFO )
	     ENDIF
	    ENDIF 
         IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, RWORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZGEGV_F95
