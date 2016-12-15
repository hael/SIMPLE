    SUBROUTINE DGEGS_F95( A, B, ALPHAR, ALPHAI, BETA, VSL, VSR, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO
   USE F77_LAPACK, ONLY: GEGS_F77 => LA_GEGS
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: ALPHAR(:), ALPHAI(:), BETA(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: VSL(:,:), VSR(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
!  LA_GEGS computes for a pair of n-by-n real nonsymmetric matrices
!  A, B: the generalized eigenvalues (alphar alpha, beta),
!  the Schur form (A, B), and optionally left and/or right
!  Schur vectors (VSL and VSR).
!
!  (If only the generalized eigenvalues are needed, use the driver SGEGV
!  instead.)
!
!  A generalized eigenvalue for a pair of matrices (A,B) is, roughly
!  speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B
!  is singular.  It is usually represented as the pair (alpha,beta),
!  as there is a reasonable interpretation for beta=0, and even for
!  both being zero.  A good beginning reference is the book, "Matrix
!  Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)
!
!  The (generalized) Schur form of a pair of matrices is the result of
!  multiplying both matrices on the left by one orthogonal matrix and
!  both on the right by another orthogonal matrix, these two orthogonal
!  matrices being chosen so as to bring the pair of matrices into
!  (real) Schur form.
!
!  A pair of matrices A, B is in generalized real Schur form if B is
!  upper triangular with non-negative diagonal and A is block upper
!  triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
!  to real generalized eigenvalues, while 2-by-2 blocks of A will be
!  "standardized" by making the corresponding elements of B have the
!  form:
!          [  a  0  ]
!          [  0  b  ]
!
!  and the pair of corresponding 2-by-2 blocks in A and B will
!  have a complex conjugate pair of generalized eigenvalues.
!
!  The left and right Schur vectors are the columns of VSL and VSR,
!  respectively, where VSL and VSR are the orthogonal matrices
!  which reduce A and B to Schur form:
!
!  Schur form of (A,B) = ( (VSL)**T A (VSR), (VSL)**T B (VSR) )
!
! =========
!
!   SUBROUTINE LA_GEGS( A, B, <alpha>, BETA, VSL, VSR, INFO )
!      <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
!      <type>(<wp>), INTENT(OUT), OPTIONAL :: <alpha'>, BETA(:)
!      <type>(<wp>), INTENT(OUT), OPTIONAL :: VSL(:,:), VSR(:,:)
!      INTEGER, INTENT(OUT), OPTIONAL :: INFO
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
!      On entry, the first of the pair of matrices whose generalized
!      eigenvalues and (optionally) Schur vectors are to be
!      computed.
!      On exit, the generalized Schur form of A.
!      Note: to avoid overflow, the Frobenius norm of the matrix
!      A should be less than the overflow threshold.
!
! B    (input/output) REAL / COMPLEX array, shape (:,:),
!      SIZE(rBA,1) == SIZE(rBA,2) == n.
!      On entry, the second of the pair of matrices whose
!      generalized eigenvalues and (optionally) Schur vectors are
!      to be computed.
!      On exit, the generalized Schur form of B.
!      Note: to avoid overflow, the Frobenius norm of the matrix
!      B should be less than the overflow threshold.
!
! ALPHAR  Only for the real case. Optional (output) REAL array,
! ALPHAI  shape (:), SIZE(ALPHA') == n. ALPHA ::= ALPHAR | ALPHAI
! ALPHA   Only for the complex case. Optional (output) COMPLEX
!         array, shape (:), SIZE(ALPHAI) == n.
! BETA Optional (output) REAL / COMPLEX array, shape (:),
!      SIZE(BETA) == n.
!      On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,n, will
!      be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
!      j=1,...,n  and  BETA(j),j=1,...,n  are the diagonals of the
!      complex Schur form (A,B) that would result if the 2-by-2
!      diagonal blocks of the real Schur form of (A,B) were further
!      reduced to triangular form using 2-by-2 complex unitary
!      transformations.  If ALPHAI(j) is zero, then the j-th
!      eigenvalue is real; if positive, then the j-th and (j+1)-st
!      eigenvalues are a complex conjugate pair, with ALPHAI(j+1)
!      negative.
!      Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
!      may easily over- or underflow, and BETA(j) may even be zero.
!      Thus, the user should avoid naively computing the ratio
!      alpha/beta.  However, ALPHAR and ALPHAI will be always less
!      than and usually comparable with norm(A) in magnitude, and
!      BETA always less than and usually comparable with norm(B).
!
! VSL  Optional (output) REAL / COMPLEX array, shape (:,:),
!      SIZE(VSL,1) == SIZE(VSL,2) == n.
!      VSL will contain the left Schur vectors. (See "Purpose", above.)
!
! VSR  Optional (output) REAL / COMPLEX array, shape (:,:),
!      SIZE(VSL,1) == SIZE(VSL,2) == n.
!      VSR will contain the right Schur vectors. (See "Purpose", above.)
!
! INFO (output) INTEGER
!      = 0:  successful exit
!      < 0:  if INFO = -i, the i-th argument had an illegal value.
!      = 1,...,n:
!            The QZ iteration failed.  (A,B) are not in Schur
!            form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
!            be correct for j=INFO+1,...,n.
!      > n:  errors that usually indicate LAPACK problems:
!            =n+1: error return from LA_GGBAL
!            =n+2: error return from LA_GEQRF
!            =n+3: error return from LA_ORMQR
!            =n+4: error return from LA_ORGQR
!            =n+5: error return from LA_GGHRD
!            =n+6: error return from LA_HGEQZ (other than failed
!                                            iteration)
!            =n+7: error return from LA_GGBAK (computing VSL)
!            =n+8: error return from LA_GGBAK (computing VSR)
!            =n+9: error return from LA_LASCL (various places)
!      If INFO is not present and an error occurs, then the program is
!         terminated with an error message.
!-------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GEGS'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBVSL, LJOBVSR
   INTEGER, SAVE :: LWORK = 0
   INTEGER :: N, LINFO, LD, ISTAT, ISTAT1, S1VSL, S2VSL, S1VSR, S2VSR, &
              SALPHAR, SALPHAI, SBETA
!  .. LOCAL ARRAYS ..
   REAL(WP), TARGET :: LLVSL(1,1), LLVSR(1,1)
   REAL(WP), POINTER :: WORK(:), &
     &                   LALPHAR(:), LALPHAI(:), LBETA(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
   IF( PRESENT(ALPHAR) )THEN; SALPHAR = SIZE(ALPHAR); ELSE; SALPHAR = N; ENDIF
   IF( PRESENT(ALPHAI) )THEN; SALPHAI = SIZE(ALPHAI); ELSE; SALPHAI = N; ENDIF
   IF( PRESENT(BETA) )THEN; SBETA = SIZE(BETA); ELSE; SBETA = N; ENDIF
   IF( PRESENT(VSL) )THEN; S1VSL = SIZE(VSL,1); S2VSL = SIZE(VSL,2); LJOBVSL = 'V'
   ELSE; S1VSL = 1; S2VSL = 1; LJOBVSL = 'N'; END IF
   IF( PRESENT(VSR) )THEN; S1VSR = SIZE(VSR,1); S2VSR = SIZE(VSR,2); LJOBVSR = 'V'
   ELSE; S1VSR = 1; S2VSR = 1; LJOBVSR = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
   ELSE IF( SIZE(B,1) /= N .OR. SIZE(B,2) /= N )THEN; LINFO = -2
   ELSE IF( SALPHAR /= N )THEN; LINFO = -3
   ELSE IF( SALPHAI /= N )THEN; LINFO = -4
   ELSE IF( SBETA /= N )THEN; LINFO = -5
   ELSE IF( PRESENT(VSL) .AND. ( S1VSL /= N .OR. S2VSL /= N ) )THEN; LINFO = -6
   ELSE IF( PRESENT(VSR) .AND. ( S1VSR /= N .OR. S2VSR /= N ) )THEN; LINFO = -7
   ELSE IF( N > 0 )THEN
      IF( PRESENT(ALPHAR) )THEN; LALPHAR => ALPHAR
      ELSE; ALLOCATE(LALPHAR(N),STAT=ISTAT); ENDIF
      IF( ISTAT == 0 )THEN
         IF( PRESENT(ALPHAI) )THEN; LALPHAI => ALPHAI
         ELSE; ALLOCATE(LALPHAI(N),STAT=ISTAT); ENDIF
      END IF
      IF( ISTAT == 0 )THEN
         IF( PRESENT(BETA) )THEN; LBETA => BETA
         ELSE; ALLOCATE(LBETA(N),STAT=ISTAT); ENDIF
      END IF
      IF( ISTAT == 0 )THEN
         LWORK = MAX( 1, 4*N, LWORK); ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
            LWORK = MAX( 1, 4*N ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
            IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
         END IF
      END IF
      IF( ISTAT == 0 ) THEN
         IF( PRESENT(VSL) )THEN
           IF( PRESENT(VSR) )THEN
              CALL GEGS_F77( LJOBVSL, LJOBVSR, N, A, LD, B, LD, LALPHAR, LALPHAI, &
                        LBETA, VSL, S1VSL, VSR, S1VSR, WORK, LWORK, LINFO )
	   ELSE
	      CALL GEGS_F77( LJOBVSL, LJOBVSR, N, A, LD, B, LD, LALPHAR, LALPHAI, &
                        LBETA, VSL, S1VSL, LLVSR, S1VSR, WORK, LWORK, LINFO )
	   ENDIF
	 ELSE
	   IF( PRESENT(VSR) )THEN 
	      CALL GEGS_F77( LJOBVSL, LJOBVSR, N, A, LD, B, LD, LALPHAR, LALPHAI, &
                        LBETA, LLVSL, S1VSL, VSR, S1VSR, WORK, LWORK, LINFO )
	   ELSE
	      CALL GEGS_F77( LJOBVSL, LJOBVSR, N, A, LD, B, LD, LALPHAR, LALPHAI, &
                        LBETA, LLVSL, S1VSL, LLVSR, S1VSR, WORK, LWORK, LINFO )
	   ENDIF
	 ENDIF  
         IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DGEGS_F95
