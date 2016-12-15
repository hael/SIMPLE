      SUBROUTINE DGESVX_F95( A, B, X, AF, IPIV, FACT, TRANS, EQUED, R, C, FERR, BERR, RCOND,     &
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
      REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
      REAL(WP), INTENT(OUT) :: X(:,:)
      INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: IPIV(:)
      REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: C(:), R(:)
      REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: AF(:,:)
      REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: FERR(:), BERR(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_GESVX computes the solution to a real or complex linear system of
! equations of the form A*X = B, A^T*X = B or A^H*X = B, where A is a 
! square matrix and X and B are rectangular matrices or vectors.
!    LA_GESVX can also optionally equilibrate the system if A is poorly
! scaled, estimate the condition number of (the equilibrated) A, return 
! the pivot growth factor, and compute error bounds.
! 
! =========
! 
!     SUBROUTINE LA_GESVX ( A, B, X, AF=af, IPIV=ipiv, FACT=fact, &
!                  TRANS=trans, EQUED=equed, R=r, C=c, FERR=ferr, &
!                  BERR=berr, RCOND=rcond, RPVGRW=rpvgrw, &
!                  INFO=info )
!          <type>(<wp>), INTENT(INOUT) :: A(:,:), <rhs>
!          <type>(<wp>), INTENT(OUT) :: <sol>
!          <type>(<wp>), INTENT(INOUT), OPTIONAL :: AF(:,:)
!          INTEGER, INTENT(INOUT), OPTIONAL :: IPIV(:)
!          CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: FACT, &
!                                       TRANS
!          CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: EQUED
!          REAL(<wp>), INTENT(INOUT), OPTIONAL :: R(:), C(:)
!          REAL(<wp>), INTENT(OUT), OPTIONAL :: <err>, RCOND, RPVGRW
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
!          <rhs>  ::= B(:,:) | B(:)
!          <sol>  ::= X(:,:) | X(:)
!          <err>  ::= FERR(:), BERR(:) | FERR, BERR
! 
! Arguments
! =========
! 
! A         (input/output) REAL or COMPLEX square array, shape (:,:).
!           On entry, the matrix A or its equilibration:
!           If FACT = 'F' and EQUED /= 'N' then A has been equilibrated 
!           by the scaling factors in R and/or C during a previous call
! 	    to LA_GESVX.
!           On exit, if FACT = 'E', then the equilibrated version of A
! 	    is stored in A; otherwise, A is unchanged.
! B         (input/output) REAL or COMPLEX array, shape (:,:) with 
!           size(B,1) = size(A,1) or shape (:) with size(B) = size(A,1).
!           On entry, the matrix B.
!           On exit, the scaled version of B if the system has been 
!           equilibrated; otherwise, B is unchanged.
! X         (output) REAL or COMPLEX array, shape (:,:) with size(X,1) =
!           size(A,1) and size(X,2) = size(B,2), or shape (:) with 
! 	    size(X) = size(A,1).
!           The solution matrix X .
! AF        Optional (input or output) REAL or COMPLEX square array, 
!           shape (:,:) with the same size as A.
!           If FACT = 'F' then AF is an input argument that contains the
!           factors L and U of (the equilibrated) A returned by a
! 	    previous call to LA_GESVX.
!           If FACT /= 'F' then AF is an output argument that contains 
!           the factors L and U of (the equilibrated) A.
! IPIV      Optional (input or output) INTEGER array, shape (:) with 
!           size(IPIV) = size(A,1).
!           If FACT = 'F' then IPIV is an input argument that contains 
!           the pivot indices from the factorization of (the 
! 	    equilibrated) A, returned by a previous call to LA_GESVX.
!           If FACT /= 'F' then IPIV is an output argument that contains
!           the pivot indices from the factorization of (the
! 	    equilibrated) A.
! FACT      Optional (input) CHARACTER(LEN=1).
!           Specifies whether the factored form of the matrix A is 
!           supplied on entry, and, if not, whether the matrix A should
!           be equilibrated before it is factored.
!            = 'N': The matrix A will be copied to AF and factored (no 
! 	          equilibration).
!            = 'E': The matrix A will be equilibrated, then copied to AF
! 	          and factored.
!            = 'F': AF and IPIV contain the factored form of (the 
! 	          equilibrated) A.
!           Default value: 'N'.
! TRANS     Optional (input) CHARACTER(LEN=1).
!           Specifies the form of the system of equations:
!            = 'N': A*X = B (No transpose)
!            = 'T': A^T*X = B (Transpose)
!            = 'C': A^H*X = B (Conjugate transpose)
! EQUED     Optional (input or output) CHARACTER(LEN=1).
!           Specifies the form of equilibration that was done.
!           EQUED is an input argument if FACT = 'F', otherwise it is an
!           output argument:
!            = 'N': No equilibration (always true if FACT = 'N').
!            = 'R': Row equilibration, i.e., A has been premultiplied by
! 	          diag(R).
!            = 'C': Column equilibration, i.e., A has been postmultiplied 
! 	          by diag(C).
!            = 'B': Both row and column equilibration.
!           Default value: 'N'.
! R         Optional (input or output) REAL array, shape (:) with size(R)
!           = size(A,1). The row scale factors for A.
!           R is an input argument if FACT = 'F' and EQUED = 'R' or 'B'.
!           R is an output argument if FACT = 'E' and EQUED = 'R' or 'B'.
! C         Optional (input or output) REAL array, shape (:) with size(C) 
!           = size(A,1). The column scale factors for A.
!           C is an input argument if FACT = 'F' and EQUED = 'C' or 'B'.
!           C is an output argument if FACT = 'E' and EQUED = 'C' or 'B'.
! FERR      Optional (output) REAL array of shape (:), with size(FERR) =
!           size(X,2), or REAL scalar.
!           The estimated forward error bound for each solution vector 
!           X(j) (the j-th column of the solution matrix X). If XTRUE is 
!           the true solution corresponding to X(j) , FERR(j) is an 
! 	    estimated upper bound for the magnitude of the largest 
! 	    element in (X(j)-XTRUE) divided by the magnitude of the 
! 	    largest element in X(j). The estimate is as reliable as the
!           estimate for RCOND and is almost always a slight 
! 	    overestimate of the true error.
! BERR      Optional (output) REAL array of shape (:), with size(BERR) =
!           size(X,2), or REAL scalar.
!           The componentwise relative backward error of each solution 
!           vector X(j) (i.e., the smallest relative change in any
!           element of A or B that makes X(j) an exact solution).
! RCOND     Optional (output) REAL.
!           The estimate of the reciprocal condition number of (the 
!           equilibrated) A. If RCOND is less than the machine precision,
!           the matrix is singular to working precision. This condition 
!           is indicated by a return code of INFO > 0.
! RPVGRW    Optional (output) REAL.
!           The reciprocal pivot growth factor ||A||inf = ||U||inf. If
!           RPVGRW is much less than 1, then the stability of the LU 
!           factorization of the (equilibrated) matrix A could be poor.
!           This also means that the solution X , condition estimator 
!           RCOND, and forward error bound FERR could be unreliable. If
!           the factorization fails with 0 < INFO <= size(A,1), then
!           RPVGRW contains the reciprocal pivot growth factor for the 
!           leading INFO columns of A.
! INFO      Optional (output) INTEGER
!           = 0: successful exit.
!           < 0: if INFO = -i, the i-th argument had an illegal value.
!           > 0: if INFO = i, and i is
!               <= n: U(i,i) = 0. The factorization has been completed,
!                    but the factor U is singular, so the solution could 
! 		   not be computed.
!               = n+1: U is nonsingular, but RCOND is less than machine 
! 	           precision, so the matrix is singular to working 
! 		   precision. Nevertheless, the solution and error
!                    bounds are computed because the computed solution 
! 		   can be more accurate than the value of RCOND would 
! 		   suggest.
!           If INFO is not present and an error occurs, then the program 
! 	    is terminated with an error message.
!----------------------------------------------------------------------
!     .. PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GESVX'
!     .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LFACT, LTRANS, LEQUED
      INTEGER :: ISTAT, ISTAT1, LD, LINFO, N, NRHS, S1AF, S2AF, SBERR, SC, SFERR, SIPIV, SR
      REAL(WP) :: LRCOND, MVR, MVC
!     .. LOCAL POINTERS ..
      INTEGER, POINTER :: IWORK(:), LPIV(:)
      REAL(WP),  POINTER :: LC(:), LR(:), LFERR(:), LBERR(:)
      REAL(WP),  POINTER :: WORK(:), LAF(:, :)
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX, PRESENT, SIZE, MINVAL, TINY
!     .. EXECUTABLE STATEMENTS ..
      LINFO = 0; ISTAT = 0; N = SIZE(A, 1); NRHS = SIZE(B, 2)
      LD = MAX(1,N)
      IF( PRESENT(RCOND) ) RCOND = 1.0_WP
      IF( PRESENT(RPVGRW) ) RPVGRW = 1.0_WP
      IF( PRESENT(FACT) )THEN
         LFACT = FACT
      ELSE
         LFACT='N'
      END IF
      IF( PRESENT(EQUED) .AND. LSAME(LFACT,'F') )THEN
         LEQUED = EQUED
      ELSE
         LEQUED='N'
      END IF
      IF( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = N
      END IF
      IF( PRESENT(AF) )THEN
         S1AF = SIZE(AF,1); S2AF = SIZE(AF,2)
      ELSE
         S1AF = N; S2AF = N
      END IF
      IF( ( PRESENT(C) ) )THEN
         SC = SIZE(C)
      ELSE
         SC = N
      END IF
      IF( ( PRESENT(C) .AND. LSAME(LFACT,'F') ) .AND.                   &
     &    ( LSAME(LEQUED,'C') .OR. LSAME(LEQUED,'B') ) )THEN
         MVC = MINVAL(C)
      ELSE
         MVC = TINY(1.0_WP)
      END IF
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
      IF( PRESENT(FERR) )THEN
         SFERR = SIZE(FERR)
      ELSE
         SFERR = NRHS
      END IF
      IF( PRESENT(BERR) )THEN
         SBERR = SIZE(BERR)
      ELSE
         SBERR = NRHS
      END IF
      IF(PRESENT(TRANS))THEN
         LTRANS = TRANS
      ELSE
         LTRANS='N'
      END IF
!     .. TEST THE ARGUMENTS
      IF( SIZE(A, 2) /= N .OR. N < 0 )THEN
         LINFO = -1
      ELSE IF( SIZE(B, 1) /= N .OR. NRHS < 0 )THEN
         LINFO = -2
      ELSE IF( SIZE(X, 1) /= N .OR. SIZE(X, 2) /= NRHS )THEN
         LINFO = -3
      ELSE IF( S1AF /= N .OR. S2AF /= N ) THEN
         LINFO = -4
      ELSE IF( SIPIV /= N )THEN
         LINFO = -5
      ELSE IF( SR /= N .OR. MVR <= 0.0_WP )THEN
         LINFO = -9
      ELSE IF( SC /= N .OR. MVC <= 0.0_WP )THEN
         LINFO = -10
      ELSE IF( SFERR /= NRHS )THEN
         LINFO = -11
      ELSE IF( SBERR /= NRHS )THEN
         LINFO = -12
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
         IF( ISTAT == 0 )THEN
            IF( .NOT.PRESENT(IPIV) )THEN
               ALLOCATE( LPIV(N), STAT=ISTAT )
            ELSE
               LPIV => IPIV
            END IF
         END IF
         IF( ISTAT == 0 )THEN
            IF( .NOT.PRESENT(R) )THEN
               ALLOCATE( LR(N), STAT=ISTAT )
            ELSE
               LR => R
            END IF
         END IF
         IF( ISTAT == 0 )THEN
            IF( .NOT.PRESENT(C) )THEN
               ALLOCATE( LC(N), STAT=ISTAT )
            ELSE
               LC => C
            END IF
         END IF
         IF( ISTAT == 0 )THEN
            IF( .NOT.PRESENT(FERR) )THEN
               ALLOCATE( LFERR(NRHS), STAT=ISTAT )
            ELSE
               LFERR => FERR
            END IF
         END IF
         IF( ISTAT == 0 )THEN
            IF( .NOT.PRESENT(BERR) )THEN
               ALLOCATE( LBERR(NRHS), STAT=ISTAT )
            ELSE
               LBERR => BERR
            END IF
         END IF
         IF( ISTAT == 0 )THEN
            ALLOCATE(WORK(4*N), IWORK(N), STAT=ISTAT )
         END IF
         IF( ISTAT == 0 )THEN
!           .. CALL LAPACK77 ROUTINE
            CALL GESVX_F77( LFACT, LTRANS, N, NRHS, A, LD, LAF, LD, LPIV, LEQUED, LR, LC, B, LD, &
                            X, LD, LRCOND, LFERR, LBERR, WORK, IWORK, LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT.PRESENT(R) ) DEALLOCATE( LR, STAT=ISTAT1 )
         IF( .NOT.PRESENT(C) ) DEALLOCATE( LC, STAT=ISTAT1 )
         IF( .NOT.PRESENT(AF) ) DEALLOCATE( LAF, STAT=ISTAT1 )
         IF( .NOT.PRESENT(IPIV) ) DEALLOCATE( LPIV, STAT=ISTAT1 )
         IF( .NOT.PRESENT(FERR) ) DEALLOCATE( LFERR, STAT=ISTAT1 )
         IF( .NOT.PRESENT(BERR) ) DEALLOCATE( LBERR, STAT=ISTAT1 )
         IF( PRESENT(RPVGRW) ) RPVGRW=WORK(1)
         IF( PRESENT(RCOND) ) RCOND=LRCOND
         IF( PRESENT(EQUED) .AND. .NOT.LSAME(LFACT,'F') ) EQUED=LEQUED
         DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
      END IF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
      END SUBROUTINE DGESVX_F95
