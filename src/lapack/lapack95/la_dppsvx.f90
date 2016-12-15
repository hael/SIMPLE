SUBROUTINE DPPSVX_F95(A, B, X, UPLO, AF, FACT, EQUED, S, FERR, BERR, RCOND, INFO)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: LSAME, ERINFO
   USE F77_LAPACK, ONLY: PPSVX_F77 => LA_PPSVX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, FACT
   CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: EQUED
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
   REAL(WP), INTENT(OUT), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:), B(:,:)
   REAL(WP), INTENT(OUT) :: X(:,:)
   REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: S(:)
   REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: AF(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: FERR(:), BERR(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_PPSVX computes the solution to a linear system of equations 
! A*X = B, where A is real symmetric or complex Hermitian and, in either
! case, positive definite, and where X and B are rectangular matrices or 
! vectors. A is stored in packed format.
!     LA_PPSVX can also optionally equilibrate the system if A is poorly
! scaled, estimate the condition number of (the equilibrated) A, and 
! compute error bounds.
! 
! =========
! 
!       SUBROUTINE LA_PPSVX( AP, B, X, UPLO=uplo, AFP=afp, &
!                  FACT=fact, EQUED=equed, S=s, FERR=ferr, &
!                  BERR=berr, RCOND=rcond, INFO=info )
!           <type>(<wp>), INTENT(INOUT) :: AP(:), <rhs>
!           <type>(<wp>), INTENT(OUT) :: <sol>
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!           <type>(<wp>), INTENT(INOUT), OPTIONAL :: AFP(:)
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: FACT
!           CHARACTER(LEN=1), INTENT(INOUT), OPTIONAL :: EQUED
!           REAL(<wp>), INTENT(INOUT), OPTIONAL :: S(:)
!           REAL(<wp>), INTENT(OUT), OPTIONAL :: <err>
!           REAL(<wp>), INTENT(OUT), OPTIONAL :: RCOND
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!           <type> ::= REAL | COMPLEX
!           <wp>   ::= KIND(1.0) | KIND(1.0D0)
!           <rhs>  ::= B(:,:) | B(:)
!           <sol>  ::= X(:,:) | X(:)
!           <err>  ::= FERR(:), BERR(:) | FERR, BERR
! 
! Arguments
! =========
! 
! AP      (input/output) REAL or COMPLEX square array, shape (:) with
!         size(AP) = n*(n + 1)=2, where n is a rank of the matrix A.
!         On entry, the upper or lower triangle of matrix A, or its 
!         equilibration, in packed storage. The elements are stored 
! 	  columnwise as follows:
! 	  if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j<=n;
! 	  if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for 1<=j<=i<=n.
!         On exit, if FACT = 'E', then the equilibrated version of A is
! 	  stored in AP; otherwise, AP is unchanged.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = n or shape (:) with size(B) = n.
!         On entry, the matrix B.
!         On exit, the scaled version of B if the system has been 
!         equilibrated; otherwise, B is unchanged.
! X       (output) REAL or COMPLEX array, shape (:,:) with size(X,1)=n
!         and size(X,2) = size(B,2), or shape (:) with size(X) = n.
!         The solution matrix X .
! UPLO    Optional (input) CHARACTER(LEN=1).
!           = 'U': Upper triangle of A is stored;
!           = 'L': Lower triangle of A is stored.
!         Default value: 'U'.
! AFP     Optional (input or output) REAL or COMPLEX array, shape (:) 
!         with the same size as AP.
!         If FACT = 'F' then AFP is an input argument that contains the
! 	  factor U or L from the Cholesky factorization of (the 
! 	  equilibrated) A, in the same storage format as A, returned by
! 	  a previous call to LA_PPSVX.
!         If FACT 6= 'F' then AFP is an output argument that contains the
!         factor U or L from the Cholesky factorization of (the
!   	  equilibrated) A in the same storage format as A.
! FACT    Optional (input) CHARACTER(LEN=1).
!         Specifies whether the factored form of the matrix A is supplied
! 	  on entry, and, if not, whether A should be equilibrated before
! 	  it is factored.
!           = 'N': The matrix A will be copied to AFP and factored 
! 	         (no equilibration).
!           = 'E': The matrix A will be equilibrated, then copied to AFP 
! 	         and factored.
!           = 'F': AFP contains the factored form of (the equilibrated) A.
!         Default value: 'N'.
! EQUED   Optional (input or output) CHARACTER(LEN=1).
!         Specifies the form of equilibration that was done.
!         EQUED is an input argument if FACT = 'F', otherwise it is an 
!         output argument:
!           = 'N': No equilibration (always true if FACT = 'N').
!           = 'Y': Equilibration, i.e., A has been premultiplied and 
!                postmultiplied by diag(S).
!         Default value: 'N'.
! S       Optional (input or output) REAL array, shape (:) with size(S)
!         = size(A,1).
!         The scaling factors for A.
!         S is an input argument if FACT = 'F' and EQUED = 'Y'.
!         S is an output argument if FACT = 'E' and EQUED = 'Y'.
! FERR    Optional (output) REAL array of shape (:), with size(FERR) = 
!         size(X,2), or REAL scalar.
!         The estimated forward error bound for each solution vector X(j)
! 	  (the j-th column of the solution matrix X). If XTRUE is the
! 	  true solution corresponding to X(j) , FERR(j) is an estimated 
! 	  upper bound for the magnitude of the largest element in 
! 	  (X(j)-XTRUE) divided by the magnitude of the largest element 
! 	  in X(j). The estimate is as reliable as the estimate for RCOND,
! 	  and is almost always a slight overestimate of the true error.
! BERR    Optional (output) REAL array of shape (:), with size(BERR) = 
!         size(X,2), or REAL scalar.
!         The componentwise relative backward error of each solution 
! 	  vector X(j) (i.e., the smallest relative change in any element
! 	  of A or B that makes X(j) an exact solution).
! RCOND   Optional (output) REAL
!         The estimate of the reciprocal condition number of (the 
! 	  equilibrated) A. If RCOND is less than the machine precision, 
! 	  the matrix is singular to working precision. This condition is
!         indicated by a return code of INFO > 0.
! INFO    Optional (output) INTEGER
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, and i is
!             <= n: the leading minor of order i of (the equilibrated) A
! 	        is not positive definite, so the factorization could 
! 		not be completed and the solution and error bounds 
! 		could not be computed. RCOND= 0 is returned.
!             = n+1: U or L is nonsingular, but RCOND is less than 
! 	        machine precision, so the matrix is singular to
! 		working precision. Nevertheless, the solution and 
! 		error bounds are computed because the computed solution
! 		can be more accurate than the value of RCOND would 
! 		suggest.
!         If INFO is not present and an error occurs, then the program is
! 	  terminated with an error message.
!-----------------------------------------------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_PPSVX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LFACT, LUPLO, LEQUED
   INTEGER :: LINFO, NRHS, N, NN, ISTAT, ISTAT1, SAF, SS, SFERR, SBERR
   REAL(WP) :: LRCOND, MVS
   COMPLEX(WP) :: WW
!  .. LOCAL POINTERS ..
   INTEGER, POINTER :: IWORK(:)
   REAL(WP),  POINTER :: LS(:), LFERR(:), LBERR(:)
   REAL(WP),  POINTER :: WORK(:), LAF( :)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, SIZE, MINVAL, TINY, REAL, INT, AIMAG
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; NN = SIZE(A, 1); NRHS = SIZE(B, 2)
   WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW)
   IF( PRESENT(RCOND) ) RCOND = 1.0_WP
   IF( PRESENT(FACT) )THEN; LFACT = FACT; ELSE; LFACT='N'; END IF
   IF( PRESENT(EQUED) .AND. LSAME(LFACT,'F') )THEN; LEQUED = EQUED
   ELSE; LEQUED='N'; END IF
   IF( PRESENT(AF) )THEN; SAF = SIZE(AF); ELSE; SAF = N*(N+1)/2; END IF
   IF( ( PRESENT(S) ) )THEN; SS = SIZE(S); ELSE; SS = N; END IF
   IF( PRESENT(S) .AND. LSAME(LFACT,'F') .AND. LSAME(LEQUED,'Y') ) THEN
       MVS = MINVAL(S); ELSE; MVS = TINY(1.0_WP); ENDIF
   IF( PRESENT(FERR) )THEN; SFERR = SIZE(FERR); ELSE; SFERR = NRHS; END IF
   IF( PRESENT(BERR) )THEN; SBERR = SIZE(BERR); ELSE; SBERR = NRHS; END IF
   IF(PRESENT(UPLO))THEN; LUPLO = UPLO; ELSE; LUPLO='U'; END IF
!  .. TEST THE ARGUMENTS
   IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
   ELSE IF( SIZE(B, 1) /= N .OR. NRHS < 0 )THEN; LINFO = -2
   ELSE IF( SIZE(X, 1) /= N .OR. SIZE(X, 2) /= NRHS )THEN; LINFO = -3
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -4
   ELSE IF( SAF /= N*(N+1)/2 ) THEN; LINFO = -5
   ELSE IF( ( .NOT. ( LSAME(LFACT,'F') .OR. LSAME(LFACT,'N') .OR. &
                    LSAME(LFACT,'E') ) ) .OR. &
      ( LSAME(LFACT,'F') .AND. .NOT.PRESENT(AF) ) )THEN; LINFO = -6
   ELSE IF( .NOT.LSAME(LEQUED,'N') .AND. .NOT.LSAME(LEQUED,'Y') )THEN; LINFO = -7
   ELSE IF( SS /= N .OR. LSAME(LFACT,'F') .AND. LSAME(LEQUED,'Y') &
      .AND. MVS  <= 0.0_WP )THEN; LINFO = -8
   ELSE IF( SFERR /= NRHS )THEN; LINFO = -9
   ELSE IF( SBERR /= NRHS )THEN; LINFO = -10
   ELSE IF ( N > 0 )THEN
      IF( .NOT.PRESENT(AF) ) THEN; ALLOCATE( LAF(N*(N+1)/2), STAT=ISTAT )
      ELSE; LAF => AF; END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(S) )THEN; ALLOCATE( LS(N), STAT=ISTAT )
         ELSE; LS => S; END IF
      END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(FERR) )THEN; ALLOCATE( LFERR(NRHS), STAT=ISTAT )
         ELSE; LFERR => FERR; END IF
      END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(BERR) )THEN; ALLOCATE( LBERR(NRHS), STAT=ISTAT )
         ELSE; LBERR => BERR; END IF
      END IF
      IF( ISTAT == 0 ) ALLOCATE(WORK(3*N), IWORK(N), STAT=ISTAT )
      IF( ISTAT == 0 )THEN
         CALL PPSVX_F77( LFACT, LUPLO, N, NRHS, A, LAF, LEQUED, LS, &
                  B, N, X, N, LRCOND, LFERR, LBERR, WORK, IWORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(S) ) DEALLOCATE( LS, STAT=ISTAT1 )
      IF( .NOT.PRESENT(AF) ) DEALLOCATE( LAF, STAT=ISTAT1 )
      IF( .NOT.PRESENT(FERR) ) DEALLOCATE( LFERR, STAT=ISTAT1 )
      IF( .NOT.PRESENT(BERR) ) DEALLOCATE( LBERR, STAT=ISTAT1 )
      IF( PRESENT(RCOND) ) RCOND=LRCOND
      IF( PRESENT(EQUED) .AND. .NOT.LSAME(LFACT,'F') ) EQUED=LEQUED
      DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
END SUBROUTINE DPPSVX_F95
