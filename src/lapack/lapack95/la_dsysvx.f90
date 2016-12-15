SUBROUTINE DSYSVX_F95(A, B, X, UPLO, AF, IPIV, FACT, &
                      FERR, BERR, RCOND, INFO)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: LSAME, ERINFO
   USE F77_LAPACK, ONLY: SYSVX_F77 => LA_SYSVX, ILAENV_F77 => ILAENV
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO, FACT
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
   REAL(WP), INTENT(OUT), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(IN) :: A(:,:), B(:,:)
   REAL(WP), INTENT(OUT) :: X(:,:)
   INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: IPIV(:)
   REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: AF(:,:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: FERR(:), BERR(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_SYSVX computes the solution to a linear system of equations 
! A*X = B, where A is a real or complex symmetric matrix and X and B are
! rectangular matrices or vectors.
!      LA_HESVX computes the solution to a linear system of equations 
! A*X = B, where A is a complex Hermitian matrix and X and B are
! rectangular matrices or vectors.
!      LA_SYSVX and LA_HESVX can also optionally estimate the condition 
! number of A and compute error bounds.
! 
! =========
! 
!         SUBROUTINE LA_SYSVX / LA HESVX( A, B, X, UPLO=uplo, AF=af, &
!                        IPIV=ipiv, FACT=fact, FERR=ferr, BERR=berr, &
!                        RCOND=rcond, INFO=info )
!              <type>(<wp>), INTENT(IN) :: A(:,:), <rhs>
!              <type>(<wp>), INTENT(OUT) :: <sol>
!              CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!              <type>(<wp>), INTENT(INOUT), OPTIONAL :: AF(:,:)
!              INTEGER, INTENT(INOUT), OPTIONAL :: IPIV(:)
!              CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: FACT
!              REAL(<wp>), INTENT(OUT), OPTIONAL :: <err>, RCOND
!              INTEGER, INTENT(OUT), OPTIONAL :: INFO
!         where
!              <type> ::= REAL | COMPLEX
!              <wp>   ::= KIND(1.0) | KIND(1.0D0)
!              <rhs>  ::= B(:,:) | B(:)
!              <sol>  ::= X(:,:) | X(:)
!              <err>  ::= FERR(:), BERR(:) | FERR, BERR
! 
! Arguments
! =========
! 
! A       (input) REAL or COMPLEX square array, shape (:,:).
!         The symmetric or Hermitian matrix A.
!         If UPLO = 'U', the upper triangular part of A contains the 
! 	  upper triangular part of the matrix A, and the strictly lower
! 	  triangular part of A is not referenced. If UPLO = 'L', the 
! 	  lower triangular part of A contains the lower triangular part
! 	  of the matrix A, and the strictly upper triangular part of A is
! 	  not referenced.
! B       (input) REAL or COMPLEX array, shape (:,:) with size(B,1) = 
!         size(A,1) or shape (:) with size(B) = size(A,1).
!         The matrix B.
! X       (output) REAL or COMPLEX array, shape (:,:) with size(X,1) = 
!         size(A,1) and size(X,2) = size(B,2), or shape (:) with size(X)
! 	  = size(A,1).
!         The solution matrix X.
! UPLO    Optional (input) CHARACTER(LEN=1).
!            = 'U': Upper triangle of A is stored;
!            = 'L': Lower triangle of A is stored.
!         Default value: 'U'.
! AF      Optional (input or output) REAL or COMPLEX array, shape (:,:)
!         with the same size as A.
!         If FACT = 'F', then AF is an input argument that contains the 
! 	  block diagonal matrix D and the multipliers used to obtain the
! 	  factor L or U from the factorization of A, returned by a
!         previous call to LA_SYSVX or LA_HESVX.
!         If FACT = 'N', then AF is an output argument that contains the
! 	  block diagonal matrix D and the multipliers used to obtain the
! 	  factor L or U from the factorization of A.
! IPIV    Optional (input or output) INTEGER array, shape (:) with 
!         size(IPIV) = size(A,1).
!         If FACT = 'F', then IPIV is an input argument that contains 
! 	  details of the row and column interchanges and the block 
! 	  structure of D.
!         If IPIV(k) > 0 , then rows and columns k and IPIV(k) were 
! 	  interchanged and D(k,k) is a 1 by 1 diagonal block.
!         If IPIV(k) < 0 , then there are two cases:
!           1. If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
! 	     columns k-1 and -IPIV(k) were interchanged and 
! 	     D(k-1:k,k-1:k) is a 2 by 2 diagonal block.
!           2. If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0, then rows and
! 	     columns k+1 and -IPIV(k) were interchanged and 
! 	     D(k:k+1,k:k+1) is a 2 by 2 diagonal block.
!         If FACT = 'N', then IPIV is an output argument that contains 
! 	  details of the row and column interchanges and the block 
! 	  structure of D; as described above.
! FACT    Optional (input) CHARACTER(LEN=1).
!         Specifies whether the factored form of the matrix A has been 
! 	  supplied on entry.
!           = 'N': The matrix A will be copied to AF and factored.
!           = 'F': AF and IPIV contain the factored form of A.
!         Default value: 'N'.
! FERR    Optional (output) REAL array of shape (:), with 
!         size(FERR) = size(X,2), or REAL scalar.
!         The estimated forward error bound for each solution vector 
! 	  X(j) (the j-th column of the solution matrix X). If XTRUE is
! 	  the true solution corresponding to X(j), FERR(j) is an
!         estimated upper bound for the magnitude of the largest element
! 	  in (X(j)-XTRUE) divided by the magnitude of the largest element
! 	  in X(j). The estimate is as reliable as the estimate for RCOND, 
! 	  and is almost always a slight overestimate of the true error.
! BERR    Optional (output) REAL array of shape (:), with size(BERR) = 
!         size(X,2), or REAL scalar.
!         The componentwise relative backward error of each solution 
! 	  vector X(j) (i.e., the smallest relative change in any element
! 	  of A or B that makes X(j) an exact solution).
! RCOND   Optional (output) REAL
!         The estimate of the reciprocal condition number of A. If RCOND 
! 	  is less than the machine
!         precision, the matrix is singular to working precision. This 
! 	  condition is indicated by a return code of INFO > 0.
! INFO    (output) INTEGER
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, and i is
!             <= n: D(i,i) = 0. The factorization has been completed, but
! 	         the block diagonal matrix D is singular, so the 
! 		 solution could not be computed.
!             = n+1: D is nonsingular, but RCOND is less than machine 
! 	         precision, so the matrix is singular to working 
! 		 precision. Nevertheless, the solution and error bounds
! 		 are computed because the computed solution can be more
! 		 accurate than the value of RCOND would suggest.
!             n is the order of A.
!         If INFO is not present and an error occurs, then the program is
! 	  terminated with an error message.
!------------------------------------------------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SYSVX'
   CHARACTER(LEN=6), PARAMETER :: BSNAME = 'DSYTRF'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LFACT, LUPLO
   INTEGER :: LINFO, NRHS, N, NB, LWORK, ISTAT, ISTAT1, SIPIV, S1AF, S2AF, SFERR, SBERR
   REAL(WP) :: LRCOND
!  .. LOCAL POINTERS ..
   INTEGER, POINTER :: IWORK(:), LPIV(:)
   REAL(WP),  POINTER :: LFERR(:), LBERR(:)
   REAL(WP),  POINTER :: WORK(:), LAF(:, :)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, SIZE, MAX
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A, 1); NRHS = SIZE(B, 2)
   IF( PRESENT(RCOND) ) RCOND = 1.0_WP
   IF( PRESENT(FACT) )THEN; LFACT = FACT; ELSE; LFACT='N'; END IF
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = N; END IF
   IF( PRESENT(AF) )THEN; S1AF = SIZE(AF,1); S2AF = SIZE(AF,2)
   ELSE; S1AF = N; S2AF = N; END IF
   IF( PRESENT(FERR) )THEN; SFERR = SIZE(FERR); ELSE; SFERR = NRHS; END IF
   IF( PRESENT(BERR) )THEN; SBERR = SIZE(BERR); ELSE; SBERR = NRHS; END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE(A, 2) /= N .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE(B, 1) /= N .OR. NRHS < 0 )THEN; LINFO = -2
   ELSE IF( SIZE(X, 1) /= N .OR. SIZE(X, 2) /= NRHS )THEN; LINFO = -3
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -4
   ELSE IF( S1AF /= N .OR. S2AF /= N ) THEN; LINFO = -5
   ELSE IF( SIPIV /= N )THEN; LINFO = -6
   ELSE IF( ( .NOT. LSAME(LFACT,'F') .AND. .NOT. LSAME(LFACT,'N') ) .OR. &
     ( LSAME(LFACT,'F') .AND. .NOT.( PRESENT(AF) .AND. PRESENT(IPIV) ) ) )THEN; LINFO = -7
   ELSE IF( SFERR /= NRHS )THEN; LINFO = -8
   ELSE IF( SBERR /= NRHS )THEN; LINFO = -9
   ELSE IF ( N > 0 )THEN
      IF( .NOT.PRESENT(AF) ) THEN; ALLOCATE( LAF(N,N), STAT=ISTAT )
      ELSE; LAF => AF; END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(IPIV) )THEN; ALLOCATE( LPIV(N), STAT=ISTAT )
         ELSE; LPIV => IPIV; END IF
      END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(FERR) )THEN; ALLOCATE( LFERR(NRHS), STAT=ISTAT )
         ELSE; LFERR => FERR; END IF
      END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(BERR) )THEN; ALLOCATE( LBERR(NRHS), STAT=ISTAT )
         ELSE; LBERR => BERR; END IF
      END IF
      IF( ISTAT == 0 )THEN
         NB = ILAENV_F77( 1, BSNAME, LUPLO, N, -1, -1, -1 )
         IF( NB <= 1 .OR. NB >= N ) NB = 1; LWORK = MAX(1,3*N,N*NB)
         ALLOCATE(WORK(LWORK), IWORK(N), STAT=ISTAT)
         IF( ISTAT /= 0 )THEN
            DEALLOCATE(WORK, IWORK, STAT=ISTAT1); LWORK = MAX(1,3*N)
            ALLOCATE(WORK(LWORK), IWORK(N), STAT=ISTAT)
            IF( ISTAT /= 0 ) THEN; LINFO = - 100
            ELSE; CALL ERINFO( -200, SRNAME, LINFO ); ENDIF
         ENDIF
      END IF
      IF( ISTAT == 0 )THEN
!        .. CALL LAPACK77 ROUTINE
         CALL SYSVX_F77( LFACT, LUPLO, N, NRHS, A, N, LAF, N, LPIV, B, N, X, N, &
                         LRCOND, LFERR, LBERR, WORK, LWORK, IWORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(AF) ) DEALLOCATE( LAF, STAT=ISTAT1 )
      IF( .NOT.PRESENT(IPIV) ) DEALLOCATE( LPIV, STAT=ISTAT1 )
      IF( .NOT.PRESENT(FERR) ) DEALLOCATE( LFERR, STAT=ISTAT1 )
      IF( .NOT.PRESENT(BERR) ) DEALLOCATE( LBERR, STAT=ISTAT1 )
      IF( PRESENT(RCOND) ) RCOND=LRCOND
      DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
END SUBROUTINE DSYSVX_F95
