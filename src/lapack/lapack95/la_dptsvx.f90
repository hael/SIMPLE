SUBROUTINE DPTSVX_F95(D, E, B, X, DF, EF, FACT, FERR, BERR, RCOND, INFO)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: LSAME, ERINFO
   USE F77_LAPACK, ONLY: PTSVX_F77 => LA_PTSVX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: FACT
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
   REAL(WP), INTENT(OUT), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(IN) :: D(:)
   REAL(WP), INTENT(IN) :: E(:), B(:,:)
   REAL(WP), INTENT(OUT) :: X(:,:)
   REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: DF(:)
   REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: EF(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: FERR(:), BERR(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_PTSVX computes the solution to a linear system of equations 
! A*X = B, where A has tridiagonal form and is real symmetric or complex
! Hermitian and, in either case, positive definite, and where X and B are
! rectangular matrices or vectors.
!    LA_PTSVX can also optionally estimate the condition number of A and
! compute error bounds.
! 
! =========
! 
!     SUBROUTINE LA_PTSVX( D, E, B, X, DF=df, EF=ef, FACT=fact, &
!                  FERR=ferr, BERR=berr, RCOND=rcond, INFO=info )
!          REAL(<wp>), INTENT(IN) :: D(:)
!          <type>(<wp>), INTENT(IN) :: E(:), <rhs>
!          <type>(<wp>), INTENT(OUT) :: <sol>
!          REAL(<wp>), INTENT(INOUT), OPTIONAL :: DF(:)
!          <type>(<wp>), INTENT(INOUT), OPTIONAL :: EF(:)
!          CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: FACT
!          REAL(<wp>), INTENT(OUT), OPTIONAL :: <err>, RCOND
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
! D       (input) REAL array, shape (:) with size(D) = n, where n is the 
!         order of A.
!         The diagonal of A.
! E       (input) REAL or COMPLEX array, shape (:) with size(E) = n-1.
!         The subdiagonal of A.
! B       (input) REAL or COMPLEX array, shape (:,:) with size(B,1) = n 
!         or shape (:) with size(B) = n.
!         The matrix B.
! X       (output) REAL or COMPLEX array, shape (:,:) with size(X,1) = n
!         and size(X,2) = size(B,2), or shape (:) with size(X) = n.
!         The solution matrix X .
! DF      Optional (input or output) REAL array, shape (:) with the same 
!         size as D.
!         If FACT = 'F', then DF is an input argument that contains the 
!         diagonal of D from the L*D*L^H factorization of A.
!         If FACT = 'N', then DF is an output argument that contains the
!         diagonal of D from the L*D*L^H factorization of A.
! EF      Optional (input or output) REAL or COMPLEX array, shape (:) with
!         the same size as E.
!         If FACT = 'F', then EF is an input argument that contains the 
!         subdiagonal of L from the L*D*L^H factorization of A.
!         If FACT = 'N', then EF is an output argument that contains the 
!         subdiagonal of L from the L*D*L^H factorization of A.
! FACT    Optional (input) CHARACTER(LEN=1).
!         Specifies whether the factored form of A has been supplied on
! 	  entry.
!           = 'N': The matrix A will be copied to DF and EF and factored.
!           = 'F': DF and EF contain the factored form of A.
!         Default value: 'N'.
! FERR    Optional (output) REAL array of shape (:), with 
!         size(FERR) = size(X,2), or REAL scalar.
!         The estimated forward error bound for each solution vector X(j)
! 	  (the j-th column of the solution matrix X). If XTRUE is the 
! 	  true solution corresponding to X(j), FERR(j) is an estimated 
!         upper bound for the magnitude of the largest element in 
! 	  (X(j)-XTRUE) divided by the magnitude of the largest element in
! 	  X(j).
! BERR    Optional (output) REAL array of shape (:), with size(BERR) = 
!         size(X,2), or REAL scalar.
!         The componentwise relative backward error of each solution 
!         vector X(j) (i.e., the smallest relative change in any element 
!         of A or B that makes X(j) an exact solution).
! RCOND   Optional (output) REAL.
!         The estimate of the reciprocal condition number of the matrix 
!         A. If RCOND is less than the machine precision, the matrix is 
!         singular to working precision. This condition is indicated by
!         a return code of INFO > 0.
! INFO    Optional (output) INTEGER
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, and i is
!             <= n: the leading minor of order i of A is not positive 
! 	         definite, so the factorization could not be completed
! 		 unless i = n, and the solution and error bounds could 
! 		 not be computed. RCOND = 0 is returned.
!             = n+1: L is nonsingular, but RCOND is less than machine 
! 	         precision, so the matrix is singular to working 
! 		 precision. Nevertheless, the solution and error
!                  bounds are computed because the computed solution can
!                  be more accurate than the value of RCOND would suggest.
!         If INFO is not present and an error occurs, then the program is
! 	  terminated with an error message.
!-----------------------------------------------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_PTSVX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LFACT
   INTEGER :: LINFO, NRHS, N, ISTAT, ISTAT1, SDF, SEF, SFERR, SBERR
   REAL(WP) :: LRCOND
!  .. LOCAL POINTERS ..
   REAL(WP),  POINTER :: LDF(:), LFERR(:), LBERR(:)
   REAL(WP),  POINTER :: WORK(:), LEF(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0
   N = SIZE(D); NRHS = SIZE(B,2)
   IF( PRESENT(RCOND) ) RCOND = 1.0_WP
   IF( PRESENT(FACT) )THEN; LFACT = FACT; ELSE; LFACT='N'; END IF
   IF( PRESENT(DF) )THEN; SDF = SIZE(DF); ELSE; SDF = N; END IF
   IF( PRESENT(EF) )THEN; SEF = SIZE(EF); ELSE; SEF = N-1; END IF
   IF( PRESENT(FERR) )THEN; SFERR = SIZE(FERR); ELSE; SFERR = NRHS; END IF
   IF( PRESENT(BERR) )THEN; SBERR = SIZE(BERR); ELSE; SBERR = NRHS; END IF
!  .. TEST THE ARGUMENTS
   IF( N < 0 ) THEN; LINFO = -1
   ELSE IF( SIZE( E ) /= N-1 .AND. N /= 0 ) THEN; LINFO = -2
   ELSE IF( SIZE(B, 1) /= N .OR. NRHS < 0 )THEN; LINFO = -3
   ELSE IF( SIZE(X, 1) /= N .OR. SIZE(X, 2) /= NRHS )THEN; LINFO = -4
   ELSE IF( SDF /= N ) THEN; LINFO = -5
   ELSE IF( .NOT.( PRESENT(DF).AND.PRESENT(EF) ) &
       .AND.( PRESENT(DF).OR.PRESENT(EF) ) )THEN; LINFO = -5
   ELSE IF( SEF /= N-1 .AND. N>0 ) THEN; LINFO = -6
   ELSE IF( ( .NOT.LSAME(LFACT,'F') .AND. .NOT.LSAME(LFACT,'N') ) .OR. &
            ( LSAME(LFACT,'F') .AND. .NOT.PRESENT(DF) ) )THEN; LINFO = -7
   ELSE IF( SFERR /= NRHS )THEN; LINFO = -8
   ELSE IF( SBERR /= NRHS )THEN; LINFO = -9
   ELSE IF ( N > 0 )THEN
      IF( .NOT.PRESENT(DF) ) THEN; ALLOCATE( LDF(N), LEF(N-1), STAT=ISTAT )
      ELSE; LDF => DF; LEF => EF; END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(FERR) )THEN; ALLOCATE( LFERR(NRHS), STAT=ISTAT )
         ELSE; LFERR => FERR; END IF
      END IF
      IF( ISTAT == 0 )THEN
         IF( .NOT.PRESENT(BERR) )THEN; ALLOCATE( LBERR(NRHS), STAT=ISTAT )
         ELSE; LBERR => BERR; END IF
      END IF
      IF( ISTAT == 0 ) ALLOCATE(WORK(2*N), STAT=ISTAT )
      IF( ISTAT == 0 )THEN
         CALL PTSVX_F77( LFACT, N, NRHS, D, E, LDF, LEF, B, N, X, N, LRCOND, &
                         LFERR, LBERR, WORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(DF) ) DEALLOCATE( LDF, LEF, STAT=ISTAT1 )
      IF( .NOT.PRESENT(FERR) ) DEALLOCATE( LFERR, STAT=ISTAT1 )
      IF( .NOT.PRESENT(BERR) ) DEALLOCATE( LBERR, STAT=ISTAT1 )
      IF( PRESENT(RCOND) ) RCOND=LRCOND
      DEALLOCATE( WORK, STAT=ISTAT1 )
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
END SUBROUTINE DPTSVX_F95
