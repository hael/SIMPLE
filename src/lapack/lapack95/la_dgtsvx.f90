SUBROUTINE DGTSVX_F95(DL, D, DU, B, X, DLF, DF, DUF, DU2, &
                           IPIV, FACT, TRANS, FERR, BERR, RCOND, INFO)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: LSAME, ERINFO
   USE F77_LAPACK, ONLY: GTSVX_F77 => LA_GTSVX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS, FACT
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
   REAL(WP), INTENT(OUT), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(IN) :: DL(:), D(:), DU(:), B(:,:)
   INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: IPIV(:)
   REAL(WP), INTENT(INOUT), OPTIONAL, TARGET :: DLF(:), DF(:), DUF(:), DU2(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: FERR(:), BERR(:)
   REAL(WP), INTENT(OUT) :: X(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_GTSVX computes the solution to a real or complex linear system 
! of equations of the form A*X = B, A^T*X = B or A^H*X = B, where A is a
! square tridiagonal matrix and X and B are rectangular matrices or
! vectors.
!     LA_GTSVX can also optionally estimate the condition number of A and 
! compute error bounds.
! 
! =========
! 
!       SUBROUTINE LA_GTSVX( DL, D, DU, B, X, DLF=dlf, DF=df, DUF=duf, &
!               DU2=du2, IPIV=ipiv, FACT=fact, TRANS=trans, FERR=ferr, &
!               BERR=berr, RCOND=rcond, INFO=info )
!           <type>(<wp>), INTENT(IN) :: DL(:), D(:), DU(:), <rhs>
!           <type>(<wp>), INTENT(OUT) :: <sol>
!           <type>(<wp>), INTENT(INOUT), OPTIONAL :: DLF(:), DF(:), &
!                                                   DUF(:), DU2(:)
!           INTEGER, INTENT(INOUT), OPTIONAL :: IPIV(:)
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: FACT, TRANS
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
! DL     (input) REAL or COMPLEX array, shape (:) with size(DL) = n-1.
!        The subdiagonal of A.
! D      (input) REAL or COMPLEX array, shape (:) with size(D) = n.
!        The diagonal of A.
! DU     (input) REAL or COMPLEX array, shape (:) with size(DU) = n-1.
!        The superdiagonal of A.
! B      (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(B,1) = n or shape (:) with size(B) = n.
!        The matrix B.
! X      (output) REAL or COMPLEX array, shape (:,:) with size(X,1) = n
!        and size(X,2) = size(B,2), or shape (:) with size(X) = n.
!        The solution matrix X .
! DLF    Optional (input or output) REAL or COMPLEX array, shape (:) with
!        size(DLF)= n-1.
!        If FACT = 'F' then DLF is an input argument that contains the 
!        multipliers that define the matrix L from the LU factorization
!        of A.
!        If FACT = 'N' then DLF is an output argument that contains the 
!        multipliers that define the matrix L from the LU factorization
!        of A.
! DF     Optional (input or output) REAL or COMPLEX array, shape (:) with
!        size(DF)= n.
!        If FACT = 'F' then DF is an input argument that contains the 
!        diagonal of the matrix U .
!        If FACT = 'N' then DF is an output argument that contains the 
!        diagonal of the matrix U .
! DUF    Optional (input or output) REAL or COMPLEX array, shape (:) with
!        size(DUF) = n-1.
!        If FACT = 'F' then DUF is an input argument that contains the
!        first superdiagonal of U.
!        If FACT = 'N' then DUF is an output argument that contains the
!        first superdiagonal of U.
! DU2    Optional (input or output) REAL or COMPLEX array, shape (:) with
!        size(DU2) = n-2.
!        If FACT = 'F', then DU2 is an input argument that contains the 
!        second superdiagonal of U.
!        If FACT = 'N', then DU2 is an output argument that contains the
!        second superdiagonal of U.
! IPIV   Optional (input or output) INTEGER array, shape (:) with 
!        size(IPIV) = n.
!        If FACT = 'F' then IPIV is an input argument that contains the
!        pivot indices from the LU factorization of A.
!        If FACT = 'N', then IPIV is an output argument that contains the
!        pivot indices from the LU factorization of A; row i of the 
!        matrix was interchanged with row IPIV(i). IPIV(i) will always
!        be either i or i+1; IPIV(i) = i indicates a row interchange was
!        not required.
! FACT   Optional (input) CHARACTER(LEN=1).
!        Specifies whether the factored form of A is supplied on entry.
!            = 'N': The matrix will be copied to DLF, DF and DUF and
! 	          factored.
!            = 'F': DLF, DF, DUF, DU2 and IPIV contain the factored form
! 	          of A.
!        Default value: 'N'.
! TRANS  Optional (input) CHARACTER(LEN=1).
!        Specifies the form of the system of equations:
!            = 'N': A*X = B (No transpose)
!            = 'T': A^T*X = B (Transpose)
!            = 'C': A^H*X = B (Conjugate transpose)
!        Default value: 'N'.
! FERR   Optional (output) REAL array of shape (:), with size(FERR) =
!        size(X,2), or REAL scalar.
!        The estimated forward error bound for each solution vector X(j)
!        (the j-th column of the solution matrix X). If XTRUE is the true
!        solution corresponding to X(j) , FERR(j) is an estimated upper
!        bound for the magnitude of the largest element in (X(j)-XTRUE)
!        divided by the magnitude of the largest element in X(j). The 
!        estimate is as reliable as the estimate for RCOND and is almost
!        always a slight overestimate of the true error.
! BERR   Optional (output) REAL array of shape (:), with size(BERR) = 
!        size(X,2), or REAL scalar.
!        The componentwise relative backward error of each solution 
!        vector X(j) (i.e.,the smallest relative change in any element of
!        A or B that makes X(j) an exact solution).
! RCOND  Optional (output) REAL.
!        The estimate of the reciprocal condition number of the matrix A.
!        If RCOND is less than the machine precision, the matrix is 
!        singular to working precision. This condition is indicated by
!        a return code of INFO > 0.
! INFO   Optional (output) INTEGER
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: if INFO = i, and i is 
!           <= n: U(i,i) = 0. The factorization has not been completed 
! 	        unless i = n. The factor U is singular, so the solution
! 		could not be computed.
!           = n+1: U is nonsingular, but RCOND is less than machine 
! 	        precision, meaning that the matrix is singular to 
! 		working precision. Nevertheless, the solution and
!                 error bounds are computed because the computed solution
! 		can be more accurate than the value of RCOND would 
! 		suggest.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!----------------------------------------------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GTSVX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LFACT, LTRANS
   INTEGER :: LINFO, NRHS, N, ISTAT, ISTAT1, SIPIV, SDLF, SDF, SDUF, SDU2, &
              SFERR, SBERR
   REAL(WP) :: LRCOND
!  .. LOCAL POINTERS ..
   INTEGER, POINTER :: IWORK(:), LPIV(:)
   REAL(WP),  POINTER :: LFERR(:), LBERR(:)
   REAL(WP),  POINTER :: WORK(:), LDLF(:), LDF(:), LDUF(:), LDU2(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0
   N = SIZE(D); NRHS = SIZE(B,2)
   IF( PRESENT(RCOND) ) RCOND = 1.0_WP
   IF( PRESENT(FACT) )THEN; LFACT = FACT; ELSE; LFACT='N'; END IF
   IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = N; END IF
   IF( PRESENT(DLF) )THEN; SDLF = SIZE(DLF); ELSE; SDLF = N-1; END IF
   IF( PRESENT(DF) )THEN; SDF = SIZE(DF); ELSE; SDF = N; END IF
   IF( PRESENT(DUF) )THEN; SDUF = SIZE(DUF); ELSE; SDUF = N-1; END IF
   IF( PRESENT(DU2) )THEN; SDU2 = SIZE(DU2); ELSE; SDU2 = N-2; END IF
   IF( PRESENT(FERR) )THEN; SFERR = SIZE(FERR); ELSE; SFERR = NRHS; END IF
   IF( PRESENT(BERR) )THEN; SBERR = SIZE(BERR); ELSE; SBERR = NRHS; END IF
   IF(PRESENT(TRANS))THEN; LTRANS = TRANS; ELSE; LTRANS='N'; END IF
!  PRINT *, LINFO, ISTAT, N, NRHS, LFACT, SIPIV, SDLF, SDF, SDUF, SDU2, SFERR, SBERR, LTRANS
!  .. TEST THE ARGUMENTS
   IF( SIZE( DL ) /= N-1 .AND. N/=0 ) THEN; LINFO = -1
   ELSE IF( N < 0 ) THEN; LINFO = -2
   ELSE IF( SIZE( DU ) /= N-1 .AND. N/=0) THEN; LINFO = -3
   ELSE IF( SIZE(B, 1) /= N .OR. NRHS < 0 )THEN; LINFO = -4
   ELSE IF( SIZE(X, 1) /= N .OR. SIZE(X, 2) /= NRHS )THEN; LINFO = -5
   ELSE IF( SDLF /= N-1 .AND. N/=0) THEN; LINFO = -6
   ELSE IF( SDF /= N ) THEN; LINFO = -7
   ELSE IF( SDUF /= N-1 .AND. N/=0 ) THEN; LINFO = -8
   ELSE IF( SDU2 /= N-2 .AND. N>1 ) THEN; LINFO = -9
   ELSE IF( SIPIV /= N )THEN; LINFO = -10
   ELSE IF( SFERR /= NRHS )THEN; LINFO = -13
   ELSE IF( SBERR /= NRHS )THEN; LINFO = -14
   ELSE IF( ( .NOT. ( LSAME(LFACT,'F') .OR. LSAME(LFACT,'N') ) ) .OR. &
       ( LSAME(LFACT,'F') .AND. .NOT.( PRESENT(DF) .AND. PRESENT(IPIV) ) ) )THEN
      LINFO = -11
   ELSE IF( .NOT.( LSAME(LTRANS,'N') .OR.  LSAME(LTRANS,'T') .OR. &
                  LSAME(LTRANS,'C') ) )THEN; LINFO = -12
   ELSE IF ( N > 0 )THEN
      IF( .NOT.PRESENT(DF) ) THEN
         ALLOCATE( LDLF(N-1),LDF(N),LDUF(N-1),LDU2(N-2), STAT=ISTAT )
      ELSE; LDLF => DLF; LDF => DF; LDUF => DUF; LDU2 => DU2; END IF
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
      IF( ISTAT == 0 ) ALLOCATE( WORK(3*N), IWORK(N), STAT=ISTAT )
      IF( ISTAT == 0 )THEN
         CALL GTSVX_F77( LFACT, LTRANS, N, NRHS, DL, D, DU, LDLF, LDF, LDUF, &
                         LDU2, LPIV, B, N, X, N, LRCOND, LFERR, LBERR, &
                         WORK, IWORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(DLF) ) DEALLOCATE( LDLF, LDF, LDUF, LDU2, STAT=ISTAT1 )
      IF( .NOT.PRESENT(IPIV) ) DEALLOCATE( LPIV, STAT=ISTAT1 )
      IF( .NOT.PRESENT(FERR) ) DEALLOCATE( LFERR, STAT=ISTAT1 )
      IF( .NOT.PRESENT(BERR) ) DEALLOCATE( LBERR, STAT=ISTAT1 )
      IF( PRESENT(RCOND) ) RCOND=LRCOND
      DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
END SUBROUTINE DGTSVX_F95
