SUBROUTINE DSYGVX_F95( A, B, W, ITYPE, JOBZ, UPLO, VL, VU, IL, IU, &
     &  M, IFAIL, ABSTOL, INFO ) 
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: LAMCH_F77 => DLAMCH
      USE F77_LAPACK, ONLY: SYGVX_F77 => LA_SYGVX
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!  .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(IN), OPTIONAL :: IL, IU, ITYPE
      INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
      REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAIL(:)
      REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
      REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_SYGVX and LA_HEGVX compute selected eigenvalues and, optionally,
! the corresponding eigenvectors of generalized eigenvalue problems of 
! the form 
!      A*z = lambda*B*z, A*B*z = lambda*z,  and  B*A*z = lambda*z,
! where A and B are real symmetric in the case of LA_SYGVX and complex
! Hermitian in the case of LA_HEGVX. In both cases B is positive 
! definite. Eigenvalues and eigenvectors can be selected by specifying
! either a range of values or a range of indices for the desired 
! eigenvalues.
! 
! =========
! 
!        SUBROUTINE LA_SYGVX / LA_HEGVX (A, B, W, ITYPE= itype, &
!             JOBZ= jobz, UPLO= uplo, VL= vl, VU= vu, IL= il, &
!             IU= iu, M= m, IFAIL= ifail, ABSTOL= abstol, INFO= info )
!          <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
!          REAL(<wp>), INTENT(OUT) :: W(:)
!          INTEGER, INTENT(IN), OPTIONAL :: ITYPE
!          CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!          REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!          INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!          INTEGER, INTENT(OUT), OPTIONAL :: M
!          INTEGER, INTENT(OUT), OPTIONAL :: IFAIL(:)
!          REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX square array, shape (:,:).
!        On entry, the matrix A.
!        If UPLO = 'U', the upper triangular part of A contains the 
!        upper triangular part of matrix A. If UPLO = 'L', the lower
!        triangular part of A contains the lower triangular part of 
!        matrix A.
!        On exit, if JOBZ = 'V', the first M columns of A contain the
!        orthonormal eigenvectors corresponding to the selected 
!        eigenvalues, with the i-th column of A holding the eigenvector
!        associated with the eigenvalue in W(i).
!        The eigenvectors are normalized as follows:
!          if ITYPE = 1 or 2: Z^H * B * Z = I ,
!          if ITYPE = 3: Z^H * B^-1 * Z = I .
!        If an eigenvector fails to converge, then that column of A 
!        contains the latest approximation to the eigenvector and the 
!        index of the eigenvector is returned in IFAIL.
!        If JOBZ = 'N', then the upper triangle (if UPLO = 'U') or the 
!        lower triangle (if UPLO = 'L') of A, including the diagonal, is
!        destroyed.
! B      (input/output) REAL or COMPLEX square array, shape (:,:) with 
!        size(B,1) = size(A,1).
!        On entry, the matrix B.
!        If UPLO = 'U', the upper triangular part of B contains the 
!        upper triangular part of matrix B. If UPLO = 'L', the lower
!        triangular part of B contains the lower triangular part of 
!        matrix B.
!        On exit, the part of B containing the matrix is overwritten by 
!        the triangular factor U or L of the Cholesky factorization 
!              B = U^H*U or B = L*L^H.
! W      (output) REAL array, shape (:) with size(W) = size(A,1).
!        The first M elements contain the selected eigenvalues in 
!        ascending order.
! ITYPE  Optional (input) INTEGER.
!        Specifies the problem type to be solved:
!           = 1: A*z = lambda*B*z
!           = 2: A*B*z = lambda*z
!           = 3: B*A*z = lambda*z
!        Default value: 1.
! JOBZ   Optional (input) CHARACTER(LEN=1).
!           = 'N': Computes eigenvalues only;
!           = 'V': Computes eigenvalues and eigenvectors.
!        Default value: 'N'.
! UPLO   Optional (input) CHARACTER(LEN=1).
!           = 'U': Upper triangles of A and B are stored;
!           = 'L': Lower triangles of A and B are stored.
!        Default value: 'U'.
! VL,VU  Optional (input) REAL.
!        The lower and upper bounds of the interval to be searched for
!        eigenvalues. VL < VU.
!        Default values: VL = -HUGE(<wp>) and VU = HUGE(<wp>), where 
!        <wp> ::= KIND(1.0) | KIND(1.0D0).
!        Note: Neither VL nor VU may be present if IL and/or IU is
!        present.
! IL,IU  Optional (input) INTEGER.
!        The indices of the smallest and largest eigenvalues to be 
!        returned. The IL-th through IU-th eigenvalues will be found. 
!        1<=IL<=IU<=size(A,1).
!        Default values: IL = 1 and IU = size(A,1).
!        Note: Neither IL nor IU may be present if VL and/or VU is
!        present.
!        Note: All eigenvalues are calculated if none of the arguments
!        VL, VU, IL and IU are present.
! M      Optional (output) INTEGER.
!        The total number of eigenvalues found. 0 <= M <= size(A,1).
!        Note: If IL and IU are present then M = IU - IL + 1.
! IFAIL  Optional (output) INTEGER array, shape (:) with size(IFAIL) = 
!        size(A,1).
!        If INFO = 0, the first M elements of IFAIL are zero.
!        If INFO > 0, then IFAIL contains the indices of the 
!        eigenvectors that failed to converge.
!        Note: IFAIL should be present if JOBZ = 'V'.
! ABSTOL Optional (input) REAL.
!        The absolute error tolerance for the eigenvalues. An approximate
!        eigenvalue is accepted as converged when it is determined to lie
!        in an interval [a,b] of width less than or equal to
!             ABSTOL + EPSILON(1.0_<wp>) * max(| a |, | b |),
!        where <wp> is the working precision. If ABSTOL <= 0, then
!        EPSILON(1.0_<wp>)* ||T||1 will be used in its place, where 
!        ||T||1 is the l1 norm of the tridiagonal matrix obtained by 
!        reducing the generalized eigenvalue problem to tridiagonal form. 
!        Eigenvalues will be computed most accurately when ABSTOL is set
!        to twice the underflow threshold 2 * LA_LAMCH(1.0_<wp>, 'S'),
!        not zero.
!        Default value: 0.0_<wp>.
!        Note: If this routine returns with 0 < INFO <= n, then some 
!        eigenvectors did not converge.
!        Try setting ABSTOL to 2 * LA_LAMCH(1.0_<wp>, 'S').
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: the algorithm failed to converge or matrix B is not
!        positive definite:
!           <= n: the algorithm failed to converge; if INFO = i, then i
! 	      eigenvectors failed to converge. Their indices are stored
! 	      in array IFAIL.
!           > n: if INFO = n+i, for 1 <= i <= n, then the leading minor 
! 	      of order i of B is not positive definite. The 
! 	      factorization of B could not be completed and no 
! 	      eigenvalues or eigenvectors were computed.
! 	  n is the order of A.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!-----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SYGVX'
      CHARACTER(LEN=6), PARAMETER :: BSNAME = 'DSYTRD'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBZ, LUPLO, LRANGE
      INTEGER :: N, LINFO, LDA, LDZ, LZ, LIL, LIU, LM, LWORK, ISTAT, &
      &     SIFAIL, LDB, LITYPE
      INTEGER, TARGET :: ISTAT1(1)
      REAL(WP) :: LABSTOL, LVL, LVU
!  .. LOCAL ARRAYS ..
      INTEGER, POINTER :: IWORK(:), LIFAIL(:)
      REAL(WP), POINTER :: WORK(:), Z(:,:)
      REAL(WP) :: WORKMIN(1)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC HUGE, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   IF (PRESENT(M)) M = 0
   N = SIZE(A,1); LDA = MAX(1,N); LDB=MAX(1,SIZE(B,1)); LINFO = 0; ISTAT = 0
   IF( PRESENT(ITYPE) )THEN; LITYPE = ITYPE; ELSE; LITYPE = 1; END IF 
   IF( PRESENT(IFAIL) )THEN; SIFAIL = SIZE(IFAIL);ELSE; SIFAIL = N; END IF
   IF (PRESENT (JOBZ)) THEN; LJOBZ=JOBZ; ELSE; LJOBZ = 'N'; ENDIF
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(VL) )THEN; LVL = VL; ELSE; LVL = -HUGE(LVL); ENDIF
   IF( PRESENT(VU) )THEN;  LVU = VU; ELSE; LVU = HUGE(LVU); ENDIF
   IF( PRESENT(IL) )THEN; LIL = IL; ELSE; LIL = 1; ENDIF
   IF( PRESENT(IU) )THEN; LIU = IU; ELSE; LIU = N; ENDIF
!  .. TEST THE ARGUMENTS
   IF( PRESENT(VL) .OR. PRESENT(VU) )THEN ; LRANGE = 'V'; LM=N
   ELSE IF( PRESENT(IL) .OR. PRESENT(IU) )THEN ; LRANGE = 'I'; LM=LIU-LIL+1
   ELSE ; LRANGE = 'A'; LM=N; END IF
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN; LINFO = -1
   ELSE IF (SIZE (B, 2) /= N ) THEN; LINFO = -2
   ELSE IF( SIZE( W ) /= N )THEN;  LINFO = -3
   ELSE IF( LITYPE < 1 .OR. LITYPE > 3 )THEN; LINFO = -4
   ELSE IF( .NOT.LSAME(LJOBZ,'V') .AND. .NOT.LSAME(LJOBZ,'N') )THEN; LINFO = -5
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -6
   ELSE IF( LVU < LVL )THEN ; LINFO = -7
   ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. &
            (PRESENT(IL) .OR. PRESENT(IU)) )THEN; LINFO = -8
   ELSE IF( LSAME(LRANGE, 'I') .AND. ( LIU < MIN( N, LIL ) .OR. LIU>N))THEN; LINFO = -9
   ELSE IF( N < LIU )THEN; LINFO = -10
   ELSE IF( SIFAIL /= N )THEN; LINFO = -12
   ELSE IF( N > 0 )THEN
     IF(LSAME(LJOBZ, 'V')) THEN
         LDZ = MAX(1,N); LZ=LM
     ELSE
         LDZ = 1; LZ=1
     ENDIF
       IF( PRESENT(IFAIL) )THEN; LIFAIL => IFAIL
       ELSE; ALLOCATE( LIFAIL(N), STAT=ISTAT )
         IF (ISTAT /= 0) THEN; LINFO = -100; GOTO 200; ENDIF         
         END IF
!     .. DETERMINE THE WORKSPACE
         ALLOCATE(IWORK(5*N), STAT=ISTAT)
         IF (ISTAT /= 0) THEN; LINFO = -100
           GOTO 300
         ENDIF
         ALLOCATE(Z(LDZ, LZ), STAT=ISTAT)
	 IF (ISTAT /= 0) THEN; LINFO = -100
	    GOTO 400
	 ENDIF
         
	 LWORK = -1
             CALL SYGVX_F77( LITYPE, LJOBZ, LRANGE, LUPLO, N, A, LDA, B, &
&                  LDB, LVL, LVU, LIL, LIU, LABSTOL, LM, W, Z, LDZ, WORKMIN, &
&                  LWORK, IWORK, LIFAIL, LINFO )
! NEXT LINE SHOULD BE LWORK = WORKMIN(1) 	   
           LWORK = 2*WORKMIN(1)

           ALLOCATE (WORK(LWORK), STAT = ISTAT)
           IF( ISTAT /= 0 )THEN; LINFO = -100
             GOTO 500
           ENDIF
             IF( LINFO == 0 )THEN
               IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL
               ELSE; LABSTOL = 2*LAMCH_F77('Safe minimum');  ENDIF
!     .. CALL LAPACK77 ROUTINE
                 CALL SYGVX_F77( LITYPE, LJOBZ, LRANGE, LUPLO, N, A, LDA, B, &
&                  LDB, LVL, LVU, LIL, LIU, LABSTOL, LM, W, Z, LDZ, WORK, &
                   LWORK, IWORK, LIFAIL, LINFO )

                 IF( PRESENT(M) ) M = LM
		 IF (LSAME(LJOBZ,'V'))  A(1:LDZ, 1:LM)=Z(1:LDZ, 1:LM)
              END IF
         DEALLOCATE (WORK, STAT = ISTAT1(1))
500      DEALLOCATE(Z)	      
400      DEALLOCATE (IWORK, STAT = ISTAT)
300      IF (.NOT.PRESENT(IFAIL)) DEALLOCATE(LIFAIL, STAT=ISTAT1(1)); END IF
200       CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSYGVX_F95
