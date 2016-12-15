SUBROUTINE ZHPGVX_F95( AP, BP, W, ITYPE, UPLO, Z, VL, VU, IL, IU, &
     &  M, IFAIL, ABSTOL, INFO )
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: LAMCH_F77 => DLAMCH
      USE F77_LAPACK, ONLY: HPGVX_F77 => LA_HPGVX
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!  .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(IN), OPTIONAL :: IL, IU, ITYPE
      INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
      REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAIL(:)
      COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
      COMPLEX(WP), INTENT(INOUT) :: AP(:), BP(:)
      REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_SPGVX and LA_HPGVX compute selected eigenvalues and, optionally, 
! the corresponding eigenvectors of generalized eigenvalue problems of 
! the form
!       A*z = lambda*B*z, A*B*z = lambda*z, and B*A*z = lambda*z,
! where A and B are real symmetric in the case of LA_SPGVX and complex
! Hermitian in the case of LA_HPGVX. In both cases B is positive
! definite. Eigenvalues and eigenvectors can be selected by specifying
! either a range of values or a range of indices for the desired 
! eigenvalues. Matrices A and B are stored in a packed format.
! 
! =========
! 
!       SUBROUTINE LA_SPGVX / LA_HPGVX( AP, BP, W, ITYPE= itype, &
!             UPLO= uplo, Z= z, VL= vl, VU= vu, IL= il, IU= iu, M= m, &
!             IFAIL= ifail, ABSTOL= abstol, INFO= info )
!         <type>(<wp>), INTENT(INOUT) :: AP(:), BP(:)
!         REAL(<wp>), INTENT(OUT) :: W(:)
!         INTEGER, INTENT(IN), OPTIONAL :: ITYPE
!         CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!         <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!         REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!         INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!         INTEGER, INTENT(OUT), OPTIONAL :: M
!         INTEGER, INTENT(OUT), OPTIONAL :: IFAIL(:)
!         REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!         INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!         <type> ::= REAL | COMPLEX
!         <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! AP       (input/output) REAL or COMPLEX array, shape (:) with 
!          size(AP) = n*(n + 1)/2, where n is the order of A and B.
!          On entry, the upper or lower triangle of matrix A in packed
! 	   storage. The elements are stored columnwise as follows:
! 	   if UPLO = 'U', AP(i +(j-1)*j/2) = A(i,j) for 1<=i<=j<=n;
! 	   if UPLO = 'L', AP(i +(j-1)*(2*n-j)/2) = A(i,j) for 1<=j<=i<=n.
!          On exit, the contents of AP are destroyed.
! BP       (input/output) REAL or COMPLEX array, shape (:) with 
!          size(BP) = size(AP).
!          On entry, the upper or lower triangle of matrix B in packed 
! 	   storage. The elements are stored columnwise as follows:
! 	   if UPLO = 'U', BP(i +(j-1)*j/2) = B(i,j) for 1<=i<=j<=n;
! 	   if UPLO = 'L', BP(i +(j-1)*(2*n-j)/2) = B(i,j) for 1<=j<=i<=n.
!          On exit, the triangular factor U or L of the Cholesky 
! 	   factorization B = U^T*U or B = L*L^T, in the same storage 
! 	   format as B.
! W        (output) REAL array, shape (:) with size(W) = n.
!          The eigenvalues in ascending order.
! ITYPE    Optional (input) INTEGER.
!          Specifies the problem type to be solved:
!             = 1: A*z = lambda*B*z
!             = 2: A*B*z = lambda*z
!             = 3: B*A*z = lambda*z
!          Default value: 1.
! UPLO     Optional (input) CHARACTER(LEN=1).
!             = 'U': Upper triangles of A and B are stored;
!             = 'L': Lower triangles of A and B are stored.
!          Default value: 'U'.
! Z        Optional (output) REAL or COMPLEX rectangular array, shape 
!          (:,:) with size(Z,1) = n and size(Z,2) = M.
!          The first M columns of Z contain the orthonormal eigenvectors
! 	   corresponding to the selected eigenvalues, with the i-th 
! 	   column of Z holding the eigenvector associated with the 
! 	   eigenvalue in W(i). The eigenvectors are normalized as 
! 	   follows:
!            if ITYPE = 1 or 2: Z^H * B * Z = I ,
!            if ITYPE = 3: Z^H * B^-1 * Z = I .
!          If an eigenvector fails to converge, then that column of Z
! 	   contains the latest approximation to the eigenvector and the
! 	   index of the eigenvector is returned in IFAIL.
! VL,VU    Optional (input) REAL.
!          The lower and upper bounds of the interval to be searched for
! 	   eigenvalues. VL < VU.
!          Default values: VL = -HUGE(<wp>) and VU = HUGE(<wp>), where 
! 	   <wp> ::= KIND(1.0) | KIND(1.0D0).
!          Note: Neither VL nor VU may be present if IL and/or IU is 
! 	   present.
! IL,IU    Optional (input) INTEGER.
!          The indices of the smallest and largest eigenvalues to be 
! 	   returned. The IL-th through IU-th eigenvalues will be found. 
! 	   1 <= IL <= IU <= size(A,1).
!          Default values: IL = 1 and IU = size(A,1).
!          Note: Neither IL nor IU may be present if VL and/or VU is 
! 	   present.
!          Note: All eigenvalues are calculated if none of the arguments
! 	   VL, VU, IL and IU are present.
! M        Optional (output) INTEGER.
!          The total number of eigenvalues found. 0 <= M <= size(A,1).
!          Note: If IL and IU are present then M = IU - IL + 1.
! IFAIL    Optional (output) INTEGER array, shape (:) with size(IFAIL) =
!          size(A,1).
!          If INFO = 0, the first M elements of IFAIL are zero.
!          If INFO > 0, then IFAIL contains the indices of the 
! 	   eigenvectors that failed to converge.
!          Note: If Z is present then IFAIL should also be present.
! ABSTOL   Optional (input) REAL.
!          The absolute error tolerance for the eigenvalues. An 
! 	   approximate eigenvalue is accepted as converged when it is
! 	   determined to lie in an interval [a,b] of width less than or
! 	   equal to
!              ABSTOL + EPSILON(1.0_<wp>) * max(|a|,|b|),
!          where <wp> is the working precision. If ABSTOL <= 0, then 
! 	   EPSILON(1.0_<wp>) * ||T||1 will be used in its place, where 
! 	   ||T||1 is the l1 norm of the tridiagonal matrix obtained by
! 	   reducing the generalized eigenvalue problem to tridiagonal
! 	   form. Eigenvalues will be computed most accurately when
! 	   ABSTOL is set to twice the underflow threshold 
! 	   2 * LA_LAMCH(1.0_<wp>, 'S'), not zero.
!          Default value: 0.0_<wp>.
!          Note: If this routine returns with 0 < INFO <= n, then some 
! 	   eigenvectors did not converge.
!          Try setting ABSTOL to 2 * LA_LAMCH(1.0_<wp>, 'S').
! INFO     Optional (output) INTEGER.
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: the algorithm failed to converge or matrix B is not 
! 	     positive definite:
!             <= n: the algorithm failed to converge; if INFO = i, then 
! 	       i eigenvectors failed to converge. Their indices are 
! 	       stored in array IFAIL.
!             > n: if INFO = n + i, for 1 <= i <= n, then the leading
! 	       minor of order i of B is not positive definite. The 
! 	       factorization of B could not be completed and no
!                eigenvalues or eigenvectors were computed.
!          If INFO is not present and an error occurs, then the program
! 	   is terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_HPGVX'
      CHARACTER(LEN=1) :: LJOBZ, LUPLO, LRANGE
!  .. LOCAL SCALARS ..
      INTEGER :: N, LINFO, LIL, LIU, LM, ISTAT, &
     &  SIFAIL, S1Z, S2Z, NN, LITYPE
      INTEGER, TARGET :: ISTAT1(1)
      COMPLEX(WP), TARGET :: LLZ(1,1)
      REAL(WP) :: LABSTOL, LVL, LVU
      COMPLEX(WP) :: WW
!  .. LOCAL ARRAYS ..
      INTEGER, POINTER :: IWORK(:), LIFAIL(:)
      COMPLEX(WP), POINTER :: WORK(:)
      REAL(WP), POINTER :: RWORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC HUGE, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; NN = SIZE(AP)
   WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW)
   IF( PRESENT(ITYPE) )THEN; LITYPE = ITYPE; ELSE; LITYPE = 1; END IF 
   IF( PRESENT(IFAIL) )THEN; SIFAIL = SIZE(IFAIL); ELSE; SIFAIL = N; END IF
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(VL) )THEN; LVL = VL; ELSE; LVL = -HUGE(LVL); ENDIF
   IF( PRESENT(VU) )THEN; LVU = VU; ELSE; LVU = HUGE(LVU); ENDIF
   IF( PRESENT(IL) )THEN; LIL = IL; ELSE; LIL = 1; ENDIF
   IF( PRESENT(IU) )THEN; LIU = IU; ELSE; LIU = N; ENDIF
   IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2)
    ELSE; S1Z = 1; S2Z = 1; ENDIF
!  .. TEST THE ARGUMENTS
   IF( PRESENT(VL) .OR. PRESENT(VU) )THEN; LRANGE = 'V'
   ELSE IF( PRESENT(IL) .OR. PRESENT(IU) )THEN; LRANGE = 'I'
   ELSE ; LRANGE = 'A'; END IF                                                 

   IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
   ELSE IF( SIZE(BP) /= SIZE(AP)  )THEN; LINFO = -2 
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -3
   ELSE IF (LITYPE <1 .OR. LITYPE >3) THEN; LINFO = -4
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -5
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= N .OR. S2Z /= N ) )THEN; LINFO = -6
   ELSE IF( LVU < LVL )THEN; LINFO = -7
   ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. &
     &  (PRESENT(IL) .OR. PRESENT(IU)) )THEN; LINFO = -8
   ELSE IF( LSAME(LRANGE,'I') .AND. ( LIU < MIN(N, LIL) .OR. LIU > N ))THEN; LINFO = -9
   ELSE IF( N < LIU )THEN; LINFO = -10
   ELSE IF( SIFAIL /= N .OR. PRESENT(IFAIL).AND..NOT.PRESENT(Z) )THEN; LINFO = -12
   ELSE IF( N > 0 )THEN
      IF( PRESENT(Z) ) THEN; LJOBZ = 'V'
        IF( PRESENT(IFAIL) )THEN; LIFAIL => IFAIL
        ELSE; ALLOCATE( LIFAIL(N), STAT=ISTAT ); END IF
        ELSE; LJOBZ = 'N'; LIFAIL => ISTAT1; ENDIF
!     .. DETERMINE THE WORKSPACE
       IF( ISTAT == 0 ) THEN
	 ALLOCATE(IWORK(MAX(1,5*N)), WORK(MAX(1,8*N)), RWORK(MAX(1,7*N)), STAT=ISTAT)
         IF( ISTAT == 0 )THEN
           IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL
           ELSE; LABSTOL = 2*LAMCH_F77('Safe minimum'); ENDIF
             IF (PRESENT (Z)) THEN
   	          CALL HPGVX_F77( LITYPE, LJOBZ, LRANGE, LUPLO, N, AP, BP, LVL, &
     &                   LVU, LIL, LIU, LABSTOL, LM, W, Z, S1Z, WORK, RWORK, IWORK, LIFAIL, &
     &                   LINFO )
             ELSE
	          CALL HPGVX_F77( LITYPE, LJOBZ, LRANGE, LUPLO, N, AP, BP, LVL, &
     &                   LVU, LIL, LIU, LABSTOL, LM, W, LLZ, S1Z, WORK, RWORK, IWORK, LIFAIL, &
     &                   LINFO )
             ENDIF
             IF( PRESENT(M) ) M = LM
           ELSE; LINFO = -100; END IF
           END IF
           IF( PRESENT(Z) .AND. .NOT.PRESENT(IFAIL) ) DEALLOCATE( LIFAIL, STAT=ISTAT1(1) )
	   DEALLOCATE(IWORK, WORK, RWORK, STAT=ISTAT1(1))
         END IF
         CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZHPGVX_F95
