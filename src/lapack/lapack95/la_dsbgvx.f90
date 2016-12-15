SUBROUTINE DSBGVX_F95( AB, BB,  W, UPLO, Z, VL, VU, IL, IU, &
      & 	M, IFAIL, Q,  ABSTOL, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP =>  DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: LAMCH_F77 => DLAMCH
      USE F77_LAPACK, ONLY: SBGVX_F77 => LA_SBGVX
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
      INTEGER, INTENT(IN), OPTIONAL :: IL, IU
      INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
      REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAIL(:)
      REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:), Q(:,:)
      REAL(WP), INTENT(INOUT) :: AB(:,:), BB(:,:)
      REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_SBGVX and LA_HBGVX compute selected eigenvalues and, optionally,
! the corresponding eigenvectors of the generalized eigenvalue problem
!                      A*z = lambda*B*z,
! where A and B are real symmetric in the case of LA_SBGVX and complex
! Hermitian in the case of LA_HBGVX. In both cases B is positive 
! definite. Matrices A and B are stored in a band format. Eigenvalues 
! and eigenvectors can be selected by specifying either a range of
! values or a range of indices for the desired eigenvalues.
! 
! =========
! 
!           SUBROUTINE LA_SBGVX / LA_HBGVX( AB, BB, W, UPLO=uplo, Z=z, &
!                   VL=vl, VU=vu, IL=il, IU=iu, M=m, IFAIL=ifail, Q=q, &
!                   ABSTOL=abstol, INFO=info )
!              <type>(<wp>), INTENT(INOUT) :: AB(:,:), BB(:,:)
!              REAL(<wp>), INTENT(OUT) :: W(:)
!              CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!              <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!              REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!              INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!              INTEGER, INTENT(OUT), OPTIONAL :: M
!              INTEGER, INTENT(OUT), OPTIONAL :: IFAIL(:)
!              <type>(<wp>), INTENT(OUT), OPTIONAL :: Q(:,:)
!              REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!              INTEGER, INTENT(OUT), OPTIONAL :: INFO
!           where 
!              <type> ::= REAL | COMPLEX
!              <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! AB       (input/output) REAL or COMPLEX array, shape (:,:) with 
!          size(AB,1) = ka + 1 and size(AB,2) = n, where ka is the number
! 	 of subdiagonals or superdiagonals in the band of A and n is 
! 	 the order of A and B.
!          On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')
! 	 triangle of A in band storage. The ka + 1 diagonals of A are 
! 	 stored in the rows of AB so that the j-th column of A is 
! 	 stored in the j-th column of AB as follows:
! 	 if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j,
! 	                                            1<=j<=n
! 	 if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka),
!                  	                            1<=j<=n.
!        On exit, the contents of AB are destroyed.
! BB     (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(BB,1) = kb + 1 and size(BB,2) = n, where kb is the number
! 	 of subdiagonals or superdiagonals in the band of B.
!        On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') 
! 	 triangle of matrix B in band storage. The kb + 1 diagonals of
! 	 B are stored in the rows of BB so that the j-th column of B
!        is stored in the j-th column of BB as follows:
! 	 if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j,
! 	                                            1<=j<=n
! 	 if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb),
!                                                   1<=j<=n.
!        On exit, the factor S from the split Cholesky factorization 
! 	              B = S^H*S.
! W      (output) REAL array, shape (:) with size(W) = n.
!        The first M elements contain the selected eigenvalues in 
! 	 ascending order.
! UPLO   Optional (input) CHARACTER(LEN=1).
!            = 'U': Upper triangles of A and B are stored;
!            = 'L': Lower triangles of A and B are stored.
!        Default value: 'U'.
! Z      Optional (output) REAL or COMPLEX square array, shape (:,:) 
!        with size(Z,1) = n.
!        The first M columns of Z contain the orthonormal eigenvectors
! 	 corresponding to the selected eigenvalues, with the i-th 
!        column of Z containing the eigenvector associated with the
!        eigenvalue in W(i). The eigenvectors are normalized so that 
! 	 Z^H*B*Z = I . If an eigenvector fails to converge, then that 
! 	 column of Z contains the latest approximation to the 
! 	 eigenvector and the index of the eigenvector is returned in 
! 	 IFAIL.
! VL,VU  Optional (input) REAL.
!        The lower and upper bounds of the interval to be searched for 
! 	 eigenvalues. VL < VU.
!        Default values: VL = -HUGE(<wp>) and VU = HUGE(<wp>), where 
! 	 <wp> ::= KIND(1.0) | KIND(1.0D0).
!        Note: Neither VL nor VU may be present if IL and/or IU is 
! 	 present.
! IL,IU  Optional (input) INTEGER.
!        The indices of the smallest and largest eigenvalues to be 
! 	 returned. The IL-th through IU-th eigenvalues will be found.
! 	 1 <= IL <= IU <= size(A,1).
!        Default values: IL = 1 and IU = size(A,1).
!        Note: Neither IL nor IU may be present if VL and/or VU is 
! 	 present.
!        Note: All eigenvalues are calculated if none of the arguments
! 	 VL, VU, IL and IU are present.
! M      Optional (output) INTEGER.
!        The total number of eigenvalues found. 0 <= M <= size(A,1).
!        Note: If IL and IU are present then M = IU - IL + 1.
! IFAIL  Optional (output) INTEGER array, shape (:) with size(IFAIL)=n.
!        If INFO = 0, the first M elements of IFAIL are zero.
!        If INFO > 0, then IFAIL contains the indices of the 
! 	 eigenvectors that failed to converge.
!        Note: If Z is present then IFAIL should also be present.
! Q      Optional, (Output) REAL or COMPLEX square array, shape(:,:)
!        with size(Q,1) = n.
!        If Z is present, the matrix used in the reduction of 
! 	 A*z = lambda*B*z to tridiagonal form.
! ABSTOL Optional (input) REAL.
!        The absolute error tolerance for the eigenvalues. An 
! 	 approximate eigenvalue is accepted as converged when it is
! 	 determined to lie in an interval [a,b] of width less than or
! 	 equal to
!              ABSTOL + EPSILON(1.0_<wp>) * max(|a|,|b|),
!        where <wp> is the working precision. If ABSTOL <= 0, then 
! 	 EPSILON(1.0_<wp>) * ||T||1 will be used in its place, where
! 	 ||T||1 is the l1 norm of the tridiagonal matrix obtained by
! 	 reducing the generalized eigenvalue problem to tridiagonal 
! 	 form. Eigenvalues will be computed most accurately when ABSTOL
! 	 is set to twice the underflow threshold 
! 	            2 * LA_LAMCH(1.0_<wp>, 'S'), not zero.
!        Default value: 0.0_<wp>.
!        Note: If this routine returns with 0 < INFO <= n, then some 
! 	 eigenvectors did not converge.
!        Try setting ABSTOL to 2 * LA_LAMCH(1.0_<wp>, 'S').
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value
!        > 0: the algorithm failed to converge or matrix B is not 
! 	    positive definite:
!           <= n: the algorithm failed to converge; if INFO = i, 
! 	      then i eigenvectors failed to converge. Their indices 
! 	      are stored in array IFAIL.
!           > n: if INFO = n + i, for 1 <= i <= n, then the leading 
! 	      minor of order i of B is not positive definite. The 
! 	      factorization of B could not be completed and no
! 	      eigenvalues or eigenvectors were computed.
!        If INFO is not present and an error occurs, then the program 
! 	 is terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SBGVX'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBZ, LUPLO, LRANGE
      INTEGER :: N, LINFO, LDAB, LDBB, LIL, LIU, LM, ISTAT, &
     &  SIFAIL, S1Z, S2Z, S1Q, S2Q, KAB, KBB
      INTEGER, TARGET :: ISTAT1(1)
      REAL(WP), TARGET :: LLZ(1,1), LLQ(1,1)
      REAL(WP) :: LABSTOL, LVL, LVU
!  .. LOCAL ARRAYS ..
      INTEGER, POINTER :: IWORK(:), LIFAIL(:)
      REAL(WP), POINTER :: WORK(:), LQ(:,:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC HUGE, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0
   KAB = SIZE(AB,1)-1; KBB = SIZE(BB,1)-1
   N = SIZE(AB,2); LDAB = MAX(1,SIZE(AB,1))
   LDBB = MAX(1, SIZE(BB,1))
   IF( PRESENT(IFAIL) )THEN; SIFAIL = SIZE(IFAIL); ELSE; SIFAIL = N; END IF
   IF( PRESENT(Q) )THEN; S1Q = SIZE(Q,1); S2Q = SIZE(Q,2)
   ELSE; S1Q = N; S2Q = N; END IF
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
      ELSE; LRANGE = 'A' ; END IF
        
      IF( KAB < 0 .OR. N < 0 ) THEN; LINFO = -1
      ELSE  IF (KBB < 0) THEN; LINFO = -2
      ELSE IF( SIZE( W ) /= N )THEN; LINFO = -3
      ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -4
      ELSE IF( PRESENT(Z) .AND. ( S1Z /= N .OR. S2Z /= N ) )THEN; LINFO = -5
      ELSE IF( LVU < LVL )THEN; LINFO = -6
      ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. &
     &    (PRESENT(IL) .OR. PRESENT(IU)) )THEN; LINFO = -7
      ELSE IF( LRANGE == 'I' .AND. ( LIU.LT.MIN( N, LIL ) .OR. LIU.GT.N)) THEN; LINFO = -8
      ELSE IF( N < LIU )THEN; LINFO = -9
      ELSE IF( SIFAIL /= N .OR. PRESENT(IFAIL).AND..NOT.PRESENT(Z) )THEN; LINFO = -10
      ELSE IF( S1Q /= N .OR. S2Q /= N .OR. PRESENT(Q).AND..NOT.PRESENT(Z) )THEN; LINFO = -11
      ELSE IF( N > 0 )THEN
          IF( PRESENT(Z) ) THEN; LJOBZ = 'V'
            IF( PRESENT(IFAIL) )THEN; LIFAIL => IFAIL
            ELSE; ALLOCATE( LIFAIL(N), STAT=ISTAT ); END IF
            IF( ISTAT == 0 )THEN
              IF( PRESENT(Q) )THEN; LQ => Q
              ELSE; ALLOCATE( LQ(N,N), STAT=ISTAT ); END IF
            END IF
          ELSE; LJOBZ = 'N'; LIFAIL => ISTAT1; LQ => LLQ; S1Q =1; ENDIF
! .. DETERMINE THE WORKSPACE
          IF( ISTAT == 0 ) THEN
            ALLOCATE(IWORK(MAX(1,5*N)), WORK(MAX(1,7*N)), STAT=ISTAT)
            IF( ISTAT == 0 )THEN
              IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL
              ELSE; LABSTOL = 2*LAMCH_F77('Safe minimum'); ENDIF
	       IF (PRESENT (Z))  THEN
                CALL SBGVX_F77( LJOBZ, LRANGE, LUPLO, N, KAB, KBB, AB, LDAB, &
     &            BB, LDBB, LQ, S1Q, LVL, LVU, LIL, LIU, LABSTOL, LM, W, Z,&
     &            S1Z, WORK, IWORK, LIFAIL, LINFO )
               ELSE
	        CALL SBGVX_F77( LJOBZ, LRANGE, LUPLO, N, KAB, KBB, AB, LDAB, &
     &            BB, LDBB, LQ, S1Q, LVL, LVU, LIL, LIU, LABSTOL, LM, W, LLZ,&
     &            S1Z, WORK, IWORK, LIFAIL, LINFO )
               ENDIF
                  IF( PRESENT(M) ) M = LM
                ELSE; LINFO = -100; END IF
                END IF
        IF( PRESENT(Z) .AND. .NOT.PRESENT(IFAIL) ) DEALLOCATE( LIFAIL, STAT=ISTAT1(1) )
        IF( PRESENT(Z) .AND. .NOT.PRESENT(Q) ) DEALLOCATE( LQ, STAT=ISTAT1(1) )
        DEALLOCATE(IWORK, WORK, STAT=ISTAT1(1))
      END IF
      CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSBGVX_F95

              
