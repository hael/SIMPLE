SUBROUTINE ZHBEVX_F95( A, W, UPLO, Z, VL, VU, IL, IU, &
                       M, IFAIL, Q,  ABSTOL, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: LAMCH_F77 => DLAMCH
   USE F77_LAPACK, ONLY: HBEVX_F77 => LA_HBEVX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
   INTEGER, INTENT(IN), OPTIONAL :: IL, IU
   INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
   REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAIL(:)
   COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:), Q(:,:)
   COMPLEX(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_SBEVX / LA_HBEVX compute selected eigenvalues and, optionally, 
! the corresponding eigenvectors of a real symmetric/complex Hermitian 
! band matrix A. Eigenvalues and eigenvectors can be selected by 
! specifying either a range of values or a range of indices for the
! desired eigenvalues.
! 
! =========
! 
!        SUBROUTINE LA_SBEVX / LA_HBEVX( AB, W, UPLO=uplo, Z=z, &
!                 VL=vl, VU=vu, IL=il, IU=iu, M=m, IFAIL=ifail, &
!                 Q=q, ABSTOL=abstol, INFO=info )
!            <type>(<wp>), INTENT(INOUT) :: AB(:,:)
!            REAL(<wp>), INTENT(OUT) :: W(:)
!            CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!            <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!            REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!            INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!            INTEGER, INTENT(OUT), OPTIONAL :: M
!            INTEGER, INTENT(OUT), OPTIONAL :: IFAIL(:)
!            <type>(<wp>), INTENT(OUT), OPTIONAL :: Q(:,:)
!            REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!            INTEGER, INTENT(OUT), OPTIONAL :: INFO
!         where
!            <type> ::= REAL j COMPLEX
!            <wp> ::= KIND(1.0) j KIND(1.0D0)
! 
! Arguments
! =========
! 
! AB     (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(AB,1) = kd + 1 and size(AB,2) = n, where kd is the number
!        of subdiagonals or superdiagonals in the band and n is the order
!        of A.
!        On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') 
!        triangle of matrix A in band storage. The kd + 1 diagonals of A 
!        are stored in the rows of AB so that the j-th column of A is 
!        stored in the j-th column of AB as follows:
!        if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j
!                                                   1<=j<=n
!        if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd)
!                                                   1<=j<=n.
!        On exit, AB is overwritten by values generated during the 
!        reduction of A to a tridiagonal matrix T . If UPLO = 'U' the 
!        first superdiagonal and the diagonal of T are returned in rows
!        kd and kd + 1 of AB. If UPLO = 'L', the diagonal and first 
!        subdiagonal of T are returned in the first two rows of AB.
! W      (output) REAL array, shape (:) with size(W) = n.
!        The first M elements contain the selected eigenvalues in 
!        ascending order.
! UPLO   Optional (input) CHARACTER(LEN=1).
!        = 'U' : Upper triangle of A is stored;
!        = 'L' : Lower triangle of A is stored.
!        Default value: 'U'.
! Z      Optional (output) REAL or COMPLEX array, shape (:,:) with 
!        size(Z,1) = n and size(Z,2) = M.
!        The first M columns of Z contain the orthonormal eigenvectors of
!        the matrix A corresponding to the selected eigenvalues, with the
!        i-th column of Z containing the eigenvector associated with the 
!        eigenvalue in W(i). If an eigenvector fails to converge, then 
!        that column of Z contains the latest approximation to the 
!        eigenvector, and the index of the eigenvector is returned in 
!        IFAIL.
!        Note: The user must ensure that at least M columns are supplied 
!        in the array Z. When the exact value of M is not known in 
!        advance, an upper bound must be used. In all cases M<=n.
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
!        The total number of eigenvalues found. 0<=M<=size(A,1).
!        Note: If IL and IU are present then M = IU - IL + 1.
! IFAIL  Optional (output) INTEGER array, shape (:) with size(IFAIL) = n.
!        If INFO = 0, the first M elements of IFAIL are zero.
!        If INFO > 0, then IFAIL contains the indices of the eigenvectors
!        that failed to converge.
!        Note: If Z is present then IFAIL should also be present.
! Q      Optional (output) REAL or COMPLEX square array, shape(:,:) with
!        size(Q,1) = n.
!        The n by n unitary matrix used in the reduction to tridiagonal 
!        form. This is computed only if Z is present.
! ABSTOL Optional (input) REAL.
!        The absolute error tolerance for the eigenvalues. An approximate
!        eigenvalue is accepted as converged when it is determined to lie
!        in an interval [a,b] of width less than or equal to
!        ABSTOL + EPSILON(1.0_<wp>) * max(|a|, |b|),
!        where <wp> is the working precision. If ABSTOL<=0, then 
!        EPSILON(1.0_<wp>)* ||T||1 will be used in its place, where
!        ||T||1 is the l1 norm of the tridiagonal matrix obtained by 
!        reducing A to tridiagonal form. Eigenvalues will be computed most
!        accurately when ABSTOL is set to twice the underflow threshold 
!        2 * LA_LAMCH(1.0_<wp>, 'Safe minimum'), not zero.
!        Default value: 0.0_<wp>.
!        Note: If this routine returns with INFO > 0, then some 
!        eigenvectors did not converge. Try setting ABSTOL to 
!        2 * LA_LAMCH(1.0_<wp>, 'Safe minimum').
! INFO   Optional (output) INTEGER
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: if INFO = i, then i eigenvectors failed to converge. Their
!        indices are stored in array IFAIL.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_HBEVX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBZ, LUPLO, LRANGE
   INTEGER :: N, LINFO, LD, LIL, LIU, LM, ISTAT, &
              SIFAIL, S1Z, S2Z, S1Q, S2Q, KD
   INTEGER, TARGET :: ISTAT1(1)
   COMPLEX(WP), TARGET :: LLZ(1,1), LLQ(1,1)
   REAL(WP) :: LABSTOL, LVL, LVU
!  .. LOCAL ARRAYS ..
   INTEGER, POINTER :: IWORK(:), LIFAIL(:)
   COMPLEX(WP), POINTER :: WORK(:), LQ(:,:)
   REAL(WP), POINTER :: RWORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC HUGE, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0
   KD = SIZE(A,1)-1; N = SIZE(A,2); LD = MAX(1,SIZE(A,1))
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
   IF( KD < 0 .OR. N < 0 ) THEN; LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= N .OR. S2Z /= N ) )THEN; LINFO = -4
   ELSE IF( LVU < LVL )THEN; LINFO = -5
   ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. &
            (PRESENT(IL) .OR. PRESENT(IU)) )THEN; LINFO = -6
   ELSE IF(( LIU < LIL .OR. LIL < 1 ) .AND. N>0 )THEN; LINFO = -7
   ELSE IF( N < LIU )THEN; LINFO = -8
   ELSE IF( SIFAIL /= N .OR. PRESENT(IFAIL).AND..NOT.PRESENT(Z) )THEN; LINFO = -10
   ELSE IF( S1Q /= N .OR. S2Q /= N .OR. PRESENT(Q).AND..NOT.PRESENT(Z) )THEN; LINFO = -11
   ELSE IF( N > 0 )THEN
      IF( PRESENT(VL) .OR. PRESENT(VU) )THEN; LRANGE = 'V'; LM = N
      ELSE IF( PRESENT(IL) .OR. PRESENT(IU) )THEN; LRANGE = 'I'; LM = LIU-LIL+1
      ELSE; LRANGE = 'A'; LM = N; END IF
      IF( PRESENT(Z) ) THEN; LJOBZ = 'V'
         IF( PRESENT(IFAIL) )THEN; LIFAIL => IFAIL
         ELSE; ALLOCATE( LIFAIL(N), STAT=ISTAT ); END IF
         IF( ISTAT == 0 )THEN
            IF( PRESENT(Q) )THEN; LQ => Q
            ELSE; ALLOCATE( LQ(N,N), STAT=ISTAT ); END IF
         END IF
      ELSE; LJOBZ = 'N'; LIFAIL => ISTAT1; LQ => LLQ; S1Q =1; ENDIF
!     .. DETERMINE THE WORKSPACE
      IF( ISTAT == 0 ) THEN
         ALLOCATE(IWORK(MAX(1,5*N)), RWORK(MAX(1,7*N)), WORK(MAX(1,N)), STAT=ISTAT)
         IF( ISTAT == 0 )THEN
            IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL
            ELSE; LABSTOL = 2*LAMCH_F77('Safe minimum'); ENDIF
	    IF (PRESENT (Z)) THEN
                 CALL HBEVX_F77( LJOBZ, LRANGE, LUPLO, N, KD, A, LD, LQ, S1Q, LVL, LVU, &
     &                 LIL, LIU, LABSTOL, LM, W, Z, S1Z, WORK, &
     &                 RWORK, IWORK, LIFAIL, LINFO )
            ELSE
	         CALL HBEVX_F77( LJOBZ, LRANGE, LUPLO, N, KD, A, LD, LQ, S1Q, LVL, LVU, &
     &                 LIL, LIU, LABSTOL, LM, W, LLZ, S1Z, WORK, &
     &                 RWORK, IWORK, LIFAIL, LINFO )
           ENDIF
            IF( PRESENT(M) ) M = LM
            W(LM+1:N) = 0.0_WP
         ELSE; LINFO = -100; END IF
      END IF
      IF( PRESENT(Z) .AND. .NOT.PRESENT(IFAIL) ) DEALLOCATE( LIFAIL, STAT=ISTAT1(1) )
      IF( PRESENT(Z) .AND. .NOT.PRESENT(Q) ) DEALLOCATE( LQ, STAT=ISTAT1(1) )
      DEALLOCATE(IWORK, RWORK, WORK, STAT=ISTAT1(1))
   END IF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZHBEVX_F95
