SUBROUTINE DSPEVX_F95( A, W, UPLO, Z, VL, VU, IL, IU, &
                       M, IFAIL, ABSTOL, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: LAMCH_F77 => DLAMCH
   USE F77_LAPACK, ONLY: SPEVX_F77 => LA_SPEVX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(IN), OPTIONAL :: IL, IU
   INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
   REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAIL(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
   REAL(WP), INTENT(INOUT) :: A(:)
   REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_SPEVX / LA_HPEVX compute selected eigenvalues and, optionally,
! the corresponding eigenvectors of a real symmetric/complex hermitian
! matrix A in packed storage. Eigenvalues and eigenvectors can be
! selected by specifying either a range of values or a range of indices
! for the desired eigenvalues.
! 
! =========
! 
!        SUBROUTINE LA_SPEVX / LA_HPEVX( AP, W, UPLO=uplo, Z=z, &
!                 VL=vl, VU=vu, IL=il, IU=iu, M=m, IFAIL=ifail, &
!                 ABSTOL=abstol, INFO=info )
!          <type>(<wp>), INTENT(INOUT) :: AP(:)
!          REAL(<wp>), INTENT(OUT) :: W(:)
!          CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!          <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!          REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!          INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!          INTEGER, INTENT(OUT), OPTIONAL :: M
!          INTEGER, INTENT(OUT), OPTIONAL :: IFAIL(:)
!          REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!          <type> ::= REAL | COMPLEX
!          <wp> ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
!  AP      (input/output) REAL or COMPLEX array, shape (:) with size(AP)=
!          n*(n+1)/2, where n is the order of A.
!          On entry, the upper or lower triangle of matrix A in packed
!          storage. The elements are stored columnwise as follows:
!          if UPLO = 'U', AP(i+(j-1)*j/2)=A(i,j) for 1<=i<=j<=n;
!          if UPLO = 'L', AP(i+(j-1)*(2*n-j)/2)=A(i,j) for 1<=j<=i<=n.
!          On exit, AP is overwritten by values generated during the
!          reduction of A to a tridiagonal matrix T . If UPLO = 'U', the
!          diagonal and first superdiagonal of T overwrite the correspond-
!          ing diagonals of A. If UPLO = 'L', the diagonal and first
!          subdiagonal of T overwrite the corresponding diagonals of A.
!  W       (output) REAL array, shape (:) with size(W) = n.
!          The eigenvalues in ascending order.
!  UPLO    Optional (input) CHARACTER(LEN=1).
!          = 'U': Upper triangle of A is stored;
!          = 'L': Lower triangle of A is stored.
!          Default value: 'U'.
! Z        Optional (output) REAL or COMPLEX array, shape (:,:) with 
!          size(Z,1) = n and size(Z,2) = M.
!          The first M columns of Z contain the orthonormal eigenvectors of 
! 	   the matrix A corresponding to the selected eigenvalues, with the
! 	   i-th column of Z containing the eigenvector associated with
!          the eigenvalue in W(i) . If an eigenvector fails to converge, 
! 	   then that column of Z contains the latest approximation to the 
!    	   eigenvector, and the index of the eigenvector is returned in 
! 	   IFAIL.
!          Note: The user must ensure that at least M columns are supplied
! 	   in the array Z. When the exact value of M is not known in 
! 	   advance, an upper bound must be used. In all cases M <= n.
! VL,VU    Optional (input) REAL.
!          The lower and upper bounds of the interval to be searched for 
! 	   eigenvalues. VL < VU.
!          Default values: VL = -HUGE(<wp>) and VU = HUGE(<wp>), where 
! 	   <wp> ::= KIND(1.0) | KIND(1.0D0).
!          Note: Neither VL nor VU may be present if IL and/or IU is 
! 	   present.
! IL,IU    Optional (input) INTEGER.
!          The indices of the smallest and largest eigenvalues to be 
!          returned. The IL-th through IU-th
!          eigenvalues will be found. 1<=IL<=IU<=size(A,1).
!          Default values: IL = 1 and IU = size(A,1).
!          Note: Neither IL nor IU may be present if VL and/or VU is 
!          present.
!          Note: All eigenvalues are calculated if none of the arguments
! 	   VL, VU, IL and IU are present.
! M        Optional (output) INTEGER.
!          The total number of eigenvalues found. 0<=M<=size(A,1).
!          Note: If IL and IU are present then M = IU - IL + 1.
! IFAIL    Optional (output) INTEGER array, shape (:) with 
!          size(IFAIL) = n.
!          If INFO = 0, the first M elements of IFAIL are zero.
!          If INFO > 0, then IFAIL contains the indices of the 
! 	   eigenvectors that failed to converge.
!          Note: If Z is present then IFAIL should also be present.
! ABSTOL   Optional (input) REAL.
!          The absolute error tolerance for the eigenvalues. An 
! 	   approximate eigenvalue is accepted as converged when it is
! 	   determined to lie in an interval [a,b] of width less than or
! 	   equal to ABSTOL+EPSILON(1.0_<wp>) * max(|a|,|b|),
!          where <wp> is the working precision. If ABSTOL<=0, then 
! 	   EPSILON(1.0_<wp>)*||T||1 will be used in its place, where
! 	   ||T||1 is the l1 norm of the tridiagonal matrix obtained by 
! 	   reducing A to tridiagonal form. Eigenvalues will be computed 
! 	   most accurately when ABSTOL is set to twice the underflow 
! 	   threshold 2*LA_LAMCH(1.0_<wp>, 'Safe minimum'), not zero.
!          Default value: 0.0_<wp>.
!          Note: If this routine returns with INFO > 0, then some 
! 	   eigenvectors did not converge. Try setting ABSTOL to 
! 	   2*LA_LAMCH(1.0_<wp>, 'Safe minimum').
! INFO     Optional (output) INTEGER.
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, then i eigenvectors failed to converge. Their
! 	   indices are stored in array IFAIL.
!          If INFO is not present and an error occurs, then the program is
! 	   terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SPEVX'
   CHARACTER(LEN=1) :: LJOBZ, LUPLO, LRANGE
!  .. LOCAL SCALARS ..
   INTEGER :: N, LINFO, LD, LIL, LIU, LM, ISTAT, &
              SIFAIL, S1Z, S2Z, NN
   INTEGER, TARGET :: ISTAT1(1)
   REAL(WP), TARGET :: LLZ(1,1)
   REAL(WP) :: LABSTOL, LVL, LVU
   COMPLEX(WP) :: WW
!  .. LOCAL ARRAYS ..
   INTEGER, POINTER :: IWORK(:), LIFAIL(:)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC HUGE, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; NN = SIZE(A)
   WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW);  LD = MAX(1,N)
   IF( PRESENT(IFAIL) )THEN; SIFAIL = SIZE(IFAIL); ELSE; SIFAIL = N; END IF
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(VL) )THEN; LVL = VL; ELSE; LVL = -HUGE(LVL); ENDIF
   IF( PRESENT(VU) )THEN; LVU = VU; ELSE; LVU = HUGE(LVU); ENDIF
   IF( PRESENT(IL) )THEN; LIL = IL; ELSE; LIL = 1; ENDIF
   IF( PRESENT(IU) )THEN; LIU = IU; ELSE; LIU = N; ENDIF
   IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2)
    ELSE; S1Z = 1; S2Z = 1; ENDIF
!  .. TEST THE ARGUMENTS
   IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= LD .OR. S2Z /= N ) )THEN; LINFO = -4
   ELSE IF( LVU < LVL )THEN; LINFO = -5
   ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. &
     &      (PRESENT(IL) .OR. PRESENT(IU)) )THEN; LINFO = -6
   ELSE IF(( LIU < LIL .OR. LIL < 1) .AND. N > 0 )THEN; LINFO = -7
   ELSE IF( N < LIU )THEN; LINFO = -8
   ELSE IF( SIFAIL /= N .OR. PRESENT(IFAIL).AND..NOT.PRESENT(Z) )THEN; LINFO = -10
   ELSE IF( N > 0 )THEN
      IF( PRESENT(VL) .OR. PRESENT(VU) )THEN; LRANGE = 'V'; LM = N
      ELSE IF( PRESENT(IL) .OR. PRESENT(IU) )THEN; LRANGE = 'I'; LM = LIU-LIL+1
      ELSE; LRANGE = 'A'; LM = N; END IF
      IF( PRESENT(Z) ) THEN; LJOBZ = 'V'
         IF( PRESENT(IFAIL) )THEN; LIFAIL => IFAIL
         ELSE; ALLOCATE( LIFAIL(N), STAT=ISTAT ); END IF
      ELSE; LJOBZ = 'N'; LIFAIL => ISTAT1; ENDIF
!     .. DETERMINE THE WORKSPACE
      IF( ISTAT == 0 ) THEN
         ALLOCATE(IWORK(MAX(1,5*N)), WORK(MAX(1,8*N)), STAT=ISTAT)
         IF( ISTAT == 0 )THEN
            IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL
            ELSE; LABSTOL = 2*LAMCH_F77('Safe minimum'); ENDIF
	    IF (PRESENT (Z)) THEN
               CALL SPEVX_F77( LJOBZ, LRANGE, LUPLO, N, A, LVL, LVU, &
     &                  LIL, LIU, LABSTOL, LM, W, Z, S1Z, WORK, &
     &                  IWORK, LIFAIL, LINFO )
            ELSE
	       CALL SPEVX_F77( LJOBZ, LRANGE, LUPLO, N, A, LVL, LVU, &
     &                  LIL, LIU, LABSTOL, LM, W, LLZ, S1Z, WORK, &
     &                  IWORK, LIFAIL, LINFO )
            ENDIF
            IF( PRESENT(M) ) M = LM
            W(LM+1:N) = 0.0_WP
         ELSE; LINFO = -100; END IF
      END IF
      IF( PRESENT(Z) .AND. .NOT.PRESENT(IFAIL) ) DEALLOCATE( LIFAIL, STAT=ISTAT1(1) )
      DEALLOCATE(IWORK, WORK, STAT=ISTAT1(1))
   END IF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSPEVX_F95
