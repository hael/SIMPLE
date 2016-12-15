!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_STEVX computes selected eigenvalues and, optionally, the 
! corresponding eigenvectors of a real symmetric tridiagonal matrix A.
! Eigenvalues and eigenvectors can be selected by specifying either a 
! range of values or a range of indices for the desired eigenvalues.
! 
! =========
! 
!        SUBROUTINE LA_STEVX( D, E, W, Z=z, VL=vl, VU=vu, &
!                         IL=il, IU=iu, M=m, IFAIL=ifail, &
!                         ABSTOL=abstol, INFO=info )
!             REAL(<wp>), INTENT(INOUT) :: D(:), E(:)
!             REAL(<wp>), INTENT(OUT) :: W(:)
!             REAL(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!             REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!             INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!             INTEGER, INTENT(OUT), OPTIONAL :: M
!             INTEGER, INTENT(OUT), OPTIONAL :: IFAIL(:)
!             REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!             INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!              <wp> ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! D        (input/output) REAL array, shape (:) with size(D) = n, where n
!          is the order of A.
!          On entry, the diagonal elements of the matrix A.
!          On exit, the original contents of D possibly multiplied by a 
!          constant factor to avoid over/underflow in computing the
! 	   eigenvalues.
! E        (input/output) REAL array, shape (:) with size(E) = n.
!          On entry, the n-1 subdiagonal elements of A in E(1) to E(n-1).
! 	   E(n) need not be set.
!          On exit, the original contents of E possibly multiplied by a 
! 	   constant factor to avoid over/underflow in computing the 
! 	   eigenvalues.
! W        (output) REAL array with size(W) = n.
!          The first M elements contain the selected eigenvalues in 
! 	   ascending order.
! Z        Optional (output) REAL or COMPLEX array, shape (:,:) with 
!          size(Z,1) = n and size(Z,2) = M.
!          The first M columns of Z contain the orthonormal eigenvectors
! 	   of A corresponding to the selected eigenvalues, with the i-th
! 	   column of Z containing the eigenvector associated with the
!          eigenvalue in W(i) . If an eigenvector fails to converge, then
!          that column of Z contains the latest approximation to the 
!          eigenvector, and the index of the eigenvector is returned in 
! 	   IFAIL.
!          Note: The user must ensure that at least M columns are 
! 	   supplied in the array Z. When the exact value of M is not 
! 	   known in advance, an upper bound must be used. In all cases
! 	   M <= n.
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
! 	   1 <= IL <= IU <= n.
!          Default values: IL = 1 and IU = n.
!          Note: Neither IL nor IU may be present if VL and/or VU is 
! 	   present.
!          Note: All eigenvalues are calculated if none of the arguments
! 	   VL, VU, IL and IU are present.
! M        Optional (output) INTEGER.
!          The total number of eigenvalues found. 0 <= M <= n.
!          Note: If IL and IU are present then M = IU - IL + 1.
! IFAIL    Optional (output) INTEGER array, shape (:) with 
!          size(IFAIL) = n.
!          If INFO = 0, the first M elements of IFAIL are zero.
!          If INFO > 0, then IFAIL contains the indices of the 
! 	   eigenvectors that failed to converge.
!          Note: If Z is present then IFAIL should also be present.
! ABSTOL   Optional (input) REAL.
!          The absolute error tolerance for the eigenvalues. An 
!          approximate eigenvalue is accepted as converged when it is
!          determined to lie in an interval [a,b] of width less than or
!          equal to ABSTOL + EPSILON(1.0_<wp>) * max(|a|,|b|),
!          where <wp> is the working precision. If ABSTOL <= 0, then
! 	   EPSILON(1.0_<wp>) * ||A||1 will be used in its place. 
! 	   Eigenvalues will be computed most accurately when ABSTOL is
! 	   set to twice the underflow threshold 
! 	   2 * LA_LAMCH(1.0_<wp>, 'Safe minimum'), not zero.
!          Default value: 0.0_<wp>.
!          Note: If this routine returns with INFO > 0, then some 
! 	   eigenvectors did not converge. Try setting ABSTOL to 
! 	   2 * LA_LAMCH(1.0_<wp>, 'Safe minimum').
! INFO     Optional (output) INTEGER
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, then i eigenvectors failed to converge.
! 	      Their indices are stored in array IFAIL.
!          If INFO is not present and an error occurs, then the program
! 	   is terminated with an error message.
!----------------------------------------------------------------------
SUBROUTINE DSTEVX_F95( D, E, W, Z, VL, VU, IL, IU, M, &
                            IFAIL, ABSTOL, INFO )
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO
   USE F77_LAPACK, ONLY: STEVX_F77 => LA_STEVX, LAMCH_F77 => DLAMCH
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(IN), OPTIONAL :: IL, IU
   INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
   REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAIL(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
   REAL(WP), INTENT(INOUT) :: D(:), E(:)
   REAL(WP), INTENT(OUT) :: W(:)
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_STEVX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBZ, LRANGE
   INTEGER :: N, LD, LIL, LIU, LM, SIFAIL, S1Z, S2Z
   INTEGER :: LINFO, ISTAT
   INTEGER, TARGET :: ISTAT1(1)
   REAL(WP), TARGET :: LLZ(1,1)
   REAL(WP) :: LABSTOL, LVL, LVU
!  .. LOCAL ARRAYS ..
   INTEGER, POINTER :: IWORK(:), LIFAIL(:)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC HUGE, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(D); LD = MAX(1,N)
   IF( PRESENT(M)) M = 0 
   IF( PRESENT(IFAIL) )THEN; SIFAIL = SIZE(IFAIL); ELSE; SIFAIL = N; END IF
   IF( PRESENT(VL) )THEN; LVL = VL; ELSE; LVL = -HUGE(LVL); ENDIF
   IF( PRESENT(VU) )THEN; LVU = VU; ELSE; LVU = HUGE(LVU); ENDIF
   IF( PRESENT(IL) )THEN; LIL = IL; ELSE; LIL = 1; ENDIF
   IF( PRESENT(IU) )THEN; LIU = IU; ELSE; LIU = N; ENDIF
   IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2)
   ELSE; S1Z = 1; S2Z = 1; ENDIF
!  .. TEST THE ARGUMENTS
   IF( N < 0 ) THEN; LINFO = -1
   ELSE IF( SIZE( E ) /= N .AND. N > 0 )THEN; LINFO = -2
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -3
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= LD .OR. S2Z /= N ) )THEN; LINFO = -4
   ELSE IF( LVU < LVL )THEN; LINFO = -5
   ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. &
     &         (PRESENT(IL) .OR. PRESENT(IU)) )THEN; LINFO = -6
   ELSE IF(( LIU < LIL .OR. LIL < 1) .AND. N>0 )THEN; LINFO = -7
   ELSE IF( N < LIU )THEN; LINFO = -8
   ELSE IF( SIFAIL /= N .OR. PRESENT(IFAIL).AND..NOT.PRESENT(Z) )THEN; LINFO = -10
   ELSE IF( N > 0 )THEN
      IF( PRESENT(VL) .OR. PRESENT(VU) )THEN; LRANGE = 'V'; LM = N
      ELSE IF( PRESENT(IL) .OR. PRESENT(IU) )THEN; LRANGE = 'I'; LM = LIU-LIL+1
      ELSE; LRANGE = 'A'; LM = N; END IF
      IF( PRESENT(Z) ) THEN; LJOBZ = 'V'
         IF( PRESENT(IFAIL) )THEN; LIFAIL => IFAIL
         ELSE; ALLOCATE( LIFAIL(N), STAT=ISTAT ); END IF
      ELSE; LJOBZ = 'N'
      LIFAIL => ISTAT1; ENDIF
!     .. DETERMINE THE WORKSPACE
      IF( ISTAT == 0 ) THEN
         ALLOCATE(IWORK(5*N), WORK(5*N), STAT=ISTAT)
         IF( ISTAT == 0 )THEN
            IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL
            ELSE; LABSTOL = 2*LAMCH_F77('Safe minimum'); ENDIF
	    IF (PRESENT(Z)) THEN
               CALL STEVX_F77( LJOBZ, LRANGE, N, D, E, LVL, LVU, LIL, LIU, &
                            LABSTOL, LM, W, Z, S1Z, WORK, IWORK, LIFAIL, LINFO )
	    ELSE
	       CALL STEVX_F77( LJOBZ, LRANGE, N, D, E, LVL, LVU, LIL, LIU, &
                            LABSTOL, LM, W, LLZ, S1Z, WORK, IWORK, LIFAIL, LINFO )
	    ENDIF		    
            IF( PRESENT(M) ) M = LM
            W(LM+1:N) = 0.0_WP
         ELSE; LINFO = -100; END IF
      END IF
      IF( PRESENT(Z) .AND. .NOT.PRESENT(IFAIL) ) DEALLOCATE( LIFAIL, STAT=ISTAT1(1) )
      DEALLOCATE(IWORK, WORK, STAT=ISTAT1(1))
   END IF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSTEVX_F95
