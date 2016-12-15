SUBROUTINE DSTEVR_F95( D, E, W, Z, VL, VU, IL, IU, M, ISUPPZ, &
     &    ABSTOL, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO
      USE F77_LAPACK, ONLY: STEVR_F77 => LA_STEVR, LAMCH_F77 => DLAMCH
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(IN), OPTIONAL :: IL, IU
      INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
      REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: ISUPPZ(:)
      REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
      REAL(WP), INTENT(INOUT) :: D(:), E(:)
      REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
!       LA_STEVR computes selected eigenvalues and, optionally, the 
! corresponding eigenvectors of a real symmetric tridiagonal matrix A.
! Eigenvalues and eigenvectors can be selected by specifying either a 
! range of values or a range of indices for the desired eigenvalues.
!       LA_STEVR uses a relatively robust representation (RRR) algorithm.
! It is usually the fastest algorithm of all and uses the least
! workspace.
! 
! =========
! 
!         SUBROUTINE LA_STEVR ( D, E, W, Z=z, VL=vl, VU=vu, &
!                         IL=il, IU=iu, M=m, ISUPPZ=isuppz, &
!                         ABSTOL=abstol, INFO=info )
!                 REAL(<wp>), INTENT(INOUT) :: D(:), E(:)
!                 REAL(<wp>), INTENT(OUT) :: W(:)
!                 REAL(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!                 INTEGER, INTENT(OUT), OPTIONAL :: ISUPPZ(:)
!                 REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!                 INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!                 INTEGER, INTENT(OUT), OPTIONAL :: M
!                 REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!                 INTEGER, INTENT(OUT), OPTIONAL :: INFO
!         where
!                 <wp> ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! D      (input/output) REAL array, shape (:) with size(D) = n, where n 
!        is the order of A.
!        On entry, the diagonal elements of the matrix A.
!        On exit, the original contents of D possibly multiplied by a 
!        constant factor to avoid over/underflow in computing the 
!        eigenvalues.
! E      (input/output) REAL array, shape (:) with size(E) = n.
!        On entry, the n-1 subdiagonal elements of A in E(1) to E(n-1) .
!        E(n) need not be set.
!        On exit, the original contents of E possibly multiplied by a 
!        constant factor to avoid over/underflow in computing the 
!        eigenvalues.
! W      (output) REAL array with size(W) = n.
!        The first M elements contain the selected eigenvalues in 
!        ascending order.
! Z      Optional (output) REAL or COMPLEX array, shape (:,:) with 
!        size(Z,1) = n and size(Z,2) = M.
!        The first M columns of Z contain the orthonormal eigenvectors of
!        A corresponding to the selected eigenvalues, with the i-th column 
!        of Z containing the eigenvector associated with the eigenvalue in
!        W(i).
!        Note: The user must ensure that at least M columns are supplied 
!        in the array Z. When the exact value of M is not known in advance,
!        an upper bound must be used. In all cases M <= n.
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
!        1 <= IL <= IU <= n.
!        Default values: IL = 1 and IU = n.
!        Note: Neither IL nor IU may be present if VL and/or VU is 
!        present.
!        Note: All eigenvalues are calculated if none of the arguments 
!        VL, VU, IL and IU are present.
! M      Optional (output) INTEGER.
!        The total number of eigenvalues found. 0 <= M <=  n.
!        Note: If IL and IU are present then M = IU - IL + 1.
! ISUPPZ Optional (output) INTEGER array, shape (:) with 
!        size(ISUPPZ) = 2*max(1,M).
!        The support of the eigenvectors in A, i.e., the indices 
!        indicating the nonzero elements. The i-th eigenvector is nonzero
!        only in elements ISUPPZ(2*i-1) through ISUPPZ(2*i).
! ABSTOL Optional (input) REAL.
!        The absolute error tolerance for the eigenvalues. An approximate
!        eigenvalue is accepted as converged when it is determined to lie
!        in an interval [a, b] of width less than or equal to 
!        ABSTOL + EPSILON(1.0_<wp>) * max(| a |, | b |),
!        where <wp> is the working precision. If ABSTOL <= 0, then 
!        EPSILON(1.0_<wp>)*||A||1 will be used in its place. Eigenvalues
!        will be computed most accurately if ABSTOL is set to 
!        LA_LAMCH( 1.0_<wp>, 'Safe minimum'), not zero.
!        Default value: 0.0_<wp>.
! INFO   Optional (output) INTEGER
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: an internal error occurred.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!-----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_STEVR'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBZ, LRANGE
      INTEGER :: N, LD, LIL, LIU, LM, SISUPPZ, S1Z, S2Z, NN
      INTEGER :: LINFO, ISTAT, LWORK, LIWORK
      REAL(WP), TARGET :: LLZ(1,1)
      REAL(WP) :: LABSTOL, LVL, LVU
!  .. LOCAL ARRAYS ..
      INTEGER:: IWORKMIN(1), DUMMY(1)
      REAL(WP), TARGET :: WORKMIN(1)
      INTEGER, POINTER :: IWORK(:), LISUPPZ(:)
      REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC HUGE, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(D); LD = MAX(1,N)
   NN=2*MAX(1,N)
   IF( PRESENT(ISUPPZ) )THEN; SISUPPZ = SIZE(ISUPPZ); ELSE; SISUPPZ = NN; END IF
   IF( PRESENT(VL) )THEN; LVL = VL; ELSE; LVL = -HUGE(LVL); ENDIF
   IF( PRESENT(VU) )THEN; LVU = VU; ELSE; LVU = HUGE(LVU); ENDIF
   IF( PRESENT(IL) )THEN; LIL = IL; ELSE; LIL = 1; ENDIF
   IF( PRESENT(IU) )THEN; LIU = IU; ELSE; LIU = N; ENDIF
   IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2)
   ELSE; S1Z = 1; S2Z = 1; ENDIF
!  .. TEST THE ARGUMENTS
   IF( N < 0 ) THEN; LINFO = -1
   ELSE IF( SIZE( E ) /= N .AND. N/=0)THEN; LINFO = -2
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -3
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= LD .OR. S2Z /= MAX(1,N) ) )THEN; LINFO = -4
   ELSE IF( SISUPPZ /= NN .OR. PRESENT(ISUPPZ).AND..NOT.PRESENT(Z) )THEN; LINFO = -5
   ELSE IF( LVU < LVL .AND. N>0)THEN; LINFO = -6
   ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. &
            (PRESENT(IL) .OR. PRESENT(IU)) )THEN; LINFO = -7
   ELSE IF (((LIU<LIL) .OR. (LIL<1)) .AND. (N>0))THEN; LINFO = -7
   ELSE IF( N < LIU )THEN; LINFO = -8
   ELSE IF( N > 0 )THEN
     IF( PRESENT(VL) .OR. PRESENT(VU) )THEN; LRANGE = 'V'; LM = N
     ELSE IF( PRESENT(IL) .OR. PRESENT(IU) )THEN; LRANGE = 'I'; LM = LIU-LIL+1
     ELSE; LRANGE = 'A'; LM = N; END IF
       IF( PRESENT(Z) ) THEN; LJOBZ = 'V'
       ELSE; LJOBZ = 'N'
       ENDIF
! .. DETERMINE THE WORKSPACE ..
! .. QUERING THE SIZE OF WORKSPACE ..
       LWORK = -1
       LIWORK = -1
       IF (PRESENT(Z)) THEN
         CALL STEVR_F77( LJOBZ, LRANGE, N, D, E,  LVL, LVU, &
     &       LIL, LIU, LABSTOL, LM, W, Z, S1Z, DUMMY, &
     &       WORKMIN, LWORK, IWORKMIN, LIWORK, LINFO )
       ELSE
         CALL STEVR_F77( LJOBZ, LRANGE, N, D, E,  LVL, LVU, &
     &       LIL, LIU, LABSTOL, LM, W, LLZ, S1Z, DUMMY, &
     &       WORKMIN, LWORK, IWORKMIN, LIWORK, LINFO )
       ENDIF
       LWORK = WORKMIN(1)
       LIWORK = IWORKMIN(1)
       
       ALLOCATE(IWORK(LIWORK), STAT=ISTAT)
       IF (ISTAT /= 0) THEN ; LINFO = -100; GOTO 100; ENDIF

       ALLOCATE(LISUPPZ(NN), STAT=ISTAT)
       IF (ISTAT /= 0) THEN ; LINFO = -100; GOTO 200; ENDIF

       ALLOCATE(WORK(LWORK), STAT=ISTAT)
       IF (ISTAT /= 0) THEN ; LINFO = -100; GOTO 300; ENDIF

       IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL
       ELSE; LABSTOL = 2*LAMCH_F77('Safe minimum'); ENDIF

       IF (PRESENT (Z)) THEN
         CALL STEVR_F77( LJOBZ, LRANGE, N, D, E, LVL, LVU, LIL, LIU, &
     &     LABSTOL, LM, W, Z, S1Z, LISUPPZ, WORK, LWORK, &
     &     IWORK, LIWORK, LINFO )
       ELSE
         CALL STEVR_F77( LJOBZ, LRANGE, N, D, E, LVL, LVU, LIL, LIU, &
     &     LABSTOL, LM, W, LLZ, S1Z, LISUPPZ, WORK, LWORK, &
     &     IWORK, LIWORK, LINFO )
       ENDIF

       IF( PRESENT(M) ) M = LM

       DEALLOCATE(WORK)
300    DEALLOCATE(LISUPPZ)
200    DEALLOCATE(IWORK)
     ENDIF
100  CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSTEVR_F95
