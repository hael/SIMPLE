SUBROUTINE DGESVD_F95( A, S, U, VT, WW, JOB, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: GESVD_F77 => LA_GESVD
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), OPTIONAL, INTENT(IN) :: JOB
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(OUT) :: S(:)
   REAL(WP), INTENT(OUT), OPTIONAL :: WW(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: U(:,:), VT(:,:)
!----------------------------------------------------------------------
! 
!  Purpose
!  =======
!  
!       LA_GESVD and LA_GESDD compute the singular values and, 
! optionally, the left and/or right singular vectors from the singular
! value decomposition  (SVD) of a real or complex m by n matrix A. The
! SVD of A is written
!                    A = U * SIGMA * V^H
! where SIGMA is an  m by n matrix which is zero except for its 
! min(m, n) diagonal elements, U is an m by m orthogonal (unitary) 
! matrix, and V is an n by n orthogonal (unitary) matrix. The diagonal
! elements of SIGMA , i.e., the values 
! 
!      sigma(i)= SIGMA(i,i), i = 1, 2,..., min(m, n)
! are the singular values of A; they are real and non-negative, and are
! returned in descending order. The first min(m, n) columns of U and V 
! are the left and right singular vectors of A, respectively.
! LA_GESDD solves the same problem as LA_GESVD but uses a divide and 
! conquer method if singular vectors are desired. For large matrices it
! is usually much faster than LA_GESVD when singular vectors are 
! desired, but uses more workspace.
! 
! Note: The routine returns V^H , not V .
! 
! ========
! 
!    SUBROUTINE LA_GESVD / LA_GESDD( A, S, U=u, VT=vt, &
!              WW=ww, JOB=job, INFO=info )  
!      <type>(<wp>), INTENT(INOUT) :: A(:,:)
!      REAL(<wp>), INTENT(OUT) :: S(:)
!      <type>(<wp>), INTENT(OUT), OPTIONAL :: U(:,:), VT(:,:)
!      REAL(<wp>), INTENT(OUT), OPTIONAL :: WW(:)
!      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOB
!      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!      where
!      <type> ::= REAL | COMPLEX
!      <wp>   ::= KIND(1.0) | KIND(1.0D0)
!  
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX array, shape (:, :) with 
!        size(A, 1) = m and size(A, 2) = n.
!        On entry, the matrix A.
!        On exit, if JOB = 'U' and U is not present, then A is 
!        overwritten with the first min(m, n) columns of U (the left
!        singular vectors, stored columnwise).
!        If JOB = 'V' and VT is not present, then A is overwritten with
!        the first min(m, n) rows of V^H (the right singular vectors, 
!        stored rowwise).
!        In all cases the original contents of A are destroyed.
! S      (output) REAL array, shape (:) with size(S) = min(m, n).
!        The singular values of A, sorted so that S(i) >= S(i+1).
! U      Optional (output) REAL or COMPLEX array, shape (:, :) with 
!        size(U, 1) = m  and size(U, 2) = m or min(m, n).
!        If size(U, 2) = m, U contains the m by m matrix U .
!        If size(U; 2) = min(m, n), U contains the first min(m, n) 
!        columns of U (the left singular vectors, stored columnwise).
! VT     Optional (output) REAL or COMPLEX array, shape (:, :) with 
!        size(VT, 1) = n or min(m, n) and size(VT, 2) = n.
!        If size(VT, 1) = n , VT contains the n by n matrix V^H .
!        If size(VT, 1) = min(m, n), VT contains the first min(m, n)
!        rows of V^H (the right singular vectors, stored rowwise).
! WW     Optional (output) REAL array, shape (:) with size(WW) = 
!        min(m, n) - 1
!        If INFO > 0, WW contains the unconverged superdiagonal elements
!        of an upper bidiagonal matrix B whose diagonal is in SIGMA (not
!        necessarily sorted). B has the same singular values as A.
!        Note: WW is a dummy argument for LA_GESDD.
! JOB    Optional (input) CHARACTER(LEN=1).
!        = 'N': neither columns of U nor rows of V^H are returned in 
!          array A.
!        = 'U': if U is not present, the first min(m, n) columns of U 
!          (the left singular vectors) are returned in array A;
!        = 'V': if VT is not present, the first min(m, n) rows of V^H 
!          (the right singular vectors) are returned in array A;
!        Default value: 'N'.
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: The algorithm did not converge.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GESVD'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOB
   CHARACTER(LEN=1) :: LJOBU, LJOBVT
   INTEGER, SAVE :: LWORK = 0
   INTEGER :: N, M, LINFO, LD, ISTAT, ISTAT1, S1U, S2U, S1VT, S2VT, &
              NN, MN, SWW
!  .. LOCAL ARRAYS ..
   REAL(WP), TARGET :: LLU(1,1), LLVT(1,1)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MIN, MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2)
   LD = MAX(1,M); MN = MIN(M,N)
   IF( PRESENT(JOB) )THEN; LJOB = JOB; ELSE; LJOB = 'N'; ENDIF
   IF( PRESENT(U) )THEN; S1U = SIZE(U,1); S2U = SIZE(U,2)
   ELSE; S1U = 1; S2U = 1; END IF
   IF( PRESENT(VT) )THEN; S1VT = SIZE(VT,1); S2VT = SIZE(VT,2)
   ELSE; S1VT = 1; S2VT = 1; END IF
   IF( PRESENT(WW) )THEN; SWW = SIZE(WW); ELSE; SWW = MN-1; ENDIF
!  .. TEST THE ARGUMENTS
   IF( M < 0 .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE( S ) /= MN )THEN; LINFO = -2
   ELSE IF( PRESENT(U) .AND. ( S1U /= M .OR. &
            ( S2U /= M .AND. S2U /= MN ) ) )THEN; LINFO = -3
   ELSE IF( PRESENT(VT) .AND. ( ( S1VT /= N .AND. S1VT /= MN ) &
            .OR. S2VT /= N ) )THEN; LINFO = -4
   ELSE IF( SWW /= MN-1 .AND. MN > 0 ) THEN; LINFO = -5
   ELSE IF( PRESENT(JOB) .AND. ( .NOT. ( LSAME(LJOB,'U') .OR. &
            LSAME(LJOB,'V') .OR. LSAME(LJOB,'N') ) .OR. &
            LSAME(LJOB,'U') .AND. PRESENT(U) .OR. &
            LSAME(LJOB,'V') .AND. PRESENT(VT)) )THEN; LINFO = -6
   ELSE
      IF( PRESENT(U) )THEN
         IF( S2U == M )THEN; LJOBU = 'A'; ELSE; LJOBU = 'S'; ENDIF
      ELSE; IF( LSAME(LJOB,'U') ) THEN; LJOBU = 'O'
         ELSE; LJOBU = 'N'; ENDIF
      ENDIF
      IF( PRESENT(VT) )THEN
         IF( S1VT == N )THEN; LJOBVT = 'A'; ELSE; LJOBVT = 'S'; ENDIF
      ELSE; IF( LSAME(LJOB,'V') )THEN; LJOBVT = 'O'
         ELSE; LJOBVT = 'N'; ENDIF
      ENDIF
      IF( ISTAT == 0 )THEN
         NN = MAX( 5, 3*MN + MAX(M,N), 5*MN )
	 LWORK = MAX( 1, NN, LWORK ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
	 IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
	  LWORK = MAX( 1, NN); ALLOCATE(WORK(LWORK), STAT=ISTAT)
	  IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
         END IF
       END IF
       IF( ISTAT == 0 ) THEN 
          IF( PRESENT(U) ) THEN
            IF ( PRESENT(VT) )THEN
              CALL GESVD_F77( LJOBU, LJOBVT, M, N, A, LD, S, U, MAX(1,S1U), &
             VT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ELSE
             CALL GESVD_F77( LJOBU, LJOBVT, M, N, A, LD, S, U, MAX(1,S1U), &
 &             LLVT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ENDIF
         ELSE
           IF ( PRESENT(VT) )THEN
             CALL GESVD_F77( LJOBU, LJOBVT, M, N, A, LD, S, LLU, MAX(1,S1U), &
 &             VT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ELSE
             CALL GESVD_F77( LJOBU, LJOBVT, M, N, A, LD, S, LLU, MAX(1,S1U), &
 &             LLVT, MAX(1,S1VT), WORK, LWORK, LINFO )
           ENDIF
         ENDIF        
         LWORK = INT(WORK(1)+1)
         IF( LINFO > 0 .AND. PRESENT(WW) ) WW(1:MN-1) = WORK(2:MN)
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DGESVD_F95
