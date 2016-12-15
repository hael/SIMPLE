SUBROUTINE ZHBGVD_F95( AB, BB, W, UPLO, Z, INFO ) 
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP =>  DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: HBGVD_F77 => LA_HBGVD
!  .. IMPLICIT STATEMENT .. 
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS .. 
   COMPLEX(WP), INTENT(INOUT) :: AB(:,:), BB(:,:)
   REAL(WP), INTENT(OUT) :: W(:)
   COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_SBGV, LA_SBGVD, LA_HBGV and LA_HBGVD compute all eigenvalues and,
! optionally, all eigenvectors of the generalized eigenvalue problem
!                    A*z = lambda*B*z,
! where A and B are real symmetric in the cases of LA_SBGV and LA_SBGVD 
! and complex Hermitian in the cases of LA_HBGV and LA_HBGVD. Matrix B
! is positive definite. Matrices A and B are stored in a band format.
!    LA_SBGVD and LA_HBGVD use a divide and conquer algorithm. If 
! eigenvectors are desired, they can be much faster than LA_SBGV and 
! LA_HBGV for large matrices but use more workspace.
! 
! =========
! 
!         SUBROUTINE LA_SBGV / LA_SBGVD / LA_HBGV / LA_HBGVD( AB, BB, &
!                                        W, UPLO=uplo, Z=z, INFO=info )
!               <type>(<wp>), INTENT(INOUT) :: AB(:,:), BB(:,:)
!               REAL(<wp>), INTENT(OUT) :: W(:)
!               CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!               <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!               INTEGER, INTENT(OUT), OPTIONAL :: INFO
!         where
!               <type> ::= REAL | COMPLEX
!               <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! AB      (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(AB,1) = ka + 1 and size(AB,2) = n, where ka is the number
! 	  of subdiagonals or superdiagonals in the band and n is the 
!    	  order of A and B.
!         On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') 
! 	  triangle of matrix A in band storage. The ka + 1 diagonals of
! 	  A are stored in the rows of AB so that the j-th column of A
!         is stored in the j-th column of AB as follows:
!         if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j,
! 	                                           1<=j<=n
!         if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka),
! 	                                           1<=j<=n.
!         On exit, the contents of AB are destroyed.
! BB      (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(BB,1) = kb + 1 and size(BB,2) = n, where kb is the number
! 	  of subdiagonals or superdiagonals in the band of B.
!         On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L')
! 	  triangle of matrix B in band storage. The kb + 1 diagonals of
! 	  B are stored in the rows of BB so that the j-th column of B
!         is stored in the j-th column of BB as follows:
!         if UPLO = 'U', BB(kb+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j,
! 	                                           1<=j<=n
!         if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb),
! 	                                           1<=j<=n.
!         On exit, the factor S from the split Cholesky factorization
! 	  B = S^H*S.
! W       (output) REAL array, shape (:) with size(W) = n.
!         The eigenvalues in ascending order.
! UPLO    Optional (input) CHARACTER(LEN=1).
!           = 'U': Upper triangles of A and B are stored;
!           = 'L': Lower triangles of A and B are stored.
!         Default value: 'U'.
! Z       Optional (output) REAL or COMPLEX square array, shape (:,:)
!         with size(Z,1) = n.
!         The matrix Z of eigenvectors, normalized so that Z^H*B*Z = I.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: the algorithm failed to converge or matrix B is not
! 	    positive definite:
!            <= n: if INFO = i, i off-diagonal elements of an 
! 	      intermediate tridiagonal form did not converge to 
! 	      zero.
!            > n: if INFO = n+i, for 1<=i<=n, then the leading minor of
! 	      order i of B is not positive definite. The factorization
! 	      of B could not be completed and no eigenvalues or
! 	      eigenvectors were computed.
!         If INFO is not present and an error occurs, then the program is
! 	  terminated with an error message.
!-----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_HBGVD'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBZ, LUPLO
      INTEGER :: LINFO, N, ISTAT, ISTAT1, S1Z, S2Z, KAB, KBB, &
          LDAB, LDBB, LWORK, LIWORK, LRWORK, IWORKMIN(1)
!  .. LOCAL ARRAYS ..
      COMPLEX(WP), TARGET :: LLZ(1,1), WORKMIN(1)
      COMPLEX(WP), POINTER :: WORK(:)
      INTEGER, POINTER :: IWORK(:)
  REAL(WP), TARGET :: RWORKMIN(1)
  REAL(WP), POINTER :: RWORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC SIZE, MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
      LINFO = 0; KAB = SIZE(AB,1)-1; N = SIZE(AB,2); LDAB = MAX(SIZE(AB,1),1)
      ISTAT = 0; KBB = SIZE(BB,1)-1; LDBB = MAX(SIZE(BB,1),1)
      IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2); LJOBZ = 'V'
      ELSE; S1Z = 1; S2Z = 1; LJOBZ = 'N'; END IF
        IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
!  .. TEST THE ARGUMENTS
          IF( KAB < 0 .OR. N < 0 ) THEN; LINFO = -1
          ELSE IF( KBB < 0 .OR. SIZE(BB,2) /= N ) THEN; LINFO = -2
          ELSE IF( SIZE(W) /= N )THEN; LINFO = -3
          ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -4
          ELSE IF( PRESENT(Z) .AND. ( S1Z /= N .OR. S2Z /= N ) )THEN; LINFO = -5
          ELSE IF( N > 0 )THEN
	    LIWORK = -1
	    LRWORK = -1
	    LWORK = -1
	    IF (PRESENT (Z)) THEN
	     CALL HBGVD_F77(LJOBZ, LUPLO, N, KAB, KBB, AB, LDAB, BB, &
     &           LDBB, W, Z, S1Z, WORKMIN, LWORK, RWORKMIN, LRWORK, &
     &           IWORKMIN, LIWORK, LINFO )
             ELSE
	       CALL HBGVD_F77(LJOBZ, LUPLO, N, KAB, KBB, AB, LDAB, BB, &
     &           LDBB, W, LLZ, S1Z, WORKMIN, LWORK, RWORKMIN, LRWORK, &
     &           IWORKMIN, LIWORK, LINFO )
             ENDIF
             LWORK = WORKMIN(1)
	     LRWORK = RWORKMIN(1)
	     LIWORK = IWORKMIN(1) 
! THEN NEXT 3 LINES SHOULD BE REMOVED WHEN THE BUG IS FIXED IN LAPACK77     
            LWORK = 2 * LWORK + 1
	    LRWORK =2 * LRWORK + 1 
	    LIWORK =2 * LIWORK + 1
	    
            ALLOCATE(WORK(LWORK), RWORK(LRWORK), IWORK(LIWORK), STAT=ISTAT)
            IF( ISTAT == 0 )THEN
	      IF (PRESENT (Z)) THEN
	        CALL HBGVD_F77( LJOBZ, LUPLO, N, KAB, KBB, AB,LDAB, BB, &
     &            LDBB, W, Z, S1Z,  WORK, LWORK, RWORK, LRWORK, IWORK, &
     &            LIWORK, LINFO )
              ELSE
	        CALL HBGVD_F77( LJOBZ, LUPLO, N, KAB, KBB, AB,LDAB, BB, &
     &            LDBB, W, LLZ, S1Z,  WORK, LWORK, RWORK, LRWORK, IWORK, &
     &            LIWORK, LINFO )
              ENDIF
            ELSE; LINFO = -100; ENDIF
	    DEALLOCATE(WORK, RWORK, IWORK, STAT=ISTAT1)
	  ENDIF
	  CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZHBGVD_F95
