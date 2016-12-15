      SUBROUTINE ZGELSY_F95( A, B, RANK, JPVT, RCOND, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO
      USE F77_LAPACK, ONLY: GELSY_F77 => LA_GELSY
!   .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL :: RANK
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
      REAL(WP), INTENT(IN), OPTIONAL :: RCOND
!   .. ARRAY ARGUMENTS ..
      INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: JPVT(:)
      COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_GELSY computes the minimum-norm least squares solution to one 
! or more real or complex linear systems A*x = b using a complete 
! orthogonal factorization of A. Matrix A is rectangular and may be 
! rankdeficient. The vectors b and corresponding solution vectors x are 
! the columns of matrices denoted B and X, respectively.
!     The routine computes a QR factorization of A with column pivoting:
!           A * P = Q * [ R11 R12 ]
!                       [  0  R22 ]
! where R11 is the largest leading submatrix whose estimated condition 
! number is less than 1/RCOND. The order of R11, RANK, is the effective 
! rank of A. R22 is considered to be negligible, and R12 is annihilated 
! by orthogonal (unitary) transformations from the right, yielding the 
! complete orthogonal (unitary) factorization
!           A * P = Q * [ T11  0  ] * Z
! 	              [  0   0  ]
! The minimum-norm least squares solution is then
!           x = P * Z^H [ T11^-1 * Q1^H * b ]
! 	              [        0          ]
! where Q1 consists of the first RANK columns of Q.
! 
! =========
! 
!          SUBROUTINE LA_GELSY( A, B, RANK=rank, &
!                     JPVT= jpvt, RCOND= rcond, INFO= info )
!              <type>(<wp>), INTENT(INOUT) :: A(:,:), <rhs>
!              INTEGER, INTENT(OUT), OPTIONAL :: RANK
!              INTEGER, INTENT(INOUT), OPTIONAL :: JPVT(:)
!              REAL(<wp>), INTENT(IN), OPTIONAL :: RCOND
!              INTEGER, INTENT(OUT), OPTIONAL :: INFO
!          where
!              <type> ::= REAL | COMPLEX
!              <wp>   ::= KIND(1.0) | KIND(1.0D0)
!              <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX array, shape (:,:).
!         On entry, the matrix A.
!   	  On exit, A has been overwritten by details of its complete
!         orthogonal factorization.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = max(size(A,1),size(A,2)) or shape (:) with 
! 	  size(B) = max(size(A,1), size(A,2)).
! 	  On entry, the matrix B.
! 	  On exit, rows 1 to size(A,2) contain the solution matrix X .
! 	  If size(A,1) >= size(A,2) and RANK = size(A,2), the residual
! 	  sum-of-squares for the solution vector in a column of B is 
!   	  given by the sum of squares of elements in rows size(A,2)+1 :
!         size(A,1) of that column.
! RANK    Optional (output) INTEGER.
!         The effective rank of A, i.e., the order of the submatrix R11.
! 	  This is the same as the order of the submatrix T11 in the 
! 	  complete orthogonal factorization of A.
! JPVT    Optional (input/output) INTEGER array, shape (:) with 
!         size(JPVT) = size(A,2).
! 	  On entry, if JPVT(i) /= 0, the i-th column of A is an initial 
! 	  column, otherwise it is a free column.
! 	  Before the QR factorization of A, all initial columns are  
! 	  permuted to the leading positions; only the remaining free 
! 	  columns are moved as a result of column pivoting during the 
! 	  factorization.
! 	  On exit, if JPVT(i) = k, then the i-th column of the matrix 
! 	  product A*P was the k-th column of A.
! RCOND   Optional (input) REAL.
!         RCOND is used to determine the effective rank of A. This is 
! 	  defined as the order of the largest leading triangular 
! 	  submatrix R11 in the QR factorization of A, with pivoting, 
! 	  whose estimated condition number < 1/RCOND.
!         Default value: 10*max(size(A,1),size(A,2))*BEPSILON(1.0_<wp>),
! 	  where <wp> is the working precision.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit
! 	  < 0: if INFO = -i, the i-th argument had an illegal value
! 	  If INFO is not present and an error occurs, then the program 
! 	  is terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GELSY'
!   .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, ISTAT1, LWORK, N, M, MN, NRHS, LRANK, SJPVT
      REAL(WP) :: LRCOND
!   .. LOCAL POINTERS ..
      INTEGER, POINTER :: LJPVT(:)
      COMPLEX(WP), POINTER :: WORK(:)
      COMPLEX(WP) :: WORKMIN(1)
      REAL(WP), POINTER :: RWORK(:)
!   .. INTRINSIC FUNCTIONS ..
      INTRINSIC SIZE, PRESENT, MAX, MIN, EPSILON
!   .. EXECUTABLE STATEMENTS ..
      LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2); NRHS = SIZE(B,2)
      MN = MIN(M,N)
      IF( PRESENT(RCOND) )THEN; LRCOND = RCOND; ELSE
          LRCOND = 100*EPSILON(1.0_WP) ; ENDIF
      IF( PRESENT(JPVT) )THEN; SJPVT = SIZE(JPVT); ELSE; SJPVT = N; ENDIF
!   .. TEST THE ARGUMENTS
      IF( M < 0 .OR. N < 0 ) THEN; LINFO = -1
      ELSE IF( SIZE( B, 1 ) /= MAX(1,M,N) .OR. NRHS < 0 ) THEN; LINFO = -2
      ELSE IF( SJPVT /= N ) THEN; LINFO = -4
      ELSE IF( LRCOND <= 0.0_WP ) THEN; LINFO = -5
      ELSE
        IF( PRESENT(JPVT) )THEN; LJPVT => JPVT
        ELSE; ALLOCATE( LJPVT(N), STAT = ISTAT ); LJPVT = 0; END IF
	
        ALLOCATE(RWORK(2*N), STAT=ISTAT)
	IF( ISTAT /= 0 ) CALL ERINFO( -200, SRNAME, LINFO ) 
	
! .. DETERMINE THE WORKSPACE ..
! .. QUERING THE SIZE OF WORKSPACE ..
          LWORK = -1
	  CALL GELSY_F77( M, N, NRHS, A, MAX(1,M), B, MAX(1,M,N), &
     &      LJPVT, LRCOND, LRANK, WORKMIN, LWORK, RWORK, LINFO )
          LWORK = WORKMIN(1)
          IF( ISTAT == 0 ) THEN
            ALLOCATE( WORK(LWORK), STAT = ISTAT )
            IF( ISTAT /= 0 ) CALL ERINFO( -200, SRNAME, LINFO )
          END IF

          IF ( ISTAT == 0 ) THEN
            CALL GELSY_F77( M, N, NRHS, A, MAX(1,M), B, MAX(1,M,N), &
     &      LJPVT, LRCOND, LRANK, WORK, LWORK, RWORK, LINFO )
          ELSE; LINFO = -100; END IF
            IF( PRESENT(RANK) ) RANK = LRANK
            IF( PRESENT(JPVT) ) JPVT = LJPVT
            DEALLOCATE(WORK, RWORK, STAT = ISTAT1 )
          END IF
          CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
        END SUBROUTINE ZGELSY_F95
