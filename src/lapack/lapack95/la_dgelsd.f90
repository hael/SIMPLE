SUBROUTINE DGELSD_F95( A, B, RANK, S, RCOND, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO
      USE F77_LAPACK, ONLY: GELSD_F77 => LA_GELSD,  ILAENV_F77 => ILAENV
!   .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL :: RANK
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
      REAL(WP), INTENT(IN), OPTIONAL :: RCOND
!   .. ARRAY ARGUMENTS ..
      REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
      REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: S(:)
!----------------------------------------------------------------------
! 
! Purpose
! ======= 
!
!       LA_GELSS and LA_GELSD compute the minimum-norm least squares 
! solution to one or more real or complex linear systems A*x = b using
! the singular value decomposition of A. Matrix A is rectangular and may
! be rank-deficient. The vectors b and corresponding solution vectors x
! are the columns of matrices denoted B and X , respectively.
!       The effective rank of A is determined by treating as zero those 
! singular values which are less than RCOND times the largest singular 
! value. In addition to X , the routines also return the right singular
! vectors and, optionally, the rank and singular values of A.
!       LA_GELSD combines the singular value decomposition with a divide
! and conquer technique. For large matrices it is often much faster than 
! LA_GELSS but uses more workspace.
! 
! ==========
! 
!        SUBROUTINE LA_GELSS / LA_GELSD( A, B, RANK=rank, S=s, &
!                                          RCOND=rcond, INFO=info )
!             <type>(<wp>), INTENT( INOUT ) :: A( :, : ), <rhs>
!             INTEGER, INTENT(OUT), OPTIONAL :: RANK
!             REAL(<wp>), INTENT(OUT), OPTIONAL :: S(:)
!             REAL(<wp>), INTENT(IN), OPTIONAL :: RCOND
!             INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!             <type> ::= REAL | COMPLEX
!             <wp>   ::= KIND(1.0) | KIND(1.0D0)
!             <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX array, shape (:,:).
!         On entry, the matrix A.
!         On exit, the first min(size(A,1), size(A,2)) rows of A are 
!         overwritten with its right singular vectors, stored rowwise.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = max(size(A,1), size(A,2)) or shape (:) with 
!         size(B) = max(size(A,1), size(A,2)).
!         On entry, the matrix B.
!         On exit, the solution matrix X .
!         If size(A,1) >= size(A,2) and RANK = size(A,2), the residual 
!         sum-of-squares for the solution in a column of B is given by 
!         the sum of squares of elements in rows size(A,2)+1:size(A,1)
!         of that column.
! RANK    Optional (output) INTEGER.
!         The effective rank of A, i.e., the number of singular values 
!         of A which are greater than the product RCOND*sigma1 , where
!  	  sigma1 is the greatest singular value.
! S       Optional (output) REAL array, shape (:) with size(S) = 
!         min(size(A,1), size(A,2)).
!         The singular values of A in decreasing order.
!         The condition number of A in the 2-norm is
! 	  K2(A)= sigma1/sigma(min(size(A,1),size(A,2)) .
! RCOND   Optional (input) REAL.
!         RCOND is used to determine the effective rank of A.
!         Singular values sigma(i)<=RCOND*sigma1  are treated as zero.
!         Default value: 10*max(size(A,1), size(A,2))*EPSILON(1.0_<wp>),
!         where <wp> is the working precision.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: the algorithm for computing the SVD failed to converge; 
! 	  if INFO = i,i off-diagonal elements of an intermediate 
! 	  bidiagonal form did not converge to zero.
!         If INFO is not present and an error occurs, then the program 
! 	  is terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GELSD'
!   .. LOCAL SCALARS ..
      INTEGER :: LINFO, ISTAT, LWORK, N, M, MN, NRHS, LRANK, SS, &
     &  LIWORK, SMLSIZ, NLVL
      REAL(WP) :: LRCOND
!   .. LOCAL POINTERS ..
      REAL(WP), POINTER :: WORK(:)
      REAL(WP), POINTER :: LS(:)
      INTEGER, POINTER :: IWORK(:)
      REAL(WP) :: WORKMIN(1)
      INTEGER :: IWORKMIN(1) 
      DOUBLE PRECISION  TWO
      PARAMETER       ( TWO = 2.0D0 ) 
!   .. INTRINSIC FUNCTIONS ..
      INTRINSIC SIZE, PRESENT, MAX, MIN, EPSILON
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2); NRHS = SIZE(B,2)
    MN = MIN(M,N)
    SMLSIZ = ILAENV_F77( 9, 'DGELSD', ' ', 0, 0, 0, 0 )
    NLVL = INT( LOG( DBLE( MAX(1,MN) ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) )
      LIWORK = 3*MAX(M,N)*(3 * NLVL + 11 )
    IF( PRESENT(RCOND) )THEN; LRCOND = RCOND; ELSE
      LRCOND = 100*EPSILON(1.0_WP) ; ENDIF
    IF( PRESENT(S) )THEN; SS = SIZE(S); ELSE; SS =MN; ENDIF
!   .. TEST THE ARGUMENTS
    IF( M < 0 .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= MAX(1,M,N) .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( SS /= MN ) THEN; LINFO = -4
    ELSE IF( LRCOND <= 0.0_WP ) THEN; LINFO = -5
    ELSE
      IF( PRESENT(S) )THEN; LS => S
      ELSE; ALLOCATE( LS(MN), STAT = ISTAT );
        IF (ISTAT /= 0) THEN
          LINFO = -100
          GOTO 100
        ENDIF
      END IF
! .. DETERMINE THE WORKSPACE ..
! .. QUERING THE SIZE OF WORKSPACE .. 
      LWORK = -1
      CALL GELSD_F77( M, N, NRHS, A, MAX(1,M), B, MAX(1,M,N), &
&       LS, LRCOND, LRANK, WORKMIN, LWORK, IWORKMIN,  LINFO )
      LWORK = WORKMIN(1)
      
      ALLOCATE( WORK(LWORK), STAT = ISTAT )
      IF (ISTAT /= 0) THEN
        LINFO = -100
        GOTO 200
      ENDIF  
      ALLOCATE( IWORK(LIWORK), STAT = ISTAT )
      IF (ISTAT /= 0) THEN
        LINFO = -100
        GOTO 250
      ENDIF
      CALL GELSD_F77( M, N, NRHS, A, MAX(1,M), B, SIZE(B,1), &
&       LS, LRCOND, LRANK, WORK, LWORK, IWORK, LINFO )
      IF( PRESENT(RANK) ) RANK = LRANK
      
300     DEALLOCATE(IWORK)
250	DEALLOCATE(WORK)
200     IF (.NOT. PRESENT(S)) DEALLOCATE(LS)
ENDIF
100     CALL ERINFO( LINFO, SRNAME, INFO, ISTAT ) 
END SUBROUTINE DGELSD_F95
