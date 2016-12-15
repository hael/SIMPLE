 SUBROUTINE DGELS_F95( A, B, TRANS, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME, LA_WS_GELS
    USE F77_LAPACK, ONLY: GELS_F77 => LA_GELS
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!       LA_GELS computes the minimum-norm least squares solution to one 
! or more real or complex linear systems of the form A*x = b, A^T*x = b
! or A^H*x = b using a QR or LQ factorization of A. Matrix A is 
! rectangular assumed to be of full rank. The vectors b and correspon-
! ding solution vectors x are the columns of matrices denoted B and X,
! respectively.
! 
! ==========
! 
!       SUBROUTINE LA_GELS( A, B, TRANS=trans, INFO=info )
!          <type>(<wp>), INTENT( INOUT ) :: A( :, : ), <rhs>
! 	   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: TRANS
! 	   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!          <type> ::= REAL | COMPLEX
! 	   <wp>   ::= KIND(1.0) | KIND(1.0D0)
!          <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A        (input/output) REAL or COMPLEX rectangular array, shape (:,:).
!          On entry, the matrix A.
! 	   On exit, if size(A,1) >= size(A,2), A is overwritten by
! 	   details of its QR factorization. If size(A,1) < size(A,2), A 
! 	   is overwritten by details of its LQ factorization.
! B        (input/output) REAL or COMPLEX array, shape (:,:) with 
!          size(B,1) = max(size(A,1), size(A,2)) or shape (:) with 
! 	   size(B) = max(size(A,1); size(A,2)).
! 	   On entry, the matrix B.
! 	   On exit, the solution matrix X. There are four cases:
! 	   1. If TRANS = 'N' and size(A,1) >= size(A,2), then rows 1 
! 	      to size(A,2) of B contain, columnwise, the least squares
! 	      solution vector(s); the residual sum of squares for the 
! 	      solution vector in a column of B is given by the sum of
! 	      squares of elements in rows size(A,2)+1 to size(A,1) of
! 	      that column.
! 	   2. If TRANS = 'N' and size(A,1) < size(A,2), then rows 1 
! 	      to size(A,2) of B contain, columnwise, the minimum norm 
! 	      solution vector(s).
! 	   3. If TRANS = 'T' or TRANS = 'C', and size(A,1)>=size(A,2),
! 	      then rows 1 to size(A,1) of B contain, columnwise, the 
! 	      minimum norm solution vector(s).
! 	   4. If TRANS = 'T' or TRANS = 'C', and size(A,1) < size(A,2),
! 	      then rows 1 to size(A,1) of B contain, columnwise, the 
! 	      least squares solution vector(s); the residual sum of 
! 	      squares for the solution vector in a column of B is given
! 	      by the sum of squares of elements in rows size(A,1)+1 to 
! 	      size(A,2) of that column.
! TRANS    Optional (input) CHARACTER(LEN=1).
!          Specifies the form of the system of equations:
! 	   = 'N': Ax = b (No transpose)
!          = 'T': A^T*x = b (Transpose)
! 	   = 'C': A^H*x = b (Conjugate transpose)
! 	   Default value: 'N'.
! INFO     Optional (output) INTEGER
!          = 0: successful exit.
! 	   < 0: if INFO = -i, the i-th argument had an illegal value.
! 	   If INFO is not present and an error occurs, then the program
! 	   is terminated with an error message.
!---------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GELS'
    CHARACTER(LEN=1), PARAMETER :: VER = 'D'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LTRANS
    INTEGER :: LINFO, ISTAT, ISTAT1, LWORK, N, M, NRHS
!   .. LOCAL POINTERS ..
    REAL(WP), POINTER :: WORK(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT, MAX, MIN
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2); NRHS = SIZE(B,2)
    IF( PRESENT(TRANS) )THEN; LTRANS = TRANS; ELSE; LTRANS = 'N'; ENDIF
!   .. TEST THE ARGUMENTS
    IF( M < 0 .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= MAX(1,M,N) .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( .NOT.( LSAME(LTRANS,'N') .OR. LSAME(LTRANS,'T') ) )THEN; LINFO = -3
    ELSE
!   .. CALCULATE THE OPTIMAL WORKSPACE ..
       LWORK = LA_WS_GELS( VER, M, N, NRHS, LTRANS )
       ALLOCATE( WORK(LWORK), STAT = ISTAT )
       IF( ISTAT /= 0 ) THEN
          DEALLOCATE( WORK, STAT=ISTAT1 ); LWORK = MIN(M,N) + MAX(1,M,N,NRHS)
          ALLOCATE( WORK(LWORK), STAT = ISTAT )
          IF( ISTAT /= 0 ) CALL ERINFO( -200, SRNAME, LINFO )
       END IF
       IF ( ISTAT == 0 ) THEN
!      .. CALL LAPACK77 ROUTINE
          CALL GELS_F77( LTRANS, M, N, NRHS, A, MAX(1,M), B, MAX(1,M,N), &
                         WORK, LWORK, LINFO )
       ELSE; LINFO = -100; END IF
       DEALLOCATE(WORK, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE DGELS_F95
