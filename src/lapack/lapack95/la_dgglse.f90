SUBROUTINE DGGLSE_F95( A, B, C, D, X, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO
    USE F77_LAPACK, ONLY: GGLSE_F77 => LA_GGLSE
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:), C(:), D(:)
    REAL(WP), INTENT(OUT) :: X(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!       LA_GGLSE solves the linear equality-constrained least squares 
! (LSE) problem:
!      min || c - A*x||2 subject to B*x = d,
! where A and B are real or complex rectangular matrices and c and d are
! real or complex vectors. Further, A is m by n, B is p by n, c is m by 1
! and d is p by 1, and it is assumed that
!       p <= n <=  m + p, rank(B) = p,  rank [ A ] = n.
!                                            [ B ]
! These conditions ensure that the LSE problem has a unique solution x. 
! This is obtained using the generalized RQ factorization of the matrices
! B and A.
! 
! =========
! 
!        SUBROUTINE LA_GGLSE( A, B, C, D, X, INFO=info )
!          <type>(<wp>), INTENT( INOUT ) :: A(:,:), B(:,:), C(:), D(:)
!          <type>(<wp>), INTENT( OUT ) :: X(:)
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(A,1) = m and size(A,2) = n.
!         On entry, the matrix A.
!         On exit, the contents of A are destroyed.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = p and size(B,2) = n.
!         On entry, the matrix B.
!         On exit, the contents of B are destroyed.
! C       (input/output) REAL or COMPLEX array, shape (:) with 
!         size(C) = m.
!         On entry, the vector c.
!         On exit, the residual sum of squares for the solution is given
!  	  by the sum of squares of elements n-p+1 to m.
! D       (input/output) REAL or COMPLEX array, shape (:) with 
!         size(D) = p.
!         On entry, The vectors d.
!         On exit, the contents of D are destroyed.
! X       (output) REAL or COMPLEX array, shape (:) with size(X) = n.
!         The solution vector x.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         If INFO is not present and an error occurs, then the program
!         is terminated with an error message.
!----------------------------------------------------------------------	
!   .. PARAMETERS ..
    CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GGLSE'
!   .. LOCAL SCALARS ..
    INTEGER :: LINFO, ISTAT, ISTAT1, N, M, MN, P
!   .. LOCAL SAVE SCALARS ..
    INTEGER, SAVE :: LWORK = 0
!   .. LOCAL POINTERS ..
    REAL(WP), POINTER :: WORK(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, MAX, MIN
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2); P = SIZE(B,1)
    MN = MIN(M,N)
!   .. TEST THE ARGUMENTS
    IF( M < 0 .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( P < 0 .OR. P > N .OR. P < N-M .OR. SIZE(B,2) /= N ) THEN; LINFO = -2
    ELSE IF( SIZE(C) /= M ) THEN; LINFO = -3
    ELSE IF( SIZE(D) /= P ) THEN; LINFO = -4
    ELSE IF( SIZE(X) /= N ) THEN; LINFO = -5
    ELSE
       LWORK = MAX( 1, M + N + P , LWORK )
       ALLOCATE( WORK(LWORK), STAT = ISTAT )
       IF( ISTAT /= 0 ) THEN
          DEALLOCATE( WORK, STAT=ISTAT1 )
          LWORK = MAX( 1, M + N + P )
          ALLOCATE( WORK(LWORK), STAT = ISTAT )
          IF( ISTAT == 0 ) CALL ERINFO( -200, SRNAME, LINFO )
       END IF
       IF ( ISTAT == 0 ) THEN
!      .. CALL LAPACK77 ROUTINE
          CALL GGLSE_F77( M, N, P, A, MAX(1,M), B, MAX(1,P), &
                          C, D, X, WORK, LWORK, LINFO )
          IF( LINFO == 0 ) LWORK = INT(WORK(1))
       ELSE; LINFO = -100; END IF
       DEALLOCATE(WORK, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE DGGLSE_F95
