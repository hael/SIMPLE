SUBROUTINE DGGGLM_F95( A, B, D, X, Y, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO
    USE F77_LAPACK, ONLY: GGGLM_F77 => LA_GGGLM
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:), D(:)
    REAL(WP), INTENT(OUT) :: X(:), Y(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!        LA_GGGLM solves the general (Gauss-Markov) linear model (GLM) 
! problem:
!      min ||y||2  subject to d=A*x + B*y
!       x
! where A and B are real or complex rectangular matrices and d is a real
! or complex vector. Further, A is n by m, B is n by p, and d is n by 1,
! and it is assumed that m <= n <= m+p, rank(A) = m, rank(A, B) = n.
! These conditions ensure that the GLM problem has unique solution 
! vectors x and y. The problem is solved using the generalized QR 
! factorization of A and B.
!        If matrix B is square and nonsingular, then the GLM problem is 
! equivalent to the weighted linear least squares problem
!         min ||B^-1 * (d-A*x)||2 
! 	   x
! 
! =========
! 
!       SUBROUTINE LA_GGGLM( A, B, D, X, Y, INFO=info )
!           <type>(<wp>), INTENT( INOUT ) :: A( :, : ), B(:,:), D(:)
!           <type>(<wp>), INTENT( OUT ) :: X(:), Y(:)
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!           <type> ::= REAL | COMPLEX
!           <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! 
! Arguments
! =========
! 
! A     (input/output) REAL or COMPLEX array, shape (:,:) with 
!       size(A,1) = n and size(A,2) = m.
!       On entry, the matrix A.
!       On exit, the contents of A are destroyed.
! B     (input/output) REAL or COMPLEX array, shape (:,:) with 
!       size(B,1) = n and size(B,2) = p.
!       On entry, the matrix B.
!       On exit, the contents of B are destroyed.
! D     (input/output) REAL or COMPLEX array, shape (:) with 
!       size(D) = n.
!       On entry, the vector d.
!       On exit, the contents of D are destroyed.
! X     (output) REAL or COMPLEX array, shape (:) with size(X) = m.
!       The solution vector x.
! Y     (output) REAL or COMPLEX array, shape (:) with size(Y) = p.
!       The solution vector y.
! INFO  Optional (output) INTEGER.
!       = 0: successful exit
!       < 0: if INFO = -i, the i-th argument had an illegal value.
!       If INFO is not present and an error occurs, then the program is
!       terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GGGLM'
!   .. LOCAL SCALARS ..
    INTEGER :: LINFO, ISTAT, ISTAT1, N, M, MN, P
!   .. LOCAL SAVE SCALARS ..
    INTEGER, SAVE :: LWORK = 0
!   .. LOCAL POINTERS ..
    REAL(WP), POINTER :: WORK(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, MAX, MIN
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; N = SIZE(A,1); M = SIZE(A,2); P = SIZE(B,2)
    MN = MIN(M,N)
!   .. TEST THE ARGUMENTS
    IF( M < 0 .OR. N < M ) THEN; LINFO = -1
    ELSE IF( P < 0 .OR. P < N-M .OR. SIZE(B,1) /= N ) THEN; LINFO = -2
    ELSE IF( SIZE(D) /= N ) THEN; LINFO = -3
    ELSE IF( SIZE(X) /= M ) THEN; LINFO = -4
    ELSE IF( SIZE(Y) /= P ) THEN; LINFO = -5
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
          CALL GGGLM_F77( N, M, P, A, MAX(1,N), B, MAX(1,N), &
                          D, X, Y, WORK, LWORK, LINFO )
          IF( LINFO == 0 ) LWORK = INT( WORK(1) )
       ELSE; LINFO = -100; END IF
       DEALLOCATE(WORK, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE DGGGLM_F95
