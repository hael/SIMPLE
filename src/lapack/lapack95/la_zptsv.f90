 SUBROUTINE ZPTSV_F95( D, E, B, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO
    USE F77_LAPACK, ONLY: PTSV_F77 => LA_PTSV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    REAL(WP), INTENT(INOUT) :: D(:)
    COMPLEX(WP), INTENT(INOUT) :: E(:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_PTSV computes the solution to a linear system of equations 
! A*X = B, where A has tridiagonal form and is real symmetric or complex
! Hermitian and, in either case, positive definite, and where X and B are
! rectangular matrices or vectors. A is factored as A = L*D*L^H, where L 
! is a unit lower bidiagonal matrix and D is a diagonal matrix. The 
! factored form of A is then used to solve the above system.
! 
! =========
! 
!        SUBROUTINE LA_PTSV( D, E, B, INFO=info )
!            REAL(<wp>), INTENT(INOUT) :: D(:)
!            <type>(<wp>), INTENT(INOUT) :: E(:), <rhs>
!            INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!            <type> ::= REAL | COMPLEX
!            <wp>   ::= KIND(1.0) | KIND(1.0D0)
!            <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! D      (input/output) REAL array, shape (:) with size(D) = n, where n 
!        is the order of A.
!        On entry, the diagonal of A.
!        On exit, the diagonal of D.
! E      (input/output) REAL or COMPLEX array, shape (:), with 
!        size(E) = n-1.
!        On entry, the subdiagonal of A.
!        On exit, the subdiagonal of L.
! B      (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(B,1) = n or shape (:) with size(B) = n.
!        On entry, the matrix B.
!        On exit, the solution matrix X.
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: if INFO = i, the leading minor of order i of A is not 
!             positive definite, and the solution has not been computed.
! 	    The factorization has not been completed unless i = n.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_PTSV'
!   .. LOCAL SCALARS ..
    INTEGER :: LINFO, N, NRHS
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0
    N = SIZE(D); NRHS = SIZE(B,2)
!   .. TEST THE ARGUMENTS
    IF( N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( E ) /= N-1 .AND. N /= 0 ) THEN; LINFO = -2
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -3
    ELSE IF ( N > 0 ) THEN
       CALL PTSV_F77( N, NRHS, D, E, B, N, LINFO )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO )
 END SUBROUTINE ZPTSV_F95
