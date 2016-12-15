 SUBROUTINE DGTSV_F95( DL, D, DU, B, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO
    USE F77_LAPACK, ONLY: GTSV_F77 => LA_GTSV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    REAL(WP), INTENT(INOUT) :: DL(:), D(:), DU(:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_GTSV computes the solution to a real or complex linear system of
! equations A*X = B, where A is a square tridiagonal matrix and X and B 
! are rectangular matrices or vectors. The LU decomposition is used to
! factor the matrix A as A = L*U , where L is a product of permutation
! and unit lower bidiagonal matrices and U is upper triangular with 
! nonzeros in only the main diagonal and first two superdiagonals.
! The factored form of A is then used to solve the above system.
! 
! Note: The system A^T*X = B may be solved by interchanging the order of
! the arguments DU and DL.
! 
! =========
! 
!       SUBROUTINE LA_GTSV( DL, D, DU, B, INFO=info ) 
!           <type>(<wp>), INTENT(INOUT) :: DL(:), D(:), DU(:), <rhs>
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!           <type> ::= REAL | COMPLEX
!           <wp> ::= KIND(1.0) | KIND(1.0D0)
!           <rhs> ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! DL    (input/output) REAL or COMPLEX array, shape (:) with 
!       size(DL) = n-1, where n is the order of A.
!       On entry, the subdiagonal of A.
!       On exit, the n-2 elements of the second superdiagonal of U in
!       DL(1),..., DL(n-2).
! D     (input/output) REAL or COMPLEX array, shape (:) with size(D) = n.
!       On entry, the diagonal of A.
!       On exit, the diagonal of U .
! DU    (input/output) REAL or COMPLEX array, shape (:) with 
!       size(DL) = n-1.
!       On entry, the superdiagonal of A.
!       On exit, the first superdiagonal of U .
! B     (input/output) REAL or COMPLEX array, shape (:,:) with 
!       size(B,1) = n or shape (:) with size(B) = n.
!       On entry, the matrix B.
!       On exit, the solution matrix X .
! INFO  Optional (output) INTEGER
!       = 0: successful exit.
!       < 0: if INFO = -i, the i-th argument had an illegal value.
!       > 0: if INFO = i, then U(i,i) = 0. The factorization has not been
!       completed unless i = n. The factor U is singular, so the solution
!       could not be computed.
!       If INFO is not present and an error occurs, then the program is 
!       terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GTSV'
!   .. LOCAL SCALARS ..
    INTEGER :: LINFO, N, NRHS
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0
    N = SIZE(D); NRHS = SIZE(B,2)
!   .. TEST THE ARGUMENTS
    IF( SIZE( DL ) /= N-1 .AND. N /= 0 ) THEN; LINFO = -1
    ELSE IF( N < 0 ) THEN; LINFO = -2
    ELSE IF( SIZE( DU ) /= N-1 .AND. N/=0 ) THEN; LINFO = -3
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -4
    ELSE IF ( N > 0 ) THEN
       CALL GTSV_F77( N, NRHS, DL, D, DU, B, N, LINFO )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO )
 END SUBROUTINE DGTSV_F95
