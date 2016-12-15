 SUBROUTINE DPOSV_F95( A, B, UPLO, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME
    USE F77_LAPACK, ONLY: POSV_F77 => LA_POSV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_POSV computes the solution to a linear system of equations 
! A*X=B, where A is real symmetric or complex Hermitian and, in either
! case, positive definite, and where X and B are rectangular matrices 
! or vectors. The Cholesky decomposition is used to factor A as
! A = U^H*U if UPLO = 'U', or A = L*L^H if UPLO = 'L'
! where U is an upper triangular matrix and L is a lower triangular 
! matrix (L = U^H ). The factored form of A is then used to solve the
! above system.
! 
! =========
! 
!         SUBROUTINE LA_POSV( A, B, UPLO=uplo, INFO=info )
!            <type>(<wp>), INTENT(INOUT) :: A(:,:), <rhs>
!            CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!            INTEGER, INTENT(OUT), OPTIONAL :: INFO
!         where
!            <type> ::= REAL | COMPLEX
!            <wp>   ::= KIND(1.0) | KIND(1.0D0)
!            <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A     (input/output) REAL or COMPLEX square array, shape (:,:).
!       On entry, the matrix A.
!       If UPLO = 'U', the upper triangular part of A contains the upper 
!       triangular part of the matrix A, and the strictly lower triangular
!       part of A is not referenced. If UPLO = 'L', the lower triangular 
!       part of A contains the lower triangular part of the matrix A, and
!       the strictly upper triangular part of A is not referenced.
!       On exit, the factor U or L from the Cholesky factorization 
!       A = U^H*U = L*L^H.
! B     (input/output) REAL or COMPLEX array, shape (:,:) with 
!       size(B,1) = size(A,1) or shape (:) with size(B) = size(A,1).
!       On entry, the matrix B.
!       On exit, the solution matrix X.
! UPLO  Optional (input) CHARACTER(LEN=1)
!        = 'U': Upper triangle of A is stored;
!        = 'L': Lower triangle of A is stored.
!       Default value: 'U'.
! INFO  Optional (output) INTEGER
!       = 0: sauccessful exit.
!       < 0: if INFO = -i, the i-th argument had an illegal value.
!       > 0: if INFO = i, the leading minor of order i of A is not
!            positive definite, so the factorization could not be 
! 	   completed and the solution could not be computed.
!       If INFO is not present and an error occurs, then the program is
!       terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_POSV'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LUPLO
    INTEGER :: LINFO, N, NRHS
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; N = SIZE(A,1); NRHS = SIZE(B,2)
    IF( PRESENT(UPLO) )THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
!   .. TEST THE ARGUMENTS
    IF( SIZE( A, 2 ) /= N .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
    ELSE IF ( N > 0 ) THEN
       CALL POSV_F77( LUPLO, N, NRHS, A, N, B, N, LINFO )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO )
 END SUBROUTINE DPOSV_F95
