 SUBROUTINE DPPSV_F95( A, B, UPLO, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME
    USE F77_LAPACK, ONLY: PPSV_F77 => LA_PPSV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    REAL(WP), INTENT(INOUT) :: A(:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_PPSV computes the solution to a linear system of equations 
! A*X = B, where A is real symmetric or complex Hermitian and, in either
! case, positive definite, and where X and B are rectangular matrices or
! vectors. A is stored in packed format. The Cholesky decomposition is 
! used to factor A as
!      A = U^H*U if UPLO = 'U', or A = L*L^H if UPLO = 'L'
! where U is an upper triangular matrix and L is a lower triangular 
! matrix (L = U^H ). The factored form of A is then used to solve the 
! above system.
! 
! =========
! 
!       SUBROUTINE LA_PPSV( AP, B, UPLO=uplo, INFO=info )
!           <type>(<wp>), INTENT(INOUT) :: AP(:), <rhs>
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!           <type> ::= REAL | COMPLEX
!           <wp>   ::= KIND(1.0) | KIND(1.0D0)
!           <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! AP      (input/output) REAL or COMPLEX array, shape (:) with size(AP) 
!         = n*(n+1)/2, where n is the order of A.
!         On entry, the upper or lower triangle of matrix A in packed
!         storage. The elements are stored columnwise as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j<=n;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for 1<=j<=i<=n.
!         On exit, the factor U or L from the Cholesky factorization 
! 	  A = U^H*U or A = L*L^H , in the same storage format as A.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = n or shape (:) with size(B) = n.
!         On entry, the matrix B.
!         On exit, the solution matrix X.
! UPLO    Optional, (input) CHARACTER(LEN=1)
!         = 'U': Upper triangle of A is stored;
!         = 'L': Lower triangle of A is stored.
!         Default value: 'U'.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, the leading minor of order i of A is not 
! 	     positive definite, so the factorization could not be 
! 	     completed and the solution could not be computed.
!         If INFO is not present and an error occurs, then the program is
! 	  terminated with an error message.
!------------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_PPSV'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LUPLO
    INTEGER :: LINFO, N, NN, NRHS
    COMPLEX(WP) :: WW
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT, REAL, INT, AIMAG
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; NN = SIZE(A); NRHS = SIZE(B,2)
    WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW)
    IF( PRESENT(UPLO) )THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
!   .. TEST THE ARGUMENTS
    IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
    ELSE IF ( N > 0 ) THEN
       CALL PPSV_F77( LUPLO, N, NRHS, A, B, N, LINFO )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO )
 END SUBROUTINE DPPSV_F95
