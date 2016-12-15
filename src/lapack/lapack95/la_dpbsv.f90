 SUBROUTINE DPBSV_F95( A, B, UPLO, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME
    USE F77_LAPACK, ONLY: PBSV_F77 => LA_PBSV
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
!     LA_PBSV computes the solution to a linear system of equations
! A*X = B, where A has band form and is real symmetric or complex
! Hermitian and, in either case, positive definite, and where X and B are
! rectangular matrices or vectors. The Cholesky decomposition is used to
! factor A as A = U^H*U if UPLO = 'U', or A = L*L^H if UPLO = 'L'
! where U is an upper triangular band matrix and L is a lower triangular 
! band matrix, each with the same number of superdiagonals or subdiagonals
! as A. The factored form of A is then used to solve the above system.
! 
! =========
! 
!      SUBROUTINE LA_PBSV( AB, B, UPLO=uplo, INFO=info )
!          <type>(<wp>), INTENT(INOUT) :: AB(:,:), <rhs>
!          CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!      where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
!          <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! AB      (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(AB,1) = kd + 1 and size(AB,2) = n, where kd is the number
!         of superdiagonals or subdiagonals in the band and n is the 
!  	  order of A.
!         On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') 
! 	  triangle of matrix A in band storage. The (kd + 1) diagonals of
! 	  A are stored in the rows of AB so that the j-th column of A is
! 	  stored in the j-th column of AB as follows:
! 	  if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j
! 	                                           1<=j<=n
!         if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd)
! 	                                           1<=j<=n.
!         On exit, the factor U or L from the Cholesky factorization 
!    	  A = U^H*U = L*L^H in the same storage format as A.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = n or shape (:) with size(B) = n.
!         On entry, the matrix B.
!         On exit, the solution matrix X .
! UPLO    Optional (input) CHARACTER(LEN=1)
!           = 'U': Upper triangle of A is stored;
!           = 'L': Lower triangle of A is stored.
!         Default value: 'U'.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, the leading minor of order i of A is not 
! 	     positive definite, so the factorization could not be 
! 	     completed and the solution could not be computed.
!         If INFO is not present and an error occurs, then the program 
! 	  is terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_PBSV'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LUPLO
    INTEGER :: LINFO, KD, N, NRHS
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; KD = SIZE(A,1)-1; N = SIZE(A,2); NRHS = SIZE(B,2)
    IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; ENDIF
!   .. TEST THE ARGUMENTS
    IF( KD < 0 .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
    ELSE IF ( N > 0 ) THEN
       CALL PBSV_F77( LUPLO, N, KD, NRHS, A, KD+1, B, N, LINFO )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO )
 END SUBROUTINE DPBSV_F95
