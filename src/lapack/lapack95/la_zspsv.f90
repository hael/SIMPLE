 SUBROUTINE ZSPSV_F95( AP, B, UPLO, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME
    USE F77_LAPACK, ONLY: SPSV_F77 => LA_SPSV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    COMPLEX(WP), INTENT(INOUT) :: AP(:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_SPSV computes the solution to a linear system of equations
! A*X = B, where A is a real or complex symmetric matrix stored in packed
! format and X and B are rectangular matrices or vectors. A diagonal
! pivoting method is used to factor A as
!    A = U*D*U^T if UPLO = 'U', or A = L*D*L^T if UPLO = 'L'
! where U (or L) is a product of permutation and unit upper (or lower)
! triangular matrices, and D is a symmetric block diagonal matrix with 
! 1 by 1 and 2 by 2 diagonal blocks. The factored form of A is then used 
! to solve the above system.
!     LA_HPSV computes the solution to a linear system of equations 
! A*X = B, where A is a complex Hermitian matrix stored in packed format
! and X and B are rectangular matrices or vectors. A diagonal pivoting 
! method is used to factor A as
!     A = U*D*U^H if UPLO = 'U', or A = L*D*L^H if UPLO = 'L'
! where U (or L) is a product of permutation and unit upper (or lower)
! triangular matrices, and D is a complex Hermitian block diagonal matrix
! with 1 by 1 and 2 by 2 diagonal blocks. The factored form of A is then 
! used to solve the above system.
! 
! =========
! 
!      SUBROUTINE LA_SPSV / LA_HESV( AP, B, UPLO=uplo, &
!                                 IPIV=ipiv, INFO=info )
!          <type>(<wp>), INTENT(INOUT) :: AP(:), <rhs>
!          CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!          INTEGER, INTENT(OUT), OPTIONAL :: IPIV(:)
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!      where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
!          <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! AP       (input/output) REAL or COMPLEX array, shape (:) with size(AP)=
!          n*(n + 1)=2, where n is the order of A.
!          On entry, the upper or lower triangle of matrix A in packed 
! 	   storage. The elements are stored columnwise as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j<=n;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for 1<=j<=i<=n.
!          On exit, the block diagonal matrix D and the multipliers used
! 	   to obtain U or L from the factorization of A, stored as a
! 	   packed triangular matrix in the same storage format as A.
! B        (input/output) REAL or COMPLEX array, shape (:,:) with 
!          size(B,1) = n or shape (:) with size(B) = n.
!          On entry, the matrix B.
!          On exit, the solution matrix X .
! UPLO     Optional (input) CHARACTER(LEN=1)
!             = 'U': Upper triangle of A is stored;
!             = 'L': Lower triangle of A is stored.
!          Default value: 'U'.
! IPIV     Optional (output) INTEGER array, shape (:) with size(IPIV)=n.
!          Details of the row and column interchanges and the block 
! 	   structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were 
! 	   interchanged, and D(k,k) is a 1 by 1 diagonal block.
!          If IPIV k < 0, then there are two cases:
!            1. If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows 
! 	      and columns (k-1) and -IPIV(k) were interchanged and
! 	      D(k-1:k,k-1:k) is a 2 by 2 diagonal block.
!            2. If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0, then rows 
! 	      and columns (k + 1) and -IPIV(k) were interchanged and
! 	      D(k:k+1,k:k+1) is a 2 by 2 diagonal block.
! INFO     Optional (output) INTEGER.
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, D(i,i) = 0. The factorization has been 
! 	     completed, but the block diagonal matrix D is singular,
! 	     so the solution could not be computed.
!          If INFO is not present and an error occurs, then the program 
! 	   is terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_SPSV'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LUPLO
    INTEGER :: LINFO, N, NN, NRHS, SIPIV, ISTAT, ISTAT1
    COMPLEX(WP) :: WW
!   .. LOCAL POINTERS ..
    INTEGER, POINTER :: LPIV(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT, REAL, INT, AIMAG
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; NN = SIZE(AP); NRHS = SIZE(B,2)
    WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW)
    IF( PRESENT(UPLO) )THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
    IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = N; END IF
!   .. TEST THE ARGUMENTS
    IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
    ELSE IF( SIPIV /= N )THEN; LINFO = -4
    ELSE IF ( N > 0 ) THEN
      IF( PRESENT(IPIV) )THEN; LPIV => IPIV
      ELSE; ALLOCATE( LPIV(N), STAT = ISTAT ); END IF
      IF( ISTAT == 0 ) THEN
         CALL SPSV_F77( LUPLO, N, NRHS, AP, LPIV, B, N, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(IPIV) )DEALLOCATE(LPIV, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE ZSPSV_F95
