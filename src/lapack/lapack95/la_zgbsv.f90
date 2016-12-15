 SUBROUTINE ZGBSV_F95( A, B, KL, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO
    USE F77_LAPACK, ONLY: GBSV_F77 => LA_GBSV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    INTEGER, INTENT(IN), OPTIONAL :: KL
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_GBSV computes the solution to a real or complex linear system
! of equations A*X = B, where A is a square band matrix and X and B are
! rectangular matrices or vectors. The LU decomposition with row
! interchanges is used to factor A as A = L*U , where L is a product of
! permutation and unit lower triangular matrices with kl subdiagonals, 
! and U is upper triangular with kl + ku superdiagonals. The factored 
! form of A is then used to solve the above system.
! 
! =========
! 
!       SUBROUTINE LA_GBSV( AB, B, KL=kl, IPIV=ipiv, INFO=info )
!             <type>(<wp>), INTENT(INOUT) :: AB(:,:), <rhs>
!             INTEGER, INTENT(IN), OPTIONAL :: KL
!             INTEGER, INTENT(OUT), OPTIONAL :: IPIV(:)
!             INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!             <type> ::= REAL | COMPLEX
!             <wp>   ::= KIND(1.0) | KIND(1.0D0)
!             <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! AB      (input/output) REAL or COMPLEX rectangular array, shape (:,:)
!         with size(AB,1) = 2*kl+ku+1 and size(AB,2) = n, where kl and ku
!         are, respectively, the numbers of subdiagonals and 
!         superdiagonals in the band of A, and n is the order of A.
!         On entry, the matrix A in band storage. The (kl + ku + 1) 
! 	  diagonals of A are stored in rows (kl + 1) to (2*kl + ku + 1) 
! 	  of AB, so that the j-th column of A is stored in the j-th 
! 	  column of AB as follows:
! 	  AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl)
! 	                                 1<=j<=n
!         The remaining elements in AB need not be set.
!         On exit, details of the factorization. U is an upper triangular 
!         band matrix with (kl + ku + 1) diagonals. These are stored in
! 	  the first (kl + ku + 1) rows of AB. The multipliers that arise
!         during the factorization are stored in the remaining rows.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = n or shape (:) with size(B) = n.
!         On entry, the matrix B.
!         On exit, the solution matrix X.
! KL      Optional (input) INTEGER.
!         The number of subdiagonals in the band of A (KL = kl).
!         The number of superdiagonals in the band is given by
! 	  ku = size(AB,1) - 2 * kl - 1.
!         Default value: (size(AB,1)-1)/3.
! IPIV    Optional (output) INTEGER array, shape (:) with size(IPIV) = n.
!         The pivot indices that define the row interchanges; row i of the
!         matrix was interchanged with row IPIV(i).
! INFO    Optional (output) INTEGER
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, U(i,i) = 0. The factorization has been 
!         completed, but the factor U is singular, so the solution could 
!         not be computed. 
!         If INFO is not present and an error occurs, then the program 
!         is terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GBSV'
!   .. LOCAL SCALARS ..
    INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV, LDA, N, NRHS, LKL, KU
!   .. LOCAL POINTERS ..
    INTEGER, POINTER :: LPIV(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0
    LDA = SIZE(A,1); N = SIZE(A,2); NRHS = SIZE(B,2)
    IF( PRESENT(KL) ) THEN; LKL = KL; ELSE; LKL = (LDA-1)/3; ENDIF
    IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = N; ENDIF
!   .. TEST THE ARGUMENTS
    IF( LDA - 2*LKL -1 < 0 .OR. LDA < 0 .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( LDA - 2*LKL -1 < 0 .OR. LKL < 0 ) THEN; LINFO = -3
    ELSE IF( SIPIV /= N )THEN; LINFO = -4
    ELSE IF ( N > 0 ) THEN
       IF( PRESENT(IPIV) )THEN; LPIV => IPIV; ELSE
           ALLOCATE( LPIV(N), STAT = ISTAT ); END IF
       IF ( ISTAT == 0 ) THEN
          KU = LDA -2*LKL -1
          CALL GBSV_F77( N, LKL, KU, NRHS, A, LDA, LPIV, B, N, LINFO )
       ELSE
          LINFO = -100
       END IF
       IF( .NOT.PRESENT(IPIV) )THEN
          DEALLOCATE(LPIV, STAT = ISTAT1 )
       END IF
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE ZGBSV_F95
