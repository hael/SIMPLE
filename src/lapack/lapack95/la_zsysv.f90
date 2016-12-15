 SUBROUTINE ZSYSV_F95( A, B, UPLO, IPIV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO, LSAME
    USE F77_LAPACK, ONLY: SYSV_F77 => LA_SYSV, ILAENV_F77 => ILAENV
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!   .. ARRAY ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_SYSV computes the solution to a linear system of equations
! A*X = B, where A is a real or complex symmetric matrix and X and B are
! rectangular matrices or vectors. A diagonal pivoting method is used to
! factor A as
!      A = U*D*U^T if UPLO = 'U', or A = L*D*L^T if UPLO = 'L'
! where U (or L) is a product of permutation and unit upper (or lower)
! triangular matrices, and D is a symmetric block diagonal matrix with 
! 1 by 1 and 2 by 2 diagonal blocks. The factored form of A is then used
! to solve the above system.
!    LA_HESV computes the solution to a linear system of equations 
! A*X = B, where A is a complex Hermitian matrix and X and B are 
! rectangular matrices or vectors. A diagonal pivoting method is used to
! factor A as
!      A = U*D*U^H if UPLO = 'U', or A = L*D*L^H if UPLO = 'L'
! where U (or L) is a product of permutation and unit upper (or lower)
! triangular matrices, and D is a complex Hermitian block diagonal
! matrix with 1 by 1 and 2 by 2 diagonal blocks. The factored form of A
! is then used to solve the above system.
! 
! =========
! 
!          SUBROUTINE LA_SYSV / LA_HESV( A, B, UPLO=uplo, &
! 	                            IPIV=ipiv, INFO=info )
!                <type>(<wp>), INTENT(INOUT) :: A(:,:), <rhs>
!                CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!                INTEGER, INTENT(OUT), OPTIONAL :: IPIV(:)
!                INTEGER, INTENT(OUT), OPTIONAL :: INFO
!          where
!                <type> ::= REAL | COMPLEX
!                <wp>   ::= KIND(1.0) | KIND(1.0D0)
!                <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX square array, shape (:,:).
!        On entry, the matrix A.
!        If UPLO = 'U', the upper triangular part of A contains the upper
!        triangular part of the matrix A, and the strictly lower 
!        triangular part of A is not referenced.
!        If UPLO = 'L', the lower triangular part of A contains the lower 
!        triangular part of the matrix A, and the strictly upper 
!        triangular part of A is not referenced.
!        On exit, the block diagonal matrix D and the multipliers used to
!        obtain the factor U or L from the factorization of A.
! B      (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(B,1) = size(A,1) or shape (:) with size(B) = size(A,1).
!        On entry, the matrix B.
!        On exit, the solution matrix X.
! UPLO   Optional (input) CHARACTER(LEN=1)
!          = 'U': Upper triangle of A is stored;
!          = 'L': Lower triangle of A is stored.
!        Default value: 'U'.
! IPIV   Optional (output) INTEGER array, shape (:) with size(IPIV) = 
!        size(A,1).
!        Details of the row and column interchanges and the block 
!        structure of D.
!        If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!        interchanged, and D(k,k) is a 1 by 1 diagonal block.
!        If IPIV k < 0, then there are two cases:
!         1. If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and 
! 	   columns (k-1) and -IPIV(k) were interchanged and 
! 	   D(k-1:k,k-1:k) is a 2 by 2 diagonal block.
!         2. If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0, then rows and
! 	   columns (k + 1) and -IPIV(k) were interchanged and 
! 	   D(k:k+1,k:k+1) is a 2 by 2 diagonal block.
! INFO   Optional (output) INTEGER
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: if INFO = i, D(i,i) = 0. The factorization has been 
!             completed, but the block diagonal matrix D is singular, so
! 	    the solution could not be computed.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!-----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_SYSV'
    CHARACTER(LEN=6), PARAMETER :: BSNAME = 'ZSYTRF'
!   .. LOCAL SCALARS ..
    CHARACTER(LEN=1) :: LUPLO
    INTEGER :: LINFO, ISTAT, ISTAT1, SIPIV, N, NRHS, LWORK, NB
!   .. LOCAL POINTERS ..
    INTEGER, POINTER :: LPIV(:)
    COMPLEX(WP), POINTER :: WORK(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; N = SIZE(A,1); NRHS = SIZE(B,2)
    IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
    IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = SIZE(A,1); END IF
!   .. TEST THE ARGUMENTS
    IF( SIZE( A, 2 ) /= N .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= N .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
    ELSE IF( SIPIV /= N )THEN; LINFO = -4
    ELSE IF ( N > 0 ) THEN
!  .. DETERMINE THE WORKSPACE
      IF( PRESENT(IPIV) )THEN; LPIV => IPIV
      ELSE; ALLOCATE( LPIV(N), STAT = ISTAT ); END IF
      IF( ISTAT == 0 )THEN
         NB = ILAENV_F77( 1, BSNAME, LUPLO, N, -1, -1, -1 )
         IF( NB <= 1 .OR. NB >= N ) NB = 1; LWORK = N*NB
         ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT /= 0 )THEN
            DEALLOCATE(WORK, STAT=ISTAT1); LWORK = 3*N
            ALLOCATE(WORK(LWORK), STAT=ISTAT)
            IF( ISTAT /= 0 ) THEN; LINFO = - 100
            ELSE; CALL ERINFO( -200, SRNAME, LINFO ); ENDIF
         ENDIF
         IF ( ISTAT == 0 ) &
!           .. CALL LAPACK77 ROUTINE
            CALL SYSV_F77( LUPLO, N, NRHS, A, N, LPIV, B, N, WORK, LWORK, LINFO )
      ELSE; LINFO = -100; END IF
      IF( .NOT.PRESENT(IPIV) )DEALLOCATE(LPIV, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE ZSYSV_F95
