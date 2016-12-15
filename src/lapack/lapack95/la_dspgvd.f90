SUBROUTINE DSPGVD_F95( AP, BP, W, ITYPE, UPLO, Z, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: SPGVD_F77 => LA_SPGVD
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
   INTEGER, INTENT(IN), OPTIONAL :: ITYPE
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: AP(:), BP(:)
   REAL(WP), INTENT(OUT) :: W(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_SPGV, LA_SPGVD, LA_HPGV and LA_HPGVD compute all eigenvalues 
! and, optionally, all eigenvectors of generalized eigenvalue problems 
! of the form A*z = lambda*B*z, A*B*z = lambda*z; and B*A*z = lambda*z,
! where A and B are real symmetric in the cases of LA_SPGV and LA_SPGVD
! and complex Hermitian in the cases of LA_HPGV and LA_HPGVD. In all 
! four cases B is positive definite. Matrices A and B are stored in a 
! packed format.
!     LA_SPGVD and LA_HPGVD use a divide and conquer algorithm. If 
! eigenvectors are desired, they can be much faster than LA_SPGV and
! LA_HPGV for large matrices but use more workspace.
! 
! =========
! 
!       SUBROUTINE LA_SPGV / LA_SPGVD / LA_HPGV / LA_HPGVD( AP, BP, &
!                         W, ITYPE=itype, UPLO=uplo, Z=z, INFO=info )
!            <type>(<wp>), INTENT(INOUT) :: AP(:), BP(:)
!            REAL(<wp>), INTENT(OUT) :: W(:)
!            INTEGER, INTENT(IN), OPTIONAL :: ITYPE
!            CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!            <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!            INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!            <type> ::= REAL | COMPLEX
!            <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! AP       (input/output) REAL or COMPLEX array, shape (:) with
!          size(AP) = n*(n + 1)/2, where n is the order of A and B.
!          On entry, the upper or lower triangle of matrix A in packed 
! 	   storage. The elements are stored columnwise as follows:
! 	   if UPLO = 'U', AP(i +(j-1)*j/2) = A(i,j) for 1<=i<=j<=n;
! 	   if UPLO = 'L', AP(i +(j-1)*(2*n-j)/2) = A(i,j) for 1<=j<=i<=n.
!          On exit, the contents of AP are destroyed.
! BP       (input/output) REAL or COMPLEX array, shape (:) and 
!          size(BP) = size(AP).
!          On entry, the upper or lower triangle of matrix B in packed 
! 	   storage. The elements are stored columnwise as follows:
! 	   if UPLO = 'U', BP(i +(j-1)*j/2) = B(i,j) for 1<=i<=j<=n;
! 	   if UPLO = 'L', BP(i +(j-1)*(2*n-j)/2) = B(i,j) for 1<=j<=i<=n.
!          On exit, the triangular factor U or L of the Cholesky 
! 	   factorization B = U^H*U or B = L*L^H, in the same storage
! 	   format as B.
! W        (output) REAL array, shape (:) with size(W) = n.
!          The eigenvalues in ascending order.
! ITYPE    Optional (input) INTEGER.
!          Specifies the problem type to be solved:
!             = 1: A*z = lambda*B*z
!             = 2: A*B*z = lambda*z
!             = 3: B*A*z = lambda*z
!          Default value: 1.
! UPLO     Optional (input) CHARACTER(LEN=1).
!             = 'U': Upper triangles of A and B are stored;
!             = 'L': Lower triangles of A and B are stored.
!          Default value: 'U'.
! Z        Optional (output) REAL or COMPLEX square array, shape (:,:)
!          with size(Z,1) = n.
!          The matrix Z of eigenvectors, normalized as follows:
!             if ITYPE = 1 or 2: Z^H * B * Z = I ,
!             if ITYPE = 3: Z^H * B^-1 * Z = I .
! INFO     Optional (output) INTEGER.
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: the algorithm failed to converge or matrix B is not 
! 	      positive definite:
!              <= n: if INFO = i, i off-diagonal elements of an 
! 	          intermediate tridiagonal form did not converge to 
! 		  zero.
!              > n: if INFO = n+i, for 1<=i<=n, then the leading minor 
! 	          of order i of B is not positive definite. The 
! 		  factorization of B could not be completed and no
!                   eigenvalues or eigenvectors were computed.
!          If INFO is not present and an error occurs, then the program
! 	   is terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SPGVD'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) ::  LJOBZ, LUPLO
      INTEGER :: LINFO, N, NN, LD, LITYPE, ISTAT, S1Z, S2Z
      INTEGER :: LWORK, LIWORK
!  .. LOCAL ARRAYS ..
      REAL(WP), TARGET :: LLZ(1,1), WORKMIN(1)
      REAL(WP), POINTER :: WORK(:)
      INTEGER :: IWORKMIN(1)
      COMPLEX(WP) :: WW
      INTEGER, POINTER :: IWORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC SIZE, MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
      LINFO = 0; ISTAT = 0; NN = SIZE(AP)
      WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW);  LD = MAX(1,N)
      IF( PRESENT(ITYPE) )THEN; LITYPE = ITYPE; ELSE; LITYPE = 1; END IF
        IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2); LJOBZ = 'V'
        ELSE; S1Z = 1; S2Z = 1; LJOBZ = 'N'; END IF
          IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
!  .. TEST THE ARGUMENTS
            IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
            ELSE IF( SIZE(BP) /= SIZE(AP)  )THEN; LINFO = -2
            ELSE IF( SIZE(W) /= N )THEN; LINFO = -3
            ELSE IF( LITYPE < 1 .OR. LITYPE > 3 )THEN; LINFO = -4
            ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -5
            ELSE IF( PRESENT(Z) .AND. ( S1Z /= LD .OR. S2Z /= N ) )THEN; LINFO = -6
            ELSE IF( N > 0 )THEN
!  QUERING THE SIZE OF WORKSPACE ...
                LWORK = -1
                LIWORK = -1
		IF (PRESENT (Z)) THEN
                   CALL SPGVD_F77( LITYPE, LJOBZ, LUPLO, N, AP, BP, W, &
     &                Z, S1Z, WORKMIN, LWORK, IWORKMIN, LIWORK, LINFO )
                ELSE
		   CALL SPGVD_F77( LITYPE, LJOBZ, LUPLO, N, AP, BP, W, &
     &                LLZ, S1Z, WORKMIN, LWORK, IWORKMIN, LIWORK, LINFO )
                ENDIF
                LWORK = WORKMIN(1)
                LIWORK= IWORKMIN(1)
                ALLOCATE(WORK(LWORK), IWORK(LIWORK), STAT=ISTAT)
                IF( ISTAT == 0 )THEN
                  IF (PRESENT(Z)) THEN
                     CALL SPGVD_F77( LITYPE, LJOBZ, LUPLO, N, AP, BP, W, &
     &                    Z, S1Z, WORK, LWORK, IWORK, LIWORK, LINFO )
                  ELSE
		     CALL SPGVD_F77( LITYPE, LJOBZ, LUPLO, N, AP, BP, W, &
     &                    LLZ, S1Z, WORK, LWORK, IWORK, LIWORK, LINFO )
                  ENDIF
                ELSE; LINFO = -100; ENDIF
                ENDIF
                DEALLOCATE(WORK, IWORK, STAT=ISTAT)
                CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSPGVD_F95

                
