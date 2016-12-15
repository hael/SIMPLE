SUBROUTINE DSYGST_F95( A, B, ITYPE, UPLO, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: SYGST_F77 => LA_SYGST
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
   INTEGER, INTENT(IN), OPTIONAL :: ITYPE
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(IN) :: B(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_SYGST / LA_HEGST reduces a real symmetric-definite or complex
! Hermitian-definite generalized eigenproblem to standard form.
!
! If ITYPE = 1, the problem is A*x = lambda*B*x,
! and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
!
! If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
! B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
!
! B must have been previously factorized as U**H*U or L*L**H 
! by LA_POTRF.
!
! =======
!
!    SUBROUTINE LA_SYGST / LA_HEGST( A, B, ITYPE, UPLO, INFO )
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!       INTEGER, INTENT(IN), OPTIONAL :: ITYPE
!       INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       <type(<wp>), INTENT(IN) :: B(:,:)
!       <type(<wp>), INTENT(INOUT) :: A(:,:)
!    where
!       <wp>   ::= KIND(1.0) | KIND(1.0D0)
!       <type> ::= REAL | COMPLEX
!
! Defaults
! ========
!
! 1. If ITYPE is not present then ITYPE = 1 is assumed.
!       
! 2. If UPLO is not present then UPLO = 'U' is assumed.
!
! Arguments
! =========
!
! A       (input/output) either REAL or COMPLEX square array, 
!         shape (:,:), size(A,1) == size(A,2) >= 0.
!         On entry, the symmetric (Hermitian) matrix A.  
!            If UPLO = 'U', the upper triangular part of A contains
!               the upper triangular part of the matrix A, and the 
!               strictly lower triangular part of A is not referenced.
!            If UPLO = 'L', the lower triangular part of A contains
!               the lower triangular part of the matrix A, and the 
!               strictly upper triangular part of A is not referenced.
!         On exit, if INFO = 0, the transformed matrix, stored in the
!            same format as A.
!
! B       (input) either REAL or COMPLEX square array,
!         shape (:,:), size(B,1) == size(A,1).
!         The triangular factor from the Cholesky factorization of B,
!         as returned by LA_POTRF.
!
! ITYPE   Optional, (input) INTEGER
!         If ITYPE is present then:
!            = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
!            = 2 or 3: compute U*A*U**H or L**H*A*L.
!         otherwise ITYPE = 1 is assumed.
!
! UPLO    Optional, (input) CHARACTER*1
!         If UPLO is present then:
!            = 'U':  Upper triangle of A is stored and B is factored as
!                    U**H*U;
!            = 'L':  Lower triangle of A is stored and B is factored as
!                    L*L**H.
!         otherwise UPLO = 'U' is assumed.
!
! INFO    Optional, (output) INTEGER
!         If INFO is present:
!            = 0: successful exit
!            < 0: if INFO = -i, the i-th argument had an illegal value
!         If INFO is not present and an error occurs, then the program
!            is terminated with an error message.
!
!------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SYGST'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LUPLO
   INTEGER :: LINFO, N, LD, LITYPE
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, SIZE, MAX
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; N = SIZE(A,1); LD = MAX(1,N)
   IF( PRESENT(ITYPE) )THEN; LITYPE = ITYPE; ELSE; LITYPE = 1; END IF
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE( B, 1 ) /= N .OR. SIZE( B, 2 ) /= N )THEN; LINFO = -2
   ELSE IF( LITYPE < 1 .OR. LITYPE > 3 )THEN; LINFO = -3
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -4
   ELSE IF( N > 0 )THEN
!  .. CALL LAPACK77 ROUTINE
      CALL SYGST_F77( LITYPE, LUPLO, N, A, LD, B, LD, LINFO )
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO)
END SUBROUTINE DSYGST_F95
