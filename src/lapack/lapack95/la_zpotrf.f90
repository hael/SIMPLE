SUBROUTINE ZPOTRF_F95( A, UPLO, RCOND, NORM, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: POTRF_F77 => LA_POTRF, LANSY_F77 => LA_LANSY, &
                                  POCON_F77 => LA_POCON
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM, UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL  :: INFO
   REAL(WP), INTENT(OUT), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_POTRF computes the Cholesky factorization of a real symmetric or
! complex Hermitian positive definite matrix A.
!
! The factorization has the form
!    A = U**H * U,  if UPLO = 'U', or
!    A = L * L**H,  if UPLO = 'L',
! where U is an upper triangular matrix and L is lower triangular.
!
! This is the block version of the algorithm, calling Level 3 BLAS.
!
! LA_POTRF optionally estimates the reciprocal of the condition number
! (in the 1-norm) of a real symmetric or complex Hermitian positive 
! definite matrix A.
! An estimate is obtained for norm(inv(A)), and the reciprocal of the
! condition number is computed as RCOND = 1 / (norm(A) * norm(inv(A))).
!
! =======
!
!    SUBROUTINE LA_POTRF( A, UPLO, RCOND, NORM, INFO )
!       <type>(<wp>), INTENT(INOUT) :: A(:,:)
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!       REAL(<wp>), INTENT(OUT), OPTIONAL :: RCOND
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
!       INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!       <type> ::= REAL | COMPLEX
!       <wp>   ::= KIND(1.0) | KIND(1.0D0)
!
! Defaults
! ========
!
! 1. If UPLO is not present then UPLO = 'U' is assumed.
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
!         On exit, if INFO = 0, the factor U or L from the Cholesky
!            factorization A = U**H*U or A = L*L**H.
!
! UPLO    Optional, (input) CHARACTER*1
!         If UPLO is present then:
!            = 'U':  Upper triangle of A is stored;
!            = 'L':  Lower triangle of A is stored.
!         otherwise UPLO = 'U' is assumed.
!
! RCOND   Optional (output) REAL
!         The reciprocal of the condition number of the matrix A 
!         computed as RCOND = 1/(norm(A) * norm(inv(A))).
! NORM    Optional (input) CHARACTER*1
!         Specifies whether the 1-norm condition number or the
!         infinity-norm condition number is required:
!           If NORM is present then:
!              = '1', 'O' or 'o': 1-norm;
!              = 'I' or 'i': infinity-norm.
!           otherwise NORM = '1' is used.
!
! INFO    Optional, (output) INTEGER
!         If INFO is present:
!            = 0: successful exit
!            < 0: if INFO = -i, the i-th argument had an illegal value
!            > 0: if INFO = i, the leading minor of order i is not
!               positive definite, and the factorization could not be
!               completed.
!         If INFO is not present and an error occurs, then the program
!            is terminated with an error message.
!
! --------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_POTRF'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LNORM, LUPLO
   INTEGER :: LINFO, N, ISTAT, ISTAT1, LD
   REAL(WP) :: ANORM
!  .. LOCAL POINTERS ..
   REAL(WP), POINTER :: RWORK(:)
   COMPLEX(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC PRESENT, MAX
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; N = SIZE(A,1); LD = MAX(1,N); ISTAT = 0
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(NORM) ) THEN; LNORM = NORM; ELSE; LNORM = '1'; END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .AND. N < 0 )THEN; LINFO = -1
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -2
   ELSE IF( ( .NOT.PRESENT(RCOND) .AND. PRESENT(NORM) ) .OR. &
            ( .NOT.LSAME(LNORM,'I') .AND. .NOT.LSAME(LNORM,'O') &
              .AND. LNORM /= '1' ) ) THEN; LINFO = -4
   ELSE IF(  N > 0 )THEN
      IF( PRESENT(RCOND) ) THEN
!     .. COMPUTE THE NORM OF THE MATRIX A
         ALLOCATE(RWORK(N), STAT=ISTAT)
         IF( ISTAT == 0 )THEN; ANORM = LANSY_F77( LNORM, LUPLO, LD, A, N, RWORK )
         ELSE; LINFO = -100; END IF
         DEALLOCATE(RWORK, STAT=ISTAT1)
      END IF
!
      IF( LINFO == 0 ) THEN
!     .. COMPUTE THE CHOLESKY FACTORS OF THE MATRIX A
         CALL POTRF_F77( LUPLO, N, A, LD, LINFO )
!
         IF( PRESENT(RCOND) .AND. LINFO == 0 ) THEN
!        .. COMPUTE THE RECIPROCAL OF THE CONDITION NUMBER OF A
            IF( ANORM == 0.0_WP )THEN; RCOND = 0.0_WP
            ELSE; ALLOCATE(WORK(2*N), RWORK(N), STAT=ISTAT)
               IF( ISTAT == 0 )THEN
                  CALL POCON_F77( LUPLO, N, A, LD, ANORM, RCOND, &
                                  WORK, RWORK, LINFO )
               ELSE; LINFO = -100; END IF
               DEALLOCATE(WORK, RWORK, STAT=ISTAT1)
            END IF
         END IF
      END IF
   ELSE IF( PRESENT(RCOND) ) THEN; RCOND = 1.0_WP; ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZPOTRF_F95
