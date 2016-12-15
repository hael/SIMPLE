      SUBROUTINE DGETRF_F95( A, IPIV, RCOND, NORM, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: LSAME, ERINFO
      USE F77_LAPACK, ONLY: GETRF_F77 => LA_GETRF, LANGE_F77 => LA_LANGE, &
                                          GECON_F77 => LA_GECON
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
      REAL(WP), INTENT( OUT ), OPTIONAL :: RCOND
!  .. ARRAY ARGUMENTS ..
      INTEGER, INTENT( OUT ), OPTIONAL, TARGET :: IPIV( : )
      REAL(WP), INTENT( INOUT ) :: A( :, : )
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_GESV computes the solution to a real or complex linear system of
! equations A*X = B, where A is a square matrix and X and B are 
! rectangular matrices or vectors. Gaussian elimination with row 
! interchanges is used to factor A as A = P*L*U , where P is a permutation
! matrix, L is unit lower triangular, and U is upper triangular. The 
! factored form of A is then used to solve the above system.
! 
! =========
! 
!       SUBROUTINE LA_GESV( A, B, IPIV=ipiv, INFO=info )
!          <type>(<wp>), INTENT(INOUT) :: A(:,:), <rhs>
!          INTEGER, INTENT(OUT), OPTIONAL :: IPIV(:)
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
!          <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX square array, shape (:,:).
!        On entry, the matrix A.
!        On exit, the factors L and U from the factorization A = P*L*U; 
!        the unit diagonal elements of L are not stored.
! B      (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(B,1) = size(A,1) or shape (:) with size(B) = size(A,1).
!        On entry, the matrix B.
!        On exit, the solution matrix X .
! IPIV   Optional (output) INTEGER array, shape (:) with size(IPIV) = 
!        size(A,1).
!        The pivot indices that define the permutation matrix P; row i of
!        the matrix was interchanged with row IPIV(i).
! INFO   Optional (output) INTEGER
!        = 0 : successful exit.
!        < 0 : if INFO = -i, the i-th argument has an illegal value.
!        > 0 : if INFO = i, then U(i,i) = 0. The factorization has been 
!        completed, but the factor U is singular, so the solution could
!        not be computed.
!        If INFO is not present and an error occurs, then the program is 
!        terminated with an error message.
!----------------------------------------------------------------------
!  .. PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GETRF'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LNORM
      INTEGER :: LINFO, M, N, LD, ISTAT, ISTAT1, MINMN, LWORK, SIPIV
      REAL(WP) :: LANORM
!  .. LOCAL POINTERS ..
      INTEGER, POINTER :: LIPIV(:), IWORK(:)
      REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC PRESENT, MAX, MIN, SIZE, TINY
!  .. EXECUTABLE STATEMENTS ..
      M = SIZE(A,1); N = SIZE(A,2); LINFO = 0; ISTAT = 0; MINMN = MIN(M,N)
      LD = MAX(1,M)
      IF( PRESENT(IPIV) )THEN; SIPIV = SIZE(IPIV); ELSE; SIPIV = MINMN; ENDIF
      IF ( PRESENT(NORM) ) THEN; LNORM = NORM; ELSE; LNORM = '1'; END IF
!  .. TEST THE ARGUMENTS
      IF( M < 0 .OR. N < 0 .OR. PRESENT(RCOND) .AND. M /= N )THEN; LINFO = -1
      ELSE IF( SIPIV /= MINMN )THEN; LINFO = -2
      ELSE IF( ( .NOT.PRESENT(RCOND) .AND. PRESENT(NORM) ) .OR. &
               ( .NOT.LSAME(LNORM,'I') .AND. .NOT.LSAME(LNORM,'O') &
                                       .AND. LNORM /= '1' ) ) THEN; LINFO = -4
      ELSE IF( M > 0 .AND. N > 0 ) THEN
         IF( PRESENT(RCOND) .AND. M == N ) THEN
!        .. COMPUTE THE NORM OF THE MATRIX A
            IF( LNORM == 'I' ) THEN; LWORK = MINMN; ELSE; LWORK = 1; END IF
            ALLOCATE( WORK(LWORK), STAT=ISTAT )
            IF( ISTAT == 0 )THEN
               LANORM = LANGE_F77( LNORM, MINMN, MINMN, A, LD, WORK )
            ELSE
               LINFO = -100
            END IF
            DEALLOCATE(WORK, STAT=ISTAT1)
         END IF
         IF( LINFO == 0 ) THEN
            IF( PRESENT(IPIV) )THEN; LIPIV => IPIV
            ELSE; ALLOCATE( LIPIV( MINMN ), STAT=ISTAT ); ENDIF
            IF( ISTAT /= 0 )LINFO = -100
         END IF
         IF( LINFO == 0 ) THEN
!           .. COMPUTE THE LU FACTORS OF THE MATRIX A
            CALL GETRF_F77( M, N, A, LD, LIPIV, LINFO )
            IF( .NOT. PRESENT(IPIV) )DEALLOCATE( LIPIV, STAT=ISTAT )
            IF( PRESENT(RCOND) ) THEN
!              .. COMPUTE THE RECIPROCAL OF THE CONDITION NUMBER OF A
               IF( LANORM <= TINY(1.0_WP) .OR. M /= N .OR. LINFO /= 0 ) THEN
                  RCOND = 0.0_WP
               ELSE 
                  ALLOCATE(WORK(4*MINMN), IWORK(MINMN), STAT=ISTAT)
                  IF( ISTAT == 0 )THEN
                     CALL GECON_F77( LNORM, MINMN, A, LD, LANORM, &
                                     RCOND, WORK, IWORK, LINFO )
                  ELSE
                     LINFO = -100
                  END IF
                  DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
               END IF
            END IF
         END IF
      ELSE IF( PRESENT(RCOND) ) THEN
         IF( M == N )THEN
            RCOND = 1.0_WP
         ELSE
            RCOND = 0.0_WP
         END IF
      END IF
      CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
      END SUBROUTINE DGETRF_F95
