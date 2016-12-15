SUBROUTINE DSYEVD_F95( A, W, JOBZ, UPLO, INFO )
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: SYEVD_F77 => LA_SYEVD, ILAENV_F77 => ILAENV
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(OUT) :: W(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_SYEV and LA_SYEVD compute all eigenvalues and, optionally, all
! eigenvectors of a real symmetric matrix A.
!      LA_HEEV and LA_HEEVD compute all eigenvalues and, optionally, all
! eigenvectors of a complex Hermitian matrix A.
!      LA_SYEVD and LA_HEEVD use a divide and conquer algorithm. If 
! eigenvectors are desired, they can be much faster than LA_SYEV and 
! LA_HEEV for large matrices but use more workspace.
! 
! =========
! 
!       SUBROUTINE LA_SYEV / LA_HEEV / LA_SYEVD / LA_HEEVD( A, W, &
!                       JOBZ=jobz, UPLO=uplo, INFO=info )
!           <type>(<wp>), INTENT(INOUT) :: A(:,:)
!           REAL(<wp>), INTENT(OUT) :: W(:)
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!           <type> ::= REAL | COMPLEX
!           <wp> ::= KIND(1.0) | KIND(1.0D0)
! 
! 
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX square array, shape (:,:).
!        On entry, the matrix A.
!        If UPLO = 'U', the upper triangular part of A contains the upper
!        triangular part of the matrix A. If UPLO = 'L', the lower 
!        triangular part of A contains the lower triangular part of the
!        matrix A.
!        On exit:
!        If JOBZ = 'V', then the columns of A contain the orthonormal 
!        eigenvectors of the matrix A in the order of the eigenvalues.
!        If JOBZ = 'N', then the upper triangle (if UPLO = 'U') or the 
!        lower triangle (if UPLO = 'L') of A, including the diagonal, is
!        destroyed.
! W      (output) REAL array, shape (:) with size(W) = size(A,1).
!        The eigenvalues in ascending order.
! JOBZ   Optional (input) CHARACTER(LEN=1).
!        = 'N': Computes eigenvalues only;
!        = 'V': Computes eigenvalues and eigenvectors.
!        Default value: 'N'.
! UPLO   Optional (input) CHARACTER(LEN=1).
!        = 'U': Upper triangle of A is stored;
!        = 'L': Lower triangle of A is stored.
!        Default value: 'U'.
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value
!        > 0: if INFO = i, then i off-diagonal elements of an
!        intermediate tridiagonal form did not converge to zero.
!        If INFO is not present and an error occurs, then the program is 
!        terminated with an error message.
!-----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SYEVD'
   CHARACTER(LEN=6), PARAMETER :: BSNAME = 'DSYTRD'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBZ, LUPLO
   INTEGER :: N, LINFO, LD, ISTAT, ISTAT1, LIWORK, LMWORK, LWORK, NB
!  .. LOCAL ARRAYS ..
   REAL(WP), POINTER :: WORK(:)
   INTEGER, POINTER :: IWORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   N = SIZE( A, 1 ); LINFO = 0; LD = MAX(1,N); ISTAT = 0
   IF( PRESENT(JOBZ) ) THEN
      LJOBZ = JOBZ
   ELSE
      LJOBZ = 'N'
   END IF
   IF( PRESENT(UPLO) ) THEN
      LUPLO = UPLO
   ELSE
      LUPLO = 'U'
   END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN
      LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN
      LINFO = -2
   ELSE IF( .NOT.LSAME(LJOBZ,'N') .AND. .NOT.LSAME(LJOBZ,'V') )THEN
      LINFO = -3
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN
      LINFO = -4
   ELSE IF( N > 0 )THEN
!  .. DETERMINE THE WORKSPACE
      IF( LSAME(LJOBZ,'V') )THEN
         LMWORK = 1+6*N+2*N**2 
         LIWORK = 3+5*N
      ELSE
         LMWORK = 2*N+1
         LIWORK = 1
      END IF
      NB = ILAENV_F77( 1, BSNAME, LUPLO, N, -1, -1, -1 )
      IF( NB <= 1 .OR. NB >= N )THEN
         NB = 1
      END IF
      LWORK = MAX( LMWORK, (NB+2)*N )
      ALLOCATE(WORK(LWORK), IWORK(LIWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN
         DEALLOCATE(WORK, IWORK, STAT=ISTAT1)
         LWORK = LMWORK
         ALLOCATE(WORK(LWORK), IWORK(LIWORK), STAT=ISTAT)
         IF( ISTAT /= 0 ) THEN
            LINFO = - 100
         ELSE
            CALL ERINFO( -200, SRNAME, LINFO )
         ENDIF
      ENDIF
!
      IF( LINFO == 0 )THEN
!        .. CALL LAPACK77 ROUTINE
         CALL SYEVD_F77( LJOBZ, LUPLO, N, A, LD, W, WORK, LWORK, &
                         IWORK, LIWORK,  LINFO )
      ENDIF
      DEALLOCATE(WORK, IWORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSYEVD_F95
