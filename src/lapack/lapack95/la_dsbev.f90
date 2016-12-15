SUBROUTINE DSBEV_F95( A, W, UPLO, Z, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: SBEV_F77 => LA_SBEV
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(OUT) :: W(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_SBEV and LA_SBEVD compute all eigenvalues and, optionally, all
! eigenvectors of a real symmetric matrix A in band form.
!     LA_HBEV and LA_HBEVD compute all eigenvalues and, optionally, all
! eigenvectors of a complex Hermitian matrix A in band form.
!     LA_SBEVD and LA_HBEVD use a divide and conquer algorithm. They are
! much faster than LA_SBEV and LA_HBEV for large matrices but use more 
! workspace.
! 
! =========
! 
!         SUBROUTINE LA_SBEV / LA_HBEV / LA_SBEVD /
!                  LA_HBEVD( AB, W, UPLO=uplo, Z=z, INFO=info )
!             <type>(<wp>), INTENT(INOUT) :: AB(:,:)
!             REAL(<wp>), INTENT(OUT) :: W(:)
!             CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!             <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!             INTEGER, INTENT(OUT), OPTIONAL :: INFO
!         where
!             <type> ::= REAL | COMPLEX
!             <wp> ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! AB     (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(AB,1) = kd + 1 and
!        size(AB,2) = n, where kd is the number of subdiagonals or 
!        superdiagonals in the band and n is the order of A.
!        On entry, the upper (if UPLO = 'U') or lower (if UPLO = 'L') 
!        triangle of matrix A in band storage. The kd + 1 diagonals of A
!        are stored in the rows of AB so that the j-th column of A is 
!        stored in the j-th column of AB as follows:
!         if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j
!                                                    1<=j<=n
!         if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd)
!                                                    1<=j<=n.
!        On exit, AB is overwritten by values generated during the 
!        reduction of A to a tridiagonal matrix T . If UPLO = 'U', the
!        first superdiagonal and the diagonal of T are returned in rows
!        kd and kd + 1 of AB. If UPLO = 'L', the diagonal and first 
!        subdiagonal of T are returned in the first two rows of AB.
! W      (output) REAL array, shape (:) with size(W) = n.
!        The eigenvalues in ascending order.
! UPLO   Optional (input) CHARACTER(LEN=1).
!        = 'U': Upper triangle of A is stored;
!        = 'L': Lower triangle of A is stored.
!        Default value: 'U'.
! Z      Optional (output) REAL or COMPLEX square array, shape (:,:) with
!        size(Z,1) = n.
!        The columns of Z contain the orthonormal eigenvectors of A in 
!        the order of the eigenvalues.
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: if INFO = i, then i off-diagonal elements of an intermediate
!             tridiagonal form did not converge to zero.
!        If INFO is not present and an error occurs, then the program is
!        terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_SBEV'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LUPLO, LJOBZ
   INTEGER :: N, KD, LINFO, LD, ISTAT, ISTAT1, S1Z, S2Z
!  .. LOCAL ARRAYS ..
   REAL(WP), TARGET :: LLZ(1,1)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; KD = SIZE(A,1)-1; N = SIZE(A,2); LD = MAX(1,SIZE(A,1))
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2); LJOBZ = 'V'
   ELSE; S1Z = 1; S2Z = 1; LJOBZ = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( KD < 0 .OR. N < 0 ) THEN; LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= N .OR. S2Z /= N ) )THEN; LINFO = -4
   ELSE IF( N > 0 )THEN
!  .. DETERMINE THE WORKSPACE
      ALLOCATE(WORK(MAX(1,3*N-2)), STAT=ISTAT)
      IF( ISTAT == 0 ) THEN
         IF( PRESENT(Z) )THEN
   !     .. CALL LAPACK77 ROUTINE
            CALL SBEV_F77( LJOBZ, LUPLO, N, KD, A, LD, W, Z, S2Z, WORK, LINFO )
	 ELSE
	    CALL SBEV_F77( LJOBZ, LUPLO, N, KD, A, LD, W, LLZ, S2Z, WORK, LINFO )
	 ENDIF   
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSBEV_F95
