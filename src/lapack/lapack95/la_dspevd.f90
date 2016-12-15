SUBROUTINE DSPEVD_F95( A, W, UPLO, Z, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: SPEVD_F77 => LA_SPEVD
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:)
   REAL(WP), INTENT(OUT) :: W(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_SPEV and LA_SPEVD compute all eigenvalues and, optionally, all
! eigenvectors of a real symmetric matrix A in packed storage.
!     LA_HPEV and LA_HPEVD compute all eigenvalues and, optionally, all
! eigenvectors of a complex Hermitian matrix A in packed storage.
!     LA_SPEVD and LA_HPEVD use a divide and conquer algorithm. If 
! eigenvectors are desired, they can be much faster than LA_SPEV and 
! LA_HPEV for large matrices but use more workspace.
! 
! =========
! 
!       SUBROUTINE LA_SPEV / LA_HPEV / LA_SPEVD / LA_HPEVD( AP, W, &
!                     UPLO=uplo, Z=z, INFO=info )
!           <type>(<wp>), INTENT(INOUT) :: AP(:)
!           REAL(<wp>), INTENT(OUT) :: W(:)
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!           <type>(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!           <type> ::= REAL | COMPLEX
!           <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! AP      (input/output) REAL or COMPLEX array, shape (:) with size(AP)=
!         n*(n+1)/2, where n is the order of A.
!         On entry, the upper or lower triangle of matrix A in packed 
!  	  storage. The elements are stored columnwise as follows:
!         if UPLO = 'U', AP(i+(j-1)*j/2)=A(i,j) for 1<=i<=j<=n;
!         if UPLO = 'L', AP(i+(j-1)*(2*n-j)/2)=A(i,j) for 1<=j<=i<=n.
!         On exit, AP is overwritten by values generated during the 
!         reduction of A to a tridiagonal matrix T . If UPLO = 'U', the
!         diagonal and first superdiagonal of T overwrite the correspond-
!         ing diagonals of A. If UPLO = 'L', the diagonal and first
! 	  subdiagonal of T overwrite the corresponding diagonals of A.
! W       (output) REAL array, shape (:) with size(W) = n.
!         The eigenvalues in ascending order.
! UPLO    Optional (input) CHARACTER(LEN=1).
!         = 'U': Upper triangle of A is stored;
!         = 'L': Lower triangle of A is stored.
!         Default value: 'U'.
! Z       Optional (output) REAL or COMPLEX square array, shape (:,:)
!         with size(Z,1) = n.
!         The columns of Z contain the orthonormal eigenvectors of A in
! 	  the order of the eigenvalues.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value
!         > 0: if INFO = i, then i off-diagonal elements of an 
! 	  intermediate tridiagonal form did not converge to zero.
!         If INFO is not present and an error occurs, then the program is
! 	  terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SPEVD'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LUPLO, LJOBZ
   INTEGER :: N, NN, LINFO, LD, ISTAT, ISTAT1, S1Z, S2Z, LWORK, LIWORK
   INTEGER, SAVE :: LWORKN = 0, LIWORKN = 0, LWORKV = 0, LIWORKV = 0
   COMPLEX(WP) :: WW
!  .. LOCAL ARRAYS ..
   REAL(WP), TARGET :: LLZ(1,1)
   INTEGER, POINTER :: IWORK(:)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; NN = SIZE(A)
   WW = (-1+SQRT(1+8*REAL(NN,WP)))*0.5; N = INT(WW);  LD = MAX(1,N)
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
   IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2); LJOBZ = 'V'
   ELSE; S1Z = 1; S2Z = 1; LJOBZ = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( NN < 0 .OR. AIMAG(WW) /= 0 .OR. REAL(N,WP) /= REAL(WW) ) THEN; LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= LD .OR. S2Z /= N ) )THEN; LINFO = -4
   ELSE IF( N > 0 )THEN
!  .. DETERMINE THE WORKSPACE
      IF( LSAME(LJOBZ,'N') )THEN
         LWORK = MAX( 1, 2*N, LWORKN ); LIWORK = MAX( 1, LIWORKN )
      ELSE
         LWORK = MAX( 1+ 6*N+N**2, LWORKV )
         LIWORK = MAX( 3+5*N, LIWORKV )
      END IF
      ALLOCATE(WORK(LWORK), IWORK(LIWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN
         DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
         IF( LSAME(LJOBZ,'N') )THEN; LWORK = MAX( 1, 2*N ); LIWORK = 1
         ELSE
	    LWORK = 1+ 6*N+N**2
	    LIWORK = 3+5*N; END IF
         ALLOCATE(WORK(LWORK), IWORK(LIWORK), STAT=ISTAT)
         IF( ISTAT == 0 ) CALL ERINFO( -200, SRNAME, LINFO )
      END IF
      IF( ISTAT == 0 ) THEN
         IF( PRESENT(Z) )THEN
           CALL SPEVD_F77( LJOBZ, LUPLO, N, A, W, Z, S2Z, WORK, LWORK, &
                         IWORK, LIWORK, LINFO )
	 ELSE
  	   CALL SPEVD_F77( LJOBZ, LUPLO, N, A, W, LLZ, S2Z, WORK, LWORK, &
                         IWORK, LIWORK, LINFO )
	 ENDIF		 
         IF (LINFO == 0 ) THEN
            IF (LSAME(LJOBZ,'N')) THEN
               LWORKN = INT(WORK(1)); LIWORKN = IWORK(1)
            ELSE; LWORKV = INT(WORK(1)); LIWORKV = IWORK(1); END IF
         END IF
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK,IWORK,STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DSPEVD_F95
