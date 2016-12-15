SUBROUTINE DSTEVD_F95( D, E, Z, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: STEVD_F77 => LA_STEVD
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: D(:), E(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: Z(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!    LA_STEV and LA_STEVD compute all eigenvalues and, optionally, all
! eigenvectors of a real symmetric tridiagonal matrix A.
!    LA_STEVD uses a divide and conquer algorithm. If eigenvectors are 
! desired, they can be much faster than LA_STEV for large matrices but
! uses more workspace.
! =========
! 
!      SUBROUTINE LA_STEV / LA_STEVD( D, E, Z=z, INFO=info )
!         REAL(<wp>), INTENT(INOUT) :: D(:), E(:)
!         REAL(<wp>), INTENT(OUT), OPTIONAL :: Z(:,:)
!         INTEGER, INTENT(OUT), OPTIONAL :: INFO
!      where
!         <wp> ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! D    (input/output) REAL array shape (:) with size(D) = n, where n is 
!      the order of A.
!      On entry, the diagonal elements of the matrix A.
!      On exit, the eigenvalues in ascending order.
! E    (input/output) REAL array, shape (:) with size(E) = n.
!      On entry, the n - 1 subdiagonal elements of A in E(1) to E(n-1).
!      E(n) need not be set but is used by the routine.
!      On exit, the contents of E are destroyed.
! Z    Optional (output) REAL square array, shape(:,:) with size(Z,1)=n.
!      The columns of Z contain the orthonormal eigenvectors of A in the
!      order of the eigenvalues.
! INFO Optional (output) INTEGER.
!      = 0: successful exit.
!      < 0: if INFO = -i, the i-th argument had an illegal value.
!      > 0: if INFO = i, then i elements of E did not converge to zero.
!      If INFO is not present and an error occurs, then the program is 
!      terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_STEVD'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBZ
   INTEGER :: N, LD, ISTAT1, S1Z, S2Z, LWORK, LIWORK
   INTEGER :: LINFO, ISTAT
   INTEGER, SAVE :: LWORKN = 0, LIWORKN = 0, LWORKV = 0, LIWORKV = 0
!  .. LOCAL ARRAYS ..
   REAL(WP), TARGET :: LLZ(1,1)
   INTEGER, POINTER :: IWORK(:)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(D); LD = MAX(1,N)
   IF( PRESENT(Z) )THEN; S1Z = SIZE(Z,1); S2Z = SIZE(Z,2); LJOBZ = 'V'
   ELSE; S1Z = 1; S2Z = 1; LJOBZ = 'N'; END IF
!  .. TEST THE ARGUMENTS
    IF( N < 0 ) THEN; LINFO = -1
   ELSE IF( SIZE( E ) /= N .AND. N > 0 )THEN; LINFO = -2
   ELSE IF( PRESENT(Z) .AND. ( S1Z /= LD .OR. S2Z /= N ) )THEN; LINFO = -3
   ELSE IF( N > 0 )THEN
!  .. DETERMINE THE WORKSPACE
      IF( LSAME(LJOBZ,'N') )THEN
         LWORK = MAX( 1, LWORKN ); LIWORK = MAX( 1, LIWORKN )
      ELSE
         LWORK = MAX( 1+ 4*N + N**2, LWORKV ) 
         LIWORK = MAX( 3+5*N, LIWORKV )
      END IF
      ALLOCATE(WORK(LWORK), IWORK(LIWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN
         DEALLOCATE( WORK, IWORK, STAT=ISTAT1 )
         IF( LSAME(LJOBZ,'N') )THEN; LWORK = 1; LIWORK = 1
         ELSE
	   LWORK = 1+ 4*N + N**2
	   LIWORK = 3+5*N; END IF
         ALLOCATE(WORK(LWORK), IWORK(LIWORK), STAT=ISTAT)
         IF( ISTAT == 0 ) CALL ERINFO( -200, SRNAME, LINFO )
      END IF
      IF( ISTAT == 0 ) THEN
         IF( PRESENT(Z) )THEN
            CALL STEVD_F77( LJOBZ, N, D, E, Z, S1Z, WORK, LWORK, &
                         IWORK, LIWORK, LINFO )
	 ELSE
	    CALL STEVD_F77( LJOBZ, N, D, E, LLZ, S1Z, WORK, LWORK, &
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
END SUBROUTINE DSTEVD_F95
