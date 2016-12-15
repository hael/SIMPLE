SUBROUTINE DORGTR_F95( A, TAU, UPLO, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: ORGTR_F77 => LA_ORGTR, ILAENV_F77 => LA_ILAENV
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(IN) :: TAU(:)
   REAL(WP), INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_ORGTR / LA_UNGTR generates a real orthogonal / complex unitary
! matrix Q which is defined as the product of elementary reflectors,
! as returned by LA_SYTRD / LA_HETRD:
! 
! if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
! if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
! =======
!
!    SUBROUTINE LA_ORGTR / LA_UNGTR( A, TAU, UPLO, INFO )
!    .. Scalar Arguments ..
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!       INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    .. Array Arguments ..
!       <type>(<wp>), INTENT(IN) :: TAU(:)
!       <type>(<wp>), INTENT(INOUT) :: A(:,:)
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
!         On entry, the vectors which define the elementary
!            reflectors, as returned by LA_SYTRD or LA_HETRD.
!         On exit the orthogonal or unitary matrix Q.
!
! TAU     (input) either REAL or COMPLEX array,
!         shape (:), size(TAU) == size(A,1)-1.
!         TAU(i) must contain the scalar factor of the elementary
!         reflector H(i), as returned by LA_SYTRD or LA_HETRD.
!
! UPLO    Optional, (input) CHARACTER*1
!         If UPLO is present then:
!            = 'U':  Upper triangle of A is stored
!            = 'L':  Lower triangle of A is stored
!         otherwise UPLO = 'U' is assumed.
!
! INFO    Optional, (output) INTEGER
!         If INFO is present:
!            = 0: successful exit
!            < 0: if INFO = -i, the i-th argument had an illegal value
!         If INFO is not present and an error occurs, then the program
!            is terminated with an error message.
!
! ----------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_ORGTR'
   CHARACTER(LEN=5), PARAMETER :: BSNAM  = 'DORGQ'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LUPLO
   CHARACTER(LEN=6) :: BSNAME
   INTEGER :: LINFO, LWORK, NB, ISTAT, ISTAT1, N, LD
!  .. LOCAL ARRAYS ..
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; N = SIZE(A,1); LD = MAX(1,N); ISTAT = 0
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE( TAU ) /= N-1 )THEN; LINFO = -2
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
   ELSE IF( N > 0 )THEN
!     .. DETERMINE THE WORKSPACE
      IF( LSAME(LUPLO,'U') )THEN; BSNAME = BSNAM // 'L'
      ELSE; BSNAME = BSNAM // 'R'; ENDIF
      NB = ILAENV_F77( 1, BSNAME, ' ', N-1, N-1, N-1, -1 )
      IF( NB < 1 .OR. NB >= N ) NB = 1
      LWORK = MAX( 1, (N-1)*NB ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK, STAT=ISTAT)
         LWORK = MAX( 1, N-1 ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT == 0 ) CALL ERINFO( -200, SRNAME, LINFO )
      ENDIF
      IF( ISTAT == 0 )THEN
!     .. CALL LAPACK77 ROUTINE
         CALL ORGTR_F77( LUPLO, N, A, LD, TAU, WORK, LWORK, LINFO )
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO, SRNAME, INFO, ISTAT)
END SUBROUTINE DORGTR_F95
