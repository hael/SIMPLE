SUBROUTINE ZSYTRD_F95( A, TAU, UPLO, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: SYTRD_F77 => LA_SYTRD, ILAENV_F77 => LA_ILAENV
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT(INOUT) :: A(:,:)
   COMPLEX(WP), INTENT(OUT) :: TAU(:)
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_SYTRD'
   CHARACTER(LEN=6), PARAMETER :: BSNAME = 'ZSYTRD'
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_SYTRD / LA_HETRD reduces a real symmetric or complex Hermitian 
! matrix A to real symmetric tridiagonal form T by an orthogonal 
! or unitary similarity transformation:
! Q**H * A * Q = T.
!
! =======
!
!    SUBROUTINE LA_HETRD / LA_SYTRD|( A, TAU, UPLO, INFO )
!    .. Scalar Arguments ..
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
!       INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    .. Array Arguments ..
!       <type>(<wp>), INTENT(INOUT) :: A(:,:)
!       <type>(<wp>), INTENT(OUT) :: TAU(:)
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
!               the upper triangular part of the matrix A.
!            If UPLO = 'L', the lower triangular part of A contains
!               the lower triangular part of the matrix A.
!         On exit:
!            If UPLO = 'U', the diagonal and first superdiagonal
!               of A are overwritten by the corresponding elements of the
!               tridiagonal matrix T, and the elements above the first
!               superdiagonal, with the array TAU, represent the unitary
!               matrix Q as a product of elementary reflectors.
!            If UPLO = 'L', the diagonal and first subdiagonal of A are 
!               overwritten by the corresponding elements of the tridiagonal
!               matrix T, and the elements below the first subdiagonal, with
!               the array TAU, represent the unitary matrix Q as a product
!               of elementary reflectors.
!            See Further Details.
!
! TAU     (output) either REAL or COMPLEX array,
!         shape (:), size(TAU) == size(A,1)-1.
!         The scalar factors of the elementary reflectors.
!         See Further Details.
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
! Further Details
! ===============
!
! If UPLO = 'U', the matrix Q is represented as a product of elementary
! reflectors
!
!    Q = H(n-1) . . . H(2) H(1).
!
! Each H(i) has the form
!
!    H(i) = I - tau * v * v'
!
! where tau is a complex scalar, and v is a complex vector with
! v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
! A(1:i-1,i+1), and tau in TAU(i).
!
! If UPLO = 'L', the matrix Q is represented as a product of elementary
! reflectors
!
!    Q = H(1) H(2) . . . H(n-1).
!
! Each H(i) has the form
!
!    H(i) = I - tau * v * v'
!
! where tau is a complex scalar, and v is a complex vector with
! v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
! and tau in TAU(i).
!
! The contents of A on exit are illustrated by the following examples
! with n = 5:
!
! if UPLO = 'U':                       if UPLO = 'L':
!
!   (  d   e   v2  v3  v4 )              (  d                  )
!   (      d   e   v3  v4 )              (  e   d              )
!   (          d   e   v4 )              (  v1  e   d          )
!   (              d   e  )              (  v1  v2  e   d      )
!   (                  d  )              (  v1  v2  v3  e   d  )
!
! where d and e denote diagonal and off-diagonal elements of T, and vi
! denotes an element of the vector defining H(i).
!
! --------------------------------------
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LUPLO
   INTEGER :: LINFO, N, LD, LWORK, NB, ISTAT, ISTAT1
!  .. LOCAL ARRAYS ..
   COMPLEX(WP), POINTER :: WORK(:)
   REAL(WP), POINTER :: D(:), E(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC SIZE, MAX, PRESENT
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; N = SIZE(A,1); LD = MAX(1,N); ISTAT = 0
   IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U'; END IF
!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN; LINFO = -1
   ELSE IF( SIZE( TAU ) /= N-1 )THEN; LINFO = -2
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -3
   ELSE IF( N > 0 )THEN
!  .. DETERMINE THE WORKSPACE
      NB = ILAENV_F77( 1, BSNAME, LUPLO, N, -1, -1, -1 )
      IF( NB > 1 .AND. NB < N )THEN; LWORK = N*NB; ELSE; LWORK = 1; ENDIF
      ALLOCATE(D(N), E(N-1), WORK(LWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN; DEALLOCATE(D, E, WORK, STAT=ISTAT1)
         LWORK = 1; ALLOCATE(D(N), E(N-1), WORK(LWORK), STAT=ISTAT)
         IF( ISTAT == 0 ) CALL ERINFO( -200, SRNAME, LINFO )
      ENDIF
      IF( ISTAT == 0 )THEN
!        .. CALL LAPACK77 ROUTINE
         CALL SYTRD_F77( LUPLO, N, A, LD, D, E, TAU, WORK, LWORK, LINFO )
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(D, E, WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZSYTRD_F95
