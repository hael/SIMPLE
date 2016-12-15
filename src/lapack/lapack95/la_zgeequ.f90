SUBROUTINE ZGEEQU_F95( A, R, C, ROWCND, COLCND, AMAX, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO
   USE F77_LAPACK, ONLY: GEEQU_F77 => LA_GEEQU
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
   REAL(WP), INTENT( OUT ), OPTIONAL :: AMAX, COLCND, ROWCND
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT( IN ) :: A( :, : )
   REAL(WP), INTENT( OUT ) :: C( : ), R( : )
!---------------------------------------------------------------------
!
! Purpose
! =======
!
! LA_GEEQU computes row and column scalings intended to equilibrate a
! rectangle matrix A and reduce its condition number.  R returns the
! row scale factors and C the column scale factors, chosen to try to
! make the largest entry in each row and column of the matrix B with
! elements B(i,j) = R(i) A(i,j) C(j) have absolute value 1.
!
! R(i) and C(j) are restricted to be between SMLNUM = smallest safe
! number and BIGNUM = largest safe number. Use of these scaling
! factors is not guaranteed to reduce the condition number of A but
! works well in practice.
!
! Arguments
! =========
!
! SUBROUTINE LA_GEEQU ( A, R, C, ROWCND, COLCND, AMAX, INFO )
!    <type>(<wp>), INTENT(IN) :: A(:,:)
!    REAL(<wp>), INTENT( OUT ) :: R(:), C(:)
!    REAL(<wp>), INTENT( OUT ), OPTIONAL :: ROWCND, COLCND, AMAX
!    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!    <type> ::= REAL | COMPLEX
!    <wp>   ::= KIND(1.0) | KIND(1.0D0)
!
! ====================
!
! A       (input) either REAL or COMPLEX array, shape (:,:).
!         The matrix A, whose equilibration factors are to be computed.
!
! R       (output) REAL array, shape (:), size(R) == size(A,1).
!         If INFO = 0 or INFO > size(A,1), R contains the row
!         scale factors for A.
!
! C       (output) REAL array, shape (:), size(C) == size(A,2).
!         If INFO = 0, C contains the column scale factors for A.
!
! ROWCND  Optional (output) REAL.
!         If INFO = 0 or INFO > size(A,1), ROWCND contains the ratio
!         of the smallest R(i) to the largest R(i).  If ROWCND >= 0.1
!         and AMAX is neither too large nor too small, it is not worth
!         scaling by R.
!
! COLCND  Optional (output) REAL.
!         If INFO = 0, COLCND contains the ratio of the smallest
!         C(i) to the largest C(i).  If COLCND >= 0.1, it is not
!         worth scaling by C.
!
! AMAX    Optional (output) REAL.
!         Absolute value of largest matrix element.  If AMAX is very
!         close to overflow or very close to underflow, the matrix
!         should be scaled.
!
! INFO    Optional (output) INTEGER
!         If INFO is present
!            = 0:  successful exit
!            < 0:  if INFO = -k, the k-th argument had an illegal value
!            > 0:  if INFO = k,  and k is
!                  <= M:  the k-th row of A is exactly zero
!                  >  M:  the (k-M)-th column of A is exactly zero
!                         where M = size(A,1)
!         If INFO is not present and an error occurs, then the program is
!            terminated with an error message.
!
!-------------------------------------------------------------------------
!  .. PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GEEQU'
!  .. LOCAL SCALARS ..
   INTEGER :: LINFO, M, N
   REAL(WP) :: LAMAX, LCOLCND, LROWCND
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC SIZE, MAX
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; M = SIZE(A, 1); N = SIZE(A, 2)
!  .. TEST THE ARGUMENTS
   IF ( SIZE(R) /= M ) THEN; LINFO = -2
   ELSE IF ( SIZE(C) /= N ) THEN; LINFO = -3
   ELSE
!     .. CALL LAPACK77 ROUTINE
      CALL GEEQU_F77( M, N, A, MAX(1,M), R, C, LROWCND, LCOLCND, LAMAX, LINFO )
      IF( PRESENT( ROWCND ) ) ROWCND = LROWCND
      IF( PRESENT( COLCND ) ) COLCND = LCOLCND
      IF( PRESENT( AMAX ) ) AMAX = LAMAX
   END IF
   CALL ERINFO( LINFO, SRNAME, INFO )
END SUBROUTINE ZGEEQU_F95
