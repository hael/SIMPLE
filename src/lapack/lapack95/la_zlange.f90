      FUNCTION ZLANGE_F95( A, NORM, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!     .. "Use Statements" ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: LANGE_F77 => LA_LANGE
!     .. "Implicit Statement" ..
      IMPLICIT NONE
!     .. "Function Type" ..
      REAL(WP) :: ZLANGE_F95
!     .. "Scalar Arguments" ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     .. "Array Arguments" ..
      COMPLEX(WP), INTENT(IN) :: A(:,:)
!-----------------------------------------------------------------
!
!  Purpose
!  =======
!
!  LA_LANGE  returns the value of the one norm,  or the Frobenius norm,
!  or the  infinity norm,  or the  element of  largest absolute value
!  of a complex matrix A.
!
!  Description
!  ===========
!
!  LA_LANGE returns the value
!
!     LA_LANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  =========
!
!  FUNCTION LA_ANGE( A, NORM, INFO )
!     REAL(<wp>) :: LA_ANGE
!     <type>(<wp>), INTENT(IN) :: <a>
!     CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: NORM
!     INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     <type> ::= REAL | COMPLEX
!     <wp>   ::= KIND(1.0) | KIND(1.0D0)
!     <a>    ::= A(:,:) | A(:)
!
!  Arguments
!  =========
!
!  A       (input) COMPLEX array, shape either (:,:) or (:).
!          If shape is (:,:) then SIZE(A,1) == m, SIZE(A,2) == n.
!          If shape is (:) then SIZE(A) == n.
!          If either m or n == 0 LA_LANGE is set to zero.
!          The m by n matrix A.
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in LA_LANGE as described
!          above.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!      If INFO is not present and an error occurs, then the program is
!         terminated with an error message.
!
!  ---------------------------------------------------------------------
!     .. "Parameters" ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_LANGE'
!     .. "Local Scalars" ..
      CHARACTER(LEN=1) :: LNORM
      INTEGER :: ISTAT, ISTAT1, LDA, LINFO, M, N
      REAL(WP), TARGET ::LLWORK(1)
!     .. "Local Pointers" ..
      REAL(WP), POINTER :: WORK(:)
!     .. "Intrinsic Functions" ..
      INTRINSIC SIZE, PRESENT, MAX
!     .. "Executable Statements" ..
      LINFO = 0; M = SIZE(A,1); N = SIZE(A,2); LDA = MAX(1,M)
      ISTAT = 0
      IF( PRESENT(NORM) )THEN; LNORM = NORM; ELSE; LNORM = '1'; ENDIF
!     .. "Testing The Arguments" ..
      IF( M < 0 .OR. N < 0 )THEN; LINFO = -1
      ELSE IF( .NOT. ( LSAME(LNORM,'M') .OR. LSAME(LNORM,'1') .OR. &
         LSAME(LNORM,'I') .OR. LSAME(LNORM,'F') .OR. &
         LSAME(LNORM,'E') ) )THEN; LINFO = -2
      ELSE
         IF( LSAME(LNORM,'I') )THEN; ALLOCATE( WORK( M), STAT=ISTAT )
         ELSE; WORK => LLWORK; ENDIF
         IF( ISTAT == 0 ) &
           ZLANGE_F95 = LANGE_F77( LNORM, M, N, A, LDA, WORK )
         IF( LSAME(LNORM,'I') )DEALLOCATE( WORK, STAT=ISTAT1 )
      ENDIF
      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
      END FUNCTION ZLANGE_F95
