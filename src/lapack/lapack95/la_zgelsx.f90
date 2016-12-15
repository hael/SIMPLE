 SUBROUTINE ZGELSX_F95( A, B, RANK, JPVT, RCOND, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD, ONLY: ERINFO
    USE F77_LAPACK, ONLY: GELSX_F77 => LA_GELSX
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL :: RANK
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    REAL(WP), INTENT(IN), OPTIONAL :: RCOND
!   .. ARRAY ARGUMENTS ..
    INTEGER, INTENT(INOUT), OPTIONAL, TARGET :: JPVT(:)
    COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_GELSX computes the minimum-norm solution to a real linear least
! squares problem:
!     minimize || A * X - B ||
! using a complete orthogonal factorization of A.  A is an m-by-n
! matrix which may be rank-deficient.
! Several right hand side vectors b and solution vectors x can be
! handled in a single call; they are stored as the columns of the
! M-by-NRHS right hand side matrix B and the n-by-nrhs solution
! matrix X.
! The routine first computes a QR factorization with column pivoting:
!     A * P = Q * [ R11 R12 ]
!                 [  0  R22 ]
! with R11 defined as the largest leading submatrix whose estimated
! condition number is less than 1/RCOND.  The order of R11, RANK,
! is the effective rank of A.
! Then, R22 is considered to be negligible, and R12 is annihilated
! by orthogonal transformations from the right, arriving at the
! complete orthogonal factorization:
!    A * P = Q * [ T11 0 ] * Z
!                [  0  0 ]
! The minimum-norm solution is then
!    X = P * Z' [ inv(T11)*Q1'*B ]
!               [        0       ]
! where Q1 consists of the first RANK columns of Q.
!
!
! Arguments
! =========
!
!  SUBROUTINE LA_GELSX( A, B, RANK, JPVT, RCOND, INFO )
!    <type>(<wp>), INTENT( INOUT ) :: A( :, : ), <rhs>
!    INTEGER, INTENT(IN), OPTIONAL :: RANK
!    INTEGER, INTENT(OUT), OPTIONAL :: JPVT(:)
!    REAL(<wp>), INTENT(IN), OPTIONAL :: RCOND
!    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    <type> ::= REAL | COMPLEX
!    <wp>   ::= KIND(1.0) | KIND(1.0D0)
!    <rhs>  ::= B(:,:) | B(:)
!
! =====================
!
! A    (input/output) Deither REAL or COMPLEX array, shape (:,:),
!      SIZE(A,1) == m, SIZE(A,2) == n.
!      On entry, the m-by-n matrix A.
!      On exit, A has been overwritten by details of its
!      complete orthogonal factorization.
!      INFO = -1 if SIZE(A,1) < 0 or SIZE(A,2) < 0
!
! B    Optional (input/output) either REAL or COMPLEX array, shape either
!      (:,:) or (:), size(B,1) or size(B) == size(A,1). SIZE(B,2) == nrhs.
!      On entry, the m-by-nrhs right hand side matrix B.
!      On exit, the n-by-nrhs solution matrix X.
!      If m >= n and RANK = n, the residual sum-of-squares for
!      the solution in the i-th column is given by the sum of
!      squares of elements n+1:m in that column.
!      INFO = -2 if SIZE(B,1) /= max(SIZE(A,1), SIZE(A,2)) or SIZE(B,2) < 0
!                    and if shape of B is (:,:) or
!                if SIZE(B) /= max(SIZE(A,1), SIZE(A,2)) or SIZE(B,2) < 0
!                   and if shape of B is (:)
!
! RANK Optional (output) INTEGER
!      The effective rank of A, i.e., the order of the submatrix
!      R11.  This is the same as the order of the submatrix T11
!      in the complete orthogonal factorization of A.
!
! JPVT Optional (input/output) INTEGER array, shape (:), SIZE(JPVT) == n
!      On entry, if JPVT(i) .ne. 0, the i-th column of A is an
!      initial column, otherwise it is a free column.  Before
!      the QR factorization of A, all initial columns are
!      permuted to the leading positions; only the remaining
!      free columns are moved as a result of column pivoting
!      during the factorization.
!      On exit, if JPVT(i) = k, then the i-th column of A*P
!      was the k-th column of A.
!      INFO = -4 if SIZE(S) /= SIZE(A,2)
!
! RCOND Optional (input) REAL
!      RCOND is used to determine the effective rank of A, which
!      is defined as the order of the largest leading triangular
!      submatrix R11 in the QR factorization with pivoting of A,
!      whose estimated condition number < 1/RCOND.
!
! INFO    (output) INTEGER
!      = 0:  successful exit
!      < 0:  if INFO = -i, the i-th argument had an illegal value
!      If INFO is not present and an error occurs, then the program is
!         terminated with an error message.
!--------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GELSX'
!   .. LOCAL SCALARS ..
    INTEGER :: LINFO, ISTAT, ISTAT1, LWORK, N, M, MN, NRHS, LRANK, SJPVT
    REAL(WP) :: LRCOND
!   .. LOCAL POINTERS ..
    INTEGER, POINTER :: LJPVT(:)
    COMPLEX(WP), POINTER :: WORK(:)
    REAL(WP), POINTER :: RWORK(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT, MAX, MIN, EPSILON
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2); NRHS = SIZE(B,2)
    MN = MIN(M,N)
    IF( PRESENT(RCOND) )THEN; LRCOND = RCOND; ELSE
       LRCOND = 100*EPSILON(1.0_WP) ; ENDIF
    IF( PRESENT(JPVT) )THEN; SJPVT = SIZE(JPVT); ELSE; SJPVT = N; ENDIF
!   .. TEST THE ARGUMENTS
    IF( M < 0 .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= MAX(1,M,N) .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( SJPVT /= N ) THEN; LINFO = -4
    ELSE IF( LRCOND <= 0.0_WP ) THEN; LINFO = -5
    ELSE
        IF( PRESENT(JPVT) )THEN; LJPVT => JPVT
       ELSE
         ALLOCATE( LJPVT(N), STAT = ISTAT )
         IF( ISTAT /= 0 ) THEN; LINFO = -200; GOTO 100; END IF
         LJPVT = 0
       END IF
       LWORK = 2*MAX( 1, MIN(M,N) + MAX( N,2*MIN(M,N) + NRHS ) )
       ALLOCATE( WORK(LWORK), STAT = ISTAT )
       IF( ISTAT /= 0 ) THEN; LINFO = -200; GOTO 200; END IF
       ALLOCATE( RWORK(MAX(1,2*N)), STAT = ISTAT )
       IF( ISTAT /= 0 ) THEN; LINFO = -200; GOTO 300; END IF
           
       CALL GELSX_F77( M, N, NRHS, A, MAX(1,M), B, MAX(1,M,N), &
                      LJPVT, LRCOND, LRANK, WORK, RWORK, LINFO )

       IF( PRESENT(RANK) ) RANK = LRANK
!       IF( PRESENT(JPVT) ) JPVT = LJPVT
       DEALLOCATE(RWORK, STAT = ISTAT1 )
300    DEALLOCATE(WORK, STAT = ISTAT1 )
200    IF (.NOT. PRESENT(JPVT)) DEALLOCATE(LJPVT, STAT = ISTAT1 )
     END IF
100 CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
END SUBROUTINE ZGELSX_F95
