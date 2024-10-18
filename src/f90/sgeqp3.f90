!> \brief \b SGEQP3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEQP3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqp3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqp3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqp3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEQP3 computes a QR factorization with column pivoting of a
!> matrix A:  A*P = Q*R  using Level 3 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the upper triangle of the array contains the
!>          min(M,N)-by-N upper trapezoidal matrix R; the elements below
!>          the diagonal, together with the array TAU, represent the
!>          orthogonal matrix Q as a product of min(M,N) elementary
!>          reflectors.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
!>          to the front of A*P (a leading column); if JPVT(J)=0,
!>          the J-th column of A is a free column.
!>          On exit, if JPVT(J)=K, then the J-th column of A*P was the
!>          the K-th column of A.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= 3*N+1.
!>          For optimal performance LWORK >= 2*N+( N+1 )*NB, where NB
!>          is the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit.
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real/complex vector
!>  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
!>  A(i+1:m,i), and tau in TAU(i).
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!>
!  =====================================================================
      SUBROUTINE SGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            INB, INBMIN, IXOVER
      PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB,&
     &                   NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEQRF, SLAQP2, SLAQPS, SORMQR, SSWAP, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      REAL               SNRM2
      EXTERNAL           ILAENV, SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
!     Test input arguments
!  ====================
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
!
      IF( INFO.EQ.0 ) THEN
         MINMN = MIN( M, N )
         IF( MINMN.EQ.0 ) THEN
            IWS = 1
            LWKOPT = 1
         ELSE
            IWS = 3*N + 1
            NB = ILAENV( INB, 'SGEQRF', ' ', M, N, -1, -1 )
            LWKOPT = 2*N + ( N + 1 )*NB
         END IF
         WORK( 1 ) = LWKOPT
!
         IF( ( LWORK.LT.IWS ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQP3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Move initial columns up front.
!
      NFXD = 1
      DO 10 J = 1, N
         IF( JPVT( J ).NE.0 ) THEN
            IF( J.NE.NFXD ) THEN
               CALL SSWAP( M, A( 1, J ), 1, A( 1, NFXD ), 1 )
               JPVT( J ) = JPVT( NFXD )
               JPVT( NFXD ) = J
            ELSE
               JPVT( J ) = J
            END IF
            NFXD = NFXD + 1
         ELSE
            JPVT( J ) = J
         END IF
   10 CONTINUE
      NFXD = NFXD - 1
!
!     Factorize fixed columns
!  =======================
!
!     Compute the QR factorization of fixed columns and update
!     remaining columns.
!
      IF( NFXD.GT.0 ) THEN
         NA = MIN( M, NFXD )
!CC      CALL SGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         CALL SGEQRF( M, NA, A, LDA, TAU, WORK, LWORK, INFO )
         IWS = MAX( IWS, INT( WORK( 1 ) ) )
         IF( NA.LT.N ) THEN
!CC         CALL SORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
!CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            CALL SORMQR( 'Left', 'Transpose', M, N-NA, NA, A, LDA, TAU,&
     &                   A( 1, NA+1 ), LDA, WORK, LWORK, INFO )
            IWS = MAX( IWS, INT( WORK( 1 ) ) )
         END IF
      END IF
!
!     Factorize free columns
!  ======================
!
      IF( NFXD.LT.MINMN ) THEN
!
         SM = M - NFXD
         SN = N - NFXD
         SMINMN = MINMN - NFXD
!
!        Determine the block size.
!
         NB = ILAENV( INB, 'SGEQRF', ' ', SM, SN, -1, -1 )
         NBMIN = 2
         NX = 0
!
         IF( ( NB.GT.1 ) .AND. ( NB.LT.SMINMN ) ) THEN
!
!           Determine when to cross over from blocked to unblocked code.
!
            NX = MAX( 0, ILAENV( IXOVER, 'SGEQRF', ' ', SM, SN, -1,&
     &           -1 ) )
!
!
            IF( NX.LT.SMINMN ) THEN
!
!              Determine if workspace is large enough for blocked code.
!
               MINWS = 2*SN + ( SN+1 )*NB
               IWS = MAX( IWS, MINWS )
               IF( LWORK.LT.MINWS ) THEN
!
!                 Not enough workspace to use optimal NB: Reduce NB and
!                 determine the minimum value of NB.
!
                  NB = ( LWORK-2*SN ) / ( SN+1 )
                  NBMIN = MAX( 2, ILAENV( INBMIN, 'SGEQRF', ' ', SM, SN,&
     &                    -1, -1 ) )
!
!
               END IF
            END IF
         END IF
!
!        Initialize partial column norms. The first N elements of work
!        store the exact column norms.
!
         DO 20 J = NFXD + 1, N
            WORK( J ) = SNRM2( SM, A( NFXD+1, J ), 1 )
            WORK( N+J ) = WORK( J )
   20    CONTINUE
!
         IF( ( NB.GE.NBMIN ) .AND. ( NB.LT.SMINMN ) .AND.&
     &       ( NX.LT.SMINMN ) ) THEN
!
!           Use blocked code initially.
!
            J = NFXD + 1
!
!           Compute factorization: while loop.
!
!
            TOPBMN = MINMN - NX
   30       CONTINUE
            IF( J.LE.TOPBMN ) THEN
               JB = MIN( NB, TOPBMN-J+1 )
!
!              Factorize JB columns among columns J:N.
!
               CALL SLAQPS( M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA,&
     &                      JPVT( J ), TAU( J ), WORK( J ), WORK( N+J ),&
     &                      WORK( 2*N+1 ), WORK( 2*N+JB+1 ), N-J+1 )
!
               J = J + FJB
               GO TO 30
            END IF
         ELSE
            J = NFXD + 1
         END IF
!
!        Use unblocked code to factor the last or only block.
!
!
         IF( J.LE.MINMN )&
     &      CALL SLAQP2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ),&
     &                   TAU( J ), WORK( J ), WORK( N+J ),&
     &                   WORK( 2*N+1 ) )
!
      END IF
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of SGEQP3
!
      END
!> \brief \b SGEQR2 computes the QR factorization of a general rectangular matrix using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEQR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEQR2 computes a QR factorization of a real m-by-n matrix A:
!>
!>    A = Q * ( R ),
!>            ( 0 )
!>
!> where:
!>
!>    Q is a m-by-m orthogonal matrix;
!>    R is an upper-triangular n-by-n matrix;
!>    0 is a (m-n)-by-n zero matrix, if m > n.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(m,n) by n upper trapezoidal matrix R (R is
!>          upper triangular if m >= n); the elements below the diagonal,
!>          with the array TAU, represent the orthogonal matrix Q as a
!>          product of elementary reflectors (see Further Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2019
!
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!>  and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, K
      REAL               AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, SLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQR2', -INFO )
         RETURN
      END IF
!
      K = MIN( M, N )
!
      DO 10 I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
         CALL SLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,&
     &                TAU( I ) )
         IF( I.LT.N ) THEN
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
            AII = A( I, I )
            A( I, I ) = ONE
            CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),&
     &                  A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
!
!     End of SGEQR2
!
      END
!> \brief \b SGEQRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEQRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEQRF computes a QR factorization of a real M-by-N matrix A:
!>
!>    A = Q * ( R ),
!>            ( 0 )
!>
!> where:
!>
!>    Q is a M-by-M orthogonal matrix;
!>    R is an upper-triangular N-by-N matrix;
!>    0 is a (M-N)-by-N zero matrix, if M > N.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!>          upper triangular if m >= n); the elements below the diagonal,
!>          with the array TAU, represent the orthogonal matrix Q as a
!>          product of min(m,n) elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,N).
!>          For optimum performance LWORK >= N*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2019
!
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!>  and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,&
     &                   NBMIN, NX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEQR2, SLARFB, SLARFT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
!
!        Determine when to cross over from blocked to unblocked code.
!
         NX = MAX( 0, ILAENV( 3, 'SGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'SGEQRF', ' ', M, N, -1,&
     &                 -1 ) )
            END IF
         END IF
      END IF
!
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
!
!        Use blocked code initially
!
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
!
!           Compute the QR factorization of the current block
!           A(i:m,i:i+ib-1)
!
            CALL SGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,&
     &                   IINFO )
            IF( I+IB.LE.N ) THEN
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               CALL SLARFT( 'Forward', 'Columnwise', M-I+1, IB,&
     &                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H**T to A(i:m,i+ib:n) from the left
!
               CALL SLARFB( 'Left', 'Transpose', 'Forward',&
     &                      'Columnwise', M-I+1, N-I-IB+1, IB,&
     &                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),&
     &                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
!
!     Use unblocked code to factor the last or only block.
!
      IF( I.LE.K )&
     &   CALL SGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,&
     &                IINFO )
!
      WORK( 1 ) = IWS
      RETURN
!
!     End of SGEQRF
!
      END
!> \brief \b SLAQP2 computes a QR factorization with column pivoting of the matrix block.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAQP2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqp2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqp2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqp2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,
!                          WORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, M, N, OFFSET
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               A( LDA, * ), TAU( * ), VN1( * ), VN2( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAQP2 computes a QR factorization with column pivoting of
!> the block A(OFFSET+1:M,1:N).
!> The block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in] OFFSET
!> \verbatim
!>          OFFSET is INTEGER
!>          The number of rows of the matrix A that must be pivoted
!>          but no factorized. OFFSET >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the upper triangle of block A(OFFSET+1:M,1:N) is
!>          the triangular factor obtained; the elements in block
!>          A(OFFSET+1:M,1:N) below the diagonal, together with the
!>          array TAU, represent the orthogonal matrix Q as a product of
!>          elementary reflectors. Block A(1:OFFSET,1:N) has been
!>          accordingly pivoted, but no factorized.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
!>          to the front of A*P (a leading column); if JPVT(i) = 0,
!>          the i-th column of A is a free column.
!>          On exit, if JPVT(i) = k, then the i-th column of A*P
!>          was the k-th column of A.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[in,out] VN1
!> \verbatim
!>          VN1 is REAL array, dimension (N)
!>          The vector with the partial column norms.
!> \endverbatim
!>
!> \param[in,out] VN2
!> \verbatim
!>          VN2 is REAL array, dimension (N)
!>          The vector with the exact column norms.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup realOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!> \n
!>  Partial column norm updating strategy modified on April 2011
!>    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
!>    University of Zagreb, Croatia.
!
!> \par References:
!  ================
!>
!> LAPACK Working Note 176
!
!> \htmlonly
!> <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a>
!> \endhtmlonly
!
!  =====================================================================
      SUBROUTINE SLAQP2( M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2,&
     &                   WORK )
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER            LDA, M, N, OFFSET
!     ..
!     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), TAU( * ), VN1( * ), VN2( * ),&
     &                   WORK( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITEMP, J, MN, OFFPI, PVT
      REAL               AII, TEMP, TEMP2, TOL3Z
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, SLARFG, SSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH, SNRM2
      EXTERNAL           ISAMAX, SLAMCH, SNRM2
!     ..
!     .. Executable Statements ..
!
      MN = MIN( M-OFFSET, N )
      TOL3Z = SQRT(SLAMCH('Epsilon'))
!
!     Compute factorization.
!
      DO 20 I = 1, MN
!
         OFFPI = OFFSET + I
!
!        Determine ith pivot column and swap if necessary.
!
         PVT = ( I-1 ) + ISAMAX( N-I+1, VN1( I ), 1 )
!
         IF( PVT.NE.I ) THEN
            CALL SSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( I )
            JPVT( I ) = ITEMP
            VN1( PVT ) = VN1( I )
            VN2( PVT ) = VN2( I )
         END IF
!
!        Generate elementary reflector H(i).
!
         IF( OFFPI.LT.M ) THEN
            CALL SLARFG( M-OFFPI+1, A( OFFPI, I ), A( OFFPI+1, I ), 1,&
     &                   TAU( I ) )
         ELSE
            CALL SLARFG( 1, A( M, I ), A( M, I ), 1, TAU( I ) )
         END IF
!
         IF( I.LT.N ) THEN
!
!           Apply H(i)**T to A(offset+i:m,i+1:n) from the left.
!
            AII = A( OFFPI, I )
            A( OFFPI, I ) = ONE
            CALL SLARF( 'Left', M-OFFPI+1, N-I, A( OFFPI, I ), 1,&
     &                  TAU( I ), A( OFFPI, I+1 ), LDA, WORK( 1 ) )
            A( OFFPI, I ) = AII
         END IF
!
!        Update partial column norms.
!
         DO 10 J = I + 1, N
            IF( VN1( J ).NE.ZERO ) THEN
!
!              NOTE: The following 4 lines follow from the analysis in
!              Lapack Working Note 176.
!
               TEMP = ONE - ( ABS( A( OFFPI, J ) ) / VN1( J ) )**2
               TEMP = MAX( TEMP, ZERO )
               TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
               IF( TEMP2 .LE. TOL3Z ) THEN
                  IF( OFFPI.LT.M ) THEN
                     VN1( J ) = SNRM2( M-OFFPI, A( OFFPI+1, J ), 1 )
                     VN2( J ) = VN1( J )
                  ELSE
                     VN1( J ) = ZERO
                     VN2( J ) = ZERO
                  END IF
               ELSE
                  VN1( J ) = VN1( J )*SQRT( TEMP )
               END IF
            END IF
   10    CONTINUE
!
   20 CONTINUE
!
      RETURN
!
!     End of SLAQP2
!
      END
!> \brief \b SLAQPS computes a step of QR factorization with column pivoting of a real m-by-n matrix A by using BLAS level 3.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAQPS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaqps.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaqps.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaqps.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,
!                          VN2, AUXV, F, LDF )
!
!       .. Scalar Arguments ..
!       INTEGER            KB, LDA, LDF, M, N, NB, OFFSET
!       ..
!       .. Array Arguments ..
!       INTEGER            JPVT( * )
!       REAL               A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ),
!      $                   VN1( * ), VN2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAQPS computes a step of QR factorization with column pivoting
!> of a real M-by-N matrix A by using Blas-3.  It tries to factorize
!> NB columns from A starting from the row OFFSET+1, and updates all
!> of the matrix with Blas-3 xGEMM.
!>
!> In some cases, due to catastrophic cancellations, it cannot
!> factorize NB columns.  Hence, the actual number of factorized
!> columns is returned in KB.
!>
!> Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A. N >= 0
!> \endverbatim
!>
!> \param[in] OFFSET
!> \verbatim
!>          OFFSET is INTEGER
!>          The number of rows of A that have been factorized in
!>          previous steps.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The number of columns to factorize.
!> \endverbatim
!>
!> \param[out] KB
!> \verbatim
!>          KB is INTEGER
!>          The number of columns actually factorized.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, block A(OFFSET+1:M,1:KB) is the triangular
!>          factor obtained and block A(1:OFFSET,1:N) has been
!>          accordingly pivoted, but no factorized.
!>          The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
!>          been updated.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] JPVT
!> \verbatim
!>          JPVT is INTEGER array, dimension (N)
!>          JPVT(I) = K <==> Column K of the full matrix A has been
!>          permuted into position I in AP.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL array, dimension (KB)
!>          The scalar factors of the elementary reflectors.
!> \endverbatim
!>
!> \param[in,out] VN1
!> \verbatim
!>          VN1 is REAL array, dimension (N)
!>          The vector with the partial column norms.
!> \endverbatim
!>
!> \param[in,out] VN2
!> \verbatim
!>          VN2 is REAL array, dimension (N)
!>          The vector with the exact column norms.
!> \endverbatim
!>
!> \param[in,out] AUXV
!> \verbatim
!>          AUXV is REAL array, dimension (NB)
!>          Auxiliary vector.
!> \endverbatim
!>
!> \param[in,out] F
!> \verbatim
!>          F is REAL array, dimension (LDF,NB)
!>          Matrix F**T = L*Y**T*A.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of the array F. LDF >= max(1,N).
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup realOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>    G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!>    X. Sun, Computer Science Dept., Duke University, USA
!>
!> \n
!>  Partial column norm updating strategy modified on April 2011
!>    Z. Drmac and Z. Bujanovic, Dept. of Mathematics,
!>    University of Zagreb, Croatia.
!
!> \par References:
!  ================
!>
!> LAPACK Working Note 176
!
!> \htmlonly
!> <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">[PDF]</a>
!> \endhtmlonly
!
!  =====================================================================
      SUBROUTINE SLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1,&
     &                   VN2, AUXV, F, LDF )
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER            KB, LDA, LDF, M, N, NB, OFFSET
!     ..
!     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ),&
     &                   VN1( * ), VN2( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            ITEMP, J, K, LASTRK, LSTICC, PVT, RK
      REAL               AKK, TEMP, TEMP2, TOL3Z
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMM, SGEMV, SLARFG, SSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, NINT, REAL, SQRT
!     ..
!     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH, SNRM2
      EXTERNAL           ISAMAX, SLAMCH, SNRM2
!     ..
!     .. Executable Statements ..
!
      LASTRK = MIN( M, N+OFFSET )
      LSTICC = 0
      K = 0
      TOL3Z = SQRT(SLAMCH('Epsilon'))
!
!     Beginning of while loop.
!
   10 CONTINUE
      IF( ( K.LT.NB ) .AND. ( LSTICC.EQ.0 ) ) THEN
         K = K + 1
         RK = OFFSET + K
!
!        Determine ith pivot column and swap if necessary
!
         PVT = ( K-1 ) + ISAMAX( N-K+1, VN1( K ), 1 )
         IF( PVT.NE.K ) THEN
            CALL SSWAP( M, A( 1, PVT ), 1, A( 1, K ), 1 )
            CALL SSWAP( K-1, F( PVT, 1 ), LDF, F( K, 1 ), LDF )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( K )
            JPVT( K ) = ITEMP
            VN1( PVT ) = VN1( K )
            VN2( PVT ) = VN2( K )
         END IF
!
!        Apply previous Householder reflectors to column K:
!        A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**T.
!
         IF( K.GT.1 ) THEN
            CALL SGEMV( 'No transpose', M-RK+1, K-1, -ONE, A( RK, 1 ),&
     &                  LDA, F( K, 1 ), LDF, ONE, A( RK, K ), 1 )
         END IF
!
!        Generate elementary reflector H(k).
!
         IF( RK.LT.M ) THEN
            CALL SLARFG( M-RK+1, A( RK, K ), A( RK+1, K ), 1, TAU( K ) )
         ELSE
            CALL SLARFG( 1, A( RK, K ), A( RK, K ), 1, TAU( K ) )
         END IF
!
         AKK = A( RK, K )
         A( RK, K ) = ONE
!
!        Compute Kth column of F:
!
!        Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**T*A(RK:M,K).
!
         IF( K.LT.N ) THEN
            CALL SGEMV( 'Transpose', M-RK+1, N-K, TAU( K ),&
     &                  A( RK, K+1 ), LDA, A( RK, K ), 1, ZERO,&
     &                  F( K+1, K ), 1 )
         END IF
!
!        Padding F(1:K,K) with zeros.
!
         DO 20 J = 1, K
            F( J, K ) = ZERO
   20    CONTINUE
!
!        Incremental updating of F:
!        F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**T
!                    *A(RK:M,K).
!
         IF( K.GT.1 ) THEN
            CALL SGEMV( 'Transpose', M-RK+1, K-1, -TAU( K ), A( RK, 1 ),&
     &                  LDA, A( RK, K ), 1, ZERO, AUXV( 1 ), 1 )
!
            CALL SGEMV( 'No transpose', N, K-1, ONE, F( 1, 1 ), LDF,&
     &                  AUXV( 1 ), 1, ONE, F( 1, K ), 1 )
         END IF
!
!        Update the current row of A:
!        A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**T.
!
         IF( K.LT.N ) THEN
            CALL SGEMV( 'No transpose', N-K, K, -ONE, F( K+1, 1 ), LDF,&
     &                  A( RK, 1 ), LDA, ONE, A( RK, K+1 ), LDA )
         END IF
!
!        Update partial column norms.
!
         IF( RK.LT.LASTRK ) THEN
            DO 30 J = K + 1, N
               IF( VN1( J ).NE.ZERO ) THEN
!
!                 NOTE: The following 4 lines follow from the analysis in
!                 Lapack Working Note 176.
!
                  TEMP = ABS( A( RK, J ) ) / VN1( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                  IF( TEMP2 .LE. TOL3Z ) THEN
                     VN2( J ) = REAL( LSTICC )
                     LSTICC = J
                  ELSE
                     VN1( J ) = VN1( J )*SQRT( TEMP )
                  END IF
               END IF
   30       CONTINUE
         END IF
!
         A( RK, K ) = AKK
!
!        End of while loop.
!
         GO TO 10
      END IF
      KB = K
      RK = OFFSET + KB
!
!     Apply the block reflector to the rest of the matrix:
!     A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
!                         A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**T.
!
      IF( KB.LT.MIN( N, M-OFFSET ) ) THEN
         CALL SGEMM( 'No transpose', 'Transpose', M-RK, N-KB, KB, -ONE,&
     &               A( RK+1, 1 ), LDA, F( KB+1, 1 ), LDF, ONE,&
     &               A( RK+1, KB+1 ), LDA )
      END IF
!
!     Recomputation of difficult columns.
!
   40 CONTINUE
      IF( LSTICC.GT.0 ) THEN
         ITEMP = NINT( VN2( LSTICC ) )
         VN1( LSTICC ) = SNRM2( M-RK, A( RK+1, LSTICC ), 1 )
!
!        NOTE: The computation of VN1( LSTICC ) relies on the fact that
!        SNRM2 does not fail on vectors with norm below the value of
!        SQRT(DLAMCH('S'))
!
         VN2( LSTICC ) = VN1( LSTICC )
         LSTICC = ITEMP
         GO TO 40
      END IF
!
      RETURN
!
!     End of SLAQPS
!
      END
