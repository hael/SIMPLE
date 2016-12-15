SUBROUTINE ZGGSVD_F95( A, B, ALPHA, BETA, K, L, U, V, Q, IWORK, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO
   USE F77_LAPACK, ONLY: GGSVD_F77 => LA_GGSVD
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO, K, L
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
   REAL(WP), INTENT(OUT) :: ALPHA(:), BETA(:)
   COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: U(:,:), V(:,:), Q(:,:)
   INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IWORK(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_GGSVD computes the generalized singular values and, optionally,
! the transformation matrices from the generalized singular value
! decomposition (GSVD) of a real or complex matrix pair (A,B), where A 
! is m by n and B is p by n. The GSVD of (A,B) is written 
!       A = U * SIGMA1(0, R)*Q^H , B = V * SIGMA2(0, R)*Q^H
! where U , V and Q are orthogonal (unitary) matrices of dimensions m by m,
! p by p and n by n, respectively. Let l be the rank of B and r the rank of
! the (m + p) * n matrix ( A )
!                        ( B )
! , and let k = r-l. Then SIGMA1 and SIGMA2 are m*(k + l) and p * (k + l) 
! "diagonal" matrices, respectively, and R is a (k + l) * (k + l)
! nonsingular triangular matrix. The detailed structure of SIGMA1 ,SIGMA2
! and R depends on the sign of (m - k - l) as follows:
!       The case m-k-l>=0:
! 
!                               k   l
!                      k      ( I   0 )
!        SIGMA1 =      l      ( 0   C )
!                    m-k-l    ( 0   0 )
! 
! 
!                           k   l
! 	SIGMA2 =    l   ( 0   S )
! 	          p - l ( 0   S )
! 
! 
!                          n-k-l   k     l
!          (0, R) =   k   (  0    R11   R12  )
!  	              l   (  0     0    R22  )
! 			
! where C^2 + S^2 = I . We define
! alpha(1)=alpha(2)=...=alpha(k) = 1, alpha(k+i)=c(i i), i=1,2,...,l
! beta(1) = beta(2) = ... = beta(k) = 0, beta(k+i) = s(i i), i=1,2,...,l
!  
! The case m-k-l < 0:
! 
!                                k    m-k    k+l-m
!            SIGMA1 =    k    (  I     0       0   )
!   	                m-k   (  0     C       0   )
! 		      
! 		      
! 		               k    m-k     k+l-m
!                        m-k   ( 0     S        0  )
!            SIGMA2 =   k+l-m  ( 0     0        I  )
!  	                 p-l   ( 0     0        0  )
! 
! 
!                                  n-k-l    k     m-k   k+l-m
! 		  	  k      (   0     R11    R12    R13  )
!              (0,R) =   m-k     (   0      0     R22    R23  )
! 	                k+l-m    (   0      0      0     R33  )
! 			
! where C^2 + S^2 = I . We define
!  alpha(1)=alpha(2)=...=alpha(k) = 1, alpha(k+i)=c(i i), i =1,2,...,m-k,
!  alpha(m+1) = alpha(m+2)=...= alpha(k+l) = 0
!  beta(1)=beta(2)= ... =beta(k)=0, beta(k+i)=s(i i), i=1,2,...,m-k,
!  beta(m+1) = beta(m+2) = ... = beta(k+l) = 1
! 
! In both cases the generalized singular values of the pair (A,B) are the
! ratios
!  sigma(i) = alpha(i)/beta(i), i = 1,2, ... ,k+l
!  
! The first k singular values are infinite. The finite singular values 
! are real and nonnegative.
!     LA_GGSVD computes the real (nonnegative) scalars alpha(i), beta(i),
! i=1,2,..., k+l , the matrix R, and, optionally, the transformation 
! matrices U , V and Q.
! 
! =========
! 
!      SUBROUTINE LA_GGSVD( A, B, ALPHA, BETA, K=k, L=l, &
!                       U=u, V=v, Q=q, IWORK=iwork, INFO=info )
!          <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
!          REAL(<wp>), INTENT(OUT) :: ALPHA(:), BETA(:)
!          INTEGER, INTENT(OUT), OPTIONAL :: K, L
!          <type>(<wp>), INTENT(OUT), OPTIONAL :: U(:,:), V(:,:), Q(:,:)
!          INTEGER, INTENT(IN), OPTIONAL :: IWORK(:)
!          INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!          <type> ::= REAL | COMPLEX
!          <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! A      (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(A,1) = m and size(A,2) = n.
!        On entry, the matrix A.
!        On exit, A contains the triangular matrix R, or part of R, as
!        follows:
!        If m-k-l >= 0, then R is stored in A(1:k+l,n-k-l+1:n).
!        If m-k-l < 0, then the matrix
!                  ( R11     R12    R13 )
! 		 (  0      R22    R23 )
!        is stored in A(1:m,n-k-l+1:n).
! B      (input/output) REAL or COMPLEX array, shape (:,:) with 
!        size(B,1) = p and size(B,2) = n.
!        On entry, the matrix B.
!        On exit, if m-k-l < 0, then R33 is stored in
!        B(m-k+1:l,n+m-k-l+1:n).
! ALPHA  (output) REAL array, shape (:) with size(ALPHA) = n
!        The real scalars alpha(i) , i = 1, 2,..., k+l.
! BETA   (output) REAL array, shape (:) with size(BETA) = n.
!        The real scalars beta(i) , i = 1, 2, ..., k+l.
!        Note: The generalized singular values of the pair (A,B) are
!        sigma(i) = ALPHA(i)/BETA(i), i = 1, 2, ...,  k+l.
!        If k + l < n, then ALPHA(k+l+1:n) = BETA(k+l+1:n) = 0.
! K, L   Optional (output) INTEGER.
!        The dimension parameters k and l.
! U      Optional (output) REAL or COMPLEX square array, shape (:,:) with
!        size(U,1) = m.
!        The matrix U .
! V      Optional (output) REAL or COMPLEX square array, shape (:,:) with
!        size(V,1) = p.
!        The matrix V .
! Q      Optional (output) REAL or COMPLEX square array, shape (:,:) with
!        size(Q,1) = n.
!        The matrix Q.
! IWORK  Optional (output) INTEGER array, shape(:) with size(IWORK) = n.
!        IWORK contains sorting information. More precisely, the loop
!             for i = k + 1, min(m, k + l)
!                   swap ALPHA(i) and ALPHA(IWORK(i))
!             end
!        will sort ALPHA so that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(n).
! INFO   Optional (output) INTEGER.
!        = 0: successful exit.
!        < 0: if INFO = -i, the i-th argument had an illegal value.
!        > 0: if INFO = 1, the algorithm failed to converge.
!        If INFO is not present and an error occurs, then the program is 
!        terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GGSVD'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBU, LJOBV, LJOBQ
   INTEGER :: M, N, P, LINFO, ISTAT, ISTAT1, S1U, S2U, S1V, S2V, &
              S1Q, S2Q, LK, LL
!  .. LOCAL ARRAYS ..
   COMPLEX(WP), TARGET :: LLU(1,1), LLV(1,1), LLQ(1,1)
   INTEGER, POINTER :: LIWORK(:)
   COMPLEX(WP), POINTER :: WORK(:)
   REAL(WP), POINTER :: RWORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2); P = SIZE(B,1)
   IF( PRESENT(U) )THEN; S1U = SIZE(U,1); S2U = SIZE(U,2); LJOBU = 'U'
   ELSE; S1U = 1; S2U = 1; LJOBU = 'N'; END IF
   IF( PRESENT(V) )THEN; S1V = SIZE(V,1); S2V = SIZE(V,2); LJOBV = 'V'
   ELSE; S1V = 1; S2V = 1; LJOBV = 'N'; END IF
   IF( PRESENT(Q) )THEN; S1Q = SIZE(Q,1); S2Q = SIZE(Q,2); LJOBQ = 'Q'
   ELSE; S1Q = 1; S2Q = 1; LJOBQ = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( M < 0 .OR. N < N )THEN; LINFO = -1
   ELSE IF( SIZE(B,2) /= N .OR. P < 0 ) THEN; LINFO = -2
   ELSE IF( SIZE( ALPHA ) /= N )THEN; LINFO = -3
   ELSE IF( SIZE( BETA ) /= N )THEN; LINFO = -4
   ELSE IF( PRESENT(U) .AND. ( S1U /= MAX(1,M) .OR.  S2U /= M ) )THEN; LINFO = -7
   ELSE IF( PRESENT(V) .AND. ( S1V /= MAX(1,P) .OR.  S2V /= P ) )THEN; LINFO = -8
   ELSE IF( PRESENT(Q) .AND. ( S1Q /= MAX(1,N) .OR.  S2Q /= N ) )THEN; LINFO = -9
   ELSE
    IF (PRESENT(IWORK)) THEN
      LIWORK => IWORK
    ELSE 
      ALLOCATE( LIWORK(MAX(1,N)), STAT=ISTAT )
    ENDIF
      IF( ISTAT == 0 ) THEN
         ALLOCATE( RWORK(MAX(1,2*N)), STAT=ISTAT )
      END IF
      IF( ISTAT == 0 ) THEN
         ALLOCATE( WORK( MAX(1,MAX(3*N,M,P)+N) ), STAT=ISTAT)
      END IF
      IF( ISTAT == 0 ) THEN
         IF( PRESENT(U) )THEN
           IF( PRESENT(V) )THEN
             IF( PRESENT(Q) )THEN
                  CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, U, MAX(1,S1U), &
                        V, MAX(1,S1V), Q, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		ELSE
		  CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, U, MAX(1,S1U), &
                        V, MAX(1,S1V), LLQ, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		ENDIF
	       ELSE
	        IF( PRESENT(Q) )THEN 
		   CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, U, MAX(1,S1U), &
                        LLV, MAX(1,S1V), Q, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		ELSE
		  CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, U, MAX(1,S1U), &
                        LLV, MAX(1,S1V), LLQ, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		ENDIF
	       ENDIF
	       ELSE
	         IF( PRESENT(V) )THEN
	           IF( PRESENT(Q) )THEN
		     CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, LLU, MAX(1,S1U), &
                        V, MAX(1,S1V), Q, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		   ELSE
		     CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, LLU, MAX(1,S1U), &
                        V, MAX(1,S1V), LLQ, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		   ENDIF
		  ELSE
		   IF( PRESENT(Q) )THEN
		      CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, LLU, MAX(1,S1U), &
                        LLV, MAX(1,S1V), Q, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		   ELSE
		     CALL GGSVD_F77( LJOBU, LJOBV, LJOBQ, M, N, P, LK, LL, A, MAX(1,M), &
                         B, MAX(1,P), ALPHA, BETA, LLU, MAX(1,S1U), &
                        LLV, MAX(1,S1V), LLQ, MAX(1,S1Q), WORK, RWORK, LIWORK, LINFO )
		   ENDIF
		  ENDIF
		 ENDIF
         IF( PRESENT(K) ) K = LK
         IF( PRESENT(L) ) L = LL
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, RWORK, STAT=ISTAT1)
      IF (.NOT. PRESENT(IWORK)) DEALLOCATE(LIWORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZGGSVD_F95
