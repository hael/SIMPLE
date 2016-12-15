       SUBROUTINE ZGGES_F95( A, B, ALPHA, BETA, VSL, VSR, SELECT, SDIM, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP =>  DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: GGES_F77 => LA_GGES
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
      INTEGER, INTENT(OUT), OPTIONAL :: SDIM
!  .. ARRAY ARGUMENTS ..
      COMPLEX(WP), INTENT(INOUT) :: A(:,:), B(:,:)
      COMPLEX(WP), INTENT(OUT) :: ALPHA(:), BETA(:)
      COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: VSL(:,:), VSR(:,:)
!  .. FUNCTIONAL ARGUMENTS ..
      INTERFACE
      LOGICAL FUNCTION SELECT( ALPHA, BETA)
      USE LA_PRECISION, ONLY: WP => DP
      COMPLEX(WP), INTENT(IN) :: ALPHA, BETA
      END FUNCTION SELECT
      END INTERFACE
      OPTIONAL :: SELECT
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_GGES computes for a pair of n by n real or complex matrices
! (A, B) the (generalized) real or complex Schur form, the generalized
! eigenvalues in the form of scalar pairs (alpha,beta), and, optionally,
! the left and/or right Schur vectors.
!      If A and B are real then the real-Schur form is computed, 
! otherwise the complex-Schur form is computed. The real-Schur form is a
! pair of real matrices (S,T) such that 1) S has block upper triangular 
! form, with 1 by 1 and 2 by 2 blocks along the main diagonal, 2) T has
! upper triangular form with nonnegative elements on the main diagonal,
! and 3) S = Q^T*A*Z and T = Q^T*B*Z, where Q and Z are orthogonal 
! matrices. The 2 by 2 blocks of S are "standardized" by making the 
! corresponding elements of T have the form
!                        [ a  0 ]
! 		         [ 0  b ]
! The complex-Schur form is a pair of matrices (S,T) such that 1) S has
! upper triangular form, 2) T has upper triangular form with nonnegative 
! elements on the main diagonal, and 3) S = Q^H*A*Z and T = Q^H*B*Z,
! where Q and Z are unitary matrices.
!       In both cases the columns of Q and Z are called, respectively, 
! the left and right (generalized) Schur vectors.
! A generalized eigenvalue of the pair (A,B) is, roughly speaking, a 
! scalar of the form lambda = alpha/beta such that the matrix 
! A -lambda*B is singular. It is usually represented as the pair
! (alpha, beta), as there is a reasonable interpretation of the case 
! beta = 0 (even if alpha = 0).
! 
! =========
! 
!           SUBROUTINE LA_GGES( A, B, <alpha>, BETA, VSL=vsl, &
!                      VSR=vsr, SELECT=select, SDIM=sdim, INFO=info )
!              <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
!              <type>(<wp>), INTENT(OUT) :: <alpha(:)>, BETA(:)
!              <type>(<wp>), INTENT(OUT), OPTIONAL :: VSL(:,:), VSR(:,:)
!              INTERFACE
!                 LOGICAL FUNCTION SELECT(<alpha(j)> , BETA(j))
!                     <type>(<wp>), INTENT(IN) :: <alpha(j)> , BETA(j)
!                 END FUNCTION SELECT
!              END INTERFACE
!              OPTIONAL :: SELECT
!              INTEGER, INTENT(OUT), OPTIONAL :: SDIM
!              INTEGER, INTENT(OUT), OPTIONAL :: INFO
!           where
!              <type>     ::= REAL | COMPLEX
!              <wp>       ::= KIND(1.0) | KIND(1.0D0)
!              <alpha>    ::= ALPHAR, ALPHAI | ALPHA
!              <alpha(:)> ::= ALPHAR(:), ALPHAI(:) | ALPHA(:)
!              <alpha(j)> ::= ALPHAR(j) , ALPHAI(j) | ALPHA(j)
! 
! Arguments
! =========
! 
! A        (input/output) REAL or COMPLEX square array, shape (:,:).
!          On entry, the matrix A.
!          On exit, the matrix S.
! B        (input/output) REAL or COMPLEX square array, shape (:,:) with
!          size(B,1) = size(A,1).
!          On entry, the matrix B.
!          On exit, the matrix T .
! <alpha>  (output) REAL or COMPLEX array, shape (:) with size(alpha) =
!          size(A,1).
! 	   The values of alpha.
!          <alpha(:)> ::= ALPHAR(:), ALPHAI(:) | ALPHA(:),
!          where
!          ALPHAR(:), ALPHAI(:) are of REAL type (for the real and 
! 	   imaginary parts) and ALPHA(:) is of COMPLEX type.
! BETA     (output) REAL or COMPLEX array, shape (:) with size(BETA) =
!          size(A,1).
!          The values of beta.
!          Note: The generalized eigenvalues of the pair (A,B) are the 
! 	   scalars lambda(j) = alpha(j)/beta(j). These quotients may 
! 	   easily over- or underflow, and beta(j) may even be zero. Thus,
! 	   the user should avoid computing them naively.
!          Note: If A and B are real then complex eigenvalues occur in 
! 	   complex conjugate pairs. Each pair is stored consecutively. 
! 	   Thus a complex conjugate pair is given by
!            lambda(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
!            lambda(j+1) = (ALPHAR(j+1) + i*ALPHAI(j+1))/BETA(j+1)
!          where
!          ALPHAI(j)/BETA(j) = -( ALPHAI(j+1)/BETA(j+1))
! VSL      Optional (output) REAL or COMPLEX square array, shape (:,:)
!          with size(VSL,1) = size(A,1).
!          The left Schur vectors.
! VSR      Optional (output) REAL or COMPLEX square array, shape (:,:)
!          with size(VSR,1) = size(A,1).
!          The right Schur vectors.
! SELECT   Optional (input) LOGICAL FUNCTION
!          LOGICAL FUNCTION SELECT( <alpha(j)> , BETA(j))
!             <type>(<wp>), INTENT(IN) :: <alpha(j)> , BETA(j)
!          where
!             <type> ::= REAL | COMPLEX
!             <wp> ::= KIND(1.0) | KIND(1.0D0)
!             <alpha(j)> ::= ALPHAR(j) , ALPHAI(j) | ALPHA(j)
!          1. SELECT must be declared as EXTERNAL or as an explicit 
! 	   interface in the calling (sub)program.
!          2. SELECT is called by LA_GGES for every computed eigenvalue
! 	   (<alpha(j)> , BETA(j)) (but only once for a complex conjugate
! 	   pair when A and B are real). It is used to select the
!          eigenvalues that will be ordered to the top left of the Schur 
! 	   form. The eigenvalue (<alpha(j)>, BETA(j)) is selected if 
! 	   SELECT(<alpha(j)>, BETA(j)) has the value .TRUE.
!          3. A selected complex eigenvalue may no longer satisfy 
! 	   SELECT(<alpha(j)>, BETA(j)) = .TRUE. after ordering, since 
! 	   ordering may change the value of complex eigenvalues
! 	   (especially if the eigenvalue is ill-conditioned); in this case
! 	   INFO is set to size(A,1) + 2 (see INFO below).
!          Note: Select must be present if SDIM is desired.
! SDIM     Optional (output) INTEGER.
!          The number of eigenvalues (after sorting) for which SELECT = 
! 	   .TRUE. (If A and B are real, then complex conjugate pairs for 
! 	   which SELECT = .TRUE. for either eigenvalue count as 2).
! INFO     Optional (output) INTEGER.
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, and i is
!            <= n: the QZ iteration failed. The matrix pair (A,B) has not
! 	      been reduced to Schur form, but (<alpha(j)>,BETA(j)) 
! 	      should be correct for j = INFO + 1,..., n.
!            = n+1: another part of the algorithm failed.
!            = n+2: after reordering, roundoff changed values of some 
! 	      complex eigenvalues so that leading eigenvalues in the 
! 	      Schur form no longer satisfy SELECT = .TRUE. This can be
! 	      caused by ordinary roundoff or underflow due to scaling.
!            = n+3: the reordering failed.
!          If INFO is not present and an error occurs, then the program 
!          is terminated with an error message.
!-----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GGES'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBVSL, LJOBVSR, LSORT
      INTEGER, SAVE :: LWORK = 0
      INTEGER :: N, LINFO, LDA, LDB, ISTAT, S1VSL, S2VSL, S1VSR, S2VSR, &
     &  SALPHA, SBETA
      INTEGER :: LSDIM
!  .. LOCAL ARRAYS ..
      COMPLEX(WP), TARGET :: LLVSL(1,1), LLVSR(1,1), WORKMIN(1)
      COMPLEX(WP), POINTER :: WORK(:)
        REAL(WP), POINTER :: RWORK(:)
      LOGICAL, POINTER :: BWORK(:)
      LOGICAL, TARGET :: LLBWORK(1)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LDA = MAX(1,N); LDB = MAX(1,SIZE(B,1))
   SALPHA = SIZE(ALPHA)
   IF (PRESENT (SELECT)) THEN
      LSORT = 'S'; ELSE ; LSORT = 'N'; ENDIF 
   SBETA = SIZE(BETA)
   IF( PRESENT(VSL) )THEN; S1VSL = SIZE(VSL,1); S2VSL = SIZE(VSL,2); LJOBVSL = 'V'
   ELSE; S1VSL = 1; S2VSL = 1; LJOBVSL = 'N'; END IF
   IF( PRESENT(VSR) )THEN; S1VSR = SIZE(VSR,1); S2VSR = SIZE(VSR,2); LJOBVSR = 'V'
   ELSE; S1VSR = 1; S2VSR = 1; LJOBVSR = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
      ELSE IF( SIZE(B,1) /= N .OR. SIZE(B,2) /= N )THEN; LINFO = -2
       ELSE IF( SALPHA /= N )THEN; LINFO = -3
       ELSE IF( SBETA /= N )THEN; LINFO = -4
       ELSE IF( PRESENT(VSL) .AND. ( S1VSL /= N .OR. S2VSL /= N) )THEN; LINFO = -5
       ELSE IF( PRESENT(VSR) .AND. ( S1VSR /= N .OR. S2VSR /= N ) )THEN; LINFO = -6

       ELSE IF( N >= 0 )THEN
       IF(LSAME(LSORT,'S')) THEN;  ALLOCATE(BWORK(N),STAT=ISTAT)
        IF (ISTAT /= 0) THEN ; LINFO=-100; GOTO 500; ENDIF
       ELSE; BWORK => LLBWORK; END IF
       
       ALLOCATE(RWORK(8*N), STAT = ISTAT)
       IF (ISTAT /= 0) THEN ; LINFO=-100; GOTO 100; ENDIF


! .. DETERMINE THE WORKSPACE ..
! .. QUERING THE SIZE OF WORKSPACE ..
        LWORK = -1
        IF( PRESENT(VSL) )THEN
	   IF( PRESENT(VSR) )THEN 
               CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                  LSDIM, ALPHA, BETA, VSL, MAX(1,S1VSL), VSR, MAX(1,S1VSR), &
&                  WORKMIN, LWORK, RWORK, BWORK,  LINFO )
           ELSE
	       CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                  LSDIM, ALPHA, BETA, VSL, MAX(1,S1VSL), LLVSR, MAX(1,S1VSR), &
&                  WORKMIN, LWORK, RWORK, BWORK,  LINFO )
           ENDIF
	 ELSE
           IF( PRESENT(VSR) )THEN 
               CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                  LSDIM, ALPHA, BETA, LLVSL, MAX(1,S1VSL), VSR, MAX(1,S1VSR), &
&                  WORKMIN, LWORK, RWORK, BWORK,  LINFO )
           ELSE
	       CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                  LSDIM, ALPHA, BETA, LLVSL, MAX(1,S1VSL), LLVSR, MAX(1,S1VSR), &
&                  WORKMIN, LWORK, RWORK, BWORK,  LINFO )
           ENDIF
	 ENDIF
        LWORK = WORKMIN(1)
        
        ALLOCATE(WORK(LWORK), STAT=ISTAT)
        IF (ISTAT /= 0) THEN ; LINFO=-100; GOTO 400; ENDIF

        IF( PRESENT(VSL) )THEN
	  IF( PRESENT(VSR) )THEN  
              CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                 LSDIM, ALPHA, BETA, VSL, MAX(1,S1VSL), VSR, MAX(1,S1VSR), &
&                 WORK, LWORK, RWORK, BWORK,  LINFO )
          ELSE
	      CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                 LSDIM, ALPHA, BETA, VSL, MAX(1,S1VSL), LLVSR, MAX(1,S1VSR), &
&                 WORK, LWORK, RWORK, BWORK,  LINFO )
          ENDIF
	ELSE
	  IF( PRESENT(VSR) )THEN  
              CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                 LSDIM, ALPHA, BETA, LLVSL, MAX(1,S1VSL), VSR, MAX(1,S1VSR), &
&                 WORK, LWORK, RWORK, BWORK,  LINFO )
          ELSE
	      CALL GGES_F77( LJOBVSL, LJOBVSR, LSORT, SELECT, N, A, LDA, B, LDB, &
&                 LSDIM, ALPHA, BETA, LLVSL, MAX(1,S1VSL), LLVSR, MAX(1,S1VSR), &
&                 WORK, LWORK, RWORK, BWORK,  LINFO )
          ENDIF
	ENDIF  
          IF (PRESENT(SDIM)) SDIM = LSDIM
          
          DEALLOCATE(WORK)
400       DEALLOCATE(RWORK) 
100       IF(LSAME(LSORT,'S')) DEALLOCATE(BWORK)
        ENDIF
500     CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZGGES_F95
