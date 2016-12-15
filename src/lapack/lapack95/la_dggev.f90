      SUBROUTINE DGGEV_F95( A, B, ALPHAR, ALPHAI, BETA, VL, VR, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO
      USE F77_LAPACK, ONLY: GGEV_F77 => LA_GGEV
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
       REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
      REAL(WP), INTENT(OUT) :: ALPHAR(:), ALPHAI(:), BETA(:)
      REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
! LA_GGEV computes for a pair of n by n real or complex matrices (A,B) 
! the generalized eigenvalues in the form of scalar pairs (alpha, beta)
! and, optionally, the left and/or right generalized eigenvectors.
!       A generalized eigenvalue of the pair (A,B) is, roughly 
! speaking, a scalar of the form  lambda=alpha/beta such that the matrix
! A-lambda*B is singular. It is usually represented as the pair
! (alpha; beta), as there is a reasonable interpretation of the case 
! beta = 0 (even if alpha = 0).
!       A right generalized eigenvector corresponding to a generalized 
! eigenvalue lambda is a vector v such that (A-lambda*B)*v=0. A left 
! generalized eigenvector is a vector u such that u^H*(A-lambda*B)=0,
! where u^H is the conjugate-transpose of u.
!       The computation is based on the (generalized) real or complex 
! Schur form of (A,B). (See LA_GGES for details of this form.)
! 
! =========
! 
!     SUBROUTINE LA_GGEV( A, B, <alpha>, BETA, VL=vl, &
!                 VR=vr, INFO=info )
!         <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
! 	  <type>(<wp>), INTENT(OUT) :: <alpha(:)>, BETA(:)
! 	  <type>(<wp>), INTENT(OUT), OPTIONAL :: VL(:,:), VR(:,:)
! 	  INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     where
! 	  <type>     ::= REAL | COMPLEX
! 	  <wp>       ::= KIND(1.0) | KIND(1.0D0)
! 	  <alpha>    ::= ALPHAR, ALPHAI | ALPHA
! 	  <alpha(:)> ::= ALPHAR(:), ALPHAI(:) | ALPHA(:)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX square array, shape (:,:).
!         On entry, the matrix A.
! 	  On exit, A has been destroyed.
! B       (input/output) REAL or COMPLEX square array, shape (:,:) with 
!         size(B,1) = size(A,1).
!   	  On entry, the matrix B.
! 	  On exit, B has been destroyed.
! <alpha> (output) REAL or COMPLEX array, shape (:) with size(alpha) = 
!         size(A,1).
! 	  The values of alpha.
!         alpha(:) ::= ALPHAR(:), ALPHAI(:) | ALPHA(:),
! 	  where
! 	  ALPHAR(:), ALPHAI(:) are of REAL type (for the real and 
! 	  imaginary parts) and ALPHA(:) is of COMPLEX type.
! BETA    (output) REAL or COMPLEX array, shape (:) with size(BETA) =
!         size(A,1).
! 	  The values of beta.
! 	  Note: The generalized eigenvalues of the pair (A,B) are the 
! 	  scalars lambda(j)=alpha(j)/beta(j). These quotients may easily
! 	  over- or underflow, and beta(j) may even be zero. Thus, the 
! 	  user should avoid computing them naively.
!         Note: If A and B are real then complex eigenvalues occur in 
! 	  complex conjugate pairs. Each pair is stored consecutively. 
! 	  Thus a complex conjugate pair is given by
!             lambda(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
!    	      lambda(j+1) = (ALPHAR(j+1) + i*ALPHAI(j+1))/BETA(j+1)
!   	  where
!     	      ALPHAI(j)/BETA(j) = -(ALPHAI(j+1)/BETA(j+1))
! VL      Optional (output) REAL or COMPLEX square array, shape (:,:) 
!         with size(VL,1) = size(A,1).
! 	  The left generalized eigenvectors u(j) are stored in the 
! 	  columns of VL in the order of their eigenvalues. Each 
! 	  eigenvector is scaled so the largest component has 
! 	      |realpart| + |imag.part| = 1,
! 	  except that for eigenvalues with alpha = beta = 0, a zero 
!     	  vector is returned as the corresponding eigenvector.
! 	  Note: If A and B are real then complex eigenvectors, like 
! 	  their eigenvalues, occur in complex conjugate pairs. The real
! 	  and imaginary parts of the first eigenvector of the pair are 
! 	  stored in VL(:,j) and VL(:,j+1) . Thus a complex conjugate 
! 	  pair is given by
! 	  u(j) = VL(:,j) + i*VL(:,j+1), u(j+1) = VL(:,j) - i*VL(:,j+1)
! VR      Optional (output) REAL or COMPLEX square array, shape (:,:) 
!         with size(VR,1) = size(A,1).
! 	  The right generalized eigenvectors v(j) are stored in the 
! 	  columns of VR in the order of their eigenvalues. Each 
! 	  eigenvector is scaled so the largest component has 
! 	      |realpart| + |imag:part| = 1,
! 	  except that for eigenvalues with alpha = beta = 0, a zero 
!  	  vector is returned as the corresponding eigenvector.
! 	  Note: If A and B are real then complex eigenvectors, like 
! 	  their eigenvalues, occur in complex conjugate pairs. The real
! 	  and imaginary parts of the first eigenvector of the pair are
! 	  stored in VR(:,j) and VR(:,j+1) . Thus a complex conjugate 
!    	  pair is given by
! 	  v(j) = VR(:,j) + i*VR(:,j+1), v(j+1) = VR(:,j) - i*VR(:,j+1)
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
! 	  < 0: if INFO = -i, the i-th argument had an illegal value.
! 	  > 0: if INFO = i, and i is
! 	     <= n: The QZ iteration failed. No eigenvectors have been
! 	           calculated, but (alpha(j), BETA(j)) should be 
! 		   correct for j = INFO+1, ..., n.
! 	     = n+1: another part of the algorithm failed.
! 	     = n+2: a failure occurred during the computation of the
! 	           generalized eigenvectors.
! 	   If INFO is not present and an error occurs, then the program
! 	   is terminated with an error message.
!-----------------------------------------------------------------------
 
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GGEV'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBVL, LJOBVR
      INTEGER, SAVE :: LWORK = 0
      INTEGER :: N, LINFO, LD, ISTAT, S1VL, S2VL, S1VR, S2VR, &
     &  SALPHAR, SALPHAI, SBETA
!  .. LOCAL ARRAYS ..
      REAL(WP), TARGET :: LLVL(1,1), LLVR(1,1), WORKMIN(1)
      REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
   SALPHAR = SIZE(ALPHAR); SALPHAI = SIZE(ALPHAI)
   SBETA = SIZE(BETA)
   IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
   ELSE; S1VL = 1; S2VL = 1; LJOBVL = 'N'; END IF
   IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
   ELSE; S1VR = 1; S2VR = 1; LJOBVR = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
   ELSE IF( SIZE(B,1) /= N .OR. SIZE(B,2) /= N )THEN; LINFO = -2
   ELSE IF( SALPHAR /= N )THEN; LINFO = -3
   ELSE IF( SALPHAI /= N )THEN; LINFO = -4
   ELSE IF( SBETA /= N )THEN; LINFO = -5
   ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -6
   ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -7
   ELSE IF( N > 0 )THEN

    
! .. DETERMINE THE WORKSPACE ..
! .. QUERING THE SIZE OF WORKSPACE ..
      LWORK = -1
      IF (PRESENT (VL)) THEN
        IF (PRESENT (VR)) THEN
           CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&	       BETA, VL, S1VL, VR, S1VR, WORKMIN, LWORK, LINFO )
        ELSE
	   CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&	       BETA, VL, S1VL, LLVR, S1VR, WORKMIN, LWORK, LINFO )
        ENDIF
      ELSE
        IF (PRESENT (VR)) THEN
           CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&	       BETA, LLVL, S1VL, VR, S1VR, WORKMIN, LWORK, LINFO )
        ELSE
	   CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&	       BETA, LLVL, S1VL, LLVR, S1VR, WORKMIN, LWORK, LINFO )
        ENDIF
       ENDIF	
      LWORK = WORKMIN(1)
      ALLOCATE(WORK(LWORK), STAT=ISTAT)
      IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 100; ENDIF
      IF (PRESENT (VL)) THEN
        IF (PRESENT (VR)) THEN
           CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&                        BETA, VL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
        ELSE
	   CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&                        BETA, VL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
        ENDIF
       ELSE
        IF (PRESENT (VR)) THEN
           CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&                        BETA, LLVL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
        ELSE
	   CALL GGEV_F77( LJOBVL, LJOBVR, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&                        BETA, LLVL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
        ENDIF
       ENDIF 	
      IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)

      DEALLOCATE(WORK)
     ENDIF
100  CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DGGEV_F95
