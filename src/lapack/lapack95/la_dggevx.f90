      SUBROUTINE DGGEVX_F95( A, B, ALPHAR, ALPHAI, BETA, VL, VR, &
     &  BALANC, ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: GGEVX_F77 => LA_GGEVX
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(OUT), OPTIONAL :: INFO
      INTEGER, INTENT(OUT), OPTIONAL :: ILO,IHI
      REAL(WP), INTENT(OUT), OPTIONAL :: ABNRM, BBNRM
!  .. ARRAY ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: BALANC
      REAL(WP), INTENT(INOUT) :: A(:,:), B(:,:)
      REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: LSCALE(:), RSCALE(:), &
     &  RCONDE(:), RCONDV(:)
      REAL(WP), INTENT(OUT) :: ALPHAR(:), ALPHAI(:), BETA(:)
      REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!      LA_GGEVX computes for a pair of n-by-n real or complex matrices
! (A, B) the generalized eigenvalues in the form of scalar pairs 
! (alpha; beta) and, optionally, the left and/or right generalized 
! eigenvectors.
!      A generalized eigenvalue of the pair (A; B) is, roughly speaking, 
! a scalar of the form lambda = alpha / beta such that the matrix 
! A - lambda * B is singular. It is usually represented as the pair
! (alpha, beta), as there is a reasonable interpretation of the case 
! beta = 0 (even if alpha = 0).
!      A right generalized eigenvector corresponding to a generalized 
! eigenvalue lambda is a vector  v  such that ( A - lambda*B)* v = 0. 
! A left generalized eigenvector is a vector u such that 
! u^H*(A-lambda*B) = 0, where u^H is the conjugate-transpose of u.
!      The computation is based on the (generalized) real or complex 
! Schur form of (A, B). (See LA_GGES for details of this form.)
!      Optionally, LA_GGEVX also computes a balancing transformation 
! (to improve the conditioning of the eigenvalues and eigenvectors), 
! reciprocal condition numbers for the eigenvalues, and reciprocal 
! condition numbers for the right eigenvectors. The balancing 
! transformation consists of a permutation of rows and columns and/or a
! scaling of rows and columns.
! 
! ==========
! 
!    SUBROUTINE LA_GGEVX( A, B, <alpha>, BETA, VL=vl, &
!          VR=vr, BALANC=balanc, ILO=ilo, IHI=ihi, &
!          LSCALE=lscale, RSCALE=rscale, ABNRM=abnrm, &
!          BBNRM=bbnrm, RCONDE=rconde, RCONDV=rcondv, &
!          INFO=info )
!        <type>(<wp>), INTENT(INOUT) :: A(:,:), B(:,:)
!        <type>(<wp>), INTENT(OUT) :: <alpha(:)>, BETA(:)
!        <type>(<wp>), INTENT(OUT), OPTIONAL :: VL(:,:), VR(:,:)
!        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: BALANC
!        INTEGER, INTENT(OUT), OPTIONAL :: ILO, IHI
!        REAL(<wp>), INTENT(OUT), OPTIONAL :: LSCALE(:),
!             RSCALE(:), RCONDE(:), RCONDV(:)
!        REAL(<wp>), INTENT(OUT), OPTIONAL :: ABNRM, BBNRM
!        INTEGER, INTENT(OUT), OPTIONAL :: INFO
!     where
!        <type>     ::= REAL | COMPLEX
!        <wp>       ::= KIND(1.0) | KIND(1.0D0)
!        <alpha>    ::= ALPHAR, ALPHAI | ALPHA
!        <alpha(:)> ::= ALPHAR(:), ALPHAI(:) | ALPHA(:)
! 
! 
! Arguments
! =========
! 
! A        (input/output) REAL or COMPLEX square array, shape (:, :).
!          On entry, the matrix A.
!          On exit, A has been overwritten. If the left, the right or 
!          both generalized eigenvectors are computed, then A contains 
! 	   the first part of the real/complex Schur form of the
!          "balanced" versions of the matrix pair (A, B).
! B        (input/output) REAL or COMPLEX square array, shape (:, :) 
!          with size(B, 1) = size(A, 1).
!          On entry, the matrix B.
!          On exit, B has been overwritten. If the left, the right or 
! 	   both generalized eigenvectors are computed, then B contains 
! 	   the second part of the real/complex Schur form of the "bal-
!          anced" versions of the matrix pair (A, B).
! <alpha>  (output) REAL or COMPLEX array, shape (:) with 
!          size(<alpha>) = size(A, 1).
!          The values of alpha.
!          <alpha(:)> ::= ALPHAR(:), ALPHAI(:) |  ALPHA(:),
!          where
!          ALPHAR(:), ALPHAI(:) are of REAL type (for the real and 
! 	   imaginary parts) and ALPHA(:) is of COMPLEX type.
! BETA     (output) REAL or COMPLEX array, shape (:) with 
!          size(BETA) = size(A,1).
!          The values of beta.
!          Note: The generalized eigenvalues of the pair (A, B) are the
! 	   scalars lambda(j) = alpha(j) / beta(j) . These quotients
!          may easily over- or underflow, and beta(j) may even be zero. 
! 	   Thus, the user should avoid computing them naively.
!          Note: If A and B are real then complex eigenvalues occur in 
! 	   complex conjugate pairs. Each pair is stored consecutively. 
! 	   Thus a complex conjugate pair is given by
!                 lambda(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
!                 lambda(j+1) = (ALPHAR(j+1) + i*ALPHAI(j+1))/BETA(j+1)
!          where
!              ALPHAI(j)/BETA(j)= - (ALPHAI(j+1)/BETA(j+1))
! VL       Optional (output) REAL or COMPLEX square array, shape (:, :)
!          with  size(VL, 1) = size(A, 1).
!         The left generalized eigenvectors u(j) are stored in the
!         columns of VL in the order of their eigenvalues. Each
!         eigenvector is scaled so the largest component has
!             |realpart| + |imag.part| = 1,
!         except that for eigenvalues with alpha = beta = 0, a zero
!         vector is returned as the corresponding eigenvector.
!         Note: If A and B are real then complex eigenvectors, like
!         their eigenvalues, occur in complex conjugate pairs. The real
!         and imaginary parts of the first eigenvector of the pair are
!         stored in VL(:,j) and VL(:,j+1) . Thus a complex conjugate
!         pair is given by
!         u(j) = VL(:,j) + i*VL(:,j+1), u(j+1) = VL(:,j) - i*VL(:,j+1)
! VR      Optional (output) REAL or COMPLEX square array, shape (:,:)
!         with size(VR,1) = size(A,1).
!         The right generalized eigenvectors v(j) are stored in the 
!         columns of VR in the order of their eigenvalues. Each 
! 	  eigenvector is scaled so the largest component has 
! 	  | realpart | + | imag.part | = 1,except that for eigenvalues
! 	  with alpha = beta = 0, a zero vector is returned as the 
! 	  corresponding eigenvector.
!         Note: If A and B are real then complex eigenvectors, like 
! 	  their eigenvalues, occur in complex conjugate pairs. The real
! 	  and imaginary parts of the first eigenvector of the pair are
! 	  stored in VR(:,j) and VR(:,j+1) . Thus a complex conjugate 
! 	  pair is given by
!         v(j) = VR(:,j) + i*VR(:,j+1), v(j+1) = VR(:,j) - i*VR(:,j+1)
! BALANC  Optional (input) CHARACTER(LEN=1).
!         Specifies the balance option to be performed.
!            = 'N': do not permute or scale;
!            = 'P': permute only;
!            = 'S': scale only;
!            = 'B': both permute and scale.
!         Default value: 'N".
!         Note: Computed reciprocal condition numbers will be for the
! 	  matrices after balancing. Permuting does not change condition
! 	  numbers (in exact arithmetic), but scaling does.
! ILO,IHI Optional (output) INTEGER.
!         ILO and IHI are integer values such that on exit A(i,j) = 0
! 	  and B(i,j) = 0 if i > j and j =1,...,ILO-1 or 
! 	  i = IHI+1,...,n.
! 	  If BALANC = 'N' or 'S', then ILO = 1 and IHI = n.
! LSCALE  Optional (output) REAL array, shape (:) with size(LSCALE) = 
!         size(A, 1).
!         Details of the permutations and scaling factors applied to 
! 	  the left side of A and B. If PL(j) is the index of the row 
! 	  interchanged with row j, and DL(j) is the scaling factor
!         applied to row j, then
!                 PL(j) = LSCALE(j),  j = 1,...,ILO-1 and IHI+1,..., n
! 	  and
! 	        DL(j) = LSCALE(j),  j = ILO, ..., IHI
! RSCALE  Optional (output) REAL array, shape (:), size(RSCALE) = 
!         size(A, 1).
!         Details of the permutations and scaling factors applied to the
! 	  right side of A and B. If PR(j) is the index of the column 
! 	  interchanged with column j, and DR(j) is the scaling factor 
! 	  applied to column j, then
!                PR(j) = RSCALE(j),  j = 1,...,ILO-1 and IHI+1,..., n
!         and
!                DR(j) = RSCALE(j),  j = ILO, ..., IHI
! ABNRM   Optional (output) REAL.
!         The l1 norm of A after balancing.
! BBNRM   Optional (output) REAL.
!         The l 1 norm of B after balancing.
! RCONDE  Optional (output) REAL array, shape (:) with size(RCONDE) = 
!         size(A, 1).
!         The reciprocal condition numbers of the eigenvalues.
! RCONDV  Optional (output) REAL array, shape (:) with size(RCONDE) = 
!         size(A, 1).
!         The estimated reciprocal condition numbers of the right 
! 	  eigenvectors. If the eigenvalues cannot be reordered to 
! 	  compute RCONDV(j) then RCONDV(j) is set to 0. This can only
!         occur when the true value would be very small.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, and i is 
! 	    <= n: The QZ iteration failed. No eigenvectors have been
! 	          calculated, but (alpha(j) , BETA(j) ) should be 
! 	   correct for j = INFO + 1,..., n.
!             = n+1: another part of the algorithm failed.
!             = n+2: a failure occurred during the computation of the
! 	          generalized eigenvectors.
!         If INFO is not present and an error occurs, then the program 
! 	is terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GGEVX'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBVL, LJOBVR, LBALANC, LSENSE
      INTEGER, SAVE :: LWORK = 0
      INTEGER :: N, LINFO, LD, ISTAT, S1VL, S2VL, S1VR, S2VR, &
     &  SALPHAR, SALPHAI, SBETA, LILO, LIHI, SRCONDE, SRCONDV, SLSCALE, SRSCALE
      REAL(WP) :: LABNRM, LBBNRM
!  .. LOCAL ARRAYS ..
      REAL(WP), TARGET :: LLVL(1,1), LLVR(1,1), WORKMIN(1)
      REAL(WP), POINTER :: LRCONDE(:), LRCONDV(:), LLSCALE(:), LRSCALE(:)
      REAL(WP), POINTER :: WORK(:)
      INTEGER, POINTER :: IWORK(:)
      LOGICAL, POINTER :: BWORK(:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
   SALPHAR = SIZE(ALPHAR); SALPHAI = SIZE(ALPHAI); SBETA = SIZE(BETA)
   IF( PRESENT(BALANC) )THEN; LBALANC = BALANC; ELSE; LBALANC = 'N'; ENDIF
   IF( PRESENT(LSCALE) )THEN; SLSCALE = SIZE(LSCALE); ELSE; SLSCALE = N; ENDIF
   IF( PRESENT(RSCALE) )THEN; SRSCALE = SIZE(RSCALE); ELSE; SRSCALE = N; ENDIF
   IF( PRESENT(RCONDE) )THEN; SRCONDE = SIZE(RCONDE); ELSE; SRCONDE = N; ENDIF
   IF( PRESENT(RCONDV) )THEN; SRCONDV = SIZE(RCONDV); ELSE; SRCONDV = N; ENDIF
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
   ELSE IF( .NOT.( LSAME(LBALANC,'N') .OR. LSAME(LBALANC,'P') .OR. &
     &  LSAME(LBALANC,'S') .OR. LSAME(LBALANC,'B') ) )THEN; LINFO = -8
   ELSE IF( SLSCALE /= N )THEN; LINFO = -11
   ELSE IF( SRSCALE /= N )THEN; LINFO = -12
   ELSE IF( SRCONDE /= N )THEN; LINFO = -15
   ELSE IF( SRCONDV /= N )THEN; LINFO = -16
   ELSE IF( N > 0 )THEN
      IF( PRESENT(RCONDE).AND.PRESENT(RCONDV) )THEN; LSENSE = 'B'
      ELSE IF( PRESENT(RCONDE) )THEN; LSENSE = 'E'
      ELSE IF( PRESENT(RCONDV) )THEN; LSENSE = 'V'; ELSE; LSENSE = 'N'; ENDIF

      ALLOCATE(BWORK(N), STAT=ISTAT)
      IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 1000; ENDIF

      ALLOCATE(IWORK(N+6), STAT=ISTAT)
      IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 900; ENDIF

      IF( PRESENT(LSCALE) )THEN; LLSCALE => LSCALE
      ELSE; ALLOCATE( LLSCALE(N), STAT=ISTAT )
        IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 800; ENDIF
      END IF

      IF( PRESENT(RSCALE) )THEN; LRSCALE => RSCALE
      ELSE; ALLOCATE( LRSCALE(N), STAT=ISTAT )
        IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 700; ENDIF
      END IF

      IF( PRESENT(RCONDV) )THEN; LRCONDV => RCONDV
      ELSE; ALLOCATE( LRCONDV(N), STAT=ISTAT )
        IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 600; ENDIF
      END IF

      IF( PRESENT(RCONDE) )THEN; LRCONDE => RCONDE
      ELSE; ALLOCATE( LRCONDE(N), STAT=ISTAT )
        IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 500; ENDIF
      END IF 

! .. DETERMINE THE WORKSPACE ..
! .. QUERING THE SIZE OF WORKSPACE ..
      LWORK = -1
      IF (PRESENT (VL)) THEN
        IF (PRESENT (VR)) THEN
           CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, VL, S1VL, VR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM, &
&             LRCONDE, LRCONDV, WORKMIN, LWORK, IWORK, BWORK, LINFO )
        ELSE
	   CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, VL, S1VL, LLVR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM, &
&             LRCONDE, LRCONDV, WORKMIN, LWORK, IWORK, BWORK, LINFO )
        ENDIF
       ELSE
        IF (PRESENT (VR)) THEN
           CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, LLVL, S1VL, VR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM, &
&             LRCONDE, LRCONDV, WORKMIN, LWORK, IWORK, BWORK, LINFO )
        ELSE
	   CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, LLVL, S1VL, LLVR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM, &
&             LRCONDE, LRCONDV, WORKMIN, LWORK, IWORK, BWORK, LINFO )
        ENDIF
       ENDIF	
      LWORK = WORKMIN(1)
      
      ALLOCATE(WORK(LWORK), STAT=ISTAT)
      IF (ISTAT /= 0) THEN; LINFO=-100; GOTO 100; ENDIF
      IF (PRESENT (VL)) THEN
        IF (PRESENT (VR)) THEN
     	    CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, VL, S1VL, VR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM,&
&             LRCONDE, LRCONDV, WORK, LWORK, IWORK, BWORK, LINFO )
        ELSE
	    CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, VL, S1VL, LLVR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM,&
&             LRCONDE, LRCONDV, WORK, LWORK, IWORK, BWORK, LINFO )
        ENDIF
      ELSE
        IF (PRESENT (VR)) THEN
     	    CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, LLVL, S1VL, VR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM,&
&             LRCONDE, LRCONDV, WORK, LWORK, IWORK, BWORK, LINFO )
        ELSE
	    CALL GGEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, B, LD, ALPHAR, ALPHAI, &
&             BETA, LLVL, S1VL, LLVR, S1VR, LILO, LIHI, LLSCALE, LRSCALE, LABNRM, LBBNRM,&
&             LRCONDE, LRCONDV, WORK, LWORK, IWORK, BWORK, LINFO )
        ENDIF
      ENDIF	
      IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
      
      IF( PRESENT(ILO) ) ILO = LILO; IF( PRESENT(IHI) ) IHI = LIHI
      IF( PRESENT(ABNRM) ) ABNRM = LABNRM
      IF( PRESENT(BBNRM) ) BBNRM = LBBNRM

      DEALLOCATE(WORK)
100  IF (.NOT. PRESENT(RCONDE)) DEALLOCATE(LRCONDE) 
500  IF (.NOT. PRESENT(RCONDV)) DEALLOCATE(LRCONDV) 
600  IF (.NOT. PRESENT(RSCALE)) DEALLOCATE(LRSCALE)
700  IF (.NOT. PRESENT(LSCALE)) DEALLOCATE(LLSCALE)
800  DEALLOCATE (IWORK) 
900  DEALLOCATE (BWORK) 
     ENDIF
1000  CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DGGEVX_F95
