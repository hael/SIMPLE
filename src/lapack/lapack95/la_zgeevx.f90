SUBROUTINE ZGEEVX_F95( A, W, VL, VR, BALANC, ILO, IHI, &
                       SCALE, ABNRM, RCONDE, RCONDV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: GEEVX_F77 => LA_GEEVX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: BALANC
   INTEGER, INTENT(OUT), OPTIONAL :: INFO, ILO, IHI
   REAL(WP), INTENT(OUT), OPTIONAL :: ABNRM
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT(INOUT) :: A(:,:)
   COMPLEX(WP), INTENT(OUT) :: W(:) 
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: SCALE(:), RCONDE(:), RCONDV(:)
   COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!        LA_GEEVX computes for a real or complex square matrix A, the 
! eigenvalues and, optionally, the left and/or right eigenvectors. 
! Optionally, it also balances A and computes reciprocal condition
! numbers for the  eigenvalues and right eigenvectors.
! A right eigenvector v(j) of A satisfies
!                   A * v(j) = lambda(j) * v(j)
! where lambda(j) is its eigenvalue. A left eigenvector u(j) of A 
! satisffies
!                   u(j)^H * A = lambda(j) * u(j)^H
! where u(j)^H denotes the conjugate-transpose of u(j). The computed
! eigenvectors are normalized to have Euclidean norm equal to 1 and 
! largest component real.
!        Balancing A involves permuting its rows and columns to make 
! it more nearly upper triangular and then scaling rows and columns by
! a diagonal similarity transformation to reduce the condition numbers 
! of the eigenvalues and eigenvectors.
!        Computed reciprocal condition numbers pertain to the matrix 
! after balancing. Permuting does not change condition numbers (in 
! exact arithmetic), but scaling does.
!
! =========
! 
!    SUBROUTINE LA_GEEVX( A, <w>, VL=vl, VR=vr, BALANC=balanc, ILO=ilo, &
!                     IHI=ihi, SCALE=scale, ABNRM=abnrm, RCONDE=rconde, &
!                     RCONDV=rcondv, INFO=info )
!         <type>(<wp>), INTENT(INOUT) :: A(:,:)
!         <type>(<wp>), INTENT(OUT) :: <w(:)>
!         <type>(<wp>), INTENT(OUT), OPTIONAL :: VL(:,:), VR(:,:)
!         CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: BALANC
!         INTEGER, INTENT(OUT), OPTIONAL :: ILO, IHI
!         REAL(<wp>), INTENT(OUT), OPTIONAL :: SCALE(:), ABNRM, &
!               RCONDE(:), RCONDV(:)
!         INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!         <type> ::= REAL | COMPLEX
!         <wp>   ::= KIND(1.0) | KIND(1.0D0)
!         <w>    ::= WR, WI | W
!         <w(:)> ::= WR(:), WI(:) | W(:)
! 
! Arguments
! =========
! 
! A        (input/output) REAL or COMPLEX square array, shape (:,:).
!          On entry, the matrix A.
!          On exit, the contents of A are destroyed.
! <w>      (output) REAL or COMPLEX array, shape (:) with size(w) = 
!          size(A,1).
!          The computed eigenvalues.
!          <w(:)> ::= WR(:), WI(:) | W(:),
!          where
!          WR(:), WI(:) are of REAL type (for the real and imaginary
!          parts) and W(:) is of COMPLEX type.
!          Note: If A is real, then a complex-conjugate pair appear 
!          consecutively, with the eigenvalue having the positive 
!          imaginary part appearing first.
! VL       Optional (output) REAL or COMPLEX square array, shape (:,:)
!          with size(VL,1) = size(A,1).
!          The left eigenvectors u(j) are stored in the columns of VL in
!          the order of their eigenvalues. Each eigenvector is scaled so
!          that the Euclidean norm is 1 and the largest component is real.
!          Note: If A is real then complex eigenvectors, like their 
!          eigenvalues, occur in complex conjugate pairs. The real and 
!          imaginary parts of the first eigenvector of the pair are
!          stored in VL(:,j) and VL(:,j+1). Thus a complex conjugate pair
!          is given by 
!            u(j) = VL(:,j) + i*VL(:,j+1), u(j+1) = VL(:,j) - i*VL(:,j+1)
! VR       Optional (output) REAL or COMPLEX square array, shape (:,:) 
!          with size(VR,1) = size(A,1).
!          The right eigenvectors v(j) are stored in the columns of VR in
!          the order of their eigenvalues.
!          Each eigenvector is scaled so that the Euclidean norm is 1 and 
!          the largest component is real.
!          Note: If A is real then complex eigenvectors, like their 
!          eigenvalues, occur in complex conjugate pairs. The real and 
!          imaginary parts of the first eigenvector of the pair are stored
!          in VR(:,j) and VR(:,j+1). Thus a complex conjugate pair is 
!          given by 
! 	     v(j) = VR(:,j) + i*VR(:,j+1), v(j+1) = VR(:,j) - i*VR(:,j+1)
! BALANC   Optional (input) CHARACTER(LEN=1).
!          Indicates whether the input matrix should be permuted and/or 
!          diagonally scaled.
!             = 'N': Do not permute or scale;
!             = 'P': Permute but do not scale;
!             = 'S': Scale but do not permute;
!             = 'B': Both permute and scale.
!           Default value: 'N'.
! ILO,IHI  Optional (output) INTEGER.
!          ILO and IHI are determined when A is balanced. The balanced
!          A(i,j) = 0 if i > j and j = 1, ..., ILO-1 or 
!          i = IHI+1, ... , size(A,1).
! SCALE    Optional (output) REAL array, shape (:) with size(SCALE) = 
!          size(A,1).
!          Details of the permutations and scaling factors applied when 
!          balancing A. If P(j) is the index of the row and column 
!          interchanged with row and column j, and D(j) is the
!          scaling factor applied to row and column j, then
!          P(j) = SCALE(j), j = 1, ..., ILO-1 and j =IHI+1, ...,  n
!          D(j) = SCALE(j), j = ILO, ... , IHI.
! ABNRM    Optional (output) REAL.
!          The l1 norm of the balanced matrix (the maximum of the sum 
!          of absolute values of elements of any column).
! RCONDE   Optional (output) REAL array, shape (:) with size(RCONDE) =
!          size(A,1). RCONDE(j) is the reciprocal condition number of 
!          the j-th eigenvalue.
! RCONDV   Optional (output) REAL array, shape (:), size(RCONDV) = 
!          size(A,1). RCONDV(j) is the reciprocal condition number of 
!          the j-th right eigenvector.
! INFO     Optional (output) INTEGER.
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, the QR algorithm failed to compute all the 
!          eigenvalues and no eigenvectors or condition numbers were
!          computed; elements 1:ILO-1 and i+1:n of <w> contain
!          eigenvalues which have converged.
!          If INFO is not present and an error occurs, then the program
!          is terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GEEVX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LBALANC, LJOBVL, LJOBVR, LSENSE
   INTEGER, SAVE :: LWORK = 0
   INTEGER :: N, LINFO, LD, ISTAT, ISTAT1, S1VL, S2VL, S1VR, S2VR, NN, &
              LILO, LIHI, SSCALE, SRCONDE, SRCONDV
   REAL(WP) :: LABNRM
!  .. LOCAL ARRAYS ..
   REAL(WP), POINTER :: RWORK(:)
   COMPLEX(WP) :: WORKMIN(1)
   REAL(WP), POINTER :: LSCALE(:), LRCONDE(:), LRCONDV(:)
   COMPLEX(WP), POINTER :: LLVL(:,:), LLVR(:,:)
   COMPLEX(WP), POINTER :: WORK(:)
   
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
   IF( PRESENT(BALANC) )THEN; LBALANC = BALANC; ELSE; LBALANC = 'N'; ENDIF
   IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
   ELSE; S1VL = 1; S2VL = N; LJOBVL = 'N'; END IF
   IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
   ELSE; S1VR = 1; S2VR = N; LJOBVR = 'N'; END IF
   IF( PRESENT(SCALE) )THEN; SSCALE = SIZE(SCALE); ELSE; SSCALE = N; ENDIF
   IF( PRESENT(RCONDE) )THEN; SRCONDE = SIZE(RCONDE); ELSE; SRCONDE = N; ENDIF
   IF( PRESENT(RCONDV) )THEN; SRCONDV = SIZE(RCONDV); ELSE; SRCONDV = N; ENDIF
!  .. TEST THE ARGUMENTS
   IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
   ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -3
   ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -4
   ELSE IF( .NOT.( LSAME(LBALANC,'N') .OR. LSAME(LBALANC,'P') .OR. &
                   LSAME(LBALANC,'S') .OR. LSAME(LBALANC,'B') ) )THEN; LINFO = -5
   ELSE IF( SSCALE /= N )THEN; LINFO = -8
   ELSE IF( SRCONDE /= N )THEN; LINFO = -10
   ELSE IF( SRCONDV /= N )THEN; LINFO = -11
   ELSE IF( N > 0 )THEN
      IF( PRESENT(BALANC) )THEN; LBALANC = BALANC; ELSE; LBALANC = 'N'; ENDIF
      IF( PRESENT(RCONDE).AND.PRESENT(RCONDV) )THEN; LSENSE = 'B'
      ELSE IF( PRESENT(RCONDE) )THEN; LSENSE = 'E'
      ELSE IF( PRESENT(RCONDV) )THEN; LSENSE = 'V'; ELSE; LSENSE = 'N'; ENDIF
      IF( PRESENT(SCALE) )THEN; LSCALE => SCALE
      ELSE
        ALLOCATE( LSCALE(N), STAT=ISTAT )
        IF (ISTAT /= 0) THEN; LINFO =-100; GOTO 100; END IF
      ENDIF
      IF( PRESENT(RCONDE) )THEN; LRCONDE => RCONDE
      ELSE
        ALLOCATE( LRCONDE(N), STAT=ISTAT )
        IF (ISTAT /= 0) THEN; LINFO =-100; GOTO 200; END IF
      ENDIF
      IF( PRESENT(RCONDV) )THEN; LRCONDV => RCONDV
      ELSE
        ALLOCATE( LRCONDV(N), STAT=ISTAT )
        IF (ISTAT /= 0) THEN; LINFO =-100; GOTO 300; END IF
      END IF
      ALLOCATE( RWORK(MAX(1,2*N)), STAT=ISTAT)
      IF (ISTAT /= 0) THEN; LINFO =-100; GOTO 400; END IF
      IF( LSAME(LSENSE,'N') .OR. LSAME(LSENSE,'E') )THEN; NN = 2*N
      ELSE IF( LSAME(LSENSE,'V') .OR. LSAME(LSENSE,'B') )THEN; NN = N*(N+2)
      ELSE; NN = 2*N; ENDIF
      IF (LSAME (LSENSE,'E') .OR. LSAME(LSENSE,'B') ) THEN
         LJOBVL = 'V'; LJOBVR = 'V'; S1VR=N; S1VL=N
      ENDIF
      IF (PRESENT (VL)) THEN; LLVL => VL
      ELSE; ALLOCATE (LLVL(N,N), STAT = ISTAT); ENDIF
      IF (PRESENT (VR)) THEN; LLVR => VR
      ELSE; ALLOCATE (LLVR(N,N), STAT = ISTAT); ENDIF
      LWORK = -1
      IF (PRESENT (VL)) THEN
         IF (PRESENT(VR)) THEN
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&              VL, S1VL, VR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&              LRCONDE, LRCONDV, WORKMIN, LWORK, RWORK, LINFO )
          ELSE
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&              VL, S1VL, LLVR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&              LRCONDE, LRCONDV, WORKMIN, LWORK, RWORK, LINFO )
          ENDIF
	ELSE
         IF (PRESENT(VR)) THEN
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&              LLVL, S1VL, VR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&              LRCONDE, LRCONDV, WORKMIN, LWORK, RWORK, LINFO )
          ELSE
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&              LLVL, S1VL, LLVR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&              LRCONDE, LRCONDV, WORKMIN, LWORK, RWORK, LINFO )
          ENDIF
         ENDIF	
       LWORK = WORKMIN(1)
      ALLOCATE(WORK(LWORK), STAT=ISTAT)
      IF( ISTAT /= 0 )THEN
        DEALLOCATE(WORK, STAT=ISTAT1)
        LWORK = MAX( 1, NN); ALLOCATE(WORK(LWORK), STAT=ISTAT)
        IF (ISTAT /= 0) THEN; LINFO =-100; GOTO 500; END IF
      END IF
      IF (PRESENT (VL)) THEN
         IF (PRESENT(VR)) THEN
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&                VL, S1VL, VR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&                LRCONDE, LRCONDV, WORK, LWORK, RWORK, LINFO )
         ELSE
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&                VL, S1VL, LLVR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&                LRCONDE, LRCONDV, WORK, LWORK, RWORK, LINFO )
	 ENDIF
       ELSE
         IF (PRESENT(VR)) THEN
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&                LLVL, S1VL, VR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&                LRCONDE, LRCONDV, WORK, LWORK, RWORK, LINFO )
         ELSE
            CALL GEEVX_F77( LBALANC, LJOBVL, LJOBVR, LSENSE, N, A, LD, W, &
&                LLVL, S1VL, LLVR, S1VR, LILO, LIHI, LSCALE, LABNRM, &
&                LRCONDE, LRCONDV, WORK, LWORK, RWORK, LINFO )
	 ENDIF
        ENDIF	 
      IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)

      IF( PRESENT(ILO) ) ILO = LILO;
      IF( PRESENT(IHI) ) IHI = LIHI
      IF( PRESENT(ABNRM) ) ABNRM = LABNRM
      
      DEALLOCATE(WORK, STAT=ISTAT1)
500   DEALLOCATE(RWORK, STAT=ISTAT1)
400   IF( .NOT. PRESENT(RCONDV) ) DEALLOCATE( LRCONDV, STAT=ISTAT1 )
300   IF( .NOT. PRESENT(RCONDE) ) DEALLOCATE( LRCONDE, STAT=ISTAT1 )
200   IF( .NOT. PRESENT(SCALE) ) DEALLOCATE( LSCALE, STAT=ISTAT1 )
      IF (.NOT.PRESENT (VL)) DEALLOCATE(LLVL)
      IF (.NOT.PRESENT (VR)) DEALLOCATE(LLVR)
   ENDIF
100   CALL ERINFO(LINFO, SRNAME, INFO, ISTAT)

END SUBROUTINE ZGEEVX_F95
