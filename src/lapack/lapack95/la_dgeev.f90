SUBROUTINE DGEEV_F95( A, WR, WI, VL, VR, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: GEEV_F77 => LA_GEEV
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(INOUT) :: A(:,:)
   REAL(WP), INTENT(OUT) :: WR(:), WI(:)
   REAL(WP), INTENT(OUT), OPTIONAL, TARGET :: VL(:,:), VR(:,:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!        LA_GEEV computes for a real or complex square matrix A, the 
! eigenvalues and, optionally, the left and/or right eigenvectors. A 
! right eigenvector v(j) of A satisfies
!                   A * v(j) = lambda(j) * v(j)
! where lambda(j) is its eigenvalue. A left eigenvector u(j) of A 
! satisffies
!                   u(j)^H * A = lambda(j) * u(j)^H
! where u(j)^H denotes the conjugate-transpose of u(j).
! 
! =========
! 
!       SUBROUTINE LA_GEEV( A, <w>, VL=vl, VR=vr, INFO=info )
!            <type>(<wp>), INTENT(INOUT) :: A(:,:)
!            <type>(<wp>), INTENT(OUT) :: <w(:)>
!            <type>(<wp>), INTENT(OUT), OPTIONAL :: VL(:,:), VR(:,:)
!            INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!            <type> ::= REAL | COMPLEX
!            <wp>   ::= KIND(1.0) | KIND(1.0D0)
!            <w>    ::= WR, WI | W
!            <w(:)> ::= WR(:), WI(:) | W(:)
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
! INFO     Optional (output) INTEGER.
!          = 0: successful exit.
!          < 0: if INFO = -i, the i-th argument had an illegal value.
!          > 0: if INFO = i, the QR algorithm failed to compute all the 
!          eigenvalues and no eigenvectors were computed. Elements 
! 	   i+1 : n of <w> contain eigenvalues which have converged.
!          n is the order of A
!          If INFO is not present and an error occurs, then the program
!          is terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GEEV'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBVL, LJOBVR
   INTEGER, SAVE :: LWORK = 0
   INTEGER :: N, NN, LINFO, LD, ISTAT, ISTAT1, S1VL, S2VL, S1VR, S2VR
!  .. LOCAL ARRAYS ..
   REAL(WP), TARGET :: LLVL(1,1), LLVR(1,1)
   REAL(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
   IF( PRESENT(VL) )THEN; S1VL = SIZE(VL,1); S2VL = SIZE(VL,2); LJOBVL = 'V'
   ELSE; S1VL = 1; S2VL = 1; LJOBVL = 'N'; END IF
   IF( PRESENT(VR) )THEN; S1VR = SIZE(VR,1); S2VR = SIZE(VR,2); LJOBVR = 'V'
   ELSE; S1VR = 1; S2VR = 1; LJOBVR = 'N'; END IF
!  .. TEST THE ARGUMENTS
   IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
   ELSE IF( SIZE( WR ) /= N )THEN; LINFO = -2
   ELSE IF( SIZE( WI ) /= N )THEN; LINFO = -3
   ELSE IF( PRESENT(VL) .AND. ( S1VL /= N .OR. S2VL /= N ) )THEN; LINFO = -4
   ELSE IF( PRESENT(VR) .AND. ( S1VR /= N .OR. S2VR /= N ) )THEN; LINFO = -5
   ELSE IF( N > 0 )THEN
      NN = 3; IF( LSAME(LJOBVL,'V').OR.LSAME(LJOBVR,'V') ) NN = NN + 1
      LWORK = MAX( 1, NN*N, LWORK); ALLOCATE(WORK(LWORK), STAT=ISTAT)
         IF( ISTAT /= 0 )THEN; DEALLOCATE(WORK,STAT=ISTAT1)
            LWORK = MAX( 1, NN*N ); ALLOCATE(WORK(LWORK), STAT=ISTAT)
            IF( ISTAT == 0) CALL ERINFO( -200, SRNAME, LINFO )
         END IF
      IF( ISTAT == 0 ) THEN
        IF( PRESENT(VL) )THEN
           IF( PRESENT(VR) )THEN
               CALL GEEV_F77( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                        VL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
	   ELSE
	      CALL GEEV_F77( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                        VL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
	   ENDIF
	 ELSE
	   IF( PRESENT(VR) )THEN
	       CALL GEEV_F77( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                        LLVL, S1VL, VR, S1VR, WORK, LWORK, LINFO )
	   ELSE
	       CALL GEEV_F77( LJOBVL, LJOBVR, N, A, LD, WR, WI, &
                        LLVL, S1VL, LLVR, S1VR, WORK, LWORK, LINFO )
	   ENDIF
	 ENDIF  
         IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)
      ELSE; LINFO = -100; ENDIF
      DEALLOCATE(WORK, STAT=ISTAT1)
   ENDIF
   CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE DGEEV_F95
