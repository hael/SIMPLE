SUBROUTINE ZHEEVR_F95( A, W, JOBZ, UPLO, VL, VU, IL, IU, M, ISUPPZ,  ABSTOL, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
      USE LA_PRECISION, ONLY: WP => DP
      USE LA_AUXMOD, ONLY: ERINFO, LSAME
      USE F77_LAPACK, ONLY: HEEVR_F77 => LA_HEEVR
!  .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!  .. SCALAR ARGUMENTS ..
      INTEGER, INTENT(IN), OPTIONAL :: IL, IU
      INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
      REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
!  .. ARRAY ARGUMENTS ..
      COMPLEX(WP), INTENT(INOUT) :: A(:,:)
      REAL(WP), INTENT(OUT) :: W(:)
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: ISUPPZ(:)
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!       LA_SYEVR / LA_HEEVR compute selected eigenvalues and, optionally,
! the corresponding eigenvectors of a real symmetric/complex Hermitian 
! matrix A. Eigenvalues and eigenvectors can be selected by specifying
! either a range of values or a range of indices for the desired 
! eigenvalues.
!       LA_SYEVR and LA_HEEVR use a relatively robust representation 
! (RRR) algorithm. It is usually the fastest algorithm of all and uses 
! the least workspace.
! 
! ========
! 
!        SUBROUTINE LA SYEVR / LA HEEVR( A, W, JOBZ=jobz, &
!                 UPLO=uplo, VL=vl, VU=vu, IL=il, IU=iu, M=m, &
!                 ISUPPZ= isuppz, ABSTOL=abstol, INFO=info )
!           <type>(<wp>), INTENT(INOUT) :: A(:,:)
!           REAL(<wp>), INTENT(OUT) :: W(:)
!           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!           REAL(<wp>), INTENT(IN), OPTIONAL :: VL, VU
!           INTEGER, INTENT(IN), OPTIONAL :: IL, IU
!           INTEGER, INTENT(OUT), OPTIONAL :: M
!           INTEGER, INTENT(OUT), OPTIONAL :: ISUPPZ(:)
!           REAL(<wp>), INTENT(IN), OPTIONAL :: ABSTOL
!           INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!           <type> ::= REAL | COMPLEX
!           <wp>   ::= KIND(1.0) | KIND(1.0D0)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX square array, shape (:,:).
!         On entry, the matrix A.
!         If UPLO = 'U', the upper triangular part of A contains the 
!         upper triangular part of the matrix A. If UPLO = 'L', the 
!         lower triangular part of A contains the lower triangular part
!         of the matrix A.
!         On exit:
!         If JOBZ = 'V', then the first M columns of A contain the 
!         orthonormal eigenvectors of the matrix A corresponding to the
!         selected eigenvalues, with the i-th column of A containing the
!         eigenvector associated with the eigenvalue in W(i).
!         If JOBZ = 'N', the upper triangle (if UPLO = 'U') or the lower 
!         triangle (if UPLO = 'L') of A, including the diagonal, is 
! 	  destroyed.
! W       (output) REAL array, shape (:) with size(W) = size(A,1).
!         The first M elements contain the selected eigenvalues in 
!         ascending order.
! JOBZ    Optional (input) CHARACTER(LEN=1).
!         = 'N': Computes eigenvalues only;
!         = 'V': Computes eigenvalues and eigenvectors.
!         Default value: 'N'.
! UPLO    Optional (input) CHARACTER(LEN=1).
!         = 'U': Upper triangle of A is stored;
!         = 'L': Lower triangle of A is stored.
!         Default value: 'U'.
! VL,VU   Optional (input) REAL.
!         The lower and upper bounds of the interval to be searched for
!         eigenvalues. VL < VU.
!         Default values: VL = -HUGE(<wp>) and VU = HUGE(<wp>), where 
!         <wp> ::= KIND(1.0) | KIND(1.0D0).
!         Note: Neither VL nor VU may be present if IL and/or IU is 
!         present.
! IL,IU   Optional (input) INTEGER.
!         The indices of the smallest and largest eigenvalues to be
!         returned. The IL-th through IU-th eigenvalues will be found.
! 	     1 <= IL <= IU <= size(A,1).
!         Default values: IL = 1 and IU = size(A,1).
!         Note: Neither IL nor IU may be present if VL and/or VU is 
!         present.
!         Note: All eigenvalues are calculated if none of the arguments
!         VL, VU, IL and IU are present.
! M       Optional (output) INTEGER.
!         The total number of eigenvalues found. 0 <= M <= size(A,1).
!         Note: If IL and IU are present then M = IU-IL+1.
! ISUPPZ  Optional (output) INTEGER array, shape (:) with 
!         size(ISUPPZ) = 2*max(1,M).
!         The support of the eigenvectors in A, i.e., the indices 
!         indicating the nonzero elements. The i-th eigenvector is
!         nonzero only in elements ISUPPZ(2i-1) through ISUPPZ(2i) .
!         Note: ISUPPZ must be absent if JOBZ = 'N'.
! ABSTOL  Optional (input) REAL.
!         The absolute error tolerance for the eigenvalues. An 
!         approximate eigenvalue is accepted as converged when it is 
!         determined to lie in an interval [a, b] of width less than or
!         equal to ABSTOL + EPSILON(1.0_<wp>) * max(|a|,|b|),
!         where <wp> is the working precision. If ABSTOL <= 0, then 
! 	  EPSILON(1.0_<wp>)* ||T||1 will be used in its place, where 
! 	  ||T||1 is the l1 norm of the tridiagonal matrix obtained by 
!         reducing A to tridiagonal form.
!         Default value: 0.0_<wp>.
!         Note: Eigenvalues are computed most accurately if ABSTOL is set
!         to LA_LAMCH(1.0_<wp>,'Safe minimum'), not zero.
! INFO    Optional (output) INTEGER
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         If INFO is not present and an error occurs, then the program
!         is terminated with an error message.
!----------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
      CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_HEEVR'
!  .. LOCAL SCALARS ..
      CHARACTER(LEN=1) :: LJOBZ, LUPLO, LRANGE
      INTEGER :: N, LINFO, LD, LDZ, LZ, LIL, LIU, LM, LWORK, ISTAT, &
     &  LIWORK, SISUPPZ, NN, LRWORK
      INTEGER:: IWORKMIN(1)
      COMPLEX(WP) :: WORKMIN(1)
      REAL(WP), POINTER :: RWORK(:)
      REAL(WP) :: RWORKMIN(1) 
      REAL(WP) :: LABSTOL, LVL, LVU
!  .. LOCAL ARRAYS ..
      INTEGER, POINTER :: IWORK(:), LISUPPZ(:)
      COMPLEX(WP), POINTER :: WORK(:), Z(:,:)
!  .. INTRINSIC FUNCTIONS ..
      INTRINSIC HUGE, PRESENT, SIZE
      INTEGER ILAENV
      EXTERNAL ILAENV
!  .. EXECUTABLE STATEMENTS ..
      N = SIZE(A,1); LD = MAX(1,N); LINFO = 0; ISTAT = 0
      NN=2*MAX(1,N)
      IF( PRESENT(JOBZ) )THEN; LJOBZ = JOBZ; ELSE; LJOBZ = 'N'; ENDIF
      IF( PRESENT(ISUPPZ) )THEN; SISUPPZ = SIZE(ISUPPZ); ELSE; SISUPPZ = NN; END IF 
      IF( PRESENT(UPLO) ) THEN; LUPLO = UPLO; ELSE; LUPLO = 'U';  END IF
      IF( PRESENT(VL) )THEN; LVL = VL; ELSE; LVL = -HUGE(LVL);  ENDIF
      IF( PRESENT(VU) )THEN; LVU = VU; ELSE; LVU = HUGE(LVU); ENDIF
      IF( PRESENT(IL) )THEN; LIL = IL; ELSE; LIL = 1; ENDIF
      IF( PRESENT(IU) )THEN; LIU = IU; ELSE; LIU = N; ENDIF
!  .. TEST THE ARGUMENTS ..
        IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN; LINFO = -1
        ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
	ELSE IF( .NOT.LSAME(LJOBZ,'N') .AND. .NOT.LSAME(LJOBZ,'V') )THEN; LINFO = -3
        ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN; LINFO = -4
        ELSE IF( LVU < LVL )THEN; LINFO = -5
        ELSE IF( (PRESENT(VL) .OR. PRESENT(VU)) .AND. (PRESENT(IL) .OR. PRESENT(IU)) )THEN
          LINFO = -6
        ELSE IF(( LIU < LIL .OR. LIL < 1) .AND. N>0)THEN; LINFO = -7
        ELSE IF( N < LIU )THEN; LINFO = -8
        ELSE IF( SISUPPZ /= NN .OR. PRESENT(ISUPPZ).AND.LSAME(LJOBZ,'N'))THEN; LINFO = -10
        ELSE IF( N > 0 )THEN
          IF( PRESENT(VL) .OR. PRESENT(VU) )THEN
            LRANGE = 'V'; LM = N
          ELSE IF( PRESENT(IL) .OR. PRESENT(IU) )THEN
            LRANGE = 'I'; LM = LIU-LIL+1
          ELSE
            LRANGE = 'A'; LM = N
          END IF
	  IF ( LSAME(LJOBZ,'V') ) THEN
	      LDZ=N; LZ=LM
	  ELSE
	      LDZ=1; LZ=1
	  ENDIF
	  ALLOCATE (Z(LDZ,LZ), STAT=ISTAT)							
!  .. QUERING THE SIZE OF WORKSPACE ..
             LWORK = -1
	     LRWORK = -1
	     LIWORK = -1
	        CALL HEEVR_F77( LJOBZ, LRANGE, LUPLO, N, A, LD, LVL, LVU, &
     &             LIL, LIU, LABSTOL, LM, W, Z, LDZ, IWORKMIN, &
     &             WORKMIN, LWORK, RWORKMIN, LRWORK, IWORKMIN, LIWORK, LINFO )
! THE NEXT  LINE SHOULD BE
!           LWORK = WORKMIN(1) WHEN  THE LAPACK IS FIXED
            LWORK = 10*WORKMIN(1)+10000
	    LRWORK = RWORKMIN(1)
	    LIWORK = IWORKMIN(1)
	    DEALLOCATE(Z)
	    
	    ALLOCATE(IWORK(LIWORK), RWORK(LRWORK), STAT=ISTAT)
            IF (ISTAT /=0) THEN ; LINFO = -100; GOTO 100; ENDIF
            ALLOCATE(WORK(LWORK), STAT=ISTAT)
            IF (ISTAT /=0) THEN ; LINFO = -100; GOTO 200; ENDIF
                
            ALLOCATE(LISUPPZ(NN), STAT=ISTAT)
            IF (ISTAT /= 0) THEN ; LINFO = -100; GOTO 300; ENDIF
	    
	    ALLOCATE(Z(LDZ, LZ), STAT = ISTAT)
	    IF(ISTAT /= 0) THEN; LINFO = -100; GOTO 400; ENDIF
                  
            IF( PRESENT(ABSTOL) )THEN; LABSTOL = ABSTOL; ELSE
!              LABSTOL = 2*LAMCH_F77('Safe minimum');  ENDIF
              LABSTOL = 0.0_WP;  ENDIF
! .. CALL LAPACK77 ROUTINE ..
              CALL HEEVR_F77( LJOBZ, LRANGE, LUPLO, N, A, LD, LVL, LVU, &
     &          LIL, LIU, LABSTOL, LM, W, Z, LDZ, LISUPPZ, &
     &          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, LINFO )
              IF( LSAME(LJOBZ,'V') ) A(1:LDZ,1:LM) = Z(1:LDZ,1:LM)
              IF( PRESENT(M) ) M = LM
	      IF(PRESENT(ISUPPZ)) ISUPPZ = LISUPPZ

              DEALLOCATE(Z)
400           DEALLOCATE(LISUPPZ)
300           DEALLOCATE(WORK)
200           DEALLOCATE(IWORK, RWORK)
         ENDIF
100         CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZHEEVR_F95
