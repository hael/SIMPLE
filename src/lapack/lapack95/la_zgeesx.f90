SUBROUTINE ZGEESX_F95( A, W, VS, SELECT, SDIM, RCONDE, RCONDV, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD, ONLY: ERINFO, LSAME
   USE F77_LAPACK, ONLY: GEESX_F77 => LA_GEESX
!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL :: INFO, SDIM
   REAL(WP), INTENT(OUT), OPTIONAL :: RCONDE, RCONDV
!  .. ARRAY ARGUMENTS ..
   COMPLEX(WP), INTENT(INOUT) :: A(:,:)
   COMPLEX(WP), INTENT(OUT) :: W(:)
   COMPLEX(WP), INTENT(OUT), OPTIONAL, TARGET :: VS(:,:)
!  .. EXTERNAL ARGUMENTS ..
      INTERFACE
         FUNCTION SELECT(W)
            USE LA_PRECISION, ONLY: WP => DP
            COMPLEX(WP), INTENT(IN) :: W
            LOGICAL :: SELECT
         END FUNCTION SELECT
      END INTERFACE
   OPTIONAL :: SELECT
!----------------------------------------------------------------------
! 
! Purpose
! =======
! 
!     LA_GEESX computes for a real/complex square matrix A, the 
! eigenvalues, the real-Schur/complex-Schur form T , and, optionally, the
! matrix of Schur vectors Z, where Z is orthogonal/unitary. This gives the
! factorization
!                     A = Z*T*Z^H.
! Optionally, it also orders the eigenvalues on the diagonal of the Schur 
! form so that selected eigenvalues are at the top left, computes a 
! reciprocal condition number for the average of the selected eigenvalues
! and computes a reciprocal condition number for the right invariant 
! subspace corresponding to the selected eigenvalues. The leading columns
! of Z form an orthonormal basis for this invariant subspace.
! A real matrix is in real-Schur form if it is block upper triangular
! with 1 by 1 and 2 by 2 blocks along the main diagonal. 2 by 2 blocks 
! are standardized in the form
!                        [ a  b ]
! 	                 [ c  a ]
! where b*c < 0. The eigenvalues of such a block are a +/- Sqrt(b*c).
! 
! =========
! 
!       SUBROUTINE LA_GEESX( A, <w>, VS=vs, SELECT=select, SDIM=sdim, &
!                              RCONDE=rconde, RCONDV=rcondv, INFO=info )
!            <type>(<wp>), INTENT(INOUT) :: A(:,:)
!            <type>(<wp>), INTENT(OUT) :: <w(:)>
!            <type>(<wp>), INTENT(OUT), OPTIONAL :: VS(:,:)
!            INTERFACE
!               LOGICAL FUNCTION SELECT(<w(j)>)
!                  <type>(<wp>), INTENT(IN) :: <w(j)>
!               END FUNCTION SELECT
!            END INTERFACE
!            OPTIONAL :: SELECT
!            INTEGER, INTENT(OUT), OPTIONAL :: SDIM
!            REAL(<wp>), INTENT(OUT), OPTIONAL :: RCONDE, RCONDV
!            INTEGER, INTENT(OUT), OPTIONAL :: INFO
!       where
!            <type> ::= REAL | COMPLEX
!            <wp>   ::= KIND(1.0) | KIND(1.0D0)
!            <w>    ::= WR, WI | W
!            <w(:)> ::= WR(:), WI(:) | W(:)
!            <w(j)> ::= WR(j), WI(j) | W(j)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX square array, shape (:,:).
!         On entry, the matrix A.
!         On exit, the Schur form T .
! <w>     (output) REAL or COMPLEX array, shape (:) with size(w) = 
!         size(A,1).
!         The computed eigenvalues in the order in which they appear on
!         the diagonal of the Schur form T.
!         <w(:)> ::= WR(:), WI(:) | W(:),
!         where
!         WR(:), WI(:) are of REAL type (for the real and imaginary 
!         parts) and W(:) is of COMPLEX type.
!         Note: If A is real, then a complex-conjugate pair appear 
!         consecutively, with the eigenvalue having the positive 
! 	  imaginary part appearing first.
! VS      Optional (output) REAL or COMPLEX square array, shape (:,:) 
!         with size(VS,1) = size(A,1).
!         The matrix Z of Schur vectors.
! SELECT  Optional (input) LOGICAL FUNCTION
!         LOGICAL FUNCTION SELECT( <w(j)>)
!            <type>(<wp>), INTENT(IN) :: <w(j)>
!         where
!            <type> ::= REAL | COMPLEX
!            <wp> ::= KIND(1.0) | KIND(1.0D0)
!            <w(j)> ::= WR(j), WI(j) | W(j)
!         1. SELECT must be declared as EXTERNAL or as an explicit 
! 	   interface in the calling (sub)program.
!         2. SELECT is called by LA_GEES for every computed eigenvalue
! 	   <w(j)> (but only once for a complex conjugate pair when A 
! 	   is real). It is used to select the eigenvalues that will be
!          ordered to the top left of the Schur form. The eigenvalue 
! 	   <w(j)> is selected if SELECT(<w(j)>) has the value .TRUE.
!         3. A selected complex eigenvalue may no longer satisfy 
! 	   SELECT(<w(j)>) = .TRUE. after ordering, since ordering may
!          change the value of complex eigenvalues (especially if the
!          eigenvalue is ill-conditioned). In this case INFO is set to
!          size(A,1) + 2 (see INFO below).
!         Note: Select must be present if SDIM, RCONDE and RCONDF are 
! 	  desired.
! SDIM    Optional (output) INTEGER.
!         The number of eigenvalues (after sorting) for which 
! 	  SELECT = .TRUE. (If A is real, complex conjugate pairs for which
! 	  SELECT = .TRUE. for either eigenvalue count as 2).
! RCONDE  Optional (output) REAL.
!         The reciprocal condition number for the average of the selected
! 	  eigenvalues.
! RCONDV  Optional (output) REAL.
!         The reciprocal condition number for the selected right 
! 	  invariant subspace.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: if INFO = i, and i is
!              <= n: the QR algorithm failed to compute all the 
!                    eigenvalues; elements 1:ilo-1 and i+1:n of <w> 
!                    contain those eigenvalues which have converged. 
!                    VS contains the matrix which reduces A to its 
!                    partially converged Schur form.
!              = n+1: the eigenvalues could not be reordered because some
! 	             eigenvalues were not sufficiently separated (the 
!                    problem is very ill-conditioned).
!              = n+2: after reordering, some leading complex eigenvalues
!                    in the Schur form no longer satisfy SELECT = .TRUE. 
!                    This can be caused by ordinary roundoff or underflow
!                    due to scaling.
!              n is the order of A.
!         If INFO is not present and an error occurs, then the program is
!         terminated with an error message.
!------------------------------------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GEESX'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1) :: LJOBVS, LSORT, LSENSE
   INTEGER, SAVE :: LWORK = 0
   INTEGER :: N, LINFO, LD, ISTAT, ISTAT1, S1VS, S2VS, LSDIM
   REAL(WP) :: LRCONDE, LRCONDV
!  .. LOCAL ARRAYS ..
   REAL(WP), POINTER :: RWORK(:)
   LOGICAL, TARGET :: LLBWORK(1)
   LOGICAL, POINTER :: BWORK(:)
   COMPLEX(WP), TARGET :: LLVS(1,1)
   COMPLEX(WP), POINTER :: WORK(:)
!  .. INTRINSIC FUNCTIONS ..
   INTRINSIC MAX, PRESENT, SIZE
!  .. EXECUTABLE STATEMENTS ..
   LINFO = 0; ISTAT = 0; N = SIZE(A,1); LD = MAX(1,N)
   IF( PRESENT(VS) )THEN; S1VS = SIZE(VS,1); S2VS = SIZE(VS,2); LJOBVS = 'V'
   ELSE; S1VS = 1; S2VS = 1; LJOBVS = 'N'; END IF
   IF( PRESENT(SDIM) .OR. PRESENT(SELECT) )THEN; LSORT = 'S'; ELSE; LSORT = 'N'; END IF
   IF( PRESENT(RCONDE).AND.PRESENT(RCONDV) )THEN; LSENSE = 'B'
   ELSE IF( PRESENT(RCONDE) )THEN; LSENSE = 'E'
   ELSE IF( PRESENT(RCONDV) )THEN; LSENSE = 'V'; ELSE; LSENSE = 'N'; ENDIF
   IF( .NOT. LSAME(LSENSE,'N') ) LSORT = 'S'
!  .. TEST THE ARGUMENTS
   IF( N < 0 .OR. SIZE(A,2) /= N )THEN; LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN; LINFO = -2
   ELSE IF( PRESENT(VS) .AND. ( S1VS /= N .OR. S2VS /= N ) )THEN; LINFO = -3
   ELSE IF( LSAME(LSORT,'S') .AND. .NOT.PRESENT(SELECT) )THEN; LINFO = -4
   ELSE IF( N > 0 )THEN
      IF( LSAME(LSORT,'S') )THEN
        ALLOCATE(BWORK(N),STAT=ISTAT)
        IF (ISTAT /= 0) THEN; LINFO = -100; GOTO 100; ENDIF
      ELSE; BWORK => LLBWORK; END IF

      ALLOCATE( RWORK(MAX(1,N)), STAT=ISTAT )
      IF (ISTAT /= 0) THEN; LINFO = -100; GOTO 200; ENDIF

      LWORK = MAX(1, 2*N, N*N/2)
      ALLOCATE(WORK(LWORK), STAT=ISTAT)
      IF (ISTAT /= 0) THEN; LINFO = -100; GOTO 300; ENDIF
      IF (PRESENT (VS)) THEN
          CALL GEESX_F77( LJOBVS, LSORT, SELECT, LSENSE, N, A, LD, LSDIM, W, &
     &        VS, S1VS, LRCONDE, LRCONDV, WORK, LWORK, &
     &        RWORK, BWORK, LINFO )
       ELSE
          CALL GEESX_F77( LJOBVS, LSORT, SELECT, LSENSE, N, A, LD, LSDIM, W, &
     &        LLVS, S1VS, LRCONDE, LRCONDV, WORK, LWORK, &
     &        RWORK, BWORK, LINFO )
      ENDIF  	
      IF( LINFO == 0 ) LWORK = INT(WORK(1)+1)

      IF( PRESENT(SDIM) ) SDIM = LSDIM
      IF( PRESENT(RCONDE) ) RCONDE = LRCONDE
      IF( PRESENT(RCONDV) ) RCONDV = LRCONDV
      IF( LSAME(LSORT,'S') ) DEALLOCATE(BWORK, STAT=ISTAT1)

      DEALLOCATE(WORK, STAT=ISTAT)
300   DEALLOCATE(RWORK, STAT=ISTAT )
200   IF( LSAME(LSORT,'S')) DEALLOCATE(BWORK,STAT=ISTAT)
    ENDIF
100 CALL ERINFO(LINFO,SRNAME,INFO,ISTAT)
END SUBROUTINE ZGEESX_F95
